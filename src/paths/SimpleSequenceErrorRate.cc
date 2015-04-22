///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SimpleSequenceErrorRate.  From high-coverage data, compute a de novo estimate for
// the per-base indel error rate in sequencing M^n, where M is a motif.  This is 
// done by grouping reads of the form F1 M^n F2, where F1 and F2 are ten-base flanks
// and n is allowed to vary within the group.  
//
// One of the challenges in making this computation is that distinct loci on the 
// genome may be combined ("overcollapsed") within one group.  An approach to 
// reducing this effect would be to use the read locations on an assembly to break 
// apart read groups into subgroups that lie in the same position on the assembly.  
// In this code we do two things to mitigate the effect of overcollapsed loci:
// first we screen out loci having unusually high coverage; second, we calculate
// the error rate parameter by computing a Poisson parameter from the fraction of
// the observed distribution which occurs with low error rate.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "math/Functions.h"
#include "system/SortInPlace.h"

// Horrors.  A char20 is really 20 bases.

class char20 {
     public:
     char x[20];
     friend int compare( char20 const& i1, char20 const& i2 )
     { return memcmp(i1.x,i2.x,5); }
     friend Bool operator==( const char20& I1, const char20& I2 )
     {    return memcmp( I1.x, I2.x, 5 ) == 0;    }
     friend Bool operator!=( const char20& I1, const char20& I2 )
     {    return memcmp( I1.x, I2.x, 5 ) != 0;    }
     friend Bool operator<( const char20& I1, const char20& I2 )
     {    return memcmp( I1.x, I2.x, 5 ) < 0;    }
     friend Bool operator>( const char20& I1, const char20& I2 )
     {    return memcmp( I1.x, I2.x, 5 ) > 0;    }
     void GetFromChars( const vec<char>& b )
     {    for ( int i = 0; i < 5; i++ )
               x[i] = 0;
          for ( int i = 0; i < 20; i++ )
          {    if ( b[i] == 'A' );
               else if ( b[i] == 'C' ) x[i/4] ^= ( 1 << 2*(i%4) );
               else if ( b[i] == 'G' ) x[i/4] ^= ( 2 << 2*(i%4) );
               else if ( b[i] == 'T' ) x[i/4] ^= ( 3 << 2*(i%4) );
               else
               {    cout << "GetFromChars: illegal character " << b[i] << endl;
                    cout << "Abort." << endl;
                    exit(1);    }    }    }
               
};

// For a sequence F1 M^n F2, a quadx consists of F1F2 (concatenated) and n, together
// with a flag to indicate whether the data are from reads or from the genome.

class quadx {

     public:

     quadx( ) : third(0), fourth(False) { }

     quadx( const char20& first_second, const int third, const Bool fourth ) :
          first_second(first_second), third(third), fourth(fourth) { }

     char20 first_second;
     int third;
     Bool fourth;

     friend int compare( quadx const& q1, quadx const& q2 )
     { int result = compare(q1.first_second,q2.first_second);
       if ( !result ) result = compare(q1.third,q2.third);
       if ( !result ) result = compare(q1.fourth,q2.fourth);
       return result; }
     friend Bool operator<( const quadx& q1, const quadx& q2 )
     {    if ( q1.first_second < q2.first_second ) return True;
          if ( q1.first_second > q2.first_second ) return False;
          if ( q1.third < q2.third ) return True;
          if ( q1.third > q2.third ) return False;
          if ( q1.fourth < q2.fourth ) return True;
          return False;    }

};

String Rotate( const String s, const int r )
{    String x( s.size( ) );
     for ( int i = 0; i < s.isize( ); i++ )
          x[i] = s[ (i+r) % s.isize( ) ];
     return x;    }

double PredictedLowFrac( const long double p, const double E, 
     const vec<int>& count, const int winner )
{    int ngroups = count.size( );
     long double predicted_low_count = 0.0;
     for ( int l = 0; l < ngroups; l++ )
     {    int c = count[l];
          long double lambda = p * double(c) * double(winner);
          int max_e = int( floor( E * double(c) ) );
          long double next = exp(-lambda);
          for ( int j = 0; j <= max_e; j++ )
          {    predicted_low_count += next;
               next *= lambda / double(j+1);    }    }
     return predicted_low_count/double(ngroups);    }

int main( int argc, char *argv[] )
{
  
     RunTime( );
  
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int_OrDefault(VERBOSITY, 0 );
     CommandArgument_Int_OrDefault(MAX_MOTIF, 3 );
     CommandArgument_Bool_OrDefault(SEE_GENOME, False);
     CommandArgument_Bool_OrDefault(REQUIRE_TRUE, False);
     CommandArgument_String_OrDefault(OUT, "");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
           "Number of threads for sorting, defaults to number of processors.");
     EndCommandArguments;

     // Define directories.
  
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;

     // Set up loggin and load data.

     ostream* out = ( OUT == "" ? &cout : new ofstream( OUT.c_str( ) ) );
     if ( VERBOSITY >= 1 ) cout << Date( ) << ": loading reads" << endl;
     vecbasevector reads( run_dir + "/frag_reads_filt.fastb" );
     vecbasevector genome;
     if (SEE_GENOME)
     {    if ( VERBOSITY >= 1 ) cout << Date( ) << ": loading genome" << endl;
          // genome.ReadAll( data_dir + "/../genome_extended.fastb" );
          // FetchReads( genome, 0, data_dir + "/../genome_extended.fasta" );
          FetchReads( genome, 0, data_dir + "/../genome.fasta" );    }

     // Define heuristic constants.

     const int flank = 10; // hardwired via char20 structure
     const int min_h = 4;
     const int max_h = 30;
     const int min_to_report = 10;

     // Define motifs to process.

     if ( VERBOSITY >= 1 ) cout << Date( ) << ": finding motifs" << endl;
     vec<String> motifs;
     for ( int k = 1; k <= MAX_MOTIF; k++ )
     {    vec<String> M;
          String m(k);
          for ( int z = 0; z < IPow(4, k); z++ )
          {    int zz = z;
               for ( int j = 0; j < k; j++ )
               {    m[j] = as_base( zz % 4 );
                    zz /= 4;    }
               Bool found = False;
               for ( int r = 0; r < k; r++ )
               {    String mrotated = Rotate( m, r ), mrotated_rc;
                    if ( Member( M, mrotated ) ) found = True;
                    StringReverseComplement( mrotated, mrotated_rc );
                    if ( Member( M, mrotated_rc ) ) found = True;    }
               if (found) continue;
               Bool is_copy_of_smaller = False;
               for ( int div = 1; div < k; div++ )
               {    if ( k % div != 0 ) continue;
                    Bool is_copy_of_div = True;
                    for ( int l = 0; l < k/div; l++ )
                    {    for ( int i = 0; i < div; i++ )
                         {    if ( m[i] != m[ i + l*div ] )
                                   is_copy_of_div = False;    }    }
                    if (is_copy_of_div) is_copy_of_smaller = True;    }
               if (is_copy_of_smaller) continue;
               M.push_back(m);    }
          motifs.append(M);    }

     // Go through the motifs.
     typedef vec<quadx>::iterator ItrQ;
     typedef CompareFunctor<quadx> CompQ;
     InPlaceParallelSorter<ItrQ,CompQ > ippsQuad(NUM_THREADS);
     typedef vec< triple<int,int,int> >::iterator ItrT;
     typedef CompareFunctor<triple<int,int,int> > CompT;
     InPlaceParallelSorter<ItrT,CompT> ippsTrip(NUM_THREADS);

     for ( int mi = 0; mi < motifs.isize( ); mi++ )
     {    String motif = motifs[mi], motif_rc;
          int nm = motif.size( );
          StringReverseComplement( motif, motif_rc );
          Bool eq_rc_rotated = False;
          for ( int i = 0; i < nm; i++ )
               if ( motif_rc == Rotate(motif, i) ) eq_rc_rotated = True;
          if ( VERBOSITY >= 1 ) cout << "\nsearching for " << motif << "\n" << endl;
          vec<quadx> F;

          // Look for runs of the motif in the reads that have ten-base flanks that
          // are not 80% or more similar to the motif.

          if ( VERBOSITY >= 1 ) cout << Date( ) << ": creating Fs" << endl;
          const double max_match_rate = 0.8;
          const size_t nbatches = 5 * omp_get_max_threads( );
          vec< vec<quadx> > Fs(nbatches);
          int sources = ( SEE_GENOME ? 2 : 1 );
          for ( int source = 0; source < sources; source++ )
          {    const vecbasevector& S = ( source == 0 ? reads : genome );
               int npasses = ( eq_rc_rotated ? 1 : 2 );
               for ( int pass = 1; pass <= npasses; pass++ )
               {    String motifx = ( pass == 1 ? motif : motif_rc );
                    vec<char> motifc( motifx.size( ) );
                    for ( int j = 0; j < nm; j++ )
                         motifc[j] = as_char( motifx[j] );
                    #pragma omp parallel for
                    for ( size_t u = 0; u < nbatches; u++ )
                    {    basevector x, y;
                         vec<char> c;
                         char20 xy;
                         size_t start = ( u * S.size( ) ) / nbatches;
                         size_t stop = ( (u+1) * S.size( ) ) / nbatches;
                         for ( size_t id = start; id < stop; id++ )
                         {    const basevector& r = S[id];
                              for (int j = flank; j <= r.isize( ) - nm - flank; j++)
                              {    int R;
                                   for ( R = 0; R < nm; R++ )
                                   {    Bool at_motif_R = True;
                                        for ( int l = 0; l < nm; l++ )
                                        {    if ( r[j+l] != motifc[(l+R) % nm] )
                                             {    at_motif_R = False;
                                                  break;    }    }
                                        if (at_motif_R) break;    }
                                   if ( R == nm ) continue;
                                   if ( r[j-1] == motifc[(nm+R-1) % nm ] ) continue;
                                   int k;
                                   for ( k = j + nm; k < r.isize( ); k++ )
                                        if ( r[k] != motifc[(k-j+R)%nm] ) break;
                                   int nh = k - j;
                                   if ( nh >= min_h*nm && nh <= max_h*nm 
                                        && k + flank <= r.isize( ) )
                                   {    x.SetToSubOf( r, j - flank, flank );
                                        y.SetToSubOf( r, k, flank );
                                        Bool simple = False;
                                        for ( int bpass = 1; bpass <= 2; bpass++ )
                                        {    const basevector& b 
                                                  = ( bpass == 1 ? x : y );
                                             int o = ( bpass == 1 
                                                  ? nm-(flank % nm) : k );
                                             int match = 0;
                                             for ( int z = 0; z < flank; z++ )
                                             {    if ( b[z] == motifc[(z+o)%nm] )
                                                  {    match++;    }    }
                                             if ( double(match)/double( flank ) 
                                                  >= max_match_rate )
                                             {    simple = True;    }    }
                                        if ( !simple )
                                        {    if ( pass == 2 )
                                             {    x.ReverseComplement( ); 
                                                  y.ReverseComplement( );
                                                  swap( x, y );    }
                                             c.clear( );
                                             for ( int l = 0; l < flank; l++ )
                                                  c.push_back( as_base( x[l] ) );
                                             for ( int l = 0; l < flank; l++ )
                                                  c.push_back( as_base( y[l] ) );
                                             xy.GetFromChars(c);
                                             Fs[u].push( xy, nh, source );    }    }
                                   j = k - 1;    }    }    }    }    }
          if ( VERBOSITY >= 1 ) cout << Date( ) << ": merging Fs" << endl;
          size_t Fsize = 0;
          for ( size_t i = 0; i < Fs.size( ); i++ )
               Fsize += Fs[i].size( );
          F.resize(Fsize);
          #pragma omp parallel for
          for ( size_t i = 0; i < Fs.size( ); i++ )
          {    size_t start = 0;
               for ( size_t j = 0; j < i; j++ )
                    start += Fs[j].size( );
               for ( size_t j = 0; j < Fs[i].size( ); j++ )
                    F[ start + j ] = Fs[i][j];    }
          if ( VERBOSITY >= 1 ) cout << Date( ) << ": sorting F" << endl;
          ippsQuad.sort( F.begin( ), F.end( ) );

          // Now identify groups.  There are two passes.  The first pass just
          // computes the coverage distribution.

          if ( VERBOSITY >= 1 ) cout << Date( ) << ": identifying groups" << endl;
          vec<int> total(max_h*nm+1);
          vec< vec<int> > freqs( max_h*nm + 1 );
          vec< triple<int,int,int> > winner_count_errs;
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 )
               {    for ( int x = min_h*nm; x <= max_h*nm; x++ )
                         Sort( freqs[x] );    }
               for ( size_t i = 0; i < F.size( ); i++ )
               {    size_t j;
                    for ( j = i + 1; j < F.size( ); j++ )
                         if ( F[j].first_second != F[i].first_second ) break;
                    int cov = 0;
                    for ( size_t r = i; r < j; r++ )
                         if ( !F[r].fourth ) cov++;
                    if ( cov >= min_to_report )
                    {    int max_freq = 0, max_freq_id = -1;
                         for ( size_t r = i; r < j; r++ )
                         {    size_t s;
                              for ( s = r + 1; s < j; s++ )
                                   if ( F[s].third != F[r].third ) break;
                              int c = 0;
                              for ( size_t l = r; l < s; l++ )
                                   if ( !F[l].fourth ) c++;
                              if ( c > max_freq )
                              {    max_freq = c;
                                   max_freq_id = F[r].third;    }
                              r = s - 1;    }
                         if ( pass == 1 )
                         {    freqs[max_freq_id].push_back(cov);
                              i = j - 1;
                              continue;    }
     
                         // We exclude cases where coverage is in the top 10%.
     
                         int max_cov = 1000000000;
                         if ( pass == 2 && freqs[max_freq_id].nonempty( ) )
                         {    int n = freqs[max_freq_id].size( );
                              max_cov = freqs[max_freq_id][(9*n)/10];    }
                         if ( cov <= max_cov )
                         {    vec<int> trues;
                              for ( size_t r = i; r < j; r++ )
                                   if ( F[r].fourth ) trues.push_back( F[r].third );
                              if (REQUIRE_TRUE)
                              {    if ( !trues.solo( ) || trues[0] != max_freq_id ) 
                                   {    i = j - 1;
                                        continue;    }    }
                              int winner = max_freq_id, errs = 0;
                              for ( size_t r = i; r < j; r++ )
                              {    if ( !F[r].fourth )
                                        errs += Abs( F[r].third - max_freq_id );    }
                              winner_count_errs.push( winner, cov, errs );
                              total[max_freq_id] += cov;
     
                              if ( VERBOSITY >= 2 )
                              {    cout << "motif = " << motif
                                        << ", winner = " << max_freq_id << ", ";
                                   for ( size_t r = i; r < j; r++ )
                                   {    size_t s;
                                        for ( s = r + 1; s < j; s++ )
                                             if ( F[s].third != F[r].third ) break;
                                        int c = 0;
                                        for ( size_t l = r; l < s; l++ )
                                             if ( !F[l].fourth ) c++;
                                        if ( r > i ) cout << " ";
                                        cout << F[r].third << "[" << c << "]";
                                        r = s - 1;    }
                                   for ( size_t r = i; r < j; r++ )
                                   {    if ( F[r].fourth ) 
                                        {    cout << " " << F[r].third 
                                                  << "[true]";    }    }
                                   cout << ", cov = " << cov
                                        << ", max_cov = " << max_cov;
                                   UniqueSort(trues);
                                   if ( ( !trues.solo( ) || trues[0] != max_freq_id )
                                        && max_freq < (int) (j-i) ) 
                                   {    cout << " **********\n";   }
                                   cout << "\n\n";    }    }     }
                    i = j - 1;    }    }
          if ( VERBOSITY >= 1 ) cout << Date( ) << ": analyzing winners" << endl;
          ippsTrip.sort( winner_count_errs.begin(), winner_count_errs.end() );
          vec<double> best_p( max_h*nm+1, -1 );
          for ( int i = 0; i < winner_count_errs.isize( ); i++ )
          {    int winner = winner_count_errs[i].first, j;
               for ( j = i + 1; j < winner_count_errs.isize( ); j++ )
                    if ( winner_count_errs[j].first != winner ) break;
               vec<int> count, errs;
               for ( int r = i; r < j; r++ )
               {    count.push_back( winner_count_errs[r].second );
                    errs.push_back( winner_count_errs[r].third );    }

               // Find the smallest E such that the groups having error rate
               // <= E comprise >= 25% of the data.  We refer to these groups as
               // the 'low bin'.

               const double target_frac = 0.25;
               int ngroups = count.size( );
               const int min_groups = 10;
               if ( ngroups < min_groups )
               {    i = j - 1;
                    continue;    }
               vec<double> rates;
               for ( int l = 0; l < ngroups; l++ )
                    rates.push_back( double( errs[l] ) / double( count[l] ) );
               Sort(rates);
               int Ep = int( ceil( double(ngroups) * target_frac ) );
               if ( Ep == ngroups ) Ep--;
               double E = rates[Ep];
               while( Ep < rates.isize( ) - 1 && rates[Ep+1] == E ) Ep++;
               double low_frac = double(Ep) / double(ngroups);

               // Choose the Poisson parameter that best matches the low bin.
               // This is done by a binary search.

               int start = 1, stop = 10000;
               double div = 100000.0;
               best_p[winner] = 0.0;
               double start_val 
                    = PredictedLowFrac( double(start)/div, E, count, winner );
               double stop_val 
                    = PredictedLowFrac( double(stop)/div, E, count, winner );
               while(1)
               {    if ( start_val <= low_frac )
                    {    best_p[winner] = double(start)/div;
                         break;    }
                    if ( stop_val >= low_frac )
                    {    best_p[winner] = double(stop)/div;
                         break;    }
                    if ( stop - start == 1 )
                    {    if ( start_val - low_frac <= low_frac - stop_val )
                              best_p[winner] = double(start)/div;
                         else best_p[winner] = double(stop)/div;
                         break;    }
                    int mid = ( start + stop ) / 2;
                    double mid_val 
                         = PredictedLowFrac( double(mid)/div, E, count, winner );
                    if ( mid_val <= low_frac ) stop = mid;
                    else start = mid;    }
               i = j - 1;    }

          // Derive summary statistics.

          if ( VERBOSITY >= 1 ) 
               cout << Date( ) << ": generating summary statistics\n" << endl;
          int top = max_h*nm;
          while( top >= 0 && total[top] == 0 ) top--;
          for ( int x = min_h*nm; x <= top; x++ )
          {    double X = double(x) / double(nm);
               double p = 100.0 * best_p[x];
               if ( p >= 0 )
                    *out << motif << " " << X << " " << p << "\n";    }    }    }
