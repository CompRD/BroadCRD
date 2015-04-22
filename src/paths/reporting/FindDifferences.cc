///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FindDifferences.  Progressively modify scaffolds until they agree with the
// reference sequence, recording the changes that are made (as if on the reference).
// Certain changes are made silently.  For example, gaps in the scaffold are filled
// under certain circumstances.
//
// IN PROGRESS!

// TO DO:
//
// 1. Identify efasta choices that improve alignment.
//
// 2. What's up with insertions that end in N?
//
// 3. Allow free scaffolds that have aligning bits.
//
// 4. Implement better handling of max_freq for alignment.  Currently it is set 
// to 1, and as a consequence, after extension, we sometimes align to the wrong
// location.  An alternative approach would be to raise max_freq to a small value
// (even 2 would help), and then choose between alternate choices by locally
// Smith-Watermanning or exploiting the density of seeds.  Human scaffold 0 is a 
// good test case.
//
// 5. For insertions, show coordinate range on reference, and where possible, a
// single coordinate.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Alignment.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "PrintAlignment.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "efasta/EfastaTools.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/reporting/PerfAlign.h"
#include "paths/reporting/Sog.h"
#include "random/Shuffle.h"

// Horrors.  A char100 is really 100 bases.  Copied from AssemblyAccuracy.cc.

class char100 {
     public:
     char x[25];
     friend int compare(  char100 const& i1, char100 const& i2 )
     { return memcmp(i1.x,i2.x,25); }
     friend Bool operator==( const char100& I1, const char100& I2 )
     {    return memcmp( I1.x, I2.x, 25 ) == 0;    }
     friend Bool operator!=( const char100& I1, const char100& I2 )
     {    return memcmp( I1.x, I2.x, 25 ) != 0;    }
     friend Bool operator<( const char100& I1, const char100& I2 )
     {    return memcmp( I1.x, I2.x, 25 ) < 0;    }
     void GetFromChars( const vec<char>& b )
     {    for ( int i = 0; i < 25; i++ )
               x[i] = 0;
          for ( int i = 0; i < 100; i++ )
          {    if ( b[i] == 'A' );
               else if ( b[i] == 'C' ) x[i/4] ^= ( 1 << 2*(i%4) );
               else if ( b[i] == 'G' ) x[i/4] ^= ( 2 << 2*(i%4) );
               else if ( b[i] == 'T' ) x[i/4] ^= ( 3 << 2*(i%4) );
               else
               {    cout << "GetFromChars: illegal character " << b[i] << endl;
                    cout << "Abort." << endl;
                    exit(1);    }    }    }
               
};

typedef char100 FASTAVECTOR;

class kmer_pos_fw_id {

     public:

     kmer_pos_fw_id( ) : pos(0), fw(False), id(0) { }

     kmer_pos_fw_id( const FASTAVECTOR& kmer, const int pos, const Bool fw,
          const size_t id ) : kmer(kmer), pos(pos), fw(fw), id(id) { }

     FASTAVECTOR kmer;
     int pos;
     Bool fw;
     size_t id;

     friend int compare( kmer_pos_fw_id const& x1, kmer_pos_fw_id const& x2 )
     { int result = compare(x1.kmer,x2.kmer);
       if ( !result ) result = compare(x1.id,x2.id);
       if ( !result ) result = compare(x1.pos,x2.pos);
       if ( !result ) result = compare(x1.fw,x2.fw);
       return result; }
     friend Bool operator<( const kmer_pos_fw_id& x1, const kmer_pos_fw_id& x2 )
     {    if ( x1.kmer < x2.kmer ) return True;
          if ( x2.kmer < x1.kmer ) return False;
          if ( x1.id < x2.id ) return True;
          if ( x1.id > x2.id ) return False;
          if ( x1.pos < x2.pos ) return True;
          if ( x1.pos > x2.pos ) return False;
          if ( x1.fw < x2.fw ) return True;
          return False;    }

          // return x1.kmer < x2.kmer;    }

};

void MergeOverlapping( sog& S )
{    for ( int j = 0; j < S.AlignsCount( ) - 1; j++ )
     {    const perf_align& a1 = S.Align(j);
          if ( !a1.Fw( ) ) continue;
          int id1 = a1.Id1( ), id2 = a1.Id2( );
          perf_align& a2 = S.Align(j+1);
          if ( !Comparable( a1, a2 ) ) continue;
          if ( !( a1.Pos1( ) < a2.pos1( ) ) ) continue;
          if ( !( a1.Pos2( ) > a2.pos2( ) ) ) continue;
          Bool is_gap = True;
          for ( int l = a1.Pos1( ); l < a2.pos1( ); l++ )
          {    if ( !S.Gap(l) )
               {    is_gap = False;
                    break;    }    }
          if ( !is_gap ) continue;
          vec<char> replacement;
          int a1Pos1 = a1.Pos1( ), a1Pos2 = a1.Pos2( );
          int left_overlap = a1Pos1 - a2.pos1( );
          int right_overlap = a1Pos2 - a2.pos2( );
          int M = Max( left_overlap, right_overlap );
          left_overlap -= M;
          right_overlap -= M;
          a1Pos1 -= M, a1Pos2 -= M;
          perf_align a( id1, id2, a1.pos1( ), 
               a2.Pos1( ) + left_overlap - right_overlap,
               a1.pos2( ), a2.Pos2( ), a1.Rc( ) );
          if ( !S.ReplaceOK( a1Pos1, a2.pos1( ), j ) ) continue;
          S.Replace( a1Pos1, a2.pos1( ), replacement, j, a );
          j--;    }    }

void MakeAlignsFromKpfi( const int K, const vec<kmer_pos_fw_id>& kpfi,
     const int nref, vec<perf_align>& paligns )
{    const size_t max_freq = 1;
     for ( size_t i = 0; i < kpfi.size( ); i++ )
     {    size_t j;
          for ( j = i + 1; j < kpfi.size( ); j++ )
               if ( kpfi[j].kmer != kpfi[i].kmer ) break;
          size_t r;
          for ( r = i; r < j; r++ ) // index of reference contig
               if ( (int) kpfi[r].id >= nref ) break;
          if ( r - i > max_freq )
          {    i = j - 1;
               continue;    }
          for ( size_t z2 = i; z2 < j; z2++ ) // index of reference contig
          {    if ( (int) kpfi[z2].id >= nref ) continue;
               for ( size_t z1 = i; z1 < j; z1++ ) // index of scaffold
               {    if ( (int) kpfi[z1].id < nref ) continue;
                    Bool rc2 = kpfi[z1].fw ^ kpfi[z2].fw;
                    int id1 = kpfi[z1].id - nref, id2 = kpfi[z2].id;
                    int pos1 = kpfi[z1].pos, pos2 = kpfi[z2].pos;
                    paligns.push( id1, id2, pos1, pos1 + K, 
                         pos2, pos2 + K, rc2 );    }    }
          i = j - 1;    }    }

void CondenseAligns( vec<perf_align>& paligns )
{    vec<Bool> to_remove( paligns.size( ), False );
     for ( size_t i = 0; i < paligns.size( ); i++ )
     {    size_t j;
          for ( j = i + 1; j < paligns.size( ); j++ )
          {    if ( !Comparable( paligns[j], paligns[i] ) ) break;
               if ( paligns[i].Fw( ) )
               {    if ( paligns[j].pos1( ) != paligns[j-1].pos1( ) + 1 ) break;
                    if ( paligns[j].pos2( ) != paligns[j-1].pos2( ) + 1 ) break;
                    }
               else
               {    if ( paligns[j].pos1( ) != paligns[j-1].pos1( ) + 1 ) break;
                    if ( paligns[j].pos2( ) != paligns[j-1].pos2( ) - 1 ) break;
                    }    }
          if ( paligns[i].Fw( ) )
          {    paligns[i].Pos1( ) += j - i - 1;
               paligns[i].Pos2( ) += j - i - 1;    }
          else
          {    paligns[i].Pos1( ) += j - i - 1;
               paligns[i].pos2( ) -= j - i - 1;    }
          for ( size_t k = i + 1; k < j; k++ )
               to_remove[k] = True;
          i = j - 1;    }
     EraseIf( paligns, to_remove );    }

// Track events.

vec<String> event_types;
vec< vec<int> > event_count;
void RecordEvent( const String& tname, const int n )
{    
     #pragma omp critical
     {    static Bool first_call = True;
          if (first_call)
          {    event_types.push_back( "substitution", "insertion", "deletion", 
                    "replacement", "free", "transposition", "hanging" );
               event_count.resize( event_types.size( ), vec<int>(7) );
               event_count[0].resize(1);
               first_call = False;    }
          int t = Position( event_types, tname );
          ForceAssertGe( t, 0 );
          if ( n <= 1 )           event_count[t][0]++;
          else if ( n <= 10 )     event_count[t][1]++;
          else if ( n <= 100 )    event_count[t][2]++;
          else if ( n <= 1000 )   event_count[t][3]++;
          else if ( n <= 10000 )  event_count[t][4]++;
          else if ( n <= 100000 ) event_count[t][5]++;
          else                    event_count[t][6]++;    }    }

// Print alignment blocks.

void PrintBlocks( const vec<int>& scaffx, const vec<sog>& sogs )
{    int aligns_count = 0;
     for ( int si = 0; si < scaffx.isize( ); si++ )
     {    int s = scaffx[si];
          if ( sogs[si].Tigs( ).size( ) == 0 ) continue; // entirely gap, temp?
          cout << "\n";
          sogs[si].Print(s);    
          aligns_count += sogs[si].AlignsCount( );    }
     cout << "\n#aligns = " << aligns_count << "\n\n";    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault(ROOT, "");
     CommandArgument_String_Doc(ASSEMBLY, "assembly, in efasta format - note "
          "that gaps should be indicated by a sequence of Ns (upper case)");
     CommandArgument_String_Doc(GENOME, "looks for genome.{fastb,fastamb}");
     CommandArgument_String(ALIGNS);
     CommandArgument_Bool_OrDefault(CACHE, True);
     CommandArgument_Bool_OrDefault(SHOW_EVENTS, False);
     CommandArgument_String_OrDefault_Doc(DUMP_INSERTIONS, "", "if of the form "
          "{fn,x}, then dump insertions of size at least n to the file fn.");
     CommandArgument_String_OrDefault_Doc(SCAFF_OUT, "", "if specified, generate "
          "fasta file of this name containing the blocks in the final scaffolds");
     CommandArgument_String_OrDefault_Doc(DUMP_FREE, "", "if specified, filename "
          "to dump free scaffolds into");
     CommandArgument_String_OrDefault_Doc(SCAFF, "", "if specified in ParseIntSet "
          "format, after alignment, process only these scaffolds");
     CommandArgument_Bool_OrDefault(PRINT_INTERMEDIATES, False);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
           "number of threads used for sorting, defaults to no. of processors");
     EndCommandArguments;

     // Set parallelization level.
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads(SHOW_EVENTS?1:NUM_THREADS);

     // Parse arguments.

     double clock = WallClockTime( );
     String dump_insertions_fn;
     int dump_insertions_min = 1000000000;
     if ( DUMP_INSERTIONS != "" )
     {    vec<String> x;
          ParseStringSet( DUMP_INSERTIONS, x );
          ForceAssertEq( x.isize( ), 2 );
          dump_insertions_fn = x[0];
          dump_insertions_min = x[1].Int( );    }
     vec<int> scaffx;
     if ( SCAFF != "" ) ParseIntSet( SCAFF, scaffx );

     // Append ROOT.

     if ( ROOT != "" )
     {    ASSEMBLY = ROOT + "/" + ASSEMBLY;
          GENOME = ROOT + "/" + GENOME;
          ALIGNS = ROOT + "/" + ALIGNS;    }

     // Test for existence of required files.

     ForceAssert( IsRegularFile(ASSEMBLY) );
     ForceAssert( IsRegularFile( GENOME + ".fastb" ) );
     ForceAssert( IsRegularFile( GENOME + ".fastamb" ) );

     // Load assembly and flatten it.  Since in this version we don't actually
     // use escaffolds, we can cache to avoid loading it every time.  Note funky
     // file naming.

     cout << Date( ) << ": loading assembly" << endl;
     vecbasevector scaffolds;
     vecbitvector scaffolds_gaps;
     String scaffolds_fn = ALIGNS + ".scaffolds";
     String scaffolds_gaps_fn = ALIGNS + ".scaffolds_gaps";
     if ( CACHE && IsRegularFile(scaffolds_fn) && IsRegularFile(scaffolds_gaps_fn) )
     {    scaffolds.ReadAll(scaffolds_fn);
          scaffolds_gaps.ReadAll(scaffolds_gaps_fn);    }
     else
     {    VecEFasta escaffolds;
          LoadEfastaIntoStrings( ASSEMBLY, escaffolds );
          scaffolds.resize( escaffolds.size( ) );
          scaffolds_gaps.resize( escaffolds.size( ) );
          for ( size_t i = 0; i != escaffolds.size( ); i++ )
               escaffolds[i].FlattenTo( scaffolds[i], scaffolds_gaps[i] );
          scaffolds.WriteAll(scaffolds_fn);
          scaffolds_gaps.WriteAll(scaffolds_gaps_fn);    }
     int nscaff = scaffolds.size( );
     if ( SCAFF == "" )
     {    for ( int j = 0; j < nscaff; j++ )
               scaffx.push_back(j);    }

     // Load bases.

     cout << Date( ) << ": loading genome" << endl;
     vecbasevector all( GENOME + ".fastb" );
     int nref = all.size( );
     vecbitvector all_amb( GENOME + ".fastamb" );
     sog mr_sog;
     mr_sog.SetGenome( &all, &all_amb ); // Not quite right, since we append to all!
     for ( int i = 0; i < nscaff; i++ )
     {    all.push_back_reserve( scaffolds[i] );
          all_amb.push_back_reserve( scaffolds_gaps[i] );    }

     // Find all perfect 100-mer matches between the scaffolds and the reference, 
     // for which the 100-mer has a unique placement on the reference.  Then merge
     // these together to form perf_aligns.

     const int K = 100;
     if ( !CACHE || !IsRegularFile(ALIGNS) )
     {    cout << Date() << ": finding matches between chunks and reference" << endl;
          vec<kmer_pos_fw_id> kpfi;
          vec<size_t> kcount( all.size( ), 0 ), koffset( all.size( ) + 1, 0 );
          vec<int> shuffle;
          Shuffle( all.size( ), shuffle );
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 )
               {    for ( size_t i = 1; i <= all.size( ); i++ )
                         koffset[i] = koffset[i-1] + kcount[i-1];
                    kpfi.resize( koffset.back( ) );
                    for ( size_t i = 0; i < all.size( ); i++ )
                         kcount[i] = 0;    }
               #pragma omp parallel for
               for ( size_t ip = 0; ip < all.size( ); ip++ )
               {    int i = shuffle[ip];
                    if ( all[i].size( ) == 0 ) continue;
                    vec<char> akmer(100);
                    FASTAVECTOR kmer;
                    unsigned int nextamb;
                    for ( nextamb = 0; nextamb < all[i].size( ); nextamb++ )
                         if ( all_amb[i][nextamb] ) break;
                    for ( unsigned int j = 0; j <= all[i].size( ) - K; j++ )
                    {    if ( j + K > nextamb )
                         {    for ( ; nextamb < all[i].size( ); nextamb++ )
                                   if ( !all_amb[i][nextamb] ) break;
                              j = nextamb;
                              if ( j == all[i].size( ) ) break;
                              for ( ; nextamb < all[i].size( ); nextamb++ )
                                   if ( all_amb[i][nextamb] ) break;
                              j--;
                              continue;    }
                         basevector b, brc;
                         if ( pass == 2 )
                         {    b.SetToSubOf( all[i], j, K );
                              brc = b;
                              brc.ReverseComplement( );
                              basevector& ab = ( b < brc ? b : brc );
                              for ( int l = 0; l < K; l++ )
                                   akmer[l] = as_base( ab[l] );
                              kmer.GetFromChars(akmer);
                              kpfi[ koffset[i] + kcount[i] ]
                                   = kmer_pos_fw_id(
                                   kmer, j, b < brc ? True : False, i );    }
                         kcount[i]++;    }    }    }
          cout << Date( ) << ": starting parallel sort, mem usage = "
               << ToStringAddCommas( MemUsageBytes( ) ) << endl;
          sortInPlaceParallel( kpfi.begin( ), kpfi.end( ), NUM_THREADS );
          cout << Date( ) << ": sort complete, mem usage = "
               << ToStringAddCommas( MemUsageBytes( ) ) << endl;
          vec<perf_align> paligns;
          MakeAlignsFromKpfi( K, kpfi, nref, paligns );
          Destroy(kpfi);
          cout << Date( ) << ": sorting paligns, mem usage = "
               << ToStringAddCommas( MemUsageBytes( ) ) << endl;
          sortInPlaceParallel( paligns.begin( ), paligns.end( ), NUM_THREADS );
          CondenseAligns(paligns);
          cout << Date( ) << ": done, found " << paligns.size( ) << " alignments" 
               << endl;

          // Extend alignments.

          for ( size_t i = 0; i < paligns.size( ); i++ )
          {    perf_align& a = paligns[i];
               const basevector& scaff = scaffolds[ a.Id1( ) ];
               const basevector& ref = all[ a.Id2( ) ];
               if ( a.Fw( ) )
               {    while( a.pos1( ) > 0 && a.pos2( ) > 0 )
                    {    if ( scaff[ a.pos1( ) - 1 ] != ref[ a.pos2( ) - 1 ] ) break;
                         if ( scaffolds_gaps[ a.Id1( ) ][ a.pos1( ) - 1 ] ) break;
                         if ( all_amb[ a.Id2( ) ][ a.pos2( ) - 1 ] ) break;
                         a.pos1( ) = a.pos1( ) - 1;
                         a.pos2( ) = a.pos2( ) - 1;    }
                    while( a.Pos1( ) < scaff.isize( ) && a.Pos2( ) < ref.isize( ) )
                    {    if ( scaff[ a.Pos1( ) ] != ref[ a.Pos2( ) ] ) break;
                         if ( scaffolds_gaps[ a.Id1( ) ][ a.Pos1( ) ] ) break;
                         if ( all_amb[ a.Id2( ) ][ a.Pos2( ) ] ) break;
                         a.Pos1( ) = a.Pos1( ) + 1;
                         a.Pos2( ) = a.Pos2( ) + 1;    }    }
               else
               {    while( a.pos1( ) > 0 && a.Pos2( ) < ref.isize( ) )
                    {    if ( scaff[ a.pos1( ) - 1 ] != 3 - ref[ a.Pos2( ) ] ) break;
                         if ( scaffolds_gaps[ a.Id1( ) ][ a.pos1( ) - 1 ] ) break;
                         if ( all_amb[ a.Id2( ) ][ a.Pos2( ) ] ) break;
                         a.pos1( ) = a.pos1( ) - 1;
                         a.Pos2( ) = a.Pos2( ) + 1;    }
                    while( a.Pos1( ) < scaff.isize( ) && a.pos2( ) > 0 )
                    {    if ( scaff[ a.Pos1( ) ] != 3 - ref[ a.pos2( ) - 1 ] ) break;
                         if ( scaffolds_gaps[ a.Id1( ) ][ a.Pos1( ) ] ) break;
                         if ( all_amb[ a.Id2( ) ][ a.pos2( ) - 1 ] ) break;
                         a.Pos1( ) = a.Pos1( ) + 1;
                         a.pos2( ) = a.pos2( ) - 1;    }    }    }

          // Merge.

          cout << Date( ) << ": merging alignments" << endl;
          vec<Bool> to_remove2( paligns.size( ), False );
          for ( size_t i = 0; i < paligns.size( ); i++ )
          {    size_t j;
               for ( j = i + 1; j < paligns.size( ); j++ )
               {    if ( !Comparable( paligns[j], paligns[i] ) ) break;
                    if ( paligns[j].pos1( ) > paligns[j-1].Pos1( ) ) break;
                    if ( paligns[j].Offset( ) != paligns[i].Offset( ) ) break;    }
               paligns[i].Pos1( ) = paligns[j-1].Pos1( );
               paligns[i].Pos2( ) = paligns[j-1].Pos2( );
               for ( size_t k = i + 1; k < j; k++ )
                    to_remove2[k] = True;
               i = j - 1;    }
          EraseIf( paligns, to_remove2 );

          // Write alignments.

          Ofstream( out, ALIGNS );
          for ( size_t i = 0; i < paligns.size( ); i++ )
          {    const perf_align& a = paligns[i];
               out << a.Id1( ) << " " << a.pos1( ) << " " << a.Pos1( ) << " + "
                    << a.Id2( ) << " " << a.pos2( ) << " " << a.Pos2( ) 
                    << " " << ( a.Rc( ) ? "-" : "+" ) << "\n";    }
          cout << Date( ) << ": alignments written" << endl;    }

     // Load alignments and index them.

     cout << Date( ) << ": loading alignments" << endl;
     vec<perf_align> aligns;
     vec< vec<int> > aligns_index(nscaff), aligns_index_pos1(nscaff);
     fast_ifstream in(ALIGNS);
     String line, orient1, orient2;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          istringstream iline( line.c_str( ) );
          int id1, id2, pos1, Pos1, pos2, Pos2;
          iline >> id1 >> pos1 >> Pos1 >> orient1 >> id2 >> pos2 >> Pos2 >> orient2;
          aligns_index[id1].push_back( aligns.size( ) );
          aligns_index_pos1[id1].push_back(pos1);
          aligns.push( id1, id2, pos1, Pos1, pos2, Pos2,
               ( orient2 == "+" ? False : True ) );    }
     for ( int j = 0; j < nscaff; j++ )
          SortSync( aligns_index_pos1[j], aligns_index[j] );

     // Convert to sogs.

     cout << Date( ) << ": converting to sogs" << endl;
     vec<sog> sogs( scaffx.size( ) );
     #pragma omp parallel for
     for ( int si = 0; si < scaffx.isize( ); si++ )
     {    int s = scaffx[si];
          vec<perf_align> sal( aligns_index[s].size( ) );
          for ( int j = 0; j < aligns_index[s].isize( ); j++ )
               sal[j] = aligns[ aligns_index[s][j] ];
          sogs[si] = sog( scaffolds[s], scaffolds_gaps[s], sal );    }
     #pragma omp parallel for
     for ( int si = 0; si < scaffx.isize( ); si++ )
          sogs[si].Validate( );

     // Remove tiny contigs.  Should be temporary.  Note that this creates
     // scaffolds that are all Ns.

     const int min_contig = 1000;
     #pragma omp parallel for
     for ( int si = 0; si < scaffx.isize( ); si++ )
     {    int s = scaffx[si];
          sog& S = sogs[si];
          vec<ho_interval> tigs = S.Tigs( ), tiny;
          for ( int j = 0; j < tigs.isize( ); j++ )
          {    if ( tigs[j].Length( ) < min_contig )
               {    for ( int l = tigs[j].Start( ); l < tigs[j].Stop( ); l++ )
                         S.SetBase( l, 'N' );
                    tiny.push_back( tigs[j] );    }    }
          vec<Bool> to_delete( S.AlignsCount( ), False );
          for ( int j = 0; j < S.AlignsCount( ); j++ )
          {    const perf_align& a = S.Align(j);
               for ( int t = 0; t < tiny.isize( ); t++ )
                    if ( Subset( a.Extent1( ), tiny[t] ) ) to_delete[j] = True;    }
          EraseIf( S.Aligns( ), to_delete );    }

     // Reverse scaffolds as needed to make them mostly forward.

     #pragma omp parallel for
     for ( int si = 0; si < scaffx.isize( ); si++ )
     {    int s = scaffx[si];
          sog& S = sogs[si];
          int fw = 0, rc = 0;
          for ( int j = 0; j < S.AlignsCount( ); j++ )
          {    if ( S.Align(j).Fw( ) ) fw += S.Align(j).Len( );
               else rc += S.Align(j).Len( );    }
          if ( rc > fw ) S.Reverse( );    }
     if (PRINT_INTERMEDIATES) PrintBlocks( scaffx, sogs );

     // Look for local anchors.

     const int max_local_gap = 20000;
     cout << Date( ) << ": looking for local anchors" << endl;
     #pragma omp parallel for
     for ( int si = 0; si < scaffx.isize( ); si++ )
     {    int s = scaffx[si];
          sog& S = sogs[si];
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<perf_align> palignsx;
               for ( int j = 0; j < S.AlignsCount( ) - 1; j++ )
               {    const perf_align& a1 = S.Align(j);
                    if ( !a1.Fw( ) ) continue;
                    int id1 = a1.Id1( ), id2 = a1.Id2( );
                    perf_align& a2 = S.Align(j+1);
                    if ( !Comparable( a1, a2 ) ) continue;
                    if ( !( a1.Pos1( ) < a2.pos1( ) ) ) continue;
                    if ( !( a1.Pos2( ) < a2.pos2( ) ) ) continue;
                    if ( a2.pos1( ) - a1.Pos1( ) > max_local_gap ) continue;
                    if ( a2.pos2( ) - a1.Pos2( ) > max_local_gap ) continue;
                    vec<kmer_pos_fw_id> kpfi;
                    for ( int xpass = 1; xpass <= 2; xpass++ )
                    {    int start = ( xpass == 1 ? a1.Pos1( ) : a1.Pos2( ) );
                         int stop = ( xpass == 1 ? a2.pos1( ) : a2.pos2( ) );
                         vec<char> akmer(100);
                         FASTAVECTOR kmer;
                         int nextamb;
                         for ( nextamb = start; nextamb < stop; nextamb++ )
                         {    if ( xpass == 1 )
                              {    if ( S.Gap(nextamb) ) break;    }
                              else
                              {    if ( all_amb[id2][nextamb] ) break;    }    }
                         for ( int j = start; j <= stop - K; j++ )
                         {    if ( j + K > nextamb )
                              {    for ( ; nextamb < stop; nextamb++ )
                                   {    if ( xpass == 1 )
                                        {    if ( !S.Gap(nextamb) ) break;    }
                                        else
                                        {    if ( !all_amb[id2][nextamb] ) 
                                                  break;   }   }
                                   j = nextamb;
                                   if ( j == stop ) break;
                                   for ( ; nextamb < stop; nextamb++ )
                                   {    if ( xpass == 1 )
                                        {    if ( S.Gap(nextamb) ) break;    }
                                        else
                                        {    if ( all_amb[id2][nextamb] ) 
                                                  break;    }    }
                                   j--;
                                   continue;    }
                              basevector b(K), brc;
                              if ( xpass == 1 )
                              {    for ( int l = 0; l < K; l++ )
                                        b.Set( l, as_char( S.Base(j+l) ) );    }
                              else 
                                   b.SetToSubOf( all[id2], j, K );
                              brc = b;
                              brc.ReverseComplement( );
                              basevector& ab = ( b < brc ? b : brc );
                              for ( int l = 0; l < K; l++ )
                                   akmer[l] = as_base( ab[l] );
                              kmer.GetFromChars(akmer);
                              kpfi.push( kmer, j, b < brc ? True : False, 
                                   ( xpass == 1 ? id1 : id2 ) );    }    }
                    Sort(kpfi);
                    vec<perf_align> paligns2;
                    MakeAlignsFromKpfi( K, kpfi, nref, paligns2 );
                    Sort(paligns2);
                    CondenseAligns(paligns2);
                    palignsx.append(paligns2);    }
               S.Aligns( ).append(palignsx);
               Sort( S.Aligns( ) );
               S.Reverse( );    }    }
     #pragma omp parallel for
     for ( int si = 0; si < scaffx.isize( ); si++ )
          sogs[si].Validate( );

     // Identify and remove duff.  Note that we should probably not do this if
     // the duff is at the very end of a scaffold.

     int max_duff = 200;
     #pragma omp parallel for
     for ( int si = 0; si < scaffx.isize( ); si++ )
     {    int s = scaffx[si];
          sog& S = sogs[si];
          vec< vec<ho_interval> > scov(nref), scov2(nref);
          vec<Bool> to_delete( S.AlignsCount( ), False );
          for ( int j = 0; j < S.AlignsCount( ); j++ )
          {    const perf_align& a = S.Align(j);
               scov[ a.Id2( ) ].push_back( a.Extent2( ) );    }
          for ( int i = 0; i < nref; i++ )
               ExtractGivenCoverage( all[i].size( ), 1, scov[i], scov2[i] );
          for ( int j = 0; j < S.AlignsCount( ); j++ )
          {    const perf_align& a = S.Align(j);
               int id2 = a.Id2( );
               if ( a.Len( ) > max_duff ) continue;
               Bool solo = True;
               for ( int r = 0; r < scov2[id2].isize( ); r++ )
               {    if ( scov2[id2][r] == a.Extent2( ) ) continue;
                    if ( Distance( scov2[id2][r], a.Extent2( ) ) < 50000 )
                         solo = False;    }
               if (solo) to_delete[j] = True;    }
          EraseIf( S.Aligns( ), to_delete );    }

     // Define coverage of genome.  I've tried to speed this up, without success.

     cout << Date( ) << ": computing coverage" << endl;
     vecbitvector gcov;
     Mimic( all, gcov ); // NO NO NO: all is bigger than cov.
     for ( int s = 0; s < nscaff; s++ )
     {    for ( int j = 0; j < aligns_index[s].isize( ); j++ )
          {    const perf_align& a = aligns[ aligns_index[s][j] ];
               for ( int l = a.pos2( ); l < a.Pos2( ); l++ )
                    gcov[ a.Id2( ) ].Set( l, True );    }    }

     // Strip ends: if the distance between a gap and the nearest aligned thing
     // is less than 100 bases, extend the gap to cover it.  We record these as
     // hanging ends.

     cout << Date( ) << ": stripping ends" << endl;
     #pragma omp parallel for
     for ( int si = 0; si < scaffx.isize( ); si++ )
     {    int s = scaffx[si];
          sog& S = sogs[si];
          bitvector cov( S.BaseCount( ), False );
          bitvector extra_gap( S.BaseCount( ), False );
          for ( int j = 0; j < S.AlignsCount( ); j++ )
          {    const perf_align& a = S.Align(j);
               for ( int k = a.pos1( ); k < a.Pos1( ); k++ )
                    cov.Set( k, True );    }
          int last_cov = -1;
          for ( int j = 0; j < S.BaseCount( ); j++ )
          {    if ( cov[j] ) last_cov = j;
               else if ( S.Gap(j) )
               {    if ( j - last_cov - 1 < 100 )
                    {    if ( j - last_cov - 1 > 0 )
                              RecordEvent( "hanging", j - last_cov - 1 );
                         for ( int k = last_cov + 1; k < j; k++ )
                              extra_gap.Set( k, True );    }
                    while( j < S.BaseCount( ) && S.Gap(j) ) j++;
                    for ( int k = j; k < Min( j + 100 - 1, S.BaseCount( ) ); k++ )
                    {    if ( cov[k] )
                         {    if ( k - j > 1 ) RecordEvent( "hanging", k-j );
                              for ( int r = j; r < k; r++ )
                                   extra_gap.Set( r, True );
                              break;    }    }
                    j--;    }    }
          for ( int j = 0; j < S.BaseCount( ); j++ )
               if ( extra_gap[j] ) S.SetBase( j, 'N' );    }
     if (PRINT_INTERMEDIATES) PrintBlocks( scaffx, sogs );

     // Identify contig transpositions.  Implementation is quite restrictive.

     cout << Date( ) << ": identifying contig transpositions" << endl;
     #pragma omp parallel for
     for ( int si = 0; si < scaffx.isize( ); si++ )
     {    int s = scaffx[si];
          sog& S = sogs[si];
          vec<ho_interval> tigs = S.Tigs( );
          int ntigs = tigs.size( );
          if ( ntigs == 0 ) continue;
          vec< vec<int> > places(ntigs);
          for ( int j = 0; j < S.AlignsCount( ); j++ )
          {    const perf_align& a = S.Align(j);
               for ( int t = 0; t < ntigs; t++ )
               {    if ( Subset( a.Extent1( ), tigs[t] ) )
                    {    places[t].push_back(j);
                         break;    }    }    }
          // for ( int t = 0; t < ntigs; t++ )
          // {    cout << "\nplacements of tig " << t << endl;
          //      for ( int j = 0; j < places[t].isize( ); j++ )
          //           cout << S.Align( places[t][j] ) << endl;    }
          vec<Bool> swap_next( ntigs - 1, False );
          const int max_dist = 5000;
          for ( int t = 0; t < ntigs - 1; t++ )
          {    // Note restrictive implementation defined by the next line:
               if ( !places[t].solo( ) || !places[t+1].solo( ) ) continue;
               const perf_align& a1 = S.Align( places[t][0] );
               const perf_align& a2 = S.Align( places[t+1][0] );
               if ( !Comparable( a1, a2 ) ) continue;
               if ( ( a1.Fw( ) && a1.pos2( ) >= a2.Pos2( )
                    && a1.pos2( ) - a2.Pos2( ) <= max_dist )
                    || ( a1.Rc( ) && a2.pos2( ) >= a1.Pos2( )
                    && a2.pos2( ) - a1.Pos2( ) <= max_dist ) )
               {    swap_next[t] = True;    }    }
          for ( int t = 0; t < ntigs - 1; t++ )
          {    if ( !swap_next[t] ) continue;
               if ( t > 0 && swap_next[t-1] ) continue;
               if ( t < ntigs - 2 && swap_next[t+1] ) continue;
               int n = Max( tigs[t].Length( ), tigs[t+1].Length( ) );
               RecordEvent( "transposition", n );
               vec<char> newseq;
               for ( int j = tigs[t+1].Start( ); j < tigs[t+1].Stop( ); j++ )
                    newseq.push_back( S.Base(j) );
               for ( int j = tigs[t].Stop( ); j < tigs[t+1].Start( ); j++ )
                    newseq.push_back( 'N' );
               for ( int j = tigs[t].Start( ); j < tigs[t].Stop( ); j++ )
                    newseq.push_back( S.Base(j) );
               for ( int j = 0; j < newseq.isize( ); j++ )
                    S.SetBase( tigs[t].Start( ) + j, newseq[j] );
               perf_align& a1 = S.Align( places[t][0] );
               perf_align& a2 = S.Align( places[t+1][0] );
               if (SHOW_EVENTS) 
               {    cout << "\nswapping contigs " << s << "." << t << " and "
                         << s << "." << t+1 << "\n";
                    PRINT(a1); PRINT(a2);    }
               int len1 = a1.Len( ), len2 = a2.Len( );
               int a1pos1 = a1.pos1( ), a1pos2 = a1.pos2( );
               a1.pos1( ) = tigs[t].Start( ) + a2.pos1( ) - tigs[t+1].Start( );
               a1.Pos1( ) = a1.pos1( ) + len2;
               a1.pos2( ) = a2.pos2( );
               a1.Pos2( ) = a2.Pos2( );
               a2.pos1( ) = tigs[t+1].Stop( ) - tigs[t].Stop( ) + a1pos1;
               a2.Pos1( ) = a2.pos1( ) + len1;
               a2.pos2( ) = a1pos2;
               a2.Pos2( ) = a2.pos2( ) + len1;    }    }
     if (PRINT_INTERMEDIATES) PrintBlocks( scaffx, sogs );
     for ( int si = 0; si < scaffx.isize( ); si++ )
          sogs[si].Validate( );

     // Do two passes, one after reversing.

     vec<String> insertions(nscaff);
     for ( int pass = 1; pass <= 2; pass++ )
     {    
          // Close clean gaps.

          cout << Date( ) << ": closing clean gaps, pass " << pass << endl;
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
          {    int s = scaffx[si];
               sog& S = sogs[si];
               for ( int j = 0; j < S.AlignsCount( ) - 1; j++ )
               {    const perf_align& a1 = S.Align(j);
                    if ( !a1.Fw( ) ) continue;
                    int id1 = a1.Id1( ), id2 = a1.Id2( );
                    perf_align& a2 = S.Align(j+1);
                    if ( !Comparable( a1, a2 ) ) continue;
                    if ( !( a1.Pos1( ) < a2.pos1( ) ) ) continue;
                    if ( !( a1.Pos2( ) < a2.pos2( ) ) ) continue;
                    Bool is_gap = True;
                    for ( int l = a1.Pos1( ); l < a2.pos1( ); l++ )
                    {    if ( !S.Gap(l) )
                         {    is_gap = False;
                              break;    }    }
                    if ( !is_gap ) continue;
                    vec<char> replacement;
                    for ( int l = a1.Pos2( ); l < a2.pos2( ); l++ )
                         replacement.push_back( S.GenomeBase( id2, l ) );
                    int start = a1.Pos1( ), stop = a2.pos1( );
                    perf_align a( id1, id2, a1.pos1( ), 
                         a2.Pos1( ) + replacement.isize( ) - ( stop - start ), 
                         a1.pos2( ), a2.Pos2( ), a1.Rc( ) );
                    S.Replace( start, stop, replacement, j, a );
                    j--;    }    }
          if (PRINT_INTERMEDIATES) PrintBlocks( scaffx, sogs );
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
               sogs[si].Validate( );

          // Handle cases of overlapping adjacent contigs.

          cout << Date( ) << ": overlapping adjacent contigs, pass " << pass << endl;
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
               MergeOverlapping( sogs[si] );
          if (PRINT_INTERMEDIATES) PrintBlocks( scaffx, sogs );
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
               sogs[si].Validate( );

          // Close indels.

          cout << Date( ) << ": closing indels, pass " << pass << endl;
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
          {    int s = scaffx[si];
               sog& S = sogs[si];
               for ( int j = 0; j < S.AlignsCount( ) - 1; j++ )
               {    const perf_align& a1 = S.Align(j);
                    if ( !a1.Fw( ) ) continue;
                    int id1 = a1.Id1( ), id2 = a1.Id2( );
                    perf_align& a2 = S.Align(j+1);
                    if ( !Comparable( a1, a2 ) ) continue;
                    int a1Pos1 = a1.Pos1( ), a1Pos2 = a1.Pos2( );
                    int left_overlap = a1Pos1 - a2.pos1( );
                    int right_overlap = a1Pos2 - a2.pos2( );
                    const int max_overlap = 100;
                    if ( left_overlap > max_overlap ) continue;
                    if ( right_overlap > max_overlap ) continue;
                    if ( left_overlap < 0 && right_overlap < 0 ) continue;
                    if ( left_overlap > 0 || right_overlap > 0 )
                    {    int M = Max( left_overlap, right_overlap );
                         left_overlap -= M, right_overlap -= M;
                         a1Pos1 -= M, a1Pos2 -= M;
                         if ( a1Pos1 < a1.pos1( ) ) continue;
                         if ( a1Pos2 < a1.pos2( ) ) continue;    }
                    if ( !( a1Pos1 <= a2.pos1( ) ) ) continue;
                    // The following wasn't supposed to happen.
                    if ( !( left_overlap < 0 ^ right_overlap < 0 ) ) continue;
                    perf_align a( id1, id2, a1.pos1( ), 
                         a2.Pos1( ) + left_overlap - right_overlap,
                         a1.pos2( ), a2.Pos2( ), a1.Rc( ) );
                    if ( !S.ReplaceOK( a1Pos1, a2.pos1( ), j ) ) continue;
                   
                    // For deletions, make sure we're not deleting something
                    // that's already covered.  We may not have the right
                    // statistic for this.

                    if ( right_overlap < 0 )
                    {    int cov = 0, start = a1Pos2, stop = a2.pos2( );
                         for ( int l = start; l < stop; l++ )
                              if ( gcov[id2][l] ) cov++;
                         int len = stop - start;
                         if ( cov > len/2 && len-cov >= 10 ) continue;    }

                    // Make the join.

                    vec<char> replacement;
                    if ( left_overlap < 0 )
                    {    int start, stop;
                         for ( start = a1Pos1; start < a2.pos1( ); start++ )
                              if ( !S.Gap(start) ) break;
                         for ( stop = a2.pos1( ) - 1; stop >= a1Pos1; stop-- )
                              if ( !S.Gap(stop) ) break;
                         if (SHOW_EVENTS)
                         {    cout << "insertion type A of\n"; 
                              S.PrintFasta( cout, start, stop );
                              cout << "at " << id2 << "." << a1Pos2 << "\n";    }
                         int n = stop - start;
                         if ( n >= dump_insertions_min )
                         {    insertions[s] += ">" + ToString(id2) + "." 
                                   + ToString(a1Pos2) + " (" + ToString(n) + ")\n";
                              int count = 0;
                              for ( int l = start; l < stop; l++ )
                              {    if ( count++ == 80 ) 
                                   {    insertions[s] += "\n";
                                        count = 0;    }
                                   insertions[s] += S.Base(l);    }
                              insertions[s] += "\n";    }
                         RecordEvent( "insertion", n );    }
                    if ( right_overlap < 0 )
                    {    if (SHOW_EVENTS)
                         {    cout << "deletion of " << id2 << "."
                                   << a1Pos2 << "-" << a2.pos2( ) 
                                   << " (" << a2.pos2( ) - a1Pos2 
                                   << " bases)\n";    }
                         for ( int l = a1Pos2; l < a2.pos2( ); l++ )
                              replacement.push_back( S.GenomeBase( id2,l ) );
                         int n = a2.pos2( ) - a1Pos2;
                         RecordEvent( "deletion", n );    }
                    S.Replace( a1Pos1, a2.pos1( ), replacement, j, a );
                    j--;    }    }
          if (PRINT_INTERMEDIATES) PrintBlocks( scaffx, sogs );
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
               sogs[si].Validate( );

          // Identify single mismatches and combine alignments across them.

          cout << Date( ) << ": finding mismatches, pass " << pass << endl;
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
          {    int s = scaffx[si];
               sog& S = sogs[si];
               vec<Bool> to_delete( S.AlignsCount( ), False );
               const int spread = 2;
               for ( int j = 0; j < S.AlignsCount( ) - spread; j++ )
               {    if ( to_delete[j] ) continue;
                    for ( int k = 1; k <= spread; k++ )
                    {    if ( to_delete[j+k] ) continue;
                         const perf_align& a1 = S.Align(j);
                         if ( !a1.Fw( ) ) continue;
                         int id1 = a1.Id1( ), id2 = a1.Id2( );
                         perf_align& a2 = S.Align(j+k);
                         if ( !Comparable( a1, a2 ) ) continue;
                         if ( a2.pos1( ) - a1.Pos1( ) != 1 ) continue;
                         if ( a2.pos2( ) - a1.Pos2( ) != 1 ) continue;

                         // Check for other small alignments that would be subsumed
                         // by the join.
     
                         int left = a1.pos1( ), right = a2.Pos1( );
                         vec<int> junk;
                         Bool fail = False;
                         for ( int z = 0; z < S.AlignsCount( ); z++ )
                         {    if ( to_delete[z] ) continue;
                              if ( z == j || z == j + k ) continue;
                              perf_align& a = S.Align(z);
                              if ( a.pos1( ) >= left && a.Pos1( ) <= right )
                              {    if ( a.Len( ) > a1.Len( ) / 5 ||
                                        a.Len( ) > a2.Len( ) / 5 )
                                   {    fail = True;
                                        break;    }
                                   junk.push_back(z);    }    
                              else if ( a.pos1( ) <= a1.Pos1( ) 
                                   && a1.Pos1( ) <= a.Pos1( ) )
                              {    fail = True;
                                   break;    }    }
                         if (fail) continue;
                         for ( int l = 0; l < junk.isize( ); l++ )
                              to_delete[ junk[l] ] = True;
     
                         // Make the join.
     
                         char ref_base = S.GenomeBase( id2, a1.Pos2( ) );
                         char ass_base = S.Base( a1.Pos1( ) );
                         if (SHOW_EVENTS)
                         {    cout << "mismatch " << ref_base << " --> " << ass_base 
                                   << " at " << a1.Id2( ) << "." << a1.Pos2( ) 
                                   << "\n";    }
                         RecordEvent( "substitution", 1 );
                         S.SetBase( a1.Pos1( ), ref_base );
                         a2.pos1( ) = a1.pos1( ), a2.pos2( ) = a1.pos2( );
                         to_delete[j] = True;    }    }
               EraseIf( S.Aligns( ), to_delete );    }
          if (PRINT_INTERMEDIATES) PrintBlocks( scaffx, sogs );
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
               sogs[si].Validate( );

          // Close gaps between alignments.

          cout << Date( ) << ": closing alignment gaps, pass " << pass << endl;
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
          {    int s = scaffx[si];
               sog& S = sogs[si];
               for ( int j = 0; j < S.AlignsCount( ) - 1; j++ )
               {    const perf_align& a1 = S.Align(j);
                    if ( !a1.Fw( ) ) continue;
                    int id1 = a1.Id1( ), id2 = a1.Id2( );
                    perf_align& a2 = S.Align(j+1);
                    if ( !Comparable( a1, a2 ) ) continue;
                    int a1Pos1 = a1.Pos1( ), a1Pos2 = a1.Pos2( );
                    int left_overlap = a1Pos1 - a2.pos1( );
                    int right_overlap = a1Pos2 - a2.pos2( );
                    const int max_gap = 100000;
                    if ( left_overlap >= 0 || right_overlap >= 0 ) continue;
                    if ( -left_overlap > max_gap || -right_overlap > max_gap ) 
                         continue;
                    perf_align a( id1, id2, a1.pos1( ), 
                         a2.Pos1( ) + left_overlap - right_overlap, a1.pos2( ), 
                         a2.Pos2( ), a1.Rc( ) );
                    if ( !S.ReplaceOK( a1Pos1, a2.pos1( ), j ) ) continue;
                    vec<char> replacement, orig;
                    int n = Max( a2.pos1( ) - a1Pos1, a2.pos2( ) - a1Pos2 );
                    for ( int l = a1Pos2; l < a2.pos2( ); l++ )
                         replacement.push_back( S.GenomeBase( id2, l ) );
                    if (SHOW_EVENTS)
                    {    cout << "\nreplacing " << id2 << "." << a1Pos2 << "-" 
                              << a2.pos2( ) << " (" << a2.pos2( ) - a1Pos2 
                              << " bases)\n";
                         int count = 0;
                         for ( int l = a1Pos2; l < a2.pos2( ); l++ )
                         {    if ( count++ == 80 ) 
                              {    cout << "\n";
                                   count = 1;    }
                              cout << S.GenomeBase( id2, l );    }
                         cout << "\nby\n";
                         S.PrintFasta( cout, a1Pos1, a2.pos1( ) );   
                         cout << "\n";    }
                    for ( int l = a1Pos1; l < a2.pos1( ); l++ )
                         orig.push_back( S.Base(l) );
                    Bool leftNfree = True, rightNfree = True;
                    Bool left_all_N = True, right_all_N = True;
                    for ( int l = 0; l < orig.isize( ); l++ )
                    {    if ( orig[l] == 'N' ) leftNfree = False;
                         else left_all_N = False;    }
                    for ( int l = 0; l < replacement.isize( ); l++ )
                    {    if ( replacement[l] == 'N' ) rightNfree = False;
                         else right_all_N = False;    }
                    Bool Nfree = ( leftNfree && rightNfree );
                    const int flank = 50;
                    Bool decomposable = False;
                    vec<char> insert;
                    int start, stop;
                    for ( start = 0; start < orig.isize( ); start++ )
                         if ( orig[start] != 'N' ) break;
                    for ( stop = orig.isize( ) - 1; stop >= 0; stop-- )
                         if ( orig[stop] != 'N' ) break;
                    for ( int l = start; l <= stop; l++ )
                         insert.push_back( orig[l] );
                    if ( !left_all_N && right_all_N )
                    {    decomposable = True;
                         int n = insert.size( );
                         if (SHOW_EVENTS) 
                              cout << "insertion type B of size " << n << "\n";
                         if ( n >= dump_insertions_min )
                         {    String header = ToString(id2) + "." 
                                   + ToString(a1Pos2) + "-" + ToString( a2.pos2( ) )
                                   + " (" + ToString(n) + ")";
                              insertions[s] += ">" + header + "\n";
                              int count = 0;
                              for ( int z = 0; z < insert.isize( ); z++ )
                              {    if ( count++ == 80 )
                                   {    insertions[s].push_back( '\n' );
                                        count = 1;    }
                                   insertions[s].push_back( insert[z] );    }
                              insertions[s].push_back( '\n' );    }
                         RecordEvent( "insertion", n );    }
                    else if ( !leftNfree && rightNfree )
                    {    decomposable = True;
                         for ( int l = 0; l < orig.isize( ); l++ )
                         {    if ( orig[l] == 'N' ) continue;
                              int m;
                              for ( m = l + 1; m < orig.isize( ); m++ )
                                   if ( orig[m] == 'N' ) break;
                              int n = m - l;
                              if (SHOW_EVENTS)
                                   cout << "insertion type C of size " << n << "\n";
                              if ( n >= dump_insertions_min )
                              {    String header = ToString(id2) + "." 
                                        + ToString(a1Pos2) + "-" 
                                        + ToString( a2.pos2( ) )
                                        + " (" + ToString(n) + ")";
                                   insertions[s] += ">" + header + "\n";
                                   int count = 0;
                                   for ( int z = l; z < m; z++ )
                                   {    if ( count++ == 80 )
                                        {    insertions[s].push_back( '\n' );
                                             count = 1;    }
                                        insertions[s].push_back( orig[z] );    }
                                   insertions[s].push_back( '\n' );    }
                              RecordEvent( "insertion", n );
                              l = m;    }    }
                    else if ( !leftNfree && !rightNfree )
                    {    decomposable = True;
                         int n = insert.size( );
                         if (SHOW_EVENTS) 
                              cout << "insertion type D of size " << n << "\n";
                         if ( n >= dump_insertions_min )
                         {    String header = ToString(id2) + "." 
                                   + ToString(a1Pos2) + "-" + ToString( a2.pos2( ) )
                                   + " (" + ToString(n) + ")";
                              insertions[s] += ">" + header + "\n";
                              int count = 0;
                              for ( int z = 0; z < insert.isize( ); z++ )
                              {    if ( count++ == 80 )
                                   {    insertions[s].push_back( '\n' );
                                        count = 1;    }
                                   insertions[s].push_back( insert[z] );    }
                              insertions[s].push_back( '\n' );    }
                         RecordEvent( "insertion", n );
                         int gstart, gstop;
                         for ( gstart = a1Pos2; gstart < a2.pos2( ); gstart++ )
                              if ( S.GenomeBase( id2, gstart ) != 'N' ) break;
                         for ( gstop = a2.pos2( ) - 1; gstop >= a1Pos2; gstop-- )
                              if ( S.GenomeBase( id2, gstop ) != 'N' ) break;
                         int nd = gstop - gstart + 1;
                         if (SHOW_EVENTS) cout << "deletion of size " << nd << "\n";
                         RecordEvent( "deletion", nd );    }
                    else if ( Nfree && a1Pos1 >= flank && a1Pos2 >= flank 
                         && S.BaseCount( ) - a2.pos1( ) >= flank
                         && all[id2].isize( ) - a2.pos2( ) >= flank )
                    {    decomposable = True;
                         basevector U( orig.size( ) + 2*flank ); 
                         basevector T( replacement.size( ) + 2*flank );
                         for ( int l = 0; l < U.isize( ); l++ )
                              U.Set( l, as_char( S.Base( a1Pos1 + l - flank ) ) );
                         for ( int l = 0; l < T.isize( ); l++ )
                              T.Set( l, all[id2][ a1Pos2 + l - flank ] );
                         alignment a;
                         SmithWatAffine( U, T, a );
                         if ( a.pos1( ) > 0 || a.Pos1( ) < U.isize( ) )
                              decomposable = False;
                         if (decomposable)
                         {    if (SHOW_EVENTS) 
                                   PrintVisualAlignment( True, cout, U, T, a );    
                              avector<int> gaps, lengths;
                              int pos1, pos2, errors;
                              a.Unpack( pos1, pos2, errors, gaps, lengths );
                              int p1 = pos1, p2 = pos2;
                              int mismatches = 0;
                              for ( unsigned int j = 0; j < gaps.length; j++ )
                              {    if ( gaps(j) > 0 )
                                   {    if (SHOW_EVENTS)
                                        {    cout << "deletion of size " << gaps(j) 
                                                  << "\n";    }
                                        int n = gaps(j);
                                        RecordEvent( "deletion", n );
                                        p2 += gaps(j);    }
                                   if ( gaps(j) < 0 )
                                   {    int n = -gaps(j);
                                        if (SHOW_EVENTS)
                                        {    cout << "insertion type E of size " 
                                                  << n << "\n";    }
                                        if ( n >= dump_insertions_min )
                                        {    insertions[s] += ">" + ToString(id2) 
                                                  + "." + ToString(a1Pos2+p2-flank) 
                                                  + " (" + ToString(n) + ")\n";
                                             int count = 0;
                                             for ( int l = p1; l < p1+n; l++ )
                                             {    if ( count++ == 80 ) 
                                                  {    insertions[s] += "\n";
                                                       count = 0;    }
                                                  insertions[s] += S.Base(l);    }
                                             insertions[s] += "\n";    }
                                        RecordEvent( "insertion", n );
                                        p1 -= gaps(j);    }
                                   for ( int x = 0; x < lengths(j); x++ )
                                   {    if ( U[p1] != T[p2] ) mismatches++;
                                        ++p1;
                                        ++p2;    }    }
                              if (SHOW_EVENTS) cout << mismatches << " mismatches\n";
                              for ( int z = 0; z < mismatches; z++ )
                                   RecordEvent( "substitution", 1 );    }    }
                    if ( !decomposable )
                    {    if (SHOW_EVENTS) cout << "not decomposable\n";
                         RecordEvent( "replacement", n );    }
                    S.Replace( a1Pos1, a2.pos1( ), replacement, j, a );
                    j--;    }    }
          if (PRINT_INTERMEDIATES) PrintBlocks( scaffx, sogs );
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
               sogs[si].Validate( );

          // Again merge overlapping adjacent contigs.

          cout << Date( ) << ": overlapping again, pass " << pass << endl;
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
               MergeOverlapping( sogs[si] );
          if (PRINT_INTERMEDIATES) PrintBlocks( scaffx, sogs );
          #pragma omp parallel for
          for ( int si = 0; si < scaffx.isize( ); si++ )
               sogs[si].Validate( );

          // Now reverse.

          for ( int si = 0; si < scaffx.isize( ); si++ )
               sogs[si].Reverse( );    }

     // Find free scaffolds.  Because we use alignments based on 100-mers that are
     // unique in the reference, it's in theory possible that scaffolds will be
     // declared free, but in fact have a good home somewhere on the genome.  This
     // may happen, but probably at a low rate.  On a test of 7 human scaffolds
     // (containing 37 contigs), there were no such instances.
     //
     // For now we only consider the cases where there is no alignment at all.
     // This isn't right.  Moreover, if we want to declare something free, even
     // though it has some alignments, we should proceed by *deleting* those
     // alignments.  This is consistent with how the sog::Print function labels
     // things.

     ofstream* doutp = NULL;
     if ( DUMP_FREE != "" )
     {    doutp = new ofstream;
          OpenOfstream( *doutp, DUMP_FREE );   }
     for ( int si = 0; si < scaffx.isize( ); si++ )
     {    int s = scaffx[si];
          const sog& S = sogs[si];
          if ( S.AlignsCount( ) > 0 ) continue;
          int n = S.BaseCount( );
          RecordEvent( "free", n );
          if ( DUMP_FREE != "" )
          {    int contig = 1;
               for ( int j = 0; j < S.BaseCount( ); j++ )
               {    if ( S.Gap(j) ) continue;
                    int l;
                    for ( l = j + 1; l < S.BaseCount( ); l++ )
                         if ( S.Gap(l) ) break;
                    S.PrintFasta( *doutp, j, l,
                         "scaffold_" + ToString(s) + "." + ToString(contig++)
                         + " (" + ToString(j) + "-" + ToString(l) + ")" );
                    j = l;    }    }    }
               
     // Print insertions.

     if ( DUMP_INSERTIONS != "" )
     {    Ofstream( dout, dump_insertions_fn );
          for ( int s = 0; s < nscaff; s++ )
               dout << insertions[s];    }
     
     // Print alignment blocks.

     PrintBlocks( scaffx, sogs );
     if ( SCAFF_OUT != "" )
     {    Ofstream( sout, SCAFF_OUT );
          for ( int si = 0; si < scaffx.isize( ); si++ )
          {    int s = scaffx[si];
               const sog& S = sogs[si];
               for ( int j = 0; j < S.AlignsCount( ); j++ )
               {    const perf_align& a1 = S.Align(j);
                    S.PrintFasta( sout, a1.pos1( ), a1.Pos1( ),
                         ToString(s) + "." + ToString( a1.pos1( ) ) + "-"
                         + ToString( a1.Pos1( ) ) );
                    if ( j < S.AlignsCount( ) - 1 )
                    {    const perf_align& a2 = S.Align(j+1);
                         if ( Comparable( a1, a2 ) && a1.Pos1( ) < a2.pos1( )
                              && a2.pos1( ) - a1.Pos1( ) <= 100000 )
                         {    S.PrintFasta( sout, a1.Pos1( ), a2.pos1( ),
                                   ToString(s) + "." + ToString( a1.Pos1( ) ) + "-"
                                   + ToString( a2.pos1( ) ) );   }   }   }   }   }

     // Print edit stats.

     cout << "SUMMARY OF EDITS\n";
     vec< vec<String> > rows;
     vec<String> row;
     row.push_back( "type", "<= 1", "<= 10", "<= 100", "<= 1000", "<= 10000", 
          "<= 100000", "> 100000" );
     rows.push_back(row); 
     for ( int t = 0; t < event_types.isize( ); t++ )
     {    row.clear( );
          row.push_back( event_types[t] );
          for ( int j = 0; j < event_count[t].isize( ); j++ )
               row.push_back( ToString( event_count[t][j] ) );
          rows.push_back(row);    }
     PrintTabular( cout, rows, 3, "lrrrrrrr" );
     cout << "\nTime used = " << TimeSince(clock) << "." << endl;    }
