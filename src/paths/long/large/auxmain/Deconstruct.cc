///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Deconstruct a scaffolded assembly made from a DISCOVAR de novo assembly
// from a.lines.fasta.  Hardwired for the moment.
//
// Probably has a bunch of bugs.
//
// Need to run on a terabyte box.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "paths/RemodelGapTools.h"

String GetOp( const String& lid, const vec<int>& inv, 
     const vec<vec<vec<vec<int>>>>& lines )
{    int l = lid.substr( 1, lid.isize( ) - 1 ).Int( );
     int lr = 1000000;
     if ( l > 0 && inv[ lines[l-1].front( )[0][0] ] ==
          lines[l].back( )[0][0] )
     {    lr = l-1;    }
     if ( l < lines.isize( ) - 1 && inv[ lines[l+1].front( )[0][0] ] ==
          lines[l].back( )[0][0] )
     {    lr = l+1;    }
     if ( lid[0] == '+' ) return '-' + ToString(lr);
     return '+' + ToString(lr);    }

int main( )
{
     RunTime( );

     // Define hardwired targets.

     String dir = "/wga/scr4/jaffe/GapToy/51400.newchem/a.final";
     String scaff = "/wga/scr4/jaffe/dovetail/discovar_hirise_i3.fasta";

     // Load assembly data.

     vec<int> inv;
     BinaryReader::readFile( dir + "/a.inv", &inv );
     vec<vec<vec<vec<int>>>> linesx;
     BinaryReader::readFile( dir + "/a.lines", &linesx );

     // Load a.lines.fasta.

     cout << Date( ) << ": loading lines" << endl;
     vec<String> lines, lid;
     const int min_line = 1000;
     String line, current;
     current.reserve(300000000);
     {    fast_ifstream in( dir + "/a.lines.fasta" );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( ">", 0 ) )
               {    if ( lid.nonempty( ) ) 
                    {    if ( current.isize( ) >= min_line )
                              lines.push_back(current);
                         else lid.pop_back( );
                         current.clear( );    }
                    lid.push_back( "+" + line.After( "line_" ) );
                    if ( lid.back( ).Contains( " " ) )
                         lid.back( ) = lid.back( ).Before( " " );    }
               else current += line;    }    }
     lines.push_back(current);
     PRINT( lines.size( ) );

     // Add their reverse complements.

     cout << Date( ) << ": adding in reverse complements" << endl;
     int L = lines.size( );
     lines.resize(2*L), lid.resize(2*L);
     #pragma omp parallel for
     for ( int i = 0; i < L; i++ )
     {    String s;
          StringReverseComplement( lines[i], s );
          lines[i+L] = s;
          lid[i+L] = "-" + lid[i].After( "+" );    }

     // Load scaffolds.

     cout << Date( ) << ": loading scaffolds" << endl;
     current.clear( );
     vec<String> scaffolds;
     {    fast_ifstream in(scaff);
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( ">", 0 ) )
               {    if ( current.nonempty( ) ) 
                    {    scaffolds.push_back(current);
                         current.clear( );    }    }
               else current += line;    }    }
     scaffolds.push_back(current);

     // Get first kmer from each line.

     cout << Date( ) << ": getting first kmers from lines" << endl;
     const int K = 200;
     vec< pair< kmer<K>, int > > X( lines.size( ) );
     for ( int i = 0; i < lines.isize( ); i++ )
     {    X[i].second = i;
          ForceAssertGe( lines[i].isize( ), K );
          String s = lines[i].substr( 0, K );
          basevector b(s);
          X[i].first = b;    }
     cout << Date( ) << ": sorting first kmers" << endl;
     ParallelSort(X);

     // Get all kmers from lines.

     cout << Date( ) << ": getting all kmers from lines" << endl;
     vecbasevector xlines( lines.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < lines.isize( ); i++ )
     {    String s = lines[i];
          for ( int j = 0; j < s.isize( ); j++ )
               if ( s[j] == 'N' ) s[j] = 'A';
          xlines[i] = basevector(s);    }
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup0( xlines, kmers_plus );

     // Match lines to scaffolds.

     cout << Date( ) << ": deconstructing scaffolds" << endl;
     for ( int s = 0; s < scaffolds.isize( ); s++ )
     {    const String& S = scaffolds[s];
          int si = s;
          cout << "\nscaffold " << s << "[l=" << S.size( ) << "]" << endl;
          int pos = 0;
          top: while( pos < S.isize( ) )
          {   
               int n = S.isize( ) - pos;
               String s = S.substr( pos, Min( n, K ) );

               // Check for an N in s.

               for ( int i = 0; i < s.isize( ); i++ )
               {    if ( s[i] == 'N' )
                    {    cout << "Some bases: ";
                         for ( int j = 0; j < i; j++ )
                              cout << s[j];
                         cout << "\n";
                         int j;
                         for ( j = i + 1; j < s.isize( ); j++ )
                              if ( s[j] != 'N' ) break;
                         cout << "see gap of size " << j - i << endl;
                         pos += j;
                         goto top;    }    }

               // Check for too short.

               if ( n < K )
               {    cout << "Some bases: " << s << "\n";
                    pos += n;
                    goto top;    }

               basevector b(s);
               kmer<K> x(b);
               int64_t low = LowerBound1( X, x ), high = UpperBound1( X, x );
               int id = -1, lpos = -1;

               if ( high - low >= 1 )
               {    vec<int> stops( high - low, 0 );
                    for ( int64_t u = low; u < high; u++ )
                    {    int id = X[u].second, mpos = pos;
                         const String& t = lines[id];
                         int j;
                         for ( j = 0; j < t.isize( ); j++ )
                         {    if ( t[j] != S[mpos] ) break;
                              mpos++;    }
                         stops[u-low] = j;    }
                    int M = Max(stops);
                    vec<int> at;
                    for ( int j = 0; j < stops.isize( ); j++ )
                         if ( stops[j] == M ) at.push_back(j);
                    cout << "see: ";
                    for ( int r = 0; r < at.isize( ); r++ )
                    {    if ( r > 0 ) cout << "; ";
                         int j = at[r];
                         int id = X[ low + j ].second;
                         cout << "line " << lid[id] 
                              << "/" << GetOp( lid[id], inv, linesx )
                              << "[l=" << lines[id].size( ) << "]";
                         if( stops[j] < lines[id].isize( ) )
                              cout << ", stop = " << stops[j];    }
                    cout << "\n";
                    pos += M;    }

               else
               {    if ( low < high )
                    {    id = X[low].second;
                         lpos = 0;    }
                    if ( low == high )
                    {    int64_t lowz = LowerBound1( kmers_plus, x );
                         int64_t highz = UpperBound1( kmers_plus, x );
                         if ( lowz == highz )
                         {    cout << "Nope, not there at all." << endl;
                              Scram(1);    }
                         vec<int> stops( highz - lowz, 0 );
                         for ( int64_t u = lowz; u < highz; u++ )
                         {    int id = kmers_plus[u].second, mpos = pos;
                              int lpos = kmers_plus[u].third;
                              const String& t = lines[id];
                              int j;
                              for ( j = lpos; j < t.isize( ); j++ )
                              {    if ( t[j] != S[mpos] ) break;
                                   mpos++;    }
                              stops[u-lowz] = j;    }
                         int M = Max(stops);
                         vec<int> at;
                         for ( int j = 0; j < stops.isize( ); j++ )
                              if ( stops[j] == M ) at.push_back(j);
                         cout << "see: ";
                         for ( int r = 0; r < at.isize( ); r++ )
                         {    if ( r > 0 ) cout << "; ";
                              int j = at[r];
                              int id = kmers_plus[ lowz + j ].second;
                              cout << "line " << lid[id] 
                                   << "/" << GetOp( lid[id], inv, linesx )
                                   << "[l=" << lines[id].size( ) << "]" 
                                   << ", start = " << lpos;
                              if ( stops[j] < lines[id].isize( ) )
                                   cout << ", stop = " << stops[j];    }
                         cout << "\n";
                         pos += M;    }

                    else
                    {    const String& t = lines[id];
                         int j;
                         for ( j = lpos; j < t.isize( ); j++ )
                         {    if ( t[j] != S[pos] ) break;
                              pos++;    }
                         int stop = j;
                         cout << "see line " << lid[id] 
                              << "/" << GetOp( lid[id], inv, linesx )
                              << "[l=" 
                              << lines[id].size( ) << "]" << ", start = " << lpos;
                         if ( stop < lines[id].isize( ) )
                              cout << ", stop = " << stop; 
                         cout << endl;    }    }
               int ns = 0;
               {    while( pos < S.isize( ) && S[pos] == 'N' )
                    {    ns++;
                         pos++;    }    }
               if ( ns > 0 ) cout << "see gap of size " << ns << endl;    }    }
     cout << "\n" << Date( ) << ": done" << endl;    }
