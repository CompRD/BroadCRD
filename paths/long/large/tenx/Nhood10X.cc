///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Find neighbors of a given line.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: library GMPXX

// require <cstddef> to use gmp in GCC-4.9
#include <cstddef>
#include <gmpxx.h>

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/tenx/TenxDirs.h"
#include "paths/long/large/tenx/TenxTools.h"

typedef mpf_class big_float;

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(L, "line ids to test");
     CommandArgument_Bool_OrDefault_Doc(S, False, "show surprise");
     CommandArgument_Int_Doc(N, "1 or 2 or 3 or 4 or 5 or 6");
     CommandArgument_Double_OrDefault_Doc(MIN_DEVS, 25,
          "minimum surprise to call proximity; the default is appropriate for "
          "the NA12878 dataset but probably not for others");
     CommandArgument_Int_OrDefault_Doc(LF, -1, "check this line for proximity");
     EndCommandArguments;

     // Hardcoded directories.

     String dir, odir, tdir;
     SetTenxDirs( N, dir, odir, tdir );

     // Control.

     const Bool avoid_ends = False;

     // Load counts and assembly.

     cout << "\n" << Date( ) << ": loading assembly" << endl;
     HyperBasevectorX hb;
     vec<int> inv, tol, lens;
     vec<vec<vec<vec<int>>>> lines;
     vec< vec< pair<int,int> > > aligns;
     #pragma omp parallel sections
     {    
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.hbx", &hb );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.inv", &inv );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.lines", &lines );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.aligns", &aligns );    }    }
     GetTol( hb, lines, tol );
     GetLineLengths( hb, lines, lens );
     vec<String> genome_names;
     fast_ifstream in( dir + "/../genome.names" );
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          genome_names.push_back(line);    }

     // Load alignments and barcodes.

     cout << Date( ) << ": loading alignments" << endl;
     vec< pair<int,int> > places;
     BinaryReader::readFile( odir + "/10X.aligns", &places );
     cout << Date( ) << ": loading line hits" << endl;
     vec<double> bc_frac;
     BinaryReader::readFile( odir + "/10X.bc_frac", &bc_frac );

     // Load special structures.

     cout << Date( ) << ": loading special data structures" << endl;
     vec< vec<int> > lhitsb;   // line to barcodes
     vec<vec<int>> lbc;        // barcode to lines
     BinaryReader::readFile( odir + "/10X.lhitsb", &lhitsb );
     BinaryReader::readFile( odir + "/10X.lbc", &lbc );

     // Go through the lines.

     cout << Date( ) << ": finding the lines" << endl;
     vec<int> ls;
     ParseIntSet( L, ls );
     vec<Bool> to_delete( ls.size( ), False );
     for ( int j = 1; j < ls.isize( ); j++ )
     {    if ( lines[ ls[j] ].front( )[0][0] 
               == inv[ lines[ ls[j-1] ].back( )[0][0] ] )
          {    to_delete[j] = True;    }    }
     EraseIf( ls, to_delete );
     vec<String> reports( ls.size( ) );
     #pragma omp parallel for
     for ( int li = 0; li < ls.isize( ); li++ )
     {    int L = ls[li];

          // Find the lines.

          const vec<int>& xb = lhitsb[L];
          vec<int> M;
          for ( int i = 0; i < xb.isize( ); i++ )
               M.append( lbc[ xb[i] ] );
          Sort(M);
          vec< pair<int,int> > lhn;
          for ( int i = 0; i < M.isize( ); i++ )
          {    int j = M.NextDiff(i);
               lhn.push( j - i, M[i] );
               i = j - 1;    }
     
          // Now do the math.

          ostringstream out;
          out << "\n";
          int prints = 0;
          for ( int i = 0; i < lhn.isize( ); i++ )
          {    int l = lhn[i].second, l2 = -1;
               if ( lens[l] < 1000 ) continue;
               Bool pair = False;
               if ( i < lhn.isize( ) - 1 )
               {    l2 = lhn[i+1].second;
                    if ( lines[l].front( )[0][0] == inv[ lines[l2].back( )[0][0] ] )
                         pair = True;    }
               vec<String> chrs;
               for ( int j = 0; j < lines[l].isize( ); j++ )
               {    for ( int r = 0; r < lines[l][j].isize( ); r++ )
                    for ( int s = 0; s < lines[l][j][r].isize( ); s++ )
                    {    int e = lines[l][j][r][s];
                         for ( int m = 0; m < aligns[e].isize( ); m++ )
                              chrs.push_back( genome_names[aligns[e][m].first] );
                         e = inv[e];
                         for ( int m = 0; m < aligns[e].isize( ); m++ )
                         {    chrs.push_back( genome_names[aligns[e][m].first] );
                                  }    }    }
               UniqueSort(chrs);
               int obs = lhn[i].first;
               // if ( obs <= 1 ) continue;

               int nbc = xb.size( );
               double expect = nbc * bc_frac[l];

               double devs = (obs-expect) / sqrt(expect);
               if ( devs < MIN_DEVS && l != LF && l2 != LF ) continue;

               // Let sum = P( X >= obs ) for a a Poisson random variable with
               // lambda = expect.
          
               long double surprise;
               if (S)
               {
               int precision = 100;
               big_float sum = 0;
               for ( int x = 0; x < obs; x++ )
               {    
                    // Let s1 = expect^x.
     
                    big_float s1( 1, precision );
                    for ( int j = 0; j < x; j++ )
                         s1 *= expect;
     
                    // Let s2 = exp(-expect).  Note funny use of 100 term sum.
     
                    big_float s2( 1, precision );
                    big_float fact( 1, precision );
                    big_float pow( 1, precision );
                    for ( int j = 1; j < 100; j++ )
                    {    pow *= -expect;
                         fact *= j;
                         s2 += pow/fact;    }

                    // Let s = s1 * s2 / x!.

                    big_float s = s1 * s2;
                    for ( int j = 1; j <= x; j++ ) 
                         s /= j;
                    sum += s;    }
     
               // Let surprise = -log10(1-sum).
     
               big_float msum = 1 - sum;
               surprise = -log10l( msum.get_d( ) );
               // if ( surprise < 15.0 ) continue;
               }

               out << "[#" << ++prints << ", " << obs << " barcodes, " 
                    << expect << " expect, "
                    << devs << " devs";
               if (S) out << setprecision(3) << ", " << surprise << " surprise";
               out << "] L" << l;
               if ( !pair ) out << " (l=" << lens[l] << ")";
               else
               {    out << "/" << l2 << " (l=" << lens[l] << ")";
                    i++;    }
               if ( chrs.nonempty( ) ) out << " chrs = " << printSeq(chrs);
               out << "\n";    }
          reports[li] = out.str( );    }

     // Print reports.

     for ( int i = 0; i < reports.isize( ); i++ )
          cout << reports[i];
     cout << "\n" << Date( ) << ": done\n" << endl;
     Scram(0);    }
