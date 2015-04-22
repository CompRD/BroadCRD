///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FilterRegionalData.  Given data for a region of a genome, obtained by
// pulling the alignments to that region, filter the data in an attempt to 
// remove ectopically aligned reads.  This necessarily entails a very careful
// compromise: deleting all reads having high quality mismatches with the
// reference would on the one hand make the data too "easy", and on the other
// hand might create holes where the reference has an error (or the reference
// differs from the sample).

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "lookup/LookAlign.h"

int main( int argc, char** argv ) 
{
     RunTime();

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     EndCommandArguments;

     // Define directory.
  
     String data_dir = PRE + "/" + DATA;

     // Define region.

     vecbasevector genome( data_dir + "/../genome.fastb" );
     String region = FirstLineOfFile( data_dir + "/../region" );
     int tig_id = region.Before( ":" ).Int( );
     int origin = region.Between( ":", "-" ).Int( );

     // Load data, alignments, and mapping qualities.

     vecbasevector bases( data_dir + "/frag_reads_orig.fastb" );
     size_t nreads = bases.size( );
     vecqualvector quals( data_dir + "/frag_reads_orig.qualb" );
     PairsManager pairs( data_dir + "/frag_reads_orig.pairs" );
     vec<look_align> aligns;
     LoadLookAligns( data_dir + "/frag_reads_orig.qltout", aligns );
     vec<int> mapq;
     fast_ifstream mapin( data_dir + "/frag_reads_orig.mapq" );
     String line;
     while(1)
     {    getline( mapin, line );
          if ( mapin.fail( ) ) break;
          mapq.push_back( line.Int( ) );    }

     // Determine which reads are aligned.

     vec<Bool> aligned( nreads, False );
     for ( int i = 0; i < aligns.isize( ); i++ )
          aligned[ aligns[i].query_id ] = True;

     // Define reads to be removed:
     // (1) Reads whose partners are unaligned or aligned outside the region
     //     are deleted.
     // (2) Amongst mapping quality zero reads, and ignoring the first and last
     //     ten bases of each read, we find all Q30 mismatches that occur between
     //     2 and 9 times on the reference.  Reads having such bases are deleted.

     vec<Bool> to_remove( nreads, False );
     vec<int> mis( genome[0].size( ), 0 );
     int stop = origin + genome[0].isize( );
     const basevector& rd2 = genome[0];
     for ( int xpass = 1; xpass <= 2; xpass++ )
     {    for ( int i = 0; i < aligns.isize( ); i++ )
          {    const look_align& la = aligns[i];
               int id = la.query_id;
               int pid = pairs.getPartnerID(id);
               if ( pid < 0 || !aligned[pid] ) to_remove[id] = True;
               if ( mapq[i] > 0 ) continue;
               basevector rd1 = bases[id];
               qualvector q = quals[id];
               if ( la.Rc1( ) )
               {    rd1.ReverseComplement( );
                    q.ReverseMe( );    }
               const align& a = la.a;
               if ( a.pos2( ) < origin || a.Pos2( ) > stop ) continue;
               vec<int> dqs;
               int p1 = a.pos1( ), p2 = a.pos2( ) - origin;
               for ( int j = 0; j < a.Nblocks( ); j++ ) 
               {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
                    if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
                    for ( int x = 0; x < a.Lengths(j); x++ ) 
                    {    if ( rd1[p1] != rd2[p2] )
                         {    if ( p1 >= 10 && p1 < rd1.isize( ) - 10 
                                   && q[p1] >= 30 ) 
                              {    if ( xpass == 1 || mis[p2] >= 2 && mis[p2] < 10 )
                                        dqs.push_back( p2 );    }    }
                         ++p1; ++p2;    }    }
               Sort(dqs);
               if ( xpass == 1 )
               {    for ( int r = 0; r < dqs.isize( ); r++ )
                         mis[ dqs[r] ]++;    }
               if ( xpass == 2 && dqs.nonempty( ) )
               {    cout << "\nread " << id << ", mapq = " << mapq[i];
                    cout << "\nerrs at:";
                    for ( int j = 0; j < dqs.isize( ); j++ )
                         cout << " " << dqs[j];
                    cout << "\n";    
                    to_remove[id] = True;
                    if ( pid >= 0 ) to_remove[pid] = True;    }    }    }

     // Reads whose partner is unaligned will be removed.

     for ( int i = 0; i < aligns.isize( ); i++ )
     {    int id = aligns[i].query_id;
          int pid = pairs.getPartnerID(id);
          if ( pid < 0 || !aligned[pid] ) to_remove[id] = True;    }

     // Report results.

     cout << "\n" << PERCENT_RATIO( 3, Sum(to_remove), (int) bases.size( ) ) 
          << " of reads removed" << endl;

     // Remove reads.

     for ( size_t i = 0; i < nreads; i++ )
     {    if ( to_remove[i] )
          {    bases[i].resize(0), quals[i].resize(0);    }    }
     bases.WriteAll( data_dir + "/frag_reads_orig.fastb" );
     quals.WriteAll( data_dir + "/frag_reads_orig.qualb" );
     vec<Bool> pairs_to_remove( pairs.nPairs( ), False );
     for ( size_t i = 0; i < pairs.nPairs( ); i++ )
     {    if ( to_remove[ pairs.ID1(i) ] || to_remove[ pairs.ID2(i) ] )
               pairs_to_remove[i] = True;    }
     pairs.removePairs(pairs_to_remove);
     pairs.Write( data_dir + "/frag_reads_orig.pairs" );    }
