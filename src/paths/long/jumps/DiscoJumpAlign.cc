///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// DiscoJumpAlign.  Adapted from HicAlign:

// (via HicAlign)  Align Hi-C reads to an assembly graph.  To do so, we first
// truncate edges in the graph so that they overlap by only K = 40 bases.
// A read is declared placed if its first kmer appears uniquely in the
// assembly, after this truncation.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(JUMPHEAD,
             "reads jumps separately from <JUMPHEAD>_R{1,2}.{fastb,qualb}");
     CommandArgument_String_Doc(ADIR,
             "Specify a.fin(al) directory for input");
     CommandArgument_String_OrDefault_Doc(OUTPUTDIR, "",
             "specify output directory for (basename <JUMPHEAD>).aligns (default: ADIR)");
     CommandArgument_Int_OrDefault_Doc(OFFSET, 0,
             "offset from start of read for single kmer");
     EndCommandArguments;

     ForceAssertGe(OFFSET,0);

     if ( OUTPUTDIR == "" ) OUTPUTDIR = ADIR;

     // HARDCODED STUFF.

     // Heuristics.

     const int max_freq = 1;
     const int K = 40;
     const int Kbig = 200;

     cout << Date() << ": using K=" << K << ", assembly K=" << Kbig << ", max_freq=" << max_freq << endl;

     // Load data.

     cout << Date( ) << ": loading data" << endl;
     vecbasevector tigs( ADIR + "/a.fastb" );
     vecbasevector bases( JUMPHEAD+ "_R1.fastb" );
     int64_t npairs = bases.size( );
     bases.ReadAll( JUMPHEAD + "_R2.fastb", True );

     ForceAssertEq( bases.size(), 2*npairs );

     // Build kmers.

     cout << Date( ) << ": building kmers" << endl;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     double kclock = WallClockTime( );
     vec<int64_t> xstarts;
     xstarts.reserve(tigs.size()+1);
     xstarts.push_back(0);
     for ( size_t i = 0; i < tigs.size( ); i++ )
     {    const basevector& u = tigs[i];
          xstarts.push_back( 
               xstarts.back( ) + Max( 0, u.isize( ) - Kbig + 1 ) );    }
     kmers_plus.resize( xstarts.back( ) + bases.size( ) );
     cout << Date() << ": ...kmers from contigs" << endl;
     #pragma omp parallel for schedule(dynamic,1)
     for ( size_t i = 0; i < tigs.size( ); i++ )
     {    const basevector& u = tigs[i];
          kmer<K> x;
          for ( int j = 0; j <= u.isize( ) - Kbig; j++ )
          {    int64_t r = xstarts[i] + j;
               x.SetToSubOf( u, j ); 
               kmers_plus[r].first = x;
               kmers_plus[r].second = i; 
               kmers_plus[r].third = j;    }    }
     cout << Date() << ": ...kmers from reads" << endl;
     #pragma omp parallel for
     for ( size_t i = 0; i < bases.size( ); i++ )
     {    kmer<K> x;
          x.SetToSubOf( bases[i], OFFSET );
          int64_t r = xstarts.back( ) + i;
          kmers_plus[r].first = x;
          kmers_plus[r].second = (int64_t) tigs.size( ) + i; 
          kmers_plus[r].third = OFFSET;    }
     ParallelSort(kmers_plus);
     cout << TimeSince(kclock) << " used building kmers" << endl;

     // Find matches.

     cout << Date( ) << ": finding matches" << endl;

     // batch kmers for parallel work, but ensure that entries matching a
     // given kmers all lie within a single block
     const int64_t batches = 1000;
     vec<int64_t> bstart(batches+1);
     for ( int64_t i = 0; i <= batches; i++ )
          bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / batches;
     for ( int64_t i = 1; i < batches; i++ )
     {    int64_t& s = bstart[i];
          while( s > bstart[i-1] && kmers_plus[s].first == kmers_plus[s-1].first )
          {    s--;    }    }

     double mclock = WallClockTime( );
     vec< vec< triple<int,int,int> > > resultsx(batches);

     // process kmer blocks
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t bi = 0; bi < batches; bi++ )
     for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ )
     {    int64_t j, m;
          // after this block [i,j) will represent a given kmer
          for ( j = i + 1; j < bstart[bi+1]; j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          // find m is index of first non-contig (i.e. read)-based kmer, if any
          for ( m = i; m < j; m++ )
               if ( kmers_plus[m].second >= (int) tigs.size( ) ) break;
          // check that there's <= max_freq hits in the contigs
          // and that there's at least one hit in the reads
          if ( m - i <= max_freq && j - m >= 1 )
          {
              // i..m-1 for each contig hit
              for ( int64_t k1 = i; k1 < m; k1++ )
               {
                    int tid = kmers_plus[k1].second;
                    // m..j-1 for each read hit
                    for ( int64_t k2 = m; k2 < j; k2++ )
                    {    int rid = kmers_plus[k2].second - (int) tigs.size( );
                         // adjust read it to account for pairs and non-paired
                         // reads
                         if ( rid < npairs ) rid = 2 * rid;
                         else rid = 2 * ( rid - npairs ) + 1;
                         int offset = kmers_plus[k1].third - kmers_plus[k2].third;
                         // push read id, contig id, offset
                         resultsx[bi].push( rid, tid, offset );    }    }    }
          i = j - 1;    }
     cout << TimeSince(mclock) << " used matching kmers" << endl;

     // Collate results.

     cout << Date( ) << ": combining results" << endl;
     vec<int64_t> starts( batches+1, 0 );
     for ( int j = 0; j < batches; j++ )
          starts[j+1] = starts[j] + (int64_t) resultsx[j].size( );
     vec< triple<int,int,int> > results( starts.back( ) );
     #pragma omp parallel for
     for ( int j = 0; j < batches; j++ )
     {    memcpy( &results[ starts[j] ], &resultsx[j][0],
               resultsx[j].size( ) * sizeof( triple<int,int,int> ) );    }
     cout << Date( ) << ": sorting" << endl;
     ParallelSort(results);

     // results gets (read_id, tig_id, offset)
     // note that offset = tig_pos - read_pos = tig_pos - 0 = tig_pos

     // Done.

     cout << Date( ) << ": finding places" << endl;
     // places gets (tig_id, offset) for each read
     // we offset by OFFSET so that this is a location on the tig,
     // rather than a relative offset of the read position and the tig position
     vec< pair<int,int> > places( bases.size( ), make_pair( -1, -1 ) );
     for ( int i = 0; i < results.isize( ); i++ )
     {  ForceAssertLt(results[i].first, places.isize());
         places[ results[i].first ]
               = make_pair( results[i].second, results[i].third + OFFSET );    }
     cout << Date( ) << ": summarizing results" << endl;
     vec< pair< pair<int,int>, pair<int,int> > > PLACES;
     // PLACES gets pair< (tig_id,offset), (tig_id, offset)> implied
     // by each read pair
     for ( int pid = 0; pid < (int) bases.size( )/2; pid++ )
     {    int id1 = 2*pid, id2 = 2*pid+1;
          if ( places[id1].first < 0 || places[id2].first < 0 ) continue;
          PLACES.push( places[id1], places[id2] );
          /*
          {    cout << "[" << pid << "] " << places[id1]
                    << " --> " << places[id2] << endl;    }
          */
               }
     cout << Date( ) << ": sorting places" << endl;
     ParallelUniqueSort(PLACES);
     cout << "yield = " << PERCENT_RATIO( 3, PLACES.isize( ), npairs ) << endl;
     PRINT2( npairs, PLACES.size( ) );
     String output = OUTPUTDIR + "/" + Basename(JUMPHEAD) + ".aligns";
     cout << Date( ) << ": writing " << output << endl;
     BinaryWriter::writeFile( output, PLACES );
     cout << Date( ) << ": done" << endl;
     cout << "peak memory usage = " << ToStringAddCommas( PeakMemUsageBytes( ) )
          << endl;
     Scram(0);    }
