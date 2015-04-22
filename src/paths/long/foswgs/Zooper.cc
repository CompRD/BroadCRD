#include "MainTools.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"

int main( )
{    RunTime( );

     // Load shb.

     SupportedHyperBasevector shb;
     BinaryReader::readFile( "sss.27.1.shbv", &shb );

     // Set up logging etc.

     long_logging logc( "" );
     long_heuristics heur( "" );
     ref_data ref;
     vec<ref_loc> readlocs;
     String VERB = "", OUT_INT_HEAD = "";
     long_logging_control log_control( ref, &readlocs, OUT_INT_HEAD, VERB );

          {    cout << Date( ) << ": implementing REQUIRE_EDGE_MATCH" << endl;
               vecbasevector rbases( "/wga/scr4/jaffe/fos_filter/tmp.fos/"
                    "frag_reads_orig.fastb" );
               vecbasevector edges;
               for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
                    edges.push_back( shb.EdgeObject(e) );
               vec<Bool> marked( edges.size( ), False );
               const int K = 100;
               ForceAssertEq( K, shb.K( ) );
               vecbasevector all(rbases);
               all.Append(edges);
               vec< triple<kmer<K>,int,int> > kmers_plus;
               cout << Date( ) << ": making kmer lookup" << endl;
               MakeKmerLookup2( all, kmers_plus );
               cout << Date( ) << ": marking edges" << endl;
               for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
               {    int64_t j;
                    for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                         if ( kmers_plus[j].first != kmers_plus[i].first ) break;
                    Bool valid = False;
                    for ( int64_t k = i; k < j; k++ )
                    {    if ( kmers_plus[k].second < (int) rbases.size( ) ) 
                              valid = True;    }
                    if (valid)
                    {    for ( int64_t k = i; k < j; k++ )
                         {    if ( kmers_plus[k].second >= (int) rbases.size( ) ) 
                              {    int id = kmers_plus[k].second 
                                        - (int) rbases.size( );
                                   marked[id] = True;    }    }    }
                    i = j - 1;    }
               PRINT2( shb.EdgeObjectCount( ), Sum(marked) );
               vec<int> dels;
               const int min_kmers_to_del = 1;
               for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
               {    if ( !marked[e] && shb.EdgeLengthKmers(e) >= min_kmers_to_del ) 
                         dels.push_back(e);    }
               // cout << "deleting " << printSeq(dels) << endl;
               cout << "deleting " << dels.size( ) << " of "
                    << shb.EdgeObjectCount( ) << " edges" << endl;
               shb.DeleteEdges(dels);
               shb.RemoveUnneededVertices( );
               shb.RemoveDeadEdgeObjects( );
               shb.RemoveEdgelessVertices( );    }

     logc.USE_GENOME_FOR_DUMP = False;
     shb.DumpFiles( "alpha", log_control, logc );    }
