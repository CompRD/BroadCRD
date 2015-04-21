///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "paths/long/EvalAssembly.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/fosmid/Fosmids.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(N);
     CommandArgument_String_OrDefault(INSTANCE, "");
     EndCommandArguments;

     uint NUM_THREADS = 0;
     SetThreads(NUM_THREADS);

     if ( INSTANCE != "" ) INSTANCE += ".";

     vec<int> fos;
     if ( N == "all" ) fos = AllFosmids( );
     else ParseIntSet( N, fos );
 
     vecbasevector G;
     for ( int i = 0; i < fos.isize( ); i++ )
     {    int n = fos[i];
          vecbasevector g;
          FetchReads( g, 0, "/wga/scr4/jaffe/CompareVars/" + ToString(n) + "/fos."
               + ToString(n) + ".fasta" );
          G.Append(g);    }

     SupportedHyperBasevector shb;
     BinaryReader::readFile( "belch." + INSTANCE + "shbv", &shb );

     vec<int> inv2;
     FixInversion( shb, inv2 );

     vec<int> gaps, indels;
     int subs = 0;
     cout << "start loop" << endl;
     #pragma omp parallel for
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    
          // Use only graph components having a 100-mer match to the reference
          // sequence.  This is a kludge.

          vecbasevector rbases;
          rbases.push_back( G[g] );
          vecbasevector edges;
          for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
               edges.push_back( shb.EdgeObject(e) );
          vec<Bool> marked( edges.size( ), False );
          const int K = 100;
          ForceAssertEq( K, shb.K( ) );
          vecbasevector all(rbases);
          all.Append(edges);
          vec< triple<kmer<K>,int,int> > kmers_plus;
          MakeKmerLookup2( all, kmers_plus );
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
          vec<vec<int>> comp;
          shb.ComponentEdges(comp);
          vec<int> dels;
          for ( int i = 0; i < comp.isize( ); i++ )
          {    vec<int> c = comp[i];
               Bool ok = False;
               for ( int j = 0; j < c.isize( ); j++ )
                    if ( marked[ c[j] ] ) ok = True;
               if ( !ok )
               {    for ( int j = 0; j < c.isize( ); j++ )
                         dels.push_back( c[j] );    }    }
          SupportedHyperBasevector shb2(shb);
          shb2.DeleteEdges(dels);
          shb2.RemoveDeadEdgeObjects( );
          shb2.RemoveEdgelessVertices( );
          vec<int> inv2;
          FixInversion( shb2, inv2 );

          // Now do the evaluation.

          ref_data ref;
          ref.G.push_back( G[g] );
          ref.G3 = ref.G;
          ref.is_circular.resize( 1, False );
          CreateGlocs( ref.G, ref.LG, ref.Glocs );
          CreateGlocs( ref.G3, ref.LG, ref.G3locs );
          CreateGlocs( ref.G3plus, ref.LG, ref.G3pluslocs );
          vec<basevector> p;
          p.push_back( G[g] );
          const int KX = 80; // no justification for this
          ref.GH.push_back( HyperBasevector( KX, p ) );
          Ofstream( out, "some.notes." + INSTANCE + ToString( fos[g] ));
          AssemblyError err = EvalAssembly( ref, shb2, inv2, out, 0, 0 );
          #pragma omp critical
          {    gaps.append(err.GetGaps()), indels.append(err.GetIndels());
               subs += err.GetSubs();
               cout << "done with " << g+1 << endl;    }    }
     ReverseSort(gaps), ReverseSort(indels);
     cout << "gaps = {" << printSeq(gaps) << "}, indels = {" << printSeq(indels)
          << "}, subs = " << subs << endl;    }
