///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "paths/long/SupportedHyperBasevector.h"

int main( )
{
     RunTime( );

     SupportedHyperBasevector shb;

     BinaryReader::readFile( "/wga/scr4/macro/scs.shbv", &shb );
     vecbasevector bases( "/wga/scr4/macro/ecoli.Feb7.1.fastb" );

     /*
     // r48688:LongProto SAMPLE=scardovia READS=#picard TMP=tmp.xxx OUT_INT_HEAD=...
     BinaryReader::readFile( "/wga/scr4/macro/scardo.shbv", &shb );
     vecbasevector bases( "/wga/scr4/macro/scardo.Feb7.1.fastb" );
     */

     bases.WriteAll( "thereads.fastb" );
     int nb = bases.size( );
     vecbasevector bases2(bases);
     bases2.Append(bases);
     for ( int i = nb; i < 2*nb; i++ )
          bases2[i].ReverseComplement( );
     bases = bases2;

     // Heuristics.

     const int max_dups = 4;

     // Define paths.

     vec<int> to_left, to_right;
     shb.ToLeft(to_left), shb.ToRight(to_right);
     int maxread = 0;
     for ( int i = 0; i < nb; i++ )
          maxread = Max( maxread, bases[i].isize( ) );
     vec<vec<int>> allpaths;
     for ( int v = 0; v < shb.N( ); v++ )
     {    vec<vec<int>> paths;
          for ( int j = 0; j < shb.From(v).isize( ); j++ )
          {    vec<int> x = { shb.EdgeObjectIndexByIndexFrom( v, j ) };
               paths.push_back(x);    }
          while(1)
          {    vec<vec<int>> paths2;
               for ( int i = 0; i < paths.isize( ); i++ )
               {    int v = to_right[ paths[i].back( ) ];
                    int tail_len = 0;
                    for ( int j = 1; j < paths[i].isize( ); j++ )
                         tail_len += shb.EdgeLengthKmers( paths[i][j] );
                    if ( tail_len >= 2*maxread || shb.From(v).empty( ) )
                    {    paths2.push_back( paths[i] );
                         continue;    }
                    for ( int l = 0; l < shb.From(v).isize( ); l++ )
                    {    vec<int> p = paths[i];
                         p.push_back( shb.EdgeObjectIndexByIndexFrom( v, l ) );
                         vec<int> q = p;
                         Sort(q);
                         int mc = 0;
                         for ( int r = 0; r < q.isize( ); r++ )
                         {    int s = q.NextDiff(r);
                              mc = Max( mc, s - r );
                              r = s - 1;    }
                         if ( mc <= max_dups ) paths2.push_back(p);    }    }
               if ( paths2 == paths ) break;
               paths = paths2;    }
          allpaths.append(paths);    }
     PRINT( allpaths.size( ) );
     vecbasevector basepaths;
     for ( int i = 0; i < allpaths.isize( ); i++ )
          basepaths.push_back( shb.Cat( allpaths[i] ) );
     int nall = 0;
     for ( int i = 0; i < allpaths.isize( ); i++ )
          nall += basepaths[i].isize( );
     PRINT(nall);

     // Align the reads to the paths.

     basepaths.WriteAll( "uber.fastb" );
     SystemSucceed( 
          "MakeLookupTable SOURCE=uber.fastb OUT_HEAD=uber LO=True NH=True" );
     const int nbatches = 100;
     SystemSucceed( "/bin/rm -f threads.*" );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int batch = 0; batch < nbatches; batch++ )
     {    int start = (batch*nb)/nbatches, stop = ((batch+1)*nb)/nbatches; 
          SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.005 SEQS=thereads.fastb "
               "L=uber.lookup SMITH_WAT=True PARSEABLE=True REQUIRE_FULL1=True "
               "SEQS_TO_PROCESS=\"[" + ToString(start) + "," + ToString(stop) 
               + ")\" > thereads.aligns." + ToString(batch) );    }
     SystemSucceed( "cat thereads.aligns.* > thereads.aligns" );
     vec<look_align> aligns;
     vec<vec<int>> aligns_index;
     LoadLookAligns( "thereads.aligns", aligns, aligns_index, nb );

     // Go through the reads.

     for ( int id = 0; id < nb; id++ )
     {    
          // Sort the paths by their alignment errors, and take only the best.

          vec<int> errs, pids;
          for ( int j = 0; j < aligns_index[id].isize( ); j++ )
          {    const look_align& la = aligns[ aligns_index[id][j] ];
               errs.push_back( la.Errors( ) );
               pids.push_back( la.target_id );    }
          SortSync( errs, pids );
          vec<vec<int>> paths;
          for ( int j = 0; j < errs.isize( ); j++ )
          {    if ( j > 0 && errs[j] > errs[j-1] ) break;
               paths.push_back( allpaths[ pids[j] ] );    }
          if ( paths.empty( ) ) continue;

          // Find the longest subpath(s) in common to the best paths.

          vec<vec<int>> commons;
          vec<int> lens;
          for ( int start = 0; start < paths[0].isize( ); start++ )
          for ( int stop = start+1; stop <= paths[0].isize( ); stop++ )
          {    vec<int> p;
               for ( int j = start; j < stop; j++ )
                    p.push_back( paths[0][j] );
               Bool shared = True;
               for ( int j = 1; j < paths.isize( ); j++ )
                    if ( !paths[j].Contains(p) ) shared = False;
               if (shared) 
               {    commons.push_back(p);
                    int len = 0;
                    for ( int l = 0; l < p.isize( ); l++ )
                         len += shb.EdgeLengthKmers( p[l] );
                    lens.push_back(len);    }    }
          ReverseSortSync( lens, commons );
          cout << "\nbest paths for read " << id << "(len=" << bases[id].size( ) << "):\n";
          for ( int i = 0; i < lens.isize( ); i++ )
          {    if ( i > 0 && lens[i] < lens[i-1] ) break;
               cout << printSeq( commons[i] ) << endl;    }    }    }
