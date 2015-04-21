///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Test code for finding read paths that might need improvement.
//
// Messes with /wga/scr4/$user/GapToy/1.
//
// Scratch file: PathsTest60.aligns.

// MakeDepend: dependency GapToy 
// MakeDepend: dependency MakeLookupTable 
// MakeDepend: dependency QueryLookupTable 

#include "FastIfstream.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "paths/HyperBasevector.h"
#include "paths/long/PlaceReads1.h"
#include "paths/long/ReadPath.h"

String ToString( const ReadPath& p )
{    if ( p.size( ) == 0 ) return "";
     else 
     {    ostringstream s;
          String o = ( p.getOffset( ) > 0 ? "+" : "" ) + ToString( p.getOffset( ) );
          while( o.size( ) < 7 ) o += " ";
          s << o << printSeq(p);
          return s.str( );    }    }

int main( )
{
     RunTime( );

     // Number of reads to check.

     int n = 100;

     // Hardcoded path.

     String dir = "/wga/scr4/" + Getenv("USER") + "/GapToy/1";

     // Generate three sets of paths.

     SystemSucceed( "GapToy X=tiny EXIT=PATHS > /dev/null" );
     ReadPathVec paths1( dir + "/a.60/a.paths" );
     SystemSucceed( "GapToy X=tiny EXIT=PATHS NEW_ALIGNER=True > /dev/null" );
     ReadPathVec paths2( dir + "/a.60/a.paths" );
     vecbasevector bases( dir + "/data/frag_reads_orig.fastb" );
     cout << "\nThere are " << ToStringAddCommas( bases.size( ) ) << " reads." 
          << endl;
     vecqualvector quals( dir + "/data/frag_reads_orig.qualb" );
     HyperBasevector hb;
     BinaryReader::readFile( dir + "/a.60/a.hbv", &hb );
     ReadPathVec paths3;
     PlaceReads1( hb, bases, quals, paths3 );

     // Find reads for which the paths disagree.

     vec<int> ids0;
     for ( int id = 0; id < n; id++ )
     {    const ReadPath &p1 = paths1[id], &p2 = paths2[id], &p3 = paths3[id];
          if ( p1 != p2 || p2 != p3 ) ids0.push_back(id);    }

     // Get QLT alignments of the first 60-mer in each read to the assembly.

     SystemSucceed( "cd " + dir + "/a.60; MakeLookupTable SOURCE=a.fastb "
          "OUT_HEAD=a LO=True > /dev/null" );
     ostringstream vout0;
     vout0 << printSeq(ids0);
     SystemSucceed( "( cd " + dir + "/data; QueryLookupTable K=12 MM=12 "
          "MC=0.15 SEQS=frag_reads_orig.fastb QUALS=frag_reads_orig.qualb "
          "L=../a.60/a.lookup TRUNCATE_TO=60 SMITH_WAT=True REQUIRE_FULL1=True "
          "FW_ONLY=True SEQS_TO_PROCESS=\"{" + vout0.str( ) + "}\" PARSEABLE=True "
          " ) > PathsTest60.aligns" );
     vec<look_align> aligns;
     vec<vec<int>> aligns_index;
     LoadLookAligns( "PathsTest60.aligns", aligns, aligns_index, bases.size( ) );

     // Exclude some reads.

     vec<int> ids;
     for ( int i = 0; i < ids0.isize( ); i++ )
     {    int id = ids0[i];
          const ReadPath &p1 = paths1[id], &p2 = paths2[id], &p3 = paths3[id];

          // Exclude a case where we think PlaceReads1 is just insensitive.

          if ( p1 == p2 && p1.size( ) >= 1 && p3.size( ) == 0 )
          {    if ( aligns_index[id].solo( ) )
               {    const look_align& la = aligns[ aligns_index[id][0] ];
                    if ( la.target_id == p1[0] 
                         && p1.getOffset( ) == la.pos2( ) )
                    {    continue;    }    }    }

          // Exclude a case where we think Default is just insensitive.

          if ( p2 == p3 && p2.size( ) >= 1 && p1.size( ) == 0 )
          {    if ( aligns_index[id].solo( ) )
               {    const look_align& la = aligns[ aligns_index[id][0] ];
                    if ( la.target_id == p2[0] 
                         && p2.getOffset( ) == la.pos2( ) )
                    {    continue;    }    }    }

          // Otherwise keep the read.

          ids.push_back(id);    }

     // Print table of diagreements.

     cout << "\nDisagreements for reads = 0.." << n-1 << ":\n" << endl;
     vec<vec<String>> rows;
     vec<String> row;
     row.push_back( "ID", "DEFAULT", "NEW_ALIGNER", "PLACE1" );
     rows.push_back(row);
     for ( int i = 0; i < ids.isize( ); i++ )
     {    int id = ids[i];
          const ReadPath &p1 = paths1[id], &p2 = paths2[id], &p3 = paths3[id];
          vec<String> row;
          row.push_back( ToString(id), ToString(p1), ToString(p2), ToString(p3) );
          rows.push_back(row);    }
     PrintTabular( cout, rows, 4 );
     cout << "\n" << ids.size( ) << " aberrant reads" << endl;

     // Truncate the reads to length 60 and report their alignments versus the
     // assembly.

     ostringstream vout;
     vout << printSeq(ids);
     cout << "\nTruncate aberrant reads to length 60 and align to assembly:\n";
     SystemSucceed( "cd " + dir + "/data; QueryLookupTable K=12 MM=12 MC=0.15 "
          "SEQS=frag_reads_orig.fastb QUALS=frag_reads_orig.qualb "
          "L=../a.60/a.lookup VISUAL=True TRUNCATE_TO=60 SMITH_WAT=True NH=True "
          "QUIET=True REQUIRE_FULL1=True FW_ONLY=True SEQS_TO_PROCESS=\"{"
          + vout.str( ) + "}\"" );    }
