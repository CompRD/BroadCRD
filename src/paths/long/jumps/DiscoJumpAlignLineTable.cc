///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "system/System.h"
#include <omp.h>

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(JUMPHEAD,
             "reads aligns from <JUMPHEAD>.aligns");
     CommandArgument_String_Doc(ADIR,
             "Specify a.fin(al) directory for input");
     CommandArgument_String_OrDefault_Doc(OUTPUTDIR, "",
             "specify output directory for Basename(<JUMPHEAD>).line_aligns (default: ADIR)");
     EndCommandArguments;

     if ( OUTPUTDIR == "" ) OUTPUTDIR = ADIR;

     String input = JUMPHEAD + ".aligns";
     cout << Date() << input << endl;
     vec<pair<pair<int,int>, pair<int,int>>> places;
     BinaryReader::readFile(JUMPHEAD+".aligns", &places);
     cout << "read " << places.size() << " edge aligns " << endl;

     vec<int> inv;
     BinaryReader::readFile(ADIR+"/a.inv", &inv);
     cout << "involution covers " << inv.size() << " entries" << endl;;

     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile(ADIR+"/a.lines", &lines);
     cout << "read " << lines.size() << " lines" << endl;

     cout << Date() << ": computing line involution" << endl;
     vec<int> lineInv(lines.size(), -1);
     size_t goods = 0;
     size_t counter = 0;
#pragma omp parallel for schedule(static,10000)
     for ( int i = 0; i < lines.isize()-1; ++i ) {
         if ( lineInv[i] != -1 || lines[i].size() == 0 ) continue;
         for ( int j = i+1; j < lines.isize(); ++j ) {
             if ( lineInv[j] != -1 || lines[j].size() == 0 ) continue;
             if ( lines[i].front()[0][0] == inv[lines[j].back()[0][0]] ) {
                 lineInv[i] = j;
                 lineInv[j] = i;
                 goods++;
                 break;
             }
         }
     }

#if 0
     // check that all lines with size > 0 are covered
     for ( int i = 0; i < lines.isize(); ++i )
         if ( lines[i].size() > 0 ) ForceAssertGt(lineInv[i],-1);
#endif

     cout << Date() << ": found " << 2*goods << " line involution matches " << endl;

     cout << Date() << ": computing edge to CANONICAL line map" << endl;
     vec<int> edgeToLine(inv.size(), -1);

#pragma omp parallel for
     for ( int i = 0; i < lines.isize(); ++i ) {
         if ( lineInv[i] == -1 ) continue;
         for ( auto const& seg : lines[i] )
             for ( auto const& path : seg )
                 for ( auto const edge : path ) {
                     ForceAssertLt(edge, edgeToLine.isize());
                     edgeToLine[edge] = i;
                 }
     }


     cout << Date() << ": computing hits " << endl;
     vec<std::map<int, size_t>> lineHits(lines.size());

     for ( auto const& hit : places ) {
         auto e1 = hit.first.first;
         auto e2 = hit.second.first;
         ForceAssertLt(e1, edgeToLine.isize() );
         ForceAssertLt(e2, edgeToLine.isize() );
         if ( edgeToLine[e1] != -1 && edgeToLine[e2] != -1 ) {
             auto l1 = edgeToLine[e1];
             auto l2 = edgeToLine[e2];
             auto il1 = lineInv[l1];
             auto il2 = lineInv[l2];
             ForceAssertLt(l1, lineHits.isize());
             ForceAssertLt(l2, lineHits.isize());
             ForceAssertLt(il1, lineHits.isize());
             ForceAssertLt(il2, lineHits.isize());
             lineHits[l1][il2]++;
             lineHits[il2][l1]++;
             lineHits[il1][l2]++;
             lineHits[l2][il1]++;
         }
     }

     String output = OUTPUTDIR + "/" + Basename(JUMPHEAD) + ".line_hits";
     cout << Date() << ": writing table " << endl;
     Ofstream(TABLE,output);

     for ( int i = 0; i < lines.isize(); ++i ) {
         if ( lineInv[i] == -1 ) continue;
         TABLE << i << ": ";
         for ( auto const& hit : lineHits[i] )
             TABLE << hit.first << "[" << hit.second << "]  ";
         TABLE << endl;
     }

     return 0;
}
