///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "paths/long/large/PullAparter.h"

int main( int argc, char* argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(IN_HEAD,
               "use IN_HEAD + {.hbv,.paths,.inv,.paths.inv}"  );
     CommandArgument_String_Doc(OUT_HEAD,
               "write OUT_HEAD.hbv" );
     CommandArgument_Int_OrDefault_Doc(EDGEID, -1, "specify (center) edge id of putative canonical repeat");
     EndCommandArguments;

     const String& hbvfn = IN_HEAD + ".hbv";
     const String& pathsfn = IN_HEAD + ".paths";
     const String& invfn = IN_HEAD + ".inv";
     const String& pathsinvfn = IN_HEAD + ".paths.inv";

     cout << "loading " << hbvfn << endl;
     HyperBasevector hbv(hbvfn);

     cout << "loading " << pathsfn << endl;
     ReadPathVec paths(pathsfn);

     cout << "loading path index " << pathsinvfn << endl;
     VecULongVec pathsinv(pathsinvfn);

     cout << "loading inversion " << invfn << endl;
     vec<int> inv;
     BinaryReader::readFile(invfn, &inv);

     vec<Bool> used;
     hbv.Used(used);
     int nused = Sum(used);
     int total = hbv.EdgeObjectCount();
     if ( total - nused > 0 ) {
          cout << "WARNING: there are " << total - nused << " dead edges" << endl;
          for ( size_t i = 0; i < used.size(); ++i )
               if ( !used[i] ) cout << i << endl;
     } else
          cout << "ALL " << total << " edges are used." << endl;

     PullAparter pa( hbv, inv, paths, pathsinv, vec<int>{}, (EDGEID >= 0) /* debug */ );

     vec<vec<int>> to_separate;
     if ( EDGEID >= 0 ) {
          if ( pa.isSeparable(EDGEID, &to_separate) ) {
               cout << "we *can* pull apart " << EDGEID << endl;
               for ( auto const& sep : to_separate ) cout << printSeq(sep) << endl;
          }
          else
               cout << "we *cannot* pull apart " << EDGEID << endl;
     } else {
          for ( int edgeid = 0; edgeid < hbv.EdgeObjectCount(); ++edgeid ) {
               if ( edgeid < inv[edgeid] ) pa.isSeparable(edgeid, &to_separate);
          }
          cout << Date() << ": done searching ALL edges" << endl;
     }
     cout << Date() << ": separating " << to_separate.size() << " paths " << endl;
     pa.SeparateAll(to_separate);
     cout << Date() << ": done.  There were " <<
               pa.getRemovedReadPaths() << " nullified read paths" << endl;

     cout << Date() << ": writing " << OUT_HEAD+".hbv" << endl;
     BinaryWriter::writeFile( OUT_HEAD + ".hbv", hbv);
     cout << Date() << ": done" << endl;

     return 0;
}
