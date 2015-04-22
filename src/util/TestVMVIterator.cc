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
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

#include "STLExtensions.h"


typedef VirtualMasterVec<ReadPath> VirtualReadPathVec;

int main( int argc, char* argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(PATHS, "readpaths filename");
     EndCommandArguments;

     VirtualReadPathVec paths(PATHS);
#if 0
     for ( size_t i = 0; i < paths.size(); ++i ) {
          ReadPath const& p = paths[i];
          cout << "direct iterator vec size=(" << paths[i].size() << ")" << endl;
          std::copy(paths[i].begin(), paths[i].end(), std::ostream_iterator<int>(cout,",") );
          cout << endl;
          cout << "reference iterator vec size=(" << paths[i].size() << ")" << endl;
          std::copy(p.begin(), p.end(), std::ostream_iterator<int>(cout,","));
          cout << endl;
     }
#endif

     for ( auto itr=paths.begin(),end=paths.end(); itr != end; ++itr )
     {
         ReadPath const& p = *itr;
         cout << "direct iterator vec size=(" << itr->size() << ")" << endl;
         std::copy(itr->begin(), itr->end(), std::ostream_iterator<int>(cout,",") );
         cout << endl;
         cout << "reference iterator vec size=(" << p.size() << ")" << endl;
         std::copy(p.begin(), p.end(), std::ostream_iterator<int>(cout,","));
         cout << endl;
     }

     return 0;
}
