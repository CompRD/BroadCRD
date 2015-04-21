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
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/ImprovePath.h"

int main( int argc, char** argv )
{    RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(IDS, "read ids (ParseIntSet)");
     CommandArgument_Bool_OrDefault(SHOW_OLD_BETTER, False);
     CommandArgument_Bool_OrDefault(SHOW_NEW_BETTER, False);
     CommandArgument_Bool_OrDefault(IMPROVE_PATHS_LARGE, True);
     EndCommandArguments;

     // Directories.

     String work_dir = "/wga/scr4/jaffe/GapToy/50718";
     String fin_dir = work_dir + "/a.final";

     // Load assembly.

     // cout << Date( ) << ": loading assembly" << endl;
     HyperBasevector hb;
     BinaryReader::readFile( fin_dir + "/a.hbv", &hb );
     vec<int> inv;
     BinaryReader::readFile( fin_dir + "/a.inv", &inv );

     // Expand read set.

     vec<int64_t> ids;
     ParseLongLongSet( IDS, ids );
     int nids = ids.size( );
     for ( int i = 0; i < nids; i++ )
     {    int pid = ids[i]/2;
          ids.push_back( 2*pid, 2*pid + 1 );    }
     UniqueSort(ids);

     // Load subset of data.

     cout << Date( ) << ": loading some reads" << endl;
     ReadPathVec paths;
     vecbasevector bases;
     VecPQVec quals;
     paths.Read( fin_dir + "/a.paths", ids );
     bases.Read( work_dir + "/data/frag_reads_orig.fastb", ids );
     quals.Read( work_dir + "/data/frag_reads_orig.qualp", ids );

     // Map reads.

     cout << Date( ) << ": mapping reads" << endl;
     double iclock = WallClockTime( );
     int count_old_better = 0, count_new_better = 0;
     int count_same = 0, count_indet = 0;
     Bool track_results = False;
     path_improver pimp;
     if (SHOW_OLD_BETTER) pimp.show_old_better = True;
     if (SHOW_NEW_BETTER) pimp.show_new_better = True;
     ImprovePaths( paths, hb, inv, bases, quals, ids, pimp,
          IMPROVE_PATHS_LARGE, False );
     cout << "\n" << TimeSince(iclock) << " used in ImprovePath" << endl << endl;
     Scram(0);    }
