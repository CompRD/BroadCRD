///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Special test code to run after a local assembly crash in GapToy.

#include "Basevector.h"
#include "MainTools.h"
#include "efasta/EfastaTools.h"
#include "paths/long/Heuristics.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/Logging.h"
#include "paths/long/LongHyper.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/RefTrace.h"
#include "paths/long/SupportedHyperBasevector.h"

int main( )
{
     RunTime( );
     uint NUM_THREADS = 1;
     SetThreads(NUM_THREADS);
     bool USE_OLD_LRP_METHOD = true;

     // Hardcoded paths.

     String run_head, tmp_dir1, tmp_dir2, out_dir;
     String user = Getenv( "USER" );
     if ( user == "jaffe" )
     {    run_head = "/wga/scr4/jaffe/tmp";
          tmp_dir1 = "/wga/dev/jaffe/BroadCRD/tmp.xxx";
          tmp_dir2 = "/wga/dev/jaffe/BroadCRD/tmp.micro";
          out_dir = "/wga/dev/jaffe/BroadCRD";    }
     else if ( user == "iainm" )
     {    run_head = "/wga/scr4/iainm/gaptoy/test";
          tmp_dir1 = "tmp.xxx";
          tmp_dir2 = "tmp.micro";
          out_dir = "/wga/scr4/iainm/gaptoy.out";    }
     else
     {    cout << "Unknown user." << endl;
          Scram(1);    }

     // Rerun assemblies.

     for ( int mu = 0; mu < 100; mu++ )
     {    PRINT(mu);
          String TMP = tmp_dir2 + "/" + ToString(mu);
          if ( !IsDirectory(TMP) ) continue;
          long_logging logc( "", "" );
          logc.STATUS_LOGGING = False;
          logc.MIN_LOGGING = False;
          ref_data ref;
          vec<ref_loc> readlocs;
          long_logging_control log_control( ref, &readlocs, "", "" );
          long_heuristics heur( "" );
          VecEFasta corrected;
          vec<pairing_info> cpartner;
          vec<int> cid;
          vecbasevector creads;
          double clock = WallClockTime( );
          CorrectionSuite( TMP, heur, logc, log_control, creads, corrected,
               cid, cpartner, NUM_THREADS, "", clock, USE_OLD_LRP_METHOD );
          int count = 0;
          for ( int l = 0; l < (int) corrected.size( ); l++ )
               if ( corrected[l].size( ) > 0 ) count++;
          if ( count == 0 ) continue;
          vec<Bool> to_delete( corrected.size( ), False );
          DefinePairingInfo(TMP, creads, to_delete, cid, corrected, cpartner, logc);
          SupportedHyperBasevector shb;
          if ( LongHyper( "", corrected, cpartner, shb, heur,
               log_control, logc, TMP, USE_OLD_LRP_METHOD ) )
          {    shb.DeleteLowCoverage( heur, log_control, logc );
               if ( !shb.NPaths( ) != 0 ) shb.TestValid(logc);    }    }    }
