///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Temporary program to track run time in main loop of LongProto.

#include "MainTools.h"

int main( )
{
     String stage = "/wga/scr4/dexter/stage";
     vec<String> all = AllFiles(stage);
     vec<Bool> to_delete( all.size( ), False );
     for ( int i = 0; i < all.isize( ); i++ )
     {    int month = all[i].Before( "-" ).Int( );
          if ( month < 10 ) to_delete[i] = True;
          int day = all[i].Between( "-", "-" ).Int( );
          if ( month == 10 && day < 12 ) to_delete[i] = True;    }
     EraseIf( all, to_delete );
     for ( int i = 0; i < all.isize( ); i++ )
     {    String dir = stage + "/" + all[i] + "/longread";
          int month = all[i].Before( "-" ).Int( );
          int day = all[i].Between( "-", "-" ).Int( );
          String plasmo, mouse;
          if ( month < 11 || ( month == 11 && day <= 23 ) )
          {    plasmo = "Assembly_of_Plasmodium_from_simulated_long_reads";
               mouse = "Assembly_of_~230_Mb_mouse_region_from_simulated_long_reads";
                    }
          else
          {    plasmo = "Test_1,_Plasmodium_from_simulated_long_reads";
               mouse = "Test_2,_230_Mb_mouse_region_from_simulated_long_reads";    }
          String pfile = dir + "/" + plasmo + "/statistics0.txt";
          String mfile = dir + "/" + mouse + "/statistics0.txt";
          String pline = LineOfOutput( 
               "cat " + pfile + " | grep 'loop complete' | Col 12" );
          String mline = LineOfOutput( 
               "cat " + mfile + " | grep 'loop complete' | Col 12" );
          double phours = pline.Double( ) / 60.0, mhours = mline.Double( );
          cout << all[i] << "  plasmo = " << pline << " minutes"
               << "   mouse = " << mline << " hours" 
               << "   ratio = " << mhours/phours << endl;    }    }
