// Copyright (c) 2000, 2001 Whitehead Institute for Biomedical Research

// ATWindows: complete fraction of genome in AT-rich windows.

#include "Basevector.h"
#include "math/Functions.h"
#include "MainTools.h"
#include <fstream>

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(SUBDIR);
     CommandArgument_UnsignedInt_OrDefault(GAP_FLOOR, 100);
     CommandArgument_UnsignedInt_OrDefault(WINDOW_SIZE, 100);
     CommandArgument_UnsignedInt_OrDefault(BUFFER, 0);
     CommandArgument_UnsignedInt_OrDefault(MIN_AT_PERCENT, 70);
     CommandArgument_String_OrDefault(ANNOT, "at_rich.ano");
     EndCommandArguments;

     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/" + SUBDIR;

     int window_size = WINDOW_SIZE;
     int buffer = BUFFER;

     // Set up standard data structures.

     vecbasevector tigs;
     tigs.ReadAll( sub_dir + "/mergedcontigs.fastb" );

     longlong windows = 0, AT_rich_windows = 0;

     String annot = sub_dir + "/" + ANNOT;
   
     std::ofstream os(annot.c_str());

     for ( int i = 0; i < (int) tigs.size( ); i++ )
       {    
	 int start = -1;
	 int end = -1;
	 for ( int j = buffer; j < (int) tigs[i].size( ) - buffer - window_size; 
	       j += window_size/10 )
	   {    
	     static basevector b;
	     b.SetToSubOf( tigs[i], j, window_size );
	     int AT = 0;
	     for ( int k = 0; k < window_size; k++ )
	       {
		 unsigned char c = b[k];
		 if ( as_base(c) == 'A' || as_base(c) == 'T' ) ++AT;    
	       }
	     if ( float(AT)/float(window_size) >= float(MIN_AT_PERCENT)/100.0 ) {
	       ++AT_rich_windows;
	       if (start == -1)
		 start = j;
	       end = j + window_size;
	       //cout << i << "." << j << "-" << j+1000 << "\n";    
	     } else {
	       if (start != -1) {
	           os << i << '\t' << start << '\t' << end << "\tat_rich\n";
	       }
	       start = -1;
	       end = -1;
	     }

	     

	     ++windows;    
	   }    

	 if (end != -1) {
             os << i << '\t' << start << '\t' << end << "\tat_rich\n";
	 }
       }
     os.close();
     ForceAssert(os);
     cout << "\n" << AT_rich_windows << " (" << setprecision(3) 
          << 100.0 * float(AT_rich_windows)
       / float(windows) << "%) of " << window_size 
          << "bp windows are >= " << MIN_AT_PERCENT << "% AT\n";    
}
