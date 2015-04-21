/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

const char *DOC =

"Creates GC plots for the specified alignments. One separate curve will be created "
"for each comma-separated element of the HEAD list. The element of the HEAD list can "
"be either a separate file head, or a '+' separated list of file heads (\"file pile\"); "
"in the latter case, all the alignments from all the files in the pile will be merged "
"to plot a single GC curve. This tool expects to find <head>.qltout file with alignments "
"and <head>.READ_EXT (default READ_EXT=fastb) with reads for each head specified. "
"Generated GC plots show the ratios of frequencies of different GC counts found in windows "
"anchored at alignment positions on the reference to the total frequencies of corresponding "
"GC counts found across all windows on the reference; the plots are normalized so that the "
"ratio at GC count=50% of the window size is equal to 1."

;

#include "bias/BiasGC.h"
#include "BasevectorTools.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"
#include "lookup/AlignCollector.h"
#include "lookup/PerfectCount.h"
#include "math/Functions.h"
#include "random/Random.h"
#include "solexa/SolexaTools.h"
#include "solexa/SolexaMetrics.h"
#include "math/SparseArray.h"
// #include <math.h>

#define BRACKFAIL                             \
{    cout << "Illegal bracketing in HEAD.\n"; \
     cout << "Abort.\n";                      \
     exit(1);    }


int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandDoc(DOC);
     CommandArgument_String_Doc(HEAD, 
          "Comma and plus-sign separated list of file name heads (without an extention) - "
	  "each continuous group of files separated by '+' signs will be used to draw "
	  "a single curve (data accumulated across such files); a separate curve will be "
          "drawn for each comma separated file/group of files. The names provided here "
	  "are interpreted as relativ to DIR if they do not start with '/'");
     CommandArgument_String_OrDefault_Doc(DIR,".",
		       "Working directory the file names are specified relative to. If "
		       "a file name is fully qualified (starts with '/'), then DIR is ignored "
					  "for that file.");
     CommandArgument_StringSet_OrDefault_Doc(REF,"",
				"Reference fastb file(s), relative to DIR or absolute path."
				"The number of specified reference files should be either one, "
				"in which case the same reference is used for all curves, or "
				"it should be the same as the number of curves (e.g. files/groups "
				"separayted by ',' in the HEAD argument. If REF is not specified, "
                                "the references are sought as <head>.reference.fastb for each "
				"specified file");
     CommandArgument_StringSet_OrDefault_Doc(TITLES,"",
				"Allows to manually set the title for each curve; if not "
				"specified, the default will be the base part of each file head "
				"(the part after the last found '/' symbol), the title for a "
				"pile will be <bh1>+<bh2>+... etc, where bh_i are corresponding "
					     "base names of the files in the pile");
     CommandArgument_String_OrDefault_Doc(READ_EXT,"fastb",
					  "Extension of binary file storing the reads, defaults to '.fastb'");
     CommandArgument_Int_OrDefault_Doc(W, 50, "window size");
     CommandArgument_String_OrDefault_Doc(OUT, "", "output file, may end with .png");
     CommandArgument_Double_OrDefault_Doc(MAX_SLOP, 0.05, "spread of 95% "
          "confidence interval must be less than this to plot point" );
     CommandArgument_Bool_OrDefault_Doc(DOGRAPH, True,
          "If True, output graph." );
     CommandArgument_Bool_OrDefault_Doc(PRINTPOINTS, False,
          "If True, print to stdout the plot points." );
     CommandArgument_Bool_OrDefault_Doc(USE_ALIGNS, False,
          "If True, instead of using reads perfect to N bases, use aligned reads.");
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     CommandArgument_Bool_OrDefault(SAVE_MEMORY,False);
     CommandArgument_Int_OrDefault_Doc(MAX_ERRORS,infinitely_many,
				       "Only the aligns with at most MAX_ERRORS errors will"
				       "be used if USE_ALIGNS=True. If USE_ALIGNS=False, this "
				       "option is ignored.");
     EndCommandArguments;



     // Check arguments.
     ForceAssertEq( W % 2, 0 );
     if (DOGRAPH) CheckValidGraphicsSuffix( OUT, True );
     if ( ! USE_ALIGNS && (unsigned int)MAX_ERRORS < infinitely_many ) {
       cerr << "WARNING: MAX_ERRORS option is specified but will be ignored since USE_ALIGNS=False" << endl;
     }

     if ( DIR[DIR.size() - 1] != '/' ) DIR+="/"; //make sure DIR is '/'-terminated


     // process HEAD argument; the argument is expected to have the form
     // head{+head}{,head{+head}}; '+'-separated heads form a "pile" and
     // each pile will be drawn as one curve.
     vec<vec<String> > HEADS;
     vec<String> HEADPILE;
     vec<String> actual_refs;
     vec<String> actual_titles;
     String filehead;
     int j = 0; // remembers previous position
     HEAD+=","; // make HEAD ','-terminated, then everything can be processed
                // in one loop below, without a special case when the string runs out:
     for ( int i = 0 ; i < HEAD.isize() ; i++ ) {
       if ( HEAD[i] == ',' || HEAD[i] == '+' ) { // file name ended
	 if  ( i==j ) {
	   cerr << "ERROR:: file name expected at position " << i 
		<< " in the HEAD argument, empty string found" << endl;
	   exit(1);
	 }
	 filehead = HEAD.substr(j, i - j) ; // store current file head
	 if ( filehead[0] != '/' ) { // is the path fully qualified?
	   filehead = DIR+filehead;
	 }
	 // if default reference is requested, check that it is 
	 // available:
	 if ( REF.size() == 0 ) {
	   if ( ! IsRegularFile( filehead+".reference.fastb" ) ) {
	     cerr << "ERROR::default reference " << filehead << ".reference.fastb" 
		  << "does not exist. Try specifying the reference(s) explicitly." << endl;
	     exit(1);
	   }
	 }

	 // if default titles are requested, add this filehead's base:
	 if ( TITLES.size() == 0 ) {
	   if ( HEADPILE.size() == 0 ) { //first file in a pile
	     actual_titles.push_back( filehead.RevAfter("/") );
	   } else { // pile continues, update the title:
	     actual_titles.back() += "+";
	     actual_titles.back() += filehead.RevAfter("/");
	   }
	 }

	 HEADPILE.push_back( filehead );

	 if ( HEAD[i] == ',' ) { // if file pile ended
	   HEADS.push_back( HEADPILE );
	   // if default references are used, make sure that all the files
	   // in the pile have the same reference, otherwise grouping them together
	   // makes no sense!
	   if ( REF.size() == 0 ) {
	     RequireSameRef( HEADPILE );
	     actual_refs.push_back(HEADPILE[0]+".reference.fastb");
	   }
	   HEADPILE.clear();
	 }
	 j=i+1; // set start position to next symbol after the ',' or '+'
	 continue; // go get next symbol
       } 
     }
     
     // check ref argument if it was specified (if default is used it was already checked above):
     if ( 0 != REF.size() && 1 != REF.size() && REF.size() != HEADS.size() ) {
       cerr << "ERROR::wrong number of references specified. " << endl
	    << "use either one reference total (will be used for all alignments) " << endl
	    << "or one reference for each file/file pile in the HEAD argument" << endl;
       exit(1);
     } 
     if ( REF.size() == 1 ) {
       actual_refs.resize(HEADS.size(),REF[0]+".fastb"); // copy same ref for all heads
     } else {
       if (REF.size() != 0 ) {
	 actual_refs.resize(REF.size());
	 for ( int k = 0 ; k < REF.isize() ; k++ ) {
	   actual_refs[k] = REF[k]+".fastb"; // we got one ref for each file explicitly specified, just copy
	 }
       }
     }

     // check the titles argument (if default is used, it was already taken care of above)
     if ( TITLES.size() != 0 ) {
       if ( TITLES.size() != HEADS.size() ) {
	 cerr << "ERROR::wrong number of titles specified. If titles are set explicitly, " << endl
	      <<"their number must match the number of curves (files/file piles in the HEAD)" << endl;
	 exit(1);
       } else {
	 actual_titles = TITLES; // we got every title explicitly specified, just copy
       }
     }

     // computed GC curves will be first stored here
     vec<Curve> gc_curves(HEADS.isize());

     for ( int hi = 0; hi < HEADS.isize( ); hi++ )  {    
          // Get list of file heads in this group.

          HEADPILE = HEADS[hi];

          if (PRINTPOINTS) {
             cout << "Head #\tPercent\tGC Bias\n";   
          }

	  vecbasevector ref(actual_refs[hi]);
          ForceAssertGt( ref.size( ), 0u );

          // If reads have already been placed on the genome, load the placements.
          // Otherwise, compute from scratch.

          vec<placement_mark> places;
          for ( int u = 0; u < HEADPILE.isize( ); u++ ) {    
	      String HEAD = HEADPILE[u];
              if (USE_ALIGNS) {    
		   int nreads = MastervecFileObjectCount( HEAD + "." + READ_EXT );
		   if ( ! IsRegularFile(HEAD+".qltout") ) {
		     FatalErr(HEAD+".qltout file is missing");
		   }
		   MaxErrDiffAlignCollector aligns(0,MAX_ERRORS,nreads);
		   //aligns.resize(nreads);
		   LoadLookAligns( HEAD + ".qltout", aligns);
		   aligns.Consolidate();
		   for ( int id = 0 ; id < nreads ; id++ ) {
		     int n_aligns = aligns.Size(id);
		     if ( n_aligns == 0 ) continue;
		     const look_align & la = aligns.Aligns(id)[randomx( ) % n_aligns];
		     places.push( la.target_id, la.Fw1(), la.pos2() );   
		   }
		   cerr << places.size() << " alignments found" << endl;
	      }
	      else {
		if ( IsRegularFile( HEAD + ".perfect_marks_27" ) )  {
		   places.appendFromBinaryFile(HEAD + ".perfect_marks_27");
		}
	      } 
	  } // end for ( u < HEADPILE.isize ) - loop w/r to HEADPILE to merge alignments

	  // compute GC curve and associated GC bias metrics
	  double biasGCValL2 = ComputeGCCurve(ref, 
					     places.begin(),
					     places.end(),
					     W,
					     gc_curves[hi],
					     SAVE_MEMORY,
					     MAX_SLOP);
	  gc_curves[hi].SetTitle(actual_titles[hi]);

	  if (PRINTPOINTS) {
	      for ( unsigned int i = 0; i < gc_curves[hi].NPoints(); i++ ) {
	          cout << hi << "\t" << gc_curves[hi].x[i] << "\t" << gc_curves[hi].y[i] << "\n";
	      }
	  }

	  String GC_bias_str =  ToString( biasGCValL2, 2 ) ;
	  cout << "bias_GC=" << GC_bias_str << endl;
     }  // end of FOR (int hi=0 ; hi < HEADS.isize()...)


     if (DOGRAPH) {
       String x_label = "GC content of " + ToString(W) + "bp window";
       GenerateGCPlot(OUT, gc_curves, x_label.c_str());
     }     
     
}
