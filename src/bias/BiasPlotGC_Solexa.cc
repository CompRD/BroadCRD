/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

const char *DOC =

"For each read perfect to N bases, associate the GC content "
"of the W-base window on the reference, starting at the beginning of the read.  "
"For each possible GC content, find the number of windows on the reference "
"having that GC content.  Then to each such GC content, associate "
"ratio (# of reads) / (# of windows on reference).  Normalize so that GC = 50% "
"has height 1."

;

#include "bias/BiasGC.h"
#include "BasevectorTools.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "feudal/BinaryStream.h"
#include "graphics/BasicGraphics.h"
#include "lookup/AlignCollector.h"
#include "lookup/PerfectCount.h"
#include "math/Functions.h"
#include "random/Random.h"
#include "solexa/SolexaTools.h"
#include "solexa/SolexaMetrics.h"
#include "math/SparseArray.h"
#include "solexa/SolexaPipeline.h"
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
          "flowcell.lane, or date_flowcell.lane, or list of such, separated by "
          "commas to keep separate or plus signs to merge; expressions like "
          "4321.{1,2,3} or 4321.{1+2} are formally expanded; at most 8 lanes in "
          "total are allowed; do not place entire list in brackets { }" );
     CommandArgument_String_OrDefault_Doc(REFHEAD,"","if specified, HEAD.REFHEAD.reference.* will "
					  "be used as reference file(s) for each lane found in HEAD");
     CommandArgument_Int_OrDefault_Doc(N, 27, "length to truncate reads to");
     CommandArgument_Int_OrDefault_Doc(W, 50, "window size");
     CommandArgument_String_OrDefault_Doc(OUT, "", "output file, may end with .png");
     CommandArgument_Double_OrDefault_Doc(MAX_SLOP, 0.05, "spread of 95% "
          "confidence interval must be less than this to plot point" );
     CommandArgument_Double_OrDefault(TITLE_FONTSIZE, 15);
     CommandArgument_Bool_OrDefault_Doc(FILTER, False,
          "If True, filter out reads whose signal looks bad." );
     CommandArgument_Bool_OrDefault_Doc(DOGRAPH, True,
          "If True, output graph." );
     CommandArgument_Bool_OrDefault_Doc(ADDMETRIC, False,
          "If True, add gc window flatness to the metrics table." );
     CommandArgument_Bool_OrDefault_Doc(PRINTPOINTS, False,
          "If True, print to stdout the plot points." );
     CommandArgument_Bool_OrDefault_Doc(USE_ALIGNS, False,
          "If True, instead of using reads perfect to N bases, use aligned reads.");
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_String_OrDefault(PIPELINE, SOLEXAPIPELINE);
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

     // Parse HEAD into HEADS.

     vec< vec<String> > HEADS;
     String HEAD2, CURHEAD;
     for ( int i = 0; i < HEAD.isize( ); i++ ) 
     {    
          // Capture up to first left bracket.

          if ( HEAD[i] != '{' ) 
          {    HEAD2 += HEAD[i];
	       continue;    }

          // Find matching right bracket.

	 int j;
	 for ( j = i + 1; j < HEAD.isize( ); j++ )
              if ( HEAD[j] == '}' ) break;
	 if ( i < 5 || HEAD[i-1] != '.' || j == HEAD.isize( ) ) BRACKFAIL
	 String brack = HEAD.substr( i, j - i + 1 );

         // Find the lanes defined by the bracketed part.

	 vec<int> lanes;
	 Bool plus = False;
	 if ( brack.Contains( "+" ) ) 
         {    if ( brack.Contains( "," ) ) BRACKFAIL
	      plus = True;
	      brack.GlobalReplaceBy( "+", "," );    }
	 ParseIntSet( brack, lanes );
	 if ( lanes.empty( ) ) BRACKFAIL

         // Walk backward until we hit the beginning or a comma.  This should give
         // us the flowcell name plus the dot that follows, including the date that 
         // comes before, if present.

         int fcdotlen = 0;
         for ( int k = HEAD2.isize( ) - 1; k >= 0; k-- )
         {    if ( HEAD2[k] == ',' ) break;
              ++fcdotlen;    }

	 String fcdot = HEAD2.substr( HEAD2.isize( ) - fcdotlen, fcdotlen );
	 HEAD2.resize( HEAD2.isize( ) - fcdotlen );
	 for ( int k = 0; k < lanes.isize( ); k++ ) 
         {    if ( k > 0 ) HEAD2 += ( plus ? '+' : ',' );
	      HEAD2 += fcdot + ToString( lanes[k] );    }
	 i = j;    }
     vec<String> HEADPILE;
     for ( int i = 0; i < HEAD2.isize( ); i++ ) 
     {    if ( HEAD2[i] == ',' ) 
          {    HEADPILE.push_back(CURHEAD), HEADS.push_back(HEADPILE);
	       HEADPILE.clear( ), CURHEAD.clear( );    } 
          else if ( HEAD2[i] == '+' ) 
          {    HEADPILE.push_back(CURHEAD), CURHEAD.clear( );    }  
          else CURHEAD += HEAD2[i];    }
     HEADPILE.push_back(CURHEAD), HEADS.push_back(HEADPILE);
     for ( int i = 0; i < HEADS.isize( ); i++ )
          Sort( HEADS[i] );

     // Create compressed version of each HEADS[i], for title.

     vec<String> TITLES;
     for ( int hi = 0; hi < HEADS.isize( ); hi++ ) {    
         vec<String> H = HEADS[hi];
	 String TITLE;
	 for ( int i = 0; i < H.isize( ); i++ )  {    
	     int j;
	     for ( j = i + 1; j < H.isize( ); j++ )
	        if ( H[j].Before( "." ) != H[i].Before( "." ) ) break;
	     if ( i > 0 ) TITLE += "+";
	     if ( j - i == 1 ) TITLE += H[i];
	     else {    
	         vec<int> lanes;
		 for ( int k = i; k < j; k++ )
                         lanes.push_back( H[k].After( "." ).Int( ) );
		 ForceAssert( lanes.UniqueOrdered( ) );
		 TITLE += H[i].Before( "." ) + ".{";
		 if ( lanes.isize( ) == lanes.back( ) - lanes.front( ) + 1
		      && lanes.size( ) > 3 ) {    
		     TITLE += ToString( lanes.front( ) ) + "+...+"
		            + ToString( lanes.back( ) );    
		 } else {    
		     for ( int k = i; k < j; k++ ) {    
		         if ( k > i ) TITLE += "+";
			 TITLE += ToString( lanes[k-i] );    
		     }    
		 }
		 TITLE += "}";    
	     }
	     i = j - 1;    
	 }
	 TITLES.push_back(TITLE);    
     }
     // Set up to plot points.  We define a fixed list of eight colors.

     vec<Curve> gc_curves(HEADS.isize());

     for ( int hi = 0; hi < HEADS.isize( ); hi++ )  {    
          // Get list of lanes in this group.

          vec<String> HEAD = HEADS[hi];

          // Create HEADXS = full paths to HEADS.

          vec<String> HEADXS;
          String HEADX = "{";
          for ( int i = 0; i < HEAD.isize( ); i++ )  {    
	      if ( i > 0 ) HEADX += ",";
	      HEADX += HEAD[i];    
	  }
          HEADX += "}";
          ExpandHead( HEADX, HEADXS, 0, PIPELINE );

          // Print out stdout header
          if (PRINTPOINTS) {
             cout << "Head #\tPercent\tGC Bias\n";   
          }

          // Load reference.
          RequireSameRef(HEADXS);
          vecbasevector ref( HEADXS[0] + ( REFHEAD=="" ? "" : ("."+REFHEAD ) ) + ".reference.fastb" );
          ForceAssertGt( ref.size( ), 0u );

          // If reads have already been placed on the genome, load the placements.
          // Otherwise, compute from scratch.

          vec<placement_mark> places;
          for ( int u = 0; u < HEADXS.isize( ); u++ ) {    
	      String HEAD = HEADXS[u];
              if (USE_ALIGNS) {    
		   int nreads = MastervecFileObjectCount( HEAD + ".fastb" );
		   if ( ! IsRegularFile(HEAD+".qltout") ) {
		     FatalErr(HEAD+".qltout file is missing");
		   }
		   MaxErrDiffAlignCollector aligns(0,MAX_ERRORS,nreads);
		   //aligns.resize(nreads);
		   LoadLookAligns( HEAD + ".qltout", aligns );
		   aligns.Consolidate();
		   for ( int id = 0; id < nreads; id++ ) {  
		       int n_aligns = aligns.Size(id);
		       if ( n_aligns == 0 ) continue;
		       const look_align& la 
			    = aligns.Aligns(id)[ randomx() % n_aligns ];
		       places.push( la.target_id, la.Fw1(), la.pos2() );    
		   }
	      }
	      else if ( !FILTER && N == 27 && IsRegularFile( HEAD + ".perfect_marks_27" ) )  {
	          places.appendFromBinaryFile(HEAD + ".perfect_marks_27");
	      } else {    
		   vecbasevector bases( HEAD + ".fastb" );
                    for ( size_t i = 0; i < bases.size( ); i++ )
                         bases[i].resize(N);
                    if (FILTER)
                    {    VEC_FROM_PARAMS_FILE(HEAD+".paramsByRead", minq10);
                         VEC_FROM_PARAMS_FILE(HEAD+".paramsByRead", lowestHeight10);
                         for ( size_t id = 0; id < bases.size( ); id++ )
                         {    if ( lowestHeight10[id] < 1000 || minq10[id]  < 2.0 )
                              bases[id].resize(0);    }    }
                    vec<placement_mark> placesi;
                    // Conditional added else interfers with points output
                    if (!PRINTPOINTS) {
                      cout << "computing perfect aligns for " << HEAD << endl;
                    }
                    PerfectPick( bases, HEAD + ( REFHEAD=="" ? "" : ("."+REFHEAD ) ) + 
				 ".reference.lookup", 
				 FW_OR_RC, placesi );
                    places.append(placesi);    
	      }    
	  } // end for ( u < HEADXS.isize ) - loop w/r to HEADXS to merge alignments

          // Compute read start points.

	  
	  double biasGCValL2 = ComputeGCCurve(ref, 
					     places.begin(),
					     places.end(),
					     W,
					     gc_curves[hi],
					     SAVE_MEMORY,
					     MAX_SLOP);
	  gc_curves[hi].SetTitle(TITLES[hi]);

	  if (PRINTPOINTS) {
	      for ( unsigned int i = 0; i < gc_curves[hi].NPoints(); i++ ) {
	          cout << hi << "\t" << gc_curves[hi].x[i] << "\t" << gc_curves[hi].y[i] << "\n";
	      }
	  }

	  String GC_bias_str =  ToString( biasGCValL2, 2 ) ;
	  cout << "bias_GC=" << GC_bias_str << endl;
	  if (ADDMETRIC) {
              for ( int u = 0; u < HEADXS.isize( ); u++ ) {
                  solexa_metric_db db(HEADXS[u] + ".metrics");
                  db.SetValue("bias_GC", GC_bias_str);
		  if (WRITE) db.WriteMetrics( HEADXS[u] + ".metrics" );
		  else db.WriteMetrics(cout);
              }
	  }
     }  // end of FOR (int hi=0 ; hi < HEADS.isize()...)


     if (DOGRAPH) {
       String x_label = "GC content of " + ToString(W) + "bp window";
       GenerateGCPlot(OUT, gc_curves, x_label.c_str(), 0, TITLE_FONTSIZE);
     }     
     
}
