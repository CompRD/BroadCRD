/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: UnibaseCopyNumberIdeal

   Determines the true copy number of each unibase by aligning
   the unibase to the reference.  Thus, it tests the limits of UnibaseCopyNumber.

   @file
 */

#include <algorithm>
#include <list>
#include <strstream>
#include "Basevector.h"
#include "BinsVec.h"
#include "Bitvector.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "PredictionStats.h"
#include "TokenizeString.h"
#include "lookup/ImperfectLookup.h"
#include "lookup/LookAlign.h"
#include "lookup/PerfectLookup.h"
#include "math/Functions.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/PdfEntry.h"
#include "paths/simulation/GenomePlacement.h"
#include "graph/Digraph.h"
#include "paths/Unipath.h"
#include "paths/UnibaseUtils.h"
#include "paths/reporting/ReftigUtils.h"


Bool cmp_ali( const look_align& la1, const look_align& la2 ){    
  if ( la1.query_id < la2.query_id ) return True;
  if ( la1.query_id > la2.query_id ) return False;
  if ( la1.mutations + la1.indels < la2.mutations + la2.indels ) return True;
  if ( la1.mutations + la1.indels > la2.mutations + la2.indels ) return False;
  if ( la1.a.pos2( ) < la2.a.pos2( ) ) return True;
  if ( la1.a.pos2( ) > la2.a.pos2( ) ) return False;
  if ( la1.target_id < la2.target_id ) return True;
  if ( la1.target_id > la2.target_id ) return False;
  return False;    
}

int main( int argc, char** argv ) {
  RunTime();

  BeginCommandArguments;
  CommandDoc("Compute unibase copy number based on alignment to a reference.");

  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int(K);
  CommandArgument_String(READS);
  CommandArgument_String(UNIBASES);
  CommandArgument_Bool_OrDefault_Doc(PERFECT, True,
	 "If set, use only perfect alignments in computing copy numbers");
  CommandArgument_Bool_OrDefault_Doc(STRICT, False,
	 "If set (and PERFECT is set False), count multiple imperfect alignment as equivalent only if they have exactly the same errors. Otherwise they are equivalent if they have the same number of errors");
  
  CommandArgument_String_OrDefault_Doc( UNIPATHS, "unipaths", "Only affects output file name" );
  CommandArgument_String_OrDefault(GENOME,"genome");
  CommandArgument_String(COUNT_SUFFIX);

  EndCommandArguments;


  dirname_t predata = PRE + "/" + DATA;
  dirname_t run_dir = predata + "/" + RUN;
  
  String unibasesFile =  run_dir + "/" + READS + "." + UNIBASES + ".k" + ToString(K);
  if (!IsRegularFile(unibasesFile) && IsRegularFile(unibasesFile + ".fastb") )
    unibasesFile += ".fastb";
  vecbasevector unibases(unibasesFile);



  // Align each unibase to the reference to determine its copy number.
  
  vec<look_align> aligns;
  cout << " Aligning " << unibases.size() << " unibases to the reference to determine true copy #" << endl;
  String alignmentsFile = run_dir + "/" + READS + "." + UNIBASES + ".k" + ToString(K) + ".qltout";
  
  if ( PERFECT ){
    // use only perfect alignments
    PerfectLookup( 12 /* seed alignments to the reference on kmers of this size */, 
		   unibases, predata + "/" + GENOME + ".lookup", aligns, FW_OR_RC );
    Ofstream(pltout, alignmentsFile);
    for (int i=0; i != aligns.isize(); ++i) {
      aligns[i].PrintParseable(pltout);
      aligns[i].PrintReadableBrief(pltout);
    }
  }else{
    // use perfect and imperfect alignments
    GetAlignsFast( 12, /* seed alignments to the reference on kmers of this size */
		   unibasesFile,
		   predata + "/" + GENOME + ".lookup",
		   alignmentsFile,
		   aligns, false, run_dir + "/tmp" );
  }
  
  sort( aligns.begin(), aligns.end(), cmp_ali );
  

  // Find best alignments

  vecbasevector ref( predata + "/" + GENOME + ".fastb" );
  cout << " Looking at " << aligns.isize() << " total alignments." << endl;
  vec<int> bestMismatches( unibases.size(), -1 );
  vec< vec<align_id_t> > bestAlignsIndex( unibases.size() );
  for ( align_id_t i = 0; i < aligns.isize( ); i++ ) {
    if ( !( i % 10000 ) )
      PRINT( i );
    const look_align& la = aligns[i];
    int uid = la.query_id;
    int mismatches = la.mutations + la.indels;
    if ( la.FullLength() ){
      if ( bestMismatches[uid] == -1 )
	bestMismatches[uid] = mismatches;
      
      if ( mismatches == bestMismatches[uid] )
	bestAlignsIndex[uid].push_back( i );
    }
  }

  
  ForceAssertEq( bestAlignsIndex.size(), unibases.size() );
  vec<copy_num_t> trueCopyNumber( unibases.size( ), 0 );
  if ( ! STRICT ){
     for ( size_t uid = 0; uid < bestAlignsIndex.size(); uid++ )
       trueCopyNumber[uid] = bestAlignsIndex[uid].size();
  }else{
    for ( size_t uid = 0; uid < bestAlignsIndex.size(); uid++ ){
      if ( bestAlignsIndex[uid].size() <= 1 )
	trueCopyNumber[uid] = bestAlignsIndex[uid].size();
      else{
	const look_align& la0 = aligns[ bestAlignsIndex[uid][0] ];
	basevector tseq0( ref[ la0.target_id ], la0.pos2(), la0.Pos2() - la0.pos2() );
	if ( la0.Rc1() )
	  tseq0.ReverseComplement();
	//cout << "\n0: ";tseq0.Print(cout);
	trueCopyNumber[uid] = 1;
	for ( size_t j = 1; j < bestAlignsIndex[uid].size(); j++ ){
	  const look_align& la = aligns[ bestAlignsIndex[uid][j] ];
	  basevector tseq( ref[ la.target_id ], la.pos2(), la.Pos2() - la.pos2() );
	  if ( la.Rc1() )
	    tseq.ReverseComplement();
	  //cout << j << ": "; tseq.Print( cout );
	  if ( tseq != tseq0 ){
	    trueCopyNumber[uid] = 0;
	    break;
	  }else
	    trueCopyNumber[uid]++;
	}      
	//PRINT2( uid, trueCopyNumber[uid] ); cout << "<<<<<<<\n";
      }
    }
  }
  

  PRINT( unibases.size() );
  //trueCopyNumber.Print(cout);
  int maxCn = Max( trueCopyNumber );
  vec<int> cnHisto( maxCn +1, 0 );
  for ( size_t uid = 0; uid < trueCopyNumber.size(); uid++ )
    cnHisto.at( trueCopyNumber[uid] )++;
  
  cout << "\ncn\tcount\n";
  for ( size_t i = 0; i < cnHisto.size(); i++ )
    if ( cnHisto[i] > 0 )
      cout << i << "\t" << cnHisto[i] << endl;
  cout << endl;

  VecPdfEntryVec cn_pdfs( unibases.size() );

  for ( size_t u = 0; u < unibases.size(); u++ ) {
    cn_pdfs[ u ].push_back( pdf_entry( trueCopyNumber[ u ], 1.0 /* this copy # with probability 1 */ ) );
  }

  cn_pdfs.WriteAll( (run_dir + "/" + READS + "." + UNIPATHS + "." + COUNT_SUFFIX + ".k" + ToString(K)).c_str() );
  
}



