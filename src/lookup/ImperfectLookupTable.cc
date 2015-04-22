/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/** Find good unique alignments with no indels quickly.

\file ImperfectLookupTable.cc

Use ImperfectLookup() to find alignments. Save two files: 

- OUT_PREFIX.qltout contains the unique alignments in parseable 
and readable brief
format.
- OUT_PREFIX.minAlignErrors.txt: one line per read, with the number of errors 
in the best alignment for each read, or -1 if there are no 
alignments.

Parameters:
- SEQS: fasta or fastb file
- QUALS: fasta qual or qualb file. If absent, ignore quals
- K: kmer size, must match lookup table
- LOOKUP_TABLE, L: file name for lookup table
- OUT_PREFIX, O: prefix for output files, by default same as fasta file.
- OUT_SUFFIX: (default qltout) suffix for output alignment file.
- FWRC (default True): align in both direction, if False FW only.
- ERR_DIFF (default 2): if the second best alignment has ERR_DIFF or fewer
errors more than the best alignment, do not save any alignments for this
read (but do save the number of errors in the best alignment).
- CHUNK: read in the sequences in chunks of this size, to
limit memory usage.

*/

#include "MainTools.h"
#include "lookup/ImperfectLookup.h"
#include "lookup/LookupTable.h"
#include "FastaFileset.h"
#include "PrintAlignment.h"

class Printer {
  vecbasevector genome, genome_rc;
  String bar;

public:
  Printer(const String & GENOME);
  void Print(vec<look_align> & aligns, const vecbasevector & seqs);  
  template <class AlignCollector> void Print(const AlignCollector & aligns, const vecbasevector & seqs);  
  void Print(ostream & os, look_align la, const basevector & seq,
	     const qualvector * qual);
  void Print(ostream & os, look_align la, const basevector & seq);
  ~Printer() { cout << bar;   }
};


int main( int argc, char *argv[] )
{    
  RunTime();

  BeginCommandArguments;

  CommandArgument_String(SEQS);
  CommandArgument_String_OrDefault(QUALS, "");
  CommandArgument_UnsignedInt_OrDefault_Doc(K,12,"Unused; deduced from lookup table");
  CommandArgument_String_Abbr(LOOKUP_TABLE, L);
  CommandArgument_String_OrDefault_Doc(MODE,"UNIQ","Mode for collecting alignments: UNIQ - uses MAX_ERR and ERR_DIFF and"
				   " outputs only alignments that pass the filter; BEST - always outputs the best alignment(s) "
				   " with <= MAX_ERRS if they exist,  ERR_DIFF is ignored; BEST_NEXTBEST - outputs "
				   " best and next best alignment(s) for each read as long as the best number of errors "
				   " <= MAX_ERRS, ERR_DIFF is ignored; ALL - ouputs all (multiple) alignments with number "
				   " of errors <= (best) + ERR_DIFF for a read as long as best number of errors (best)<= MAX_ERRS. "
				   " Note: UNIQ is the fastest, ALL is the slowest"); 
  CommandArgument_String_Abbr_OrDefault_Doc(LOOKUP_TABLE_ALT, L_ALT, "",
       "Alternate lookup table: use if LOOKUP_TABLE does not exist.");
  CommandArgument_String_Abbr_OrDefault(OUT_PREFIX, O, "");
  CommandArgument_String_OrDefault(MIN_ERRORS_PREFIX, "");
  CommandArgument_String_OrDefault(OUT_SUFFIX, ".ilt.qltout");
  CommandArgument_Bool_OrDefault(FWRC, True);
  CommandArgument_UnsignedInt_OrDefault(ERR_DIFF,2);
  CommandArgument_Int_OrDefault_Doc(CHUNK,500000,"Align CHUNK reads at a time (saves memory)");
  CommandArgument_Bool_OrDefault(WRITE_MINALIGNS, True);
  CommandArgument_UnsignedInt_OrDefault(MAX_ERRS,4);
  CommandArgument_Bool_OrDefault_Doc(PRINT, False, "Print visual alignments to stdout");
  CommandArgument_Bool_OrDefault_Doc(WRITE, False, "Write visual alignments to the output file "
				     "in addition to the textual alignment information");
  CommandArgument_String_OrDefault(GENOME, "");
  //  CommandArgument_Bool_OrDefault_Doc(READ_BY
  CommandArgument_Int_OrDefault_Doc(START,0, "First read to process.");
  CommandArgument_Int_OrDefault_Doc(END,-1, "One past last read to align (if -1, end of file).");
  CommandArgument_Int_OrDefault_Doc(TRIM_LEFT,0,"Align left-trimmed reads. Alignment coordinates will be returned "
				    "in terms of full read (i.e. if we trim with TRIM_LEFT=5 and trimmed "
                                    "read aligns at position Z on the ref, alignment will be saved as having "
				    "target start position Z and query start position 5, not 0). "
			    "NOTE: alignments with up to TRIM_LEFT bases overhang beyond the contig end may be returned");
  CommandArgument_Int_OrDefault_Doc(TRIM_RIGHT,0,"Trim bases from the end of each read. Alignment coordinates will be "
				    "returned in terms of full read. Alignments with up to TRIM_RIGHT bases overhang "
				    "beyond the contig end may be returned");
  CommandArgument_Int_OrDefault_Doc(MAX_FREQ,0,"Ignore Kmers with frequencies above this threshold");
  CommandArgument_Bool_OrDefault_Doc(TIMER,False,"Report cumulative I/O and CPU time usage")

  EndCommandArguments;

  ForceAssertGe(TRIM_LEFT,0);
  ForceAssertGe(TRIM_RIGHT,0);

  if ( TRIM_LEFT < 0 || TRIM_RIGHT < 0 ) {
    cout << "Negative values for TRIM_LEFT or TRIM_RIGHT are not allowed" << endl;
    exit(1);
  }

  if ( LOOKUP_TABLE_ALT != "" && !IsRegularFile(LOOKUP_TABLE) )
  {    LOOKUP_TABLE = LOOKUP_TABLE_ALT;
       cout << "Using alternate lookup table.\n";    }

  ForceAssertEq( PRINT || WRITE , !GENOME.empty() );

  if (OUT_PREFIX.empty()) {
    OUT_PREFIX = SEQS.SafeBefore(".fa");
  }

  if ( TRIM_LEFT > 0 || TRIM_RIGHT > 0 ) {
    cout << "Read trimming requested; trimmed read ends will not be used in alignments/mismatch counts." << endl;
  }

  AlignCollectorBase *aligns ;

  if ( MODE == "UNIQ" ) {
      aligns = new  UniqueByErrDiffAlignCollector(ERR_DIFF,MAX_ERRS); 
  } else {
    if ( MODE == "BEST" ) {
      aligns = new BestAlignCollector(ERR_DIFF,MAX_ERRS); 
    } else {
      if ( MODE == "BEST_NEXTBEST") {
	aligns = new BestNextBestAlignCollector(ERR_DIFF,MAX_ERRS); 
      } else {
	if ( MODE == "ALL" ) {
	  aligns = new MaxErrDiffAlignCollector(ERR_DIFF,MAX_ERRS); 
	} else {
	  cout << "Unrecognized value specified for MODE" << endl;
	  exit(1);
	}
      }
    }
  }

  vec<TaskTimer *> timers;

  if ( TIMER ) {
    timers.push_back(new TaskTimer());
    timers.push_back(new TaskTimer());
  }    

  // if reads are given as fastb file, we can read them in chunks;
  // for fasta files chunks are not implemented yet...
  bool using_read_chunks = false; 
  bool using_qual_chunks = false; 
  ulonglong nreads = 0; // the number of reads we have read/are going to read:

  vecbasevector allseqs;
  if (SEQS.Contains("fasta") ) {
      FastFetchReads(allseqs, 0, SEQS);
      nreads = allseqs.size(); // we've just read them all
  } else {
    using_read_chunks = true; // if we got fastb we will read it by chunks
    nreads = MastervecFileObjectCount(SEQS); // have not read in the reads yet, but can find out the size!
  }
  if ( -1==END || (unsigned int)END > nreads ) {
      END = nreads;
  }
  PRINT3( nreads, START, END );
  
  vecqualvector * allquals = 0;
  if (!QUALS.empty()) {
    allquals = new vecqualvector;
    if (QUALS.Contains("qualb",-1)) {
      //      allquals->ReadAll(QUALS);
      using_qual_chunks = true;
    } else {
      // fetch all from quala - we can not read ranges yet:
      FastFetchQuals(*allquals, 0, QUALS);
    }
  }
  
  if (OUT_PREFIX.empty()) {
    OUT_PREFIX = SEQS.SafeBefore(".fa");
  }

  lookup_table look(LOOKUP_TABLE);
  
  cout << Date() << ": input read.\n";


  Ofstream(impltout, OUT_PREFIX + OUT_SUFFIX);
  PrintCommandPretty(impltout);

  Printer * printer = GENOME.empty() ? 0 : new Printer(GENOME);

  ostream * errStream = 0;
  String errFname = OUT_PREFIX + MIN_ERRORS_PREFIX + ".minAlignErrors.txt";
  if (WRITE_MINALIGNS) errStream = new ofstream(errFname.c_str());

  vecqualvector * quals = (0 == allquals ? 0 : new vecqualvector);

  longlong reads_processed = 0;
  longlong uniquealigned = 0;
  longlong aligned_well = 0;

  for (int start = START; start < END; start += CHUNK) {

    int end = min(start + CHUNK, END);
    PRINT2(start, end);
    vecbasevector seqs;
    if ( using_read_chunks ) {
      // we work with fastb file, can read by chunks:
      seqs.ReadRange(SEQS,start, end);
    } else {
      seqs.Append(allseqs, start, end);
    }
    if (quals) {
      quals->clear();
      if ( using_qual_chunks ) {
	// we got qualb, read range directly from disk!
	quals->ReadRange(QUALS,start,end);
      } else {
	// we can not read ranges directly from quala files yet; just copy for now:
	quals->Append(*allquals, start, end);
      }
    }
    
    if ( TRIM_LEFT > 0 || TRIM_RIGHT > 0 ) {
      // not very efficient: causes re-writing the vectors that were
      // already copied once from the original allseds and allquals vectors
      for ( size_t i = 0 ; i < seqs.size() ; i++ ) {
	seqs[i].SetToSubOf(seqs[i],TRIM_LEFT,seqs[i].isize() - TRIM_LEFT - TRIM_RIGHT);
	if ( quals ) {
	  (*quals)[i].SetToSubOf((*quals)[i],TRIM_LEFT,(*quals)[i].size() - TRIM_LEFT - TRIM_RIGHT);
	}
      }
    }

    aligns->clear();

    ImperfectLookup(look, seqs, 
		    *aligns,
		    FWRC ? FW_OR_RC : FW,
		    timers.size() == 2 ? &timers : 0  /* timers */, 
		    quals,
		    false /* use reverse*/,
		    MAX_FREQ);

    //adjust query_ids.
    aligns->AdjustQueryIds(start);

    if ( TRIM_LEFT > 0 || TRIM_RIGHT > 0 ) {
      // adjust alignment start positions:
      for ( unsigned int a = 0 ; a != aligns->size(); ++a ) {
	if ( ! aligns->UniquelyAligned(a) ) continue; // not uniquely aligned, don't bother
	if ( aligns->Align(a).IsQueryRC() ) aligns->MutableAlign(a).a.AddToStartOnQuery(TRIM_RIGHT);
	else aligns->MutableAlign(a).a.AddToStartOnQuery(TRIM_LEFT);
	aligns->MutableAlign(a).query_length += ( TRIM_LEFT+TRIM_RIGHT );
      }
    }
    
    unsigned int uniquely_aligned_in_chunk = UniquelyAlignedCount(*aligns);
    unsigned int aligned_in_chunk = AlignedCount(*aligns);
    cout << Date() << ": " << aligned_in_chunk << " reads aligned in chunk.\n";
    cout << Date() << ": " << uniquely_aligned_in_chunk << " reads aligned uniquely in chunk.\n";
    uniquealigned += uniquely_aligned_in_chunk;
    aligned_well += aligned_in_chunk;
    reads_processed += seqs.size();

    cout << "Running counts:" << endl;
    PRINT3( reads_processed, uniquealigned, aligned_well );

    if (WRITE_MINALIGNS) {
      for ( unsigned int a = 0 ; a < aligns->size() ; a++ ) {
	(*errStream) << aligns->MinErrors(a) << endl;
      }
    }

    aligns->Print(impltout,true); // print parseable and readable brief
    impltout.flush();

  } // END for ( int start = START; start < END ; start += CHUNK ) - the loop processing reads chunk by chunk

  delete printer;      
  delete errStream;
  delete allquals;
  delete quals;

  if ( timers.size() == 2 ) {
    cout << "Loading time: " << *(timers[0]) << endl;
    cout << "Aligning time: " << *(timers[1]) << endl;
    delete timers[0];
    delete timers[1];
  }
  
  const longlong N = END-START;
  cout << "Processed " << N << " reads.\n";
  cout << "Aligned " << VALUE_AND_RATIO(3, aligned_well, N) << "\n";
  cout << "Uniquely aligned " << VALUE_AND_RATIO(3, uniquealigned, N) << "\n";
  cout << Date( ) << ": done.\n";
  return 0;
}

Printer::Printer(const String & GENOME): genome(GENOME), genome_rc(genome) {
  for ( size_t i = 0; i < genome.size( ); i++ ) {
    genome_rc[i].ReverseComplement( );
  }
  bar = "----------------------------------------------------------"
    "--------------------------\n";
  cout << "\n";
}

void Printer::Print(vec<look_align> & aligns, 
		    const vecbasevector & seqs) {
  for ( int x = 0; x < aligns.isize( ); x++ ) {   
    cout << bar << "\n";
    look_align& la = aligns[x];
    int id = la.query_id;
    int tig = la.target_id;
    const basevector& g = ( la.rc1 ? genome_rc[tig] : genome[tig] );
    if (la.rc1) la.a.ReverseThis( seqs[id].size( ), g.size( ) );
    cout << id << ( la.rc1 ? "rc" : "fw" ) << " vs " << tig << ", " 
	 << la.mutations << " mismatches/" << la.indels 
	 << " indels (of " << la.query_length << "), from "
	 << la.pos1( ) << "-" << la.Pos1( ) << " to " << la.pos2( ) 
	 << "-" << la.Pos2( ) << " (of " << la.target_length << ")\n";
    PrintVisualAlignment( False, cout, seqs[id], g, la.a );    
  }
}

template <class AlignCollector>
void Printer::Print(const AlignCollector & aligns, 
		    const vecbasevector & seqs) {
  for ( unsigned int x = 0; x < aligns.size( ); x++ ) {   
    if ( !aligns.UniquelyAligned(x) ) continue;
    cout << bar << "\n";
    look_align la = aligns.Align(x);
    int id = la.query_id;
    int tig = la.target_id;
    const basevector& g = ( la.rc1 ? genome_rc[tig] : genome[tig] );
    if (la.rc1) la.a.ReverseThis( seqs[id].size( ), g.size( ) );
    cout << id << ( la.rc1 ? "rc" : "fw" ) << " vs " << tig << ", " 
	 << la.mutations << " mismatches/" << la.indels 
	 << " indels (of " << la.query_length << "), from "
	 << la.pos1( ) << "-" << la.Pos1( ) << " to " << la.pos2( ) 
	 << "-" << la.Pos2( ) << " (of " << la.target_length << ")\n";
    PrintVisualAlignment( False, cout, seqs[id], g, la.a );    
  }
}
  
void Printer::Print(ostream & os, look_align la, 
		    const basevector & seq,
		    const qualvector * qual) {
  int tig = la.target_id;
  //const basevector& g = genome[tig];
  const basevector& g = ( la.rc1 ? genome_rc[tig] : genome[tig] );
  //if (la.rc1) return;
  if (la.rc1) la.a.ReverseThis( seq.size( ), g.size( ) );
  PrintVisualAlignment( False, os, seq, g, la.a, *qual );  
}

void Printer::Print(ostream & os, look_align la, 
		    const basevector & seq) {
  int tig = la.target_id;
  //const basevector& g = genome[tig];
  const basevector& g = ( la.rc1 ? genome_rc[tig] : genome[tig] );
  //if (la.rc1) return;
  if (la.rc1) la.a.ReverseThis( seq.size( ), g.size( ) );
  PrintVisualAlignment( False, os, seq, g, la.a);  
}
  
