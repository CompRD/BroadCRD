///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Reverse complements each basevector in a fastb file. Can also reverse the "
  "quality scores in associated qualb and/or qlt files. Has an option to just "
  "reverse the base sequence rather than compute the reverse complement. "
  "Supports files containing paired reads and allows the reversal of either "
  "read of a pair, or both. Pairing information may be supplied in a pairs/pairto "
  "file or implied by the read order (interleaved where read A is followed by "
  "read B)";

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "PairsManager.h"
#include "lookup/LookAlign.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(READS_IN,
    "Fastb file containing base sequences to reverse complement");
  CommandArgument_String_OrDefault_Doc(READS_OUT, "",
    "Output fastb file. Default is READ_IN.rc.fastb");
  CommandArgument_Bool_OrDefault_Doc(QUALS, False,
    "Reverse associated qualb files");
  CommandArgument_Bool_OrDefault_Doc(QLTS, False,
    "Update associated qltout files");
  CommandArgument_Bool_OrDefault_Doc(RC, True,
    "If false, then just reverse sequences, don't reverse complement");
  CommandArgument_Bool_OrDefault_Doc(PAIRED, False,
    "Reads are paired, either interleaved or with an associated pairs/pairto file");
  CommandArgument_String_OrDefault_Valid_Doc(REVERSE_PAIRED, "both",
    "{both,first,second}",
    "For paired reads, which read of the pair to reverse");
  CommandArgument_String_OrDefault_Valid_Doc(PAIRING_INFO, "file",
    "{file,interleaved}",
    "Use pairs/pairto file or assume paired reads are interleaved");

  CommandArgument_UnsignedInt_OrDefault_Doc( TRIMMED_SIZE, 0, 
    "Final size after optional trimming" );
  CommandArgument_String_OrDefault_Valid_Doc(TRIM_FIRST_SIDE, "begin",
    "{begin,end}",
    "for paired reads, which side (AFTER REVERSAL) to trim the first read of the pair");
  CommandArgument_String_OrDefault_Valid_Doc(TRIM_SECOND_SIDE, "begin",
    "{begin,end}",
    "for paired reads, which side (AFTER REVERSAL) to trim the first read of the pair");

  EndCommandArguments;

  if (TRIMMED_SIZE > 0 && PAIRED && PAIRING_INFO == "file") {
    cout << "trimming not implemented for non-interleaved paired reads. No trimming will be performed." 
         << endl;
    TRIMMED_SIZE = 0;
  }
  if ( TRIMMED_SIZE > 0 && QLTS ) {
    cout << "updating qlts with trimming not implemented: no trimming will be performed." 
         << endl;
    TRIMMED_SIZE = 0;
  }

  // Sanity check args

  if (QLTS && (REVERSE_PAIRED != "both" || TRIMMED_SIZE !=0 ) )
      InputErr("Trimming or reversing only one read of a pair cannot be combined with QLTS=True");

  // Strip .fastb from filenames

  READS_IN = READS_IN.SafeBefore(".fastb");
  READS_OUT = READS_OUT.SafeBefore(".fastb");

  if (READS_OUT == "")
    READS_OUT = READS_IN + (RC ? ".rc" : ".reversed");

  // Load reads

  cout << Date() << " Loading reads..." << endl;
  vecbasevector reads(READS_IN + ".fastb");
  size_t nreads = reads.size();

  // Load quals

  vecqualvector quals;
  if (QUALS && IsRegularFile( READS_IN + ".qualb") ) {
    cout << Date() << " Loading quals..." << endl;
    quals.ReadAll(READS_IN + ".qualb");
    ForceAssertEq(reads.size(), quals.size());
  } else if (QUALS) {
    cout << "Warning: Unable to find quals. (qualb)" << endl;
    QUALS = false;
  }

  cout << Date() << " Found " << nreads << " reads" << endl;

  // Reverse Reads (and quals)

  if (!PAIRED || REVERSE_PAIRED=="both") {

    // Unpaired or both paired reads

    cout << Date() << " Reversing reads" << endl;
    if (RC)
      ReverseComplement(reads);
    else
      for ( size_t i = 0; i < nreads; ++i )
	reads[i].Reverse();

    if (QUALS) {
      cout << Date() << " Reversing quals" << endl;
      for ( size_t i = 0; i < nreads; ++i )
	quals[i].ReverseMe();
    }

  } 
  else if (PAIRING_INFO == "interleaved") {

    // reverse only one of each pair of interleaved reads 

    size_t start;
    if (REVERSE_PAIRED == "first") {
      start = 0;
      cout << Date() << " Reversing only first pair of interleaved paired reads (and quals)" << endl;
    }
    else {
      start = 1;
      cout << Date() << " Reversing only second pair of interleaved paired reads (and quals)" << endl;
    }

    for (size_t i = start; i < nreads; i+=2) {
      if (RC)
	reads[i].ReverseComplement();
      else
	reads[i].Reverse();
      if (QUALS) 
	quals[i].ReverseMe();
    }
  } 
  else {

    // reverse only one of each pair of reads as defined by a pairs/pairto file

    cout << Date() << " Loading pairing information..." << endl;
    PairsManager pairs;
    if (IsRegularFile(READS_IN + ".pairs"))
      pairs.Read(READS_IN + ".pairs");
    else
      pairs.ReadFromPairtoFile(READS_IN + ".pairto", nreads);
    size_t npairs = pairs.nPairs();

    cout << Date() << " Reversing paired reads (and quals)" << endl;
    for (size_t pairId = 0; pairId < npairs; pairId++ ) {
      size_t readId = (REVERSE_PAIRED == "first" ? 
                       pairs.ID1(pairId) : pairs.ID2(pairId));
      if (RC)
	reads[readId].ReverseComplement();
      else
	reads[readId].Reverse();
      if (QUALS) 
	quals[readId].ReverseMe();
    }
    
  }
  
  // --------------- trimming and saving --------------

  if (TRIMMED_SIZE > 0) {

    IncrementalWriter<bvec> bases_writer((READS_OUT + ".fastb").c_str());
    IncrementalWriter<qvec> quals_writer((READS_OUT + ".qualb").c_str());
    qvec qual_new(TRIMMED_SIZE);

    for (size_t i = 0; i < nreads; i++) {
      const bvec & read_orig = reads[i];

      ForceAssertGe(read_orig.size(), TRIMMED_SIZE);


      size_t ibase0;
      if (!PAIRED) {
        ibase0 = (TRIM_FIRST_SIDE == "end") ? 0 : read_orig.size() - TRIMMED_SIZE;
        if (i == 0)  cout << Date() << " Trimming reads - keeping bases " 
				   << ibase0 << " to " << ibase0 + TRIMMED_SIZE - 1 << endl;
      }
      else { // paired reads
        if ((i % 2) == 0)  // even read 
          ibase0 = (TRIM_FIRST_SIDE == "end") ? 0 : read_orig.size() - TRIMMED_SIZE;
        else  // odd read
          ibase0 = (TRIM_SECOND_SIDE == "end") ? 0 : read_orig.size() - TRIMMED_SIZE;

        if (i == 0 || i == 1) cout << Date() << " Trimming " << (i == 0 ? "first" : "second") 
				   << " read of pair - keeping bases " 
				   << ibase0 << " to " << ibase0 + TRIMMED_SIZE - 1 << endl;
      }

      bases_writer.add(bvec(read_orig, ibase0, TRIMMED_SIZE));
      if (QUALS) {
        //quals_writer.add(qvec(quals[i], ibase0, TRIMMED_SIZE));
        CopyQuals(quals[i], ibase0, qual_new, 0, TRIMMED_SIZE);
        quals_writer.add(qual_new);
      }
    }
    bases_writer.close();
    quals_writer.close();


  }
  else { // --------------- saving without trimming -----------------
 
    // Write reversed reads (and quals)
    
    cout << Date() << " Writing files." << endl;
    reads.WriteAll(READS_OUT + ".fastb");
    if (QUALS)
      quals.WriteAll(READS_OUT + ".qualb");
    
  }

  // Update qlts.
  if ( QLTS ) {
    String qlt_in = READS_IN + ".qltout";
    String qlt_out = READS_OUT + ".qltout";

    if ( ! IsRegularFile( qlt_in ) ) {
      cout << "Warning: Unable to find aligns. (qltout)" << endl;
    }
    else {
      cout << Date( ) << " Loading and updating aligns" << endl;
      vec<look_align_plus> qlts;
      LoadLookAlignPlus( qlt_in, qlts );
      
      ofstream qout( qlt_out.c_str( ) );
      for (size_t ii=0; ii<qlts.size( ); ii++) {
	qlts[ii].rc1 = ! qlts[ii].rc1;
	qlts[ii].WriteParseable( qout );
      }
      qout.close( );
    }
  }

  cout << Date() << " Finished" << endl;

}
