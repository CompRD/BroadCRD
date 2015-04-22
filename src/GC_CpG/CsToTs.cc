/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


/** Create 2 new FASTB files (fw and rc) with all Cs replaced by Ts
\file CsToTs.cc if a QUALB is given, create also the twin qualb files.

This is useful if you are doing a bisulfite experiment, which changes
all Cs to Ts, and then sequencing and trying to align.

Parameters:
   - FASTB: input fastb file
   - QUALB: if passed, do also qualb
   - OUT_PREFIX (O): prefix for new .fastb file, which will have two 
sequences for each sequence in the original file, the first being the
fw sequence with all Cs replaced by Ts, and the second the rc sequence
with all Cs similarly replaced by Ts. There will be two output files,
labeled OUT_PREFIX.TtoC.fw.fastb and OUT_PREFIX.TtoC.rc.fastb. If a
QUALB is passed, then there will be twin qualb files for the output
fastb. By default the prefix will be deduced from the input fastb
name.
   .
*/

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "FastaFileset.h"


int main( int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String_OrDefault_Doc(FASTB,"","Name or prefix of fastb file");
  CommandArgument_String_OrDefault_Doc(FASTA,"","Name or prefix of fasta file (mutually exclusive with FASTB)");
  CommandArgument_String_OrDefault_Doc(QUALB, "","Take this qualb file and create new qualb files synchronised "
                                                  "with the converted output sequences");
  CommandArgument_String_Abbr_OrDefault_Doc(OUT_PREFIX, O, "","If not specified, prefix of the input FASTA or FASTB file"
                                                              "will be used");
  EndCommandArguments;

  if ( FASTA.empty() && FASTB.empty() ) {
    cout << "Either FASTA or FASTB source must be specified" << endl;
    exit(1);
  }

  if ( !FASTA.empty() && ! FASTB.empty() ) { 
    cout << "FASTA and FASTB can not be simultaneously specified. Pick one!" << endl;
    exit(1);
  }

  if (OUT_PREFIX.empty()) {
    if ( !FASTB.empty() ) OUT_PREFIX = FASTB.SafeBefore(".fastb");
    else OUT_PREFIX = FASTA.SafeBefore(".fasta"); // we know by now that if FASTB is empty then FASTA is not
  }

  vecbasevector in;
  vecString inNames;

  if ( ! FASTB.empty() ) { // if we have to read FASTB, get the right name
    if ( ! IsRegularFile( FASTB ) ) {
      if ( IsRegularFile( FASTB + ".fastb" ) ) {
	FASTB+=".fastb";
      } else {
	cout << "Neither " << FASTB << " nor " << FASTB << ".fastb file exists. Aborting." << endl;
	exit(1);
      }
    }
  }

  if ( ! FASTA.empty() ) { // if we have to read FASTA, get the right name
    if ( ! IsRegularFile( FASTA ) ) {
      if ( IsRegularFile( FASTA + ".fasta" ) ) {
	FASTA+=".fasta";
      } else {
	cout << "Neither " << FASTA << " nor " << FASTA << ".fasta file exists. Aborting. " << endl;
	exit(1);
      }
    }
  }

  // ready to read; we know by now that input file exists

  if ( ! FASTB.empty() ) {
    in.ReadAll(FASTB); // read sequence data
    if ( IsRegularFile( FASTB + ".names" ) ) { // try to read names
      inNames.ReadAll( FASTB+".names" );
    } else {
      cout << "Warning: no sequence names (*.fastb.names) are found. Conversion will proceed." << endl;
    }
  } 

  if ( !FASTA.empty() ) { // if we have to read fasta
    FastFetchReads(in, &inNames, FASTA);
  }

  if ( ! inNames.empty() ) { // double check that numbers of sequences and sequence names match...
    ForceAssertEq( in.size(), inNames.size() );
  }

  cout << "Sequences imported successfully" << endl;
  flush(cout);

  // let us first process qualb (if requested) and forget about them:
  if ( QUALB != "" ) {
    if ( ! IsRegularFile(QUALB) ) {
    	if ( IsRegularFile(QUALB+".qualb") ) QUALB+=".qualb";
	else {
	    cout << "Neither " << QUALB << " nor " << QUALB << ".qualb files exist. Aborting. " << endl;
	    exit(1);
	}
    }
    vecqvec qin;
    qin.ReadAll( QUALB );
    qin.WriteAll(OUT_PREFIX+".fw.CtoT.qualb"); // fw qual is the same as the original one!

    vecqvec::iterator end(qin.end());
    for ( vecqvec::iterator itr(qin.begin()); itr != end; ++itr )
        itr->ReverseMe();
    qin.WriteAll(OUT_PREFIX+".rc.CtoT.qualb"); // save reverse complement qualities
  }
  // done with quality files, let's now revert sequences:


  if ( ! FASTA.empty() ) {
     // if we work with fasta, we can convert and write on the fly:
     // (this will save us some memory)
     Ofstream( fwfasta , OUT_PREFIX+".fw.CtoT.fasta");
     Ofstream( rcfasta , OUT_PREFIX+".rc.CtoT.fasta");

     basevector contig;
     for ( size_t i = 0 ; i < in.size() ; i++ ) {
         contig = in[i];
	 for ( unsigned int j = 0 ; j < contig.size() ; j++ ) {
	 	if ( contig[j] == BASE_C ) contig.Set(j,BASE_T); // convert C to T
         }
	 contig.Print(fwfasta, inNames[i]+" [C-->T]");

	 // now convert reverse complement:
	 contig.ReverseComplement(in[i]);
	 for ( unsigned int j = 0 ; j < contig.size() ; j++ ) {
	 	if ( contig[j] == BASE_C ) contig.Set(j,BASE_T); // convert C to T
         }
	 contig.Print(rcfasta, inNames[i]+" [rc] [C-->T]");
     }
     // WE ARE DONE NOW! 
     fwfasta.close();
     rcfasta.close();
  }

  if ( ! FASTB.empty() ) {
     // uh-oh, we are working with fastb, got to convert whole vecbasevec, then write it all:

     vecbasevector fw, rc;
     vecString fwNames, rcNames;

     fw.reserve(in.size());
     rc.reserve(in.size());

     basevector contig;
     for (size_t i = 0; i != in.size(); ++i) {

         contig = in[i]; // retrieve sequence (contig)

         for (unsigned int j = 0; j != contig.size() ; ++j ) {
             if (contig[j] == BASE_C ) contig.Set(j, BASE_T) ; // change C to T
         }

         fw.push_back(contig); // save converted forward seq

         // if we have names, create name for the fw seq:
         if ( inNames.size() != 0 ) fwNames.push_back( inNames[i] + " [C-->T]" );
 
         // now convert reverse complement of the sequence:
         contig.ReverseComplement(in[i]);

         for (unsigned int j = 0; j != contig.size() ; ++j ) {
             if (contig[j] == BASE_C ) contig.Set(j, BASE_T) ; // change C to T
         }

         rc.push_back(contig); // save converted revese complement seq

         // if we have names, create name for the rc seq:
         if ( inNames.size() != 0 ) rcNames.push_back(inNames[i] + " [rc] [C-->T]" );
     }
     fw.WriteAll(OUT_PREFIX+".fw.CtoT.fastb");
     rc.WriteAll(OUT_PREFIX+".rc.CtoT.fastb");
     if ( inNames.size() != 0 ) {
        fwNames.WriteAll(OUT_PREFIX+".fw.CtoT.fastb.names");
        rcNames.WriteAll(OUT_PREFIX+".rc.CtoT.fastb.names");
     }
     // AND WE ARE DONE WITH FASTB!!
  }

  return 0;
}

