/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "FastaFileset.h"

/**
 * SimpleACGTContent
 *
 * Determine frequency of A, C, G, and T from a vector of bases. Optionally,
 * it can determine such frequency by looking only at bases with quality
 * bigger than a given threshold. If OUTFILE is not assigned, output will
 * be sent to cout.
 *
 * BASES: full path name for fastb file
 * SELECT: use these ids only (full path name, file loaded with macro READ)
 * OUTFILE: optional output file (full path name)
 * QUALS: optional full path name for qualb file
 * MIN_QUAL: used only if QUALS is given (use base iff qual >= MIN_QUAL)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String_Doc( BASES,"Input: collection of sequences/contigs/reads (fasta or fastb file)" );
  CommandArgument_String_OrDefault_Doc( SELECT, "","If specified, must be a full qualified path to the file "
					"cotaining sequence ids; only bases from sequences with specified ids in the "
					"BASES file will be counted");
  CommandArgument_String_OrDefault_Doc( OUTFILE, "","Output file; if not specified results are printed to stdout" );
  CommandArgument_String_OrDefault_Doc( QUALS, "","Base quality file" );
  CommandArgument_Bool_OrDefault_Doc(PER_SEQUENCE, False,"If True, print table <seq_id> <#A> <#C> <#G> <#T> "
				     "(one line per sequence/contig in BASES file); otherwise print cumulative counts across all sequences");
  CommandArgument_UnsignedInt_OrDefault_Doc(TRIM_LEFT,0,"Ignore TRIM_LEFT leftmost bases in each sequence");
  CommandArgument_UnsignedInt_OrDefault_Doc(TRIM_RIGHT,0,"Ignore TRIM_RIGHT rightmost bases in each sequence");
  CommandArgument_UnsignedInt_OrDefault_Doc( MIN_QUAL, 20,"If specified (and QUALS is given), then count only bases with "
					       "qualities better than MIN_QUAL" );
  EndCommandArguments;

  // Output stream.
  ostream *outp = 0;
  if ( OUTFILE == "" )
    outp = (ostream*) &cout;
  else
    outp = new ofstream( OUTFILE.c_str( ) );
  ostream &out = *outp;
  
  // Load.
  //  cout << Date( ) << ": loading bases" << endl;
  vecbasevector bases;
  LoadReads(bases,BASES);

  vecqualvector quals;
  if ( QUALS != "" ) {
    cout << Date( ) << ": loading quals" << endl;
    quals.ReadAll( QUALS );
  }

  // Select reads.
  vec<Bool> select;
  if ( SELECT == "" )
    select.resize( bases.size( ), True );
  else {
    select.resize( bases.size( ), False );
    READ( SELECT, vec<int>, ids );
    for (int ii=0; ii<(int)ids.size( ); ii++)
      select[ ids[ii] ] = True;
  }

  // Total tally.
  longlong tot_valid = 0;
  longlong tot = 0;
  longlong totA = 0;
  longlong totC = 0;
  longlong totG = 0;
  longlong totT = 0;

  // Parse bases.
  for (int ii=0; ii<(int)bases.size( ); ii++) {
    if ( ! select[ii] )
      continue;
    const basevector &cg_bases = bases[ii];
    tot += cg_bases.size( );
    for (unsigned int jj=TRIM_LEFT; jj<cg_bases.size( )-TRIM_RIGHT; jj++) {
      if ( QUALS != "" && (int)quals[ii][jj] < MIN_QUAL )
	continue;
      tot_valid++;
      char base = as_base( cg_bases[jj] );
      if ( 'A' == base )
	totA++;
      else if ( 'C' == base )
	totC++;
      else if ( 'G' == base )
	totG++;
      else if ( 'T' == base )
	totT++;
      else
	ForceAssert( 1 == 0 );
    }
    if ( PER_SEQUENCE ) {
      out << ii << ' ' << totA << ' ' << totC << ' ' << totG << ' ' << totT << endl;
      totA = totC = totG = totT = 0;
    }
  }

  if ( ! PER_SEQUENCE ) {
    // Print result.
    out << "\nArguments:\n"
	<< "BASES=" << BASES << "\n";
    if ( SELECT != "" )
      out << "SELECT=" << SELECT << "\n";
    if ( QUALS != "" )
      out << "QUALS=" << QUALS << "\n"
	  << "MIN_QUAL=" << MIN_QUAL << "\n";
    
    String str_of;
    if ( QUALS != "" ) {
      double ratio_valid = 100.0 * SafeQuotient( tot_valid, tot );
      out << "\nTotal bases found: " << tot << "\n"
	  << "Total bases with qual >= " << MIN_QUAL
	  << ": " << tot_valid
	  << " (" << ToString( ratio_valid, 2 )
	  << "% of the total)\n";
      str_of = "% of the bases with qual >= " + ToString( MIN_QUAL );
    }
    else {
      out << "\nTotal bases found: " << tot_valid << "\n";
      str_of = "% of the total";
    }

    String str_ratioA = ToString( 100.0 * SafeQuotient( totA, tot_valid ), 2 );
    String str_ratioC = ToString( 100.0 * SafeQuotient( totC, tot_valid ), 2 );
    String str_ratioG = ToString( 100.0 * SafeQuotient( totG, tot_valid ), 2 );
    String str_ratioT = ToString( 100.0 * SafeQuotient( totT, tot_valid ), 2 );
    out << "As: " << totA << " (" << str_ratioA << str_of << ")\n"
	<< "Cs: " << totC << " (" << str_ratioC << str_of << ")\n"
	<< "Gs: " << totG << " (" << str_ratioG << str_of << ")\n"
	<< "Ts: " << totT << " (" << str_ratioT << str_of << ")\n"
	<< endl;
    out << "GCPercent: " << 100.0 * SafeQuotient(totC+totG, tot_valid) << endl;
    
    // Done.
    cout << Date( ) << ": done" << endl;
  }

}
