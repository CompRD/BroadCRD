///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "Vec.h"
#include "FeudalMimic.h"
#include "lookup/LookAlign.h"
#include "random/Random.h"
#include "PairsManager.h"
#include "feudal/BinaryStream.h"

/**
 * SimulatePolymorphicReads
 *
 * Simulate diploid reads from haploid reads.
 * 
 * Input consists of a polymorphic reference (the output of
 * MutateReference), plus reads (the output from build_micro, run with
 * DATA_ONLY set to true). SimulatePolymorphicReads will randomly
 * assign an haplotype to each pair of reads, and change single bases
 * on the reads to match the haplotype on the diploid reference.
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( REFERENCE_NAME )
  CommandArgument_String( DATA_IN );
  CommandArgument_String( DATA_OUT );
  CommandArgument_String( RANGE )
  CommandArgument_String_OrDefault( REF_HEAD, "genome_diploid" );
  CommandArgument_String_OrDefault( FRAG, "frag_reads_orig" );
  CommandArgument_String_OrDefault( JUMP, "jump_reads_orig" );
  CommandArgument_String_OrDefault( LONG_JUMP, "long_jump_reads_orig" );
  CommandArgument_UnsignedInt_OrDefault( SEED, 666 );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  EndCommandArguments;
  

  srandomx(SEED);

  // Dir and file names.
  String full_base     = PRE + "/" + REFERENCE_NAME;
  String full_data_in  = full_base + "/" + DATA_IN;
  String full_data_out = full_base + "/" + DATA_OUT;

  String dip_genome_head       = full_base + "/genome_diploid";
  String dip_genome_fastb_file = dip_genome_head + ".fastb";
  String dip_genome_snps_file  = dip_genome_fastb_file + ".snp_locs";
  String genome_fastb_file     = full_base + "/genome.fastb";

  String log_file = full_data_out + "/SimulatePolymorphicReads.log";

  Mkpath( full_data_out );
  //ofstream log( log_file.c_str() );
  ostream& log = cout;

  PrintCommandPretty( log );
  //cout << "Sending log to " << log_file << "\n" << endl;

  // Three sets of reads (not all need to exist).
  int n_found = 0;
  vec<String> file_heads;
  vec<String> qlt_files;
  for (int ii=0; ii<3; ii++) {
    String head = full_data_in + "/";
    if ( ii == 0 )      head += FRAG;
    else if ( ii == 1 ) head += JUMP;
    else                head += LONG_JUMP;
    file_heads.push_back( head );
    String qlt_file = head + ".qltout";
    if ( IsRegularFile( qlt_file ) ) {
      n_found++;
      qlt_files.push_back( qlt_file );
      String short_file = qlt_file;
      while ( short_file.Contains( "/" ) )
	short_file = short_file.After( "/" );
      log << Date( ) << ": found " << short_file << endl;
    }
    else
      qlt_files.push_back( "" );
  }
  if ( n_found < 1 ) {
    String str_err = "Fatal error, no aligns found. Stopping now.";
    log << str_err << endl;
    cout << str_err << endl;
    return 1;
  }
  
  // Load diploid genome.
  vecbvec dip_genome( dip_genome_fastb_file );
  vecbvec genome( genome_fastb_file );

  size_t genomeSize = 0;
  for ( size_t i = 0; i < genome.size(); i++ )
    genomeSize += genome[i].size();
  cout << "Haploid size of genomic region: " << genomeSize << endl;

  vec< vec<int> > to_snps;
  Mimic( dip_genome, to_snps );
  for ( size_t i = 0; i < to_snps.size(); i++ )
    for ( size_t j = 0; j < to_snps[i].size(); j++ )
      to_snps[i][j] = -1;  // initialize, negative indicates no snp at this position

  vec< pair<int,int> > snp_locs;
  BinaryReader::readFile( dip_genome_snps_file, &snp_locs );
  int nSNPs = round( snp_locs.size()/ 2.0 );
  cout << "Number of SNPs in genomic region: " << nSNPs << endl;

  for (int ii=0; ii<snp_locs.isize( ); ii++){
    if (VERBOSE )
      cout << snp_locs[ii].first << "\t"
	   << snp_locs[ii].second << "\t"
	   << as_base( dip_genome[0][snp_locs[ii].second] ) << " "
	   << as_base( dip_genome[1][snp_locs[ii].second] ) << "\n";
    int chr = snp_locs[ii].first, pos = snp_locs[ii].second;
    to_snps.at( chr ).at( pos ) = ii;
  }

  
  String srange = RANGE.After(":");
  int gStart = 0;
  if ( srange != "" )
    gStart = atoi(srange.Before("-").c_str());

  size_t nmutated = 0, nModReads = 0;
  ulonglong nAlignedBases = 0;
  for (int ii=0; ii<3; ii++) {
    String head = "";
    if ( ii == 0 )      head = FRAG;
    else if ( ii == 1 ) head = JUMP;
    else                head = LONG_JUMP;
 
    String qlt_file       = full_data_in + "/" + head + ".qltout";
    String in_fastb_file  = full_data_in + "/" + head + ".fastb";
    String out_fastb_file = full_data_out + "/" + head + ".fastb";
    String in_pairs_file  = full_data_in + "/" + head + ".pairs";

    cout << Date() << ": reading sequence file " << in_fastb_file << endl;
    vecbvec seqs( in_fastb_file );

    if ( IsRegularFile( qlt_file ) ) { 
      vec<look_align> aligns;
      cout << Date() << ": reading alignment file " << qlt_file << endl;
      LoadLookAligns( qlt_file, aligns );
      int n_reads = seqs.size();
      vec<int> aligns_index( n_reads, -1 );
      for ( size_t i = 0; i < aligns.size(); i++ ){
	ForceAssertLt( aligns[i].QueryId(), n_reads );
	aligns_index.at( aligns[i].QueryId() ) = i;
      }
      
      PairsManager pairs( in_pairs_file );
      // introduce SNPs into paired reads
      for ( size_t pid = 0; pid < pairs.nPairs(); pid++ ){
	if ( randomx() % 2 == 0 )
	  continue;  // randomly mutate only half of pairs

	size_t qid = pairs.ID1( pid );
	for ( size_t it = 0; it <= 1; it++ ){
	  if ( it == 1 ) qid = pairs.ID2( pid );
	  
	  int nsame = 0, naligned = 0;
	  if ( aligns_index.at(qid) >= 0 ){
	    look_align la = aligns.at( aligns_index.at( qid ) );
	    if ( la.indels > 0 )
	      continue; // do not consider alignments with indels
	    int tstart = la.StartOnTarget() - gStart;
	    int tend   = la.EndOnTarget() - gStart;
	    
	    bvec qseq = seqs[qid];
	    if ( la.Rc1() ) qseq.ReverseComplement();
	    
	    for ( int tpos  = tstart; tpos < tend; tpos++ ){
	      if ( tpos < 0 || tpos >= dip_genome[0].isize() )
		continue; 
	      naligned++;
	      if ( dip_genome[0][tpos] == qseq[ tpos - tstart ] ) 
		nsame++;
	      if ( to_snps[0][tpos] >= 0 && dip_genome[0][tpos] == qseq[ tpos - tstart ] ){
		qseq.Set( tpos - tstart, dip_genome[1][tpos] ); 
		nmutated++;
	      }      
	    }
	    nAlignedBases += naligned;
	    float percIde = 100.0;
	    if ( naligned > 0 ) percIde *= (float)nsame/(float)naligned;
	    //PRINT6( qid, seqs[qid].size(), nsame, naligned, percIde, ToStringBool(la.Fw1()) );
	    if ( la.Rc1() ) qseq.ReverseComplement();
	    if ( seqs[qid] != qseq) 
	      nModReads++;
	    seqs[qid] = qseq;
	  }
	}
      }
    }
    seqs.WriteAll( out_fastb_file ); 
  }
  cout << "Number of bases aligned: " << nAlignedBases << endl;
  cout << "Number of bases mutated: " << nmutated << endl; 

  int nexpected = round( nAlignedBases * nSNPs / (double)genomeSize );
  cout << "Approx. expected number of mutated bases: " << nexpected << " (assuming no read errors)" << endl;
  cout << "Number of modified reads: " << nModReads << endl;
  // Done.
  cout << Date( ) << ": done" << endl;
  log << Date( ) << ": done" << endl;
  //log.close( );

}
