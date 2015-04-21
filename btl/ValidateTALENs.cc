///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Basevector.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "btl/CMonoRep.h"
#include "btl/LoaderTALENs.h"
#include "btl/MonomerUtilsTALENs.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/AlignTwoBasevectors.h"
#include "pairwise_aligners/PerfectAlignment.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "tiled/AlignsOnConsensus.h"
#include "util/RunCommand.h"
// MakeDepend: dependency Fastb
// MakeDepend: dependency Qualb
// MakeDepend: dependency FastAlignShortReads

/**
 * Validate
 */
void Validate( const String &AB1_DIR,
	       const String &OUT_DIR,
	       const String &parts_file,
	       const String &prod_name,
	       const String &reference,
	       const vecbvec &primers,
	       const vecbvec &parts,
	       const vecString &parts_ids,
	       const bool &R2 )
{
  // File names.
  String strRC = ( R2 ? String( "R2" ) : String( "R1" ) );
  String ab1F1_file = AB1_DIR + "/" + prod_name + "-F1.ab1";
  String ab1F2_file = AB1_DIR + "/" + prod_name + "-F2.ab1";
  String ab1R_file = AB1_DIR + "/" + prod_name + "-" + strRC + ".ab1";
  
  String out_dir = OUT_DIR + "/" + prod_name;
  String log_file = out_dir + "/main.log";
  String csvline_file = out_dir + "/csvline.txt";

  String ref_head = out_dir + "/ref";

  String tmp_file =  out_dir + "/reads.tmp";
  String fasta_file = out_dir + "/reads.fasta";
  String qual_file = out_dir + "/reads.qual";
  String fastb_file = out_dir + "/reads.fastb";
  String qualb_file = out_dir + "/reads.qualb";

  String consensus_head = out_dir + "/consensus";
  String cbases_file = out_dir + "/consensus.fastb";
  String cquals_file = out_dir + "/consensus.qualb";
  
  String onreads_qlt_file = out_dir + "/on_reads.qlt";
  String onconsensus_qlt_file = out_dir + "/on_consensus.qlt";
  String visual_file = out_dir + "/aligns.visual";
  
  String devnull = "/dev/null";

  Mkpath( out_dir );
  ofstream log( log_file.c_str( ) );
  log << "Validation of " << prod_name << "\n";

  // Clean from old runs, if needed.
  if ( IsRegularFile( tmp_file ) ) Remove ( tmp_file );
  if ( IsRegularFile( fasta_file ) ) Remove( fasta_file );
  if ( IsRegularFile( qual_file ) ) Remove( qual_file );
  
  // Fasta first.
  RunCommand( "echo \">\"" + prod_name + "-F1 >> " + tmp_file );
  RunCommand( "extract_seq " + ab1F1_file + " >> " + tmp_file );

  RunCommand( "echo \">\"" + prod_name + "-" + strRC + " >> " + tmp_file );
  RunCommand( "extract_seq " + ab1R_file + " >> " + tmp_file );

  if ( IsRegularFile( ab1F2_file ) ) {
    RunCommand( "echo \">\"" + prod_name + "-F2 >> " + tmp_file );
    RunCommand( "extract_seq " + ab1F2_file + " >> " + tmp_file );
  }

  // Replace "-" with "N" (extract_seq uses "-" instead of "N").
  String line;
  ofstream out( fasta_file.c_str( ) );
  ifstream in( tmp_file.c_str( ) );
  while ( in ) {
    getline( in, line );
    if ( ! in ) break;
    if ( ! line.Contains( ">", 0 ) )
      for (int ii=0; ii<line.isize( ); ii++)
	if ( line[ii] == '-' ) line[ii] = 'N';
    out << line << "\n";
  }
  in.close( );
  out.close( );

  // Remove tmp.
  if ( IsRegularFile( tmp_file ) ) Remove( tmp_file );

  // Qual file.
  RunCommand( "echo \">\"" + prod_name + "-F1 >> " + qual_file );
  RunCommand( "extract_qual " + ab1F1_file + " >> " + qual_file );
  
  RunCommand( "echo \">\"" + prod_name + "-" + strRC + " >> " + qual_file );
  RunCommand( "extract_qual " + ab1R_file + " >> " + qual_file );
  
  if ( IsRegularFile( ab1F2_file ) ) {
    RunCommand( "echo \">\"" + prod_name + "-F2 >> " + qual_file );
    RunCommand( "extract_qual " + ab1F2_file + " >> " + qual_file );
  }

  // Fastb and qualb.
  RunCommandWithLog( "Fastb FASTAMB=True PRE= FILE=" + fasta_file, devnull );
  RunCommandWithLog( "Qualb PRE= FILE=" + qual_file, devnull );
  
  // Tag terminator parts.
  vec<bool> is_termin( parts_ids.size( ), false );
  for (size_t ii=0; ii<parts_ids.size( ); ii++)
    if ( parts_ids[ii].Contains( "termin" ) )
      is_termin[ii] = true;
  
  // Load bases and quals, sanity check.
  vecbvec rbases( fastb_file );
  vecqvec rquals( qualb_file );
  if ( rbases.size( ) != rquals.size( ) ) {
    log << "\nFATAL ERROR: bases and quals have different sizes" << endl;
    PrintErrorCsvLine( csvline_file, prod_name, reference );
    return;
  }
  if ( rbases.size( ) != 2 && rbases.size( ) != 3 ) {
    log << "\nFATAL ERROR: bases and quals must have size 2 or 3" << endl;
    PrintErrorCsvLine( csvline_file, prod_name, reference );
    return;
  }
  for (int ii=0; ii<(int)rbases.size( ); ii++) {
    if ( rbases[ii].size( ) != rquals[ii].size( ) ) {
      log << "\nFATAL ERROR: bases/quals entries have different sizes" << endl;
      PrintErrorCsvLine( csvline_file, prod_name, reference );
      return;
    }
  }

  // Visual stream (to print various aligns).
  ofstream vis_out( visual_file.c_str( ) );
  
  // Turn reference into a CMonoRep, and save bases/quals.
  CMonoRep mr_ref( "ref", &log );
  mr_ref.SetFromReference( reference, primers, parts, parts_ids );
  mr_ref.SaveBasesAndQuals( ref_head );
  mr_ref.Print( parts_ids, log );

  // Align parts to reads.
  String theCommand
    = "FastAlignShortReads TARGET=" + fastb_file
    + " QUERY=" + parts_file
    + " OUT=" + onreads_qlt_file
    + " K=8 MAX_ERRORS=64 MAX_INDELS=16";
  RunCommandWithLog( theCommand, devnull );
  
  // Reload aligns, and place monomers/terminators.
  vec<look_align_plus> hits;
  LoadLookAlignPlus( onreads_qlt_file, hits );
  
  // Filter aligns, and find chains.
  vec<look_align_plus> chains;
  AlignsToChains( hits, is_termin, chains );
  
  // Print chains.
  log << "\nAlignments of monomers/terminators to reads:\n\n";
  PrintChains( chains, is_termin, parts_ids, log );
  
  // Turn chains in CMonoReps.
  CMonoRep mr_fw1( "fw1", &log );
  mr_fw1.SetFromHits( chains, parts, rbases[0], rquals[0], 0, false );
  mr_fw1.VerbosePrint( parts_ids, parts, vis_out );

  CMonoRep mr_fw2( "fw2", &log );
  if ( rbases.size( ) > 2 ) {
    mr_fw2.SetFromHits( chains, parts, rbases[2], rquals[2], 2, false );
    mr_fw2.VerbosePrint( parts_ids, parts, vis_out );
  }

  CMonoRep mr_rc( "rc", &log );
  mr_rc.SetFromHits( chains, parts, rbases[1], rquals[1], 1, true );
  mr_rc.VerbosePrint( parts_ids, parts, vis_out );

  // Remove nasty reads, if possible.
  if ( mr_fw1.IsAwful( ) ) {
    log << "\nWarning: fw1 is awful, discarding it" << endl;
    mr_fw1.Clear( );
  }
  if ( mr_fw2.IsAwful( ) ) {
    log << "\nWarning: fw2 is awful, discarding it" << endl;
    mr_fw2.Clear( );
  }
  if ( mr_rc.IsAwful( ) ) {
    log << "\nWarning: rc is awful, discarding it" << endl;
    mr_rc.Clear( );
  }

  // Merge chains (build draft consensus).
  CMonoRep mr_consensus( "consensus", &log );
  bool ok = false;
  if ( rbases.size( ) < 3 ) {
    ok = mr_consensus.Merge( mr_fw1, mr_rc );
    mr_consensus.Print( parts_ids, log );
  }
  else {
    int off_12 = 0;
    int mat_12 = 0;
    int mis_12 = 0;
    bool ok_12 =  Offset( mr_fw1, mr_fw2, off_12, 0, &mat_12, &mis_12 );

    int off_1r = 0;
    int mat_1r = 0;
    int mis_1r = 0;
    bool ok_1r =  Offset( mr_fw1, mr_rc, off_1r, 0, &mat_1r, &mis_1r );
    
    // Decide if we first merge fw1-fw, or fw1-rc. HEURISTICS embedded here!
    int len_12 = mat_12 + mis_12;
    bool poor_12 = ( (!ok_12) || ( len_12 < 3 && mis_12 > 0 ) );
    bool good_1r = ( ok_1r && mat_1r >=2 && mis_1r < 1 )
      || ( ok_1r && mat_1r >=3 && mis_1r < 2 );
    
    // Merge fw1-rc first, and then add fw2.
    if ( poor_12 && good_1r ) {
      CMonoRep mr_inter( "inter", &log );
      ok = mr_inter.Merge( mr_fw1, mr_rc );
      mr_inter.Print( parts_ids, log );
      
      if ( ok ) ok = mr_consensus.Merge( mr_inter, mr_fw2 );
      mr_consensus.Print( parts_ids, log);
    }

    // Merge fw1-fw2 first, and then add rc.
    else {
      CMonoRep mr_inter( "inter", &log );
      ok = mr_inter.Merge( mr_fw1, mr_fw2 );
      mr_inter.Print( parts_ids, log );
      
      if ( ok ) ok = mr_consensus.Merge( mr_inter, mr_rc );
      mr_consensus.Print( parts_ids, log);
    }
    
  }
  if ( ! ok ) {
    log << "\nFATAL ERROR: merging failed" << endl;
    PrintErrorCsvLine( csvline_file, prod_name, reference );
    return;
  }
  
  // Refine consensus.
  {
    mr_consensus.Refine( mr_fw1, mr_fw2, mr_rc, consensus_head ); 
    mr_consensus.SaveBasesAndQuals( consensus_head );

    theCommand
      = "FastAlignShortReads TARGET=" + cbases_file
      + " QUERY=" + parts_file
      + " OUT=" + onconsensus_qlt_file
      + " K=8 MAX_ERRORS=64 MAX_INDELS=16";
    RunCommandWithLog( theCommand, devnull );
    
    bvec cbases = mr_consensus.Bases( );
    qvec cquals = mr_consensus.Quals( );
    
    vec<look_align_plus> hits;
    LoadLookAlignPlus( onconsensus_qlt_file, hits );
    
    vec<look_align_plus> chains;
    AlignsToChains( hits, is_termin, chains );
    
    mr_consensus.SetFromHits( chains, parts, cbases, cquals, 0, false );
  }

  // Evaluate on reference.
  log << "\n";
  mr_consensus.EvalOnRef( mr_ref, parts_ids, prod_name, csvline_file, vis_out );
  
  // Close streams.
  vis_out.close( );
  log.close( );
  
}



/**
 * ValidateTALENs
 *
 * Merge the two (or three) reads sequenced from a TALEN product, and
 * validate the consensus against a set of given parts.
 *
 * INPUT
 *
 *   <AB1_DIR>: dir with the .ab1 files (traces). These are converted
 *      into .fasta and .qual by the io-lib utilities extract_seq and
 *      extract_qual.  NOTE: the assumption is that there are two
 *      paired reads for each product, with names <name>-F1.ab1 and
 *      <name>-R1.ab1.  <name> must be unique (two different products
 *      must have different <name>s), and it is always assumed that F1
 *      is the fw read, and that R1 is the rc read. Optionally, a
 *      third (fw) read can be added to the set, and this must be
 *      called <name>-F2.ab1.
 * 
 *   <PARTS_DIR>: dir with the fasta of the monomers (and terminators)
 *      used to construct (assemble) the TALENs. These must end with
 *      ".fasta".
 *
 *   <PRIMERS_FASTA>: fasta file with the 4 bases primers.
 *
 *   <REF>: flat file with references (in monomer space). The file
 *      must contain a list of lines (one per product), each
 *      containing comma separated tag name and monomers, as in:
 *          tale7,NG,HD,NI,NN,NN,HDcseq,NG,HD,NG,NN,NN,NI,NI
 *      NB: the last monomer is assumed to be a terminator, and the
 *          name of the product must match the given <name>.
 *
 *   <PRODUCT>: if not empty, validate this product only.
 *
 *   <SKIP>: if not empty, skip selected product.
 *
 *   <R2>: if true, look for <name>-R2.ab1 reads, rather than
 *      <name>-R1.ab1 (for rc reads only).
 *
 * OUTPUT
 *
 *   <OUT_DIR>: where output will be stored.
 *
 * SANTEMP - TODO:
 *
 *   Improve monomer-based aligner (SmithWaterman) in Validate( ).
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;
  CommandArgument_String( AB1_DIR );
  CommandArgument_String( PARTS_DIR );
  CommandArgument_String( PRIMERS_FASTA );
  CommandArgument_String( OUT_DIR );
  CommandArgument_String( REF );
  CommandArgument_String_OrDefault( PRODUCT, "" );
  CommandArgument_String_OrDefault( SKIP, "" );
  CommandArgument_Bool_OrDefault( R2, False );
  EndCommandArguments;

  // Make sure extract_qual and extract_seq are in the path.
  TestExecutableByWhich( "extract_qual" );
  TestExecutableByWhich( "extract_seq" );
  
  // Make output dir and log stream.
  Mkpath( OUT_DIR );
  
  // Load products and references (the latter only if REF != "skip").
  vec<String> prod_names;
  LoadProducts( AB1_DIR, R2, prod_names, cout );
  if ( PRODUCT != "" ) {
    prod_names.clear( );
    prod_names.resize( 1, PRODUCT );
  }
  
  vec<String> references;
  if ( ! LoadReferences( REF, prod_names, references, &cout ) )
    return 1;   // Exit on failure!
  
  // Load parts, and parts ids (used to tag terminators).
  String parts_file = LoadParts( PARTS_DIR, &cout );
  String parts_ids_file = parts_file.SafeBefore( ".fastb" ) + ".ids";
  if ( parts_file == "" ) return 1;   // Exit on failure!

  vecbvec parts( parts_file );
  vecString parts_ids( parts_ids_file );
  {
    for (size_t ii=0; ii<parts_ids.size( ); ii++) {
      bool monomer = parts_ids[ii].Contains( "monomer" );
      bool termin = parts_ids[ii].Contains( "termin" );
      ForceAssert( monomer || termin );
    }
  }

  // Load primers.
  vecbvec primers;
  LoadPrimers( PRIMERS_FASTA, primers );

  // Loop over all products.
  for (size_t prod_id=0; prod_id<prod_names.size( ); prod_id++) {
    const String &name = prod_names[prod_id];
    const String &ref = references[prod_id];
    if ( PRODUCT != "" && name != PRODUCT ) continue;
    if ( SKIP != "" && name == SKIP ) continue;
    
    cout << Date( ) << ": validating " << name << endl;
    Validate(AB1_DIR,OUT_DIR,parts_file,name,ref,primers,parts,parts_ids,R2);
  }
  cout << endl;

  // Done.
  cout << Date( ) << ": ValidateTALENs done" << endl;
  
}

