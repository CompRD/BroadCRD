///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Intvector.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "paths/ReadLoc.h"

/**
 * StringPileup
 *
 * Print pileup, rc-ing contig if needed.
 */
String StringPileup( const bool rc1,
		     const int clen,
		     const int pos,
		     const vec<dumbcall> &calls,
		     const bool brief = false,
		     const IntVec *pann = 0 )
{
  const dumbcall &column = rc1 ? calls[clen-pos-1] : calls[pos];

  String pileup;
  pileup.reserve( 1 + column.CountAll( ) );

  if ( pann ) {
    int as_int = rc1 ? (*pann)[clen-pos-1] : (*pann)[pos];
    if ( as_int == 0 ) pileup += "x ";
    else if ( as_int == 1 ) pileup += ". ";
    else if ( as_int == 2 ) pileup += "s ";
    else ( ForceAssert( 1 == 0 ) );
  }
  
  int As = rc1 ? column.base[3] : column.base[0];
  int Cs = rc1 ? column.base[2] : column.base[1];
  int Gs = rc1 ? column.base[1] : column.base[2];
  int Ts = rc1 ? column.base[0] : column.base[3];
  int Ds = column.base[4];
  int Is = column.base[5];

  if ( brief ) {
    int ntot = column.CountAll( );
    double ratioA = 100.0 * SafeQuotient( As, ntot );
    double ratioC = 100.0 * SafeQuotient( Cs, ntot );
    double ratioG = 100.0 * SafeQuotient( Gs, ntot );
    double ratioT = 100.0 * SafeQuotient( Ts, ntot );
    double ratioD = 100.0 * SafeQuotient( Ds, ntot );
    double ratioI = 100.0 * SafeQuotient( Is, ntot );

    if ( As > 0 )
      pileup += ToString( As ) + "A (" + ToString( ratioA, 1 ) + "%), "; 
    if ( Cs > 0 )
      pileup += ToString( Cs ) + "C (" + ToString( ratioC, 1 ) + "%), "; 
    if ( Gs > 0 )
      pileup += ToString( Gs ) + "G (" + ToString( ratioG, 1 ) + "%), "; 
    if ( Ts > 0 )
      pileup += ToString( Ts ) + "T (" + ToString( ratioT, 1 ) + "%), "; 
    if ( Ds > 0 )
      pileup += ToString( Ds ) + "D (" + ToString( ratioD, 1 ) + "%), "; 
    if ( Is > 0 )
      pileup += ToString( Is ) + "I (" + ToString( ratioI, 1 ) + "%), "; 
    
    if ( column.CountAll( ) > 0 )
      pileup.resize( pileup.size( ) - 2 );   // remove trailing ", "
  }
  else {
    for (int ii=0; ii<As; ii++) pileup += "A";
    for (int ii=0; ii<Cs; ii++) pileup += "C";
    for (int ii=0; ii<Gs; ii++) pileup += "G";
    for (int ii=0; ii<Ts; ii++) pileup += "T";
    for (int ii=0; ii<column.base[4]; ii++) pileup += "D";
    for (int ii=0; ii<column.base[5]; ii++) pileup += "I";
  }
  
  return pileup;
}

/**
 * DisplayPileupsOnReference
 *
 * It prints pilups of reads onto the portion of contig aligning the
 * selected range, in the same format as DisplayLocs (run with
 * SHOW_PILEUPS=True). It expects the output from EvalScaffolds.
 *
 * REF_HEAD: relative to DATA
 * RANGE: eg "T:beg-end" (show bases [beg, end) on target T)
 * ANNOTATIONS: optionally load and show annotations (VecIntVec)
 * BRIEF: print pile ups in compact form
 */
int main(int argc, char *argv[])
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( RANGE );
  CommandArgument_String_OrDefault( REF_HEAD, "genome" );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
  CommandArgument_String_OrDefault( ANNOTATIONS, "" );
  CommandArgument_Bool_OrDefault( BRIEF, False );
  EndCommandArguments;

  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

  String genome_fastb_file = data_dir + "/" + REF_HEAD + ".fastb";

  String assembly_head = sub_dir + "/" + ASSEMBLY;
  String contigs_fastb_file = assembly_head + ".contigs.fastb";
  String qlt_file = assembly_head + ".Eval/aligns.qlt";
  
  // Parse range.
  int t_id = RANGE.Before( ":" ).Int( );
  int t_beg = RANGE.After( ":" ).Before( "-" ).Int( );
  int t_end = RANGE.After( "-" ).Int( );

  // Load aligns.
  if ( ! IsRegularFile( qlt_file ) ) {
    cout << "Fatal error, did you run EvalScaffolds? Leaving.\n" << endl;
    return 1;
  }
  vec<int> t_sel( 1, t_id );
  vec<look_align_plus> hits;
  LoadLookAlignPlus( qlt_file, hits, 0, &t_sel );

  VecIntVec annotations;
  if ( ANNOTATIONS != "" )
    annotations.ReadAll( ANNOTATIONS );
  
  // Find align containing selected range.
  int hit_id = -1;
  for (int ii=0; ii<hits.isize( ); ii++) {
    if ( hits[ii].a.pos2( ) > t_beg ) continue;
    if ( hits[ii].a.Pos2( ) < t_end ) continue;
    if ( hit_id > -1 ) {
      cout << "Fatal error: two candidates found for RANGE. Leaving." << endl;
      return 1;
    }
    hit_id = ii;
  }
  const look_align_plus &hit = hits[hit_id];
  const align &al = hit.a;
  vec<int> q_sel( 1, hit.query_id );

  // Load remaining data.
  read_locs_on_disk locs_loader( assembly_head, run_dir );

  vecbvec genome;
  genome.SparseRead( genome_fastb_file, t_sel, 0 );

  vecbvec contigs;
  contigs.SparseRead( contigs_fastb_file, q_sel, 0 );
  
  // Log some info.
  cout << "c" << hit.query_id << " aligns " << RANGE << "\n" << endl;
  
  // Generate pileups.
  vec< vec<String> > table;
  String sIndel = "-";
  String sMismatch = "*";
  String s0 = "";

  bvec contig_rc;
  if ( hit.Rc1( ) ) {
    contig_rc = contigs[hit.query_id];
    contig_rc.ReverseComplement( );
  }
  const bvec &contig = hit.Rc1( ) ? contig_rc : contigs[hit.query_id];
  const bvec &ref = genome[t_id];
  const IntVec *pann = ANNOTATIONS != "" ? &( annotations[hit.query_id] ): 0;
  
  vec<read_loc> locs;
  locs_loader.LoadContig( hit.query_id, locs );
  
  vec<dumbcall> calls;
  Pileup( contig, locs, run_dir, calls );
  
  bool rc1 = hit.Rc1( );
  int pos1 = al.pos1( );
  int pos2 = al.pos2( );
  int clen = contig.size( );
  for (int block_id=0; block_id<al.Nblocks( ); block_id++) {
    if ( pos2 > t_end ) break;
  
    // Move on gap.
    int gap = al.Gaps( block_id );
    if ( gap < 0 ) {
      for (int ii=0; ii<-gap; ii++) {
	if ( pos2 >= t_beg && pos2 < t_end ) {
	  String sPos = ToString( rc1 ? clen - pos1 - 1 : pos1 );
	  if ( rc1 ) sPos += "rc";
	  String sTig = String( as_base( contig[pos1] ) );
	  String pileups = StringPileup( rc1, clen, pos1, calls, BRIEF, pann );
	  table.push_back( MkVec( s0, sPos, sIndel, s0, sTig, pileups ) );
	}
	pos1++;
      }
    }
    if ( gap > 0 ) {
      for (int ii=0; ii<gap; ii++) {
	if ( pos2 >= t_beg && pos2 < t_end ) {
	  String gPos = ToString( pos2 );
	  String sRef = String( as_base( ref[pos2] ) );
	  table.push_back( MkVec( gPos, s0, sIndel, sRef, s0, s0 ) );
	}
	pos2++;
      }
    }
    
    // Move on length.
    int length = al.Lengths( block_id );      
    for (int ii=0; ii<length; ii++) {
      if ( pos2 >= t_beg && pos2 < t_end ) {
	String gPos = ToString( pos2 );
	String sPos = ToString( rc1 ? clen - pos1 - 1 : pos1 );
	if ( rc1 ) sPos += "rc";
	String sTig = String( as_base( contig[pos1] ) );
	String sRef = String( as_base( ref[pos2] ) );
	String sDecor = ( sTig == sRef ? s0 : sMismatch );
	String pileups = StringPileup( rc1, clen, pos1, calls, BRIEF, pann );
	table.push_back( MkVec( gPos, sPos, sDecor, sRef, sTig, pileups ) );
      }
      pos1++;
      pos2++;
    }

  }
  
  // Print pileups.
  PrintTabular( cout, table, 2, "rllll" );
  
}
