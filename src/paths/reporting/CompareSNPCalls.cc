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
#include "util/RunCommand.h"
// MakeDepend: dependency DisplayPileupsOnReference

/**
 * MapToContig
 */
int MapToContig( const look_align_plus &hit, const int refpos )
{
  const align &al = hit.a;
  bool rc1 = hit.Rc1( );
  int pos1 = al.pos1( );
  int pos2 = al.pos2( );
  int clen = hit.query_length;
  int contig_id = hit.query_id;

  // Move on the align.
  for (int block_id=0; block_id<al.Nblocks( ); block_id++) {
    
    // Move on gap.
    int gap = al.Gaps( block_id );
    if ( gap < 0 )
      pos1 += -gap;
    if ( gap > 0 ) {
      for (int ii=0; ii<gap; ii++) {
	if ( pos2 == refpos ) 
	  return ( rc1 ? clen - pos1 - 1 : pos1 );
	pos2++;
      }
    }
    
    // Move on length.
    int length = al.Lengths( block_id );      
    for (int ii=0; ii<length; ii++) {
      int realpos = rc1 ? clen - pos1 - 1 : pos1;
      if ( pos2 == refpos )
	return ( rc1 ? clen - pos1 - 1 : pos1 );
      pos1++;
      pos2++;
    }
  }
  
  // An error.
  return -1;

}

/**
 * MapToRef
 */
int MapToRef( const look_align_plus &hit, const int cgpos )
{
  const align &al = hit.a;
  bool rc1 = hit.Rc1( );
  int pos1 = al.pos1( );
  int pos2 = al.pos2( );
  int clen = hit.query_length;

  // Move on align.
  for (int block_id=0; block_id<al.Nblocks( ); block_id++) {
      
    // Move on gap.
    int gap = al.Gaps( block_id );
    if ( gap < 0 ) {
      for (int ii=0; ii<-gap; ii++) {
	if ( cgpos == ( rc1 ? clen - pos1 - 1 : pos1 ) )
	  return pos2;
	pos1++;
      }
    }
    if ( gap > 0 )
      pos2 += gap;
      
    // Move on length.
    int length = al.Lengths( block_id );      
    for (int ii=0; ii<length; ii++) {
      if ( cgpos == ( rc1 ? clen - pos1 - 1 : pos1 ) )
	return pos2;
      pos1++;
      pos2++;
    }
  }  
  
  // An error.
  return -1;

}

/**
 * CompareSNPCalls
 *
 * Compare SNPs called by EstimatePolymorphism against a given SNPs
 * database (saved as a plain text file).
 *
 * SNPS_DB: the SNPs data base (as a plain text file)
 * REF_HEAD: if given, show pileups against ref (huge output)
 */
int main(int argc, char *argv[])
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SNPS_DB );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
  CommandArgument_String_OrDefault( REF_HEAD, "" );
  EndCommandArguments;

  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

  String assembly_head = sub_dir + "/" + ASSEMBLY;
  String contigs_fastb_file = assembly_head + ".contigs.fastb";
  String qlt_file = assembly_head + ".Eval/aligns.qlt";
  String annotations_file = assembly_head + ".annotations";
  
  // Load.
  if ( ! IsRegularFile( qlt_file ) ) {
    cout << "Fatal error, did you run EvalScaffolds? Leaving.\n" << endl;
    return 1;
  }
  vec<look_align_plus> hits;
  LoadLookAlignPlus( qlt_file, hits );

  if ( ! IsRegularFile( annotations_file ) ) {
    cout << "Fatal error, annotations not found. Leaving.\n" << endl;
    return 1;
  }
  VecIntVec annotations;
  annotations.ReadAll( annotations_file );

  if ( ! IsRegularFile( SNPS_DB ) ) {
    cout << "Fatal error, SNPs database not found. Leaving.\n" << endl;
    return 1;
  }
  vec<int> snps_orig;
  ifstream snps_in ( SNPS_DB.c_str( ) );
  while ( 1 ) {
    int pos;
    snps_in >> pos;
    if ( ! snps_in ) break;
    snps_orig.push_back( pos );
  }
  snps_in.close( );

  // SNPs from the data base as a VecIntVec.
  VecIntVec annotations_db;
  annotations_db.resize( annotations.size( ) );
  for (size_t ii=0; ii<annotations.size( ); ii++)
    annotations_db[ii].resize( annotations[ii].size( ) );
  
  // Map contig id to hit for the contig.
  vec<int> to_hit( annotations.size( ), -1 );
  
  // Counters.
  int n_good = 0;      // snp on ref mapped on a unique contig
  int n_missing = 0;   // snp not captured by any contig
  int n_multiple = 0;  // snp mapped on two overlapping contigs
  
  // Convert snps_orig into contig coordinates.
  for (int snp_id=0; snp_id<snps_orig.isize( ); snp_id++) {
    int dbsnp = snps_orig[snp_id];
    
    // Find align containing dbsnp.
    int hit_id = -1;
    for (int ii=0; ii<hits.isize( ); ii++) {
      if ( hits[ii].a.pos2( ) > dbsnp ) continue;
      if ( hits[ii].a.Pos2( ) <= dbsnp ) continue;
      if ( hit_id > -1 ) {
	hit_id = -2;
	break;
      }
      hit_id = ii;
    }
    
    // No contig contains this snp.
    if ( hit_id == -1 ) {
      n_missing++;
      continue;
    }

    // More than one contig contains this snp.
    if ( hit_id == -2 ) {
      n_multiple++;
      continue;
    }

    // Found contig containing snp.
    const look_align_plus &hit = hits[hit_id];

    // Update counters.
    n_good++;
    int cgpos = MapToContig( hit, dbsnp );
    annotations_db[hit.query_id][cgpos] = 2;
    to_hit[hit.query_id] = hit_id;
  }
  
  // Compare annotations.
  int n_confirmed = 0;
  int n_missing2 = 0;
  int n_unverified = 0;
  int n_unmapped = 0;
  for (int contig_id=0; contig_id<(int)annotations.size( ); contig_id++) {
    for (int pos=0; pos<(int)annotations[contig_id].size( ); pos++) {
      int annot_db = annotations_db[contig_id][pos];
      int annot_ass = annotations[contig_id][pos];
      if ( annot_db == annot_ass ) {
	if ( annot_db == 2 ) n_confirmed++;
	continue;
      }
      if ( annot_db != 2 && annot_ass != 2 ) continue;
      int hit_id = to_hit[contig_id];
      if ( hit_id < 0 ) {
	n_unmapped++;
	continue;
      }

      // A discrepancy with the given db.
      int on_ref = MapToRef( hits[hit_id], pos );
      if ( annot_db == 2 ) {
	n_missing2++;
	cout << " missing snp at ";
      }
      else {
	n_unverified++;
	cout << " unverified snp at ";
      }
      cout << on_ref << "  =  c" << contig_id << " at " << pos << "\n";

      // Print pileup.
      if ( REF_HEAD == "" ) continue;
      
      const align &al = hits[hit_id].a;
      int beg = Max( al.pos1( ), pos - 5 );
      int end = Min( al.Pos1( ), pos + 5 );
      int range1 = on_ref - ( pos - beg );
      int range2 = on_ref + ( end - pos );
      String pdr = "PRE=" + PRE + " DATA=" + DATA + " RUN=" + RUN;
      String theCommand
	= " DisplayPileupsOnReference " + pdr
	+ " REF_HEAD=" + REF_HEAD
	+ " ANNOTATIONS=" + annotations_file
	+ " RANGE=0:" + ToString( range1 ) + "-" + ToString( range2 )
	+ " BRIEF=True";
      RunCommand( theCommand );

      cout << "\n";
      for (int ii=0; ii<80; ii++) cout << "=";
      cout << "\n" << endl;
    }
  }
  cout << endl;
  
  // This should not happen.
  if ( n_unmapped > 0 )
    cout << "Fatal error: unmapped SNPs reported.\n" << endl;
  
  // Final stats.
  double ratio_confirmed = 100.0 *  SafeQuotient( n_confirmed, n_good );
  double ratio_missing2 = 100.0 *  SafeQuotient( n_missing2, n_good );
  cout << "OVERALL STATISTICS\n"
       << "\n"
       << "      SNPs in the given database: " << snps_orig.size( ) << "\n"
       << "                in assembly gaps: " << n_missing << "\n"
       << "      mapped to multiple contigs: " << n_multiple << "\n"
       << " mapped to unique contigs (good): " << n_good << "\n"
       << "\n"
       << "   confirmed SNPS (in both sets): " << n_confirmed
       << " (" << ToString( ratio_confirmed, 1 ) << "% of good)\n"
       << " missing SNPS (only in given db): " << n_missing2
       << " (" << ToString( ratio_missing2, 1 ) << "% of good)\n"
       << "     unverified SNPS (not in db): " << n_unverified << "\n"
       << endl;
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}
