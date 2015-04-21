/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "ApplyCorrections.h"
#include "Basevector.h"
#include "Corrector.h"
#include "Qualvector.h"
#include "ReadLocation.h"
#include "STLExtensions.h"
#include "String.h"

/**
 * ApplyCorrections
 */
void ApplyCorrections( const vec<corrector> &fixes,
		       const vec<int> & first_locs,
		       vec<read_location> &locs,
		       vecbasevector &bases,
		       vecqualvector &quals,
		       ostream *snp_info_out )
{
  // Snp stream.
  ofstream devnull ( "/dev/null" );
  ostream &snpout = snp_info_out ? *snp_info_out : devnull;

  // Determine bases_out/quals_out raw size.
  longlong n_bases = 0;
  for (int ii=0; ii<(int)bases.size( ); ii++)
    n_bases += bases[ii].size( );

  for (int ii=0; ii<(int)fixes.size( ); ii++) {
    if ( IsEmpty( fixes[ii].old_base_ ) ) n_bases += 1;
    if ( IsGap( fixes[ii].old_base_ ) ) n_bases += 1;
    if ( IsGap( fixes[ii].new_base0_ ) ) n_bases -= 1;
  }

  int n_objects = (int)bases.size( );

  // Reserve memory.
  vecbasevector bases_out;
  vecqualvector quals_out;
  bases_out.Reserve( n_bases / 16 + n_objects, n_objects );
  quals_out.Reserve( n_bases + n_objects, n_objects );

  // Apply changes.
  int fix_id = 0;
  int last_contig_id = 0;
  int n_contigs = (int)bases.size( );
  while ( fix_id < (int)fixes.size( ) ) {
    int contig_id = fixes[fix_id].id_;

    // Copy unmodified contigs as they are.
    while ( last_contig_id < n_contigs && last_contig_id < contig_id ) {
      bases_out.push_back( bases[last_contig_id] );
      quals_out.push_back( quals[last_contig_id] );
      last_contig_id++;
    }

    // Fixes for this contig.
    int first_id = fix_id;
    int last_id = fix_id + 1;
    while ( last_id < (int)fixes.size( ) && contig_id == fixes[last_id].id_ )
      last_id++;

    // The actual bases and quals that will be modified.
    basevector local_bases = bases[contig_id];
    qualvector local_quals = quals[contig_id];
    vec<String> local_snps( bases[contig_id].size( ), "" );

    // Do corrections.
    for (int ii=last_id-1; ii>=first_id; ii--) {
      const corrector &fix = fixes[ii];
      char old_base = fix.old_base_;
      int old_qual = fix.old_qual_;

      char new_base0 = fix.new_base0_;
      int new_qual0 = fix.new_qual0_;

      char new_base1 = fix.new_base1_;
      int new_qual1 = fix.new_qual1_;

      // Old base is the same as first haplotype of new base (record SNP).
      if ( old_base == new_base0 && old_qual == new_qual0 ) {
	local_snps[fix.pos_] = IsEmpty( new_base1 ) ? "" : fix.NewBase( );
	continue;
      }

      // Insertion.
      if ( IsGap( old_base ) || IsEmpty( old_base ) ) {

	// Resize (fix.pos_ < 0 if bases need to be added at the beginning).
	if ( fix.pos_ < 0 ) {
	  local_bases.resize( local_bases.size( ) + 1 );
	  local_quals.resize( local_quals.size( ) + 1 );
	  local_snps.resize( local_snps.size( ) + 1 );
	}
	else {
	  if ( IsEmpty( old_base ) ) {
	    if ( (int)local_bases.size( ) < fix.pos_ + 1 ) {
	      local_bases.resize( fix.pos_ + 1 );
	      local_quals.resize( fix.pos_ + 1 );
	      local_snps.resize( fix.pos_ + 1 );
	    }
	  }
	  else {
	    local_bases.resize( local_bases.size( ) + 1 );
	    local_quals.resize( local_quals.size( ) + 1 );
	    local_snps.resize( local_snps.size( ) + 1 );
	  }
	}

	// Shift bases after insertion.
	if ( fix.pos_ < 0 || ! IsEmpty( old_base ) ) {
	  for (int jj=local_bases.size( )-2; jj>=Max(0, fix.pos_); jj--) {
	    local_bases.Set( jj + 1, local_bases[jj] );
	    local_quals[jj+1] = local_quals[jj];
	    local_snps[jj+1] = local_snps[jj];
	  }
	}

	// Add base.
	int real_pos = Max( 0, fix.pos_ );
	local_bases.Set( real_pos, as_char( new_base0 ) );
	local_quals[real_pos] = char( new_qual0 );
	local_snps[real_pos] = IsEmpty( new_base1 ) ? "" : fix.NewBase( );

	// Adjust locations.
	for (int jj=first_locs[contig_id]; jj<(int)locs.size( ); jj++) {
	  if ( locs[jj].Contig( ) != contig_id )
	    break;
	  int old_start = locs[jj].StartOnContig( );
	  if ( old_start < fix.pos_ )
	    continue;
	  locs[jj].SetStartOnContig( old_start + 1 );
	}

	continue;
      }

      // Deletion.
      if ( IsGap( new_base0 ) ) {
	for (int jj=fix.pos_; jj<(int)local_bases.size( )-1; jj++) {
	  local_bases.Set( jj, local_bases[jj+1] );
	  local_quals[jj] = local_quals[jj+1];
	  local_snps[jj] = local_snps[jj+1];
	}
	local_bases.resize( local_bases.size( ) - 1 );
	local_quals.resize( local_quals.size( ) - 1 );
	local_snps.resize( local_snps.size( ) - 1 );

	for (int jj=first_locs[contig_id]; jj<(int)locs.size( ); jj++) {
	  if ( locs[jj].Contig( ) != contig_id )
	    break;
	  int old_start = locs[jj].StartOnContig( );
	  if ( old_start < fix.pos_ )
	    continue;
	  locs[jj].SetStartOnContig( old_start - 1 );
	}

	continue;
      }

      // Mutation (or change quality score only).
      local_bases.Set( fix.pos_, as_char( new_base0 ) );
      local_quals[fix.pos_] = char( new_qual0 );
      local_snps[fix.pos_] = IsEmpty( new_base1 ) ? "" : fix.NewBase( );
    }

    // Update (copy fixed contig).
    for (int kk=0; kk<(int)local_snps.size( ); kk++) {
      if ( local_snps[kk] == "" ) continue;
      snpout << last_contig_id << "\t" << kk << "\t" << local_snps[kk] << "\n";
    }

    bases_out.push_back( local_bases );
    quals_out.push_back( local_quals );
    last_contig_id++;

    // Next contig.
    fix_id = last_id;
  }

  // Copy remaining contigs.
  while ( last_contig_id < (int)bases.size( ) ) {
    bases_out.push_back( bases[last_contig_id] );
    quals_out.push_back( quals[last_contig_id] );
    last_contig_id++;
  }

  for (int jj=0; jj<(int)locs.size( ); jj++) {
    int new_len = bases_out[locs[jj].Contig( )].size( );
    locs[jj].SetLengthOfContig( new_len );
  }

  // Swap.
  bases = bases_out;
  quals = quals_out;

  // Cap quality scores.
  for (vecqvec::size_type ii=0; ii<quals.size( ); ii++)
    for (qvec::size_type jj=0; jj<quals[ii].size( ); jj++)
      quals[ii][jj] = static_cast<qual_t>(Min( (int)quals[ii][jj], fin_qual));

  // Flush stream.
  snpout << endl;
}

