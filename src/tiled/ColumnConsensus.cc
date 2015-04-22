/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Vec.h"
#include "tiled/ColumnConsensus.h"
#include "tiled/Haploqual.h"
#include "tiled/PaddedSeq.h"
#include "tiled/CSnp.h"

/**
 * column_consensus
 * Constructor
 */
column_consensus::column_consensus( ) :
  A_ ( 0 ),
  C_ ( 0 ),
  G_ ( 0 ),
  T_ ( 0 ),
  pad_ ( 0 )
{
  this->SetDefaults( );
}

/**
 * column_consensus
 * Destructor
 */
column_consensus::~column_consensus( )
{
  if ( A_ ) delete ( A_ );
  if ( C_ ) delete ( C_ );
  if ( G_ ) delete ( G_ );
  if ( T_ ) delete ( T_ );
  if ( pad_ ) delete ( pad_ );
}

/**
 * column_consensus
 * SetDefaults
 */
void column_consensus::SetDefaults( )
{
  good_qual_ = 20;
  max_to_exclude_ = 20;
  min_to_include_ = 40;
}

/**
 * column_consensus
 * Add
 */
void column_consensus::Add( char base, int qual )
{
  // Switch to upper case.
  ToUpperCase( base );
  
  if ( 'A' == base ) {
    if ( !A_ ) A_ = new haploqual( 'A' );
    A_->AddQual( qual );
  }
  else if ( 'C' == base ) {
    if ( !C_ ) C_ = new haploqual( 'C' );
    C_->AddQual( qual );
  }
  else if ( 'G' == base ) {
    if ( !G_ ) G_ = new haploqual ( 'G' );
    G_->AddQual( qual );
  }
  else if ( 'T' == base ) {
    if ( !T_ ) T_ = new haploqual ( 'T' );
    T_->AddQual( qual);
  }
  else if ( gap_base == base ) {
    if ( !pad_ ) pad_ = new haploqual( gap_base );
    pad_->AddQual( qual );
  }
  else {
    // Only A, a, C, c, G, g, T, t, or gap_base are valid bases.
    ForceAssert( 1 == 0 );
  }
}

/**
 * column_consensus
 * Coverage
 */
int column_consensus::Coverage( ) const
{
  int cov = 0;
  if ( A_ ) cov += A_->Coverage( );
  if ( C_ ) cov += C_->Coverage( );
  if ( G_ ) cov += G_->Coverage( );
  if ( T_ ) cov += T_->Coverage( );
  if ( pad_ ) cov += pad_->Coverage( );
  
  return cov;
}

/**
 * column_consensus
 * Coverage
 */
int column_consensus::Coverage( char base ) const
{
  // Switch to upper case.
  ToUpperCase( base );
  
  if ( 'A' == base )
    if ( A_ ) return A_->Coverage( );
    else return 0;
  else if ( 'C' == base )
    if ( C_ ) return C_->Coverage( );
    else return 0;
  else if ( 'G' == base )
    if ( G_ ) return G_->Coverage( );
    else return 0;
  else if ( 'T' == base )
    if ( T_ ) return T_->Coverage( );
    else return 0;
  else if ( gap_base == base )
    if ( pad_ ) return pad_->Coverage( );
    else return 0;
  
  // Not a valid base.
  ForceAssert( 1 == 0 );
  return 0;
}

/**
 * column_consensus
 * Consensus
 *
 * If tag_snp is true, then SNPs are tagged as such, and quality score will
 * be set to snp_qual. Otherwise, consensus base will be one of the two
 * haplotypes, and quality score will be set to untrusted_qual.
 */
bool column_consensus::Consensus( char &base, int &qual, bool tag_snp ) const
{
  // Fall back default values.
  bool consensus_not_found = false;
  bool consensus_found = true;
  base = empty_base;
  qual = empty_qual;
  
  // Create and sort (descending) a local vector of non-null haploquals.
  vec< haploqual* > haplo;
  this->PackageSort( haplo );
  order_haploqual_Consensus sorter;
  sort( haplo.rbegin( ), haplo.rend( ), sorter );
  
  // No haplotypes (should never be the case).
  if ( haplo.size( ) == 0 )
    return consensus_not_found;
  
  // One haplotype only.
  if ( haplo.size( ) == 1 ) {
    base = haplo[0]->Base( );
    qual = haplo[0]->ConsensusQual( );
    return consensus_found;
  }
  
  // More than two haplotypes, with high quality scores. Notice that
  //  return value is consensus_found, since we know from the quality
  //  score that this is an untrusted base.
  if ( (int)haplo.size( ) > 2 && haplo[2]->ConsensusQual( ) >= good_qual_ ) {
    base = haplo[0]->Base( );
    qual = untrusted_qual;
    return consensus_found;
  }
  
  // Two or more haplotypes: an SNP.
  if ( this->IsSNP( *haplo[0], *haplo[1] ) ) {
    if ( tag_snp ) {
      base = this->BaseSNP( haplo[0]->Base( ), haplo[1]->Base( ) );
      qual = snp_qual;
    }
    else {
      base = haplo[0]->Base( );
      qual = untrusted_qual;
    }
    return consensus_found;
  }

  // Two or more haplotypes: only the first is significant.
  base = haplo[0]->Base( );
  qual = haplo[0]->ConsensusQual( );
  return consensus_found;

}

/**
 * column_consensus
 * SnpConsensus
 *
 * Generate SNP consensus. If an inconsistency is found, return false.
 * A true SNP is returned iff no consistencies are found.
 *
 * nhap: either 1 or 2 (haplotypes allowed)
 * trusted: tag for untrusted regions
 */
bool column_consensus::SnpConsensus( CSnp &snp, int nhap, bool trusted ) const
{
  snp.Set( empty_base, empty_qual );
  ForceAssert( nhap == 1 || nhap == 2 );

  // Scale heuristics if trusted = false.
  int multiplier = trusted ? 1 : 2;
  int max_exclude = max_to_exclude_ / multiplier;
  int min_include = min_to_include_ * multiplier;

  // Sort haploquals.
  vec< haploqual* > haplo;
  this->PackageSort( haplo );
  
  // Set base and qual to the best (highest qual) value.
  if ( haplo.size( ) < 1 )
    return false;
  char base1 = (haplo[0])->Base( );
  int qual1 = (haplo[0])->SumQual( );
  snp.Set( base1, qual1 );
  
  // There is only one haplotype.
  if ( haplo.size( ) == 1 )
    return true;
  
  // Only one haplotype allowed.
  if ( nhap == 1 ) {
    snp.Set( base1, 1 );
    return ( ( haplo[1] )->SumQual( ) <= max_exclude );
  }

  // Two haplotypes allowed.
  if ( haplo.size( ) > 2 && ( haplo[2] )->SumQual( ) > max_exclude ) {
    snp.Set( base1, 1 );
    return false;
  }
  
  // From here on we only look at the first two haplotypes.
  char base2 = (haplo[1])->Base( );
  int qual2 = (haplo[1])->SumQual( );
  if ( qual2 <= max_exclude ) {
    snp.Set( base1, qual1 );
    return true;
  }
  if ( qual2 < min_include ) {
    snp.Set( base1, 1 );
    return false;
  }
  
  // A true SNP.
  snp.Set( base1, qual1, base2, qual2 );
  return true;
}

/**
 * column_consensus
 * PrintInfo
 */
void column_consensus::PrintInfo( ostream &out, bool newline ) const
{
  if ( ! ( A_ || C_ || G_ || T_ || pad_ ) ) out << "empty_col  ";
  if ( A_ ) out << "A." << A_->SumQual( ) << "  ";
  if ( C_ ) out << "C." << C_->SumQual( ) << "  ";
  if ( G_ ) out << "G." << G_->SumQual( ) << "  ";
  if ( T_ ) out << "T." << T_->SumQual( ) << "  ";
  if ( pad_ ) out << gap_base << "." << pad_->SumQual( ) << "  ";
  out << ( newline ? "\n" : "" );
}

/**
 * column_consensus
 * IsSNP
 */
bool column_consensus::IsSNP( const haploqual &haplo1,
			      const haploqual &haplo2 ) const
{
  int qual1 = haplo1.ConsensusQual( );
  int qual2 = haplo2.ConsensusQual( );
  ForceAssert( qual1 >= qual2 );
  
  // haplo2 has poor quality score.
  if ( qual2 < good_qual_ )
    return false;
  
  // Regardless of quality score, haplo1 appears many more times.
  if ( haplo1.Coverage( ) > 3 && haplo2.Coverage( ) < 2 )
    return false;
  
  // Probably an SNP.
  return true;
}

/**
 * column_consensus
 * BaseSNP
 */
char column_consensus::BaseSNP( char base1, char base2 ) const
{
  // Switch to upper case.
  ToUpperCase( base1 );
  ToUpperCase( base2 );
  
  // Get to base1 < base2 (where pad is always last).
  ForceAssert( base1 != base2 );
  if ( base1 == gap_base )
    swap( base1, base2 );
  if ( base2 != gap_base && base2 < base1 )
    swap( base1, base2 );
  
  // Find proper code.
  if ( base1 == 'A' ) {
    if ( base2 == 'C' ) return snp_AC_base;
    else if ( base2 == 'G' ) return snp_AG_base;
    else if ( base2 == 'T' ) return snp_AT_base;
    else if ( base2 == gap_base ) return snp_Apad_base;
  }
  else if ( base1 == 'C' ) {
    if ( base2 == 'G' ) return snp_CG_base;
    else if ( base2 == 'T' ) return snp_CT_base;
    else if ( base2 == gap_base ) return snp_Cpad_base;
  }
  else if ( base1 == 'G' ) {
    if ( base2 == 'T' ) return snp_GT_base;
    else if ( base2 == gap_base ) return snp_Gpad_base;
  }
  else if ( base1 == 'T' ) {
    if ( base2 == gap_base ) return snp_Tpad_base;
  }
  
  // We should never get here.
  ForceAssert( 1 == 0 );
  return empty_base;
}

/**
 * column_consensus
 * PackageSort
 *
 * Put haploquals in a sorted vector, so 
 */
void column_consensus::PackageSort( vec< haploqual* > &haplo ) const
{
  haplo.clear( );
  if ( A_ ) haplo.push_back( A_ );
  if ( C_ ) haplo.push_back( C_ );
  if ( G_ ) haplo.push_back( G_ );
  if ( T_ ) haplo.push_back( T_ );
  if ( pad_ ) haplo.push_back( pad_ );
  
  order_haploqual_SumQual sorter;
  sort( haplo.rbegin( ), haplo.rend( ), sorter );
}


