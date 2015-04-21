///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "VecUtilities.h"
#include "paths/CMasterBase.h"

/**
 * classCMasterBase
 * Constructor
 */
CMasterBase::CMasterBase( )
{
  this->ClearAll( );
  this->SetDefaults( );
}

/**
 * classCMasterBase
 * ClearAll
 */
void CMasterBase::ClearAll( )
{
  idsA_.clear( );
  idsC_.clear( );
  idsG_.clear( );
  idsT_.clear( );
  idsDels_.clear( );
  idsIns_.clear( );
  basesIns_.clear( );
}

/**
 * class CMasterBase
 * SetMinCovs
 */
void CMasterBase::SetMinCovs( pair<double,double> min_covs )
{
  min_covs_ = min_covs;
  min_coverage_ = -1;
}

/**
 * class CMasterBase
 * AddBase
 */
void CMasterBase::AddBase( const size_t &id, const char base )
{
  if ( base == 'A' ) idsA_.push_back( id );
  else if ( base == 'C' ) idsC_.push_back( id );
  else if ( base == 'G' ) idsG_.push_back( id );
  else if ( base == 'T' ) idsT_.push_back( id );
  else ForceAssert( 1 == 0 );
}

/**
 * class CMasterBase
 * AddInsertion
 */
void CMasterBase::AddInsertion( const size_t &id, const bvec &ins )
{
  int insertion_id = -1;
  for (int ii=0; ii<(int)basesIns_.size( ); ii++) {
    if ( basesIns_[ii] == ins ) {
      insertion_id = ii;
      break;
    }
  }
  if ( insertion_id < 0 ) {
    vec<size_t> vec_id( 1, id );
    basesIns_.push_back( ins );
    idsIns_.push_back( vec_id );
  }
  else
    idsIns_[insertion_id].push_back( id );
}

/**
 * class CMasterBase
 * AddInsertion
 */
void CMasterBase::AddDeletion( const size_t &id )
{
  idsDels_.push_back( id );
}

/**
 * class CMasterBase
 * InsCoverage
 */
size_t CMasterBase::InsCoverage( ) const
{
  size_t n_ins = 0;
  for (int ii=0; ii<idsIns_.isize( ); ii++)
    n_ins += idsIns_[ii].size( );
  return n_ins;
}

/**
 * class CMasterBase
 * DelCoverage
 */
size_t CMasterBase::DelCoverage( ) const
{
  return idsDels_.size( );
}

/**
 * class CMasterBase
 * BaseCoverage
 */
size_t CMasterBase::BaseCoverage( ) const
{
  return idsA_.size( ) + idsC_.size( ) + idsG_.size( ) + idsT_.size( );
}

/**
 * class CMasterBase
 * TotCoverage
 */
size_t CMasterBase::TotCoverage( ) const
{
  return this->InsCoverage( ) + this->DelCoverage( ) + this->BaseCoverage( );
}

/**
 * class CMasterBase
 * MostCommonBase
 */
char CMasterBase::MostCommonBase( ) const
{
  vec< pair<int,int> > cov2base
    = MkVec( make_pair( idsA_.isize( ), 0 ),
	     make_pair( idsC_.isize( ), 1 ),
	     make_pair( idsG_.isize( ), 2 ),
	     make_pair( idsT_.isize( ), 3 ) );
  sort( cov2base.rbegin( ), cov2base.rend( ) );
  return cov2base[0].second;
}

/**
 * class CMasterBase
 * AddInsertion
 */
void CMasterBase::PrintBrief( const String &info, ostream &out ) const
{
  const size_t tot_cov = this->TotCoverage( );
  const double min_cov = this->MinCoverage( );
  
  out << info << "  ";
  if ( tot_cov < 1 ) {
    out << "no bases\n";
    return;
  }
  out << "cov: " << tot_cov << "  n_hap: " << this->NHaplotypes( ) << "  ";
    
  double ratioA = SafeQuotient( idsA_.size( ), tot_cov );
  double ratioC = SafeQuotient( idsC_.size( ), tot_cov );
  double ratioG = SafeQuotient( idsG_.size( ), tot_cov );
  double ratioT = SafeQuotient( idsT_.size( ), tot_cov );
  double ratioD = SafeQuotient( this->DelCoverage( ), tot_cov );

  // Show haplotypes.
  if ( ratioA > min_cov )
    out << "a[" << ToString( 100.0 * ratioA, 2 ) << "%]  ";

  if ( ratioC > min_cov )
    out << "c[" << ToString( 100.0 * ratioC, 2 ) << "%]  ";
  
  if ( ratioG > min_cov )
    out << "g[" << ToString( 100.0 * ratioG, 2 ) << "%]  ";

  if ( ratioT > min_cov )
    out << "t[" << ToString( 100.0 * ratioT, 2 ) << "%]  ";

  if ( ratioD > min_cov )
    out << "del[" << ToString( 100.0 * ratioD, 2 ) << "%]  ";

  for (int ii=0; ii<idsIns_.isize( ); ii++) {
    double ratioI = SafeQuotient( idsIns_[ii].size( ), tot_cov );
    if ( ratioI > min_cov )
      out << "ins." << basesIns_[ii].ToString( )
	  << ".[" << ToString( 100.0 * ratioI, 2 ) << "%]  ";
  }
  
  // Done.
  out << "\n";
}

/**
 * class CMasterBase
 * ListHaplotypes
 *
 * Send to stream a one liner with all the variants for this event,
 * with the exception of "base" (which is assumed to be one of A, C,
 * G, or T).
 */
void CMasterBase::ListHaplotypes( const String &info,
				  const char base,
				  ostream &log ) const
{
  const double min_cov = this->MinCoverage( );
  const size_t tot_cov = this->TotCoverage( );
  if ( tot_cov < 1 ) {   // nothing found, empty event.
    log << info << "  na\n";
    return;
  }
  
  pair<int,String> mainVar;
  if ( 'A' == base ) mainVar = make_pair( idsA_.isize( ), "A" );
  else if ( 'C' == base ) mainVar = make_pair( idsC_.isize( ), "C" );
  else if ( 'G' == base ) mainVar = make_pair( idsG_.isize( ), "G" );
  else if ( 'T' == base ) mainVar = make_pair( idsT_.isize( ), "T" );
  else ForceAssert( 1 == 0 );

  vec< pair<int,String> > cov2var;
  if ( SafeQuotient( idsA_.size( ), tot_cov ) > min_cov )
    if ( 'A' != base ) cov2var.push_back( make_pair( idsA_.isize( ), "A" ) );

  if ( SafeQuotient( idsC_.size( ), tot_cov ) > min_cov )
    if ( 'C' != base ) cov2var.push_back( make_pair( idsC_.isize( ), "C" ) );
  
  if ( SafeQuotient( idsG_.size( ), tot_cov ) > min_cov )
    if ( 'G' != base ) cov2var.push_back( make_pair( idsG_.isize( ), "G" ) );

  if ( SafeQuotient( idsT_.size( ), tot_cov ) > min_cov )
    if ( 'T' != base ) cov2var.push_back( make_pair( idsT_.isize( ), "T" ) );
  
  if ( SafeQuotient( idsDels_.isize( ), tot_cov ) > min_cov )
    cov2var.push_back( make_pair( idsDels_.isize( ), "*" ) );

  for (int ii=0; ii<idsIns_.isize( ); ii++) {
    if ( SafeQuotient( idsIns_[ii].size( ), tot_cov ) > min_cov ) {
      String str_ins = "i." + basesIns_[ii].ToString( );
      cov2var.push_back( make_pair( idsIns_[ii].isize( ), str_ins ) );
    }
  }

  // Sort these, so variants with higher coverage appear first.
  sort( cov2var.rbegin( ), cov2var.rend( ) );
  
  // Print variants.
  log << info << "  [" << mainVar.second << " " << mainVar.first << "]";
  for (int ii=0; ii<cov2var.isize( ); ii++)
    log << "  (" << cov2var[ii].second << " " << cov2var[ii].first << ")";
  log << "\n";
  
}

/**
 * class CMasterBase
 * NHaplotypes
 */
int CMasterBase::NHaplotypes( ) const
{
  const double min_cov = this->MinCoverage( );
  const size_t tot_cov = this->TotCoverage( );
  if ( tot_cov < 1 ) return 0;

  int n_hap = 0;
  if ( SafeQuotient( idsA_.size( ), tot_cov ) > min_cov ) n_hap++;
  if ( SafeQuotient( idsC_.size( ), tot_cov ) > min_cov ) n_hap++;
  if ( SafeQuotient( idsG_.size( ), tot_cov ) > min_cov ) n_hap++;
  if ( SafeQuotient( idsT_.size( ), tot_cov ) > min_cov ) n_hap++;
  if ( SafeQuotient( idsDels_.isize( ), tot_cov ) > min_cov ) n_hap++;
  for (int ii=0; ii<idsIns_.isize( ); ii++)
    if ( SafeQuotient( idsIns_[ii].size( ), tot_cov ) > min_cov ) n_hap++;
    
  return n_hap;
}

/**
 * class CMasterBase
 * SetDefaults
 * private
 */
void CMasterBase::SetDefaults( )
{
  // Remove noise - Do not consider an haplotype if the ratio between
  //  its coverage and the total base coverage is smaller than a given
  //  threshold: this threshold depends on the coverage of the base,
  //  but is always contained in the interval specified by min_covs_.
  //  See also MinCoverage( ) below.

  min_covs_ = make_pair( 0.02, 0.06 );
  min_coverage_ = -1;
  
}

/**
 * class CMasterBase
 * MinCoverage
 * private
 *
 * min_coverage is interpolated linearly as:
 *
 *   min_covs_.first        if coverage <= 2000
 *   a line                 2000 <= coverage <= 8000
 *   min_covs_.second       if coverage >= 8000
 *
 */
double CMasterBase::MinCoverage( ) const
{
  if ( min_coverage_ < 0 ) {
    const double a = 2000.0;
    const double b = 6000.0;
    const double m = min_covs_.first;
    const double M = min_covs_.second;
    double cov = this->TotCoverage( );

    if ( cov < a )      min_coverage_ = m;
    else if ( cov < b ) min_coverage_ = ( ( M-m )*cov + m*b - M*a ) / ( b-a );
    else                min_coverage_ = M;
  }
  
  return min_coverage_;
}
