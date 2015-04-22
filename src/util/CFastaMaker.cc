/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"
#include "SupersHandler.h"
#include "math/Functions.h"
#include "util/CFastaMaker.h"

/**
 * CFastaMaker
 * Constructor
 */
CFastaMaker::CFastaMaker( const String &supers_file,
			  const vecbasevector &bases, 
			  const vecqualvector &quals ) :
  bases_ ( bases ),
  quals_ ( quals ),
  fout_ ( 0 ),
  qout_ ( 0 ),
  log_ ( 0 )
{
  this->SetDefaults( supers_file );
}

/**
 * CFastaMaker
 * Destructor
 */
CFastaMaker::~CFastaMaker( )
{
  if ( fout_ ) *fout_ << flush;
  if ( qout_ ) *qout_ << flush;
  delete( supers_ );
}

/**
 * CFastaMaker
 * SetGapFloor
 */
void CFastaMaker::SetGapFloor( int min_gap )
{
  min_gap_ = min_gap;
}

/**
 * CFastaMaker
 * SetStreams
 */
void CFastaMaker::SetStreams( ostream *fout, ostream *qout, ostream *log )
{
  fout_ = fout;
  qout_ = qout;
  log_ = log;
}

/**
 * CFastaMaker
 * SuperChunk
 */
void CFastaMaker::SuperChunk( int super_id, int *begin, int *end )
{
  ForceAssert( fout_ );
  ofstream devnull ( "/dev/null" );
  ostream &fout = *fout_;
  ostream &qout = qout_ ? *qout_ : devnull;

  const superb &sup = (*supers_)[super_id];
  int super_len = sup.TrueLength( );
  int n_supers = sup.Ntigs( );

  String descriptor
    =  ">super" + ToString( super_id )
    + "_[" + ToString( begin ? *begin : 0 )
    + "," + ToString( end ? *end : super_len )
    + ")_" + ToString( super_len );
  
  fout << descriptor << "\n";  
  qout << descriptor << "\n";

  for (int cgpos=0; cgpos<n_supers; cgpos++) {
    this->SetToPrint( sup.Tig( cgpos ), begin, end );
    this->PrintFasta( );
    this->PrintQual( );
  }
  this->EndFasta( );
  this->EndQual( );
}

/**
 * CFastaMaker
 * SetDefaults
 */
void CFastaMaker::SetDefaults( const String &supers_file )
{
  fwidth_ = 80;
  qwidth_ = 25;
  min_gap_ = 100;

  fpos_ = 1;
  qpos_ = 1;
  contig_id_ = -1;
  contig_beg_ = -1;
  contig_end_ = -1;
  contig_Ns_ = 0;

  supers_ = new shandler( bases_.size( ) );
  supers_->LoadFromFile( supers_file, &min_gap_ );
}

/**
 * CFastaMaker
 * SetToPrint
 *
 * Find details on the window to be printed.
 */
void CFastaMaker::SetToPrint( int contig_id, int *begin, int *end )
{
  int super_id = supers_->ToSuper( contig_id );
  int contig_pos = supers_->PosOnSuper( contig_id );
  int contig_start = supers_->StartOnSuper( contig_id );
  int contig_len = (int)bases_[contig_id].size( );
  int sbegin = begin ? *begin : 0;
  int send = end ? *end : (*supers_)[super_id].TrueLength( );
  
  contig_id_ = contig_id;
  contig_beg_ = Min( Max( 0, sbegin - contig_start ), contig_len );
  contig_end_ = Max( 0, Min( contig_len, send - contig_start ) );
  
  bool is_last = ( contig_pos == (*supers_)[super_id].Ntigs( ) - 1 );
  int gap_start = contig_start + contig_len;
  int gap_len = is_last ? 0 : (*supers_)[super_id].Gap( contig_pos );
  int gap_beg = Min( Max( 0, sbegin - gap_start ), gap_len );
  int gap_end = Max( 0, Min( gap_len, send - gap_start ) );

  contig_Ns_ = Max( 0, gap_end - gap_beg );

  // Log event.
  if ( ! log_ )
    return;
  
  if ( contig_beg_ != contig_end_ )
    *log_ << "c_" << contig_id
	  << "=s" << super_id
	  << "_" << contig_pos
	  << "." << (*supers_)[super_id].Ntigs( ) - 1
	  << "\t[" << contig_beg_
	  << "," << contig_end_
	  << ")_" << contig_len
	  << "\n";

  if ( contig_Ns_ > 0 )
    *log_ << "gap s" << super_id
	  << "_" << contig_pos
	  << "." << (*supers_)[super_id].Ntigs( ) - 1
	  << "\t[" << gap_beg
	  << "," << gap_end
	  << ")_" << gap_len
	  << "\n";
  
}

/**
 * CFastaMaker
 * PrintFasta
 *
 * Print bases in contig_id_ as specified by contig_*.
 */
void CFastaMaker::PrintFasta( )
{
  ForceAssert( fout_ );
  ostream &fout = *fout_;

  const basevector &base = bases_[contig_id_];
  for (int pos=contig_beg_; pos<contig_end_; pos++) {
    if ( fpos_ > fwidth_ ) {
      fout << "\n";
      fpos_ = 1;
    }
    fout << as_base( base[pos] );
    fpos_++;
  }
  
  for (int ii=0; ii<contig_Ns_; ii++) {
    if ( fpos_ > fwidth_ ) {
      fout << "\n";
      fpos_ = 1;
    }
    fout << "N";
    fpos_++;
  }
}

/**
 * CFastaMaker
 * PrintQual
 *
 * Print quals in contig_id_ as specified by contig_*.
 */
void CFastaMaker::PrintQual( )
{
  if ( !qout_ ) return;
  ostream &qout = *qout_;

  const qualvector &qual = quals_[contig_id_];
  for (int pos=contig_beg_; pos<contig_end_; pos++) {
    if ( qpos_ > qwidth_ ) {
      qout << "\n";
      qpos_ = 1;
    }
    qout << (int)qual[pos] << " ";
    qpos_++;
  }
  
  for (int ii=0; ii<contig_Ns_; ii++) {
    if ( qpos_ > qwidth_ ) {
      qout << "\n";
      qpos_ = 1;
    }
    qout << "0 ";
    qpos_++;
  }
}

/**
 * CFastaMaker
 * EndFasta
 */
void CFastaMaker::EndFasta( )
{
  ForceAssert( fout_ );
  *fout_ << ( fpos_ == 1 ? "" : "\n" );
}

/**
 * CFastaMaker
 * EndQual
 */
void CFastaMaker::EndQual( )
{
  if ( !qout_ ) return;
  *qout_ << ( qpos_ == 1 ? "" : "\n" );
}

