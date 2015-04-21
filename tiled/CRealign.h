/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_REALIGN_H
#define C_REALIGN_H

#include "Basevector.h"
#include "LocsHandler.h"
#include "pairwise_aligners/KmerAligner.h"
#include "tiled/TAlign.h"

/**
 * class CRealign
 *
 * Realign reads to consensus.
 */
class CRealign {

public:
  
  CRealign( const vecbasevector &cbases,
	    const vecbasevector &rbases,
	    const lhandler &locs,
	    ostream *log = 0,
	    ostream *dbfastalog = 0 );
  
  ~CRealign( ) { if ( plog_ ) plog_->flush( ); }

  void SetMinMaxBand( int min, int max );

  void SetMaxDiscrep ( int discrep ) { max_discrep_ = discrep; }

  void SetMaxError( float error ) { max_error_ = error; }

  void SetVerbose( bool verbose ) { verbose_ = verbose; }
  
  // Set contig and create a (default) KmerAligner.
  void SetContig( int contig_id );

  // Try to align the given read to consensus (return false on failure).
  bool PlaceLoc( int locpos, alignment &al );

  // Same as above, but fill a t_align instead.
  bool PlaceLoc( int locpos, t_align &tal );

  // Set locpos_ and find seeds.
  void FindSeeds( int locpos );
  
  // Send to log informations about the seeds.
  void ReportSeeds( ) const;
  
  
private:

  void SetDefaults( );

  // Dump fasta of failed alignment onto dbfastalog_ (if defined).
  void DumpFasta( );

  // Determine offset and bandwidth. Return false if not found.
  bool GetOffset( int &offset, int &band ) const;

  // Return the seed closest to the locs initial guess (-1 if not found).
  int ClosestSeed( ) const;

  // Returns ( placement_by_seed - original_placement ).
  int Discrepancy( int seedpos ) const;

  // Set al to be the "bad" one-base alignment (at the end of contig).
  void SetToOneBaseAlign( int locpos, alignment &al ) const;

  // Return false on failure.
  bool RunAlign2Bvecs( alignment &al, bool K8 = false ) const;
  
  // Decide if given alignment is valid.
  bool IsValid( const alignment &al ) const;
  
  
private:
  
  const vecbasevector &cbases_;  // contigs
  const vecbasevector &rbases_;  // reads
  const lhandler &locs_;         // read locations
  ostream *plog_;                // log stream pointer (may be null)
  ostream *dbfastalog_;          // debug log (outputs fasta, may be null)

  int min_band_;                 // used when chopping contig to localize align
  int max_band_;                 // used when chopping contig to localize align
  int max_discrep_;              // max placement discrepancy (if >=0)
  float max_error_;              // max error rate (if <1.0)
  bool verbose_;                 // verbose log

  int contig_id_;                // contig id
  bool set12_;                   // flag for the setup of aligner12_
  KmerAligner<12> aligner12_;    // fine grain seeding tool
  KmerAligner<24> aligner24_;    // large grain (default) seeding tool
  
  int locpos_;                   // selected read (as pos in locs_)
  basevector the_read_;          // bases of selected read (rc-ed if needed)
  basevector the_contig_;        // bases of chopped contig
  pair<int,int> cgwin_;          // window defining how contig was chopped
  vec< pair<int,int> > seeds_;   // seeds/weights for selected read (sorted)
  
};

#endif

