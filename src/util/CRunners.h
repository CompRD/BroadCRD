/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_RUNNERS_H
#define C_RUNNERS_H

#include "String.h"
#include "Vec.h"

/**
 * class CRunners
 *
 * Class to help store informations about modules to run, subdir level
 * and such. The concept is that a core SUBDIR is given in input, and
 * each time a runner (a module to be run) is added, then the module
 * will be run from SUBDIR.level to SUBDIR.(level+1).
 */
class CRunners {

public:

  CRunners( const String &PRE, const String &DATA, const String &RUN );

  // Set initial SUBDIR and final OUTDIR.
  void SetSubdirOutdir( const String &SUBDIR, const String &OUTDIR );

  // Initial setup, reset level_, empty existing data.
  void Setup( const String &SUBDIR, const String &OUTDIR );

  // For example in_type=SUBDIR, or in_type=SUB_DIR (out_type may be empty).
  void SetInOut( const String &in_type = "", const String &out_type = "" );
  
  // Add module to list of runners.
  void Add( const String &module, const String &args = "", bool last = false );
  
  // Get all runners.
  vec<String> GetRunners( ) const;

  // The usual "PRE=... DATA=... RUN=..." thing.
  String PreDataRun( ) const;
  
  // Full run dir.
  String FullRun( ) const;

  // Full sub dir (at this level_, or at OUTDIR if last = true).
  String FullSub( bool last = false ) const;
  
  
private:

  String PRE_;           // Arachne's PRE
  String DATA_;          // DATA
  String RUN_;           // RUN
  String SUBDIR_;        // initial SUBDIR
  String OUTDIR_;        // final OUTDIR

  int level_;            // subdir level
  vec<String> runners_;  // the runners
  
  String in_type_;       // current type for input dir (e.g. SUBDIR, or INDIR)
  String out_type_;      // same, for output
  
};

#endif
