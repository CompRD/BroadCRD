///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Oct 5, 2012 
//

#include "simulation/ReadSimulatorSimpleCore.h"
#include "Basevector.h"


struct ReferenceTools {
  typedef enum _FileType { FASTB, FASTA } FileType;
  static void fromFile( BaseVecVec* bvv_p, const String& name, FileType type, int ref_id = -1, int ref_start = -1, int ref_stop = -1 ) {
    if (type == FASTB)
      reference_get(bvv_p, name, "", ref_id, ref_start, ref_stop, 0, "", 0);
    else {
      ForceAssertEq(static_cast<int>(type), static_cast<int>(FASTA));  // TODO: bvv_p check is whack... sorry.
      reference_get(bvv_p, "", name, ref_id, ref_start, ref_stop, 0, "", 0);
    }
  }
  static void fromFastb( BaseVecVec* bvv_p, const String& name, int ref_id = -1, int ref_start = -1, int ref_stop = -1 ) {
    fromFile( bvv_p, name, FASTB, ref_id, ref_start, ref_stop );
  }
  static void fromFasta( BaseVecVec* bvv_p, const String& name, int ref_id = -1, int ref_start = -1, int ref_stop = -1 ) {
    fromFile( bvv_p, name, FASTA, ref_id, ref_start, ref_stop);
  }
  static void fromRandom( BaseVecVec* bvv_p, const size_t g_size, const unsigned seed = 999 ) {
    reference_get(bvv_p, "", "", -1, -1, -1, g_size, "", seed );
  }

  static void toFastb( BaseVecVec* bvv_p, const String& head_name ) {
    bvv_p->WriteAll( head_name + ".fastb" );
  }
private:
  ReferenceTools();	// don't actually try to instantiate this
};


