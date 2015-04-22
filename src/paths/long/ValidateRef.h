///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef VALIDATE_REF
#define VALIDATE_REF

#include "CoreTools.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector.h"

void ValidateRef( const SupportedHyperBasevector& shb, const String& TMP,
     const long_logging_control& log_control, const long_logging& logc );

#endif
