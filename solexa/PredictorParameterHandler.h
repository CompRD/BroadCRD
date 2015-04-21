///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PREDICTOR_PARAMETER_HANDLER_H
#define PREDICTOR_PARAMETER_HANDLER_H



/// Interface for getting/putting additional parameters in a PhredTable.
///
/// \class PredictorParameterHandler

struct PredictorParameterHandler {
  virtual ~PredictorParameterHandler() {}
  virtual String GetPredictorParameters() const = 0;
  virtual void SetPredictorParameters(const String & params) = 0;
};

#endif //PREDICTOR_PARAMETER_HANDLER_H
