///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SOLEXA_TOOLS_H
#define SOLEXA_TOOLS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "FetchReads.h"
#include "lookup/LookAlign.h"
#include "random/Shuffle.h"
#include "solexa/FourBase.h"
#include "solexa/PredictorParameterHandler.h"
#include "Map.h"
#include "solexa/SolexaPipeline.h"
#include "Boolvector.h"

/** Parameters needed for producing error predictors.
\class SolexaPredictorParameters
*/
class SolexaPredictorParameters: public PredictorParameterHandler {
  //parameters used in SffRead::GetPredictors:
private:
  typedef map<String,String> smap;
  smap params_;

public:
  virtual ~SolexaPredictorParameters() {}
  virtual String GetPredictorParameters() const;
  virtual void SetPredictorParameters(const String & params);

  int Get(const String & s) const {
    if (!IsKey(params_,s)) return -1;
    else return params_.find(s)->second.Int();
  }

  String GetString(const String & s) const {
    if (!IsKey(params_,s)) return "";
    else return params_.find(s)->second;
  }

  void SetFromFile(const String & fname);
};


/// Function: LoadSolexaData
///
/// Load Solexa data from location HEAD, generating:
///   bases - the reads
///   basesrc - the reverse complements of the reads
///   I - the intensities
///   ref - the reference
///   rref - the reverse complement of the reference
///   aligns - alignments of the reads to the reference
///   aligns_index - index by read of the alignments
void LoadSolexaData( const String& HEAD, vecbasevector& bases,
     vecbasevector& basesrc, VecFourBaseVec& I, vecbasevector& ref,
     vecbasevector& rref, vec<look_align>& aligns, vec< vec<int> >& aligns_index, Bool loadIntensities = True );

/// Load only bases and intensities.

void LoadSolexaData( const String& HEAD, vecbasevector& bases,
     vecbasevector& basesrc, VecFourBaseVec& I, Bool loadIntensities = True );

/*
   Function: ExpandHead

   Given a file prefix representing a <string set> -- most commonly based on
   a set of <lane> numbers -- return a list of strings of individual prefixes,
   one per lane.
*/
void ExpandHead( const String& HEAD, vec<String>& HEADS, int trunc = 0,
     const String& pipeline = SOLEXAPIPELINE, Bool paired = False );

void RequireSameRef( const vec<String>& HEADS );

void SetRef( const String& REF, String& REFX, const String& HEAD,
     const vec<String>& HEADS );

/// Load from multiple heads.  The argument HEADS is in ParseStringSet format.
/// See ExpandHead in SolexaTools.cc for how the list is expanded.
///
/// trunc: use subdirectory trunc_n (only works for certain forms)

void LoadSolexaData( const String& HEADS, const String& QLTOUT,
     const String& REF, vecbasevector& bases, VecFourBaseVec& I,
     vecbasevector& ref, vecbasevector& rref, vec<look_align>& aligns,
     vec< vec<int> >& aligns_index, Bool loadIntensities = True );

void LoadSolexaData( const String& HEADS, const String& HEADS_TAIL, 
     const int trunc, const String& REF, vecbasevector& bases, VecFourBaseVec& I,
     vecbasevector& ref, Bool loadIntensities = True );

void LoadSolexaData( const String& HEAD, const String& HEAD_TAIL, const int trunc,
     const String& REF, vecbasevector& bases, vecbasevector& ref );

void LoadSolexaData( const String& HEADS, vecbasevector& bases, 
        VecFourBaseVec& I, Bool loadIntensities = True );

void LoadSolexaData( const String& HEADS, const String& REF, 
     vecbasevector& ref, vecbasevector& rref );

/// FilterPileups: filter out read pileups: coverage > 2 * mean coverage.  Return
// number of reads filtered.  Changes best entry to -1 to signify filtering out.

int FilterPileups(
     const vecbasevector& bases,             // the reads
     const vecbasevector& ref,               // the reference
     const vec<look_align>& aligns,          // alignments of reads to reference
     const vec< vec<int> >& aligns_index,    // index by read of alignments
     const vec<float>& minq10,               // min quality of first 10 bases
     vec<int>& best,                         // index of best alignment (modified)
     const double MIN_RATIO = 2.0,           // defines filtering
     Bool verbose = False
          );

/// GetMinQ10: for each read, find the minimum quality of its first 10 bases.
/// Templatized to allow use of serial or mapped feudal access.
/// dummyCall = True - do not really check call quality.  Used when did not
/// load intensities.
// TODO: potentially dangerous truncation of index by nreads
template< class MasterVecLike >
void GetMinQ10( const MasterVecLike & I, int nreads, vec<float>& minq10, Bool dummyCall = False )
{    minq10.resize(nreads);
     for ( int id = 0; id < nreads; id++ ) 
     {    minq10[id] = 1000000000;
          if (dummyCall) continue;
          FourBaseVec const& fbv = I[id];
          if ( fbv.size( ) < 10 ) minq10[id] = 1.0;
          else
          {    for ( int j = 0; j < 10; j++ )
               {    minq10[id] = Min( 
                         minq10[id], fbv[j].CallQuality( ) );
		}    
	  }    
     }    
}

/// Select reads that are not homopolymers and measure decay/height.
/// Takes in a filename and uses ReadOne, since it picks reads at
/// random from the file.
void GetDecayHeight( const String & intFname, const vecbasevector & reads,
		     NormalDistribution & decay, NormalDistribution & height,
		     const int N=5000);

/// Function: FindRun
///
/// Given a <flowcell>, find the corresponding <run> in the <Broad pipeline directory>.
/// Return its path as _run_.  Also return the date of the run.
void FindRun( const String& flowcell, String& run, String& run_date,
     const String& pipeline = SOLEXAPIPELINE );

/// Function: FindPairedRun
///
/// Given a <flowcell>, find the corresponding paired <run>s in the <Broad pipeline directory>.
/// Return its path as _run_.
void FindPairedRun( const String& flowcell, String& runA, String& runB,
		    const String& pipeline );

/// Read the indicated column of paramsbyread file into vec v.
void LoadVecFromParamsFile( const String & file, int column, vec<float> & v);

/// Read the indicated column of paramsbyread file into vec v.
void LoadVecFromParamsFile( const String & file, const String &  columnName, 
			    vec<float> & v);

#define VEC_FROM_PARAMS_FILE(fname, vecname) \
vec<float> vecname; \
LoadVecFromParamsFile(fname, #vecname, vecname);


/**
   Function: LoadReadAligns

   For a set of reads in a <.fastb> file, load their alignments to the reference
   from a <.qltout> file.

   Parameters:

      qltoutFileName - the <.qltout> file from which to load the alignments.
*/
void LoadReadAligns( const String& qltoutFileName, read_aligns_plus_t& readAligns );

// Creates a BaitMap to do this, so it you already have one created
// use the BaitMap method instead.
void LoadBaitMapAsMask(vec<boolvector> *mask, const String &file_name, int padding=0);

   
#endif
