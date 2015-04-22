/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// UnipathCoverage
//
// Given a read set and unipaths created from them, predict the
// number of copies of each unipath in the genome.  Write this prediction to a 
// file reads.unipaths.predicted_count.k*, where * is K.  It is a
// VecPdfEntryVec, where each pdf_entry gives a number of copies and a
// probability for that number, just as if it were a double vector of
// pair<int,double> (but that isn't cross-platform byte compatible).
//
// pdf here refers to "probability density function", although the distribution
// here is discrete ("what is the probability that the copy number is N" so
// this is actually a probability mass function.

#include "Basevector.h"
#include "MainTools.h"
#include "paths/KmerPath.h"
#include "paths/UnipathCoverageCore.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN); 
     CommandArgument_Int(K);
     CommandArgument_String_OrDefault(READS, "reads");
     CommandArgument_Int_OrDefault(UNIPATH_TO_TRACE, -1);
     CommandArgument_Double_OrDefault(THRESH, 0.0001);
     CommandArgument_Double_OrDefault(ERROR_RATE, 0.0);
     CommandArgument_UnsignedInt_OrDefault(USE_THIS_GENOME_SIZE, 0);
     CommandArgument_Bool_OrDefault(GC_BIASED, False);
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Int_OrDefault(PLOIDY, 1);
     CommandArgument_String_OrDefault(CN_SUFFIX, "");
     EndCommandArguments;

     String datadir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;

     vecKmerPath paths( run_dir + "/" + READS + ".paths.k" + ToString(K) );
     vecKmerPath paths_rc( run_dir + "/" + READS + ".paths_rc.k" + ToString(K) );
     BREAD2( run_dir + "/" + READS + ".pathsdb.k" + ToString(K),
	     vec<tagged_rpint>, pathsdb );
     vecKmerPath unipaths( run_dir + "/" + READS + ".unipaths.k" + ToString(K) );
     BREAD2( run_dir + "/" + READS + ".unipathsdb.k" + ToString(K),
          vec<tagged_rpint>, unipathsdb );
     vec<nbases_t> lengths( paths.size( ) );
     
     // Read reads.lengths file if it exists; otherwise calculate lengths
     if ( IsRegularFile( run_dir + "/" + READS + ".lengths" ) )
       BinaryReader::readFile( run_dir + "/" + READS + ".lengths", &lengths );
     else
       for ( size_t i = 0; i < paths.size( ); i++ )
	 lengths[i] = paths[i].KmerCount( ) + K - 1;
     
     VecPdfEntryVec p;
     vec<double> unipath_bias;
     if (GC_BIASED)
     {    READX( run_dir + "/" + READS + ".unipaths.gc_bias.k" + ToString(K), 
               unipath_bias );    }
     UnipathCoverageCore( K, paths, paths_rc, pathsdb, unipaths, unipathsdb, 
          lengths, p, UNIPATH_TO_TRACE, THRESH, ERROR_RATE, USE_THIS_GENOME_SIZE, 
          PLOIDY, ( GC_BIASED ? &unipath_bias : 0 ) );
     if (WRITE)
     {    if ( CN_SUFFIX != "" ) CN_SUFFIX = "." + CN_SUFFIX;
          p.WriteAll( (run_dir + "/" + READS + ".unipaths.predicted_count.k"
               + ToString(K) + CN_SUFFIX).c_str() );    }    }
