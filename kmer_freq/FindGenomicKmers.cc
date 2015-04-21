/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/*
   Program: FindGenomicKmers

   Given a reference genome, for a given <kmer shape>, finds the frequency of each kmer of that shape
   in the reference genome.

   Program parameters:

      PRE - the <WGA data dir>
      DATA - the data directory for our ALLPATHS project
      Ks - the <kmer shapes> to use; everything is done for each kmer shape in turn.
      REBUILD -  if REBUILD is True, always rebuild the kmer frequency histograms
          from the data, i.e. don't load them in from disk if they exist.

   Input files:

      <genome.fastb>  -  reference genome

   Output files:

      <genome.fastb.freq_table.kS> - genomic kmers with frequency

   Part of <reference genome analysis>.
*/

#include "MainTools.h"

#include "math/Functions.h"
#include "kmers/KmerShape.h"

#include "kmer_freq/WriteKmerFrequencies.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA); 
  CommandArgument_KShapes2(K, Ks);
  CommandArgument_Bool_OrDefault(REBUILD, False);
  CommandArgument_String_OrDefault_Doc(GENOME, "genome",
    "Genome filename name minus .fastb extension"); 

  EndCommandArguments; 

  String data_dir = PRE + "/" + DATA;

  String genomeFile = data_dir + "/" + GENOME + ".fastb";
  String kmersFileBase = data_dir + "/" + GENOME + ".freq_table.k";

  vecbasevector genome( genomeFile );

  for ( unsigned int i = 0; i < Ks.size(); ++i ) {
#define CASE(_KSHAPE) \
        if ( REBUILD || ! IsRegularFile( kmersFileBase + ToString(_KSHAPE::getId()) ) ) \
          WriteKmerFrequencies<_KSHAPE>(genome, kmersFileBase + ToString(_KSHAPE::getId()), true )
    
    DISPATCH_ON_KSHAPE(Ks[i], CASE);
  }

  return 0;
}
