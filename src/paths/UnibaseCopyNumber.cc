/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: UnibaseCopyNumber

   Determines copy number of each unibase using read alignments.
   Adjusts for GC bias.  Currently this is only implemented for K=20.  To implement
   other values, add another COUNT_HITS line.

   INPUT FILES:
     reads.fastb
     reads.gc_bias.k*
     reads.trusted_unibases.corrected.k*.fastb
     reads.trusted_unipaths.k*.fastb
     reads.trusted_unipaths.k*
     reads.trusted_unipathsdb.k*

   OUTPUT FILES:
     reads.<UNIPATHS>.predicted_count.k*
   
   CACHED FILES:
     reads.ilt.qltout.query_is_aligned (NO, NAME CHANGED.)
     reads.ilt.qltout.num_aligns_per_target (NO, NAME CHANGED.)

   @file
*/

#include "Basevector.h"
#include "kmers/KmerRecord.h"
#include "MainTools.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"
#include "paths/KmerPathDatabase.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"
#include "paths/UnipathCoverageCore.h"

int main( int argc, char** argv ) {
  
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String(UNIBASES);
  CommandArgument_Int(K);
  CommandArgument_String_OrDefault(READS, "reads");
  CommandArgument_String_OrDefault(READS_EC, "reads_ec");
  CommandArgument_String_OrDefault_Doc( UNIPATHS, "unipaths", "Only affects output file name" ); 
  CommandArgument_String_OrDefault(BIAS_CURVE, "");
  CommandArgument_Double_OrDefault(THRESH, 0.01);
  CommandArgument_Double_OrDefault(ERR_RATE, 0.0);
  CommandArgument_Bool_OrDefault(WRITE, True);
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  EndCommandArguments;

     // Define directories.
  
     String data_dir = PRE + "/" + DATA;
     String run_dir = data_dir + "/" + RUN;

     // Get genome size.  Note that this is cheating.

     longlong genomeSize = StringOfFile( run_dir + "/genome.size", 1 ).Int( );
     ForceAssert( genomeSize > 0 );
     // Load bias curve.
  
     vec<double> biasCurve;
     if ( ! BIAS_CURVE.empty() )
          READX( run_dir + "/" + BIAS_CURVE + ".k" + ToString(K), biasCurve );

     // Load unibases and find involution of it.

     String unibasesFile =  run_dir + "/" + READS + "." + UNIBASES 
       + ".k" + ToString(K);
     if (!IsRegularFile(unibasesFile) && IsRegularFile(unibasesFile + ".fastb"))
       unibasesFile += ".fastb";
     vecbasevector unibases(unibasesFile);
     vec<unipath_id_t> toRc;
     UnibaseInvolution( unibases, toRc );

     // Get the read lengths.

     String ReadsFile = run_dir + "/" + READS_EC + ".fastb";
     size_t nreads = MastervecFileObjectCount(ReadsFile);
     vec<int> readlen(nreads);
     int batchsize = 100000;
     {    vecbasevector reads;                                                
          for ( size_t i = 0; i < nreads; i += batchsize )
          {    int start = i, stop = Min( nreads, i + batchsize );            
               reads.ReadRange( ReadsFile, start, stop );
               for ( size_t j = 0; j < reads.size( ); j++ )
                    readlen[i+j] = reads[j].size( );
               reads.clear( );    }    }

     // Align the reads to the unibases.  We do this by finding the start point of
     // each read in the unibases.

     #define COUNT_HITS(KVAL)                                                      \
          else if ( K == KVAL )                                                    \
          {                                                                        \
               /* Form kmer records associated with first kmers in the reads. */   \
                                                                                   \
               const int K(KVAL);                                                  \
               vec< kmer_record<K,1> > readstarts, unistarts;                      \
               readstarts.reserve(nreads);                                         \
               vecbasevector reads;                                                \
               for ( size_t i = 0; i < nreads; i += batchsize )                    \
               {    size_t start = i, stop = Min( nreads, i + batchsize );         \
                    reads.ReadRange( ReadsFile, start, stop );                     \
                    for ( size_t j = 0; j < reads.size( ); j++ )                   \
                    {    if ( reads[j].size( ) >= static_cast<unsigned>(K) )       \
                         {    static basevector b;                                 \
                              b.SetToSubOf( reads[j], 0, K );                      \
                              static kmer_record<K,1> rec;                         \
                              rec.Set( b, i+j, 0 );                                \
                              readstarts.push_back(rec);    }    }                 \
                    reads.clear( );    }                                           \
                                                                                   \
               /* Form kmer records associated with kmers in the unibases. */      \
                                                                                   \
               size_t nrecs = 0;                                                   \
               for ( size_t i = 0; i < unibases.size( ); i++ )                     \
               {    if ( unibases[i].size( ) >= static_cast<unsigned>(K) )         \
                         nrecs += ( unibases[i].size( ) - K + 1 );    }            \
               unistarts.reserve(nrecs);                                           \
               for ( size_t i = 0; i < unibases.size( ); i++ )                     \
               {    for ( unsigned j = 0; j <= unibases[i].size( ) - K; j++ )      \
                    {    static basevector b;                                      \
                         b.SetToSubOf( unibases[i], j, K );                        \
                         static kmer_record<K,1> rec;                              \
                         rec.Set( b, i, j );                                       \
                         unistarts.push_back(rec);    }    }                       \
                                                                                   \
               /* Match them up. */                                                \
                                                                                   \
               Sort(readstarts), Sort(unistarts);                                  \
               size_t ulast = 0;                                                   \
               for ( size_t i = 0; i < readstarts.size( ); i++ )                   \
               {    size_t j;                                                      \
                    Bool eq = False;                                               \
                    for ( j = ulast; j < unistarts.size( ); j++ )                  \
                    {    if ( unistarts[j] > readstarts[i] ) break;                \
                         if ( unistarts[j].EqualKmers( readstarts[i] ) )           \
                         {    eq = True;                                           \
                              break;    }    }                                     \
                    if ( j == unistarts.size( ) ) break;                           \
                    if (eq)                                                        \
                    {    int id = readstarts[i].GetId( );                          \
                         int uid = unistarts[j].GetId( );                          \
                         if ( static_cast<unsigned>(unistarts[j].GetPos( ) + readlen[id]) \
                              <= unibases[uid].size( ) )                           \
                         {    numAligns[uid]++;                                    \
                              numAligns[ toRc[uid] ]++;                            \
                              readIsAligned[id] = True;    }    }                  \
                    ulast = j;    }    }
     vec<int> numAligns( unibases.size( ), 0 );
     String ALIGNS = READS_EC + ":" + UNIBASES;
     String numAlignsFile = run_dir + "/" + ALIGNS + ".num_aligns_per_target";
     vec<Bool> readIsAligned( nreads, False );
     String readIsAlignedFile = run_dir + "/" + ALIGNS + ".query_is_aligned";
     if ( IsRegularFile(numAlignsFile) && IsOlder( ReadsFile, numAlignsFile )
          && IsOlder( unibasesFile, numAlignsFile ) )
     {    BinaryReader::readFile( numAlignsFile, &numAligns );
          BinaryReader::readFile( readIsAlignedFile, &readIsAligned );    }
     else
     {    if ( 0 == 1 );
          COUNT_HITS(19)
          COUNT_HITS(20)
          COUNT_HITS(64)
          COUNT_HITS(80)
          COUNT_HITS(96)
          COUNT_HITS(100)
          else
          {    cout << "COUNT_HITS not implemented for K = " << K << "." << endl;
               cout << "Abort." << endl;
               exit(1);    }
          BinaryWriter::writeFile( numAlignsFile, numAligns );
          BinaryWriter::writeFile( readIsAlignedFile, readIsAligned );    }

  int numReadsAligned = Sum( readIsAligned );

  longlong sumLens = 0;
  for ( size_t r = 0; r < nreads; ++r )
    if ( readIsAligned[r] ) sumLens += readlen[r];
  double aveAlignedReadLen = (double) sumLens / (double) numReadsAligned;
  PRINT( aveAlignedReadLen );

  PRINT( (double)sumLens / (double)genomeSize );
  PRINT( (double)numReadsAligned / (double)genomeSize );

  VecPdfEntryVec cn_pdfs( unibases.size() );
  for ( size_t u = 0; u < unibases.size(); ++u ) {
    if ( !( u % 200 ) )
      DPRINT( u );
    int numKmers = unibases[u].size() - K + 1;
    ForceAssertGe( numKmers, 0 );
    if ( ! biasCurve.empty() ) {
      int gc = unibases[u].GcBases( 0, K );
      double biasAve = 0;
      for ( int kstart = 0; kstart < numKmers; ++kstart ) {
        biasAve += biasCurve[gc];
        gc -= IsGC( unibases[u][kstart] );
        if ( kstart + K < (int) unibases[u].size() )
          gc += IsGC( unibases[u][kstart+K] );
      }
      biasAve /= (double) numKmers;
      numAligns[u] = (int) round( (double)numAligns[u] / biasAve );
    }
    PdfEntryVec pdf;
    double expectedNumAligns;
    // ImperfectLookup only produced alignments for reads that are subsumed by the
    // unipath, so we need to adjust the number of positions used for the coverage
    // calculation.  CopyNumber() will add aveReadLen to whatever we pass in, and there
    // are (unipath_length - average_read_length) places for a fully contained alignment
    // to start, which is the value we want CopyNumber() to end up with.
    int numAlignStartPositions = unibases[u].size() - (int)(round(2.*aveAlignedReadLen));
    int minPositions = 1 - (int)(floor(aveAlignedReadLen-1));
    if ( numAlignStartPositions < minPositions ) numAlignStartPositions = minPositions;
    CopyNumber( numAligns[u], numAlignStartPositions, aveAlignedReadLen, numReadsAligned,
                genomeSize, K, pdf, THRESH, ERR_RATE, &expectedNumAligns );
    if (VERBOSE) {
      PRINT5( u, unibases[ u ].size(), numAligns[ u ], pdf, expectedNumAligns );
      double maxProb = -1.0;
      int maxProbCN = -1;
      for ( unsigned int i = 0; i < pdf.size(); ++i )
        if ( pdf[i].Prob() > maxProb ) {
          maxProb = pdf[i].Prob();
          maxProbCN = pdf[i].NumCopies();
        }
      PRINT2( maxProbCN, maxProb );
    }
    cn_pdfs[u] = pdf;
  }

  if (WRITE)
    cn_pdfs.WriteAll( (run_dir + "/" + READS + "." + UNIPATHS + ".predicted_count.k" + ToString(K)).c_str() );
}



