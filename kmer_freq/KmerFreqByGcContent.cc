/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: KmerFreqByGcContent
   
   Generates histograms of kmer frequency verus kmer count for all possible
   GC content given either a set of reads or precomputed kmer frequency tables.
   It will only compute the kmer frequency tables for each K if these tables
   do not already exist.
   Histogram filename is given in HIST_OUT.
   Values of K to compute histogram are given in K. To specify multiple
   values of K at once use the form: K="{16,20}"
   Output histogram files for each value of K are formatted:
   Row = freq (0->MAX_FREQ),  Column = GC content (0->K)

   If COMPARE_GENOMIC is True (default is False) then also calculates two additional
   histograms. One from those kmers that are genomic, and one from those kmers that
   aren't. Requires that kmer frequency tables have already been generated for
   genome.fastb in the DATA directory. Output histogram files are formatted as
   described above and are named HIST_OUT_genomic.Kx and HIST_OUT_nongenomic.Kx

   If REBUILD is True, always rebuild the histograms from the data, i.e. don't load 
   them in from disk if they exist.

   \ingroup grp_gc
   \file
*/


#include "MainTools.h"

#include "kmers/KmerShape.h"
#include "kmer_freq/KmerFrequencyTable.h"
#include "kmer_freq/WriteKmerFrequencies.h"
#include "math/Functions.h"

// Write gc histogram file with the following format:
// Row = freq (0->maxFreq), Column = GC content (0->K)
void writeHistogram( vec<vec<longlong> >& gc_hist, String filename, const int maxFreq ) {
  Ofstream( out, filename);
  for (int freq = 0; freq <= maxFreq; freq++) {
    for (int gc = 0; gc < gc_hist.isize(); gc++) 
      out << gc_hist[gc][freq] << " ";
    out << "\n";
  }
  out.close();
}


int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA); 
  CommandArgument_String(RUN); 
  CommandArgument_String(K);
  CommandArgument_String(READS_IN);
  CommandArgument_Int_OrDefault(MAX_FREQ,200);
  CommandArgument_Bool_OrDefault(COMPARE_GENOMIC, False);
  CommandArgument_Double_OrDefault(THRES, 1.25);
  CommandArgument_String_OrDefault(HIST_OUT, "gc_hist");
  CommandArgument_Bool_OrDefault(REBUILD, False);
  EndCommandArguments; 

  vec<KmerShapeId> Ks;
  ParseKmerShapeIdSet( K, Ks, false );
  ForceAssertSupportedKShapes( Ks );

  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;

  String kmersFileBase = run_dir + "/" + READS_IN + ".nonunique.k";
  String histFileBase = run_dir + "/" + HIST_OUT;
  String genomicFileBase = data_dir + "/genome.fastb.freq_table.k";

  vecbasevector allReads;
  vec<KmerFrequencyTable*> tablePtrs;
  
  for ( unsigned int i = 0; i < Ks.size(); ++i ) {
#define CASE(K) \
          if (REBUILD || !IsRegularFile(kmersFileBase + ToString(K::getId()))) { \
            if (allReads.empty()) \
	      allReads.ReadAll(run_dir + "/" + READS_IN); \
	    WriteKmerFrequencies<K>( allReads, kmersFileBase + ToString(K::getId()) ); \
          } \
          tablePtrs.push_back( new KmerFrequencyTable( K::getId(), \
            kmersFileBase+ToString(K::getId()) ) )
    DISPATCH_ON_KSHAPE( Ks[i], CASE );
  }

  // Generate Histogram for each vale of K
  for ( unsigned int i = 0; i < Ks.size(); ++i ) {
    cout << "Generating Histogram for K=" << Ks[i] << "\n";

    int localMin = tablePtrs[i]->GetFirstLocalMin(THRES);
    int peak = tablePtrs[i]->GetPeak();
    cout << "All GC values local min = "  << localMin << ", peak = " << peak << "\n";

    // Write GC summary file
    Ofstream( summary, run_dir + "/" + HIST_OUT + ".k"+ToString(Ks[i]) + ".summary" );
    for (int gc = 0; gc < GetKmerSize(Ks[i]) + 1; gc++) {
      localMin = tablePtrs[i]->GetFirstLocalMin(THRES, gc);
      peak = tablePtrs[i]->GetPeak(gc);
      PRINT3(gc, localMin, peak);
      summary << gc << "\t" << localMin << "\t" << peak << "\n";
    }
    summary.close();

    // Create gc histogram file for this value of K
    vec<vec<longlong> > gc_hist;
    tablePtrs[i]->GetGCHistograms( gc_hist );
    writeHistogram(gc_hist, histFileBase + ".k"+ ToString(Ks[i]), MAX_FREQ);

    // Write genomic/non genomic histograms
    if (COMPARE_GENOMIC) {
      if (!IsRegularFile(genomicFileBase + ToString(Ks[i])))
	cout << "Warning: Could not find file containing genomic kmers: \n"
	     << genomicFileBase + ToString(Ks[i]) << "\n";
      else {	
	// Load KmerShortMap of genomic kmers
	KmerShortMap genomicKmers(Ks[i], genomicFileBase + ToString(Ks[i]));
      
	vec<vec<longlong> > gc_hist_genomic;
	vec<vec<longlong> > gc_hist_nongenomic;

	// Generate and write genomic and nongenomic histograms
	tablePtrs[i]->GetGCHistograms( gc_hist_genomic, gc_hist_nongenomic, genomicKmers );
	writeHistogram(gc_hist_genomic, histFileBase + "_genomic.k"+ ToString(Ks[i]), MAX_FREQ);
	writeHistogram(gc_hist_nongenomic, histFileBase + "_nongenomic.k"+ ToString(Ks[i]), MAX_FREQ);
      }
      
      cout << "\n";
    }
  }
  
  return 0;
}
