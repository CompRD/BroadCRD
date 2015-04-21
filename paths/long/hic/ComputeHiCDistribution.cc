///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "paths/long/hic/HiCDefs.h"
#include "lookup/SAM.h"
#include "ParallelVecUtilities.h"
#include <cmath>

static const String base_dir = "/wga/scr4/HiC/bwa-mem";
static const String suffix = ".paired.sorted.bam";
static const vec<String> all_libs = {
        "HIC089","HIC117","HIC118","HIC152","HIC156","HIC162","HIC163",
        "HIC164","HIC165","HIC167","HIC168",
        "HIC169","HIC170","HIC171","HIC176","HIC177","HIC179_comb",
        "HIC183","HIC184","HIC187","HIC198a",
        "HIC198b","HIC207","HIC228","HIC229","HIC231","HIC236","HIC237"
};


int main( int argc, char* argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_UnsignedInt_OrDefault_Doc(BINS,0,"number of bins");
     CommandArgument_UnsignedInt_OrDefault_Doc(BINSIZE,0,"size of bins");
     CommandArgument_UnsignedInt_OrDefault_Doc(LO,1000,"low value to include");
     CommandArgument_UnsignedInt_OrDefault_Doc(HI,0,"hi value to include");
     CommandArgument_String_OrDefault_Doc(REGION,"","region to extract from indexed bam(s)");
     CommandArgument_String_OrDefault_Doc(OUTHEAD,"","write <OUT_HEAD>.sepdist");
     CommandArgument_StringSet_Doc(LIBS, "set of libraries (e.g. HIC117) or ALL" );
     CommandArgument_Bool_OrDefault_Doc(SAMPLES,False,
             "if True, prints ALL samples as lines starting with 'DIST'")
     CommandArgument_Bool_OrDefault_Doc(HIST,False,
             "if True, prints histogram buckets as lines starting with 'HIST'");
     EndCommandArguments;


     vec<String> const& libs = ( LIBS.size() == 1 && LIBS[0] == "ALL" ) ? all_libs : LIBS;
     vec<unsigned> dists;
     dists.reserve(100000);

     unsigned max=0;
     size_t total=0, below1Mb=0;
     size_t raw_idx=0;
     size_t dot = 100000;
     cout << "counting in [" << LO << "," << ( !HI ? "Inf":ToString(HI) )
             << ")" << ( REGION != ""  ? " in " + REGION : "" ) << endl;
     for ( auto const& lib : libs ) {
         SAM::BAMFile bam(base_dir + "/" + lib + suffix, REGION);
         SAM::Record rec;
         while ( bam.nextRecord(rec) ) {
             ++raw_idx;
             if ( raw_idx % dot == 0 ) cout << "." << std::flush;
             if ( raw_idx % (4*5*dot) == 0 ) cout << endl;
             else if ( raw_idx % (5*dot) == 0 ) cout << " ";
             if ( raw_idx % (10*4*5*dot) == 0) cout << endl;
             if ( rec.isMapped() && rec.isMateMapped() &&       // both are mapped
                     rec.isFirstReadOfPair() &&                 // arbitrarily do this just for the first of the pair
                     rec.getRefName() == rec.getMateRefName() ) {       // only consider the same chromosome
                 // other things to think about -- correct distances for innies vs. outies
                 // (does the length of the read matter) passing filters, duplicates.

                 auto p1 = rec.getRefPos();
                 auto p2 = rec.getMateRefPos();
                 unsigned dist = ( p1 < p2 ) ? ( p2 - p1 ) : ( p1 - p2 );
                 if ( dist >= LO and (!HI || dist <= HI) ) {
                         dists.push_back(dist);
                         max = std::max(dist, max);
                         total++;
                         if ( dist < 1000000 ) below1Mb++;
                         if ( SAMPLES ) cout << "DIST " << dist << endl;
                 }
             }
         }
     }
     PRINT2(below1Mb, total);
     cout << "fraction in that range below 1Mb sep=" << static_cast<double>(below1Mb) / total << endl;

     // ideal would be non-uniform binning, but we're fine with a lot of bins
     // we're going to smooth the distribution a bit
     cout << Date() << ": sorting dists" << endl;
     ParallelSort(dists);
     cout << Date() << ": done sorting dists" << endl;

     size_t n = dists.size();
     size_t iq1 = n*0.25;
     size_t iq3 = n*0.75;
     unsigned q1 = dists[iq1];
     unsigned q3 = dists[iq3];
     size_t iqr = q3 - q1;

     cout << "number of samples is " << n << endl;
     cout << "quartiles 1,3=" << q1 << "," << q3 << ", inter-quartile range=" << iqr << endl;
     if ( BINS == 0 && BINSIZE == 0 ) {
         cout << "computing bins based on F-D rule" << endl;
         size_t fd_bins = 2*iqr / pow(n, (1./3.));          // freedman-diaconis rule
         BINS = fd_bins;
     } else if ( BINS == 0 ) {
         cout << "computing bins based on BINSIZE" << endl;
         BINS = max / BINSIZE;
     }

     cout << "number of bins=" << BINS << endl;

     if ( HI == 0 ) HI = max;

     SepHistogram hist( LO, HI+1., BINS );
     for ( auto const dist : dists ) hist.Add(dist);
     hist.Print(std::cout, "HIST ");

     if ( OUTHEAD != "" ) {
         cout << Date() << ": writing " << OUTHEAD+".sepdist" << endl;
         BinaryWriter::writeFile(OUTHEAD+".sepdist", hist);
     }

     return 0;
}

