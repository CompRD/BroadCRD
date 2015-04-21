/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "BasevectorTools.h"
#include "graphics/BasicGraphics.h"
#include "math/SparseArray.h"
#include "lookup/LookAlign.h"
#include <math.h>



/// Simple short-hand wrapper that computes counts of alignment start 
/// positions on the reference. \c ref is a reference genome (fw strand only).
/// PlacementInputIterator: an iterator
/// into a container (e.g. vector) of \em any placement- or alignment-like
/// objects that describe placements of some queries onto the reference
/// and support accessor methods ::StartOnTarget(), ::TargetId(),
/// ::IsQueryRC(). Note: it is assumed that the alignments describe
/// positions of either forward or reverse complement queries on the forward
/// strand of the genome only. CountContainer: a two-way mutable random 
/// access container (e.g. int[][], but can be vector<vector<int> > or 
/// vec <SparseArray <int> > etc) - anything that supports indexing operation 
/// \c counter[i][j] and the value indexed in this way can be incremented. 
/// \c counter must be pre-allocated in such a way that any addressing 
/// \c counter[i][j] with i in [0, 2*ref.size()-1] and j in 
/// [0, ref[i/2].size()-1] is legal. 
/// 
/// Upon return, for each i, j, d (d=0 or 1) the \c counter[2*i+d][j] 
/// holds the total number of placements/alignments found in 
/// the collection spanned by [from, to) that start at position j on the
/// forward (d=0) or reverse (d=1) strand of the target (reference) contig i.
/// [Note that the representation is "inverted" in a sense that we do not count
/// both forward and rc query alignments on the forwad reference strand, 
/// but map the alignments of original, unmodified queries onto the appropriate
/// reference strand, either fw or rc, and count the start 
/// positions of these original query sequences on both strands of the reference 
/// independently.


template <class PlacementInputIterator, class CountContainer>
void CountStartsOnRef(const vecbasevector & ref,
		 PlacementInputIterator from, 
		 PlacementInputIterator to, 
		 CountContainer & counter) {

  for ( PlacementInputIterator it = from; it != to ; it++ ) {
      unsigned int pos = it->StartOnTarget();
      unsigned int tig = it->TargetId();
      if ( tig >= ref.size() ) {
	FatalErr("ERROR::Inconsistent target index: " << tig);
      }
      if ( pos > ref[tig].size() ) {
	FatalErr("ERROR::Illegal position "<< pos << " on target " << tig);
      }
      if ( it->IsQueryRC( ) ) pos = ref[ tig ].size( ) - pos - 1;
      ++counter[ 2*tig + ( it->IsQueryRC() ? 1 : 0 ) ][pos];    
  }
}


/// Same as CountStartsOnRef(const vecbasevector &, PlacementInputIterator,
/// PlacementInputIterator, CountContainer &) (see docs), but accepts alignments
/// stored in AlignCollector instead of begin()/end() iterators into an STL-style
/// container. Only uniquely aligned reads are counted.
template <class AlignCollector, class CountContainer>
void CountUniqueStartsOnRef(const vecbasevector & ref,
		            const AlignCollector & aligns,
			    CountContainer & counter) {

  for ( unsigned int i = 0 ; i < aligns.size() ; i++ ) {
      if ( ! aligns.UniquelyAligned(i) ) continue;
      const look_align & la = aligns.Align(i);

      unsigned int pos = la.StartOnTarget();
      int tig = la.TargetId();
      if ( tig >= ref.size() ) {
	FatalErr("ERROR::Inconsistent target index: " << tig);
      }
      if ( pos > ref[tig].size() ) {
	FatalErr("ERROR::Illegal position "<< pos << " on target " << tig);
      }
      if ( la.IsQueryRC( ) ) pos = ref[ tig ].size( ) - pos - 1;
      ++counter[ 2*tig + ( la.IsQueryRC() ? 1 : 0 ) ][pos];    
  }
}



struct Curve {
  vector<double> x;
  vector<double> y;
  String title;
  void push(double x_, double y_) { x.push_back(x_); y.push_back(y_); }
  unsigned int size() const { return x.size() ; }
  unsigned int NPoints() const { return x.size(); } // alias for size
  void SetTitle(const String & s) { title = s; }
};


/// Generates GC curve for the placements from the
/// [from, to) range on the reference \c ref and returns GC
/// bias metrics. 
/// PlacementInputIterator must reference a type that
/// has methods ::StartOnTarget(), ::TargetId(), and ::IsQueryRC()
/// defined
///
/// @param ref Reference genome, on which the placements
/// are specified
/// @param from Iterator (start) into the collection of placements 
/// @param to Iterator (end) into the collection of placements
/// @param W size of the sliding window, within which GC content is to be computed
/// @param c[out] GC curve will be saved into this object
/// @param SaveMem if true, the method will use slower but much more memory
/// efficient internal data representation. For small reference sizes, this option
/// is ignored (set to false automatically regardless of the passed value)
/// @param MaxSlop threshold parameter used to exclude from GC curve the points
/// that seem to fluctuate too wildly

template <class PlacementInputIterator>
double ComputeGCCurve(
	       const vecbasevector & ref,
	       PlacementInputIterator from,
	       PlacementInputIterator to,
	       int W,
 	       Curve & c,
	       bool SaveMem = false,
	       double MaxSlop = 0.05) {

     vec<int> GC_ref( W+1 , 0 ) , GC( W+1, 0 );

     unsigned int genome_size = 0;
     for ( size_t i = 0 ; i < ref.size() ; i++) {
         genome_size+= (ref[i].size() << 1);
     }
     // ignore savemem if genome is really small
     if ( genome_size < 50000000u ) SaveMem = false; 
     
     if ( SaveMem ) {
	      
         // compact counter, slower but saves memory: will physically
         // store the counts only for positions on the reference (both
         // strands), where these counts are non-zero
         vec< SparseArray<int,hash_map<unsigned int, int> > > starts( 2 * ref.size( ) );
	 for ( unsigned int i = 0; i < starts.size( ); i++ ) {
	      // set the total size of the "array" of counts
	      starts[i].resize( ref[i/2].size( ), 0 );

	      // give a hint about the required storage size N
	      // (physical storage, not the total "size" of the sparse array!).
	      // Different specializations of SparseArray know
	      // how to deal with the hint (or ignore it if they don't);
	      // for instance, the implementation backed with hash_map
	      // will try to reserve N buckets for most efficient access.
	      starts[i].ReserveStorage( static_cast<unsigned int>
			 (   ((double)ref[i/2].size()/(double)genome_size)*
			     distance(from,to)
			  )
		      );
	 }
	 // for each base on the fw and rc strands of the reference
	 // genome, count the number of placements/alignments that
	 // start there
	 CountStartsOnRef(ref, from, to, starts);

	 // Generate stats for GC content of windows starting at every
	 // position on the reference and at the read start 
	 // points.
	 ComputeWindowGCCounts(GC_ref, GC, W, ref, starts);
	    
     } else {
         // flat representation (fast but requires huge memory
         // for large genomes): reserve an individual int counter
         // for each and every position on both strands of reference genome
         vec< vec<int> > starts( 2 * ref.size( ) );
	 for ( int i = 0; i < starts.isize( ); i++ ) {
	     starts[i].resize( ref[i/2].size( ), 0 ); // allocate!
	 }

	 // for each base on the fw and rc strands of the reference
	 // genome, count the number of placements/alignments that
	 // start there
	 CountStartsOnRef(ref, from, to, starts);

	 // Generate stats for GC content of windows starting at every
	 // position on the reference and at the read start 
	 // points.
	 ComputeWindowGCCounts(GC_ref, GC, W, ref, starts);
	    
     }

     // Generate GC representation ratios GC[i]/GC_ref[i].
     vec<double> GC_ratio(W+1,0.0);
     for ( int i = 0; i <= W; i++ ) {    
	 if ( GC_ref[i] > 0 ) {
	      GC_ratio[i] = double( GC[i] ) / double( GC_ref[i] );    
	 }
     }
     if ( GC_ref[W/2] == 0 ) {    
         cerr << "Not a single window with GC content 50% found on the reference.\n" 
	      << "Giving up.\n";
	 exit(1);    
     }
     if ( GC[W/2] == 0 ) {    
         cerr << "Not a single read with associated window's GC content 50% found.\n"
	      << "Giving up.\n";

	 for ( int i = 0 ; i <=W ; i++ ) {
	   cerr << ( ((double)i/W)*100.0 ) << "\t" << GC[i] << "\t" << GC_ratio[i] << endl;
	 }
	 exit(1);    
     }
     double mid = GC_ratio[W/2];
     //     if (VERBOSE) PRINT(mid);

     // normalize the ratio of (GC=n for reads)/(GC=n on the reference)
     // frequencies (make ratio = 1 at n=W/2):
     for ( int i = 0; i <= W; i++ ) {    
       if ( GC_ref[i] > 0 ) {    
	   GC_ratio[i] /= mid;
	   //	   if (VERBOSE) PRINT2( i, GC_ratio[i] );    
       }    
     }

     // Our GC counts are a "sample" - treat them as random
     // Poisson events, approximate Poisson with Gauss - both
     // mean (i.e. observed count) and (expected) std. dev.
     // equal to Poisson's lambda parameter. Estimate
     // (gaussian)  95% confidence interval - even though we 
     // actually have just one measurement.
     double p = 0.95;
     double z = InverseNormalCDF( 1.0 - (1.0-p)/2.0 );
     for ( int i = 0; i <= W; i++ ) {    
          double u = GC[i];
	  if ( GC_ref[i] > 0 && GC[i] > 0 ) {    
	      double low = u - z * sqrt(u), high = u + z * sqrt(u);
	      low = Max( low, 0.0 );
	      // If our reads sample the genome uniformly,
	      // we expect GC_U[i]=GC_ref[i]*mid reads with GC count i.
	      // Instead, we counted GC[i] reads with GC count i
	      // and the confidence interval of the measurement 
	      // is estimated as (low,high).
	      // If this estimated confidence interval is larger
	      // than the MAX_SLOP fraction of GC_U[i] (not the
	      // actually measured GC[i]!!), we ignore the measurement.
	      // Otherwise, save graph points:
	      if ( (high-low) / ( GC_ref[i] * mid ) < MaxSlop ) {
                    // if ( low > 0 && high/low <= 1.05 )
		c.push( 100.0 * double(i)/double(W) , // x - gc percentage
			GC_ratio[i] );    // y - normalized ratio
	      }    
	      
	  }    
     }


     if ( c.size() < 2 ) {
          cerr << "ERROR:: No or just one point on the curve. "
                 << "Curve/metric make no sense.  Exiting.\n";
	  exit(1);
     } 

     double biasGCValL2 = 0.0;
     double h, fx1, fx0;
     double length = c.x[c.size()-1] - c.x[0]; // length of the covered x-interval
  
     // compute the area under ( 1-y)^2 curve, where y(z) is our
     // normalized curve of the ratio of number of reads with GC percentage z
     // to the number of windows on the reference with the same GC percentage.
     // Thus what we measure here is the overall deviation of y from 1, i.e.
     // the degree of sample non-uniformity, with which reads represent the 
     // reference.
     for ( unsigned int i = 0; i < c.size( ) - 1; i++ ) {
         h = c.x[i+1] - c.x[i];
	 fx0 = (1.0 - c.y[i])*(1.0 - c.y[i]);
	 fx1 = (1.0 - c.y[i+1])*(1.0 - c.y[i+1]);
	 biasGCValL2 += (h/2.0)*(fx0 + fx1);
     }
     
     // compute the final metric as interval length (sqrt of the area)
     // normalized by the actual length available for computation
     // (since we do not always integrate in the [0%, 100%] range but
     // rather over available points, which usually span smaller,
     // variable interval - the result should not depend on how many
     // points we were lucky to get). Scale the final result by 100
     // - that's arbitrary, of course.
     biasGCValL2 = 100.0*sqrt(biasGCValL2)/length;

     return biasGCValL2;
}


/// Generates GC plot and saves it into a file.
///
/// @param file_name name of the file to save the plot into
/// @param curves vector of precomputed curves to plot (x=GC percentage,
/// y=normalized frequency of reads); each curve's title (if it is set)
/// will be shown in the plot's legend
/// @param x_label_text Text to print for the X axis label (if specified)
/// @param colors if specified, colors[i] color will be used for curves[i]
/// curve; if not specified up to 8 default colors will be used (can not
/// pass more than 8 curves in this case!)
void GenerateGCPlot(String & file_name,
	     const vec<Curve> & curves,
	     const char * x_label_text = 0,
	     const vec<color> * colors = 0,
             const double TITLE_FONTSIZE = 15 );

