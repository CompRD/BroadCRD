///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Program: UnibaseCopyNumber2

   Determines copy number of each unibase using k-mer counts (experimental).
   Adjusts for GC bias.

   INPUT FILES (in run dir):
     <READS_EC>.fastb
     <READS>.<UNIBASES>.k<K>

   OUTPUT FILES:
     reads.<UNIPATHS>.predicted_count.k<K>
   
   CACHED FILES:
     <READS_EC>.<K>merParcels
     <READS>.<UNIBASES>.<K>merParcels
   
   @file
*/

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/PdfEntry.h"
#include "paths/UnibaseUtils.h"
#include "paths/Ulink.h"
#include "paths/UnibaseCopyNumber2Core.h"
#include "paths/UnipathFixerTools.h"

Bool cmp_ridx( const segalign& a1, const segalign& a2 ){    
  return a1.rid < a2.rid;    
}

void MakeUlinksX( const String& run_dir, const String& JUMP_READS, const int K, 
     const vec<int>& to_rc, vec<segalign>& JSEGS, const vecbasevector& unibases, 
     const vec<int>& innie_sep, const vec<int>& innie_dev, 
     const vec<double>& innie_percent, const double min_innie_percent,
     const int min_kmers, const int max_devs, vec<ulink_with_uids>& ulinks ){
  
  cout << Date( ) << ": set up for link computation" << endl;
  String head = run_dir + "/" + JUMP_READS;
  vecbasevector reads( head + ".fastb" );
  uint64_t nreads = reads.size( );
  vec<int> readlen( reads.size( ) );
  for ( size_t i = 0; i < reads.size( ); i++ )
    readlen[i] = reads[i].size( );
  Destroy(reads);
  PairsManager pairs( head + ".pairs" );
  vec<segalign>& SUGS = JSEGS;
  ParallelSort( SUGS, cmp_ridx );
  vec<size_t> S_START_RID(nreads+1);{   
    size_t SEG_POS = 0;
    for ( uint64_t rid = 0; rid <= nreads; rid++ ){    
      while( SEG_POS < SUGS.size( ) && SUGS[SEG_POS].rid < rid ) ++SEG_POS;
      S_START_RID[rid] = SEG_POS;    
    }    
  }
  cout << Date( ) << ": compute links from "
       << ToStringAddCommas( pairs.nPairs( ) ) << " pairs" << endl;
#pragma omp parallel for
  for ( size_t i = 0; i < pairs.nPairs( ); i++ ){    
    longlong id1 = pairs.ID1(i), id2 = pairs.ID2(i);
    
    // Require that jump reads have only one placement.  This may be
    // too stringent.
    
    if ( S_START_RID[id1+1] - S_START_RID[id1] > 1 ) continue;
    if ( S_START_RID[id2+1] - S_START_RID[id2] > 1 ) continue;
    
    for (size_t j1 = S_START_RID[id1]; j1 < S_START_RID[id1+1]; j1++){    
      for (size_t j2 = S_START_RID[id2]; j2 < S_START_RID[id2+1]; j2++){   
	int u1 = SUGS[j1].u, u2 = to_rc[ SUGS[j2].u ];
	if ( u1 == u2 || u1 == to_rc[u2] ) continue;
	if ( unibases[u1].isize( ) - K + 1 < min_kmers ) continue;
	if ( unibases[u2].isize( ) - K + 1 < min_kmers ) continue;
	int nu1 = unibases[u1].size( ), nu2 = unibases[u2].size( );
	int nr1 = readlen[id1], nr2 = readlen[id2];
	int upos1 = SUGS[j1].upos, upos2 = SUGS[j2].upos;
	int lib = pairs.libraryID(i);
	// Library -> separation/deviation
	int sep_in = innie_sep[lib] + nr1 + nr2
	  - ( nu1 - upos1 + SUGS[j1].rpos + nu2 
	      - upos2 + SUGS[j2].rpos );
	int dev_in = innie_dev[lib];
	int sep_out = pairs.getLibrarySep(pairs.libraryID(i))
	  - ( upos1 - SUGS[j1].rpos
	      + upos2 - SUGS[j2].rpos ) + nr1 + nr2;
	int dev_out = pairs.getLibrarySD(pairs.libraryID(i));
	
	// Add this link to the ulinks list - in both fw
	// and rc form, as both an innie and an outie.
	
#pragma omp critical
	{    
	  if ( innie_percent[lib] >= min_innie_percent ){    
	    ulinks.push( u1, u2, sep_in, dev_in, 
			 Min( nu1, upos1 + nr1 ), 
			 Max( 0, nu2 - upos2 - nr2 ), 1 );
	    ulinks.push( to_rc[u2], to_rc[u1], sep_in, 
			 dev_in, Min( nu2, upos2 + nr2 ), 
			 Max( 0, nu1 - upos1 - nr1 ), 1 );    
	  }
	  ulinks.push( u2, u1, sep_out, dev_out,
		       Min( nu2, nu2 - upos2 ), Max(0, upos1), 1 );
	  ulinks.push( to_rc[u1], to_rc[u2], sep_out, 
		       dev_out, Min( nu1, nu1 - upos1 ),
		       Max(0, upos2), 1 );    }    }    }    }
  cout << Date( ) << ": sorting " << ToStringAddCommas( ulinks.size( ) ) 
       << " unipath links" << endl;
  ParallelSort(ulinks);    
}

void GapStatsAlt( vec<int> gap, vec<int> gapdev, int& gap_ave, int& gapdev_ave ){
     // If there are less than six gaps, we directly compute their mean.
     // Otherwise, we attempt to remove outliers, as follows.  We sort 
     // the gaps and extract the middle half.  From this middle half, we compute the 
     // mean and standard deviation.  Then we select those gaps lying withing 5 
     // standard deviations of the mean, and form their mean.  

     vec<NormalDistribution> S;
     if ( gap.size( ) >= 6 ){    
       vec<int> mid_gaps;
       SortSync( gap, gapdev );
       for ( unsigned int i = gap.size( )/4; i < 3*(1+gap.size( ))/4; i++ )
	 mid_gaps.push_back( gap[i] );
       float sum1 = 0, sum2 = 0;
       for ( unsigned int i = 0; i < mid_gaps.size( ); i++ ){    
	 sum1 += mid_gaps[i];
	 sum2 += float(mid_gaps[i]) * float(mid_gaps[i]);    
       }
       float n = mid_gaps.size( );
       float mean = sum1/n;
       float sd = sqrt(sum2/n - mean * mean);
       float start = mean - 5 * sd, stop = mean + 5 * sd;
       for ( unsigned int i = 0; i < gap.size( ); i++ ){    
	 if ( start <= gap[i] && gap[i] <= stop )
	   S.push( gap[i], gapdev[i] );    
       }    
     }
     else{    
       for ( int l = 0; l < gap.isize( ); l++ )
	 S.push( gap[l], gapdev[l] );    
     }
     NormalDistribution s = CombineNormalDistributions(S);
     gap_ave = int(round(s.mu_));
     gapdev_ave = int(round(s.sigma_));    
}

int main( int argc, char** argv ) {
  
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int(K);
  CommandArgument_Int_OrDefault( PLOIDY, -1 );
  
  // Infixes for input/output file names.
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_String_OrDefault( READS_EC, "reads_ec" );
  CommandArgument_String_OrDefault( UNIPATHS, "unipaths" );
  CommandArgument_String_OrDefault( UNIBASES, "unibases" );
  
  // Heuristics for bias computation.
  CommandArgument_Bool_OrDefault(CORRECT_BIAS, CORRECT_BIAS_DEFAULT);
  CommandArgument_Double_OrDefault(THRESH, THRESH_DEFAULT);
  CommandArgument_Double_OrDefault(ERR_RATE, ERR_RATE_DEFAULT);
  CommandArgument_Bool_OrDefault(LOWER, False);
  CommandArgument_Bool_OrDefault(LOWER_VERBOSE, False);
  CommandArgument_Bool_OrDefault(CN2_VERBOSE, False);
  
  // Runtime control.
  CommandArgument_Int_OrDefault( NUM_THREADS, 16 );
  CommandArgument_Bool_OrDefault(WRITE, True);
  CommandArgument_Bool_OrDefault(VERBOSE, False);

  CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_filt_cpd");
  CommandArgument_String_OrDefault(FRAG_READS, "frag_reads_filt_cpd");
  CommandArgument_String_OrDefault(FRAG_READS_EDIT, "frag_reads_edit");
  CommandArgument_Int_OrDefault(MAX_PLACEMENTS, 50);
  
  EndCommandArguments;


  // Define directories.
  
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;

  
  // Make sure ploidy file (or values) are available
  if ( PLOIDY <= 0 ){
    if ( IsRegularFile( data_dir + "/ploidy" ) )
      PLOIDY = StringOfFile( data_dir + "/ploidy", 1 ).Int( );
    else 
      FatalErr("Require genome ploidy. Please supply either a ploidy "
	       "file in the data, or a value for the PLOIDY option");
    ForceAssertGt( PLOIDY, 0 );
  }

  // Load unibases.
  String reads_head = run_dir + "/" + READS_EC;
  String unibases_head = run_dir + "/" + READS + "." + UNIBASES;
  
  cout << Date() << ": Loading unibases" << endl;
  String unifile = unibases_head + ".k" + ToString(K);
  vecbasevector unibases(unifile);
  
  cout << Date() << ": Loading reads" << endl;
  vecbasevector reads( reads_head + ".fastb" );

  VecPdfEntryVec cn_pdfs;
  vec<int64_t> n_kmer_hits;
  UnibaseCopyNumber2Core( K, reads, unibases, NUM_THREADS, reads_head, unibases_head, 
	   PLOIDY, VERBOSE, cout, cn_pdfs, n_kmer_hits, THRESH, ERR_RATE, CORRECT_BIAS );

     // Lower copy numbers.

  if ( PLOIDY > 1 ) LOWER = False;
  if ( !LOWER ){    
    if (WRITE){    
      cn_pdfs.WriteAll( (run_dir + "/" + READS + "." + UNIPATHS 
			 + ".predicted_count.k" + ToString(K)).c_str() );
      BinaryWriter::writeFile( run_dir + "/" + READS + "." + UNIPATHS + ".kmer_hits.k"
		    + ToString(K), n_kmer_hits );    
    }
    cout << Date( ) << ": UnibaseCopyNumber2 Done!" << endl;
    return 0;    
  }
  
  // Load innie stats.
  
  vec<int> innie_sep, innie_dev;
  vec<double> innie_percent;
  fast_ifstream iin( run_dir + "/" + JUMP_READS + ".outies" );
  String line;
  while(1){    
    getline( iin, line );
    if ( iin.fail( ) ) break;
    istrstream iline( line.c_str( ) );
    int sep, dev;
    double iper;
    iline >> sep >> dev >> iper;
    innie_sep.push_back(sep), innie_dev.push_back(dev);
    innie_percent.push_back(iper);
    cout << Date( ) << ": jump innies " << sep << " +/- " << dev << " (" 
	 << setiosflags(ios::fixed) << setprecision(1)
	 << iper << "%)" << endl;    
  }
  
  // Align reads.
  
  vec< triple<int64_t,int64_t,int> > ALIGNS, JALIGNS;
  AlignReadsToUnipaths( run_dir, JUMP_READS, FRAG_READS, FRAG_READS_EDIT, 
			True, MAX_PLACEMENTS, unifile, ALIGNS, JALIGNS );
  
  // Set up ancillary data structures for unibases.
  
  size_t nuni = unibases.size( );
  vec<int> to_rc;
  UnibaseInvolution( unibases, to_rc );
  
  // Directly convert alignments into segments, then get rid of alignments.
  
  vec<segalign> SEGS, JSEGS;
  SEGS.resize( ALIGNS.size( ) ), JSEGS.resize( JALIGNS.size( ) );
  for ( size_t i = 0; i < ALIGNS.size( ); i++ ){    
    SEGS[i] = segalign( True, ALIGNS[i].first, 0, ALIGNS[i].second,
			ALIGNS[i].third );    
  }
  Destroy(ALIGNS);
  ParallelSort(SEGS);
  for ( size_t i = 0; i < JALIGNS.size( ); i++ ){    
    JSEGS[i] = segalign( True, JALIGNS[i].first, 0, JALIGNS[i].second,
			 JALIGNS[i].third );    
  }
  Destroy(JALIGNS);
  ParallelSort(JSEGS);
  
  // Compute unibase copy numbers.  
  
  vec<double> CN_raw(nuni), GC(nuni);
  vec<int> hits( nuni, 0 ), ulen(nuni), ids( nuni, vec<int>::IDENTITY );
  for ( size_t i = 0; i < SEGS.size( ); i++ ){    
    hits[ SEGS[i].u ]++;
    hits[ to_rc[ SEGS[i].u ] ]++;    
  }
  for ( size_t i = 0; i < nuni; i++ )
    ulen[i] = unibases[i].size( );
  ParallelReverseSortSync( ulen, ids );
  int sample = Min( 100, (int) nuni );
  int64_t total_bases = 0, total_hits = 0;
  for ( int i = 0; i < sample; i++ ){    
    total_bases += ulen[i]; total_hits += hits[ ids[i] ];    
  }
  for ( size_t i = 0; i < nuni; i++ )
    ulen[i] = unibases[i].size( );
  double hits_per_base = double(total_hits) / double(total_bases);
  for ( size_t u = 0; u < nuni; u++ ){    
    if ( (int) u > to_rc[u] ) continue;
    double cov = ( double(hits[u]) / double(ulen[u]) ) / hits_per_base;
    CN_raw[u] = double(PLOIDY) * cov;
    CN_raw[ to_rc[u] ] = double(PLOIDY) * cov;
    int gc = 0;
    for ( int j = 0; j < unibases[u].isize( ); j++ )
      if ( unibases[u][j] == 1 || unibases[u][j] == 2 ) gc++;
    GC[u] = GC[ to_rc[u] ] = double(gc) / double( unibases[u].isize( ) );
    if (CN2_VERBOSE) {    
      cout << "initially setting copy number of unipath " << u 
	   << "[l=" << unibases[u].isize( ) - K + 1 << ",hits=" << hits[u] 
	   << ",gc=" << setiosflags(ios::fixed) << setprecision(1) 
	   << 100.0 * GC[u] << "%] to " << CN_raw[u] << "\n";    
    }    
  }
  
  // Heuristics.
  
  const double max_devs = 6.0;
  const int min_kmers = 120;
  const int min_links_initial = 2;
  const double min_innie_percent = 5.0;
  
  // For each pair of unipaths that are connected by two or more links, predict
  // their order and separation.  Both orders may be possible.  Note that a more 
  // wholistic approach may be needed.
  
  vec<ulink_with_uids> ulinks;
  MakeUlinksX( run_dir, JUMP_READS, K, to_rc, JSEGS, unibases, innie_sep, 
	       innie_dev, innie_percent, min_innie_percent, min_kmers, max_devs, ulinks );
  
  // Make the list of putative unipath links into a set of condensed links.
  // Here we use GapStats to merge links, and we apply a min_links_initial 
  // threshold.  A condensed links is accepted if its predicted overlap is no 
  // more than K-1, plus slop.
  
  cout << Date( ) << ": condense links" << endl;
  vec<ulink_with_uids> condensed_links;
  vec<int> most_links_from( nuni, 0 ), most_links_to( nuni, 0 );
  vec< vec< pair<ho_interval,int> > > cov_from(nuni), cov_to(nuni);
  vec<double> cn_from(nuni, 0), cn_to(nuni, 0);
  for ( size_t i = 0; i < ulinks.size( ); i++ ){    
    size_t j;
    int u1 = ulinks[i].u1, u2 = ulinks[i].u2;
    for ( j = i + 1; j < ulinks.size( ); j++ )
      if ( ulinks[j].u1 != u1 || ulinks[j].u2 != u2 ) break;
    int nlinks = j - i;
    if ( nlinks >= min_links_initial ){    
      vec<int> seps, devs;
      int sep, dev, start1 = unibases[u1].size( ), stop2 = 0;
      for ( size_t m = i; m < j; m++ ){    
	const ulink_with_uids& f = ulinks[m];
	seps.push_back(f.sep), devs.push_back(f.dev);    
	start1 = Min( start1, f.start1 );
	stop2 = Max( stop2, f.stop2 );    }
      GapStatsAlt( seps, devs, sep, dev );
      if ( sep + max_devs * dev >= -(K-1) ){    
	condensed_links.push( u1, u2, sep, dev, start1, stop2, j-i );
	int dev_mult = 3;
	int Start2 = sep + dev_mult * dev;
	int Stop2 = sep + unibases[u2].isize( ) - dev_mult * dev;
	if ( Start2 < Stop2 ) {    
	  ho_interval h2( Start2, Stop2 );
	  cov_from[u1].push( h2, nlinks );    
	}
	int Start1 = -sep - unibases[u1].isize( ) + dev_mult * dev;
	int Stop1 = -sep - dev_mult * dev;
	if ( Start1 < Stop1 ) {    
	  ho_interval h1( Start1, Stop1 );
	  cov_to[u2].push( h1, nlinks );    
	}
	if ( nlinks > most_links_from[u1] ){    
	  most_links_from[u1] = j-i;
	  cn_from[u1] = CN_raw[u2];    
	}
	if ( nlinks > most_links_to[u2] ){    
	  most_links_to[u2] = j-i;
	  cn_to[u2] = CN_raw[u1];    
	}    
      }    
    }
    i = j - 1;    
  }
  /*
    for( size_t i = 0; i < condensed_links.size( ); i++ )
    {    const ulink_with_uids& l = condensed_links[i];
    cout << "\n" << l.u1 << "[CN_raw=" << CN_raw[l.u1] << ",len=" 
    << unibases[l.u1].size( ) << "] --> " << l.u2 << "[CN_raw=" 
    << CN_raw[l.u2] << ",len=" << unibases[l.u2].size( ) << "]: ";
    cout << l.sep << " +/- " << l.dev << ", start1 = " << l.start1 
    << ", stop2 = " << l.stop2 << ", nlinks = " << l.nlinks
    << endl;    }
  */
  
  const double cn_thresh = 0.7;
  const int min_links_mult = 3;
  for ( size_t u = 0; u < nuni; u++ ){    
    int cnx;
    GetMostLikelyValue( cnx, cn_pdfs[u] );
    if ( cnx == 1 ) continue;
    if ( cn_from[u] >= cn_thresh * CN_raw[u]
	 && cn_to[u] >= cn_thresh * CN_raw[u] ){    
      int max_from = 0, max_to = 0;
      for ( int j = 0; j < cov_from[u].isize( ); j++ )
	max_from = Max( max_from, cov_from[u][j].second );
      for ( int j = 0; j < cov_to[u].isize( ); j++ )
	max_to = Max( max_to, cov_to[u][j].second );
      Bool overlap = False;
      for ( int j1 = 0; j1 < cov_from[u].isize( ); j1++ ){    
	if ( cov_from[u][j1].second < max_from/min_links_mult ) continue;
	for ( int j2 = j1 + 1; j2 < cov_from[u].isize( ); j2++ ){    
	  if ( cov_from[u][j2].second < max_from/min_links_mult ) 
	    continue;
	  if ( Meets( cov_from[u][j1].first, cov_from[u][j2].first ) )
	    overlap = True;    
	}    
      }
      for ( int j1 = 0; j1 < cov_to[u].isize( ); j1++ ){    
	if ( cov_to[u][j1].second < max_to/min_links_mult ) continue;
	for ( int j2 = j1 + 1; j2 < cov_to[u].isize( ); j2++ ){    
	  if ( cov_to[u][j2].second < max_to/min_links_mult ) 
	    continue;
	  if ( Meets( cov_to[u][j1].first, cov_to[u][j2].first ) )
	    overlap = True;    
	}    
      }
      if (overlap) continue;
      if ( LOWER_VERBOSE && (int) u <= to_rc[u] ) {   
	cout << "\nsetting copy number of unipath " << u 
	     << "[l=" << unibases[u].isize( ) - K + 1 << ",hits=" 
	     << hits[u] << ",gc=" << setiosflags(ios::fixed) 
	     << setprecision(1) << 100.0 * GC[u] << "%] to 1\n";
	PRINT3( CN_raw[u], cn_from[u], cn_to[u] );    
	cout << "from:\n";
	for ( int j = 0; j < cov_from[u].isize( ); j++ ){   
	  cout << cov_from[u][j].first << "[" 
	       << cov_from[u][j].second << "]\n";    
	}
	cout << "to:\n";
	for ( int j = 0; j < cov_to[u].isize( ); j++ ){    
	  cout << cov_to[u][j].first << "[" 
	       << cov_to[u][j].second << "]\n";    
	}    
      }
      PdfEntryVec copyno_prob;
      copyno_prob.push_back( pdf_entry( 1, 1 ) );
      cn_pdfs[u] = copyno_prob;    
    }    
  }
  
  if (WRITE)
    {
      cn_pdfs.WriteAll( (run_dir + "/" + READS + "." + UNIPATHS + ".predicted_count.k" + ToString(K)).c_str() );
      BinaryWriter::writeFile( run_dir + "/" + READS + "." + UNIPATHS + ".kmer_hits.k"
		    + ToString(K), n_kmer_hits );
    }
  
  cout << Date( ) << ": UnibaseCopyNumber2 Done!" << endl;
  return 0;
}
