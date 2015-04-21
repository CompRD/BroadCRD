///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "paths/AssemblyCleanupTools.h"
#include "lookup/LookAlign.h"
#include "math/NStatsTools.h"
#include "paths/reporting/ReftigUtils.h"
#include "FastIfstream.h"

/**
 * EvalScaffoldsAlt
 *
 */

Bool cmp_best_ali( const look_align& la1, const look_align& la2 ){    
  if ( la1.GaplessLength() > la2.GaplessLength() ) return True;
  if ( la1.GaplessLength() < la2.GaplessLength() ) return False;
  if ( la1.Errors() < la2.Errors() ) return True;
  if ( la1.Errors() > la2.Errors() ) return False;
  if ( la1.indels < la2.indels ) return True;
  if ( la1.indels > la2.indels ) return False;
  return False;    
}

void get_gapped_lengths( const vec<superb>& scaffolds, vec<int>& glengths ){
  glengths.resize( scaffolds.size() );
  for ( size_t si = 0; si < scaffolds.size(); si++ )
    glengths[si] = scaffolds[si].FullLength();
  return;
}

void get_ungapped_lengths( const vec<superb>& scaffolds, vec<int>& ulengths ){
  ulengths.resize( scaffolds.size() );
  for ( size_t si = 0; si < scaffolds.size(); si++ )
    ulengths[si] = scaffolds[si].ReducedLength();
  return;
}

void get_nstats( vec<int>& NSels, vec<int>& lengths, vec<int>& stats ){
  vec<int> Nidx;
  stats.resize( NSels.size() );
  BasicNStats( lengths, Nidx, &NSels ); 
  for ( size_t i = 0; i < NSels.size(); i++ )
    stats[i] = lengths.at( Nidx.at(i) );
}

int get_CN( vec<look_align> las, const vecbvec& ref ){
  int true_CN = -1;
  vecbvec tseqs;
  if ( las.size() > 1 ){
    sort( las.begin(), las.end(), cmp_best_ali );
    int minErr = las.begin()->Errors(), maxAliLen = las.begin()->GaplessLength();
    for ( size_t i = 0; i < las.size(); i++ ){
      if ( las[i].Errors() <= minErr && las[i].GaplessLength() >= maxAliLen ){
	tseqs.resize( i+1 );
	tseqs.back().SetToSubOf( ref.at(las[i].TargetId()), las[i].StartOnTarget(), 
				 las[i].EndOnTarget() - las[i].StartOnTarget()       );
      }
    }
    int cn = 1;
    bvec seq0 = tseqs[0];
    for ( size_t i = 1; i < tseqs.size(); i++ ){
      if ( tseqs[i] == seq0 ) {
	cn++;
      }
    }
    true_CN = cn;
  }else{
    true_CN = las.isize();
  }
  ForceAssertGe( true_CN, 0 );
   return true_CN;
}

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( ASSEMBLY_HEAD );
  CommandArgument_String_OrDefault( REF, "genome" );
  CommandArgument_Int_OrDefault( MAX_OVERLAP, 95 );
  CommandArgument_Int_OrDefault( MAX_SEP, 20000 );
  CommandArgument_Bool_OrDefault( FORCE, True );
  CommandArgument_String_OrDefault( OUT_DIR, "" );
  EndCommandArguments;


  vec<int> NSels{5,10,25,50,75,90,99};

  int max_overlap = MAX_OVERLAP;
  int max_sep     = MAX_SEP;
  // Dir and file names.
  String scaffolds_file  = ASSEMBLY_HEAD + ".superb";
  String contigs_file    = ASSEMBLY_HEAD + ".contigs.fastb";

  String qlt_file        = ASSEMBLY_HEAD + ".contigs.qlt";
  String ref_file        = REF + ".fastb";
  String ref_lookup_file = REF + ".lookup";

  String misjoin_file    = ASSEMBLY_HEAD + ".misjoins";
  
  String out_dir = (OUT_DIR == "" ? ASSEMBLY_HEAD + ".EvalAlt" : OUT_DIR);
  Mkpath( out_dir );

  vecbvec ref( ref_file );

  vec<superb> scaffolds;
  ReadSuperbs( scaffolds_file, scaffolds );
  cout << "Nscaffolds= " << scaffolds.size() << endl;

  vecbvec contigs( contigs_file );
  cout << "Ncontigs= " << contigs.size() << endl;

  vec<look_align> aligns;
  if ( FORCE || ! IsRegularFile( qlt_file ) ) {
    cout << Date( ) << ": aligning contigs" << endl;
    // Note that this function uses OMP (badly)
    GetAlignsFast( 96, contigs_file, ref_lookup_file, qlt_file, aligns, !FORCE, out_dir ); 
  }else{
    cout << Date() << ": loading contig alignments" << endl;
    LoadLookAligns( qlt_file, aligns );
  }
  cout << Date() << ": found " << aligns.size() << " alignments" << endl;

  vec< vec< size_t > > aligns_index( contigs.size() );
  for ( size_t i = 0; i < aligns.size(); i++ ){
    size_t qid = aligns[i].QueryId();
    aligns_index.at(qid).push_back(i);
  }
  
 
  // find high copy number unipaths
  
  vec<int> trueCN( contigs.size(), 0 );
  for ( unsigned si = 0; si < scaffolds.size(); si++ ){
    for ( int ti = 0; ti < scaffolds[si].Ntigs(); ti++ ){
      int tid = scaffolds[si].Tig(ti);
      if ( aligns_index.at(tid).size() > 1 ){
	vec<look_align> las;
	for ( unsigned i = 0; i < aligns_index[tid].size(); i++ )
	  las.push_back( aligns[ aligns_index[tid][i] ] );
	trueCN[tid] = get_CN( las, ref );
	
      }else{
	trueCN[tid] = aligns_index[tid].isize();
      }
    }
  }

  size_t unaligned_sum = 0;
  for ( size_t tid = 0; tid < trueCN.size(); tid++ )
    if ( trueCN[tid] == 0 )
      unaligned_sum += contigs[tid].size();
  PRINT( unaligned_sum );


  size_t nCNg1 = 0;
  map<size_t,size_t> cnHisto;
  for ( size_t tid = 0; tid < trueCN.size(); tid++ ){
    cnHisto[ trueCN[tid] ]++;
    if ( trueCN[tid] > 1 ) nCNg1++;
  }
  cout << " --------------- Copy Number evaluation\n";
  cout.setf( ios::right);
  cout << setw(4) << "cn" << setw(10) << "ntigs" << "\n";
  for ( map<size_t,size_t>::const_iterator im = cnHisto.begin(); im != cnHisto.end(); im++ )
    cout << setw(4) << im->first << setw(10) << im->second << "\n";
  cout << "\n";
  cout << " total number CN>1 = " << nCNg1 << "\n";
  cout << " total number of sequences = " << contigs.size() << "\n";
  cout << "----------------\n\n";



  // evaluate scaffold separations (gaps). Count and list misjoins and inversions.
  size_t nMisjoins = 0, nCloseBadMisjoins = 0, nCloseAcceptMisjoins = 0, nDistantMisjoins = 0; 
  size_t nTested = 0, nInversions = 0;
  ofstream ofm( misjoin_file.c_str() );
  vec<double> zs;
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    for ( int ti = 1; ti < scaffolds[si].Ntigs(); ti++ ){
      int tid1 = scaffolds[si].Tig(ti -1);
      int tid2 = scaffolds[si].Tig(ti);
      int gap  = scaffolds[si].Gap(ti -1);
      ForceAssert( scaffolds[si].Dev(ti -1) >= 0 );
      int dev  = scaffolds[si].Dev(ti -1) > 0 ?  scaffolds[si].Dev(ti -1) : 1;
      
      if ( aligns_index.at(tid1).size() == 1 && aligns_index.at(tid2).size() == 1 ){
	look_align & la1 = aligns.at( aligns_index[tid1][0] );
	look_align & la2 = aligns.at( aligns_index[tid2][0] );
	if ( la1.TargetId() != la2.TargetId() ) continue;
	
	nTested++;
	
  	double sep = 0.31415, extension = 0.31415; // initialization
	if ( la1.Fw1() && la2.Fw1() ){
	  sep       = la2.StartOnTarget() - la1.EndOnTarget();
	  extension = la2.EndOnTarget() - la1.EndOnTarget();
	}
	else if ( la1.Rc1() && la2.Rc1() ){
	  sep = la1.StartOnTarget() - la2.EndOnTarget();
	  extension = la1.StartOnTarget() - la2.StartOnTarget();
	}
	
	if ( sep == 0.31415 ){
	  ofm << "-----------INVERSION\n";
	  nInversions++;
	  la1.PrintReadableBrief(ofm); 
	  la1.PrintVisual(ofm, contigs.at( la1.QueryId() ), ref.at( la1.TargetId() ) );
	  la2.PrintReadableBrief(ofm); 
	  la2.PrintVisual(ofm, contigs.at( la2.QueryId() ), ref.at( la2.TargetId() ) );
	  ofm << "---------------------------------------------------\n\n";
	}else if ( sep < -max_overlap && extension <= 0 ){
	  ofm << "-----------OVERLAPPING BAD MISJOIN\n";
	  nMisjoins++;
	  nCloseBadMisjoins++;
	  PRINT_TO(ofm,sep);
	  la1.PrintReadableBrief(ofm); 
	  la1.PrintVisual(ofm, contigs.at( la1.QueryId() ), ref.at( la1.TargetId() ) );
	  la2.PrintReadableBrief(ofm); 
	  la2.PrintVisual(ofm, contigs.at( la2.QueryId() ), ref.at( la2.TargetId() ) );
	  ofm << "---------------------------------------------------\n\n";
	}else if ( sep < -max_overlap && extension > 0 ){
	  ofm << "-----------OVERLAPPING ACCEPTABLE MISJOIN\n";
	  nMisjoins++;
	  nCloseAcceptMisjoins++;
	  PRINT_TO(ofm,sep);
	  la1.PrintReadableBrief(ofm); 
	  la1.PrintVisual(ofm, contigs.at( la1.QueryId() ), ref.at( la1.TargetId() ) );
	  la2.PrintReadableBrief(ofm); 
	  la2.PrintVisual(ofm, contigs.at( la2.QueryId() ), ref.at( la2.TargetId() ) );
	  ofm << "---------------------------------------------------\n\n";
	}else if ( sep > max_sep ){
	  ofm << "-----------DISTANT MISJOIN\n";
	  nMisjoins++;
	  nDistantMisjoins++;
	  PRINT_TO(ofm,sep);
   	  la1.PrintReadableBrief(ofm); 
	  la1.PrintVisual(ofm, contigs.at( la1.QueryId() ), ref.at( la1.TargetId() ) );
	  la2.PrintReadableBrief(ofm); 
	  la2.PrintVisual(ofm, contigs.at( la2.QueryId() ), ref.at( la2.TargetId() ) );
	  ofm << "---------------------------------------------------\n\n";
	}

	double z = 10000000000.0;
	if ( sep != 0.31415 ) z = (double) abs( gap -sep ) / (double)dev;
	zs.push_back( z );
      }
    }
  }
  ofm.close();

  cout << "\n----------- Misjoin evaluation\n";
  PRINT6( nTested, nInversions, nMisjoins, nCloseBadMisjoins, nCloseAcceptMisjoins, nDistantMisjoins );
  cout << " ------------\n\n";
  

  vec<int> zts{1,2,3,5,10};
  map<int,size_t> gapHisto;
  for ( size_t i = 0; i < zts.size(); i++ ){
    int zthresh = zts[i];
    for ( size_t j = 0; j < zs.size(); j++ ){
      if ( zs[j] > zthresh )
	gapHisto[zthresh]++;
    }
  }

  cout << "------------ Gap prediction evaluation\n";
  cout << setw(6) << "gap_z" << setw(10) << "n_gaps" << setw(14) << "perc_gaps_out" << endl;
  for ( size_t i = 0; i < zts.size(); i++ ){
    int zthresh = zts[i];
    double perc_gaps = (double)gapHisto[zthresh] / (double)zs.size() * 100.0;
    cout << setw(6) << zthresh << setw(10) << gapHisto[zthresh] << setw(14) << perc_gaps << "\n";
  }
  cout << " -----------\n\n";

  
  vec<int> gapped_lengths, ungapped_lengths;
  get_gapped_lengths( scaffolds, gapped_lengths);
  get_ungapped_lengths( scaffolds, ungapped_lengths );

  vec<int> gapped_nstats, ungapped_nstats;
  get_nstats( NSels, gapped_lengths, gapped_nstats );
  get_nstats( NSels, ungapped_lengths, ungapped_nstats );
  


  // compute good gapped segments of scaffolds

  vec< vec<look_align> > good_segs;
  vec<look_align> current_seg;
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    current_seg.clear(); // each scaffold produces new segment
    if ( scaffolds[si].Ntigs() <= 0 ){
      FatalErr( "require at least one contig in a scaffold!" );
    }
    else if ( scaffolds[si].Ntigs() == 1 ){
      int tid = scaffolds[si].Tig(0);
      if ( aligns_index.at(tid).size() == 1 ){
	look_align la = aligns.at( aligns_index[tid][0] ); 
	current_seg.push_back(la);
	good_segs.push_back( current_seg ); // push single new
      }
    }else{
      PRINT2( si, scaffolds[si].Ntigs() );
      vec<Bool> used( scaffolds[si].Ntigs(), False );
      // find begin candidate
      for ( int tbi = 0; tbi < scaffolds[si].Ntigs(); tbi++ ){
	if ( used[tbi] ) continue;
	int tid0 = scaffolds[si].Tig(tbi);
	if ( aligns_index.at(tid0).size() != 1 ){
	  used[tbi] = True;
	  continue;
	}
	else{
	  used[tbi] = True;
	  look_align la0 = aligns.at( aligns_index[tid0][0] );
	  int targetId0 = la0.TargetId();
	  current_seg.clear();
	  current_seg.push_back( la0 );
	  // find end;
	  int tei = tbi;
	  for ( int ci = tbi + 1; ci < scaffolds[si].Ntigs(); ci++ ){
	    int tid2 = scaffolds[si].Tig(ci);
	    if ( aligns_index.at(tid2).size() != 1 ){
	      tei = ci -1;
	      cout << " because multiple alingments\n";
	      break;
	    }else{
	      int tid1 = scaffolds[si].Tig(ci-1);
	      look_align la1 = aligns.at( aligns_index.at(tid1)[0] );
	      look_align la2 = aligns.at( aligns_index.at(tid2)[0] );
	      if ( la2.TargetId() != targetId0 || la1.Fw1() != la2.Fw1() ){
		tei = ci - 1;
		cout << " because inversion\n";
		break;
	      }else{
		int sep = 0, extension = 0;
		if ( la1.Fw1() && la2.Fw1() ){
		  sep = la2.StartOnTarget() - la1.EndOnTarget();
		  extension = la2.EndOnTarget() - la1.EndOnTarget();
		}else if ( la1.Rc1() && la2.Rc1() ){
		  sep = la1.StartOnTarget() - la2.EndOnTarget();
		  extension = la1.StartOnTarget() - la2.StartOnTarget();
		}
		
		if ( (sep < -max_overlap && extension < 0 ) || sep > max_sep ){
		  tei = ci -1;
		  cout << " because misjoin\n";
		  break;
		}else tei = ci;
	      }
	    }
	  }  
	  PRINT3( tbi, tei, tei - tbi + 1 );
	  if ( tei > tbi ){
	    for ( int j = tbi + 1; j <= tei; j++  ){
	      int tid = scaffolds[si].Tig(j);
	      look_align la = aligns.at( aligns_index[tid][0] );
	      current_seg.push_back( la );
	      used[j] = True;
	    }
	  } 
	}
	good_segs.push_back( current_seg );
	current_seg.clear();
      }
    }  
  }
  
  cout << "Ngood_segments= " << good_segs.size() << endl;

  vec<int> good_segs_lengths;
  vec< vec<ho_interval> > targetHos( ref.size() );
  for ( size_t gsi = 0; gsi < good_segs.size(); gsi++ ){
    if ( good_segs[gsi].size() < 1 ) continue;
    int targetId = good_segs[gsi][0].TargetId();
    int targetLen = ref[targetId].isize();
    vec<ssize_t> tags;
    for ( size_t i = 0; i < good_segs[gsi].size(); i++ ){
      const look_align & la = good_segs[gsi][i];
      ForceAssertEq( la.TargetId(), targetId );
      tags.push_back( la.StartOnTarget() );
      tags.push_back( la.EndOnTarget() );
    }
    int max = Max(tags) < targetLen ? Max(tags) : targetLen;
    int min = Min(tags) >= 0 ? Min(tags) : 0;
    int len = max - min;
    ForceAssertGt( len, 0 );
    good_segs_lengths.push_back(len);
    targetHos.at(targetId).push_back( ho_interval(min, max) );
  }


  size_t good_overlap_sum = 0, good_total_covered = 0;
  for ( size_t ri = 0; ri < targetHos.size(); ri++ ){
    if ( targetHos[ri].size() < 1 ) continue;
    good_total_covered += TotalCovered( targetHos[ri] );
    for ( int hi1 = 0; hi1 < targetHos[ri].isize()-1; hi1++ ){
      for ( int hi2 = hi1 +1; hi2 < targetHos[ri].isize(); hi2++ ){
	if ( Overlap( targetHos[ri][hi1], targetHos[ri][hi2] ) > max_overlap )
	  good_overlap_sum += Overlap( targetHos[ri][hi1], targetHos[ri][hi2] ) - max_overlap;
      }
    }
  }
  PRINT( good_overlap_sum );
  PRINT( good_total_covered );
  cout << endl;

  ReverseSort( good_segs_lengths );
  //PRINT3( good_segs_lengths.size(), good_segs_lengths.front(), good_segs_lengths.back() );
  //cout << "good_segs_lengths: "; good_segs_lengths.Print( cout ); cout << endl;
  size_t total_good_gapped_length = Sum( good_segs_lengths );
  PRINT( total_good_gapped_length );

  vec<int> good_gapped_nstats;
  get_nstats( NSels, good_segs_lengths, good_gapped_nstats );
 
  cout << "Number of scaffolds = " << scaffolds.size() << endl;
  cout << "Number of contigs   = " << contigs.size() << endl;

  cout << "\n----------- N-stats\n";
  cout.setf( ios::right);
  cout << setw(6) << "     " << setw(16) << "ungapped" << setw(16) << "gapped" << setw(16) << "good_gapped" << endl;
  cout << setw(6) << "total" << setw(16) << Sum(ungapped_lengths) << setw(16) << Sum(gapped_lengths) << setw(16) << Sum(good_segs_lengths) << endl;
  for ( size_t isel = 0; isel < NSels.size(); isel++ )
    cout << setw(6) << NSels[isel] << setw(16) << ungapped_nstats[isel] << setw(16) << gapped_nstats[isel] << setw(16) << good_gapped_nstats[isel] << endl;
  cout << " ------------\n";

  cout << Date( ) << ": done" << endl;
  
}

