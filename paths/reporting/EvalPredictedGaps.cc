///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "math/NStatsTools.h"
#include "FastIfstream.h"
#include "paths/reporting/ReftigUtils.h"

/**
 * EvalPredictedGaps 
 *
 * evaluate gaps predicted by FindUnipathGaps
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
  CommandArgument_String( HEAD );
  CommandArgument_String_OrDefault( REF, "genome" );
  CommandArgument_Bool_OrDefault( FORCE, True );
  CommandArgument_String_OrDefault( OUT_DIR, "" );
  CommandArgument_Int_OrDefault_Doc( SEARCH_K, 40, 
      "K-mer size used for seeding perfect alignments");
  EndCommandArguments;


  vec<int> NSels{25,50,75,90,99};

  // Dir and file names.


  String seqs_file       = HEAD + ".fastb";     
  String qlt_file        = HEAD + ".qlt";
  String ref_file        = REF + ".fastb";
  String ref_lookup_file = REF + ".lookup";

  String pred_gap_file   = HEAD + ".predicted_gaps.txt";
  
  String out_dir = (OUT_DIR == "" ? HEAD + ".EvalPredGaps" : OUT_DIR);
  Mkpath( out_dir );

  vecbvec seqs( seqs_file );
  vecbvec ref( ref_file );

  // computing sequence length stats
  

  vec<int> seqsLens( seqs.size(), 0 );
  for ( size_t is = 0; is < seqs.size(); is++ )
    seqsLens[is] = seqs[is].size();
  vec<int> seqsLens_nstats;
  get_nstats( NSels, seqsLens, seqsLens_nstats);
  cout << "Input sequence lengths stats ----------" << endl;
  cout.setf( ios::right);
  cout << setw(6) << "Nstat" << setw(12) << "value" << endl;
  cout << setw(6) << "total" << setw(12) << Sum(seqsLens) << endl;
  for ( size_t isel = 0; isel < NSels.size(); isel++ )
    cout << setw(6) << NSels[isel] << setw(12) << seqsLens_nstats[isel] << endl;
  cout << " ------------\n";


  vec<look_align> aligns;
  if ( FORCE || ! IsRegularFile( qlt_file ) ) {
    cout << Date( ) << ": aligning sequences" << endl;
    // Note that this function uses OMP (badly)
    GetAlignsFast( SEARCH_K, seqs_file, ref_lookup_file, qlt_file, aligns, !FORCE, out_dir ); 
  }else{
    cout << Date() << ": loading sequence alignments" << endl;
    LoadLookAligns( qlt_file, aligns );
  }
  cout << Date() << ": found " << aligns.size() << " alignments" << endl;

  vec< vec< size_t > > aligns_index( seqs.size() );
  for ( size_t i = 0; i < aligns.size(); i++ ){
    size_t qid = aligns[i].QueryId();
    aligns_index.at(qid).push_back(i);
  }
  
 
  // find high copy number unipaths

  
  vec<int> trueCN( seqs.size(), 0 );
  for ( unsigned si = 0; si < seqs.size(); si++ ){
    if ( aligns_index.at(si).size() > 1 ){
      vec<look_align> las;
      for ( unsigned i = 0; i < aligns_index[si].size(); i++ )
	  las.push_back( aligns[ aligns_index[si][i] ] );
      trueCN[si] = get_CN( las, ref );
      
    }else{
      trueCN[si] = aligns_index[si].isize();
    }
  }

  size_t nCNg1 = 0;
  map<size_t,size_t> cnHisto;
  for ( size_t si = 0; si < trueCN.size(); si++ )
    cnHisto[ trueCN[si] ]++;
  cout << " --------------- Copy Number evaluation\n";
  cout.setf( ios::right);
  cout << setw(4) << "cn" << setw(10) << "nseqs" << "\n";
  for ( map<size_t,size_t>::const_iterator im = cnHisto.begin(); im != cnHisto.end(); im++ ){
    cout << setw(4) << im->first << setw(10) << im->second << "\n";
    if ( im->first > 1 ) nCNg1 += im->second;
  }
  cout << "\n";
  cout << " total number CN>1 = " << nCNg1 << "\n";
  cout << " total number of sequences = " << seqs.size() << "\n";
  cout << "----------------\n\n";



  
  // evalulate predicted gaps between unipaths

  vec< vec<ho_interval> > targetHos( ref.size() );
  vec<double> zs;
  fast_ifstream in( pred_gap_file);
  String line;
  while(1) {    
    getline( in, line );
    if ( in.fail( ) ) break;
    if ( line.Contains( "#", 0 ) ) continue;
    int u1, u2, psep, pdev, nlinks;
    istringstream iline( line.c_str( ) );
    iline >> u1 >> u2 >> psep >> pdev >> nlinks;
    if ( trueCN.at(u1) > 1 || trueCN.at(u2) > 1 ) continue;
    if ( trueCN[u1] == 0 || trueCN[u2] == 0 ) continue;
    look_align & la1 = aligns.at( aligns_index[u1][0] );
    look_align & la2 = aligns.at( aligns_index[u2][0] );
    if ( la1.TargetId() != la2.TargetId() ) continue;
    vec<int> points;
    points.push_back( la1.StartOnTarget(), la1.EndOnTarget(), 
		      la2.StartOnTarget(), la2.EndOnTarget() );
    Sort( points );
    targetHos.at( la1.TargetId() ).push_back( ho_interval( points.front() +1, points.back() -1 ) );
  
    double tsep = 0.31415, extension = 0.31415; // initialization
    if ( la1.Fw1() && la2.Fw1() ){
      tsep       = la2.StartOnTarget() - la1.EndOnTarget();
      extension  = la2.EndOnTarget() - la1.EndOnTarget();
    }
    else if ( la1.Rc1() && la2.Rc1() ){
      tsep      = la1.StartOnTarget() - la2.EndOnTarget();
      extension = la1.StartOnTarget() - la2.StartOnTarget();
    }
    
    //double z = 1000000.0;
    if ( tsep != 0.31415 ){
      double z = (double) abs( tsep -psep ) / (double)pdev;
      zs.push_back( z );
    }
  }  
  vec<int> un_segs, cov_segs;
  for ( int ti = 0; ti < targetHos.isize(); ti++ ){
    vec<ho_interval> un, cov;
    Uncovered( 1, targetHos[ti], un );
    ExtractGivenCoverage( ref.at(ti).size(), 1, targetHos[ti], cov );
    for ( int ui = 0; ui < un.isize(); ui++ )
      un_segs.push_back( un[ui].Stop() - un[ui].Start() );
    for ( int ci = 0; ci < cov.isize(); ci++ )
      cov_segs.push_back( cov[ci].Stop() - cov[ci].Start() );
  }
  PRINT3( un_segs.size(), Min(un_segs), Max(un_segs) );
  PRINT3( cov_segs.size(), Min(cov_segs), Max(cov_segs) );
  vec<int> un_nstats, cov_nstats;
  get_nstats( NSels, un_segs, un_nstats );
  get_nstats( NSels, cov_segs, cov_nstats );
  cout << "Coverage in predicted links between sequences:" << endl;
  cout.setf( ios::right);
  cout << setw(6) << "Nstat" << setw(12) << "uncovered" << setw(12) << "covered" << endl;
  cout << setw(6) << "total" << setw(12) << Sum(un_segs) << setw(12) << Sum(cov_segs) << endl;
  for ( size_t isel = 0; isel < NSels.size(); isel++ )
    cout << setw(6) << NSels[isel] << setw(12) << un_nstats[isel] << setw(12) << cov_nstats[isel] << endl;
  cout << " ------------\n";

  

  
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

  
}

