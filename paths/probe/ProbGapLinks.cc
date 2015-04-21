///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// /ProbGapLinks
// extract the scaffold read pairs that link the gaps in scaffolds
// 
#include "Equiv.h"
#include "MainTools.h"
#include "CoreTools.h"
#include "Superb.h"
#include "math/Functions.h"
#include "system/ParsedArgs.h"
#include "util/ReadTracker.h"
#include "paths/Alignlet.h"
#include "paths/MapAlignletToRef.h"
#include "feudal/VirtualMasterVec.h"
#include "paths/Alignlet.h"
#include "SupersHandler.h"
#include "PairsManager.h"
#include "paths/PairDistCorrection.h"
#include "paths/PairDistFitting.h"
#include "math/Array.h"
#include "paths/LinkingPairs.h"
#include "feudal/BinaryStream.h"


void DumpLinkInfo( ostream& out, const vec<ProbFuncIntDist>& distrs, int gap_real, 
    const vec<int> &cutoffs, const vec< vec<int> > & x1s, const vec< vec<int> > & x2s,
   const vec< vec< pair<int,int> > > & links ) 
{
  int nlibs = links.size();
  ForceAssertEq( nlibs, x1s.isize() );
  ForceAssertEq( nlibs, x2s.isize() );
  ForceAssertEq( nlibs, distrs.isize() );
  out << "gap_real= " << gap_real << endl;
  out << "cutoffs= " << cutoffs << endl;
  vec< vec<double> > x1s_smooth(nlibs);
  vec< vec<double> > x2s_smooth(nlibs);
  for ( int libid = 0; libid < nlibs; libid++ ) {
    x1s_smooth[libid].resize( x1s[libid].size() );
    x2s_smooth[libid].resize( x2s[libid].size() );
    for( int i = 0; i < x1s[libid].isize(); i++ ) x1s_smooth[libid][i] = x1s[libid][i];
    for( int i = 0; i < x2s[libid].isize(); i++ ) x2s_smooth[libid][i] = x2s[libid][i];
    int delta = 50;
    SmoothArrayGaussian( x1s_smooth[libid], delta );
    SmoothArrayGaussian( x2s_smooth[libid], delta );
  }

  for ( int libid = 0; libid < nlibs; libid++ )
    out << "pairs"<<libid<< " " <<  distrs[libid].Avg() * sqrt(2.0) << " " 
      << distrs[libid].SD() * sqrt(2.0) << endl;
}


void PrepairLinks( const PairsManager& pairs, const vec<alignlet>& aligns, const vec<int>& index, 
    const vec<int>& cutoffs,
      map< int, vec< vec< int > > > &xs, 
      map<int, map<int, vec< vec< pair<int, int> > > > > &link_maps,
      map<int, map<int, vec< vec< longlong > > > > &pair_maps, bool verbose = false)
{
  int nLib = cutoffs.size();
  int counter=0;
  for( size_t iPair=0; iPair < pairs.nPairs(); iPair++ ) 
  {
    int aid1 = index[pairs.ID1(iPair)];
    int aid2 = index[pairs.ID2(iPair)];
    int libID = pairs.libraryID(iPair);
    // the two reads of the pair
    vec <int> aids;
    aids.push_back(aid1);
    aids.push_back(aid2);
    vec <int> tigs(2, -1);
    vec <int> tindex(2, -1);
    vec <int> dist_to_end(2, -1);
    // pairs sorting
    for( size_t kk = 0; kk < 2; kk++ ){
      int aid = aids[kk];
      if( aid < 0) continue;
      const alignlet & a = aligns[aid];
      int cg = a.TargetId();
      int len = a.TargetLength();
      int cutoff_len = Min( len, cutoffs[libID] ); // cutoff length
      int d_to_end = a.Fw1() ? len - a.Pos2(): a.pos2(); // distance to the end, in [0, len)
      if ( d_to_end >= cutoff_len ) continue;
      // for read located within cutoff distance from the edge only
      dist_to_end[kk] = d_to_end;
      tindex[kk] =  a.Fw1() ? 2 * cg : 2 * cg + 1;
      xs[ tindex[kk] ].resize(nLib);
      xs[ tindex[kk] ][libID].resize( cutoff_len, 0 );
      xs[ tindex[kk] ][libID][ d_to_end ]++;
      tigs[kk] = cg; // save the information
    }
    // Only deal with the valid pair
    if ( tigs[0] < 0 || tigs[1] < 0 || tigs[0] == tigs[1] ) continue;
    if ( dist_to_end[0] + dist_to_end[1] >= cutoffs[libID] ) continue;
    // save links from index1 to index2, where index1 < index2
    if ( tindex[0] < tindex[1] ) {
      link_maps[tindex[0]][tindex[1]].resize(nLib);
      link_maps[tindex[0]][tindex[1]][libID].push( 
	  make_pair(dist_to_end[0], dist_to_end[1]) );
      pair_maps[tindex[0]][tindex[1]].resize(nLib);
      pair_maps[tindex[0]][tindex[1]][libID].push(iPair);
    }else{
      link_maps[tindex[1]][tindex[0]].resize(nLib);
      link_maps[tindex[1]][tindex[0]][libID].push( 
	  make_pair(dist_to_end[1], dist_to_end[0]) );
      pair_maps[tindex[1]][tindex[0]].resize(nLib);
      pair_maps[tindex[1]][tindex[0]][libID].push(iPair);
    }
    counter++;
    if (verbose) cout << "            number of scaffold links= " << counter<< endl;
  }
}

int main(int argc, char *argv[])
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(FULL_RUN);
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String(ASSEMBLY);
  CommandArgument_String(ALIGNS);
  CommandArgument_String_OrDefault(GALIGNS, ""); // reads aligned to reference genome. 
  CommandArgument_String(READS);
  CommandArgument_Int(C1);
  CommandArgument_Int_OrDefault(C2, -1); // if not set, find the next contig following C1 in the scaffolds
					 // if set, don't assume C2 and C1 connected in the scaffold
  CommandArgument_Bool_OrDefault(FW1, True);
  CommandArgument_Bool_OrDefault(FW2, False);
  CommandArgument_Int_OrDefault(REAL_GAP, -10000); // real gap size
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  EndCommandArguments;

  // The input files
  String sub_dir       = FULL_RUN + "/ASSEMBLIES/" + SUBDIR;
  String superb_file   = sub_dir + "/" + ASSEMBLY + ".superb" ;
  String contig_file   = sub_dir + "/" + ASSEMBLY + ".contigs.fastb" ;
  String alignlet_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
  String alignlet_index_file = alignlet_file + ".index";
  String galignlet_file = sub_dir + "/" + GALIGNS + ".qltoutlet";
  String galignlet_index_file = galignlet_file + ".index";
  String reads_file    = FULL_RUN + "/"+ READS + ".fastb";
  String pairs_file    = FULL_RUN + "/"+ READS + ".pairs";

  // read the scaffolds and derive the next contigs if necessary
  shandler supers( -1, superb_file );
  int cg1 = C1, cg2 = C2;
  if ( cg2 < 0 ) {
    for( int s=0; s<supers.Size(); s++) {
      superb super = supers[s];
      for(int p=0; p<super.Ntigs(); p++) {
	if ( super.Tig(p) == cg1 ) { 
	  cg2 = super.Tig( p+1 );
	  break;
	}
      }
    }
    int s1 = supers.ToSuper(cg1);
    int s2 = supers.ToSuper(cg2);
    int p1 = supers.PosOnSuper(cg1);
    int p2 = supers.PosOnSuper(cg2);
    superb super1 = supers[s1];
    superb super2 = supers[s1];

    if ( C2 < 0 ) {
      int gap = super1.Gap(p1);
      cout << "estimated gap between tigs " << cg1 << " " << cg2 << " is " << gap << endl;
    }
  }

  int len1 = -1;
  int len2 = -1;
  PairsManager pairs(pairs_file);

  vec<alignlet> aligns;
  BinaryReader::readFile(alignlet_file,&aligns);
  vec<int> index;
  BinaryReader::readFile(alignlet_index_file,&index);

  vec<alignlet> galigns;
  vec<int> gindex;
  if ( GALIGNS != "" ) {
    BinaryReader::readFile(galignlet_file,&galigns);
    BinaryReader::readFile(galignlet_index_file,&gindex);
  }

  uint nLib = pairs.nLibraries();
  vec<int> libSeps(nLib, 0);
  vec<int> libSDs(nLib, 0);
  for(size_t i=0;i<nLib;i++)
  {
    libSeps[i] = pairs.getLibrarySep(i);
    libSDs[i] = pairs.getLibrarySD(i);
  }

  // the distribution file
  String reads_head = FULL_RUN + "/" + READS;
  vec<PM_LibraryStats> stats = pairs.getLibraryStats(reads_head + ".fastb");
  String file =  reads_head + ".distribs" ; //  for distribution
  vec< IntDistribution > distribs;
  BinaryReader::readFile( file.c_str(), &distribs );
  vec<ProbFuncIntDist> pfids;  // vec of probability function of IntDistributions
  for( size_t i=0; i<distribs.size(); i++ )
    pfids.push(ProbFuncIntDist(distribs[i], distribs[i].x_min(), distribs[i].x_max(), false, -2*stats[i].mean_len)); // reads end-to-end separation
    //pfids.push( ProbFuncIntDist( distribs[i] ) ); // Invaraint pair separation

  int MIN_N_LINKS = 5; // minimal number of links required to apply the correction
  int MAX_OVERLAP = 2000;

  vec< int > cutoffs(nLib); // cutoff for dist_to_end1 + dist_to_end2;
  for( size_t libId = 0; libId < nLib; libId++ ) {
    int libSep = pairs.getLibrarySep(libId);
    int libSD = pairs.getLibrarySD(libId);
    cutoffs[libId] = libSep + libSD * 5 + MAX_OVERLAP; 
    cout << "pairs" << libId << " " << pfids[libId].Avg()*sqrt(2.0) << " " <<  pfids[libId].SD() * sqrt(2.0)  << endl;
  }

  // xs[index1][libId][pos] coverage counter
  map< int, vec< vec< int > > > xs;  
  // link_maps[index1][index2][libID] is an array of pairs of links
  typedef vec< pair<int, int> > Pairs;
  map<int, map<int, vec< vec< pair<int, int> > > > > link_maps;  
  map<int, map<int, vec< vec< longlong > > > > pair_maps;  

  cout << Date() << ": Start sorting" << endl;
  PrepairLinks( pairs, aligns, index, cutoffs, xs, link_maps, pair_maps, false );
  cout << Date() << ": End sorting" << endl;

  // let look at the selected pairs between cg1 and cg2, fw cg1 and rc on cg2
  int index1 = FW1? cg1 * 2 : cg1 * 2 + 1; 
  int index2 = FW2? cg2 * 2 : cg2 * 2 + 1;
  if ( index1 > index2 ) swap( index1, index2 ); // force index1 smaller than index2, the way links are saved
  vec< vec< pair<int,int> > >& all_links = link_maps[index1][index2];
  vec< vec< longlong > >& pair_ids = pair_maps[index1][index2];

  for( size_t libId = 0; libId < pair_ids.size(); libId++ ) {
    for( size_t kk = 0; kk < pair_ids[libId].size(); kk ++) {
      longlong pairId = pair_ids[libId][kk];
      longlong rid1 = pairs.ID1(pairId);
      longlong rid2 = pairs.ID2(pairId);
      // check the original alignment
      alignlet& a1 = aligns[index[rid1]];
      alignlet& a2 = aligns[index[rid2]];
      // do some convertion so that rid1 is from cg1 and rid2 is from cg2
      if ( a1.TargetId() == cg2 && a2.TargetId() == cg1 )
	swap(rid1, rid2);
      else if ( a1.TargetId() == cg1 && a2.TargetId() == cg2 )
	; // do nothing
      else
	exit(1);
      if ( cg1 > cg2 ) swap( rid1, rid2 );
      // check alignments on the genome
      int real_sep = -1000;
      int gaid1 = gindex[rid1];
      int gaid2 = gindex[rid2];
      if ( gaid1 < 0 || gaid2 < 0 ) real_sep = -10000; 
      else {
	String dir1 = galigns[gaid1].Fw1() ? "+" : "-" ; 
	String dir2 = galigns[gaid2].Fw1() ? "+" : "-" ; 
	// invariant seprations
	int linkPos1 =  galigns[gaid1].Fw1() ? galigns[gaid1].Pos2() : galigns[gaid1].pos2();
	int linkPos2 =  galigns[gaid2].Fw1() ? galigns[gaid2].Pos2() : galigns[gaid2].pos2();
	if( galigns[gaid1].Fw1())
	  real_sep = linkPos2 - linkPos1;
	else
	  real_sep = linkPos1 - linkPos2;
	//real_sep = min(10000, real_sep);
	//real_sep = max(-1000, real_sep);
      }
      int dist_to_end1 = all_links[libId][kk].first;
      int dist_to_end2 = all_links[libId][kk].second;
      // do some convertion so that rid1 is from cg1 and rid2 is from cg2
      if ( cg1 > cg2 ) swap( dist_to_end1, dist_to_end2 );
      cout <<"link "<< libId  << " " << dist_to_end1 <<" "<<dist_to_end2<<" "<< real_sep << " " << real_sep - dist_to_end2 - dist_to_end1<< endl;
    }
  }
  // coverage over the two contigs , smoothed using 100-base window
  cout << "len1= " << len1 << endl;
  cout << "len2= " << len2 << endl;
  vec< vec<int> >& x1 = xs[ index1 ];
  vec< vec<int> >& x2 = xs[ index2 ];

  int gap_real = REAL_GAP;
  DumpLinkInfo(cout, pfids, gap_real, cutoffs, x1, x2, all_links);
}
