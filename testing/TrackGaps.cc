///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// TrackGaps
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
#include "feudal/BinaryStream.h"
#include "feudal/VirtualMasterVec.h"
#include "paths/Alignlet.h"
#include "SupersHandler.h"
#include "PairsManager.h"
#include "paths/PairDistCorrection.h"


int main(int argc, char *argv[])
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(FULL_RUN);
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String(ASSEMBLY);
  CommandArgument_String(ALIGNS);
  CommandArgument_String(GALIGNS);
  CommandArgument_String(READS);
  CommandArgument_Int(C1);
  CommandArgument_Int_OrDefault(C2, -1);
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  CommandArgument_String_OrDefault( DIST, "" );
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
  }

  int s1 = supers.ToSuper(cg1);
  int s2 = supers.ToSuper(cg2);
  ForceAssertEq(s1,s2);
  
  int p1 = supers.PosOnSuper(cg1);
  int p2 = supers.PosOnSuper(cg2);
  ForceAssertEq(p2-p1, 1);

  superb super = supers[s1];
  int len1 = super.Len(p1);
  int len2 = super.Len(p2);

  // tig length and connection on supers 
  map<int,int> next_tigs;

  int gap = super.Gap(p1);
  cout << "gap between tigs " << cg1 << " " << cg2 << " is " << gap << endl;

  PairsManager pairs(pairs_file);

  vec<alignlet> aligns;
  BinaryReader::readFile(alignlet_file,&aligns);
  vec<int> index;
  BinaryReader::readFile(alignlet_index_file,&index);

  vec<alignlet> galigns;
  BinaryReader::readFile(galignlet_file,&galigns);
  vec<int> gindex;
  BinaryReader::readFile(galignlet_index_file,&gindex);


  uint nLibs = pairs.nLibraries();
  vec<int> libSeps(nLibs, 0);
  vec<int> libSDs(nLibs, 0);
  for(size_t i=0;i<nLibs;i++)
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
    //pfids.push( ProbFuncIntDist( distribs[i] ) ); // Invaraint pair separation
    pfids.push(ProbFuncIntDist(distribs[i],  -2*stats[i].mean_len)); // reads end-to-end separation

  vec< vec< pair<int,int> > > all_links(nLibs);
  for(size_t i=0;i<pairs.nPairs();i++)
  {
    int sep = pairs.sep(i);
    int sd = pairs.sd(i);
    int libID = pairs.libraryID(i);
    size_t rid1 = pairs.ID1(i);
    size_t rid2 = pairs.ID2(i);
    int aid1 = index[rid1];
    int aid2 = index[rid2];
    if (aid1 < 0 || aid2 < 0) continue;
    int cg1Target = aligns[aid1].TargetId();  
    int cg2Target = aligns[aid2].TargetId(); 
    if ( ! ( cg1Target == cg1 && cg2Target == cg2  ||
	     cg1Target == cg2 && cg2Target == cg1 ) ) continue;
    if (cg1Target != cg1 ) swap(cg1Target, cg2Target), swap(aid1,aid2);   // make sure links from cg1 to cg2
    const alignlet & a1 = aligns[aid1];
    const alignlet & a2 = aligns[aid2];

    int dist_to_end1, dist_to_end2;
    if (a1.Fw1()) 
      dist_to_end1 = len1 - a1.Pos2();
    else
      dist_to_end1 = a1.pos2();
    if (a2.Fw1()) 
      dist_to_end2 = len2 - a2.Pos2();
    else 
      dist_to_end2 = a2.pos2();

    /* !!!!!!!!!!! filters !!!!!!!!!!!!!!!!!!!!*/
    if ( dist_to_end1 + dist_to_end2 > sep + sd * 5 + 1000) continue;
    if ( ! a1.Fw1() || a2.Fw1() ) continue;

    all_links[libID].push( make_pair( dist_to_end1, dist_to_end2) );

    // check alignments on the genome
    int real_sep = -1000;
    int gaid1 = gindex[rid1];
    int gaid2 = gindex[rid2];
    if ( gaid1 < 0 || gaid2 < 0 ) real_sep = -500; 
    else {
      String dir1 = galigns[gaid1].Fw1() ? "+" : "-" ; 
      String dir2 = galigns[gaid2].Fw1() ? "+" : "-" ; 
      int linkPos1 =  galigns[gaid1].Fw1() ? galigns[gaid1].Pos2() : galigns[gaid1].pos2();
      int linkPos2 =  galigns[gaid2].Fw1() ? galigns[gaid2].Pos2() : galigns[gaid2].pos2();
      if( galigns[gaid1].Fw1())
	real_sep = linkPos2 - linkPos1;
      else
	real_sep = linkPos1 - linkPos2;
      real_sep = min(10000, real_sep);
      real_sep = max(-1000, real_sep);
    }
    cout <<"link "<< libID  << " " << dist_to_end1 <<" "<<dist_to_end2<<" "<< real_sep << " " << real_sep - dist_to_end2 - dist_to_end1<< endl;
  }
  cout << "len1= " << len1 << endl;
  cout << "len2= " << len2 << endl;

  int gap_eval, std_eval;
  GapByFitKSCombine( pfids, len1, len2, all_links, gap_eval, std_eval );
  cout << "gap_eval= " << gap_eval << endl;
  for(size_t i=0;i<nLibs;i++)
  {
    int gap_i, std_i;
    GapByFitKS( pfids[i], len1, len2, all_links[i], gap_i, std_i );
    cout << "lib" << i << "#" << all_links[i].size() << " " << gap_i << endl;
  }
}
