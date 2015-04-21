///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "MainTools.h"
#include "FetchReads.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Alignlet.h"
#include "FastaFileset.h"

#include <map>

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable

Bool cmp_target( const look_align& la1, const look_align& la2 ){    
  if ( la1.target_id < la2.target_id ) return True;
  if ( la1.target_id > la2.target_id ) return False;
  if ( la1.pos2( ) < la2.pos2( ) ) return True;
  if ( la1.pos2( ) > la2.pos2( ) ) return False;
  if ( la1.query_id < la2.query_id ) return True;
  return False;    
}


void GetAligns( const int K, const String& fastb_file, const String& lookup_file, 
		const String& aligns_file, vec<look_align>& aligns,
		const Bool USE_CACHE, const String& tmp_dir ){    
  if ( !USE_CACHE || !IsRegularFile(aligns_file) ) {    
    vec< triple<int64_t,int64_t,int> > ALIGNS;
    String gfastb = lookup_file.Before( ".lookup" ) + ".fastb";
    if ( !IsRegularFile(fastb_file) ){  
      vecbasevector seqs;
      FetchReads( seqs, 0, fastb_file.Before( ".fastb" ) + ".fasta" );
      seqs.WriteAll(fastb_file);    
    }
    vecbasevector seqs(fastb_file);
    int nproc = omp_get_max_threads( ), N = seqs.size( );
    for ( int j = 0; j < nproc; j++ ){    
      Ofstream( idsout, aligns_file + ".ids." + ToString(j) );
      for ( int i = 0; i < N; i++ )
	if ( i % nproc == j ) idsout << i << "\n";    
    }
    #pragma omp parallel for
    for ( int j = 0; j < nproc; j++ ){    
      SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.3 MF=5000 MAX_ERROR_PERCENT=2 " 
		     + ARG(SEQS, fastb_file) 
		     + ARG(SEQS_IS_FASTB, True) 
		     + ARG(L, lookup_file) + ARG(PARSEABLE, True) 
		     + ARG(SMITH_WAT, True) + ARG(TMP_DIR,tmp_dir)
		     + ARG(SEQS_TO_PROCESS, "@" + aligns_file + ".ids." + ToString(j))
		     + ARG(FW_ONLY, True )
		     + " > " + aligns_file + "." + ToString(j) );  

      SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.3 MF=5000 MAX_ERROR_PERCENT=2 " 
		     + ARG(SEQS, fastb_file) 
		     + ARG(SEQS_IS_FASTB, True) 
		     + ARG(L, lookup_file) + ARG(PARSEABLE, True) 
		     + ARG(SMITH_WAT, True) + ARG(TMP_DIR,tmp_dir)
		     + ARG(SEQS_TO_PROCESS, "@" + aligns_file + ".ids." + ToString(j))
		     + ARG(RC_ONLY, True )
		     + " >> " + aligns_file + "." + ToString(j) );  
    }
    {    
      Ofstream( out, aligns_file );
      for ( int j = 0; j < nproc; j++ ){    
	vec<look_align> A;
	LoadLookAligns( aligns_file + "." + ToString(j), A );
	for ( int i = 0; i < A.isize( ); i++ )
	  A[i].PrintParseable(out);
	Remove( aligns_file + ".ids." + ToString(j) );
	Remove( aligns_file + "." + ToString(j) );    
      }    
    }
  }
  LoadLookAligns( aligns_file, aligns );    
  sort( aligns.begin( ), aligns.end( ), cmp_target );    
}




/**
 * ProbeInversions
 *
 
*/


int main( int argc, char *argv[] ){
  
  RunTime( );
  
  BeginCommandArguments;

  CommandArgument_String_OrDefault( PRE, "" );
  CommandArgument_String( ASSEMBLY );
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_String_OrDefault( REF, "" );
  CommandArgument_Int_OrDefault_Doc( WINDOW, 100000, 
    "maximum size of segment in which to check for local inversion occurance" );
  CommandArgument_Bool_OrDefault( REBUILD_QUERIES, True );
  CommandArgument_Bool_OrDefault( REBUILD_LOOKUP, False );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  EndCommandArguments;

  // Thread control (OMP used in FixScaffoldsCore)
    
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );


  if ( PRE != "" )
    PRE += "/";

  String assembly_head = ASSEMBLY.SafeAfterLast("/").SafeBefore(".fast");
  PRINT(assembly_head);

  String assembly_file = PRE + ASSEMBLY;

  String work_dir = PRE + "tmpProbeInv4";
  if ( ! IsRegularFile(work_dir) )
    Mkdir777(work_dir);

  // read assembly file 


  cout << Date( ) << ": loading assembly file" << endl;
  vecbasevector aseqs;
  vec<fastavector> afastaseqs;
  if ( ASSEMBLY.Contains(".fastb") ){
    aseqs.ReadAll(assembly_file );
    for ( size_t i = 0; i < aseqs.size(); i++ )
      afastaseqs.push_back( fastavector(aseqs[i]) );
    ForceAssertEq( aseqs.size(), afastaseqs.size() );
  }
  else if ( ASSEMBLY.Contains(".fasta") ){
    vecString readNames;
    FastFetchReads(aseqs, &readNames, assembly_file );
    LoadFromFastaFile( assembly_file, afastaseqs );
    ForceAssertEq( aseqs.size(), afastaseqs.size() );
  }else
    FatalErr("I don't know how to read ASSEMBLY file=" + ASSEMBLY );      
  cout << "Found " << aseqs.size() << " assembly sequences" << endl; 


  // find local inversions in assembly
  String inversions_file = work_dir + "/" + assembly_head + ".inversions.fastb";
  String invregs_file = work_dir + "/" + assembly_head + ".invregs.fasta";
  vecbasevector invseqs, invnhoods;
  vec<fastavector> invseqsregs;
  if ( ! IsRegularFile(inversions_file) || REBUILD_QUERIES ){
    cout << Date() << ": pathing" << endl;
    vecKmerPath paths, paths_rc;
    vec<tagged_rpint> pathsdb;
    ReadsToPathsCoreY( aseqs, K, paths, paths_rc, pathsdb, "", NUM_THREADS);
    
    cout << Date() << ": unipathing" << endl;
    vecKmerPath unipaths;
    vec<tagged_rpint> unipathsdb;
    Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb, True);
    
    PRINT4( aseqs.size(), paths.size(), paths_rc.size(), unipaths.size() );
    
    
    map< int, map< int, pair<int,vec<alignlet> > > > loci;  
    
    for ( size_t u = 0; u < unipaths.size(); u++ ){
      vec<alignlet> faligns;
      KmerPath upath = unipaths[u];
      //PRINT(upath);
      kmer_id_t kf0 = upath.Start();
      //PRINT( kf0 );
      vec<path_interval_id_t> flocs;
      //PRINT( pathsdb.size() );
      Contains( pathsdb, kf0, flocs ); 
      if ( flocs.size() <= 0 )
	cout << "very strange, unipath segment not found in database" << endl;
      // cout << "flocs: "; flocs.Print(cout); cout << endl;
      
      for ( size_t l = 0; l < flocs.size(); l++ ){
	const tagged_rpint& t = pathsdb[flocs[l]];
	int tid = t.PathId();
	if ( tid < 0 ) continue;
	read_id_t rid = t.ReadId();
	longlong offset = kf0 - t.Start();
	//PRINT3( u, tid, offset );
	for ( size_t ikpi = 0; ikpi < t.PathPos(); ikpi++ )
	  offset += paths[tid][ikpi].Length();
	alignlet la( offset, offset + upath.TotalLength() + K -1, rid, paths[tid].TotalLength() + K -1, True );
	faligns.push_back( la );
      }
      
      
      KmerPath urpath = upath;
      urpath.Reverse();
      if ( upath == urpath )
	cout << "FOUND PALINDROMIC UNIPATH of length " << upath.TotalLength() + K -1 << endl;
      kmer_id_t kr0 = urpath.Start();
      vec<path_interval_id_t> rlocs;
      Contains( pathsdb, kr0, rlocs );
      
      vec<alignlet> raligns;
      for ( size_t l = 0; l < rlocs.size(); l++ ){
	const tagged_rpint& t = pathsdb[rlocs[l]];
	int tid = t.PathId();
	if ( tid < 0 ) continue;
	read_id_t rid = t.ReadId();
	longlong offset = kr0 - t.Start();
	for ( size_t ikpi = 0; ikpi < t.PathPos(); ikpi++ )
	  offset += paths[rid][ikpi].Length();
	raligns.push_back( alignlet( offset, offset + urpath.TotalLength() + K -1, rid, paths[tid].TotalLength() + K -1, False ) );
      }
      
      if ( faligns.size() == 0 || raligns.size() == 0 )
	continue;
      
      for ( size_t i = 0; i < faligns.size(); i++ ){
	alignlet& la1 = faligns[i];
	for ( size_t j = 0; j < raligns.size(); j++ ){
	  alignlet& la2 = raligns[j];
	  if ( la1.TargetId() == la2.TargetId() && la1.Fw1() == !la2.Fw1() ){
	    int tid = la1.TargetId();
	    int q0  = la1.pos2(), q1 = la1.Pos2();
	    int t0  = la2.pos2(), t1 = la2.Pos2();
	    int dist = q0 <= t0 ? abs(t1 - q0) : abs(q1 - t0);
	    int center = q0 <= t0 ? q0 + abs(round((double)dist/2.0)) : t0 + abs(round((double)dist/2.0));
	    if ( abs(dist) < WINDOW ){
	      loci[tid][center].first++; 
	      loci[tid][center].second.push_back( la1, la2 );
	    }
	  }  
	}
      }
      
    }
    
    
    for ( map<int,map<int,pair<int,vec<alignlet> > > >::iterator i = loci.begin(); i != loci.end(); i++ ){
      int tid = i->first;
      size_t tlen = aseqs[tid].size();
      vec<ho_interval> intervs;
      for ( map<int,pair<int,vec<alignlet> > >::iterator j = i->second.begin(); j != i->second.end(); j++ ){
	if ( j->second.first  > 0 ){
	  //PRINT4( i->first, j->first, j->second.first, j->second.second.size() );
	  vec<alignlet> las = j->second.second;
	  for ( size_t k = 1; k < las.size(); k++ ){
	    ForceAssertEq( (int)tlen, las[k].TargetLength() );
	    intervs.push_back( ho_interval( las[k].pos2(), las[k].Pos2() ) );
	  }
	}
      
	
	vec<ho_interval> cov;
	ExtractGivenCoverage( tlen, 1, intervs, cov );
	// find minimum edge distance
	vec<int> eds;
	for ( size_t ic = 0; ic < cov.size(); ic++ ){
	  eds.push_back( cov[ic].Start() );
	  eds.push_back( aseqs[tid].size() - cov[ic].Stop() );
	}
	int minEdgeDist = Min(eds);
	//ForceAssertGt( minEdgeDist, 10 );
	

	for ( size_t ic = 0; ic < cov.size(); ic++ ){
	  if ( abs(cov[ic].Stop() - cov[ic].Start()) > K ){
	    invseqs.push_back( basevector( aseqs[tid], cov[ic].Start(), cov[ic].Stop() - cov[ic].Start() ) );
	    invseqs.push_back( invseqs.back() );
	    invseqs.back().ReverseComplement();

	  }
	}		
	vec<int> poss;
	for ( size_t ci = 0; ci < cov.size(); ci++ )
	  if ( abs(cov[ci].Stop() - cov[ci].Start()) > K ){
	    poss.push_back( cov[ci].Start() );
	    poss.push_back( cov[ci].Stop() );
	  }
	if ( poss.size() > 0 ){
	  basevector fv;
	  PRINT2( Min(poss), Max(poss) );
	  fv.SetToSubOf( aseqs[tid], Min(poss), Max(poss)-Min(poss) );
	  invseqsregs.push_back( fastavector(fv) );
	}
      }
    }
    invseqs.UniqueSort();
    Sort( invseqsregs );

    vec<int> toRc;
    UnibaseInvolution( invseqs, toRc );
    
    for ( size_t u = 0; u < invseqs.size(); u++ )
      if ( ! invseqs[u].empty() && toRc[u] != (int)u )
	invseqs[ toRc[u] ].resize(0);
    
    vecbasevector newtmp;
    for ( size_t u = 0; u < invseqs.size(); u++ )
      if ( invseqs[u].size() != 0 )
	newtmp.push_back( invseqs[u] );
    
    invseqs = newtmp;   
  
    invseqs.WriteAll( inversions_file );
    ofstream aout( invregs_file.c_str() );
    for ( size_t regid = 0; regid < invseqsregs.size(); regid++ ){
      invseqsregs[regid].Print(aout, ToString(regid));
    }
    aout.close();
  }else{
    invseqs.ReadAll( inversions_file );
    LoadFromFastaFile( invregs_file, invseqsregs );
  }
  cout << "There are " << invseqs.size() << " candidate inversion sequences" << endl;


  // prepare assembly files for alignments

  String lookup_file1 = work_dir + "/" + assembly_head + ".lookup";
  if ( ! IsRegularFile( lookup_file1 ) || REBUILD_LOOKUP ){
    cout << Date() << ": Creating a lookup table for assembly file" << endl;
    SystemSucceed( "MakeLookupTable" + ARG(SOURCE, assembly_file )
		   + ARG(OUT_HEAD, lookup_file1.SafeBefore(".lookup")) + ARG(LOOKUP_ONLY, True)
		   + ARG( NH, True ) );
  }

  if ( ! assembly_file.Contains(".fastb") && ! IsRegularFile( assembly_file.SafeBefore(".fast") + ".fastb" ) ){
    cout << Date() << ": creating assembly fastb file" << endl;
    SystemSucceed( "Fasta2Fastb" + ARG(IN, assembly_file )
		   + ARG(OUT, lookup_file1.SafeBefore(".lookup") + ".fastb" ) + ARG( NAMES, False ) +
			   ARG( NH, True ) );
  }
  
  cout << Date() << ": aligning candidate inversions to assembly" << endl;
  String inversions_aligns1_file = inversions_file.SafeBefore(".fast") + ".assembly.qltout";
  vec<look_align> inversions_aligns1;
  GetAligns( 12, inversions_file, lookup_file1, inversions_aligns1_file, inversions_aligns1, false, work_dir );  

  vec< vec<int> > ia1_index( invseqs.size() );
  for ( int ai = 0; ai < inversions_aligns1.isize(); ai++ ){
    size_t qid = inversions_aligns1[ai].QueryId();
    ia1_index.at(qid).push_back( ai );
  }


  // if reference specified align candidate inversions and see if they can be confirmed

  if ( REF.nonempty() ){
    String wronginvs_file = work_dir + "/" + assembly_head + ".inversions.wrong.fastb";
    String ref_file = PRE + REF;
    String ref_head = ref_file.SafeAfterLast("/").SafeBefore(".fast");
    
    // prepare reference files for alignment

    String lookup_file2 = work_dir + "/" + ref_head + ".lookup";
    if ( ! IsRegularFile( lookup_file2 ) || REBUILD_LOOKUP ){
      cout << Date() << ": Creating a lookup table for reference file" << endl;
      SystemSucceed( "MakeLookupTable" + ARG(SOURCE, ref_file )
		     + ARG(OUT_HEAD, lookup_file2.SafeBefore(".lookup")) + ARG(LOOKUP_ONLY, True)
		     + ARG( NH, True ) );
    }
    
    if ( ! ref_file.Contains(".fastb") && ! IsRegularFile( ref_file.SafeBefore(".fast") + ".fastb" ) ){
      cout << Date() << ": creating reference fastb file" << endl;
      SystemSucceed( "Fasta2Fastb" + ARG(IN, ref_file )
		     + ARG(OUT, lookup_file2.SafeBefore(".lookup") + ".fastb" ) + ARG( NAMES, False ) +
		     ARG( NH, True ) );
    }

    cout << Date() << ": aligning candidate inversions to reference " << endl;
    String inversions_aligns2_file = inversions_file.SafeBefore(".fast") + ".reference.qltout";
    vec<look_align> inversions_aligns2;

    GetAligns( 12, inversions_file, lookup_file2, inversions_aligns2_file, inversions_aligns2, false, work_dir );
    

    int ninvers = MastervecFileObjectCount( inversions_file );

    vec< vec<int> > ia2_index( ninvers );
    for ( int ai = 0; ai < inversions_aligns2.isize(); ai++ ){
      //inversions_aligns2[ai].PrintReadableBrief( cout );
      size_t qid = inversions_aligns2[ai].QueryId();
      ia2_index.at(qid).push_back( ai );
    }

    vec<Bool> isWrongInver( ninvers, True );
    for ( size_t qid = 0; qid < ia2_index.size(); qid++ ){
      for ( int i1 = 0; i1 < ia2_index[qid].isize(); i1++ ){
	const look_align& la1 = inversions_aligns2[ ia2_index[qid][i1] ];
	for ( int i2 = i1 + 1; i2 < ia2_index[qid].isize(); i2++ ){
	  const look_align& la2 = inversions_aligns2[ ia2_index[qid][i2] ];
	  
	  int len1 = abs( la1.EndOnQuery() - la1.StartOnQuery() );
	  int len2 = abs( la2.EndOnQuery() - la2.StartOnQuery() );
	  //PRINT4( len1, len2, la1.QueryLength(), Max(len1,len2) );
	  if ( (double)Min(len1,len2)/(double)la1.QueryLength() > 0.9 ){
	    isWrongInver[qid] = False;
	  }
	  
	  //int dist = la1.pos2() < la2.pos2() ? abs( la2.Pos2() - la1.pos2() ) : abs( la1.Pos2() - la2.pos2() );
	  //if ( dist < 2.0 * (double)WINDOW ) // 100% distance error allowed
	    //isWrongInver[qid] = False;
	}
      }
    }

    if ( Sum( isWrongInver ) == 0 )
      cout << "Havn't found wrong local inversions" << endl;
    
    cout << "Found " << Sum( isWrongInver ) << " wrong inversions" << endl;

    vecbasevector wronginvs_seqs;
    for ( size_t qid = 0; qid < isWrongInver.size(); qid++ ){
      if ( isWrongInver[qid] ){
	cout << "\nWrong inversion: " << endl;
	PRINT3(qid, invseqs[qid].size(), ToStringBool( isWrongInver[qid] ) );
	wronginvs_seqs.push_back( invseqs[qid] );
	cout << "  assembly alignments:" << endl;
	for ( vec<int>::iterator ir = ia1_index[qid].begin(); ir != ia1_index[qid].end(); ir++ )
	  inversions_aligns1[*ir].PrintReadableBrief( cout );
	cout << "\n  reference alignments:" << endl;
	for ( vec<int>::iterator ir = ia2_index[qid].begin(); ir != ia2_index[qid].end(); ir++ )
	  inversions_aligns2[*ir].PrintReadableBrief( cout );
	cout << "                  ---------\n";
      }
    }
    
    wronginvs_seqs.WriteAll( wronginvs_file );

  }
  
  
  cout << Date( ) << ": done" << endl;

}

