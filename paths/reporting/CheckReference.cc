///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable

#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "util/SearchFastb2Core.h"
#include "lookup/LookAlign.h"
#include "Map.h"

/**
 * CheckReference
 *
 * Check for possible errors in reference genome with respect to sequencing data. 
 * This is done by editing a reference according to alignments of sequence sections 
 * (unibases for example) that contain differences between assembly and a reference genome.
 * Then we align jumping reads to a reference and edited reference and count support
 * specific unique support for both cases.
 *
*/

// this is used for sorting alternative alignments of the query according to their quality (as defined here)
Bool cmp_best_ali( const look_align& la1, const look_align& la2 ){    
  if ( la1.GaplessLength() > la2.GaplessLength() ) return True;
  if ( la1.GaplessLength() < la2.GaplessLength() ) return False;
  if ( la1.Errors() < la2.Errors() ) return True;
  if ( la1.Errors() > la2.Errors() ) return False;
  if ( la1.indels < la2.indels ) return True;
  if ( la1.indels > la2.indels ) return False;
  return False;    
}

// this is used for sorting alignments on the reference
Bool cmp_pos( const look_align& la1, const look_align& la2 ){    
  if ( la1.TargetId() < la2.TargetId() ) return True;
  if ( la1.TargetId() > la2.TargetId() ) return False;
  if ( la1.StartOnTarget() < la2.StartOnTarget() ) return True;
  if ( la1.StartOnTarget() > la2.StartOnTarget() ) return False;
  if ( la1.EndOnTarget() < la2.EndOnTarget() ) return True;
  if ( la1.EndOnTarget() > la2.EndOnTarget() ) return False;
  return False;    
}

int main( int argc, char *argv[] ){

  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  
  
  CommandArgument_String_OrDefault_Doc(REF,"genome", 
    "looks for DATA_DIR/REF.fastb");
  CommandArgument_String_Doc(SEQUENCES, 
    "looks for RUN_DIR/SEQUENCES.fastb with sequences containing assembly/reference differences");
  CommandArgument_String_OrDefault_Doc(SEQS_TO_PROCESS, "all",
    "ids of sequences used for alignment (default: take all)");
  CommandArgument_String_OrDefault_Doc(REF_EDIT,"genome.edit", 
    "looks for DATA_DIR/REF_EDIT.fastb");
  CommandArgument_String_OrDefault_Doc(JUMPS,"jump_reads_filt", 
    "looks for RUN_DIR/JUMPS.fastb");
  CommandArgument_UnsignedInt_OrDefault_Doc(TRIM_LEN, 60, 
    "length to trim out");
  CommandArgument_UnsignedInt_OrDefault_Doc(SEARCH_K, 40, 
    "K-mer size used in alignemnts");

  EndCommandArguments;


  String data_dir = PRE + "/" + DATA;
  String run_dir  = data_dir + "/" + RUN;
  
  String ref_file =  data_dir + "/" + REF + ".fastb";
  if ( ! IsRegularFile(ref_file) )
    SystemSucceed( "Fastb PRE=" + data_dir + " FILE=" + REF + ".fasta NH=True" );


  vecbvec ref( ref_file );
  

  // Align sequences to reference


  String sequences_file = run_dir + "/" + SEQUENCES;
  if ( ! IsRegularFile( sequences_file ) )
    sequences_file += ".fastb";
  vecbvec sequences( sequences_file );

  String seqs_to_process = SEQS_TO_PROCESS;
  PRINT( seqs_to_process );

  String ref_lookup_file = ref_file.SafeBefore(".fastb") + ".lookup";
  if ( ! IsRegularFile(ref_lookup_file) ){
    cout << Date() << ": Creating a lookup table for reference file" << endl;
    SystemSucceed( "MakeLookupTable" + ARG(SOURCE, ref_file )
		   + ARG(OUT_HEAD, ref_lookup_file.SafeBefore(".lookup")) 
		   + ARG(LOOKUP_ONLY, True)
		   + ARG( NH, True ) );
  }
  String tmp_dir = run_dir + "/checkRefTmp";
  Mkpath( tmp_dir );

  String aligns_file = tmp_dir + "/" + SEQUENCES + "." + ( seqs_to_process != "all" ? seqs_to_process.Between("{","}") : "all" ) + ".qlt";
  SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 " 
		 + ARG(SEQS,sequences_file) + ARG(SEQS_IS_FASTB, True) 
		 + ARG(L, ref_lookup_file) + ARG(PARSEABLE, True) 
		 + ARG(SMITH_WAT, True) + ARG(TMP_DIR,tmp_dir)
		 + ARG(SEQS_TO_PROCESS,"\"" + seqs_to_process + "\"")
		 + ARG(VISUAL,True)
		 + ARG(OUTPUT,aligns_file) );    

  vec<look_align> aligns;
  LoadLookAligns( aligns_file, aligns );
  for ( int i = 0; i < aligns.isize( ); i++ ){
    aligns[i].PrintReadableBrief(cout);
    //aligns[i].PrintVisual( cout, sequences[ aligns[i].QueryId()], ref[ aligns[i].TargetId()] );
  }
  
  map<size_t,vec<size_t> > aligns_index;
  for ( size_t i = 0; i < aligns.size(); i++ ){
    size_t qid = aligns[i].QueryId();
    aligns_index[qid].push_back(i);
  }
  
  vec<look_align> valid_aligns;
  for ( map<size_t,vec<size_t> >::iterator it = aligns_index.begin(); it != aligns_index.end(); it++ ){
    size_t qid = it->first;
    vec<size_t>& indices = it->second;
    vec<look_align> alts;
    for ( size_t ia = 0; ia < indices.size(); ia++ )
      alts.push_back( aligns[ indices[ia] ] );
    sort( alts.begin(), alts.end(), cmp_best_ali );
    vecbvec tseqs(1);
    tseqs.back().SetToSubOf( ref.at(alts[0].TargetId()), alts[0].StartOnTarget(), 
				 alts[0].EndOnTarget() - alts[0].StartOnTarget()       );
    valid_aligns.push_back( alts[0] );
    for ( size_t ia = 1; ia < alts.size(); ia++ ){
      basevector tseq;
      tseq.SetToSubOf( ref.at(alts[ia].TargetId()), alts[ia].StartOnTarget(), 
		       alts[ia].EndOnTarget() - alts[ia].StartOnTarget()       );
      if ( tseq == tseqs[0] )
	valid_aligns.push_back( alts[ia] );
      else break;
    }
  }

  sort( valid_aligns.begin(), valid_aligns.end(), cmp_pos );

  


  // xxdp todo: make sure that all alignments are found
  vec< StdMap<int,char> > mutations( ref.size() );
  vec< StdMap<int,basevector> > insertions( ref.size() );
  vec< StdMap<int,bool> > deletions( ref.size() );
  int allmuts = 0, allindels=0;

  for ( size_t ia = 0; ia < valid_aligns.size(); ia++){
    const look_align & la = valid_aligns[ia];
    int qid     = la.QueryId();
    int tid     = la.TargetId();
    int nmuts   = la.mutations;
    int nindels = la.indels;
    basevector query  = sequences[qid];
    if ( la.Rc1() ) query.ReverseComplement();
    const basevector& target = ref[tid];
    align ar = la.a;
    int qstart = ar.pos1();
    int tstart = ar.pos2();
    int qpos = qstart;
    int tpos = tstart;
    int checkindels = 0, checkmuts = 0;
    //PRINT(ar.Nblocks());
    for ( int bi = 0; bi < ar.Nblocks(); bi++ ){
      int length = ar.Lengths(bi);
      int gap = bi < ar.Nblocks() -1 ? ar.Gaps(bi+1) : 0;
      for ( int il = 0; il < length; il++ ){
	//PRINT2( qpos, tpos );
	if ( query[qpos] != target[tpos] ){
	  checkmuts++;
	  //cout << "mutation: " << as_base(query[qpos]) << " vs " << as_base(target[tpos]) << " "; PRINT5( qpos, tpos, bi, length, il );
	}
	if ( deletions[tid].find(tpos) != deletions[tid].end() ){
	  // inconsistent alignments. Deletion was suggested in this position before.
	  FatalErr("inconsistent alignment with deletion");
	}
	if ( mutations[tid].find(tpos) != mutations[tid].end() && query[qpos] != mutations[tid][tpos] ){
	  FatalErr("inconsistent alignments");
	}
	mutations[tid][tpos] = query[qpos];
	qpos++;
	tpos++;
      }
      if ( gap < 0 ){
	//PRINT3( gap, qpos, tpos );
	basevector insertion( query,qpos+1,-gap);
	if ( insertions[tid].find(qpos) != insertions[tid].end() && insertions[tid][qpos] != insertion ){
	  FatalErr("inconsistent insertions");
	}
	insertions[tid][tpos] = insertion;
	qpos += -gap;
	checkindels += -gap;
	//PRINT3( "endgap", qpos, tpos );
      }else if ( gap > 0 ){
	//PRINT3( gap, qpos, tpos );
	for ( int gi = 0; gi < gap; gi++ ){
	  deletions[tid][tpos + gi] = true;
	}
	tpos += gap;
	checkindels += gap;
	//PRINT3( "endgap", qpos, tpos );
      }
    }
    ForceAssertEq( nmuts, checkmuts );
    ForceAssertEq( nindels, checkindels );
  }
  
  // build edited reference
  vecbvec ref_edit( ref.size() );
  for ( size_t tid = 0; tid < ref.size(); tid++ ){
    basevector contig(ref[tid]);
    int indelSum = 0;
    if ( mutations[tid].size() > 0 )
      for ( map<int,char>::const_iterator it = mutations[tid].begin(); it != mutations[tid].end(); it++ )
	contig.Set( it->first, it->second );
    
    if ( insertions[tid].size() > 0 || deletions[tid].size() > 0 ){
      size_t cpos = 0;
      basevector new_contig;
      for ( size_t cpos = 0; cpos < contig.size(); cpos++ ){
	if ( insertions[tid].find(cpos) == insertions[tid].end() && deletions[tid].find(cpos) == deletions[tid].end() )
	  new_contig.push_back( contig[cpos] );
	else if ( deletions[tid].find(cpos) != deletions[tid].end() ){
	  indelSum--;
	}
	else if ( insertions[tid].find(cpos) != insertions[tid].end() ){
	  new_contig.push_back( contig[cpos] );
	  new_contig = Cat( new_contig, insertions[tid][cpos] );
	  indelSum += insertions[tid][cpos].isize();
	}
      }
      PRINT( new_contig.size() );
      contig = new_contig;
    }
    ref_edit[tid] = contig;
    //PRINT5( ref.size(), tid, ref[tid].isize(), indelSum, ref_edit[tid].isize() );
    ForceAssertEq( ref[tid].isize() + indelSum, ref_edit[tid].isize() );
  }
  

  String ref_edit_file =  tmp_dir + "/" + REF_EDIT + ".fastb";
  ref_edit.WriteAll( ref_edit_file );

  
  // align pairs 


  String pairs_file = run_dir + "/" + JUMPS + ".pairs";
  PairsManager pairs( pairs_file);
  PRINT( pairs.nPairs() );
  
  String jumps_file = run_dir + "/" + JUMPS + ".fastb";  
  vecbasevector reads( jumps_file );
  int64_t nreads = reads.size( );
  int len = TRIM_LEN;
  for ( size_t id = 0; id < reads.size( ); id++ ){    
    ForceAssertGe( reads[id].isize( ), len );
    reads[id].SetToSubOf( reads[id], reads[id].isize( ) - len, len );    
  }
  
  String jumps_trim_file = tmp_dir + "/" + JUMPS + "_trimmed.fastb"; 
  reads.WriteAll( jumps_trim_file );

  String F1     = jumps_trim_file;
  String F2_old = ref_file;
  String F2_new = ref_edit_file;
  
  int K = SEARCH_K;
  vec< triple<int64_t,int64_t,int> > aligns_old, aligns_new;
  SearchFastb2( F1, F2_old, K, &aligns_old, 0, -1, 0.90, False );
  SearchFastb2( F1, F2_new, K, &aligns_new, 0, -1, 0.90, False );
  
  vec< vec<int> > aligns_old_index(nreads), aligns_new_index(nreads);
  for ( int i = 0; i < aligns_old.isize( ); i++ )
    aligns_old_index[ aligns_old[i].first ].push_back(i);
  PRINT( aligns_old.size() );
  for ( int i = 0; i < aligns_new.isize( ); i++ )
    aligns_new_index[ aligns_new[i].first ].push_back(i);
  PRINT( aligns_new.size() );

  vec<int64_t> good_old, good_new;
  for ( int pass = 1; pass <= 2; pass++ ){    
    const vec< triple<int64_t,int64_t,int> >& aligns
      = ( pass == 1 ? aligns_old : aligns_new );
    const vec< vec<int> >& aligns_index 
      = ( pass == 1 ? aligns_old_index : aligns_new_index );
    for ( size_t pid = 0; pid < pairs.nPairs( ); pid++ ){    
      int64_t id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
      Bool good = False;
      for ( int j1 = 0; j1 < aligns_index[id1].isize( ); j1++ )
	for ( int j2 = 0; j2 < aligns_index[id2].isize( ); j2++ ){    
	  int p1 = aligns_index[id1][j1], p2 = aligns_index[id2][j2];
	  if ( aligns[p1].second != aligns[p2].second ) continue;
	  if ( aligns[p1].third < 0 ) swap( p1, p2 );
	  if ( aligns[p1].third < 0 ) continue;
	  if ( aligns[p2].third >= 0 ) continue;
	  int pos1 = aligns[p1].third, pos2 = -aligns[p2].third-1;
	  int Pos1 = pos1 + reads[ aligns[p1].first ].isize( );
	  int sep = pos2 - Pos1;
	  if ( sep < 0 || sep > 10000 ) continue;
	  good = True;    
	}
      if ( !good ) continue;
      if ( pass == 1 ) good_old.push_back(pid);
      else good_new.push_back(pid);    
    }    
  }
  PRINT2( good_old.size(), good_new.size() );
  
  int delta1 = 0, delta2 = 0;
  for ( int j = 0; j < good_old.isize( ); j++ )
    if ( !BinMember( good_new, good_old[j] ) ) delta1++;
  for ( int j = 0; j < good_new.isize( ); j++ )
    if ( !BinMember( good_old, good_new[j] ) ) delta2++;
  
  cout << delta1 << " pairs favor " + REF + "\n";
  cout << delta2 << " pairs favor " + REF_EDIT + "\n";    
}
