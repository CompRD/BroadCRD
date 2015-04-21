///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeScaffolds.  Build scaffolds from assembly graph.
// Run this after AlignPairsToHyperLG.

#include "Equiv.h"
#include "Fastavector.h"
#include "Charvector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "feudal/FeudalTools.h"
#include "graph/Digraph.h"
#include "lookup/LookAlign.h"
#include "paths/Alignlet.h"
#include "paths/HyperFastavector.h"
#include "paths/PseudoUnipatherLib.h"
#include "paths/SInsertion.h"
#include "paths/SuckScaffolds.h"


// class that allows merging and breaking contigs and updating fastavectors and alignments at the same time
class ContigsManager {
  
private:
  vec<fastavector> * _contigs_fv; //fastavectors of contigs
  vec<alignlet> * _aligns;  //alignment information
  vec<int> * _aligns_index; //alignments index
  vec< vec<int> > _contigs2reads; //map of reads aligned to contigs;

public:
  // constructor
  ContigsManager( vec<fastavector> & contigs_fv, vec<alignlet> & aligns,
		  vec<int> & aligns_index ){
   

    _contigs_fv   = & contigs_fv;
    _aligns       = & aligns;
    _aligns_index = & aligns_index;
    // build the map of reads aligned to contigs
    _contigs2reads.resize( _contigs_fv->size() );
    for ( size_t rid = 0; rid < _aligns_index->size(); rid++ ){
      if ( (*_aligns_index)[rid] < 0 ) continue;
      alignlet& la = (*_aligns)[ (*_aligns_index)[rid] ];
      size_t t = la.TargetId();
      ForceAssertLt( t, _contigs2reads.size() );
      _contigs2reads[t].push_back( rid );
    }
  }
  // split contig
  void SplitContig( const size_t & t, const size_t & s1, const size_t & s2 ){
    fastavector tig_fv = (*_contigs_fv)[t];
    ForceAssert( s1 <= s2 ); ForceAssert( s2 < tig_fv.size() );
    
    fastavector tg0; tg0.SetToSubOf( tig_fv, 0, s1 -1 );
    (*_contigs_fv).push_back( tg0 ); 
    size_t tg0_id = (*_contigs_fv).size() -1;
    
    fastavector tg1; 
    size_t tg1_id;
    if ( s2 > s1 ){
      tg1.SetToSubOf( tig_fv, s1, s2 -1 );
      (*_contigs_fv).push_back( tg1 ); 
      tg1_id = (*_contigs_fv).size() -1;
    }
    
    fastavector tg2; tg2.SetToSubOf( tig_fv, s2, tig_fv.size() - s2 + 2 );
    (*_contigs_fv).push_back( tg2 ); 
    size_t tg2_id = (*_contigs_fv).size() -1;
  
    // update alignments
    _contigs2reads.resize( _contigs_fv->size() );
    for ( vec<int>::iterator it = _contigs2reads[t].begin(); it != _contigs2reads[t].end(); it++ ){
      size_t rid   = *it;
      alignlet& la = (*_aligns)[rid];
      if ( la.Pos2() < (int) s1 ){
	alignlet aloc( la.pos2(), la.Pos2(), tg0_id, 
		       (*_contigs_fv).back().size(), la.TargetLength() > 0 ? True : False );
	(*_aligns)[rid] = aloc;
	_contigs2reads[tg0_id].push_back( rid );
      }else if ( s2 > s1 && la.pos2() >= (int)s1 && la.Pos2() < (int)s2 ){
	alignlet aloc( la.pos2(), la.Pos2(), tg1_id, 
		       (*_contigs_fv).back().size(), la.TargetLength() > 0 ? True : False );
	(*_aligns)[rid] = aloc;
	_contigs2reads[tg1_id].push_back( rid );
      }
      else if ( la.pos2() >= (int) s2 && la.Pos2() < (int) tig_fv.size() ){
	alignlet aloc( la.pos2(), la.Pos2(), tg2_id, 
		       (*_contigs_fv).back().size(), la.TargetLength() > 0 ? True : False );
	(*_aligns)[rid] = aloc;
	_contigs2reads[tg2_id].push_back( rid );
      }
      else{
	(*_aligns_index)[rid] = -1;
      }
    }
    
    (*_contigs_fv)[t].clear();
    _contigs2reads[t].clear();
  }
  
  void SplitContig( const size_t & ts, const size_t & s1 ){
    SplitContig( ts, s1, s1 );
  }
  
  void MergeContigs( const size_t & t1, const size_t & t2, const align & a ){
    
    fastavector f1 = (*_contigs_fv)[t1];
    fastavector f2 = (*_contigs_fv)[t2];
    
    ForceAssertLt( a.pos2(), (int) f2.size() ); 

    int nblocks = a.Nblocks( );

    vecbasevector vb1 = f1.AllBasevectors();
    vecbasevector vb2 = f2.AllBasevectors();
    
    basevector b1 = vb1[0];
    basevector b2 = vb2[0];

    int errors = ActualErrors( b1, b2, a );

    int pos1 = a.pos1( ), pos2 = a.pos2( ), Pos1 = a.Pos1( ), Pos2 = a.Pos2( );
    const avector<int> &gaps = a.Gaps( ), &lengths = a.Lengths( );
    if ( pos1 > 0 && pos2 > 0 ){    
      cout << "MergeContigs: I've got pos1 = " << pos1 <<
	" and pos2 = " << pos2 << ", which is not what I expected.\n";
      ForceAssert( 0 == 1 );    
    }
    if ( Pos1 < (int)f1.size() -1 && Pos2 < (int)f2.size() -1 ){
      cout << "MergeContigs: I've got:"; 
      PRINT4( Pos1, f1.size(), Pos2, f2.size() );
      cout << " which is not what I expected.\n";
      ForceAssert( 0 == 1 );
    }
    
    fastavector c;
    c.resize( b1.size( ) + b2.size( ) + errors );
    vec<int> mapf1Toc( f1.size(), -1 ), mapf2Toc( f2.size(), -1 );
    
    int i = 0, p1 = pos1, p2 = pos2;
    if ( pos2 == 0 )
      for ( ; i < pos1; i++ ){
	c[i] = f1[i];
	mapf1Toc[i] = i;
      }
    else 
      for ( ; i < pos2; i++ ){
	c[i] = f2[i];
	mapf2Toc[i] = i;
      }

    for ( int j = 0; j < nblocks; j++ ){    
      if ( gaps(j) > 0 )
	for ( int x = 0; x < gaps(j); x++ ){    
	  c[i] = f2[p2];
	  mapf2Toc[p2] = i;
	  ++p2;
	  ++i;    
	}
      if ( gaps(j) < 0 )
	for ( int x = 0; x < -gaps(j); x++ ){    
	  c[i] = f1[p1];
	  mapf1Toc[p1] = i;
	  ++p1;
	  ++i;    
	}
      for ( int x = 0; x < lengths(j); x++ ){    
	if ( f1[p1] == f2[p2] ) c[i] = f1[p1];   
	else c[i]= 'N';   
	mapf1Toc[p1] = mapf2Toc[p2] = i;
	++p1;
	++p2;
	++i; 
      }   
    }
    
    if ( p1 < (int) f1.size() && p2 < (int) f2.size() ){
      cout << "MergeContigs: I've got:";
      PRINT4( p1, f1.size(), p2, f2.size() );
      cout << ", which is not what I expected.\n";
      ForceAssert( 0 == 1 );
    }
    
    if ( p1 < (int) f1.size() ){
      for ( int j = 0; j < (int) f1.size() - Pos1; j++ ){
	c[i] = f1[p1];
	mapf1Toc[p1] = i;
	p1++;
	i++;
      }
    }else{
      for ( int j = 0; j < (int) f2.size() - Pos2; j++ ){
	c[i] = f2[p2];
	mapf2Toc[p2] = i;
	p2++;
	i++;
      }
    }
    c.resize(i);
    
    (*_contigs_fv)[t1].clear();
    (*_contigs_fv)[t2].clear();
    (*_contigs_fv).push_back(c);
    int new_tg_id = (*_contigs_fv).size() -1;
    _contigs2reads.resize( _contigs_fv->size() );
    
    for ( vec<int>::iterator it = _contigs2reads[t1].begin(); 
	  it != _contigs2reads[t1].end(); it++ ){
      size_t rid   = *it;
      alignlet& la = (*_aligns)[rid];
      int pos2 = la.pos2();
      alignlet aloc( mapf1Toc[pos2], la.Pos2() + (mapf1Toc[pos2] - pos2 ), new_tg_id, 
		     (*_contigs_fv).back().size(), la.TargetLength() > 0 ? True : False );
      (*_aligns)[rid] = aloc;
      _contigs2reads[new_tg_id].push_back( rid );
    }
    _contigs2reads[t1].clear();
    for ( vec<int>::iterator it = _contigs2reads[t2].begin(); it != _contigs2reads[t2].end(); it++ ){
      size_t rid   = *it;
      alignlet& la = (*_aligns)[rid];
      int pos2 = la.pos2();
      alignlet aloc( mapf2Toc[pos2], la.Pos2() + (mapf2Toc[pos2] - pos2 ), new_tg_id, 
		     (*_contigs_fv).back().size(), la.TargetLength() > 0 ? True : False );
      (*_aligns)[rid] = aloc;
      _contigs2reads[new_tg_id].push_back( rid );
    }
    _contigs2reads[t2].clear();
  }

  // overloading
  void MergeContigs( const size_t & t1, const size_t & t2, const int & offset ){
    int pos1, pos2;
    int length;
    if ( offset >= 0 ){
      ForceAssertLt( offset, (int) (*_contigs_fv)[t1].size() );
      pos1 = offset; pos2 = 0;
      length = (*_contigs_fv)[t1].size() - offset;
    }else{
      ForceAssertLt( -offset, (int) (*_contigs_fv)[t2].size() );
      pos1 = 0; pos2 = (*_contigs_fv)[t2].size() + offset;
      length = -offset;
    }
    align pseudo_a;
    pseudo_a.Setpos1( pos1 );
    pseudo_a.Setpos2( pos2 );
    pseudo_a.SetGap( 0, 0 );
    pseudo_a.SetLength( 0, length ); 
    MergeContigs( t1, t2, pseudo_a );
  }

};




// An slink is a link between scaffolds.  The first scaffold is not represented
// in the object.

class slink {

     public:

     slink( ) { }
     slink( const int s2, const Bool fw2, const int sep, const int dev )
          : s2(s2), fw2(fw2), sep(sep), dev(dev) { }

     int s2;
     Bool fw2;
     int sep, dev;

     friend Bool operator<( const slink& l1, const slink& l2 ){    
       if ( l1.s2 < l2.s2 ) return True;
       if ( l1.s2 > l2.s2 ) return False;
       return l1.fw2 < l2.fw2;    
     }

};

int main( int argc, char *argv[] ){
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  
  // Input (from <PRE>/<DATA>/<RUN>):
  CommandArgument_String_OrDefault( READS, "all_reads" );
  
  // Input (from <PRE>/<DATA>/<RUN>/<SUBDIR>):
  CommandArgument_String_OrDefault( ALIGNS, "" );
  CommandArgument_String_OrDefault( HYPER, "hyper_plus_edited" );
  
  // Either input or output. INITIAL_SCAFFOLDS also lives in
  //  <PRE>/<DATA>/<RUN>/<SUBDIR>, but with a twist. If
  //  READ_INITIAL_SCAFFOLDS is True, then INITIAL_SCAFFOLDS is
  //  loaded from file.  Otherwise, the initial scaffolds are
  //  generated with ConvertToSuperb, and saved to INITIAL
  //  SCAFFLODS (if WRITE = True, see below).
  CommandArgument_String_OrDefault( INITIAL_SCAFFOLDS,
				    "linear_scaffolds.initial" );
  
  CommandArgument_String_OrDefault( SCAFFOLDS_IN, "initial_scaffolds" );
  
  // Output saved in <PRE>/<DATA>/<RUN>/<SUBDIR>/<ASSEMBLIES>:
  //   <SCAFFOLDS_OUT>.{contigs.fasta,assembly.fasta,superb}
  
  CommandArgument_String_OrDefault( SCAFFOLDS_OUT, "scaffolds" );
  
  // Output (in <PRE>/<DATA>/<RUN>/<SUBDIR>):
  CommandArgument_String_OrDefault( SCAFFOLDS, "linear_scaffolds" );
  
  // The following three arguments must be in sync. For each
  // iteration, accept links between two scaffolds if and only if
  // there are at least MIN_LINKS links of size in the interval
  // specified by MIN_PAIR_SEPS and MAX_PAIR_SEPS.
  CommandArgument_String_OrDefault(MIN_PAIR_SEPS, "");
  CommandArgument_String_OrDefault(MAX_PAIR_SEPS, ""); 
  CommandArgument_String_OrDefault( MIN_LINKS, "{20,6,4,2}" );
  
  // If true, run SuckScaffold to absorb small scaffolds into the
  // gaps of larger ones.
  CommandArgument_Bool_OrDefault( SUCK_SCAFFOLDS, False );
  
  // If true, save both initial super structure and result.
  CommandArgument_Bool_OrDefault( WRITE, True );
  
  // Various other arguments to tune the algorithm:
  CommandArgument_Int_OrDefault_Doc(MIN_EDGE_TO_SAVE, 1000,
      "if scaffold can't be linearized, convert edges of this size or "
				    "larger into separate scaffolds");
  CommandArgument_Int_OrDefault_Doc(MAX_EDGE_TO_DISCARD, 10000,
      "don't accept a linearization of a scaffold if it omits an edge "
				    "of this size or larger");
  CommandArgument_Int_OrDefault(VERBOSITY, 0);
  CommandArgument_Int_OrDefault( MIN_LINKS_TO_PRINT, 1000000000 );
  CommandArgument_Int_OrDefault(MAX_THREADS, 16);
  CommandArgument_Bool_OrDefault(JOIN, True);
  CommandArgument_Double_OrDefault_Doc(CONNECTIVITY_LIMIT, 1000000,
      "scaffolds with more than this many links/kb will not be linked");

  EndCommandArguments;
  
  
  if ( ALIGNS == "" )
    ALIGNS = READS;
  
  
  // Define heuristic constants. These could be made into arguments.
  const int dev_mult = 3;
  const int pair_seps_bin_size = 100000; // large value => effectively no binning 
  const int mean_jump_construct_size = 300; // mean size of sheared fragments
  // Set up directories.
 

   // Dir and file names.

  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  
  String reads_file = run_dir + "/" + READS + ".fastb";
  String pairs_file = run_dir + "/" + READS + ".pairs";
  String contigs_in_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
  String scaffolds_in_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
  String aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
  String index_file = sub_dir + "/" + ALIGNS + ".qltoutlet.index";
  
  String out_suffix = SCAFFOLDS_OUT + ( JOIN ? "" : ".unjoined" );
  String contigs_out_file = sub_dir + "/" + out_suffix + ".contigs.fasta";
  String assembly_out_file = sub_dir + "/" + out_suffix + ".assembly.fasta";
  String scaffolds_out_file = sub_dir + "/" + out_suffix + ".superb";
  
  // Load scaffold files (vec<fastavector> and superb.)
  cout << Date( ) << ": loading scaffolds" << endl;
  vec<superb> scaffolds;
  ReadSuperbs( scaffolds_in_file, scaffolds );
  vec<fastavector> scaffold_contigs;
  LoadFromFastaFile( contigs_in_file, scaffold_contigs );
  size_t ntigs = scaffold_contigs.size();
  
  // Check that the files are consistent.
  size_t n_tigs_in_scaffolds = 0;
  for ( size_t i = 0; i < scaffolds.size(); i++ ) {
    n_tigs_in_scaffolds += scaffolds[i].Ntigs();
    for ( int j = 0; j < scaffolds[i].Ntigs(); j++ ) {
      int tig = scaffolds[i].Tig(j);
      ForceAssertEq( (int)scaffold_contigs[tig].size(), scaffolds[i].Len(j) );
    }
  }
  ForceAssertGe( ntigs, n_tigs_in_scaffolds );
  
  // Load alignments and convert to lightweight form.
  // aligns0_index works as this: say aligns0_index[id] = ii, then:
  //  ii = -3 :   look_aligns not informative for this read 
  //              (e.g. both reads align to the same scaffold)
  //  ii = -2 :   no look_aligns for this read
  //  ii = -1 :   more than one look_align for this read (id)
  //  ii >= 0 :   aligns0[ii] is the only look_align for this read
  
  cout << Date( ) << ": loading aligns" << flush;
  vec<alignlet> aligns0;
  BinaryReader::readFile( aligns_file, &aligns0 );

  vec<int> aligns0_index;
  BinaryReader::readFile( index_file, &aligns0_index );
  longlong nreads = MastervecFileObjectCount( reads_file );
  cout << " (" << aligns0.size( ) << " aligns found)" << endl;
  
  ContigsManager contigs_info( scaffold_contigs, aligns0, aligns0_index );
  
  // Load read pairs.
  
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( run_dir + "/" + READS + ".pairs" );
  
  // Read pair separation limits:  With long insert read libraries,
  // we may wish to consider them in sequence to find nearby joins
  // before matches which might jump over other scaffolds. Default
  // uses everthing in one iteration, as before. --bruce 30 Sep 09
  
  vec<int> min_links_set;
  ParseIntSet(MIN_LINKS, min_links_set, False);

  vec<int> min_pair_seps;
  vec<int> max_pair_seps;
  if ( ! MIN_PAIR_SEPS.empty() ){   
    // use separation limits from arguments
    ForceAssert( ! MAX_PAIR_SEPS.empty() );
    ParseIntSet(MIN_PAIR_SEPS, min_pair_seps, False);
    ParseIntSet(MAX_PAIR_SEPS, max_pair_seps, False);
  }else{
    // derive separation limits from library information
    vec<int> libID2Sep( pairs.nLibraries() );
    for (unsigned libID = 0; libID < pairs.nLibraries(); libID++ ){
      libID2Sep[libID] = pairs.getLibrarySep( libID );
    }
    vec<int> sortedLibSeps = libID2Sep; 
    UniqueSort( sortedLibSeps );
    int minimum_sep = sortedLibSeps[0] < 0 ? sortedLibSeps[0] : 0;
    if ( VERBOSITY >= 2 ){
      cout << "sortedLibSeps: "; sortedLibSeps.Print( cout ); cout << endl;
    }
    min_pair_seps.push_back( minimum_sep );
    max_pair_seps.push_back( sortedLibSeps[0] );
    int bin_beg = sortedLibSeps[0];
    for ( unsigned i = 1; i < sortedLibSeps.size(); i++ ){      
      if ( sortedLibSeps[i] > bin_beg + pair_seps_bin_size ){
	min_pair_seps.push_back( minimum_sep );
	max_pair_seps.push_back( sortedLibSeps[i] );
	bin_beg = sortedLibSeps[i];
      }else{
	min_pair_seps.back() = minimum_sep;
	max_pair_seps.back() = sortedLibSeps[i];
      }
    }
    ForceAssertEq( max_pair_seps.back(), sortedLibSeps.back() );
  }
  if ( VERBOSITY >= 1 ){
    cout << "min_pair_seps: "; min_pair_seps.Print(cout); cout << endl;
    cout << "max_pair_seps: "; max_pair_seps.Print(cout); cout << endl;
  }
  ForceAssertEq( min_pair_seps.size(), max_pair_seps.size() );



  // Track orientations of contigs.
  
  vec<Bool> rctig( ntigs, False );
  
  // Now begin the scaffold joining operation.  This is repeated until there
  // no more joins are made.
  
  cout << Date( ) << ": Joining scaffolds... " << endl;
  
  if ( VERBOSITY >= 1 ){    
    cout << "\ninitial scaffolds:\n\n";
    for ( size_t i = 0; i < scaffolds.size(); i++ )
      scaffolds[i].Print( cout, "scaffold_" + ToString(i) );    
  }
  int iterationCt = 0, allJoinsCt = 0; // help track global progress
    
  for ( unsigned mls = 0; mls < min_links_set.size(); mls++ ){ 
    size_t min_links_size = min_links_set.size();
    int min_links         = min_links_set[ mls ];
   
    for (size_t mps = 0; mps < max_pair_seps.size(); ++mps){
      int min_pair_sep = min_pair_seps[mps];
      int max_pair_sep = max_pair_seps[mps];
      int joins = 0; int localJoinsCt =  0; // help track local progress
      
      if ( VERBOSITY >= 1 )
	cout << Date() << ": Using min_links=" << min_links
	     << " separations=[" << min_pair_sep
	     << "," << max_pair_sep
	     << "]" << endl;
      do {
	
	int nscaff = scaffolds.size( );
	iterationCt++;
	if ( VERBOSITY >= 1 ) 
	  PRINT( iterationCt );
	joins = 0;
	
	// Build links between scaffolds.  There is one pass for each orientation.
	// If s1 is a scaffold id, then fw_fw[s1] consists of triples
	// (s2, sep, dev) where s2 is another scaffold id such that there are
	// "enough" read pairs having one read forward on s1 and the other forward
	// on s2.  The construction is symmetric: fw_fw[s2] includes
	// (s1, sep, dev).
	//
	// Note that the link is actually from s1 to s2rc (or from s2 to s1rc).
	//
	// The vectors fw_rc, rc_fw, and rc_rc are similarly defined.
	
	longlong outConflictCt = 0, inConflictCt = 0;
	longlong propAlgndPairsCt = 0, pairsCt = 0;
	
	vec< vec< triple<int,int,int> > >
	  fw_fw(nscaff), fw_rc(nscaff), rc_fw(nscaff), rc_rc(nscaff);
	for ( int pass = 1; pass <= 2; pass++ ){
	  // Reverse stuff.
	  for ( int i = 0; i < nscaff; i++ )
	    scaffolds[i].Reverse( );
	  for ( size_t i = 0; i < aligns0.size(); i++ )
	    aligns0[i].Reverse( );
	  // Build indices that allow one find the scaffold and position on the
	  // scaffold for a given contig.
	  
	  vec<int> tig_to_scaffold( ntigs, -1 );
	  vec<int> tig_to_scaffold_pos( ntigs, -1 );
	  for ( int i = 0; i < nscaff; i++ ){    
	    for ( int j = 0; j < scaffolds[i].Ntigs( ); j++ ){    
	      tig_to_scaffold[ scaffolds[i].Tig(j) ] = i;
	      tig_to_scaffold_pos[ scaffolds[i].Tig(j) ] = j;    
	    }    
	  }
	  
	  // Find forward links from scaffolds.
	  
	  vec< vec<slink> > slinks(nscaff);
	  for ( ulonglong pair_ID = 0; pair_ID < pairs.nPairs(); pair_ID++ ){
	    longlong id1 = pairs.ID1( pair_ID );
	    longlong id2 = pairs.ID2( pair_ID );
	    // We only consider read pairs for which both ends are uniquely
	    // placed.
	    if ( aligns0_index[id1] < 0 || aligns0_index[id2] < 0 )
	      continue;
	    if ( pairs.sep(pair_ID) < min_pair_sep ||
		 pairs.sep(pair_ID) > max_pair_sep)
	      continue;
	    pairsCt++;
	    
	    // Figure out what the link is.  Ignore intrascaffold links, as
	    // well as those for which id1 is not forward on the contig.
	    
	    const alignlet& la1 = aligns0[ aligns0_index[id1] ];
	    if ( ! la1.Fw1() ) continue;
	    const alignlet& la2 = aligns0[ aligns0_index[id2] ];
	    int t1 = la1.TargetId( ), t2 = la2.TargetId( );
	    int s1 = tig_to_scaffold[t1], s2 = tig_to_scaffold[t2];
	    if ( s1 < 0 || s2 < 0 ) {
	      continue;
	      
	    }
	    if ( s1 == s2 ) {
	      // speedup, set negative index to skip these read pairs in next iterations
	      aligns0_index[id1] = -3;  
	      aligns0_index[id2] = -3;
	      continue;
	    }
	    
	    
	    
	    
	    int p1 = tig_to_scaffold_pos[t1], p2 = tig_to_scaffold_pos[t2];
	    const superb &S1 = scaffolds[s1], &S2 = scaffolds[s2]; 
	    
	    propAlgndPairsCt++;
	    // Now do some math to decide the gap and its deviation, and
	    // whether the gap is within bounds.  There are lots of
	    // heuristics here.
	    
	    int dist_to_end1 = scaffold_contigs[t1].size() - la1.Pos2( )
	      + S1.SubSuperLength( p1+1, S1.Ntigs( ) - 1 );
	    int dist_to_end1_dev
	      = S1.SubSuperLengthDev( p1+1, S1.Ntigs( ) - 1 );

	    int dist_to_oend1 = la1.Pos2( ) + S1.SubSuperLength( 0, p1 -1 );

	    int dist_to_end2, dist_to_end2_dev;
	    int dist_to_oend2;
	    if ( la2.Fw1( ) ){    
	      dist_to_end2 = scaffold_contigs[t2].size() - la2.Pos2( )
		+ S2.SubSuperLength( p2+1, S2.Ntigs( ) - 1 );
	      dist_to_end2_dev
		= S2.SubSuperLengthDev( p2+1, S2.Ntigs( ) - 1 );

	      dist_to_oend2 = la2.Pos2( ) + S2.SubSuperLength( 0, p2 -1 );
	    }
	    else{    
	      dist_to_end2 = la2.pos2( ) + S2.SubSuperLength( 0, p2-1 );
	      dist_to_end2_dev = S2.SubSuperLengthDev( 0, p2-1 );
	      
	      dist_to_oend2 = scaffold_contigs[t2].size() - la2.pos2( ) 
		+ S2.SubSuperLength( p2+1, S2.Ntigs( ) - 1 );
	    }
	    
	    // do not consider near other edge alignments, possibly avoid some 'innies' 
	    // from sheared libraries
	    if ( dist_to_oend1 + dist_to_oend2 < mean_jump_construct_size ) continue;
	    
	    // disregard alignments that are too far away from edges, they are likely
	    // to be just random noise. (They could be used to indicate problems with scaffolds)
	    if ( dist_to_end1 - dev_mult * dist_to_end1_dev > 
		 pairs.sep(pair_ID) + dev_mult * pairs.sd(pair_ID) ){
	      continue;
	    } 
	    if ( dist_to_end2 - dev_mult * dist_to_end2_dev  > 
		 pairs.sep(pair_ID) + dev_mult * pairs.sd(pair_ID) ){
	      continue;
	    } 
	    
	    int sep = pairs.sep( pair_ID ) - dist_to_end1 - dist_to_end2;
	    int dev = dist_to_end1_dev + dist_to_end2_dev + pairs.sd( pair_ID );
	    int min_gap = sep - dev_mult * dev;
	    int max_gap = sep + dev_mult * dev;
	    
	    slinks[s1].push( s2, la2.Fw1( ), sep, dev );    
	  }
	  
	  // Condense the links.  We look at all the links between two
	  // scaffolds and decide if there are enough.
	  
	  vec< vec<slink> > slinks2(nscaff), slinks2_to(nscaff);
	  for ( int s1 = 0; s1 < nscaff; s1++ ){    
	    Sort( slinks[s1] );
	    for ( size_t i = 0; i < slinks[s1].size(); i++ ){    
	      int s2 = slinks[s1][i].s2;
	      Bool fw2 = True;
	      size_t j;
	      vec<double> seps, devs, seps_rc, devs_rc;
	      for ( j = i; j < slinks[s1].size(); j++ ){ 
		if ( slinks[s1][j].s2 != s2 ) break;
		if (slinks[s1][j].fw2) {
		  seps.push_back( slinks[s1][j].sep );
		  devs.push_back( slinks[s1][j].dev );
		} else {
		  seps_rc.push_back( slinks[s1][j].sep );
		  devs_rc.push_back( slinks[s1][j].dev );
		}
	      }
	      if (seps.size() < seps_rc.size()) {
		fw2 = False;
		seps = seps_rc;
		devs = devs_rc;
	      }
	      
	      int count = (int)seps.size();
	      if ( count >= MIN_LINKS_TO_PRINT ) 
		PRINT3( s1, s2, count );
	      
	      // Decide if there are enough links.
	      
	      if ( count < min_links ){   
		i = j - 1;
		continue;    
	      }
	      
	      // Convert to a single link.  The mechanism for defining
	      // the separation and deviation is heuristic.
	      
	      SortSync( devs, seps );
	      double Sep, Dev;
	      CombineStats( seps, devs, Sep, Dev );
	      
	      // Save the link.

	      slinks2[s1].push(s2, fw2, int(round(Sep)), int(round(Dev)) );
	      // Re-using this class to record the pointers to each scaffold
	      slinks2_to[s2].push(s1, fw2, int(round(Sep)), int(round(Dev)) );
	      i = j - 1;    
	    }    
	  }
	  
	  // Identify small scaffolds with high connectivity.
	  // We'll define connectivity as scaffold links in and
	  // out per kb. Anything with a connectivity metric
	  // larger than CONNECTIVITY_LIMIT will not be
	  // considered for joining. --bruce 6 Oct 09
	  for (int s = 0; s < nscaff; s++) {
	    const superb &S = scaffolds[s];
	    double length = S.SubSuperLength( 0, S.Ntigs( ) - 1 );
	    double connectivity
	      = 1000.0 * (slinks2[s].size() + slinks2_to[s].size()) / length;
	    
	    if (connectivity > CONNECTIVITY_LIMIT) {
	      if (VERBOSITY>=3)
		cout << "connectivity limit exceeded: s=" << s
		     << " c=" << connectivity << " l=" << length << endl;
	      for (size_t i = 0; i < slinks2_to[s].size(); ++i) {
		int s1 = slinks2_to[s][i].s2;
		for (size_t j = 0; j < slinks2[s1].size(); ++j) {
		  if (slinks2[s1][j].s2 == s) {
		    slinks2[s1].Erase(j, j+1);
		    break;
		  }
		}
	      }
	    }
	  }
	  
	  // Identify the cases where of the scaffolds linked forward to by a
	  // given scaffold, one is clearly closest.
	  
	  for ( int s1 = 0; s1 < nscaff; s1++ ){    
	    if ( slinks2[s1].empty( ) ) continue;
	    int sep2 = slinks2[s1][0].sep;
	    size_t closest = 0;
	    for ( size_t i = 1; i < slinks2[s1].size(); i++ ){    
	      if ( slinks2[s1][i].sep < sep2 ){    
		sep2 = slinks2[s1][i].sep;
		closest = i;    
	      }    
	    }
	    int s2 = slinks2[s1][closest].s2;
	    int dev2 = slinks2[s1][closest].dev;
	    Bool fw2 = slinks2[s1][closest].fw2, OK = True;
	    for ( size_t i = 0; i < slinks2[s1].size(); i++ ){    
	      if ( i == closest ) continue;
	      int s2x = slinks2[s1][i].s2;
	      int sep2x = slinks2[s1][i].sep;
	      int dev2x = slinks2[s1][i].dev;
	      
	      // Make sure s2x could be to the right of s2, but could not
	      // be to the left of s2.
	      
	      const superb &S2 = scaffolds[s2], &S2x = scaffolds[s2x];
	      int s2_len = S2.SubSuperLength( 0, S2.Ntigs( ) - 1 );
	      int s2_dev = S2.SubSuperLengthDev( 0, S2.Ntigs( ) - 1 );
	      int s2x_len = S2x.SubSuperLength( 0, S2x.Ntigs( ) - 1 );
	      int s2x_dev = S2x.SubSuperLengthDev( 0, S2x.Ntigs( ) - 1 );
	      if ( sep2x + s2x_len - dev_mult * (s2x_dev + dev2x)
		   <= sep2 + s2_len){    
		OK = False;
		if (VERBOSITY >= 2 ){ 
		  cout << "Outgoing conflict on " << s1 << " between "; PRINT2(s2,s2x); 
		}
		if ( VERBOSITY < 3 )
		  break;
    	      }
	      if (VERBOSITY >= 3)
		cout << "    link s1=" << s1
		     << " s2x=" << s2x
		     << " sep2x=" << sep2x
		     << " dev2x=" << dev2x
		     << " OK=" << (OK?1:0) << endl;
	    }
	    if (VERBOSITY >= 3)
	      cout << "  link s1=" << s1
		   << " s2=" << s2
		   << " sep2=" << sep2
		   << " dev2=" << dev2
		   << " OK=" << (OK?1:0)
		   << " (closest)" << endl;
	    
	    if ( !OK ){
	      outConflictCt++;
	      continue;
	    }
	    // Record the links.
	    
	    if ( pass == 2 ){    
	      if (fw2) fw_fw[s1].push( s2, sep2, dev2 );
	      else fw_rc[s1].push( s2, sep2, dev2 );    
	    } else{    
	      if (fw2) rc_rc[s1].push( s2, sep2, dev2 );
	      else rc_fw[s1].push( s2, sep2, dev2 );    
	    }    
	  }    
	}
	
	if (VERBOSITY >= 2 )
	  PRINT2( pairsCt, propAlgndPairsCt );
	
	if ( VERBOSITY >= 2 ){    
	  for ( int i = 0; i < nscaff; i++ ){    
	    for ( size_t j = 0; j < fw_fw[i].size(); j++ )
	      cout << "fw_fw: " << i << " --> " << fw_fw[i][j] 
		   << "\tlengths: " << scaffolds[i].SubSuperLength(0, scaffolds[i].Ntigs() -1 ) 
		   << " " << scaffolds[fw_fw[i][j].first].SubSuperLength(0, scaffolds[fw_fw[i][j].first].Ntigs() -1 ) << "\n";
	    for ( size_t j = 0; j < fw_rc[i].size(); j++ )
	      cout << "fw_rc: " << i << " --> " << fw_rc[i][j] 
		   << "\tlengths: " << scaffolds[i].SubSuperLength(0, scaffolds[i].Ntigs() -1 ) 
		   << " " << scaffolds[fw_rc[i][j].first].SubSuperLength(0, scaffolds[fw_rc[i][j].first].Ntigs() -1 ) << "\n";
	    for ( size_t j = 0; j < rc_fw[i].size(); j++ )
	      cout << "rc_fw: " << i << " --> " << rc_fw[i][j] 
		   << "\tlengths: " << scaffolds[i].SubSuperLength(0, scaffolds[i].Ntigs() -1 ) 
		   << " " << scaffolds[rc_fw[i][j].first].SubSuperLength(0, scaffolds[rc_fw[i][j].first].Ntigs() -1 ) << "\n";
	    for ( size_t j = 0; j < rc_rc[i].size(); j++ )
	      cout << "rc_rc: " << i << " --> " << rc_rc[i][j] 
		   << "\tlengths: " << scaffolds[i].SubSuperLength(0, scaffolds[i].Ntigs() -1 ) 
		   << " " << scaffolds[rc_rc[i][j].first].SubSuperLength(0, scaffolds[rc_rc[i][j].first].Ntigs() -1 ) << "\n";   
	  }  
	}
	if ( !JOIN ) break;
	
	if ( VERBOSITY >= 2 ) 
	  PRINT( outConflictCt );
	// Make the scaffold links into a graph.  The vertices are 0fw, 0rc, 1fw,
	// 1rc, etc.  This is painful but purely mechanical.
	//if (VERBOSITY >=2 ) cout << Date() << " buiding graph" << endl;
	typedef pair<int,int> sepdev;
	vec< vec<int> > from(2*nscaff), to(2*nscaff);
	vec< vec<int> > from_edge_obj(2*nscaff), to_edge_obj(2*nscaff);
	vec<sepdev> edges;
	for ( int i = 0; i < nscaff; i++ ){    
	  for ( size_t j = 0; j < fw_fw[i].size(); j++ ){    
	    int nedges = edges.size( );
	    from[2*i].push_back( 2*fw_fw[i][j].first+1 );
	    to[ 2*fw_fw[i][j].first+1 ].push_back(2*i);
	    from_edge_obj[2*i].push_back( edges.size( ) );
	    to_edge_obj[ 2*fw_fw[i][j].first+1 ].push_back(nedges);
	    edges.push( fw_fw[i][j].second, fw_fw[i][j].third );    
	  }
	  for ( size_t j = 0; j < fw_rc[i].size(); j++ ){    
	    int nedges = edges.size( );
	    from[2*i].push_back( 2*fw_rc[i][j].first );
	    to[ 2*fw_rc[i][j].first ].push_back(2*i);
	    from_edge_obj[2*i].push_back( edges.size( ) );
	    to_edge_obj[ 2*fw_rc[i][j].first ].push_back(nedges);
	    edges.push( fw_rc[i][j].second, fw_rc[i][j].third );    
	  }
	  for ( size_t j = 0; j < rc_fw[i].size(); j++ ){    
	    int nedges = edges.size( );
	    from[2*i+1].push_back( 2*rc_fw[i][j].first+1 );
	    to[ 2*rc_fw[i][j].first+1 ].push_back(2*i+1);
	    from_edge_obj[2*i+1].push_back( edges.size( ) );
	    to_edge_obj[ 2*rc_fw[i][j].first+1 ].push_back(nedges);
	    edges.push( rc_fw[i][j].second, rc_fw[i][j].third );    
	  }
	  for ( size_t j = 0; j < rc_rc[i].size(); j++ ){    
	    int nedges = edges.size( );
	    from[2*i+1].push_back( 2*rc_rc[i][j].first );
	    to[ 2*rc_rc[i][j].first ].push_back(2*i+1);
	    from_edge_obj[2*i+1].push_back( edges.size( ) );
	    to_edge_obj[ 2*rc_rc[i][j].first ].push_back(nedges);
	    edges.push( rc_rc[i][j].second, rc_rc[i][j].third );    
	  }    
	}
	digraphE<sepdev> G( from, to, edges, to_edge_obj, from_edge_obj );
	
	if ( VERBOSITY >= 2 ){ 
	  cout << "number of edges in the graph "; PRINT(edges.size());
	}
	
	// Cleanup strand bidiriectionals
	for ( int s1 = 0; s1 < nscaff*2 -1; s1++ ){
	  for ( int s2 = s1 + 1; s2 < nscaff*2; s2++ ){
	    if ( G.EdgesBetween(s1,s2).size() > 0 ){
	      
	      if ( G.EdgesBetween(s2,s1).size() > 0 ){
		G.DeleteEdges(  G.EdgesBetween(s1,s2) );
		G.DeleteEdges(  G.EdgesBetween(s2,s1) );
		if (VERBOSITY >= 2){
		  PRINT( G.EdgesBetween(s1,s2).size() );
		  PRINT( G.EdgesBetween(s2,s1).size() );
		  cout << "DELETING BIDIRECTIONAL: edges between s1 = " 
		       << s1 << " and s2 = " << s2 << " were deleted" << endl;
		}
	      }
	    }	    
	  }
	}
	
	// Cleanup incoming (symmetrize): check for conflicts, chose closest link	
	if ( VERBOSITY >= 2 ) 
	  cout << Date() << " Cleaning up graph" << endl;    
	for ( int s = 0; s < nscaff*2; s++ ){
	  if ( G.To(s).size() > 1 ){
	    vec< pair<int,sepdev> > links( G.To(s).size() );
	    for ( size_t i = 0; i < G.To(s).size(); i++ ){
	      int s0 = G.To(s).at(i);
	      ForceAssert( G.EdgeObjectsBetween(s0,s).size() ==1 );
	      links[i] = pair<int,sepdev>(s0, G.EdgeObjectsBetween(s0,s).at(0)); 
	    }
	    
	    //find closest
	    ForceAssert( links.size() >= 1 );
	    size_t closest = 0;
	    int sep0 = links.at(0).second.first;
	    for ( size_t i = 1; i < links.size(); i++ ){
	      if ( links.at(i).second.first < sep0 ){
		sep0    = links.at(i).second.first;
		closest = i;
	      }
	    }
	    
	    Bool OK = True;
	    int s0 = links.at(closest).first;
	    for ( size_t i = 0; i < links.size(); i++ ){
	      if ( i == closest ) continue;
	      int s0x = links[i].first;
	      int sep0x = links[i].second.first;
	      int dev0x = links[i].second.second;
	      
	      const superb &S0 = scaffolds.at( s0/2 );
	      const superb &S0x = scaffolds.at( s0x/2 );
	      int s0_len = S0.SubSuperLength( 0, S0.Ntigs( ) - 1 );
	      int s0_dev = S0.SubSuperLengthDev( 0, S0.Ntigs( ) - 1 );
	      int s0x_len = S0x.SubSuperLength( 0, S0x.Ntigs( ) - 1 );
	      int s0x_dev = S0x.SubSuperLengthDev( 0, S0x.Ntigs( ) - 1 );
	      if ( sep0x + s0x_len - dev_mult * (s0x_dev + dev0x)
		   <= sep0 + s0_len){    
		OK = False;
		if ( VERBOSITY >= 2 ){ 
		  cout << "Incoming conflict on " << s << " between "; PRINT2(s0,s0x);
		}
		inConflictCt++;
		break; // one conflict is enough to stop
	      }
	    }
	    if ( ! OK ){
	      // remove connections
	      if ( VERBOSITY >= 2 )
		cout << " removing all connections to scaffold vertex " << s <<  endl;
	      for ( size_t i = 0; i < links.size(); i++ ){
		int s0 = links[i].first;
		vec<int> linkEdgesTo = G.EdgesBetween(s0, s);
		ForceAssert( linkEdgesTo.size() == 1 );
		G.DeleteEdges( linkEdgesTo );
	      }
	    }else{
	      // remove all except the nearest
	      if ( VERBOSITY >= 2 )
		cout << " removing distant connections to scaffold vertex " << s <<  endl;
	      for ( size_t i = 0; i < links.size(); i++ ){
		if ( i == closest ) continue;
		int s0 = links[i].first;
		vec<int> linkEdgesTo = G.EdgesBetween(s0, s);
		ForceAssertEq( 1u, linkEdgesTo.size() );
		G.DeleteEdges( linkEdgesTo );
	      }
	      
	    } 
	  }
	}

	if ( VERBOSITY >= 2 )
	  PRINT( inConflictCt );
	
	
	// Now use the links to define lines in the scaffold graph.  There are
	// heuristics embedded here.
	
	vec< vec<int> > lines;
	vec<Bool> used( nscaff, False );
	for ( int pass = 1; pass <= 2; pass++ ){    
	  for ( int s0 = 0; s0 < nscaff*2; s0++ ){
	    if ( VERBOSITY >= 2 ){
	      if ( G.From(s0).size() > 1 )
		cout << "s0: G.From(" << s0 << ").size() = " << G.From(s0).size() << endl;
	      if ( G.To(s0).size() > 1 )
		cout << "s0: G.To(" << s0 << ").size() = " << G.To(s0).size() << endl;
	    }
	    int s = s0;
	    if ( used[s/2] || ( pass == 1 && !G.Source(s) ) ) continue;
	    if ( pass == 1 && G.From(s0).size() < 1 ) continue;
	    vec<int> line;
	    line.push_back(s);
	    used[s/2] = True;
	    int extendedCt = 0;
	    while(1){    
	      if ( VERBOSITY >= 2 ){
		if ( G.From(s).size() > 1 )
		  cout << "G.From(" << s << ").size() = " << G.From(s).size() << endl;
		if ( G.To(s).size() > 1 )
		  cout << "G.To(" << s << ").size() = " << G.To(s).size() << endl;
	      }
	      if ( G.From(s).size() != 1 ) break;
	      s = G.From(s)[0];
	      if ( G.To(s).size() > 1 || used[s/2] ) break;
	      used[s/2] = True;
	      line.push_back(s); 
	      extendedCt++;
	    }
	    lines.push_back(line);    
	  }    
	}
	
	// Merge the scaffolds in the lines to form new scaffolds.  This is a
	// purely formal operation.
	
	vec<superb> new_scaffolds;
	vec<Bool> reversed( ntigs, False );
	if ( VERBOSITY >= 1 ){    
	  cout << "\nlines now:\n";
	  for ( size_t i = 0; i < lines.size(); i++ ){    
	    for ( size_t j = 0; j < lines[i].size(); j++ ){    
	      if ( j > 0 ) cout << " --> ";
	      cout << lines[i][j]/2
		   << ( lines[i][j] % 2 == 0 ? "fw" : "rc" ) << " (" << scaffolds[ lines[i][j]/2 ].SubSuperLength(0, scaffolds[ lines[i][j]/2 ].Ntigs() -1 ) << ")";    
	    }
	    cout << "\n";    
	  }    
	}
	for ( size_t i = 0; i < lines.size(); i++ ){    
	  joins += lines[i].size( ) - 1;
	  superb x;
	  int ntigs = 0;
	  for ( size_t j = 0; j < lines[i].size(); j++ )
	    ntigs += scaffolds[ lines[i][j]/2 ].Ntigs( );
	  x.SetNtigs(ntigs);
	  int count = 0;
	  for ( size_t j = 0; j < lines[i].size(); j++ ){ 
	    int y = lines[i][j];
	    int s = y/2;
	    superb S = scaffolds[s];
	    
	    // In the rc case, reverse S, and also track the change.
	    
	    if ( y % 2 == 1 ){    
	      S.Reverse( );
	      for ( int k = 0; k < S.Ntigs( ); k++ ){    
		rctig[ S.Tig(k) ] = !rctig[ S.Tig(k) ];
		reversed[ S.Tig(k) ] = True;    
	      }    
	    }
	    
	    // Add S to the scaffold.
	    
	    for ( int k = 0; k < S.Ntigs( ); k++ ){    
	      x.SetTig( count, S.Tig(k) );
	      x.SetLen( count, S.Len(k) );
	      if ( k < S.Ntigs( ) - 1 ){    
		x.SetGap( count, S.Gap(k) );
		x.SetDev( count, S.Dev(k) );
		count++;    
	      }    
	    }
	    
	    // Add the between scaffold gap.
	    
	    if ( j+1 < lines[i].size() ){    
	      x.SetGap( count, G.EdgeObjectByIndexFrom( y, 0 ).first );
	      x.SetDev( count, G.EdgeObjectByIndexFrom( y, 0 ).second );    
	    }
	    count++;    
	  }
	  new_scaffolds.push_back(x);    
	}
	
	// Reverse the alignments lying over contigs that have just
	// been reversed.
	
	for ( size_t i = 0; i < aligns0.size(); i++ )
	  if ( reversed[ aligns0[i].TargetId( ) ] ) aligns0[i].Reverse( );
	
	// Check to see if done.
	
	scaffolds = new_scaffolds;
	if ( VERBOSITY >= 1 )
	  cout << Date( ) << ": made " << joins << " joins" << endl;
	allJoinsCt += joins;
	localJoinsCt += joins;
      } while (joins > 0);
      
      // Absorb small scaffolds into larger ones.
      
      if ( SUCK_SCAFFOLDS ) {
	ofstream devnull ( "/dev/null" );
	ostream *plog = ( VERBOSITY > 0 ? &cout : &devnull );
	int n_sucked = SuckScaffolds( pairs, aligns0_index, aligns0,
				      scaffolds, rctig, plog );
	cout << Date( ) << ": " << n_sucked << " scaffolds sucked" << endl;
      }
      
    }
    
  }
  if ( VERBOSITY >= 1 )
    cout << "total number of joins made = " << allJoinsCt << endl;
  
  // Flip rc contigs. 
  
  cout << Date() << ": flipping rc contigs in the fastavector" << endl;
  for ( size_t i = 0; i < ntigs; i++ )
    if ( rctig[i] )
      scaffold_contigs[i].ReverseComplement();
  
  // Report some stats.
  
  cout << "\nEdges in HyperFastavector: " << ntigs << endl;
  vec<int> contig_sizes, scaffold_sizes;
  for ( size_t i = 0; i < scaffolds.size(); i++ ){    
    for ( int j = 0; j < scaffolds[i].Ntigs( ); j++ )
      contig_sizes.push_back(  scaffolds[i].Len(j) );
    scaffold_sizes.push_back( scaffolds[i].ReducedLength( ) );    
  }
  Sort(contig_sizes), Sort(scaffold_sizes);
  contig_sizes.EraseValue( 0 );
  scaffold_sizes.EraseValue( 0 );
  
  cout << "N50 contig size = " << N50(contig_sizes) 
       << " (" << contig_sizes.size() << " contigs)" << endl;
  cout << "N50 scaffold size (without gaps) = " << N50(scaffold_sizes) 
       << " (" << scaffold_sizes.size() << " scaffolds)" << endl;
  cout << "total bases in contigs = " << BigSum(scaffold_sizes) << endl;
  
  // Save results.  We create a separate vector to track orientation of the
  // contigs.  It would be better to encapsulate this in class superb.
  
  if ( !WRITE ) exit(0);
  
  cout << "\n" << Date( ) << ": writing output files" << endl;
     
  WriteSuperbs( scaffolds_out_file , scaffolds );
  Ofstream( contig_fasta, contigs_out_file );
  for ( size_t i = 0; i < ntigs; i++ )
    scaffold_contigs[i].Print( contig_fasta, "contig_" + ToString(i) );
  contig_fasta.close();
  
  WriteScaffoldedFasta( assembly_out_file, scaffold_contigs, scaffolds );
  
  cout << Date( ) << ": done" << endl;
  
}
