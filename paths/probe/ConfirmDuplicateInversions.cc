///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "CoreTools.h"
#include "MainTools.h"
#include "FetchReads.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "FastaFileset.h"
#include "paths/ReadLoc.h"

#include <map>
#include <omp.h>
// MakeDepend: library OMP
// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable

class pdata{

public:
  size_t id1_, id2_;
  bool fw1_, fw2_;
  short read_class_;

  pdata( uint64_t id1, bool fw1, uint64_t id2, bool fw2 ){
    id1_ = id1; id2_ = id2;
    fw1_ = fw1; fw2_ = fw2;
  }
  
};

Bool cmp_target( const look_align& la1, const look_align& la2 ){    
  if ( la1.target_id < la2.target_id ) return True;
  if ( la1.target_id > la2.target_id ) return False;
  if ( la1.pos2( ) < la2.pos2( ) ) return True;
  if ( la1.pos2( ) > la2.pos2( ) ) return False;
  if ( la1.query_id < la2.query_id ) return True;
  return False;    
}

struct cmp_errors : public std::binary_function<look_align,look_align,bool>
{
    bool operator()( const look_align& la1, const look_align& la2 ) const
    { return  la1.Errors() < la2.Errors(); }
};

/**
 * ConfirmDuplicateInversions
 *
 
*/


int main( int argc, char *argv[] ){
  
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  
  CommandArgument_String( ASSEMBLY );
  CommandArgument_String( REF);
  CommandArgument_String( INVERSIONS);

  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Bool_OrDefault( REBUILD_LOOKUP, False );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  EndCommandArguments;

  // Thread control (OMP used in FixScaffoldsCore)
  
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  String data_dir = PRE + "/" + DATA + "/";
  String run_dir  = data_dir + RUN + "/";
  String subdir   = run_dir + "/ASSEMBLIES/" + SUBDIR + "/";


  cout << Date() << ": reading contigs file" << endl;
  String contigs_file = subdir + ASSEMBLY.SafeBeforeLast(".assembly") + ".contigs.fastb";
  ForceAssert( IsRegularFile( contigs_file ));
  vecbasevector contigs( contigs_file );

  cout << Date() << ": reading superbs" << endl;
  String superb_file = subdir +  ASSEMBLY.SafeBeforeLast(".assembly") + ".superb";
  ForceAssert( IsRegularFile( superb_file ));
  vec<superb> supers;
  ReadSuperbs( superb_file, supers);

  const int n_contigs = contigs.size( );
  
  // Maps.
  cout << Date() << ": computing maps" << endl;
  vec<int> super_len( supers.size( ), 0 );
  vec<int> start_on_super( n_contigs );
  vec<int> to_super( n_contigs );
  vec<int> tig_len( n_contigs );
  for (int ii=0; ii<supers.isize( ); ii++) {
    int pos = 0;
    super_len[ii] = supers[ii].TrueLength( );
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++) {
      start_on_super[ supers[ii].Tig( jj ) ] = pos;
      to_super[ supers[ii].Tig( jj ) ] = ii;
      tig_len[ supers[ii].Tig( jj ) ] = supers[ii].Len( jj );
      pos += supers[ii].Len( jj );
      if ( jj < supers[ii].Ntigs( ) - 1 ) pos += supers[ii].Gap( jj );
    }
  }
  
  // Read in inversion coordinates
  
  String inversions_file = subdir + "tmpProbeInv4/" + INVERSIONS;
  vec< vec<int> > inversions_in;
  {
    ifstream ifin( inversions_file.c_str() );
    String line;
    while( getline(ifin,line) ){
      vec<String> fields;
      Tokenize( line, '\t', fields ); //tig1 \t start1 \t end1 \t tig2 \t start2 \t end2\n";
      //PRINT(line);
      vec<int> values;
      for ( size_t i = 0; i < fields.size(); i++ ){
	values.push_back( atoi(fields[i].c_str()) );
      }
      inversions_in.push_back( values );
    }
  }



  // Get read data associated with inversions

  String assembly_head = subdir + ASSEMBLY.SafeAfterLast("/").SafeBefore(".fast");
  PRINT(assembly_head);
  
  cout << Date() << ": preparing readloc reading" << endl;
  read_locs_on_disk locs_file( assembly_head.SafeBeforeLast(".assembly"), run_dir );

  PRINT( inversions_in.size() );
  vec<size_t> rids;
  vec< vec< vec<pdata> > > assem_data( inversions_in.size(), vec< vec<pdata> > (2) );
  for ( size_t i = 0; i < inversions_in.size(); i++ ){
    for ( int j = 0; j <= 1; j++ ){
      int s1 = inversions_in[i][0+j*3], start1 = inversions_in[i][1+j*3], end1 = inversions_in[i][2+j*3];
      int tig1 = -1;
      for ( size_t ti = 0; ti < to_super.size(); ti++ )
	if ( to_super[ti] == s1 && start_on_super[ti] <= start1 && start_on_super[ti] + tig_len[ti] <= end1 )
	  tig1 = ti;
      ForceAssertGe( tig1, 0 );
      vec<read_loc> locs;
      locs_file.LoadContig( tig1, locs );
      for (int loc_id=0; loc_id<locs.isize( ); loc_id++) {
	const read_loc &rloc = locs[loc_id];
	if ( rloc.ReadClass() != 1 ) continue;
	if ( ! rloc.PartnerPlaced() ) continue;
	int tig2 = rloc.PartnerContigId();
	if ( to_super.at(tig1) != to_super.at(tig2) ) continue;
	int r1s   = start_on_super[tig1] + rloc.Start();
	int r1e   = start_on_super[tig1] + rloc.Stop();
	int r2s   = start_on_super[tig2] + rloc.PartnerStart();
	int r2e   = start_on_super[tig2] + rloc.PartnerStop();
	if ( ( r1s >= start1 && r1e <= end1 && ( r2s > end1 || r2e <= start1 ) ) ||
	     ( r2s >= start1 && r2e <= end1 && ( r1s > end1 || r1e <= start1 ) )  
	     ){
	  //PRINT6( start1, end1, r1s, r1e, r2s, r2e );
	  size_t id1 = rloc.ReadId();
	  size_t id2 = rloc.PartnerReadId();
	  bool fw1   = rloc.Fw();
	  bool fw2   = rloc.PartnerFw();
	  pdata data_loc( id1, fw1, id2, fw2 );
	  assem_data[i][j].push_back( data_loc );
	  rids.push_back( id1, id2 );
	}
      }
    }
  }
  
  UniqueSort(rids);
  String ids_file = subdir + "tmpProbeInv4/" + assembly_head.SafeAfterLast("/") + ".ids_to_process";
  ofstream iout( ids_file.c_str() );
  for ( size_t i = 0; i < rids.size(); i++ )
    iout << rids[i] << "\n";
  iout.close();
  


  // get read alignments

  String fastb_file = run_dir + "jump_reads_filt_cpd.fastb";
  String aligns_file = subdir + "tmpProbeInv4/selectReads.qltout";
  String lookup_file = data_dir + REF.SafeBefore(".fast") + ".lookup";
  
  if ( ! IsRegularFile( lookup_file ) ){
    cout << Date() << ": Creating a lookup table for reference file" << endl;
    SystemSucceed( "MakeLookupTable" + ARG(SOURCE, REF )
		   + ARG(OUT_HEAD, lookup_file.SafeBefore(".lookup")) + ARG(LOOKUP_ONLY, True)
		   + ARG( NH, True ) );
  }

  String tmp_dir = "./";
  cout << Date() << ": aligning reads" << endl;
  SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.3 MF=5000 " 
		 + ARG(SEQS, fastb_file) 
		 + ARG(SEQS_IS_FASTB, True) 
		 + ARG(L, lookup_file) + ARG(PARSEABLE, True) 
		 + ARG(SMITH_WAT, True) + ARG(TMP_DIR,tmp_dir)
		 + ARG(KEEP_BEST,2)
		 + ARG(SEQS_TO_PROCESS, "@" + ids_file )
		 + ARG(OUTPUT, aligns_file) );
		 
  
  vec<look_align> aligns;
  LoadLookAligns( aligns_file, aligns );
  sort( aligns.begin( ), aligns.end( ), cmp_target );  
  
  // build index
  map<size_t,vec<size_t> > aligns_index;
  for ( size_t i = 0; i < aligns.size(); i++ ){
    size_t qid = aligns[i].QueryId();
    aligns_index[qid].push_back(i);
  }
    
  for ( map<size_t,vec<size_t> >::iterator it = aligns_index.begin(); it != aligns_index.end(); it++ ){
    size_t qid = it->first;
    vec<size_t>& indices = it->second;

    vec<look_align> alignsLoc( indices.size() );
    for ( size_t i = 0; i < indices.size(); i++ ){
      alignsLoc[i] = aligns[ indices[i] ];
    }
    SortSync( alignsLoc, indices, cmp_errors() );
    int minErr = alignsLoc[0].Errors();
    size_t i_cut = alignsLoc.size();
    for ( size_t i = 0; i < alignsLoc.size(); i++ )
      if ( alignsLoc[i].Errors() != minErr ){
	indices.resize(i);
	break;
      }
  }


  for ( size_t ii = 0; ii < assem_data.size(); ii++ ){
    cout << "\nInversion number = " << ii << endl;
    inversions_in[ii].Print(cout); cout << endl;
    int opp_a1 = 0, along_a1 = 0;
    int opp_r1 = 0, along_r1 = 0;
    PRINT(  assem_data[ii][0].size() );
    for ( size_t k = 0; k < assem_data[ii][0].size(); k++ ){
      pdata& d = assem_data[ii][0][k];
      if ( d.fw1_ != d.fw2_ ) opp_a1++;
      else along_a1++;

      for ( size_t m = 0; m < aligns_index[d.id1_].size(); m++ ){
	look_align& la1 = aligns[ aligns_index[d.id1_][m] ];
	for ( size_t n = 0; n < aligns_index[d.id2_].size(); n++ ){
	  look_align& la2 = aligns[ aligns_index[d.id2_][n] ];
	  if ( la1.Fw1() != la2.Fw1() ) opp_r1++;
	  else along_r1++;
	}
      }
    }
    int opp_a2 = 0, along_a2 = 0;
    int opp_r2 = 0, along_r2 = 0;
    PRINT(  assem_data[ii][1].size() );
    for ( size_t k = 0; k < assem_data[ii][1].size(); k++ ){
      pdata& d = assem_data[ii][1][k];
      if ( d.fw1_ != d.fw2_ ) opp_a2++;
      else along_a2++;

      for ( size_t m = 0; m < aligns_index[d.id1_].size(); m++ ){
	look_align& la1 = aligns[ aligns_index[d.id1_][m] ];
	for ( size_t n = 0; n < aligns_index[d.id2_].size(); n++ ){
	  look_align& la2 = aligns[ aligns_index[d.id2_][n] ];
	  if ( la1.Fw1() != la2.Fw1() ) opp_r2++;
	  else along_r2++;
	}
      }
    }


    PRINT8( opp_a1, along_a1, opp_r1, along_r1, opp_a2, along_a2, opp_r2, along_r2 );

  }
  
  cout << Date( ) << ": done" << endl;

}


