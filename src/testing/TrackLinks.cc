///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// TrackLinks.  Display the links from scaffold S to other scaffolds.  The argument
// MAX_OVERLAP determines how large the implies overlap between S and the other
// scaffolds can be - in order to display the links.
// 
// Derived from ShowLinks. Now we check if the linking information also exist 
// in the initial scaffolds

#include "Equiv.h"
#include "MainTools.h"
#include "CoreTools.h"
#include "Superb.h"
#include "math/Functions.h"
#include "paths/ReadLoc.h"
#include "system/ParsedArgs.h"
#include "util/ReadTracker.h"
#include "paths/Alignlet.h"
#include "paths/MapAlignletToRef.h"
#include "feudal/VirtualMasterVec.h"
#include "feudal/BinaryStream.h"

// the datastructure to represent how a read is aligned to refrence genome
// multiple alignment is possible because repeats
// status code to check if the align cannot found
struct AlignToRef{
  int status;
  vec<alignlet> aligns;
};

class Link {
  public:
    int64_t readloc_id;
    longlong read_id1;
    longlong read_id2;
    int m1;
    int s1;
    int m2;
    int s2;
    int sep;
    int dev;
    int read_class;
    int dist_to_end1;
    int dist_to_end2;
    Bool rc1;
    Bool rc2;

    Link( ) { }
    Link(const longlong read_id1, const longlong read_id2,
       	const int m1, const int s1, const int m2, const int s2, const int sep, 
	const int dev, const int read_class, const int dist_to_end1, 
	const int dist_to_end2, const Bool rc1, const Bool rc2) : 
       read_id1(read_id1), read_id2(read_id2), m1(m1), s1(s1),
    m2(m2), s2(s2), sep(sep), dev(dev), read_class(read_class), 
    dist_to_end1(dist_to_end1), dist_to_end2(dist_to_end2), rc1(rc1), 
    rc2(rc2) { }

    friend Bool operator<( const Link& l1, const Link& l2 )
    {    return l1.s2 < l2.s2;     }


};

void PrintRefAlign(ostream& out, vec<alignlet> aligns) {
  if (aligns.empty()) {
    out << "NA";
    return;
  }
  String sign = aligns[0].Fw1()? "+" : "-";
  int pos =     aligns[0].Fw1()? aligns[0].pos2() : aligns[0].Pos2();
  out<< aligns[0].pos2()<< sign <<aligns[0].TargetId();
  for(size_t ii=1;ii<aligns.size();ii++){
    String sign = aligns[ii].Fw1()? "+" : "-";
    int pos =     aligns[ii].Fw1()? aligns[ii].pos2() : aligns[ii].Pos2();
    out<< "&"<< pos << sign <<aligns[ii].TargetId();
  }
}

int main(int argc, char *argv[])
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String_OrDefault(SUBDIR, "test");
  CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
  CommandArgument_Int(SC1);
  CommandArgument_Int_OrDefault(SC2,-1);
  CommandArgument_String_OrDefault(EC_READS, "_ec");
  CommandArgument_String_OrDefault(INITIAL_ASSEMBLY, "initial_scaffolds");
  CommandArgument_String_OrDefault(ALIGNLET, "scaffold_reads");
  CommandArgument_String_OrDefault_Doc(READLOCS_PREFIX, "",
      "if specified, file extension is .READLOCS_PREFIX.readlocs "
      "instead of .readlocs");
  CommandArgument_Int_OrDefault(MAX_OVERLAP, 10000);
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  EndCommandArguments;

  // Define directories.

  String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 
  String eval_dir = sub_dir + "/" + ASSEMBLY + ".Eval";
  String contigs_aligns_file = eval_dir + "/aligns.qlt";
  String eval_dir0 = sub_dir + "/" + INITIAL_ASSEMBLY+ ".Eval";
  String contigs_aligns_file0 = eval_dir0 + "/aligns.qlt";

  // The contigs
  String contig_fn = sub_dir + "/" + ASSEMBLY + ".contigs.fastb" ;
  String contig_fn0 = sub_dir + "/" + INITIAL_ASSEMBLY + ".contigs.fastb" ;
  int ntigs = MastervecFileObjectCount( contig_fn );
  int ntigs0 = MastervecFileObjectCount( contig_fn0 );
  VirtualMasterVec<bvec> contigs(contig_fn.c_str());
  VirtualMasterVec<bvec> contigs0(contig_fn0.c_str());

  // load the alignment of scaffold to reference
  cout << Date( ) << ": loading aligns of contigs onto reference" << endl;
  vec<look_align_plus> cg_hits;
  LoadLookAlignPlus( contigs_aligns_file, cg_hits );
  // Maps.
  vec< vec<int> > to_hits( ntigs );   // contig_id -> ids in cg_hits
  for (int ii=0; ii<cg_hits.isize( ); ii++)
    to_hits[ cg_hits[ii].query_id ].push_back( ii );

  cout << Date( ) << ": loading aligns of contigs0 onto reference" << endl;
  vec<look_align_plus> cg_hits0;
  LoadLookAlignPlus( contigs_aligns_file0, cg_hits0 );
  // Maps.
  vec< vec<int> > to_hits0( ntigs0 );   // contig_id -> ids in cg_hits
  for (int ii=0; ii<cg_hits0.isize( ); ii++)
    to_hits0[ cg_hits0[ii].query_id ].push_back( ii );

  // Load the scaffolds
  vec<superb> scaffolds;
  ReadSuperbs( sub_dir + "/" + ASSEMBLY + ".superb", scaffolds );
  vec<int> to_super( ntigs, -1 ), to_super_pos( ntigs, -1 );
  vec<int> to_super_posr( ntigs, - 1 );
  for ( int i = 0; i < scaffolds.isize( ); i++ )
  {    int n = scaffolds[i].Ntigs( );
    for ( int j = 0; j < n; j++ )
    {    to_super[ scaffolds[i].Tig(j) ] = i;
      to_super_pos[ scaffolds[i].Tig(j) ] = j;   
      to_super_posr[ scaffolds[i].Tig(j) ] = n - j - 1;    
    }    
  }

  // Load scaffolds0.
  vec<superb> scaffolds0;
  ReadSuperbs( sub_dir + "/" + INITIAL_ASSEMBLY + ".superb", scaffolds0 );
  vec<int> to_super0( ntigs0, -1 ), to_super_pos0( ntigs0, -1 );
  vec<int> to_super_posr0( ntigs0, - 1 );
  for ( int i = 0; i < scaffolds0.isize( ); i++ )
  {    int n = scaffolds0[i].Ntigs( );
    for ( int j = 0; j < n; j++ )
    {    to_super0[ scaffolds0[i].Tig(j) ] = i;
      to_super_pos0[ scaffolds0[i].Tig(j) ] = j;   
      to_super_posr0[ scaffolds0[i].Tig(j) ] = n - j - 1;    
    }    
  }

  // Load the filt_cpd reads
  String jump_reads_filt_fn = run_dir + "/jump_reads_filt_cpd.fastb";
  String long_jump_reads_filt_fn = run_dir + "/long_jump_reads_filt.fastb";
  longlong n_jump_reads_filt = IsRegularFile(jump_reads_filt_fn) ?
    MastervecFileObjectCount(jump_reads_filt_fn) : 0;
  cout << "n_jump_reads_filt "<< n_jump_reads_filt <<  endl;
  longlong n_long_jump_reads_filt = IsRegularFile(long_jump_reads_filt_fn) ?
    MastervecFileObjectCount(long_jump_reads_filt_fn) : 0 ;
  cout << "n_long_jump_reads_filt "<< n_long_jump_reads_filt <<  endl;

  // Load the ec reads
  String jump_reads_ec_fn = run_dir + "/jump_reads" + EC_READS + ".fastb";
  String long_jump_reads_ec_fn = run_dir + "/long_jump_reads" + EC_READS + ".fastb";
  longlong n_jump_reads_ec = IsRegularFile(jump_reads_ec_fn) ?
    MastervecFileObjectCount(jump_reads_ec_fn) : 0;
  cout << "n_jump_reads_ec "<< n_jump_reads_ec <<  endl;
  longlong n_long_jump_reads_ec = IsRegularFile(long_jump_reads_ec_fn) ?
    MastervecFileObjectCount(long_jump_reads_ec_fn) : 0;
  cout << "n_long_jump_reads_ec "<< n_long_jump_reads_ec <<  endl;

  // Load the qltoutlets
  String qlt_fn =  sub_dir + "/"+ ALIGNLET +".qltoutlet";
  String qlt_index_fn =  qlt_fn +".index";
  ForceAssert( IsRegularFile(qlt_fn) );
  ForceAssert( IsRegularFile(qlt_index_fn) );
  // Load the _filtered.qltoutlets
  String qlt_filt_fn =  sub_dir + "/"+ ALIGNLET +"_filtered.qltoutlet";
  String qlt_filt_index_fn =  qlt_fn +".index";
  ForceAssert( IsRegularFile(qlt_filt_fn) );
  ForceAssert( IsRegularFile(qlt_filt_index_fn) );

  // Built mapper from filt_cpd reads to ec reads
  // Read the mapping file if already generated
  vec<longlong> jump_mapping;
  vec<longlong> long_jump_mapping;
  String jump_mapping_file = "tmp_jump_reads_filt2ec.mapping";
  String long_jump_mapping_file = "tmp_long_jump_reads_filt2ec.mapping";
  if ( ! IsRegularFile(jump_mapping_file) ) {
    String jump_ec_tracker_fn =  run_dir + "/jump_reads" + EC_READS ;
    ReadTracker jumpTr;
    cout << "load " <<jump_ec_tracker_fn << ".readtrack" << endl;
    jumpTr.Load(jump_ec_tracker_fn);
    cout << "jump_reads Readtrack size: " << jumpTr.size() << endl;
    jump_mapping.resize(n_jump_reads_filt,-1);
    for(uint64_t i=0;i<jumpTr.size();i++)
    {
      uint64_t read_id_filt = jumpTr.GetReadIndex(i);
      jump_mapping[read_id_filt] = i;
    }
    BinaryWriter::writeFile(jump_mapping_file, jump_mapping);
  }
  if ( ! IsRegularFile(long_jump_mapping_file) ) {
    String long_jump_ec_tracker_fn =  run_dir + "/long_jump_reads" + EC_READS;
    ReadTracker long_jumpTr;
    cout << "load " <<long_jump_ec_tracker_fn << ".readtrack" << endl;
    long_jumpTr.Load(long_jump_ec_tracker_fn);
    cout << "long_jump_reads Readtrack size: " << long_jumpTr.size() << endl;
    long_jump_mapping.resize(n_long_jump_reads_filt,-1);
    for(uint64_t i=0;i<long_jumpTr.size();i++)
    {
      uint64_t read_id_filt = long_jumpTr.GetReadIndex(i);
      long_jump_mapping[read_id_filt] = i;
    }
    BinaryWriter::writeFile(long_jump_mapping_file, long_jump_mapping);
  }

  // Process contigs.
  int s1 = SC1;
  int s2 = SC2;
  const superb& S1 = scaffolds[s1];
  String head = sub_dir + "/" + ASSEMBLY;
  if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
  read_locs_on_disk locs_file( head, run_dir );
  vec<Link> Links;
  vec<read_loc> readlocs; // save the readlocs for further analysis
  for ( int p1 = 0; p1 < scaffolds[s1].Ntigs( ); p1++ )
  {    
    int m1 = scaffolds[s1].Tig(p1);
    vec<read_loc> locs;
    locs_file.LoadContig( m1, locs );
    for ( int j = 0; j < locs.isize( ); j++ )
    {
      const read_loc& rl = locs[j];
      if ( !rl.PartnerPlaced( ) ) continue;
      if ( rl.Frag() ) continue; // only jump and long jump reads
      int x1 = rl.Start( );
      int y1 = scaffolds[s1].SubSuperLength( 0, p1 - 1 ) + x1;
      int m2 = rl.PartnerContigId( );
      int s2 = to_super[m2], p2 = to_super_pos[m2];
      if ( s2 < 0 ) continue;
      const superb& S2 = scaffolds[s2];
      if ( s2 == s1 ) continue;
      // only selected links
      if ( SC2 > 0 && s2 != SC2 ) continue;

      int x2 = rl.PartnerStart( );
      int y2 = scaffolds[s2].SubSuperLength( 0, p2 - 1 ) + x2;
      // get the read_id and read_id_p
      longlong read_id = rl.ReadId();
      longlong read_id_p = rl.PartnerReadId();

      // For now just consider 1/2 of the cases:
      int sep=-1, dev=-1;
      int dist_to_end1=-1;
      int dist_to_end2=-1;
      if ( rl.Fw( ) && rl.PartnerRc( ) )
      {    
	dist_to_end1 = S1.Len(p1) - rl.Stop( ) + S1.SubSuperLength( p1 + 1, S1.Ntigs( ) - 1 );
	dist_to_end2 = rl.PartnerStart( ) + S2.SubSuperLength( 0, p2 - 1 );
	sep = rl.Sep( ) - dist_to_end1 - dist_to_end2;
	dev = rl.Dev( );
	if ( sep < -MAX_OVERLAP ) continue;
	Links.push(Link(read_id, read_id_p, m1, s1, m2, s2, sep, dev, rl.ReadClass( ),
	      dist_to_end1, dist_to_end2, False, True) );
      }
      if ( rl.Fw( ) && rl.PartnerFw( ) )
      {    
	dist_to_end1 = S1.Len(p1) - rl.Stop( ) + S1.SubSuperLength( p1 + 1, S1.Ntigs( ) - 1 );
	dist_to_end2 = S2.Len(p2) - rl.PartnerStop( ) + S2.SubSuperLength( p2 + 1, S2.Ntigs( ) - 1 );
	sep = rl.Sep( ) - dist_to_end1 - dist_to_end2;
	dev = rl.Dev( );
	if ( sep < -MAX_OVERLAP ) continue;
	Links.push(Link(read_id, read_id_p, m1, s1, m2, s2, sep, dev, rl.ReadClass( ),
	      dist_to_end1, dist_to_end2, False, False) );
      }
      if ( rl.Rc( ) && rl.PartnerRc( ) )
      {    
	dist_to_end1 = rl.Start( ) + S1.SubSuperLength( 0, p1 - 1 );
	dist_to_end2 = rl.PartnerStart( ) + S2.SubSuperLength( 0, p2 - 1 );
	sep = rl.Sep( ) - dist_to_end1 - dist_to_end2;
	dev = rl.Dev( );
	if ( sep < -MAX_OVERLAP ) continue;
	Links.push(Link(read_id, read_id_p, m1, s1, m2, s2, sep, dev, rl.ReadClass( ),
	      dist_to_end1, dist_to_end2, True, True) );
      }    
      if ( rl.Rc( ) && rl.PartnerFw( ) )
      {    
	dist_to_end1 = rl.Start( ) + S1.SubSuperLength( 0, p1 - 1 );
	dist_to_end2 = S2.Len(p2) - rl.PartnerStop( ) + S2.SubSuperLength( p2 + 1, S2.Ntigs( ) - 1 );
	sep = rl.Sep( ) - dist_to_end1 - dist_to_end2;
	dev = rl.Dev( );
	if ( sep < -MAX_OVERLAP ) continue;
	Links.push(Link(read_id, read_id_p, m1, s1, m2, s2, sep, dev, rl.ReadClass( ),
	      dist_to_end1, dist_to_end2, True, False) );
      }    
      readlocs.push(rl);
    }    
  }
  ForceAssertEq( readlocs.size(), Links.size());

  // Create the mapping of the read ids belong to the link
  set<longlong> jump_read_ids_filt;
  set<longlong> long_read_ids_filt;
  map<longlong,longlong> jump_reads_mapping;
  map<longlong,longlong> long_reads_mapping;
  for(size_t i=0;i<Links.size();i++)
  {
    if (Links[i].read_class == 1 )
    {
      jump_read_ids_filt.insert(Links[i].read_id1);
      jump_read_ids_filt.insert(Links[i].read_id2);
    }
    if (Links[i].read_class == 2 )
    {
      long_read_ids_filt.insert(Links[i].read_id1);
      long_read_ids_filt.insert(Links[i].read_id2);
    }
  }

  // create mapping from filt readid to ec readid
  // also get list of all scaffold_reads ids
  set<longlong> sr_ids;
  for(set<longlong>::iterator it=jump_read_ids_filt.begin();
      it!= jump_read_ids_filt.end(); ++it)
  {
    vec<longlong> temp;
    temp.ReadRange(jump_mapping_file,*it, (*it)+1);
    jump_reads_mapping[*it] = temp.back();
    if (temp.back() >=0)
      sr_ids.insert(temp.back());
  }
  for(set<longlong>::iterator it=long_read_ids_filt.begin();
      it!= long_read_ids_filt.end(); ++it)
  {
    vec<longlong> temp;
    temp.ReadRange(long_jump_mapping_file,*it, (*it)+1);
    long_reads_mapping[*it] = temp.back();
    if (temp.back() >=0)
      sr_ids.insert(temp.back() + n_jump_reads_ec);
  }
  if (VERBOSE){
    std::cout << "Links.size()= " << Links.size() << std::endl;
    cout << "    n jump reads"<<jump_reads_mapping.size() << endl;
    cout << "    n long reads"<< long_reads_mapping.size() << endl;
    cout << "    n scaffold reads " << sr_ids.size() << endl;
  }

  // create subset of the qltoutlets and qltoutlets_index array that 
  // contains only the reads relavent
  map<longlong,int> qlt_index_mapping; 
  map<int,alignlet> qlt_mapping; 
  set<int> alignlet_ids;
  for(set<longlong>::iterator it=sr_ids.begin();
      it!= sr_ids.end(); ++it)
  {
    vec<int> temp;
    temp.ReadRange(qlt_index_fn, *it, (*it)+1);
    if (temp.back() >=0 )
    {
      qlt_index_mapping[*it]=temp.back();
      alignlet_ids.insert(temp.back());
    }
  }
  cout << "Alignlet size " << alignlet_ids.size() << endl;
  for(set<int>::iterator it=alignlet_ids.begin();it!=alignlet_ids.end(); ++it)
  {
    vec<alignlet> temp2;
    temp2.ReadRange(qlt_fn, *it, (*it)+1);
    qlt_mapping[*it] = temp2.back();
  }


  //Sort(Links);

  // for all the links
  int n_missed = 0;
  int n_unaligned = 0;
  int n_aligned = 0;
  int n_aligned_filt = 0; // number of pairs aligned in .qltoutlets but not in _filter.qltoutlets
  int n_long = 0;
  for(size_t i=0;i<Links.size();i++)
  {
    read_loc& rl = readlocs[i];
    string dir_str1= Links[i].rc1? "-" : "+";
    string dir_str2= Links[i].rc2? "-" : "+";
    int m1 = Links[i].m1;
    int m2 = Links[i].m2;
    //int s1 = Links[i].s1;
    //int s2 = Links[i].s2;
    //int p1 = to_super_pos[m1];
    //int p2 = to_super_pos[m2];
    int len1 = contigs[m1].size();
    int len2 = contigs[m2].size();
    vec<alignlet> result1;
    vec<alignlet> result2;
    {
      int pos2 = rl.Start( );
      int Pos2 = rl.Stop( );
      int target_id = rl.ContigId( ); 
      Bool fw = rl.Fw( );
      alignlet align(pos2, Pos2, target_id, len1, fw );
      MapAlignletToRef( align, cg_hits, to_hits, result1 );
    }
    // the partner reads
    {
      int pos2 = rl.PartnerStart( );
      int Pos2 = rl.PartnerStop( );
      int target_id = rl.PartnerContigId( ); 
      Bool fw = rl.PartnerFw( );
      alignlet align(pos2, Pos2, target_id, len2, fw );
      MapAlignletToRef( align, cg_hits, to_hits, result2 );
    }

    if (VERBOSE) {
      cout <<"s"<<Links[i].s2<<" "<< Links[i].read_class <<" r/" <<Links[i].read_id1<< dir_str1
	<< "," << Links[i].read_id2 << dir_str2
	<< " c/" << Links[i].m1 <<","<<Links[i].m2<<" ";
      cout << "pos/";
      PrintRefAlign(cout,result1);
      cout << ",";
      PrintRefAlign(cout,result2);
      cout << " ";
      cout << "gap/" << Links[i].sep<<"~"<<Links[i].dev<<" ";
      cout << "sep/" << rl.Sep()<<" ";
    }
    longlong ec_read_id1 = -1, ec_read_id2 = -1;
    if  (Links[i].read_class == 1 )
    {
      ec_read_id1 = jump_reads_mapping[Links[i].read_id1];
      ec_read_id2 = jump_reads_mapping[Links[i].read_id2];
    }
    if  (Links[i].read_class == 2 )
    {
      ec_read_id1 = long_reads_mapping[Links[i].read_id1];
      ec_read_id2 = long_reads_mapping[Links[i].read_id2];
      n_long++;
    }
    if ( ec_read_id1 < 0 || ec_read_id2 <0 ) // not mapped
    {
      if (VERBOSE)
	cout << " r/"<< ec_read_id1<<","<<ec_read_id2 << " missed "<< endl;
      n_missed ++;
      continue;
    }
    // convertion to scaffold read id
    longlong scaffold_read_id1 = ec_read_id1;
    longlong scaffold_read_id2 = ec_read_id2;
    if  (Links[i].read_class == 2 ){
      scaffold_read_id1 += n_jump_reads_ec;
      scaffold_read_id2 += n_jump_reads_ec;
    }
    if (VERBOSE)
      cout << " r/"<<scaffold_read_id1 <<","<<scaffold_read_id2;
    int alignlet_id1 = -1;
    int alignlet_id2 = -1;
    bool align1Filt = False; // are they filtered out by RemoveHCN?
    bool align2Filt = False;
    if ( qlt_index_mapping.find(scaffold_read_id1) != qlt_index_mapping.end())
      alignlet_id1 = qlt_index_mapping[scaffold_read_id1];
    if ( qlt_index_mapping.find(scaffold_read_id2) != qlt_index_mapping.end())
      alignlet_id2 = qlt_index_mapping[scaffold_read_id2];
    if ( alignlet_id1 > 0) {
      vec<int> temp;
      temp.ReadRange(qlt_filt_index_fn, alignlet_id1, alignlet_id1+1);
      if(temp.back() < 0) align1Filt = True;
    }
    if ( alignlet_id2 > 0){
      vec<int> temp;
      temp.ReadRange(qlt_filt_index_fn, alignlet_id1, alignlet_id1+1);
      if(temp.back() < 0) align2Filt = True;
    }
    String strFilt1= align1Filt? "*":"";
    String strFilt2= align1Filt? "*":"";
    if (VERBOSE)
      cout << " a/"<< alignlet_id1 <<strFilt1<<","<<alignlet_id2<<strFilt2<<" " ;
    if ( alignlet_id1 < 0 || alignlet_id2 < 0 )
    {
      // not alignlet found
      if (VERBOSE)
	cout << "unaligned "<< endl;
      n_unaligned ++;
      continue;
    }
    else
    {
      alignlet& align1 = qlt_mapping[alignlet_id1];
      alignlet& align2 = qlt_mapping[alignlet_id2];
      // which scaffold in these two contigs belongs to?
      int c1 = align1.TargetId();
      int c2 = align2.TargetId();
      int s1 = to_super0[c1];
      int s2 = to_super0[c2];
      //cout << "ps/" << s1 << "," << s2 <<" ";
      if (s1 < 0 || s2 < 0) {
	cout << "not in scaffold"<< endl;
	continue;
      } 
      // calculate the distance to end and separation
      //superb &S1 = scaffolds0[s1];
      //superb &S2 = scaffolds0[s2];
      String fw1= align1.Fw1()?"+":"-";
      String fw2= align2.Fw1()?"+":"-";
      int len1 = align1.TargetLength();
      int len2 = align2.TargetLength();
      int p1= to_super_pos0[c1];
      int p2= to_super_pos0[c2];
      int dist_to_end1 = -1;
      int dist_to_end2 = -1;
      ForceAssertEq(contigs0[c1].isize(), len1);
      ForceAssertEq(contigs0[c2].isize(), len2);
      // Convert the alignment position to reference genome

      vec<alignlet> result1;
      vec<alignlet> result2;
      MapAlignletToRef( align1, cg_hits0, to_hits0, result1 );
      MapAlignletToRef( align2, cg_hits0, to_hits0, result2 );
      if (VERBOSE) {
	cout << "aligned " << "c2/"<< c1<< "," <<c2<<" ";
	cout << "pos2/";
	PrintRefAlign(cout,result1);
	cout << ",";
	PrintRefAlign(cout,result2);
	cout << " ";
      }
      const superb& S1_0 = scaffolds0[s1];
      const superb& S2_0 = scaffolds0[s2];

      // use the ShowLinks display convention
      if (align1.Fw1()) // fw1
	dist_to_end1 = S1_0.Len(p1) - align1.Pos2() + S1_0.SubSuperLength( p1 + 1, S1_0.Ntigs( ) - 1 );
      else // rc
	dist_to_end1 = align1.pos2() + S1_0.SubSuperLength( 0, p1 - 1 );
      if (align2.Fw1()) // fw1
	dist_to_end2 = S2_0.Len(p2) - align2.Pos2() + S2_0.SubSuperLength( p2 + 1, S2_0.Ntigs( ) - 1 );
      else // rc
	dist_to_end2 = align2.pos2() + S2_0.SubSuperLength( 0, p2 - 1 );
      int sep = rl.Sep() - dist_to_end1 - dist_to_end2;
      if (VERBOSE) 
      {
	//cout << "pos/";
	//cout << p1 << "(" << len1 <<")"<<fw1<<", ";
	//cout << p2 << "(" << len2<<")"<<fw2<<" ";
	//cout << "sep/"<< sep;
      }
      cout << endl;
      if (align1Filt || align2Filt ) n_aligned_filt ++;
      else n_aligned ++;
    }
  }
  cout << "n_links " << Links.size() << endl;
  std::cout << "  n_missed= " << n_missed << std::endl;
  std::cout << "  n_unaligned= " << n_unaligned << std::endl;
  std::cout << "  n_aligned= " << n_aligned << std::endl;
  std::cout << "  n_aligned_filt= " << n_aligned_filt << std::endl;
  std::cout << "  n_long= " << n_long << std::endl;

}
