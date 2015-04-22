/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2010) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "paths/Alignlet.h"
#include "Alignment.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "paths/Alignlet.h"
#include "paths/ReadLoc.h"
#include "util/ReadTracker.h"
#include "efasta/EfastaTools.h"
#include "feudal/QualNibbleVec.h"
#include "feudal/VirtualMasterVec.h"
#include "paths/KmerPath.h"
#include <getopt.h>


// To print out brief information about certain files for debugging. Also print
// first 10 lines of the records. 

int main( int argc, char *argv[] )
{
  RunTime( );

  //BeginCommandArguments;
  ////CommandArgument_StringSet(FILES);
  ////CommandArgument_String_OrDefault(CMD,"size");
  //CommandArgument_Int_OrDefault(MAX_READS, 10); //maximum 10 reads a time
  //EndCommandArguments;
  longlong START_READS = 0; 
  int MAX_READS = 10; 
  int MAX_LENGTH = 120;
  String format = ""; // force format
  String extra_options = ""; // additional options specific for each format
  int c;
  while( (c = getopt(argc,argv, "m:f:l:s:x:"))!= -1){
    switch(c){
      case 's':
	START_READS = atoi(optarg);
	break;
      case 'm':
	MAX_READS = atoi(optarg);
	break;
      case 'l':
	MAX_LENGTH = atoi(optarg);
	break;
      case 'f':
	format = optarg;
	break;
      case 'x':
	extra_options = optarg;
	break;
      default:
	cout << "Invalid options " << endl;
	exit(1);
    }
  }
  if (MAX_READS < 0) MAX_READS = std::numeric_limits<int>::max();
  if (MAX_LENGTH < 0) MAX_LENGTH = std::numeric_limits<int>::max();
  vec<String> FILES;
  for(int i=optind;i<argc;i++){
    FILES.push_back(argv[i]);
  }

  for(size_t i=0; i<FILES.size();i++)
  {
    String file_name = FILES[i];
    cout << " ------------ " << file_name << " ----------" << endl;
    size_t dot_pos = file_name.find_last_of('.');
    String ext =  file_name.substr(dot_pos, file_name.size() - dot_pos);
    cout << " ext : " << ext << endl;
    if (format != ""){
      cout << "Force format to be " << format << endl;
      ext = format;
    }
    if ( ext == ".qltoutlet") {
      vec<alignlet> aligns;
      longlong n_aligns = BinaryVecNumElements(file_name);
      aligns.ReadRange( file_name, START_READS, min(n_aligns,START_READS+MAX_READS));
      cout << "n_aligns " << n_aligns << endl;
      // find the aligned reads
      int count = 0;
      for(size_t i=0;i<aligns.size();i++)
      {
	alignlet& align = aligns[i];
	String rc= align.Fw1()? "+":"-";
	cout << i+START_READS << " ["<<align.pos2() << "," << align.Pos2() << ") c/" 
	  << align.TargetId() << rc<<" L/" << align.TargetLength() << endl;
      }
    }
    else if ( ext == ".pairs") {
      vec<String> lib_names;
      vec<int> lib_sep;
      vec<int> lib_sd;
      longlong nreads;
      ReadPairsManagerLibInfo(file_name, nreads, lib_names, lib_sep, lib_sd );
      cout << "nreads: "<< nreads << endl;
      for(size_t i=0;i<lib_names.size();i++){
	cout << "lib "<<i<<":"<<lib_names[i]<<" "<<lib_sep[i]
	  <<" "<<lib_sd[i]<< endl;
      }
      PairsManager pairs( file_name);
      size_t n_pairs = pairs.nPairs( );
      cout << "n_pairs " << n_pairs << endl;
      cout << "the pairs: "<< endl;
      for(size_t i=START_READS;i<min(size_t(START_READS+MAX_READS), n_pairs);i++)
      {
	longlong id1 = pairs.ID1(i);
	longlong id2 = pairs.ID2(i);
        int libid = pairs.libraryID(i);
	cout <<i <<" "<< id1 <<" " <<id2 << " libid= " << libid << endl;
      }
    }
    else if ( ext == ".index") {
      vec<int> index;
      longlong n_index= BinaryVecNumElements(file_name);
      int ndisplay = std::min(int(n_index), MAX_READS);
      index.ReadRange( file_name, START_READS, START_READS+ndisplay);
      cout << "n_index " << n_index<< endl;
      // find the aligned reads
      for(size_t i=0;i<index.size();i++)
      {
	cout << i+START_READS << " " << index[i] << endl;
      }
    }
    else if ( ext == ".superb") {
      vec<superb> superbs;
      ReadSuperbs( file_name, superbs );
      int n_superb = superbs.size();
      uint n_superb_tigs = 0;
      for ( int i = 0; i < superbs.isize( ); i++ ) 
	for ( int j = 0; j < superbs[i].Ntigs( ); j++ ) 
	  n_superb_tigs++;
      cout << "n_superb = "<< n_superb<<  endl;
      cout << "n_superb_tigs = "<< n_superb_tigs<<  endl;
      int ndisplay = std::min(int(n_superb), MAX_READS);
      for(int i=START_READS;i<superbs.isize();i++){
	if (i - START_READS >= ndisplay) break;
	cout << "s"<<i << "#"<< superbs[i].Ntigs()<<"("<<superbs[i].FullLength() <<"): ";
	for(int j=0;j<superbs[i].Ntigs();j++){
	  cout << superbs[i].Tig(j)<< " ";
	  if ( extra_options.Contains("g") && j< superbs[i].Ntigs() -1) cout << "("<<superbs[i].Gap(j) <<") ";
	  if ( j > MAX_LENGTH) {
	    cout << " .. ";
	    break;
	  }
	}
	cout << endl;
      }
    }
    else if ( ext == ".efasta") {
      VecEFasta scaffolds, contigs;
      LoadEfastaIntoStrings( file_name, scaffolds );
      vec<superb> superbs;
      SplitEfastaIntoContigs(scaffolds, contigs, superbs );
      vecbasevector tigs( contigs.size( ) );

      int n_superb = superbs.size();
      int n_contig = contigs.size();
      std::cout << "n_superb= " << n_superb << std::endl;
      std::cout << "n_contig= " << n_contig << std::endl;
      uint n_superb_tigs = 0;
      for ( int i = 0; i < superbs.isize( ); i++ ) 
	for ( int j = 0; j < superbs[i].Ntigs( ); j++ ) 
	  n_superb_tigs++;
      cout << "n_superb_tigs = "<< n_superb_tigs<<  endl;
      int ndisplay = std::min(int(n_superb), MAX_READS);
      for(int i=0;i<superbs.isize();i++){
	if (i >= ndisplay) break;
	cout << "s"<<i << "."<< superbs[i].Ntigs()<<"("<<superbs[i].FullLength() <<"): ";
	for(int j=0;j<superbs[i].Ntigs();j++){
	  cout << superbs[i].Tig(j)<< " ";
	  if ( j > MAX_LENGTH) {
	    cout << " .. ";
	    break;
	  }
	}
	cout << endl;
      }
    }
    else if ( ext == ".fasta") {
      vec<fastavector> reads_in;
      LoadFromFastaFile( file_name, reads_in );
      int n_reads = reads_in.size();
      cout << "n_reads "<< n_reads <<  endl;
      for(int i=START_READS;i<reads_in.isize();i++)
      {
	if ( i-START_READS >= MAX_READS) break;
	fastavector & base = reads_in[i];
	int length = base.size();
	if (MAX_LENGTH > 0){
	  base = base.SetToSubOf(base,0,std::min(MAX_LENGTH,length));
	  cout <<i <<" "<< base.ToString() << ".. "<< length <<endl;
	}else{
	  cout <<i <<" "<< base.ToString() <<endl;
	}
      }
    }
    else if ( ext == ".fastb") {
      longlong n_reads = MastervecFileObjectCount(file_name);
      cout << "n_reads "<< n_reads <<  endl;
      VirtualMasterVec<bvec> reads_in(file_name.c_str());

      int ndisplay = std::min(int(n_reads), MAX_READS);
      for(longlong i=START_READS;i< min(START_READS+MAX_READS,n_reads);i++)
      {
	bvec base = reads_in[i];
	int length = base.isize();
	if (MAX_LENGTH > 0){
	  base = base.SetToSubOf(base,0,std::min(MAX_LENGTH,length));
	  cout <<i <<" "<< base.ToString() << ".. "<< length <<endl;
	}else{
	  cout <<i <<" "<< base.ToString() <<endl;
	}
      }
    }
    else if ( ext == ".qualb") {
      longlong n_reads = MastervecFileObjectCount(file_name);
      cout << "n_reads "<< n_reads <<  endl;
      vecqualvector reads_in;
      int ndisplay = std::min(int(n_reads), MAX_READS);
      reads_in.reserve(ndisplay);
      reads_in.ReadRange(file_name, START_READS, START_READS+ndisplay);
      for(size_t i=0;i<reads_in.size();i++)
      {
	qualvector & base = reads_in[i];
	Print(cout, base, "test", MAX_LENGTH);
      }
    }
    else if ( ext == ".readlocs") {
      string scaffold_name = file_name.substr(0,dot_pos);
      read_locs_on_disk locs_file( scaffold_name, "./");
      cout << scaffold_name << endl;
      uint64_t locs_size = locs_file.getLocsSize();
      std::cout << "locs_size= " << locs_size << std::endl;
      vec<read_loc> locs;
      locs_file.LoadContig(START_READS, locs );
      int ndisplay = std::min(int(locs_size), MAX_READS);
      for (int i=0; i< locs.isize(); i++) {
	if (i >= ndisplay) continue;
	read_loc & loc = locs[i];
	String rc1 = loc.Fw()?"+":"-";
	String rc2 = loc.PartnerFw()?"+":"-";
	cout << i << " ["<< loc.Start()<<","<<loc.Stop()<<") "
	  << "c/" << loc.ContigId() << rc1;
	cout << " "<<uint(loc.ReadClass()) << " " << 
	  loc.ReadId() << " " << loc.PartnerReadId() << endl;
	cout << i << " ["<< loc.PartnerStart()<<","<<loc.PartnerStop()<<") "
	  << "c/" << loc.PartnerContigId() << rc2;
	cout << " "<<uint(loc.ReadClass()) << " " << 
	  loc.ReadId() << " " << loc.PartnerReadId() << endl;
      }
    }
    else if ( ext == ".readtrack") {
      string reads_name = file_name.substr(0,dot_pos);
      ReadTracker rt;
      rt.Load(reads_name);
      cout << "Readtrack size: " << rt.size() << endl;
      vec<String> files = rt.source_files;
      cout << "Source size: " << files.size() << endl;
      for(int i=0;i<files.isize();i++){
	cout << "source " << i <<" "<< files[i] << endl;
      }
      for (int i=START_READS; i< min(START_READS+MAX_READS, int64_t(rt.size())); i++) {
	String source = rt.GetReadSource(i);
	uint64_t index = rt.GetReadIndex(i);
	cout << i <<" "<< index << " " << source << endl;
      }
    }
    else if ( ext == ".unipaths") {
      string reads_name = file_name.substr(0,dot_pos);
      vecKmerPath paths;
      paths.ReadAll(file_name);
      cout << "paths.size()= " << paths.size() << endl;
      for(size_t i=START_READS;i<min(size_t(START_READS+MAX_READS), paths.size());i++)
      {
	KmerPath &kpath = paths[i];
	cout << i <<" "<<kpath.TotalLength() << endl;
      }
    }
    else{
      cout << "unsupported ext = " << ext << endl;
    }
  }

}
