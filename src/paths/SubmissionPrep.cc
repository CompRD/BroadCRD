///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Prepares assembly for submission to NCBI (or similar). It removes small "
  "contigs and scaffolds, removes duplicate scaffolds, sets a minimum gap size, "
  "generates AGP and tbl files, and optionally removes identified contaminants.";


#include "CoreTools.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "agp/AgpFile.h"
#include "efasta/EfastaTools.h"
#include "VecUtilities.h"
#include "paths/AssemblyCleanupTools.h"
#include "math/HoInterval.h"

/**
 * SubmissionPrep
 *
 * Prepare scaffolds for submission (order, clean, remove small contigs, etc. ).
 *
 */


size_t scaffoldsTotLen( const vec<superb>& scaffolds ){
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    len += scaffolds[is].FullLength();
  return len;
}

size_t scaffoldsRedLen( const vec<superb>& scaffolds ){
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    len += scaffolds[is].ReducedLength();
  return len;
}

size_t scaffoldsNTigs( const vec<superb>& scaffolds ){
  size_t ntigs = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    ntigs += scaffolds[is].Ntigs();
  return ntigs;
}

void scaffoldsPrintStats( const vec<superb>& scaffolds, ostream& out ){
  out << Date() << ": Nscaffolds     = " << scaffolds.size() << endl;
  out << Date() << ": Ncontigs       = " << scaffoldsNTigs( scaffolds ) << endl;
  out << Date() << ": Reduced length = " << scaffoldsRedLen( scaffolds ) << endl;
  out << Date() << ": Total length   = " << scaffoldsTotLen( scaffolds ) << endl;
  out << Date() << ":             ----------" << endl;
  return;
}


String zero_padded_int(const int value, const int length) {
  ostringstream out;
  out.fill('0');
  out.width(length);
  out << value;
  return out.str();
}

String contig_id_str(const int contig_id)
{
  return "contig" + zero_padded_int(contig_id+1, 6);
}

String scaffold_id_str(const int scaf_id)
{
  return "scaffold" + zero_padded_int(scaf_id+1, 5);
}

// restrict efasta given start position and length measured
// using the first option in the curly brackets.
void RestrictEfasta( efasta& source, int & start, int & stop ) {
  if ( stop > source.Length1() || stop < 0 ||
       start < 0 || start >= source.Length1() ){
    PRINT3( start, stop, source.Length1() );
    FatalErr("Inconsistent start and stop positions");
  }
  if ( stop < source.Length1() )
    source.Erase1( stop, source.Length1() );
  if ( start > 0 )
    source.Erase1( 0, start -1 );
  return;
}

// structure to hold contamination removal list (specified in optional tab-delimited input file)
struct contam_spec {
  int contig_id;		// contig id (0-based)
  int contig_length;		// contig length (sanity check to make sure we have right one)
  int contig_begin;		// first base to remove, 0-based coordinates
  int contig_end;		// last base to remove + 1; interval is [begin,end)
  int cut;			// whether to cut the scaffold at this point (non-zero to cut)
};

bool operator<(const contam_spec &a, const contam_spec &b) {
  return (a.contig_id < b.contig_id) ||
    ((a.contig_id == b.contig_id) && (a.contig_begin < b.contig_begin));
}


typedef struct contam_spec contam_spec;
typedef vec<contam_spec> contam_vec;

int load_contamination_list( String in_contam_file, contam_vec& v_contam, const VecEFasta& efastas, Bool& is_zero_based){
  if ( ! IsRegularFile( in_contam_file ) )
      FatalErr("contamination file " + in_contam_file + " not found");
  ifstream in( in_contam_file.c_str() );

  // flags to record the contig ID base in the contamination file
  bool zero_based = true;
  bool one_based = true;

  int contamTotLen = 0;

  uint32_t line_counter = 0;
  String line;
  // reading data in the format:
  cout << Date() << ": reading in contamination file" << endl;
  while ( getline(in,line) ){
    vec<String> tokens;
    contam_spec c;
    Tokenize( line, tokens );
    // Ignore anything beyond 4th column
    ForceAssert( tokens.size() >= 4);
    c.contig_id = tokens[0].Int();
    c.contig_length = tokens[1].Int();
    c.contig_begin = tokens[2].Int();
    c.contig_end   = tokens[3].Int();
    c.cut = 0;

    if (tokens.size() > 4 && tokens[4] == "1") {
      c.cut = tokens[4].Int();
    }
    
    // Determine if the contig ids are 0 or 1 based, and
    // Determine if contig ids and lengths are valid
    bool zero_range = (c.contig_id >= 0) && (c.contig_id < static_cast<int> (efastas.size()));
    bool one_range  = (c.contig_id > 0) && (c.contig_id <= static_cast<int> (efastas.size()));
    bool zero_length = zero_range && (c.contig_length == efastas[c.contig_id].Length1()) ;
    bool one_length  = one_range && (c.contig_length == efastas[c.contig_id - 1].Length1());

    // Contamination file appears to be zero based (or undecided)
    if (zero_based && !zero_length && (!one_based || !one_length) ) {
      cout << Date() << ": An error occurred whilst parsing the contamination file at line: " << line_counter << endl;
      cout << line << endl;
      if (zero_range) {
	cout << Date() << ": Contamination file and assembly contig lengths do not match:" << endl;
	cout << Date() << ":   Contamination contig ID= " << c.contig_id << ", length= " << c.contig_length << endl; 
	cout << Date() << ":   Assembly      contig ID= " << c.contig_id << ", length= " << efastas[c.contig_id].Length1() << endl; 
      } else if ( one_based) 
	cout << Date() << ": The contig ID may be out of range [0," << efastas.size() - 1 << "] "
	     << "if the contamination file is zero based." << endl;
      else 
	cout << Date() << ": The contig ID is out of range [0," << efastas.size() - 1 << "]" << endl;
    }
    
    // Contamination file appears to be one based (or undecided)
    if (one_based && !one_length && (!zero_based || !zero_length) ) {
      if (!zero_based)  {
	cout << Date() << ": An error occurred whilst parsing the contamination file at line: " << line_counter << endl;
	cout << line << endl;
      }
      if (one_range) {
	cout << Date() << ": Contamination file and assembly contig lengths do not match:" << endl;
	cout << Date() << ":   Contamination contig ID= " << c.contig_id     << ", length= " << c.contig_length << endl; 
	cout << Date() << ":   Assembly      contig ID= " << c.contig_id - 1 << ", length= " << efastas[c.contig_id -1].Length1() << endl; 
      } else if (zero_based)
	cout << Date() << ": The contig ID may be out of range [1," << efastas.size() << "] "
	     << "if the contamination file is one based." << endl;
      else
	cout << Date() << ": The contig ID is out of range [1," << efastas.size() << "]" << endl;

      cout << Date() << ": Note: ALLPATHS-LG starts numbering contigs IDs at zero, but NCBI starts at one." << endl;
      if (!zero_based)
	cout << Date() << ":       SubmissionPrep has determined this contamination file is one based." << endl;
      else if (zero_based) { 
	cout << Date() << ":       SubmissionPrep was unable to determine if the contamination file is zero or one based, " << endl;
	cout << Date() << ":       so error messages for both cases are given above." << endl;
      }
    }

    // Update contamination file base flag
    zero_based &= zero_length;
    one_based &= one_length;
    
    if (!zero_based && !one_based)
      FatalErr("Unable to parse contamination file - stopping");
    
  
    // Determine if the section to trim is sensible
    int contig_size = efastas[c.contig_id - (zero_based ? 0 : 1)].Length1();
    bool begin_ok = (c.contig_begin >= 0) && (c.contig_begin < contig_size);
    bool end_ok = (c.contig_end >= 0 ) && (c.contig_end <= contig_size);
    bool range_ok = (c.contig_end >= c.contig_begin);
    if (false == (begin_ok && end_ok && range_ok))  {
      cout << Date() << ": An error occurred whilst parsing the contamination file at line: " << line_counter << endl;
      cout << line << endl;
      if (!begin_ok)
	cout << Date() << ": The contamination begin position " << c.contig_begin 
	     << " is out of range [0, " << contig_size - 1 << "]" << endl;
      if (!end_ok)
	cout << Date() << ": The contamination end position " << c.contig_end 
	     << " is out of range [0, " << contig_size << "]" << endl;
      if (!range_ok)
	cout << Date() << ": The contamination end position "  << c.contig_end 
	     << " is less than the begin position " << c.contig_begin << endl;
      FatalErr("Unable to parse contamination file - stopping");
    }

    // Line parsed successfully, add to list.
    contamTotLen += c.contig_end - c.contig_begin;
    v_contam.push_back(c);
    line_counter++;
  }

  in.close();

  // Adjust contig ID base if neccessary
  if (zero_based){
    is_zero_based = True;
    cout << "\nINFO: Contamination list contig ID values are zero based\n" << endl;
  }
  else if (one_based) {
    is_zero_based = False;
    cout << "\nINFO: Contamination list contig ID values are one based\n" << endl;    
    for (size_t index = 0; index < v_contam.size(); index++) 
      --v_contam[index].contig_id;
  }

  // Order by contigID and start base
  Sort(v_contam);

  cout << Date() << ": Total contamination sequence length = " << contamTotLen << endl;
  return contamTotLen;
}

int remove_contamination_list( const contam_vec& v_contam,
			       VecEFasta& efastas, vec<superb>& scaffolds,
			       vec<String>& tigMap, const String save_contam_file ){


  vec<superb> new_tscaffolds( efastas.size() ); 
  vec<size_t> cuts_ids;
  vec<Bool> modif_contigs( efastas.size(), False );
  vec< pair<int,int> > beg_gaps( efastas.size() );
  vec< pair<int,int> > end_gaps( efastas.size() );
  vec<fastavector> contamination;
  vec<String> contamination_ids;

  vec<contam_vec> tig_contam(efastas.size());
  for ( size_t i = 0; i < v_contam.size(); i++ ) 
    tig_contam[v_contam[i].contig_id].push_back(v_contam[i]);
    

  int contamCheckLen = 0;
  for ( size_t tid = 0; tid < tig_contam.size(); tid++ ){ 
    if ( tig_contam.at(tid).size() == 0 ) 
      continue;
    else
      modif_contigs[tid] = True;
       
    basevector tbases;
    efastas[tid].FlattenTo( tbases );

    int specContamLen = 0, resContamLen = 0;
    // map of cut points
    vec<int> cut_points;
    for ( int i = tig_contam[tid].isize() -1; i >= 0; i-- ){
      int clen = tig_contam[tid][i].contig_length;
      int cbeg = tig_contam[tid][i].contig_begin;
      int cend = tig_contam[tid][i].contig_end;

      ForceAssertEq( (int)efastas[tid].Length1(), clen );
      specContamLen +=  cend - cbeg;
      if (tig_contam[tid][i].cut)
	cut_points.push_back( cend );

      basevector  contam_bv(tbases, cbeg, specContamLen);
      fastavector contam_fv(contam_bv);
      contamination.push_back(contam_fv);
      contamination_ids.push_back("contam_contig" + zero_padded_int(tid, 6) + "_" + ToString(cbeg) + "_" + ToString(cend));
    }
    cout << Date() << ": cuts_points: "; cut_points.Print(cout); cout << endl;

    // make sure there are no overlaps
    vec<ho_interval> intervs;
    for ( size_t i = 0; i < tig_contam[tid].size(); i++ )
      intervs.push_back( ho_interval( tig_contam[tid][i].contig_begin, tig_contam[tid][i].contig_end ) );
    ForceAssertEq( Sum(intervs), TotalCovered(intervs) );
    contamCheckLen += Sum(intervs);

    cout << Date() << ": introducing cut_points in contig " << tid << endl;
    vec<ho_interval> remains;
    Uncovered( efastas[tid].Length1(), intervs, remains ); 
 
    for ( size_t cpi = 0; cpi < cut_points.size(); cpi++ ){
      int cp = cut_points[cpi];
      size_t remains_orig_size = remains.size();
      for ( size_t rpi = 0; rpi < remains_orig_size; rpi++ ){	
	int start = remains[rpi].Start();
	int stop  = remains[rpi].Stop();
	if ( start <= cp && stop > cp ){
	  ho_interval lho( start, cp );
	  remains[rpi] = lho;
	  ho_interval rho( cp, stop );
	  remains.push_back( rho );
	}
      }
    }
     
    Sort( remains );
    //PRINT( remains.size() );
    //for ( int l = 0; l < remains.isize(); l++ )
    //PRINT2( remains[l].Start(), remains[l].Stop() );

    cout << Date() << ": updating sequences" << endl;
    
    efasta efasta_orig    = efastas[tid];

    new_tscaffolds[tid].SetNtigs( remains.size() );
    if ( remains.size() == 0 ){
      efastas[tid].resize(0);
      new_tscaffolds[tid].SetNtigs( 1 );
      new_tscaffolds[tid].SetTig(0, tid);
      new_tscaffolds[tid].SetLen(0, 0);
    } else {

      for ( int ri = 0; ri < remains.isize(); ri++ ){
	ho_interval rem = remains[ri];
	efasta locEfasta;
	
	efasta_orig.Extract1( rem.Start(), rem.Stop(), locEfasta );

	new_tscaffolds[tid].SetLen(ri, locEfasta.Length1() );
	
	bool cut = BinPosition( cut_points, rem.Start() ) == -1 ? false : true;

	
	// If cut is set, then we want to cut the scaffold *before* this contig.
	// We'll push the contig we want to cut before in a list of cuts.
	if (cut) cout << Date() << ": CUT AT " << rem.Start() << endl;
	
	//cout << Date() << ": updating efasta and base vecs" << endl;
	if ( ri == 0 ){
	  new_tscaffolds[tid].SetTig(ri, tid);
	  efastas[tid] = locEfasta;
	  if (cut) cuts_ids.push_back(tid);
	} else {
	  int new_tid = efastas.size();
	  new_tscaffolds[tid].SetTig(ri, new_tid );
	  efastas.push_back( locEfasta );
	  tigMap.push_back( ToString(tid) );
	  if (cut) cuts_ids.push_back(new_tid);
	}
	if ( ri != remains.isize() -1 ){
	  new_tscaffolds[tid].SetGap(ri, remains[ri+1].Start() - remains[ri].Stop() );
	  new_tscaffolds[tid].SetDev(ri, 1);
	}
	if (cut) cout << Date() << ": CUT BEFORE TIG " << cuts_ids.back() << endl;
      }
    }
  }
  cout << Date() << ": cuts_ids: "; cuts_ids.Print(cout); cout << endl;

  // update scaffolds
  size_t endGapSum = 0;
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    superb & s = scaffolds[si];
    for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ) {
      size_t tid = s.Tig(tpos);

      if ( tid < modif_contigs.size() && modif_contigs[tid] ){

	size_t lenBefore = s.FullLength();
	size_t lenAfter = 0;
	if ( tpos == 0 ){
	  lenAfter += beg_gaps[tid].first;
	  if ( s.Ntigs() > 1 && new_tscaffolds[tid].FullLength() == 0
	       && new_tscaffolds[tid].Ntigs() == 0 ){
	    lenAfter += end_gaps[tid].first;
	    lenAfter += s.Gap(tpos);
	  }
	}
	if ( tpos == s.Ntigs() -1 ) {
	  lenAfter += end_gaps[tid].first;
	  if ( s.Ntigs() > 1 && new_tscaffolds[tid].FullLength() == 0
	       && new_tscaffolds[tid].Ntigs() == 0 ){
	    lenAfter += beg_gaps[tid].first;
	    lenAfter += s.Gap(tpos -1 );
	  }
	}
	endGapSum += lenAfter;

	s.ReplaceTigBySuper( tpos, new_tscaffolds[tid], 
			     beg_gaps[tid].first, beg_gaps[tid].second,
			     end_gaps[tid].first, end_gaps[tid].second);

	modif_contigs[tid] = False;

	lenAfter += s.FullLength();

	tpos += new_tscaffolds[tid].Ntigs() -1;
      }
    }
  }  

  if (cuts_ids.size() > 0) {
    cout << Date() << ": cutting scaffolds" << endl;
    for (size_t c = 0; c < cuts_ids.size(); ++c) {
      for ( size_t si = 0; si < scaffolds.size(); si++ ) {
	superb & s = scaffolds[si];
	for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ) {
	  size_t tid = s.Tig(tpos);
	  bool cut = false;
	  if (cuts_ids[c] == tid) cut = true;
	  if (cut) {
	    // Break scaffold after this contig
	    cout << Date() << ": FOUND CUT TIG " << tid << " IN SCAF " << si << " tpos=" << tpos << " ntigs=" << s.Ntigs() << endl;
	    superb new_s;
	    new_s.SetNtigs(s.Ntigs() - tpos);
	    for (int t = tpos; t < s.Ntigs(); ++t) {
	      int i = t - tpos;
	      new_s.SetTig(i, s.Tig(t));
	      new_s.SetLen(i, s.Len(t));
	      if (t < s.Ntigs() - 1) {
		new_s.SetGap(i, s.Gap(t));
		new_s.SetDev(i, s.Dev(t));
	      }
	    }
	    s.SetNtigs(tpos);
	    scaffolds.push_back(new_s);
	    cuts_ids[c] = -1;
	  }
	}
      }
    }
  }
      
  if (save_contam_file != "") {
    cout << Date() << ": saving contamination segments to " << save_contam_file << endl;
    Ofstream(contam_stream, save_contam_file);
    for (u_int i = 0; i < contamination.size(); ++i)
      contamination[i].Print(contam_stream, contamination_ids[i]);
    contam_stream.close();
  }
  return contamCheckLen;
}

void
translate_rings_file(const String input_filename, const String output_filename, const vec<String> & scaffMap )
{
  Ifstream (in, input_filename);
  Ofstream (out, output_filename);
  String line;
  out << "#scaffold\tgap_size\tgap_dev\tn_links" << endl;
  while (getline(in,line)){
    vec<String> tokens;
    Tokenize( line, tokens );
    if (tokens.size() == 4 && tokens[0][0] >= '0' && tokens[0][0] <= '9') {
      int scaffold    = tokens[0].Int();
      int gap_size    = tokens[1].Int();
      int gap_dev     = tokens[2].Int();
      int n_links     = tokens[3].Int();

      for (size_t s = 0; s < scaffMap.size(); ++s) {
	if (scaffMap[s].Int() == scaffold) {
	  scaffold = s;
	  break;
	}
      }

      out << scaffold_id_str(scaffold);
      out << "\t" << gap_size << "\t" << gap_dev << "\t" << n_links << endl;      
    }
  }
  in.close();
  out.close();
}


void
remove_contamination(String in_contam_file,
		     Assembly& assembly,
		     int min_contig_solo,
		     int min_contig_in,
		     String out_save_contam_file)
{
  
  vec<String> scaffMap  = assembly.scaffMap();
  vec<String> tigMap    = assembly.tigMap();
  vec<superb> scaffolds = assembly.scaffolds();
  VecEFasta efastas   = assembly.efastas();
  digraphE<sepdev> SG   = assembly.SG();

  cout << Date() << ": loading contaminant information" << endl;
  contam_vec v_contam;
  Bool list_zero_based = True;
  int contamTotLen = 
    load_contamination_list(in_contam_file, v_contam, efastas, list_zero_based);

  if ( ! list_zero_based )
    cout << "\nNOTE:  contamination list is 'one' based (first contig number is 1) however the following logging output is zero based (first contig number is 0).\n" << endl;

  int contamCheckLen = 
    remove_contamination_list( v_contam, efastas, scaffolds, tigMap, out_save_contam_file);
  if ( contamTotLen != contamCheckLen ){
    cout << Date() << ":  There seems to be a problem with contamination specification, are there contaminant overlaps?" << endl;
    ForceAssertEq( contamTotLen, contamCheckLen );
  }
  
  Assembly assembly_new( scaffolds, efastas, &scaffMap, &tigMap );

  assembly_new.remove_small_contigs(min_contig_solo, min_contig_in);
  assembly_new.remove_small_scaffolds(min_contig_solo);
  assembly_new.remove_unused_contigs();
  assembly_new.check_integrity();

  assembly = assembly_new;
}


int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc( HEAD_IN,
    "Input assembly files <HEAD_IN>.*");
  CommandArgument_String_Doc( HEAD_OUT,
    "Output assembly files <HEAD_OUT>.*" );
  CommandArgument_Int_OrDefault_Doc( MIN_CONTIG_SIZE_SOLO, 1000, 
    "Remove solo contigs that are smaller than MIN_CONTIG_SIZE_SOLO");
  CommandArgument_Int_OrDefault_Doc( MIN_CONTIG_SIZE_IN, 200, 
    "Remove contigs within scaffolds that are smaller than MIN_CONTIG_SIZE_IN");
  CommandArgument_Bool_OrDefault_Doc( REORDER, False,
    "Reorder scaffolds, largest to smallest.");
  CommandArgument_Bool_OrDefault_Doc( DEDUP, True,
    "Remove duplicate scaffolds.");
  CommandArgument_Int_OrDefault_Doc( MIN_GAP, 10, 
    "Gaps smaller than this will be increased to this value (per NCBI)");
  CommandArgument_String_OrDefault_Doc( CONTAM, "", 
    "Full path name of contamination file with: contig_id contig_lenght begin end");
  CommandArgument_Bool_OrDefault_Doc( REMOVE_CONTAM_LAST, False,
    "Remove contamination information at the end of the process. This option assumes that contamination coordinates corresopond to the modified assembly created by the first part of the code.");
  CommandArgument_Bool_OrDefault_Doc( SAVE_CONTAM, True,
    "Save fasta of removed contamination segments.");
  CommandArgument_UnsignedInt_OrDefault_Doc( FASTA_BATCH_SIZE, 10000,
    "Output contig fastas won't contain more than this many elements.");
  CommandArgument_String_OrDefault_Doc( ISOLATE_CONTIGS, "",
    "Detach specified contigs from their scaffolds, but don't remove them.");
  CommandArgument_String_OrDefault_Doc( EXTRA_CONTIGS, "",
    "Fasta file with extra contigs to be merged into the assembly as singleton contigs.");
  CommandArgument_Bool_OrDefault_Doc( AGP2, True,
    "Generate AGP file using version 2.0 spec (otherwise 1.1)");
  EndCommandArguments;

  // File names.
  String in_efasta_file = HEAD_IN + ".contigs.efasta";
  String in_superb_file = HEAD_IN + ".superb";
  String in_contam_file = CONTAM;

  String out_efasta_file        = HEAD_OUT + ".contigs.efasta";
  String out_fasta_file         = HEAD_OUT + ".contigs.fasta";
  String out_fastamb_file       = HEAD_OUT + ".contigs.fastamb";
  String out_tbl_file           = HEAD_OUT + ".contigs.tbl";
  String out_superb_file        = HEAD_OUT + ".superb";
  String out_fastb_file         = HEAD_OUT + ".fastb";
  String out_contig_map_file    = HEAD_OUT + ".contigs.mapping";
  String out_scaffold_map_file  = HEAD_OUT + ".superb.mapping";
  String out_save_contam_file   = SAVE_CONTAM ? HEAD_OUT + ".contam.fasta" : "";

  vec<int> isolate_contigs;
  ParseIntSet( ISOLATE_CONTIGS, isolate_contigs, false );

  // Lodading scaffolds
  cout << Date() << ": loading superb file" << endl;
  vec<superb> scaffolds;
  ReadSuperbs( in_superb_file, scaffolds );
  cout << Date() <<  ": Input scaffolds:" << endl;
  scaffoldsPrintStats( scaffolds, cout ); 

  // reading contig information
  cout << Date( ) << ": loading contigs efasta file" << endl;
  if ( ! IsRegularFile( in_efasta_file) )
    FatalErr("input file " + in_efasta_file + " not found");
  VecEFasta efastas;
  LoadEfastaIntoStrings(in_efasta_file, efastas);
  cout << Date() <<  ": Input contigs " << efastas.size() << endl;


  if (EXTRA_CONTIGS != "") {
    cout <<  Date( ) << ": loading extra contigs" << endl;
    vec<fastavector> extra_contigs;
    LoadFromFastaFile( EXTRA_CONTIGS, extra_contigs );
    
    // Add the extra contigs in as supers
    for (size_t i = 0; i < extra_contigs.size(); ++i) {
      size_t tig = efastas.size();
      efasta contig(extra_contigs[i]);
      superb super;
      super.PlaceFirstTig(tig, contig.Length1());
      scaffolds.push_back(super);
      efastas.push_back(contig);
    }
    cout << Date( ) << ": " << extra_contigs.size() << " extra contigs appended as singleton scaffolds" << endl;
    cout << Date() << ": number of extra contig read = " << extra_contigs.size() << endl;
  }

  // make sure that there are no N's at the ends
  for ( size_t id = 0; id < efastas.size(); id++ ){
    basevector bases;
    efastas[id].FlattenTo( bases );
    String sbases = bases.ToString();
    ForceAssert( sbases.at(0) != 'N' );
    ForceAssert( sbases.at(0) != 'n' );

    ForceAssert( sbases.back() != 'N' );
    ForceAssert( sbases.back() != 'n' );
  }

  if (isolate_contigs.size() > 0) {
    int nscaffolds_orig_size = scaffolds.size();
    for (int i = 0; i < isolate_contigs.isize(); ++i) {
      int isolate = isolate_contigs[i];
      cout << Date( ) << " isolating contig " << isolate << " into it's own scaffold" << endl;
      for (int s = 0; s < nscaffolds_orig_size; ++s) {
	for (int t = 0; t < scaffolds[s].Ntigs(); ++t) {
	  int tig = scaffolds[s].Tig(t);
	  if (tig == isolate) {
	    int len = scaffolds[s].Len(t);
	    PRINT4(tig, s, t, len);
	    superb iso_super;
	    iso_super.PlaceFirstTig(tig, len);
	    scaffolds[s].RemoveTigByPos(t);
	    scaffolds.push_back(iso_super);
	    break;
	  }
	}
      }
    }
  }

  size_t origScaffoldsTotLen = scaffoldsTotLen( scaffolds );
  size_t origScaffoldsRedLen = scaffoldsRedLen( scaffolds );

  Assembly assembly( scaffolds, efastas );

  // check initial basic integrity
  assembly.check_integrity();


  // remove contamination
  vec<fastavector> contamination;
  if ( in_contam_file.nonempty() && !REMOVE_CONTAM_LAST) {
    remove_contamination(in_contam_file, assembly, MIN_CONTIG_SIZE_SOLO, MIN_CONTIG_SIZE_IN,
			 out_save_contam_file);
  }

  assembly.remove_small_contigs( MIN_CONTIG_SIZE_SOLO, MIN_CONTIG_SIZE_IN );
  assembly.remove_small_scaffolds( MIN_CONTIG_SIZE_SOLO );
  assembly.remove_unused_contigs();
  assembly.check_integrity();


  if (DEDUP) assembly.dedup_exact();
  if (REORDER) assembly.reorder();

  // check integrity
  assembly.check_integrity();

  // remove contamination
  if ( in_contam_file.nonempty() && REMOVE_CONTAM_LAST) {
    remove_contamination(in_contam_file, assembly, MIN_CONTIG_SIZE_SOLO, MIN_CONTIG_SIZE_IN,
			 out_save_contam_file);
  }

  // reset small or negative gaps
  assembly.set_min_gaps(MIN_GAP);
  

  // writing output
  cout << Date() << ": writing superb" << endl;
  WriteSuperbs( out_superb_file, assembly.scaffolds() );

  cout << Date() <<  ": Output scaffolds:" << endl;
  scaffoldsPrintStats( assembly.scaffolds(), cout ); 

  cout << Date() << ": number of contigs written = " << assembly.efastas().size() << endl;


  cout << Date() << ": writing assembly files" << endl;
  Ofstream( efout, out_efasta_file );
  Ofstream( fout, out_fasta_file );
  Ofstream( famb, out_fastamb_file );
  Ofstream( tblout, out_tbl_file );

  int batch = 0;
  String batch_suffix = (assembly.efastas().size()>FASTA_BATCH_SIZE) ? ("_" + zero_padded_int(batch, 2)) : "";
  Ofstream(fsa, HEAD_OUT + batch_suffix + ".fsa");
  Ofstream(tbl, HEAD_OUT + batch_suffix + ".tbl");

  vec<fastavector> fastas(assembly.efastas().size());
  int snp_count = 0,indel_count = 0, ambiguous_bases = 0;

  for ( size_t id = 0; id < assembly.efastas().size(); id++ ) {
    String name = contig_id_str(id);

    // output efasta record
    assembly.efastas()[id].Print(efout, name);

    // chunk into files of no more than 10,000 contigs for NCBI
    if ((id+1) % FASTA_BATCH_SIZE == 0) {
      cout << "." << flush;
      fsa.close();
      ++batch;
      batch_suffix = "_" + zero_padded_int(batch, 2);
      OpenOfstream( fsa, HEAD_OUT + batch_suffix + ".fsa");
      tbl.close();
      OpenOfstream( tbl, HEAD_OUT + batch_suffix + ".tbl");
    }

    // code fragment from EfastaToFasta. --bruce
    // flatten to fsa (fasta) without ambiguities
    fastavector v;
    vec<Ambiguity> va;
    assembly.efastas()[id].FlattenTo( v, va, False );
    fastas[id] = v;

    v.Print(fsa, name);
    v.Print(fout, name);

    if (va.size()) {
      ostringstream feature;

      feature << ">Feature " << name << "\n";

      // add annotation of alternate paths to tbl file
      for (size_t i = 0; i < va.size(); i++)
	feature << va[i].to_annotation();
      tbl << feature.str();
      tblout << feature.str();
    }

    int snps = 0, indels = 0;
    int amb =  assembly.efastas()[id].AmbCount(snps, indels);
    ambiguous_bases += amb;
    snp_count += snps;
    indel_count += indels;

    // flatten to fasta with ambiguities for possible fixup
    assembly.efastas()[id].FlattenTo( v, va, True );
    v.Print(famb, name);

  }
  efout.close();
  tbl.close();
  tblout.close();
  fsa.close();
  fout.close();

  if (assembly.efastas().size() >= FASTA_BATCH_SIZE)
    cout << endl;

  // Generate agp file...code fragment from agp/SuperbsToAgp. --bruce
  cout << Date() << ": writing agp file" << endl;
  Ofstream(agpfile, HEAD_OUT + ".agp");
  for (u_int ii=0; ii<assembly.scaffolds().size( ); ii++) {
    agp_chromosome agp( zero_padded_int(ii+1, 5 )) ;
    
    // Loop over the contigs.
    for (int jj=0; jj<assembly.scaffolds()[ii].Ntigs( ); jj++) {

      // Add contig.
      String cg_name = contig_id_str(assembly.scaffolds()[ii].Tig(jj));

      int cg_len = assembly.scaffolds()[ii].Len( jj );
      ForceAssertEq((u_int)cg_len, fastas[assembly.scaffolds()[ii].Tig( jj )].size());
      int cg_start = 0;
      int cg_stop = cg_len - 1;
      bool cg_RC = false;
      agp_contig new_contig( cg_name, cg_len, cg_start, cg_stop, cg_RC );
      new_contig.SetType( agp_contig::wgs_contig );
      agp.AddContig ( new_contig );

      // Add gap.
      if ( jj < assembly.scaffolds()[ii].Ntigs( ) - 1 ) {
	String gap_type = AGP2 ? "scaffold" : "fragment";
	bool is_bridged = true;
	int gap_len = assembly.scaffolds()[ii].Gap( jj );
	String linkage = AGP2 ? "paired-ends" : "";
	agp_gap new_gap( gap_type, gap_len, is_bridged, linkage );
	agp.AddGap( new_gap );
      }
    }
    
    // Print info for super.
    String base_name = "scaffold";
    agp.Print( agpfile, &base_name );
  }

  //vecbasevector obases( assembly.efastas().size() );
  //for ( size_t id = 0; id < assembly.efastas().size(); id++ )
  //  assembly.efastas()[id].FlattenTo( obases[id] );
  
  //obases.WriteAll( out_fastb_file );


  cout << Date() << ": writing gapped assembly efasta" << endl;
  WriteScaffoldedEFasta( HEAD_OUT + ".assembly.efasta", assembly.efastas(), assembly.scaffolds(), True);

  cout << Date() << ": writing gapped assembly fasta" << endl;
  vec<Bool> rc(fastas.size(), False);
  WriteScaffoldedFasta( HEAD_OUT + ".assembly.fasta", fastas, assembly.scaffolds(), rc, MIN_GAP, 'N', True);

  Ofstream( cmout, out_contig_map_file );
  for ( size_t id = 0; id < assembly.tigMap().size(); id++ )
    cmout << contig_id_str(id) + " from " + ToString( assembly.tigMap()[id] ) << "\n";

  Ofstream( smout, out_scaffold_map_file );
  for ( size_t is = 0; is < assembly.scaffMap().size(); is++ )
    smout << scaffold_id_str(is) + " from " + ToString( assembly.scaffMap()[is] ) << "\n";

  {
    size_t newScaffoldsTotLen = 0, newScaffoldsRedLen = 0, newContigTotLen = 0;
    vec<int> scaffolds_full(assembly.scaffolds().size());
    vec<int> scaffolds_reduced(assembly.scaffolds().size());
    for ( size_t is = 0; is < assembly.scaffolds().size(); is++ ){
      int tlen = assembly.scaffolds()[is].FullLength();
      int rlen = assembly.scaffolds()[is].ReducedLength();
      scaffolds_full[is] = tlen;
      newScaffoldsTotLen += tlen;
      scaffolds_reduced[is] = rlen;
      newScaffoldsRedLen += rlen;
    }
    vec<int> contig_sizes(fastas.size());
    for ( size_t ic = 0; ic < fastas.size(); ic++ ){
      int len = fastas[ic].size();
      newContigTotLen += len;
      contig_sizes[ic] = len;
    }

    String rings_file_in = HEAD_IN + ".rings";
    String rings_file_out = HEAD_OUT + ".rings";
    if (IsRegularFile(rings_file_in)) {
      cout << Date() << ": Translating rings file" << endl;
      translate_rings_file(rings_file_in, rings_file_out, assembly.scaffMap());
    }


    cout << Date() << ": Assembly statistics:" << endl;
    cout << Date() << ": Number of scaffolds: " << ToStringAddCommas(assembly.scaffolds().size()) << endl;
    cout << Date() << ": Total scaffold length including gaps: " << ToStringAddCommas(newScaffoldsTotLen)
	 << " (N50 = " << ToStringAddCommas(N50(scaffolds_full)) << ")" << endl;
    cout << Date() << ": Total scaffold length excluding gaps: " << ToStringAddCommas(newScaffoldsRedLen)
	 << " (N50 = " << ToStringAddCommas(N50(scaffolds_reduced)) << ")" << endl;
    cout << Date() << ": Number of contigs: " << ToStringAddCommas(fastas.size())
	 << " (N50 = " << ToStringAddCommas(N50(contig_sizes)) << ")" << endl;

    if (ambiguous_bases > 0)
      cout << Date() << ": Ambiguous bases: " << ToStringAddCommas(ambiguous_bases)
	   << " (rate = 1/" << ToStringAddCommas(newScaffoldsRedLen / ambiguous_bases) << ")" << endl;
    if (snp_count > 0)
      cout << Date() << ": SNP count: " << ToStringAddCommas(snp_count)
	   << " (rate = 1/" << ToStringAddCommas(newScaffoldsRedLen / snp_count) << ")" << endl;
    if (indel_count > 0)
      cout << Date() << ": Indel count: " << ToStringAddCommas(indel_count)
	   << " (rate = 1/" << ToStringAddCommas(newScaffoldsRedLen / indel_count) << ")" << endl;
  }

  cout << Date() << ": Done!" << endl;
}

