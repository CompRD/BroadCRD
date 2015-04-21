///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <sstream>

#include "Basevector.h"
#include "CoreTools.h"
#include "ParseSet.h"
#include "TokenizeString.h"
#include "lookup/LookAlign.h"
#include "lookup/LookAlignSort.h"
#include "lookup/LookAlignFinder.h"
#include "solexa/FourBase.h"
#include "solexa/SolexaTools.h"
#include "FastaFileset.h"
#include "util/BaitMap.h"

void SolexaPredictorParameters::SetFromFile(const String & fname) {
  Ifstream(is, fname);
  String params = Slurp(is);
  SetPredictorParameters(params);
}
  
void SolexaPredictorParameters::SetPredictorParameters(const String & params) {
  istringstream is(params);
  String parameter, value;
  while (true) {
    is >> parameter >> value;
    if (!is) break;
    params_.insert(make_pair(ToLower(parameter), ToLower(value)));
  }
#define KEY_DEFAULT(x,y) if (!IsKey(params_,String(#x)) ) params_[#x] = #y;
  //Set default values as needed.
  KEY_DEFAULT(localrange,2);
  KEY_DEFAULT(posstart ,5 );
  KEY_DEFAULT(minq10i ,1 );
  KEY_DEFAULT(locali ,1 );
  KEY_DEFAULT(o1 ,1 );
  KEY_DEFAULT(o2 ,0 );
#undef KEY_DEFAULT
#if 0
  //Set default values as needed.
  if (!IsKey(params_,"localrange") ) params_["localrange"] = "2";
  if (!IsKey(params_,"posstart") )   params_["posstart"] = "5";
  if (!IsKey(params_,"minq10i") )     params_["minq10i"] = "1";
  if (!IsKey(params_,"locali") )      params_["locali"] = "1";
  if (!IsKey(params_,"o1") )          params_["o1"] = "1";;
  if (!IsKey(params_,"o2") )          params_["o2"] = "0";
#endif
}


String SolexaPredictorParameters::GetPredictorParameters() const {
  ostringstream os;
  for (smap::const_iterator i=params_.begin(); i != params_.end(); ++i) {
    os << i->first << "\t" << i->second  << "\t";
  } 
  return os.str();
}

void LoadSolexaData( const String& HEAD, vecbasevector& bases,
		     vecbasevector& basesrc, VecFourBaseVec& I,
		     vecbasevector& ref,
		     vecbasevector& rref, vec<look_align>& aligns, 
		     vec< vec<int> >& aligns_index, Bool loadIntensities )
{
  LoadSolexaData(HEAD, bases, basesrc, I, loadIntensities);
  // Define file names.

  String ref_file = HEAD + ".reference.fastb";
  String aligns_file = HEAD + ".qltout";
     
  // Load data.
  ref.ReadAll(ref_file );
  rref = ref;
  for ( size_t i = 0; i < ref.size( ); i++ )
    rref[i].ReverseComplement( );
  LoadLookAligns( aligns_file, aligns, aligns_index, bases.size( ) );
}

void LoadSolexaData( const String& HEAD, vecbasevector& bases,
		     vecbasevector& basesrc, VecFourBaseVec& I, Bool loadIntensities ) {
  // Define file names.
  String fastb_file = HEAD + ".fastb";
  String intensity_file = HEAD + ".intensities";
     
  // Load data.

  bases.ReadAll(fastb_file);
  basesrc = bases;
  for ( size_t id = 0; id < bases.size( ); id++ )
    basesrc[id].ReverseComplement( );
  if (loadIntensities) I.ReadAll(intensity_file);
}

void ExpandHead( const String& HEAD, vec<String>& HEADS, int trunc,
     const String& pipeline, Bool paired )
{    vec<String> heads;


  // Below is a special case where the user wants a specific directory 
  // using pipeline and a specific head that may not necessarily follow the 
  // Solexa pipeline conventions.
  if ( !paired && 
       !HEAD.Contains("{") && !HEAD.Contains("}") && !HEAD.Contains(".") ) {
     if ( trunc == 0 )
        HEADS.push_back(pipeline + "/" + HEAD);
     else
        HEADS.push_back(pipeline + "/trunc_" + ToString(trunc) + "/" + HEAD);
  } else {


     // Remove outer brackets.

     String HEADX = HEAD;
     if ( HEAD.Contains( "{", 0 ) && HEAD.Contains( "}", -1 ) )
          HEADX = HEAD.substr( 1, HEAD.isize( ) - 2 );

     // Break at commas that are not inside brackets.  Parse.

     int brackets = 0;
     for ( int i = 0; i < HEADX.isize( ); i++ )
     {    int j;
          for ( j = i; j < HEADX.isize( ); j++ )
          {    if ( HEADX[j] == ',' && brackets == 0 ) break;
               if ( HEADX[j] == '{' ) ++brackets;
               if ( HEADX[j] == '}' ) --brackets;    }
          vec<String> headsi;
          ParseStringSet( HEADX.substr( i, j - i ), headsi );
          heads.append(headsi);
          i = j;    }

     // Expand the heads.

     HEADS.clear( );
     for ( int i = 0; i < heads.isize( ); i++ )
     {    ForceAssert( heads[i].Contains( "." ) );
          String predot = heads[i].Before( "." ), postdot = heads[i].After( "." );
          if ( postdot.Contains( "[" ) )
          {    trunc = postdot.Between( "[", "]" ).Int( );
               postdot = postdot.Before( "[" );    }
          if ( !predot.Contains( "_" ) && !predot.Contains( "/" ) )
	  {    // predot is the number of the flowcell.
	       if (paired) {
	         String runA, runB;
		 // find the paired runs for that flowcell.
		 FindPairedRun( predot, runA, runB, pipeline );
		 if ( trunc == 0 )
		   HEADS.push_back( runA + "/" + heads[i], runB + "/" + heads[i] );
		 else 
		   HEADS.push_back(  runA + "/trunc_" + ToString(trunc) + "/" + heads[i],
				     runB + "/trunc_" + ToString(trunc) + "/" + heads[i]);
               } else {
		 String run, run_date;
		 // find the run for that flowcell.
		 FindRun( predot, run, run_date, pipeline );
		 if ( trunc == 0 ) HEADS.push_back( run + "/" + heads[i] );
		 else 
		   {    HEADS.push_back( 
		        run + "/trunc_" + ToString(trunc) + "/" + heads[i] );    }
	       }
	       continue;    }
          if ( predot.Contains( "_" ) )
          {    String date = predot.Before( "_" );
               String FC = predot.After( "_" );
               // Removed && FC.IsInt( ) below since FCs now can have letters
               if ( date.size( ) == 6 && date.IsInt( ) )
               {    HEADS.push_back( 
                         pipeline + "/" + predot + "/" + FC + "." + postdot );
                    continue;    }    }
          HEADS.push_back( heads[i] );    }    }
     ForceAssert( HEADS.nonempty( ) );    }

void RequireSameRef( const vec<String>& HEADS ) {    
    if ( HEADS.size() < 2 ) return;

    vecbasevector ref;
    String fastb_name = HEADS[0]+".reference.fastb";

    if ( IsRegularFile( fastb_name ) ) {
        ref.ReadRange( fastb_name, 0, MastervecFileObjectCount( fastb_name ) );
    } else {
        FastFetchReads(ref, 0, HEADS[0]+".reference.fasta" );
    }

    ForceAssertGt( ref.size( ), 0u );
    for ( int i = 1; i < HEADS.isize( ); i++ ) {    
        vecbasevector refi;
	fastb_name = HEADS[i]+".reference.fastb";
	if ( IsRegularFile( fastb_name ) ) {
	  refi.ReadRange( fastb_name, 0, MastervecFileObjectCount( fastb_name ) );
	} else {
	  FastFetchReads(refi, 0, HEADS[i]+".reference.fasta" );
	}
	if ( !( ref == refi ) )  {    
	  cout << HEADS[0] + ".reference.fast*"
		   << " is different from "
		   << HEADS[i] + ".reference.fast*.\n";
	  PRINT2( ref.size( ), refi.size( ) );   
	  FatalErr("All experiments must use same reference");
	}
    }
}

void SetRef( const String& REF, String& REFX, const String& HEAD, 
     const vec<String>& HEADS )
{    REFX = REF;
     if ( REF == "" ) RequireSameRef(HEADS);
     if ( REF != "" && HEADS.solo( ) && HEAD.Contains( "." ) )
     {    String predot = HEAD.Before( "." );
          String run, run_date;
          FindRun( predot, run, run_date );
          if ( IsRegularFile( run + "/" + REF ) ) REFX = run + "/" + REF;    }    }

void LoadSolexaData( const String& HEAD, const String& QLTOUT, 
     const String& REF, vecbasevector& bases, VecFourBaseVec& I,
     vecbasevector& ref, vecbasevector& rref, vec<look_align>& aligns, 
     vec< vec<int> >& aligns_index, Bool loadIntensities )
{    vec<String> HEADS;
     ExpandHead( HEAD, HEADS );
     int nreads = 0;
     String REFX;
     SetRef( REF, REFX, HEAD, HEADS );
     FetchReads( ref, 0, REF == "" ? HEADS[0] + ".reference.fasta" : REFX );
     rref = ref;
     ReverseComplement(rref);
     for ( int i = 0; i < HEADS.isize( ); i++ )
     {    bases.ReadAll( HEADS[i] + ".fastb", i > 0 );
          if (loadIntensities) I.ReadAll( HEADS[i] + ".intensities", i > 0 );
          vec<look_align> alignsi;
          vec< vec<int> > alignsi_index;
          LoadLookAligns( HEADS[i] + "." + QLTOUT, alignsi, alignsi_index, 
               bases.size( ) - nreads );
          for ( int j = 0; j < alignsi.isize( ); j++ )
               alignsi[j].query_id += nreads;
          aligns.append(alignsi);
          nreads = bases.size( );    }
     aligns_index.clear_and_resize(nreads);
     for ( int i = 0; i < aligns.isize( ); i++ )
          aligns_index[ aligns[i].query_id ].push_back(i);    }

void LoadSolexaData( const String& HEAD, const String& HEAD_TAIL, const int trunc,
     const String& REF, vecbasevector& bases, VecFourBaseVec& I,
     vecbasevector& ref, Bool loadIntensities )
{    vec<String> HEADS;
     ExpandHead( HEAD, HEADS, trunc );
     for ( int i = 0; i < HEADS.isize( ); i++ )
     {    bases.ReadAll( HEADS[i] + HEAD_TAIL + ".fastb", i > 0 );
          if (loadIntensities) I.ReadAll( HEADS[i] + HEAD_TAIL + ".intensities", i > 0 );    }
     String REFX;
     SetRef( REF, REFX, HEAD, HEADS );
     FetchReads( ref, 0, REF == "" ? HEADS[0] + ".reference.fasta" : REFX );    }

void LoadSolexaData( const String& HEAD, const String& HEAD_TAIL, const int trunc,
     const String& REF, vecbasevector& bases, vecbasevector& ref )
{    vec<String> HEADS;
     ExpandHead( HEAD, HEADS, trunc );
     for ( int i = 0; i < HEADS.isize( ); i++ )
          bases.ReadAll( HEADS[i] + HEAD_TAIL + ".fastb", i > 0 );
     String REFX;
     SetRef( REF, REFX, HEAD, HEADS );
     FetchReads( ref, 0, REF == "" ? HEADS[0] + ".reference.fasta" : REFX );    }

void LoadSolexaData( const String& HEAD, vecbasevector& bases, VecFourBaseVec& I, Bool loadIntensities )
{    vec<String> HEADS;
     ExpandHead( HEAD, HEADS );
     for ( int i = 0; i < HEADS.isize( ); i++ )
     {    bases.ReadAll( HEADS[i] + ".fastb", i > 0 );
          if (loadIntensities) I.ReadAll( HEADS[i] + ".intensities", i > 0 );    }    }

void LoadSolexaData( const String& HEAD, const String& REF, 
     vecbasevector& ref, vecbasevector& rref )
{    vec<String> HEADS;
     ExpandHead( HEAD, HEADS );
     String REFX;
     SetRef( REF, REFX, HEAD, HEADS );
     FetchReads( ref, 0, REF == "" ? HEADS[0] + ".reference.fasta" : REFX );
     rref = ref;
     ReverseComplement(rref);    }

int FilterPileups
(
 const vecbasevector& bases,           // the reads
 const vecbasevector& ref,             // the reference
 const vec<look_align>& aligns,        // alignments of reads to reference
 const vec< vec<int> >& aligns_index,  // index by read of alignments
 const vec<float>& minq10,             // min quality of first 10 bases
 vec<int>& best,                       // index of best alignment (modified)
 const double MIN_RATIO,               // defines filtering
 const Bool verbose
 )
{
  int nreads = bases.size( );
  vec< vec<int> > coverage( ref.size( ) );
  longlong total_bases = 0, total_ref = 0;
  for ( size_t i = 0; i < ref.size( ); i++ ) {    
    coverage[i].resize( ref[i].size( ), 0 );
    total_ref += ref[i].size( );    
  }
  for ( int id = 0; id < nreads; id++ ) {    
    if ( best[id] < 0 ) continue;
    if ( minq10[id] < MIN_RATIO || aligns_index[id].empty( ) ) continue;
    total_bases += bases[id].size( );
    const look_align& la = aligns[ aligns_index[id][ best[id] ] ];
    int id2 = la.target_id;
    int pos2 = la.a.pos2( ), Pos2 = la.a.Pos2( );
    for ( int p = pos2; p < Pos2; p++ )
      coverage[id2][p]++;    
  }
  double mean_coverage = double(total_bases) / double(total_ref);
  int nfiltered = 0;
  for ( int id = 0; id < nreads; id++ ) {    
    if ( best[id] < 0 ) continue;
    if ( minq10[id] < MIN_RATIO || aligns_index[id].empty( ) ) continue;
    const look_align& la = aligns[ aligns_index[id][ best[id] ] ];
    int id2 = la.target_id;
    int pos2 = la.a.pos2( ), Pos2 = la.a.Pos2( );
    for ( int p = pos2; p < Pos2; p++ ) {    
      if ( coverage[id2][p] > 2.0 * mean_coverage ) {    
	++nfiltered;
	if (verbose) 
	  cout << "Filtered: " << bases[id].ToString( ) << "\n";
	best[id] = -1;    
	break;    
      }    
    }    
  }
  return nfiltered;    
}

void GetDecayHeight( const String & intFname, const vecbasevector & reads,
		     NormalDistribution & decay, NormalDistribution & height,
		     const int N) {
  size_t SIZE = MastervecFileObjectCount(intFname);
  AssertEq(SIZE, reads.size());
  if (0 == SIZE) return;
  
  //We use a shuffled vector of indices
  // to pick N non-homopolymer reads at random out of all the possible ones.
  vec<int> shuffled;  
  Shuffle(SIZE, shuffled, SIZE);//use SIZE as the random seed.
  VecFourBaseVec I;

  vec<float> vdecay, vheight;
  int found = 0;
  I.ReadOne(intFname,0);
  Decay<FourBaseVec > dec(I[0].size());
  for (int i=0; found != N && i != shuffled.isize(); ++i) {
    int pos = shuffled[i];
    I.clear();
    I.ReadOne(intFname, pos);
    if (!reads[pos].IsHomopolymer()) {
      float d = dec(I[0]);
      float h = FirstHeight(I[0]);
      if (d > 0 && h > 0) {
	vdecay.push_back(d);
	vheight.push_back(h);
	++found;
      }
    }
  }
  height = SafeMeanStdev(vheight);
  decay = SafeMeanStdev(vdecay);
}

#define ABORT(msg)                   \
{    cout << "\n" << msg << "\n";    \
     cout << "Abort.\n";             \
     exit(1);    }

/// Given a <flowcell>, find the corresponding <run> in the 
/// <Broad pipeline directory>.
/// Return its path as _run_.  Also return the date of the run.

void FindRun( const String& flowcell, String& run, String& run_date, 
     const String& pipeline )
{    run = "";
     vec<String> all_runs = AllFiles(pipeline);
     for ( int i = 0; i < all_runs.isize( ); i++ )
     {    if ( all_runs[i].Contains( "_" + flowcell, -1 ) )
          {    if ( run != "" ) ABORT( "Multiple plausible run directories." );
               run = all_runs[i];    }    }
     if ( run == "" ) ABORT( "Couldn't find run directory for " << flowcell << "." );
     run_date = run.Before( "_" );
     ForceAssertEq( run_date.size( ), 6u );
     run = pipeline + "/" + run;    
}

void FindPairedRun( const String& flowcell, String& runA, String& runB,
     const String& pipeline )
{    runA = "";  runB = ""; 
     vec<String> all_runs = AllFiles(pipeline);
     for ( int i = 0; i < all_runs.isize( ); i++ )
     {    if ( all_runs[i].Contains( "_" + flowcell, -1 ) )
          {    if( runA == "")  runA = all_runs[i];
	       else if (runB == "" )  runB = all_runs[i];
	       else ABORT( "More than two plausible run directories." );    }    }
     if ( runA == "" ) ABORT( "Couldn't find paired run directories for " << flowcell << "." );
     if ( runB == "" ) ABORT( "Couldn't find 2nd paired run directory for " << flowcell << "." );
     if (runA > runB) runA.Swap(runB);
     String runA_date = runA.Before( "_" );
     String runB_date = runB.Before( "_" );
     ForceAssertEq( runA_date.size( ), 6u );
     ForceAssertEq( runB_date.size( ), 6u );
     runA = pipeline + "/" + runA;
     runB = pipeline + "/" + runB;    
}

void LoadVecFromParamsFile( const String & file, const String & columnName, 
			    vec<float> & v) {
  String firstline = LineOfOutput("head -1 " + file);
  vec<String> names;
  Tokenize(firstline, names);
  String lowname = ToLower(columnName);
  int column = -1;
  for (int i=0; i != names.isize(); ++i) {
    names[i].ToLower();
    if (names[i] == lowname) { column = i+1; break; }
  }
  if (-1 == column) FatalErr( "No parameter with name " << columnName
			      << " found in file " << file);
  LoadVecFromParamsFile(file, column, v);
}

void LoadVecFromParamsFile( const String & file, int column, vec<float> & v) {
  temp_file tmp(file+"_XXXXXX");
  SystemSucceed("cut -d' ' -f " + ToString(column) + " " + file
		+ " > " + tmp);
  Ifstream(is,tmp);
  //if first line is not numeric (current format) ignore it.
  //first line might be numeric in the old format.
  char c= is.peek();
  if (!isdigit(c)) is.ignore(numeric_limits<std::streamsize>::max(),'\n');
  v.ReadFromTextStream(is);
}

/**
   For a set of reads in a <.fastb> file, load their alignments to the reference
   from a <.qltout> file.
*/
void LoadReadAligns( const String& qltoutFileName, read_aligns_plus_t& readAligns ) {
  order_lookalign_ErrorRate  sorter;
  cout << "running LookAlignFinder on " << qltoutFileName << endl;
  for (LookAlignFinder finder(qltoutFileName); !finder.empty(); ++finder) {
    int r = finder.QueryId();
    vec<look_align_plus> & aligns = finder.Aligns();
    sort(aligns.begin(), aligns.end(), sorter);
    readAligns[r] = aligns;
  }
}

void LoadBaitMapAsMask(vec<boolvector> *mask, const String &file_name, int padding)
{
    BaitMap baitMap(file_name);
    baitMap.MakeBoolVectorMask(mask, padding);
}
