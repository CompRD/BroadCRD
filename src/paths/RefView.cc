/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: RefView

   Show a window of the reference -- the reads that align there, the unipaths
   that align there, etc.   Help us see why a window of the reference was
   missing from the assembly or assembled incorrectly.

   The reads and/or the unipaths must have already been aligned to the assembly.

   The window information is written out to an HTML file.

   @file
*/

#include <sstream>
#include <algorithm>
#include "MainTools.h"
#include "ParseSet.h"
#include "ReadLocation.h"
#include "Qualvector.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"
#include "paths/simulation/GenomePlacement.h"
#include "solexa/FourBase.h"
#include "system/HTMLUtils.h"
#include "paths/DisplayUtils.h"
#include "lookup/LookAlign.h"
#include "graph/DotUtils.h"
#include "feudal/BinaryStream.h"

/**
   Function: IntervalsOverlap

   Tests whether two closed discrete intervals [x_lower,x_upper] and
   [y_lower,y_upper] share at least one point.
*/
Bool IntervalsOverlap( int x_lower, int x_upper, int y_lower, int y_upper ) {
  return x_lower <= y_lower && y_lower <= x_upper ||
         y_lower <= x_lower && x_lower <= y_upper;
}



#if 0

class HTMLDoc {
  String id;
  String title;
  vec< HTMLDoc > subs;

  Bool IsLeaf() const { return subs.empty(); }

};

/**
   Class: HTMLFrames

   Represents a set of frames with unique names, arranged as the user specifies.
   Hides the details of which frames nest in which, and
   the exact path for referencing frames from each other.

   Right now, supports only a vertical layout of a bunch of frames.
*/
class HTMLFrames {

  void AddFrame( String frameName, String frameTitle, String frameBody );

  void WriteFrames( String dir );
  
};  // class HTMLFrames

#endif  

/**
   Create a table showing the genomic unipaths.

   implementation notes:

   so, 

*/
void ShowGenomicUnipaths( nbases_t K,
			  genome_part_id_t GPART, genome_part_pos_t FROM, genome_part_pos_t TO,
			  const vecKmerPath& genome_paths, const vecKmerPath& genome_paths_rc,
			  const vec< big_tagged_rpint >& genome_pathsdb,
			  const vecKmerPath& genome_unipaths,
			  const vec< big_tagged_rpint >& genome_unipathsdb,
			  String HTML_DIR ) {

  ForceAssert( 0 <= GPART && static_cast<size_t>(GPART) < genome_paths.size() );
  ForceAssert( 0 <= FROM && FROM <= TO && TO < Kmers2Bases( genome_paths[ GPART ].KmerCount(), K ) );

  // Find the genome path interval containing our start position.
  const KmerPath& gpath = genome_paths[ GPART ];
  KmerPathLoc curLoc = gpath.Begin() + FROM;
  KmerPathLoc toLoc = gpath.Begin() + TO;

  vec< genome_placement > genome_unipath_placements;
  vec< html_t > genome_unipath_td;
  vec< html_attrs_t > genome_unipath_td_attrs;

  vec< unipath_interval_id_t > unipathIntervalId;
  genome_part_pos_t curGenomePos = FROM;
  while ( curLoc <= toLoc ) {
    Contains( genome_unipathsdb, curLoc.GetKmer(), unipathIntervalId );
    ForceAssert( unipathIntervalId.solo() );

    unipath_id_t uid = genome_unipathsdb[ unipathIntervalId.front() ].ReadId();

    cout << "uid=" << uid << " len=" << genome_unipaths[ uid ].KmerCount() << endl;

    if ( genome_unipaths[ uid ].KmerCount() > 100 ) {
      genome_unipath_placements.push( uid, genome_unipaths[ uid ].KmerCount(), GPART,
				      curGenomePos, curGenomePos + genome_unipaths[ uid ].KmerCount(),
				      False /* not rc */, copy_num_t(1) );
      genome_unipath_td.push( ToString( uid ) );
      genome_unipath_td_attrs.push( "" );
    }
    
    curLoc += genome_unipaths[ uid ].KmerCount();
    curGenomePos += genome_unipaths[ uid ].KmerCount();
  }

  cout << " computing placements..." << endl;

  Ofstream( html, HTML_DIR + "/unipaths.html" );
  html << HTMLHead( "Unipaths" );
  
  html <<  ShowGenomePlacements( genome_unipath_placements,
				 genome_unipath_td,
				 genome_unipath_td_attrs, GPART,
				 genome_paths[ GPART ].KmerCount() ) << endl;
  html << HTMLTail();
}


int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int_OrDefault(MAX_READS,3000);
  CommandArgument_String_OrDefault_Doc(GENOME, "genome", "the reference genome" );
  CommandArgument_Int_OrDefault_Doc(GPART, 0, "genome part (genome.fasta record) from"
				    " which we take the window." );
  CommandArgument_Int_Doc(FROM, "start coordinate of the window on the specified genome part");
  CommandArgument_Int_Doc(TO, "end coordinate of the window on the specified genome part");
  CommandArgument_Int_OrDefault_Doc( K, 20, "the kmer size used for read paths and unipaths" );
  CommandArgument_String_OrDefault_Doc( UNIPATHS, "",
					"the unipaths" );
  CommandArgument_String_OrDefault_Doc( UNIPATH_PLACEMENTS, "",
					"placements of unipaths on the reference,"
					" as found by PathsToLocs" );
  CommandArgument_String_OrDefault_Doc( HTML_DIR, "rv", "name of html dir for output" );
  CommandArgument_Bool_OrDefault_Doc( SHOW_UNIPATHS, True, "show unipaths approximated from the reads" );
  CommandArgument_Bool_OrDefault_Doc( SHOW_READS, True, "show reads" );
  CommandArgument_Int_OrDefault_Doc( MIN_SEP, -1, "only show reads with at least this separation" );
  CommandArgument_Int_OrDefault_Doc( MAX_SEP, -1, "only show reads with at most this separation"  );
  CommandArgument_String_OrDefault_Doc( SHOW_THESE_READS, "",
					"show only these read ids (from a vec<read_id_t> in a file)" );
  CommandArgument_String_OrDefault( READS, "reads_orig" );
  EndCommandArguments;

  const String data_dir = PRE + "/" + DATA;
  const String run_dir = PRE + "/" + DATA + "/" + RUN;
  String KS = ToString(K);

  if ( !IsDirectory( HTML_DIR ) )
    Mkdir777( HTML_DIR );

  vec< read_id_t > readsToShow;
  if ( SHOW_THESE_READS.nonempty() ) {
    BinaryReader::readFile( SHOW_THESE_READS, &readsToShow );
    UniqueSort( readsToShow );
  }

  {
    cout << Date() << " Loading read aligns..." << endl;
    vec<GaplessAlign> aligns;
    BinaryReader::readFile( data_dir + "/reads_orig.aligns_to_ref.uniq.gapless_aligns", &aligns );
    cout << Date() << " Loaded " << aligns.size() << " aligns." << endl;

    int allFrom = -1, allTo = -1;
    for ( align_id_t alignId = 0; alignId < aligns.isize(); alignId++ ) {
      const GaplessAlign& a = aligns[ alignId ];
      if ( alignId == 0 ) {
	allFrom = a.StartOnTarget();
	allTo = a.EndOnTarget();
      } else {
	update_min( allFrom, a.StartOnTarget() );
	update_max( allTo, a.EndOnTarget() );
      }
    }
    PRINT2( allFrom, allTo );

    
    cout << Date() << " Loading pair information..." << endl;
    int nreads = MastervecFileObjectCount( data_dir + "/" + READS + ".fastb" );
    PRINT( nreads ); 
    vec<read_pairing> pairs;
    ReadPairsFile( data_dir + "/" + READS + ".pairto", pairs );
    
    cout << Date() << " Loaded " << pairs.isize() << " pairs." << endl;
    ForceAssertEq( nreads, 2 * pairs.isize() );

     vec< pair_id_t > read2pair( nreads, -1 );
     for ( pair_id_t pairId = 0; pairId < pairs.isize(); pairId++ ) {
       const read_pairing& p = pairs[ pairId ];
       ForceAssertLe( 0, p.id1 );
       ForceAssertLt( p.id1, nreads );
       ForceAssertLe( 0, p.id2 );
       ForceAssertLt( p.id2, nreads );
       
       read2pair[ p.id1 ] = read2pair[ p.id2 ] = pairId;
     }

     cout << Date() << " Built read2pair" << endl;
    
    vec<GaplessAlign> selectedAligns;
    for ( align_id_t alignId = 0; alignId < aligns.isize(); alignId++ ) {
      const GaplessAlign& a = aligns[ alignId ];
      ForceAssertLe( 0, a.QueryId() );
      ForceAssertLt( a.QueryId(), nreads );
      nbases_t sep = pairs[ read2pair[ a.QueryId() ] ].sep;
      if ( ( FROM < 0 || FROM < a.StartOnTarget() ) &&
	   ( TO < 0  ||  a.EndOnTarget() < TO ) &&
	   ( MIN_SEP < 0  ||  sep > MIN_SEP ) &&
	   ( MAX_SEP < 0  ||  sep < MAX_SEP ) &&
	   ( readsToShow.empty() || BinMember( readsToShow, a.QueryId() ) ) )
	selectedAligns.push_back( a );
    }

    PRINT( selectedAligns.size() );

    if ( selectedAligns.isize() > MAX_READS ) {
      vec<GaplessAlign> selectedAlignsSample( MAX_READS );
      random_sample( selectedAligns.begin(), selectedAligns.end(),
		     selectedAlignsSample.begin(), selectedAlignsSample.end() );
      selectedAligns.swap( selectedAlignsSample );
    }

    PRINT( selectedAligns.size() );

    vec< html_t > readTd( selectedAligns.isize() );
    vec< html_attrs_t > readAttrs( selectedAligns.isize() );
    for ( align_id_t alignId = 0; alignId < selectedAligns.isize(); alignId++ ) {
      dot_color_t fg, bg;
      DotUtils::RandomFgBg( fg, bg );
      readAttrs[ alignId ] = " BGCOLOR=\"" + bg + "\" ";
      readTd[ alignId ] = "<FONT color=\"" + fg + "\"> </FONT>";
    }
    Ofstream( html, HTML_DIR + "/uniqueReads.html" );
    
    html << HTMLHead( "Unipaths" );
    
    vecKmerPath genome_paths( data_dir + "/genome.paths.k" + KS );

    cout << Date() << " placing " << selectedAligns.isize() << " reads..." << endl;

    BinaryWriter::writeFile( "falsedBridgeReadsAligns.dat", selectedAligns );
    BinaryWriter::writeFile( "falsedBridgeReadsTd.dat", readTd );
    BinaryWriter::writeFile( "falsedBridgeReadsAttrs.dat", readAttrs );

    html <<  ShowGenomePlacements( selectedAligns,
				   readTd,
				   readAttrs, GPART,
				   genome_paths[ GPART ].KmerCount() ) << endl;
    html << HTMLTail();
    
    exit( 0 );
    
  }

  
  {

    cout << Date() << " Loading genome unipaths..." << endl;
    vecKmerPath genome_paths( run_dir + "/genome.paths.k" + KS );
    vecKmerPath genome_paths_rc( run_dir + "/genome.paths_rc.k" + KS );
    BREAD2( run_dir + "/genome.pathsdb_big.k" + KS, vec<big_tagged_rpint>, genome_pathsdb );
    vecKmerPath genome_unipaths( run_dir + "/genome.unipaths.k" + KS );
    BREAD2( run_dir + "/genome.unipathsdb_big.k" + KS, vec<big_tagged_rpint>, genome_unipathsdb );
    CheckUnipathSanity( genome_unipaths, genome_unipathsdb );
    cout << Date() << " Loaded genome unipaths." << endl;
    ShowGenomicUnipaths( K, GPART, FROM, TO, genome_paths, genome_paths_rc, genome_pathsdb, genome_unipaths, genome_unipathsdb,
			 HTML_DIR );
    exit(1);
    
  }


  CpIfNewer( "auxfiles/utils.js", HTML_DIR + "/utils.js", True /* ignore errors */ );
  CpIfNewer( "auxfiles/tablecloth.js", HTML_DIR + "/tablecloth.js", True );
  CpIfNewer( "auxfiles/tablecloth.css", HTML_DIR + "/tablecloth.css", True );

  ForceAssert( IsRegularFile( run_dir + "/reads.fastb" ) );
  ForceAssert( IsRegularFile( run_dir + "/reads.pairto" ) );
  ForceAssert( IsDirectory( HTML_DIR ) );
  ForceAssertLe( FROM, TO );
  
  ForceAssert( UNIPATHS.nonempty() && UNIPATH_PLACEMENTS.nonempty() && SHOW_UNIPATHS ||
	       UNIPATHS.empty() && UNIPATH_PLACEMENTS.empty() && !SHOW_UNIPATHS );
  
  cout << Date() << ": Loading genome..." << endl;
  vecbasevector genome( data_dir + "/genome.fastb" );
  cout << Date() << ": Loading reads..." << endl;
  //vecbasevector reads_new( run_dir + "/reads.fastb" );
  vecbasevector reads( run_dir + "/reads.fastb" );

  vecqualvector quals;
  if ( IsRegularFile( run_dir + "/reads.qualb" ) )
    quals.ReadAll( run_dir + "/reads.qualb"  );
  VecFourBaseVec intensities;
  if ( IsRegularFile( run_dir + "/reads.intensities" ) )
    intensities.ReadAll( run_dir + "/reads.intensities" );

  //BREAD2( run_dir + "/reads.id_map", vec<int>, id_map );

  Ofstream( html, HTML_DIR + "/reads.html" );
  html << HTMLHead( "Window " + ToString( GPART ) + ":" +
		    ToString( FROM ) + "-" + ToString( TO ),
		    
		    "<script type=\"text/javascript\">"
		    "  function PleaseScroll( hor ) {"
		    "var curLeft = (window.pageXOffset)?(window.pageXOffset):(document.documentElement)?document.documentElement.scrollLeft:document.body.scrollLeft;"
		    "  if ( hor != curLeft ) {"
		    "      var curTop = "
		    "             (window.pageYOffset)?(window.pageYOffset):(document.documentElement)?document.documentElement.scrollTop:document.body.scrollTop;"
		    "      scrollTo( hor, curTop );"
		    "  }"
		    " }"
		    " "
		    "function scrollL()"
		    "{"
		    "var left = (window.pageXOffset)?(window.pageXOffset):(document.documentElement)?document.documentElement.scrollLeft:document.body.scrollLeft;"
		    "   window.top.ScrollAll( left, self );"
		    "}"
		    " window.onscroll = scrollL; "
		    "</script>"
		    "<style>"
		    "table, td{ font: 100% Arial, Helvetica, sans-serif; }"
		    "table { width: 100%; border-collapse: collapse; margin: 1em 0;}"
		    "th,td { text-align: center; padding: .3em; border: 1px solid #fff; }"
		    "</style>"
		    );

  nbases_t NCOLS = TO - FROM + 1;
  
  html << "<table cellspacing=\"0\" >\n";

  html << "<colgroup span=\"" << NCOLS << "\">\n";
  for ( genome_part_pos_t gpos = FROM; gpos <= TO; gpos++ ) {
    html << "<col align=\"center\" valign=\"center\" style=\"background-color:" << (gpos % 2 ? "#f8fbfc" : "#e5f1f4" )
	 << "\"></col>\n";
  }
  html << "</colgroup>\n";

  html << "<thead>\n";
  html << "<tr><td style=\"text-align: left; background:#f8fbfc\" colspan=\"" << NCOLS << "\">Reference</td></tr>\n";

  html << "<tr>";
  for ( genome_part_pos_t gpos = FROM; gpos <= TO; gpos++ ) {
    html << "<th>" << as_base( genome[ GPART ][gpos] ) << "</th>";
  }
  html << "</tr>\n";
  html << "</thead>\n";
  html << "<tbody>\n";

  int nreads = reads.size();

  vec<read_location> readlocs;
  vec<int> readlocs_index;
  cout << Date() << ": Loading true read locs..." << endl;
  READX( run_dir + "/reads.ref.locs", readlocs );
  readlocs_index.resize(nreads);

  for ( int i = 0; i < readlocs.isize( ); i++ ) {
    const read_location& rl = readlocs[i];
    ForceAssertLe( 0, rl.ReadId() );
    ForceAssertLt( rl.ReadId(), nreads );
    readlocs_index[ rl.ReadId( ) ] = i;
  }

  
  
  {
    Ofstream( fhtml, HTML_DIR + "/genomeAndReads.html" );
    fhtml << HTMLHeadFrames( "Window",
			    "<script type=\"text/javascript\">"
			    " function ScrollAll( hor, asker ) { "
			    "   if ( asker != genome ) { genome.PleaseScroll( hor ); } "
			    "   if ( asker != reads ) { reads.PleaseScroll( hor ); } "
			    "}"
			    "</script>"
			    );
    fhtml << "<frameset rows=\"5%,95%\" >" << endl;
    fhtml << "<frame name=\"genome\" id=\"genome\" src=\"reads.html\" />\n";
    fhtml << "<frame name=\"reads\" id=\"reads\" src=\"reads.html\" />\n";
    fhtml << "</frameset>"<< endl;
    
    fhtml << HTMLTailFrames();
  }

  
  if ( SHOW_UNIPATHS ) {

  vecbasevector bases( run_dir + "/reads.fastb" );
  vecKmerPath paths( run_dir + "/reads.paths.k" + KS );
  vecKmerPath paths_rc( run_dir + "/reads.paths_rc.k" + KS );
  BREAD2( run_dir + "/reads.pathsdb.k" + KS, vec<tagged_rpint>, pathsdb );
  vecKmerPath unipaths( run_dir + "/reads.unipaths.k" + KS );
  cout << Date() << " Read " << unipaths.size() << " unipaths." << endl;
  
    

    vec<genome_placement> unipath_placements;
    String PATHS = run_dir + "/reads.unipaths.k" + ToString(K);
    //  vecKmerPath unipaths;
    
    BinaryReader::readFile( UNIPATH_PLACEMENTS, &unipath_placements );
    //unipaths.ReadAll( PATHS );

  vecbasevector unibases;
  {    KmerBaseBroker kbb;
  kbb.Initialize( K, bases, paths, paths_rc, pathsdb );
  unibases.reserve(unipaths.size());
  for ( size_t i = 0; i < unipaths.size( ); i++ )
    unibases.push_back( kbb.Seq( unipaths[i] ) );    }

    
    html << "<tr><td  style=\"text-align: left; background:#f8fbfc\" colspan=\"" << NCOLS << "\">Unipaths:</td></tr>\n";
    
    ForceAssert( unipaths.size() == unibases.size() );
    
    for ( int rc = 0; rc < 2; rc++ ) {
      for ( int readLocId = 0; readLocId < unipath_placements.isize(); readLocId++ ) {
	const genome_placement& rloc = unipath_placements[ readLocId ];
	ForceAssertLe( 0, rloc.GetReadId() );
	ForceAssertLt( static_cast<size_t>(rloc.GetReadId()), unipaths.size() );
	ForceAssertEq( rloc.GetKmerCount(), unipaths[ rloc.GetReadId() ].KmerCount() );
	
	if ( rloc.GetGenomeId() == GPART  &&
	     ( rc == 0 && rloc.IsFw()  ||
	       rc == 1 && rloc.IsRc() )
	     && IntervalsOverlap( rloc.GetStartOnGenome(), rloc.GetEndOnGenome(),
				  FROM, TO ) ) {
	  
	  html << "<tr>";
	  
	  cout << "unipath: " << rloc.GetStartOnGenome() << " - " << rloc.GetEndOnGenome() << endl;
	  
	  for ( genome_part_pos_t genomePos = FROM; genomePos <= TO; genomePos++ ) {
	    html << "<td>";
	    
	    read_pos_t posInRead = rloc.IsFw() ? genomePos - rloc.GetStartOnGenome() :
	      rloc.GetStartOnGenome() + rloc.GetKmerCount() + K - 1 - genomePos;
	    
	    if ( posInRead >= 0  &&  posInRead < (read_pos_t)(rloc.GetKmerCount() + K - 1) ) {
	      Bool atReadStart = ( posInRead == 0 );
	      Bool atReadEnd = ( posInRead == (read_pos_t)(rloc.GetKmerCount() + K - 2) );
	      
	      if ( rloc.IsFw() && atReadStart  ||  rloc.IsRc() && atReadEnd )
		html << "[";
	      html << as_base( unibases[ rloc.GetReadId() ][ posInRead ] );
	      if ( rloc.IsFw() && atReadEnd  ||  rloc.IsRc() && atReadStart )
		html << "]";
	    } else {
	      html << ".";
	    }
	    
	    html << "</td>";
	  }
	  
	  
	  html << "</tr>\n";
	  
	  //	static int cnt = 0;
	  //	if ( cnt++ > 5 )
	  //	  exit(1);
	  
	}  // if this read intersects the [FROM,TO] interval of the reference
      }  // for each read placement
      
      if ( rc == 0 ) {
	html << "<tr>";
	
	for ( genome_part_pos_t genomePos = FROM; genomePos <= TO; genomePos++ ) {
	  html << "<td>";
	  html << "-";
	  html << "</td>";
	}
	
	html << "</tr>";
      }
      
    }  // for rc=0..1
    
  }
  
  if ( SHOW_READS ) {

    html << "<tr><td style=\"text-align: left; background:#f8fbfc\" colspan=\"" << NCOLS << "\">Reads:</td></tr>\n";

    Sort( readlocs );
    
    for ( int rc = 0; rc < 2; rc++ ) {
      for ( int readLocId = 0; readLocId < readlocs.isize(); readLocId++ ) {
	const read_location& rloc = readlocs[ readLocId ];
	if ( !( 0 <= rloc.ReadId() && static_cast<size_t>(rloc.ReadId()) < reads.size() &&
		rloc.LengthOfRead() == reads[ rloc.ReadId() ].size() ) ) {
	  cout << "skipping bad read: " << readLocId << ": " << rloc << endl;
	  continue;
	}
	ForceAssertLe( 0, rloc.ReadId() );
	ForceAssertLt( static_cast<size_t>(rloc.ReadId()), reads.size() );
	ForceAssertEq( rloc.LengthOfRead(), reads[ rloc.ReadId() ].size() );
	if ( rloc.Contig() == GPART  &&
	     ( rc == 0 && rloc.Fw()  ||
	       rc == 1 && rloc.Rc() )
	     && IntervalsOverlap( rloc.StartOnContig(), rloc.StopOnContig(),
				  FROM, TO ) ) {
	  
	  html << "<tr>";
	  
	  cout << "read: " << rloc.StartOnContig() << " - " << rloc.StopOnContig() << endl;
	  
	  for ( genome_part_pos_t genomePos = FROM; genomePos <= TO; genomePos++ ) {

	    read_pos_t posInRead = rloc.Fw() ? genomePos - rloc.StartOnContig() :
	      rloc.StartOnContig() + rloc.LengthOfRead() - 1 - genomePos;
	    
	    if ( posInRead >= 0  &&  posInRead < (read_pos_t)rloc.LengthOfRead() ) {

	      html << "<td ";
	      if ( !quals.empty() )
		html << " title=\"qual=" << int(quals[ rloc.ReadId() ][ posInRead ]) << "\" ";
	      html << "  >";
	      
	      Bool atReadStart = ( posInRead == 0 );
	      Bool atReadEnd = ( posInRead == (read_pos_t)(rloc.LengthOfRead()-1) );
	      
	      if ( rloc.Fw() && atReadStart  ||  rloc.Rc() && atReadEnd )
		html << "[";

	      base_t genomeBase = genome[ GPART ][genomePos];
	      if ( rloc.Rc() )
		genomeBase = GetComplementaryBase( genomeBase );
	      html << ( reads[ rloc.ReadId() ][ posInRead ] == genomeBase  ? '~' : as_base( reads[ rloc.ReadId() ][ posInRead ] ) );
	      
	      if ( rloc.Fw() && atReadEnd  ||  rloc.Rc() && atReadStart )
		html << "]";
	      
	      html << "</td>";
	    } else {
	      html << "<td>.</td>";
	    }
	    
	  }
	  
	  
	  html << "</tr>\n";
	  
	  //	static int cnt = 0;
	  //	if ( cnt++ > 5 )
	  //	  exit(1);
	  
	}  // if this read intersects the [FROM,TO] interval of the reference
      }  // for each read placement
      
      if ( rc == 0 ) {
	html << "<tr>";
	
	for ( genome_part_pos_t genomePos = FROM; genomePos <= TO; genomePos++ ) {
	  html << "<td>";
	  html << "-";
	  html << "</td>";
	}
	
	html << "</tr>";
      }
      
    }  // for rc=0..1
    
  }
  
  html << "</tbody>\n";
  html << "</table>\n";
  
  
  html << HTMLTail();
}



