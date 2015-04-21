/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "kmers/KmerRecord.h"
#include "graph/Digraph.h"
#include "PrettyPrintTable.h"
#include "kmer_freq/KAdjGraph.h"
#include "paths/Unipath.h"
#include "reporting/PrintLengthsNstats.h"
#include "feudal/BinaryStream.h"

template <int K>
void Stats( const String &base,
	    const String &unibases,
	    const String &unigraph,
	    ostream &out,
	    bool verbose );

/**
 * BasicHyperStats.cc
 *
 * Basic digraph-based assembly stats.
 *
 * Input (all files are relative to PRE/DATA/RUN):
 *    ../BASE.UNIBASES.kK.fastb: vecbasevector of unibases
 *    ../BASE.UNIGRAPH.kK (digraph)
 * Output:
 *    ../BASE.BasicHyperStats.kK.out if ARCHIVE=True (or else to cout)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_Int( K );
  CommandArgument_String_OrDefault( BASE, "reads" );
  CommandArgument_String_OrDefault( UNIBASES, "unibases" );
  CommandArgument_String_OrDefault( UNIGRAPH, "unigraph" );
  CommandArgument_Bool_OrDefault( ARCHIVE, True );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  EndCommandArguments;
  
  // File names.
  String strK = ToString( K );
  String base = PRE + "/" + DATA + "/" + RUN + "/" + BASE;
  String out_file = base + ".BasicHyperStats.k" + strK + ".out";

  // Out stream.
  ofstream archive_out;
  if ( ARCHIVE ) archive_out.open( out_file.c_str( ) );
  ostream &out = ARCHIVE ? archive_out : * (ostream *) &cout;
  if ( ARCHIVE ) PrintCommandPretty( out );

  // Run.
#define STATS( K_ ) Stats<K_> ( base, UNIBASES, UNIGRAPH, out, VERBOSE );
  DISPATCH_ON_K( K, STATS );
  
  // Close stream.
  if ( ARCHIVE ) archive_out.close( );
}

/**
 * Stats
 */
template <int K>
void Stats( const String &base,
	    const String &unibases,
	    const String &unigraph,
	    ostream &out,
	    bool verbose = false)
{
  // File names.
  String strK = "k" + ToString( K );
  String unibases_file = base + "." + unibases + "." + strK + ".fastb";
  String unigraph_file = base + "." + unigraph + "." + strK;
  
  ForceAssert( IsRegularFile( unibases_file ) );
  ForceAssert( IsRegularFile( unigraph_file ) );

  // Load.
  out << Date( ) << ": loading unibases" << endl;
  vecbvec bases;
  bases.ReadAll( unibases_file );

  out << Date( ) << ": loading unigraph" << endl;
  digraph graph;
  BinaryReader::readFile( unigraph_file, &graph );

  // Get components.
  vec< vec<int> > components;
  graph.Components( components );
  
  vec<int> counts( components.size( ), 0 );
  vec<int> lens( components.size( ), 0 );
  for (int ii=0; ii<(int)components.size( ); ii++) {
    const vec<int> &comp = components[ii];
    counts[ii] = comp.size( );
    for (int jj=0; jj<(int)comp.size( ); jj++)
      lens[ii] += bases[ comp[jj] ].size( );
  }
  
  // Print results.
  out << "\n"
      << "connected components: " << components.size( ) << "\n"
      << "total unipath length: " << BigSum( lens ) << "\n"
      << "number of vertices:   " << BigSum( counts ) << "\n"
      << "\n";

  out << "N-stats of lengths of unipaths:\n";
  PrintLengthsNstats( lens, "", out ); 
  out << "\n";

  out << "N-stats of counts of vertices in each component:\n";
  PrintLengthsNstats( counts, "", out );
  out << "\n";
  
  // Verbose mode.
  if ( ! verbose ) return;

  const int lastcomp = components.size( ) - 1;
  for (int comp_id=0; comp_id<(int)components.size( ); comp_id++) {
    const vec<int> &comp = components[comp_id];
    const int lastjj = comp.size( ) - 1;

    vec<int> complens( comp.size( ), 0 );
    for (int ii=0; ii<(int)comp.size( ); ii++)
      complens[ii] = bases[ comp[ii] ].size( );
    sort( complens.begin( ), complens.end( ) );

    out << "C" << comp_id << "/" << lastcomp
	<< "   " << BigSum( complens )
	<< " bp total   " << N50( complens )
	<< " bp N50   " << comp.size( )
	<< ( comp.size( ) == 1 ? " vertex\n" : " vertices\n" );
    
    for (int jj=0; jj<(int)comp.size( ); jj++) {
      const int vert_id = comp[jj];
      vec<int> sources = graph.From( vert_id );
      vec<int> sinks = graph.To( vert_id );
      
      out << "  c" << comp_id << "/" << lastcomp
	  << ".v" << jj << "/" << lastjj
	  << "   u" << vert_id 
	  << "." << bases[vert_id].size( )
	  << "bp\n";
      
      vec< vec<String> > table;
      vec<String> line( 2 );
      line[0] = ToString( sinks.size( ) ) + "-in";
      line[1] = ToString( sources.size( ) ) + "-out";
      table.push_back( line );
      int msize = Max( sinks.size( ), sources.size( ) );
      for (int kk=0; kk<msize; kk++) {
	line[0] = "";
	if ( kk < (int)sinks.size( ) ) {
	  int sid = sinks[kk];
	  int len = bases[sid].size( );
	  line[0] = ToString( sid ) + "." + ToString( len ) + "bp";
	}
	line[1] = "";
	if ( kk < (int)sources.size( ) ) {
	  int sid = sources[kk];
	  int len = bases[sid].size( );
	  line[1] = ToString( sid ) + "." + ToString( len ) + "bp";
	}
	table.push_back( line );
      }
      BeautifyTable( table );
      for (int kk=0; kk<(int)table.size( ); kk++) {
	out << "    ";
	for (int mm=0; mm<(int)table[kk].size( ); mm++)
	  out << table[kk][mm] << "   ";
	out << "\n";
      }
    }
    out << "\n";
  }

  // Done.
  out << Date( ) << ": done" << endl;
}

/**
 * Instantiate templates
 */
#define INSTANTIATE_STATS( K_, dummy )		\
template void Stats<K_>( const String &,	\
			 const String &,	\
			 const String &,	\
			 ostream &,             \
			 bool )

FOR_ALL_K( INSTANTIATE_STATS,unused);

