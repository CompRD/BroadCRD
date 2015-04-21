///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AssemblySummaryStats.  Generate a small file of critical statistics for an
// assembly for which there is a reference sequence.

#include "FastIfstream.h"
#include "MainTools.h"
#include <sstream>

// MakeDepend: dependency AssemblyAccuracy
// MakeDepend: dependency AssemblyCoverage
// MakeDepend: dependency ScaffoldAccuracy

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(SUBDIR);
     CommandArgument_String_OrDefault( SCAFFOLDS, "linear_scaffolds" );
     CommandArgument_String_OrDefault( OUT_HEAD, "summary_stats" );
     EndCommandArguments;

     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

     Ofstream( out, sub_dir + "/" + OUT_HEAD );

     String REFHEAD = data_dir + "/genome";
     String ASSEMBLY = sub_dir + "/" + SCAFFOLDS + ".fasta";
     String AC_HEAD = sub_dir + "/" + SCAFFOLDS;

     String REF = REFHEAD + ".fasta";
     String TAIL = " ASSEMBLY=" + ASSEMBLY + " NH=True";
     String ACH = " HEAD=" + AC_HEAD;

     fast_pipe_ifstream sa( "ScaffoldAccuracy REFHEAD=" + REFHEAD + TAIL );
     fast_pipe_ifstream ac( "AssemblyCoverage REF=" + REF + TAIL + ACH );
     fast_pipe_ifstream aa( "AssemblyAccuracy REF=" + REF + TAIL );

     String line;
     vec< vec<String> > rows;
     vec<String> row;

     while(1)
     {    getline( ac, line );
          ForceAssert( !ac.fail( ) );
          if ( !line.Contains( " covered", -1 ) ) continue;
          row.push_back( line.Before( " " ) );
          row.push_back( "of genome is covered by contigs (from AssemblyCoverage)" );
          rows.push_back(row);
          break;    }

     while(1)
     {    getline( sa, line );
          ForceAssert( !sa.fail( ) );
          if ( line.Contains( "N50 contig size = ", 0 ) ) break;    }
     row.clear( );
     double N50contig = line.Between( "= ", " " ).Int( ) / 1000.0;
     ostringstream outc;
     outc << setiosflags(ios::fixed) << setprecision(1) << N50contig;
     row.push_back( outc.str( ) );
     row.push_back( "N50 contig size in kb (from Scaffold Accuracy)" );
     rows.push_back(row);

     getline( sa, line );
     ForceAssert( !sa.fail( ) );
     ForceAssert( line.Contains( "N50 scaffold size ", 0 ) );
     row.clear( );
     double N50scaffold = line.Between( "= ", " " ).Int( ) / 1000.0;
     ostringstream outs;
     outs << setiosflags(ios::fixed) << setprecision(1) << N50scaffold;
     row.push_back( outs.str( ) );
     row.push_back( "N50 scaffold size in kb, excluding gaps "
          "(from Scaffold Accuracy)" );
     rows.push_back(row);

     getline( ac, line );
     ForceAssert( !ac.fail( ) );
     ForceAssert( line.Contains( " duplicated", -1 ) );
     row.clear( );
     row.push_back( line.Before( " " ) );
     row.push_back( "duplication fraction in contigs (from AssemblyCoverage)" );
     rows.push_back(row);

     getline( aa, line );
     ForceAssert( !aa.fail( ) );
     ForceAssert( line.Contains( " of contig bases are ambiguous", -1 ) );
     double amb = line.Before( "%" ).Double( ) * 100.0;
     ostringstream outa;
     outa << setiosflags(ios::fixed) << setprecision(2) << amb;
     row.clear( );
     row.push_back( outa.str( ) );
     row.push_back( "fraction of ambiguous bases, per 10^4 bases "
          "(from AssemblyAccuracy)" );
     rows.push_back(row);

     while(1)
     {    getline( aa, line );
          ForceAssert( !aa.fail( ) );
          if ( !line.Contains( "perfect", 0 ) ) continue;
          line = line.Between( "perfect", "%" );
          line.GlobalReplaceBy( " ", "" );
          row.clear( );
          row.push_back( line + "%" );
          row.push_back( "fraction of contig chunks that are perfect "
               "(from AssemblyAccuracy)" );
          rows.push_back(row);
          break;    }

     while(1)
     {    getline( sa, line );
          ForceAssert( !sa.fail( ) );
          if ( !line.Contains( "validation rate", 0 ) ) continue;
          row.clear( );
          row.push_back( line.After( "= " ) );
          row.push_back( "validation rate of scaffolds at 100 kb "
               "(from ScaffoldAccuracy)" );
          rows.push_back(row);
          break;    }

     PrintTabular( out, rows, 3 );
     PrintTabular( cout, rows, 3 );    }
