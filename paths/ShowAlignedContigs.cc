///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ShowAlignedContigs: a tool to display an alignment of contigs or scaffolds to
// the reference sequence.

// MakeDepend: dependency QueryLookupTable

#include "Basevector.h"
#include "Fastavector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "PrintAlignment.h"
#include "Superb.h"
#include "Vec.h"
#include "lookup/LookAlign.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_Int(K);
     CommandArgument_String_OrDefault(SCAFFOLDS_IN, "linear_scaffolds0.patched");
     CommandArgument_String_OrDefault(PATCHDIR, "patch");
     CommandArgument_String_OrDefault_Doc(TIGS, "", "if unspecified, process all "
          "contigs; otherwise it is one of the following: \n(a) a list of contig "
          "ids (in ParseIntSet format) or \n(b) the letter s followed by a list of "
          "scaffolds or \n(c) s<scaffold id>.<list of indices of contigs in the "
          "scaffold");
     CommandArgument_String_OrDefault_Doc(TARGETS, "", "for assessment, list of "
          "reference targets to use, in ParseIntSet format, to speed up alignment");
     CommandArgument_Bool_OrDefault(VERBOSE, True);
     CommandArgument_String_OrDefault(TMP_FILE, "ShowAlignedContigs.fasta");
     CommandArgument_String_OrDefault(GENOME, "genome");
     EndCommandArguments;

     // Begin.

     double clock = WallClockTime( );

     // Define directories, etc.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     cout << Date( ) << ": " << run_dir << endl;
     String pdir = run_dir + "/" + PATCHDIR;
     String ch_head = pdir + "/PostPatcher." + SCAFFOLDS_IN + ".";

     // Note some files that are needed.

     String TIGS_file = ch_head + "TIGS";

     // Get total number of contigs.

     size_t ntigs = MastervecFileObjectCount( 
          sub_dir + "/" + SCAFFOLDS_IN + ".contigs.vecfasta" );

     // Parse arguments.

     vec<int> tigs;
     if ( TIGS == "" )
     {    for ( size_t j = 0; j < ntigs; j++ )
               tigs.push_back(j);    }
     else if ( TIGS.Contains( "s", 0 ) )
     {    TIGS = TIGS.After( "s" );
          vec<superb> scaffolds;
          ReadSuperbs( sub_dir + "/" + SCAFFOLDS_IN + ".superb", scaffolds );
          if ( TIGS.Contains( "." ) )
          {    int scaffold = TIGS.Before( "." ).Int( );
               ForceAssertLt( scaffold, scaffolds.isize( ) );
               vec<int> spos;
               ParseIntSet( TIGS.After( "." ), spos );
               for ( int j = 0; j < spos.isize( ); j++ )
                    tigs.push_back( scaffolds[scaffold].Tig( spos[j] ) );    }
          else
          {    vec<int> s;
               ParseIntSet( TIGS, s );
               for ( int i = 0; i < s.isize( ); i++ )
               {    int scaffold = s[i];
                    ForceAssertLt( scaffold, scaffolds.isize( ) );
                    for ( int j = 0; j < scaffolds[scaffold].Ntigs( ); j++ )
                         tigs.push_back( scaffolds[scaffold].Tig(j) );    }    }    }
     else ParseIntSet( TIGS, tigs );
     int used_tigs = tigs.size( );

     // Fetch contigs.

     vec<fastavector> orig;
     for ( int it = 0; it < tigs.isize( ); it++ )
     {    int tig = tigs[it];
          vecfastavector contigs;
          contigs.ReadOne( 
               sub_dir + "/" + SCAFFOLDS_IN + ".contigs.vecfasta", tig );
          orig.push_back( contigs[0] );    }

     // Evaluate contigs.

     vec<fastavector> combo;
     {    Ofstream( out, TMP_FILE );
          for ( int i = 0; i < used_tigs; i++ )
               combo.push_back( orig[i] );
               for ( int i = 0; i < combo.isize( ); i++ )
                    combo[i].Print( out, ToString(i) );    }
     cout << "\n" << Date( ) << ": aligning to reference" << endl;
     String TARGETS_arg = ( TARGETS == "" ? ""
          : " TARGETS_TO_PROCESS=\"" + TARGETS + "\"" );
     fast_pipe_ifstream in( "QueryLookupTable K=12 MM=12 MC=0.15 "
          "SEQS=" + TMP_FILE + " L=" + data_dir + "/" + GENOME + ".lookup "
          "PARSEABLE=True" + TARGETS_arg );
     String line;
     vec<look_align> aligns;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( !line.Contains( "QUERY", 0 ) ) continue;
          look_align la;
          la.ReadParseable(line);
          aligns.push_back(la);    }
     cout << Date( ) << ": done" << endl;
     vec< vec<look_align> > old_aligns(used_tigs);
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    int id = aligns[i].query_id;
          old_aligns[id].push_back( aligns[i] );    }
     vec< vec<look_align> > old_aligns_best(used_tigs);
     for ( int i = 0; i < used_tigs; i++ )
          if ( old_aligns[i].nonempty( ) ) old_aligns_best[i] = old_aligns[i];
     vec<int> tids;
     for ( int i = 0; i < used_tigs; i++ )
     {    for ( int j = 0; j < old_aligns_best[i].isize( ); j++ )
               tids.push_back( old_aligns_best[i][j].target_id );    }
     UniqueSort(tids);
     vecbasevector targets;
     targets.Read( data_dir + "/" + GENOME + ".fastb", tids );
     cout << "\n=============================================================="
          << "======================\n";
     cout << "\nALIGNMENTS OF CONTIGS\n" << endl;
     int64_t total_bases = 0, total_errs = 0;
     int unaligneds = 0, incompletes = 0;
     for ( int i = 0; i < tigs.isize( ); i++ )
     {    int tig = tigs[i];
          cout << "==========================================================="
               << "=========================\n" << "\nCONTIG " << tig << "\n" 
               << endl;
          if ( old_aligns_best[i].empty( ) ) 
          {    cout << "(no alignment found)" << endl;
               unaligneds++;    }
          else
          {    const look_align& la = old_aligns_best[i][0];
               la.PrintReadableBrief(cout);
               la.PrintVisual( cout, combo[la.query_id],
                    targets[ Position( tids, la.target_id ) ] );    
               total_errs += la.Errors( );
               total_bases += la.query_length;
               if ( la.pos1( ) > 0 || la.Pos1( ) < (int) la.query_length )
                    incompletes++;    }    }

     // Print summary stats.

     cout << "==========================================================="
          << "=========================\n\n";
     cout << "total bases = " << ToStringAddCommas(total_bases) << endl;
     cout << "error rate = " << PERCENT_RATIO( 2, total_errs, total_bases ) << endl;
     cout << "unaligned = " << unaligneds << endl;
     cout << "incomplete alignments = " << incompletes << endl;
     cout << "unaligned+incomplete rate = " << setprecision(2)
          << double(unaligneds+incompletes) / ( double(total_bases)/1000000.0 )
          << " per Mb" << endl;    }
