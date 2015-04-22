/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Program: EvalHyper
//
// Carry out the HyperKmerPath evaluation that 
// occurs at the tail end of LocalizeReads.  You can run this command with all the 
// LocalizeReads arguments.  

#include "Basevector.h"
#include "Bitvector.h"
#include "MainTools.h"
#include "ReadPairing.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "paths/AlignHyperKmerPath.h"
#include "paths/EvalUtils.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(INSTANCE, "");
     CommandArgument_String(SUBDIR);
     CommandArgument_String_OrDefault(WRUN, "run");
     CommandArgument_Bool_OrDefault(USE_TRUTH, True);
     CommandArgument_String_OrDefault_Doc(TPDIR, "/showTrusted", 
          "where to put the trusted path viewer");
     CommandArgument_String_OrDefault_Doc(ALT_DATA_DIR, "", 
          "alternate data_dir, allowing one to substitute a different reference");
     CommandArgument_String_OrDefault(HYPER, "hyper");
     CommandArgument_Bool_OrDefault(FILTER_ALIGNS, True);
     CommandArgument_Bool_OrDefault(DUMP_BASES, False);
     CommandArgument_Bool_OrDefault(DIPLOID, False);
     CommandArgument_Bool_OrDefault(REORDER, False);
     CommandArgument_Int_OrDefault_Doc( K, 20,
          "kmer size (must match the ones used elsewhere)." );
     CommandArgument_Int_OrDefault_Doc( MIN_TP_LEN, -1, 
          "minimum trusted path length" );
     CommandArgument_Bool_OrDefault(LIBSTATS, False);
     CommandArgument_String_OrDefault(REPORT_SUFFIX, "");
     EndCommandArguments;

     // Set up directories.

     SUBDIR += INSTANCE;
     String data_dir = PRE + "/" + DATA;
     if ( ALT_DATA_DIR != "" ) data_dir = ALT_DATA_DIR;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String wdata_dir = sub_dir;
     String wrun_dir = sub_dir + "/" + WRUN;

     if ( FirstLineOfFile( data_dir + "/ploidy" ).Int() == 2 ) {
       cout << "Setting DIPLOID to True." << endl;
       DIPLOID = True;
     }

     // Generate library stats.

     if (LIBSTATS)
     {    Ofstream( out, sub_dir + "/report.data" );
          vecbasevector reads( run_dir + "/reads.fastb" );
          vec<read_pairing> pairs;
          ReadPairsFile( run_dir + "/reads.pairto", pairs );
          vec< pair<int,int> > len_dev;
          for ( int i = 0; i < pairs.isize( ); i++ )
          {    int id1 = pairs[i].id1, id2 = pairs[i].id2;
               len_dev.push( pairs[i].sep + reads[id1].isize( ) 
                    + reads[id2].isize( ), pairs[i].sd );    }
          Sort(len_dev);
          vec<int> len, dev, count;
          for ( int i = 0; i < len_dev.isize( ); i++ )
          {    int j = len_dev.NextDiff(i);
               len.push_back( len_dev[i].first );
               dev.push_back( len_dev[i].second );
               count.push_back( j - i );
               i = j - 1;    }
          vec< vec<String> > rows;
          vec<String> row1;
          row1.push_back( "len", "dev", "count" );
          rows.push_back( row1 );
          for ( int i = 0; i < len.isize( ); i++ )
          {    vec<String> row;
               row.push_back( 
                    ToString(len[i]), ToString(dev[i]), ToString(count[i]) );
               rows.push_back(row);    }
          PrintTabular( out, rows, 2, "rrr" );    }

     // Load data.

     vecbasevector genome; // this is empty if !USE_TRUTH, which is okay
     vecbitvector genome_amb;
     if (USE_TRUTH) genome = vecbasevector( data_dir + "/genome.fastb" );
     if ( USE_TRUTH && IsRegularFile( data_dir + "/genome.fastamb" ) )
          genome_amb.ReadAll( data_dir + "/genome.fastamb" );
     String hkpFname = sub_dir + "/" + HYPER;
     PRINT( hkpFname );
     ForceAssert( IsRegularFile( hkpFname ) );
     HyperKmerPath *h = NULL;
     h = new HyperKmerPath( hkpFname );
     ForceAssertEq( K, h->K( ) );
     KmerBaseBroker* kbb = new KmerBaseBroker( wrun_dir, K );

     // Evaluate assembly, reordering to follow reference.

     EvaluateAssembly( *h, kbb, data_dir, wrun_dir, sub_dir, genome, genome_amb,
          DIPLOID, USE_TRUTH, FILTER_ALIGNS, False, REORDER, MIN_TP_LEN,
		       HYPER, REPORT_SUFFIX );

     if (REORDER) BinaryWriter::writeFile( sub_dir + "/" + HYPER + ".reordered", *h );

     // Dump bases.

     if (DUMP_BASES)
     {    for ( int i = 0; i < h->EdgeObjectCount( ); i++ )
               kbb->Seq( h->EdgeObject(i) ).Print( cout, BaseAlpha(i) );    }    }
