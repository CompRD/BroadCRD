///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: private

#include "CoreTools.h"
#include "FastIfstream.h"
#include "ParseSet.h"
#include "paths/long/DataSpec.h"
#include "paths/long/LoadAndCorrect.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/fosmid/FosmidPool.h"
#include "paths/long/local/Setup.h"

void AddTo( const Bool direct, vec<picard_spec>& specs, const String& flowcell, 
     const String& picard_run, const int lane, const String& lib )
{    specs.push( direct, flowcell, picard_run, lane, lib );    }

void AddTo( const Bool direct, vec<picard_spec>& specs, const String& flowcell, 
     const String& picard_run, const vec<int>& lanes, const String& lib )
{    for ( int lane : lanes ) specs.push(direct, flowcell, picard_run, lane, lib);    }

void AddTo( const Bool direct, vec<picard_spec>& specs, const String& flowcell, 
     const String& picard_run, const int lane, const vec<String>& libs )
{    for ( String lib : libs ) specs.push(direct, flowcell, picard_run, lane, lib);    }

Bool IfDel( vec<String>& dataset, const String& x )
{    if ( !Member( dataset, x ) ) return False;
     dataset.EraseValue(x);
     return True;    }

int ext_hh = 400;

void SetupIlluminaData( 
     // inputs:
     const long_data_spec& spec, const String& SAMPLE, const String& DATASET, 
     const vec<String>& regions, const double COVERAGE_CAP, const Bool HPOOL_ORIG,
     // outputs:
     vec<picard_spec>& specs, vec<String>& chrlist, String& region,
     // input and output:
     double& SELECT_FRAC,
     // logging:
     const long_logging& logc )
{
     // Set up for rhody.

     if ( SAMPLE == "rhody" ) // equivalent to build_micro ... miseq-alt
     {    
          // old data:

          /*
          String lh = "Solexa-121745.Tag_IlluminaHTKit";
          AddTo( False, specs, "A1PJD", "C1-508_2012-09-24_2012-10-04", 1, 
               { lh + "2", lh + "3", lh + "4", lh + "5" } );
          */

          // new data:

          AddTo( True, specs, "A6FWA", "C1-516_2013-12-20_2013-12-22", 1,
               "Solexa-208790" );
          // AddTo( True, specs, "A6VP5", "C1-516_2013-12-20_2013-12-22", 1,
          //      "Solexa-208790" ); // extra coverage
          // SELECT_FRAC = 0.4;

          chrlist.push_back( "chr1", "plasmid1", "chr2", "plasmid2", "plasmid3", 
               "plasmid4", "plasmid5" );    }

     // Set up for arabidopsis.

     if ( SAMPLE == "arabidopsis" )
     {    AddTo( True, specs, "H75EWADXX", "C1-516_2013-09-21_2013-10-03", 1,
               "Solexa-183228" );
          AddTo( True, specs, "H75EWADXX", "C1-516_2013-09-21_2013-10-03", 2,
               "Solexa-183228" );
          chrlist.push_back( "\"gi|240254421|ref|NC_003070.9|\"" );
          chrlist.push_back( "\"gi|240254678|ref|NC_003071.7|\"" );
          chrlist.push_back( "\"gi|240255695|ref|NC_003074.8|\"" );
          chrlist.push_back( "\"gi|240256243|ref|NC_003075.7|\"" );
          chrlist.push_back( "\"gi|240256493|ref|NC_003076.8|\"" );
          chrlist.push_back( "\"gi|26556996|ref|NC_001284.2|\"" );
          chrlist.push_back( "\"gi|7525012|ref|NC_000932.1|\"" );    }
  
     // Set up for ecoli11.

     if ( SAMPLE == "ecoli11" )
     {    AddTo( True, specs, "A42CA", "C1-516_2013-05-28_2013-05-30", 1,
               "Ec_11_9941_PCR_Free" );
          for ( int j = 1; j <= 15; j++ )
               chrlist.push_back( "Supercontig_6." + ToString(j) );    }

     // Set up for scardovia.

     if ( SAMPLE == "scardovia" )
     {    specs.push( True, "", "", -1, "", "/seq/picard/A6FWA/"
              "C1-516_2013-12-20_2013-12-22//1/Solexa-208793/A6FWA.1.unmapped.bam" );
          // AddTo( True, specs, "A6FWA", "C1-516_2013-12-20_2013-12-22", 1,
          //      "Solexa-208793" );
          // AddTo( True, specs, "A6VP5", "C1-516_2013-12-20_2013-12-22", 1,
          //      "Solexa-208793" ); // extra coverage, would need to tweak
          SELECT_FRAC = 0.2;    }

     // Set up for tb148.

     if ( SAMPLE == "tb148" )
     {    AddTo( True, specs, "A6FWA", "C1-516_2013-12-20_2013-12-22", 1,
               "Solexa-208791" );
          // AddTo( True, specs, "A6VP5", "C1-516_2013-12-20_2013-12-22", 1,
          //      "Solexa-208791" ); // extra coverage
          chrlist.push_back( "supercont1.1" );
          SELECT_FRAC = 0.4;    }

     // Set up for bifi.

     if ( SAMPLE == "bifi" )
     {    AddTo( True, specs, "A6FWA", "C1-516_2013-12-20_2013-12-22", 1,
               "Solexa-208788" );
          // AddTo( True, specs, "A6VP5", "C1-516_2013-12-20_2013-12-22", 1,
          //      "Solexa-208788" ); // extra coverage
          for ( int c = 1; c <= 9; c++ )
               chrlist.push_back( "supercont2." + ToString(c) );
          SELECT_FRAC = 0.2;    }

     // Set up for ecoli_scs.

     if ( SAMPLE == "ecoli_scs" )
     {    AddTo( True, specs, "A6FWA", "C1-516_2013-12-20_2013-12-22", 1,
               "Solexa-208786" );
          // AddTo( True, specs, "A6VP5", "C1-516_2013-12-20_2013-12-22", 1,
          //      "Solexa-208786" ); // extra coverage
          // AddTo( True, specs, "A6FWA", "C1-516_2013-12-20_2013-12-22", 1,
          //      "Solexa-208787" ); // extra coverage
          // AddTo( True, specs, "A6VP5", "C1-516_2013-12-20_2013-12-22", 1,
          //      "Solexa-208787" ); // extra coverage
          chrlist.push_back( "ecoli_K12_MG12655" );
          SELECT_FRAC = 0.4;    }

     // Set up for ecoli12.

     if ( SAMPLE == "ecoli12" )
     {    AddTo( True, specs, "A6FWA", "C1-516_2013-12-20_2013-12-22", 1,
               "Solexa-208789" );
          // AddTo( True, specs, "A6VP5", "C1-516_2013-12-20_2013-12-22", 1,
          //      "Solexa-208789" ); // extra coverage
          chrlist.push_back( "ecoli_K12_MG12655" );
          SELECT_FRAC = 0.4;    }

     // Set up for entero.

     if ( SAMPLE == "entero" )
     {    AddTo( True, specs, "A6FWA", "C1-516_2013-12-20_2013-12-22", 1,
               "Solexa-208792" );
          // AddTo( True, specs, "A6VP5", "C1-516_2013-12-20_2013-12-22", 1,
          //      "Solexa-208792" ); // extra coverage
          chrlist.push_back( "Contig2141", "Contig2126", "Contig498", 
               "Contig2166", "Contig200" );
          SELECT_FRAC = 0.4;    }

     // Set up for bcereus.

     if ( SAMPLE == "bcereus" )
     {    AddTo( True, specs, "A42CA", "C1-516_2013-05-28_2013-05-30", 1,
               "VD022_PCR_Free" );
          for ( int64_t j = 7000000184643749l; j <= 7000000184643779l; j++ )
               chrlist.push_back( ToString(j) );
          SELECT_FRAC = 0.25;    }

     // Set up for tb.

     if ( SAMPLE == "tb" )
     {    AddTo( True, specs, "A42CA", "C1-516_2013-05-28_2013-05-30", 1,
               "Haarlem_PCR_Free" );
          chrlist.push_back( "supercont2.1" );
          SELECT_FRAC = 0.35;    }

     // Expand DATASET and note some abbreviations.

     vec<String> dataset0, dataset;
     ParseStringSet( "{" + DATASET + "}", dataset0 );
     for ( int i = 0; i < dataset0.isize( ); i++ )
     {    if ( dataset0[i] == "F1" )
          {    dataset.push_back( "17E_PD", "16E_MD", "15E_DD" );    }
          else if ( dataset0[i] == "F2" )
          {    dataset.push_back( "23H_LM", "24H_CM", "25H_JM" );    }
          else if ( dataset0[i] == "F3" )
          {    dataset.push_back( "65T_CR", "66T_NG", "67T_SR", 
                    "68T_DR", "69T_GG" );    }
          else dataset.push_back( dataset0[i] );    }

     // Set up for plasmo.  Two datasets.  Note that the reference sequences are
     // not the same, so you can't mix them, and coordinates for dataset 1 may
     // be wrong.

     if ( SAMPLE == "plasmo" )
     {    if ( dataset.empty( ) ) dataset.push_back( "2" );
          // First dataset made in Chesterford April 2012, not PCR-free.
          auto dataset_seq=dataset;
          for( auto& current_set : dataset_seq){
              if ( current_set=="1" && IfDel( dataset, "1" ) )
              {    AddTo( False, specs, "CHESTERFO", "NORUNNUMBER", 1, "Solexa-90948" );
                   chrlist.push_back( "Pf3D7_01", "Pf3D7_02", "Pf3D7_03", "Pf3D7_04" );
                   chrlist.push_back( "Pf3D7_05", "Pf3D7_06", "Pf3D7_07", "Pf3D7_08" );
                   chrlist.push_back( "Pf3D7_09", "Pf3D7_10", "Pf3D7_11", "Pf3D7_12" );
                   chrlist.push_back( "Pf3D7_13", "Pf3D7_14", "M76611",
                        "PFC10_API_IRAB" );    }
              if ( current_set=="2" && IfDel( dataset, "2" ) )
              {    AddTo( True, specs, "A2RD9", "C1-516_2013-05-28_2013-05-30", 1,
                        "Solexa-160369" );
                   chrlist.push_back( "M76611", "PFC10_API_IRAB", "Pf3D7_01_v3" );
                   chrlist.push_back( "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3" );
                   chrlist.push_back( "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_07_v3" );
                   chrlist.push_back( "Pf3D7_08_v3", "Pf3D7_09_v3", "Pf3D7_10_v3" );
                   chrlist.push_back( "Pf3D7_11_v3", "Pf3D7_12_v3", "Pf3D7_13_v3" );
                   chrlist.push_back( "Pf3D7_14_v3" );    }
          }
     }

     // Set up for rhino.

     if ( SAMPLE == "rhino" )
     {    specs.push( True, "", "", -1, "", "/sysman/picard/H77KJADXX/"
               "C1-516_2014-02-12_2014-02-17/1/Solexa-216256/"
               "H77KJADXX.1.aligned.duplicates_marked.bam" );
          specs.push( True, "", "", -1, "", "/sysman/picard/H77KJADXX/"
               "C1-516_2014-02-12_2014-02-17/2/Solexa-216256/"
               "H77KJADXX.2.aligned.duplicates_marked.bam" );
          for ( int64_t j = 1; j <= 3086; j++ )
          {    String id = ToString(j);
               while( id.size( ) < 5 ) id = "0" + id;
               chrlist.push_back( "scaffold" + id );    }     }

     // Set up for human.

     if ( SAMPLE == "human" || SAMPLE == "human.hpool2" || SAMPLE == "human.hpool3" )
     {    if ( dataset.empty( ) ) dataset.push_back( "1" );
          auto dataset_seq=dataset;
          for( auto& current_set : dataset_seq){
              if ( current_set=="0" 
                    && IfDel( dataset, String("0") ) ) // NA12878, PCR library
              {    AddTo( False, specs, "H01RHADXX", "C1-508_2012-10-24_2012-10-31",
                        {1,2}, "Solexa-126138" );    }

              if ( current_set=="1" 
                    && IfDel( dataset, String("1") ) ) // NA12878, PCR-free library
              {    if ( !spec.NEW_Q2 && !spec.NEW_Q2_NEW_ALIGN )
                   {    AddTo( False, specs, "H01UJADXX", 
                             "C1-508_2012-11-01_2012-11-04",
                             {1,2}, "Solexa-125532" );    }
                   else if (spec.NEW_Q2)
                   {    specs.push( True, "", "", -1, "", "/seq/tng/datasets/"
                             "20130418-H01UJADXX-no-eamss/H01UJADXX.1/"
                             "NA12878/NA12878.bam" );
                        specs.push( True, "", "", -1, "", "/seq/tng/datasets/"
                             "20130418-H01UJADXX-no-eamss/H01UJADXX.2/"
                             "NA12878/NA12878.bam" );    }
                   else if (spec.NEW_Q2_NEW_ALIGN)
                   {    specs.push( True, "", "", -1, "",
                             "/fg/reich3/heng/tmp/NA12878-new/NA12878-new.bam" );   
                              }    }

              if ( current_set=="1s" && IfDel( dataset, "1s" ) ) // NA12878, PCR, 100 base reads, 180 bp frags
              {    AddTo( True, specs, "202PBABXX", "C1-202_2010-06-26_2013-05-13",
                        {1,2,3,4,5,6,7,8}, "Solexa-23661" );
                   AddTo( True, specs, "61PHDAAXX", "C1-202_2010-07-01_2013-05-13",
                        {6,7,8}, "Solexa-23661" );    }

              if ( current_set=="1.1" && IfDel( dataset, String("1.1") ) )
              {    if ( !spec.NEW_Q2 )
                   {    AddTo( False, specs, "H01UJADXX", 
                             "C1-508_2012-11-01_2012-11-04",
                             1, "Solexa-125532" );    }
                   else
                   {    specs.push( True, "", "", -1, "", "/seq/tng/datasets/"
                              "20130418-H01UJADXX-no-eamss/H01UJADXX.1/"
                              "NA12878/NA12878.bam" );    }    }
              if ( current_set=="1.2" && IfDel( dataset, String("1.2") ) )
              {    if ( !spec.NEW_Q2 )
                   {    AddTo( False, specs, "H01UJADXX",
                              "C1-508_2012-11-01_2012-11-04", 2, 
                              "Solexa-125532" );    }
                   else
                   {    specs.push( True, "", "", -1, "", "/seq/tng/datasets/"
                             "20130418-H01UJADXX-no-eamss/H01UJADXX.2/"
                             "NA12878/NA12878.bam" );    }    }
              if ( current_set=="2" && IfDel( dataset, 
                   String("2") ) ) // NA12878, another PCR-free library
              {    AddTo( False, specs, "H06HDADXX", "C1-508_2013-01-10_2013-01-13",
                        {1,2}, "Solexa-135852" );
                   AddTo( False, specs, "H06JUADXX", 
                        "C1-508_2013-01-10_2013-01-13", 1,
                        "Solexa-135852" );    }
              if ( current_set=="2.1" && IfDel( dataset, String("2.1") ) )
              {    AddTo( False, specs, "H06HDADXX", 
                        "C1-508_2013-01-10_2013-01-13", 1,
                        "Solexa-135852" );    }
              if ( current_set=="2.2" && IfDel( dataset, String("2.2") ) )
              {    AddTo( False, specs, "H06HDADXX", 
                        "C1-508_2013-01-10_2013-01-13", 2,
                        "Solexa-135852" );    }
              if ( current_set=="2.3" && IfDel( dataset, 
                   String("2.3") ) ) // NA12878, another PCR-free lib
              {    AddTo( False, specs, "H06JUADXX", 
                        "C1-508_2013-01-10_2013-01-13", 1,
                        "Solexa-135852" );    }
              if ( current_set=="3" 
                   && IfDel( dataset, String("3") ) ) // NA12882, data from Illumina
              {    String dir400 
                         = "/seq/tng/datasets/20130201-illumina-cambridge-2x400";
                   specs.push( True, "", "", -1, "",
                         dir400 + "/NA12882/NA12882.bam" );    }
              if ( current_set=="4" 
                   && IfDel( dataset, String("4") ) ) // NA12892, PCR-free
              {    AddTo( False, specs, "H06JHADXX", "C1-508_2013-01-10_2013-01-13",
                        {1,2}, "Solexa-135853" );
                   AddTo( False, specs, "H06JUADXX", 
                        "C1-508_2013-01-10_2013-01-13", 2,
                        "Solexa-135853" );    }

              if ( current_set == "1x"
                   && IfDel( dataset, String("1x") ) ) // NA12878, PCR-free improved

              {    String iroot = "/wga/scr4/vendor/illumina/2014-07-04/Conversion";
                   String abam = ".aligned.sorted.bam";
                   specs.push( True, "", "", -1, "", iroot +
                       "/140613_HSQ1185_0757_BH9TEUADXX/1/L1_NA12878/H9TEUADXX.1" 
                       + abam );
                   specs.push( True, "", "", -1, "", iroot +
                       "/140613_HSQ1185_0757_BH9TEUADXX/2/L2_NA12878/H9TEUADXX.2"
                       + abam );    }

              if ( current_set=="5" 
                   && IfDel( dataset, String("5") ) ) // NA12891, PCR-free
              {    AddTo( False, specs, "H03N7ADXX", "C1-508_2013-01-07_2013-01-10",
                        {1,2}, "Solexa-135851" );
                   AddTo( False, specs, "H05F1ADXX", 
                        "C1-508_2013-01-15_2013-01-18", 2, "Solexa-135851" );    }
              if ( current_set=="5.1" 
                   && IfDel( dataset, String("5.1") ) )
              {    AddTo( False, specs, "H03N7ADXX", "C1-508_2013-01-07_2013-01-10",
                        1, "Solexa-135851" );    }
              if ( current_set=="5.2" 
                   && IfDel( dataset, String("5.2") ) )
              {    AddTo( False, specs, "H03N7ADXX", "C1-508_2013-01-07_2013-01-10",
                        2, "Solexa-135851" );    }
              if ( current_set=="5.3" 
                   && IfDel( dataset, String("5.3") ) ) // NA12891, PCR-free
              {    AddTo( False, specs, "H05F1ADXX", 
                        "C1-508_2013-01-15_2013-01-18", 2, "Solexa-135851" );    }
          }
          for ( int c = 1; c <= 22; c++ )
               chrlist.push_back( ToString(c) );
          chrlist.push_back( "X", "Y", "MT" );    }

     // Define region for WGS datasets.

     if ( SAMPLE == "rhody" || SAMPLE == "ecoli11" || SAMPLE == "bcereus" 
          || SAMPLE == "tb" || SAMPLE == "plasmo" || SAMPLE == "human" 
          || SAMPLE == "arabidopsis" || SAMPLE == "entero" || SAMPLE == "bifi" 
          || SAMPLE == "tb148" || SAMPLE == "ecoli12" || SAMPLE == "rhino" )
     {    for ( int i = 0; i < regions.isize( ); i++ )
          {    if ( i > 0 ) region += " ";
               if ( !regions[i].Contains( ":" ) ) 
               {    if ( SAMPLE == "human" ) region += regions[i];
                    else region += chrlist[ regions[i].Int( ) ];    }
               else
               {    Bool bad = False;
                    String chr = regions[i].Before( ":" ); 
                    String range = regions[i].After( ":" );
                    if ( !range.Contains( "-" ) ) bad = True;
                    if ( !logc.TREAT_AS_UNKNOWN )
                    {    if ( SAMPLE == "human" )
                         {    if ( !Member( chrlist, chr ) ) bad = True;    }
                         else if ( !chr.IsInt( ) ) bad = True;    }
                    String start = range.Before( "-" ), stop = range.After( "-" );
                    if ( !start.IsInt( ) || !stop.IsInt( ) ) bad = True;
                    if (bad) FAIL_MSG( "Failed to parse X argument." );
                    if ( SAMPLE != "human" && !logc.TREAT_AS_UNKNOWN ) 
                         chr = chrlist[ chr.Int( ) ];
                    region += chr + ":" + ToString( start.Int( ) ) + "-" 
                         + ToString( stop.Int( ) );    }    }    }

     // Handle hpool1, hpool2 and hpool3.

     if ( SAMPLE == "hpool1" || SAMPLE == "hpool2" || SAMPLE == "hpool3" )
     {    if ( SAMPLE == "hpool1" )
          {    const int version = 2;
               if ( version == 1 )
               {    AddTo( False, specs, "A1PEB", "C1-708_2012-10-05_2012-10-08", 1,
                         { "Solexa-122851", "Solexa-122852" } );
                    if ( SELECT_FRAC < 0 && COVERAGE_CAP < 0 ) 
                         SELECT_FRAC = 0.15;    }
               if ( version == 2 )
               {    AddTo( False, specs, "A1MUT", "C1-508_2012-10-24_2012-10-26", 1, 
                         { "Solexa-122851", "Solexa-122852" } );
                    if ( SELECT_FRAC < 0 && COVERAGE_CAP < 0 ) 
                         SELECT_FRAC = 0.125;    }    }
          if ( SAMPLE == "hpool2" )
          {    if ( !HPOOL_ORIG )
               {    specs.push( True, "", "", -1, "", 
                         "/wga/scr4/human_data/NA12878_A2925/"
                         "Solexa-127359_aligned.sorted.bam" );    }
               else
               {    AddTo( False, specs, "A2925", "C1-508_2012-11-12_2012-11-14", 1,
                         "Solexa-127359" );    }    }
          if ( SAMPLE == "hpool3" )
          {    if ( !HPOOL_ORIG )
               {    specs.push( True, "", "", -1, "", 
                         "/wga/scr4/human_data/NA12878_A2925/"
                         "Solexa-127365_aligned.sorted.bam" );    }
               else
               {    AddTo( False, specs, "A2925", "C1-508_2012-11-12_2012-11-14", 1,
                         "Solexa-127365" );    }    }    }
     if ( SAMPLE == "hpool1" || SAMPLE == "hpool2" || SAMPLE == "human.hpool2"
          || SAMPLE == "hpool3" || SAMPLE == "human.hpool3" )
     {    vec<String> regionsx;
          if ( SAMPLE != "hpool1" )
          {    vec< vec< pair<String,String> > > junctions, breaks, edits; 
               ParseFosmidPoolMetainfo( regionsx, junctions, breaks, edits );    }
          else // SAMPLE = hpool1
          {    fast_ifstream in( "/wga/dev/references/Homo_sapiens/"
                    "WIBR_Fosmid_Pool.regions" );
               String line;
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    int start = line.Between( ":", "-" ).Int( );
                    int stop = line.After( "-" ).Int( );
                    start += Fosmid_trim_back, stop -= Fosmid_trim_back;
                    regionsx.push_back( line.Before( ":" ) + ":" 
                         + ToString(start) + "-" + ToString(stop) );    }    }
          vec<int> xs;
          for ( int j = 0; j < regions.isize( ); j++ )
          {    if ( !regions[j].IsInt( ) )
               {    PRINT2( j, regions[j] );
                    FatalErr("For now, X has to be a ParseIntSet-style list of "
                              " integers.\nAbort.");    }
               xs.push_back( regions[j].Int( ) );    }
          for ( int j = 0; j < xs.isize( ); j++ )
          {    if ( j > 0 ) region += " ";
               if ( SAMPLE == "human.hpool2" || SAMPLE == "human.hpool3" )
               {    String r = regionsx[ xs[j] ];
                    int start = r.Between( ":", "-" ).Int( );
                    int stop = r.After( "-" ).Int( );
                    start -= ext_hh, stop += ext_hh;
                    region += r.Before( ":" ) + ":" 
                         + ToString(start) + "-" + ToString(stop);    }
               else region += regionsx[ xs[j] ];    }    
          if ( xs.solo( ) ) 
               cout << Date( ) << ": using region " << region << endl;    }


     String CEPH = "/wga/scr4/human_data/CEPH";
     auto dataset_seq=dataset;
     for( auto& current_set : dataset_seq){
     // Handle Illumina/Platinum datasets for CEPH/UTAH pedigree 1463.
         if ( current_set=="P_momsmom" && IfDel( dataset, "P_momsmom" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12892_S1.bam" );
         if ( current_set=="P_momsdad" && IfDel( dataset, "P_momsdad" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12890_S1.bam" );
         if ( current_set=="P_dadsmom" && IfDel( dataset, "P_dadsmom" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12891_S1.bam" );
         if ( current_set=="P_dadsdad" && IfDel( dataset, "P_dadsdad" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12889_S1.bam" );
         if ( current_set=="P_mom" && IfDel( dataset, "P_mom" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12878_S1.bam" );
         if ( current_set=="P_dad" && IfDel( dataset, "P_dad" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12877_S1.bam" );
         if ( current_set=="P_d1" && IfDel( dataset, "P_d1" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12879_S1.bam" );
         if ( current_set=="P_d2" && IfDel( dataset, "P_d2" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12880_S1.bam" );
         if ( current_set=="P_d3" && IfDel( dataset, "P_d3" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12881_S1.bam" );
         if ( current_set=="P_d4" && IfDel( dataset, "P_d4" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12885_S1.bam" );
         if ( current_set=="P_d5" && IfDel( dataset, "P_d5" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12887_S1.bam" );
         if ( current_set=="P_s1" && IfDel( dataset, "P_s1" ) )
         {    specs.push( True, "", "", -1, "", CEPH + "/NA12882_S1.bam" );
              specs.push( True, "", "", -1, "", CEPH + "/NA12882_2_S1.bam" );    }
         if ( current_set=="P_s2" && IfDel( dataset, "P_s2" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12883_S1.bam" );
         if ( current_set=="P_s3" && IfDel( dataset, "P_s3" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12884_S1.bam" );
         if ( current_set=="P_s4" && IfDel( dataset, "P_s4" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12886_S1.bam" );
         if ( current_set=="P_s5" && IfDel( dataset, "P_s5" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12888_S1.bam" );
         if ( current_set=="P_s6" && IfDel( dataset, "P_s6" ) )
              specs.push( True, "", "", -1, "", CEPH + "/NA12893_S1.bam" );

     // Handle cancer.

         if ( current_set=="HCC1954" && IfDel( dataset, "HCC1954" ) )
         {    AddTo( True, specs, "H7AH6ADXX", "C1-516_2013-09-05_2013-09-08", {1,2},
                   "Solexa-178364" );
              AddTo( True, specs, "H7AH7ADXX", "C1-516_2013-09-05_2013-09-08", {1,2},
                   "Solexa-178364" );    }

     // Handle MacArthur muscle disease samples.

         if ( current_set=="17E_PD" && IfDel( dataset, "17E_PD" ) )
         {    AddTo( True, specs, "H09P1ADXX", "C1-508_2013-04-12_2014-05-23", {1,2},
                   "Solexa-151576" );
              AddTo( True, specs, "H09PGADXX", "C1-508_2013-04-12_2014-05-23", {1,2},
                   "Solexa-151576" );    }
         if ( current_set=="23H_LM" && IfDel( dataset, "23H_LM" ) )
         {    AddTo( True, specs, "H09PDADXX", "C1-508_2013-04-12_2014-05-23", 2,
                   "Solexa-151574" );
              AddTo( True, specs, "H09NJADXX", "C1-508_2013-04-12_2014-05-23", 1,
                   "Solexa-151574" );
              AddTo( True, specs, "H0KTMADXX", "C1-508_2013-05-17_2014-05-23", 1,
                   "Solexa-151574" );    }
         if ( current_set=="65T_CR" && IfDel( dataset, "65T_CR" ) )
         {    AddTo( True, specs, "H09P8ADXX", "C1-508_2013-04-12_2014-04-24", 2,
                   "Solexa-151577" );
              AddTo( True, specs, "H09RDADXX", "C1-508_2013-04-12_2014-04-24", {1,2},
                   "Solexa-151577" );
              AddTo( True, specs, "H0KTMADXX", "C1-508_2013-05-17_2014-05-23", 2,
                   "Solexa-151577" );    }
         if ( current_set=="66T_NG" && IfDel( dataset, "66T_NG" ) )
         {    AddTo( True, specs, "H09NKADXX", "C1-508_2013-04-12_2014-05-23", 1,
                   "Solexa-151580" );
              AddTo( True, specs, "H0A1RADXX", "C1-508_2013-04-12_2014-04-24", {1,2},
                   "Solexa-151580" );    }
         if ( current_set=="68T_DR" && IfDel( dataset, "68T_DR" ) )
         {    AddTo( True, specs, "H09H9ADXX", "C1-508_2013-04-12_2014-04-24", {1,2},
                   "Solexa-151581" );    }
         if ( current_set=="69T_GG" && IfDel( dataset, "69T_GG" ) )
         {    AddTo( True, specs, "H09HFADXX", "C1-508_2013-04-11_2014-04-24", {1,2},
                   "Solexa-151579" );
              AddTo( True, specs, "H0KREADXX", "C1-508_2013-05-17_2014-05-23", 2,
                   "Solexa-151579" );    }
         if ( current_set=="67T_SR" && IfDel( dataset, "67T_SR" ) )
         {    AddTo( True, specs, "H09NJADXX", "C1-508_2013-04-12_2014-05-23", 2,
                   "Solexa-151582" );
              AddTo( True, specs, "H09P8ADXX", "C1-508_2013-04-12_2014-04-24", 1,
                   "Solexa-151582" );    }
         if ( current_set=="16E_MD" && IfDel( dataset, "16E_MD" ) )
         {    AddTo( True, specs, "H09NKADXX", "C1-508_2013-04-12_2014-05-23", 2,
                   "Solexa-151573" );
              AddTo( True, specs, "H09NHADXX", "C1-508_2013-04-12_2014-05-23", 1,
                   "Solexa-151573" );    }
         if ( current_set=="24H_CM" && IfDel( dataset, "24H_CM" ) )
         {    AddTo( True, specs, "H09PDADXX", "C1-508_2013-04-12_2014-05-23", 1,
                   "Solexa-151583" );
              AddTo( True, specs, "H09PJADXX", "C1-508_2013-04-12_2014-05-23", 2,
                   "Solexa-151583" );
              AddTo( True, specs, "H0KREADXX", "C1-508_2013-05-17_2014-05-23", 1,
                   "Solexa-151583" );    }
         if ( current_set=="15E_DD" && IfDel( dataset, "15E_DD" ) )
         {    AddTo( True, specs, "H09PJADXX", "C1-508_2013-04-12_2014-05-23", 1,
                   "Solexa-151575" );
              AddTo( True, specs, "H09NHADXX", "C1-508_2013-04-12_2014-05-23", 2,
                   "Solexa-151575" );    }
         if ( current_set=="25H_JM" && IfDel( dataset, "25H_JM" ) )
         {    AddTo( True, specs, "H09PHADXX", "C1-508_2013-04-11_2014-05-23", {1,2},
                   "Solexa-151578" );
              AddTo( True, specs, "H0JG5ADXX", "C1-508_2013-05-17_2014-05-23", {1,2},
                   "Solexa-151578" );    }
     }

     // Check for illegal datasets.

     if ( dataset.nonempty( ) )
     {    cout << "\nFound illegal DATASET(s): " << printSeq(dataset) << "\n" 
               << endl;
          _exit(1);    }    }

void Dexterize( const String& var, String& value, String& SAMPLE, String& X, 
     double& GENOME_SUB_PERCENT, String& TMP, String& READS )
{    Bool illegal = False;
     if ( !value.Contains( ":" ) ) illegal = True;
     if ( !illegal )
     {    if ( !value.After( ":" ).Contains( ":" ) ) illegal = True;    }
     if (illegal)
     {    FatalErr("Illegal use of #dexter for " << var << ": " << value);    }
     String sample = value.Between( ":", ":" ).ToLower( );
     String date = value.After( ":" ).After( ":" );
     String dir1 = "/wga/scr4/dexter/stage/" + date + "/longread";
     if ( !IsDirectory(dir1) )
     {    FatalErr(var << ": can't find directory " << dir1);    }
     vec<String> all = AllFiles(dir1), matches;
     for ( int j = 0; j < all.isize( ); j++ )
     {    String d = all[j];
          if ( d.ToLower( ).Contains(sample) ) matches.push_back( all[j] );    }
     if ( !matches.solo( ) )
     {    FatalErr(var << ": found " << matches.size( ) << " matching subdirs");   }
     String dir2 = dir1 + "/" + matches[0];
     if ( var == "IN_SHBV" ) value = dir2 + "/assembly.1.shbv";
     if ( var == "IN_SHBV_FINAL" ) value = dir2 + "/assembly.final.shbv";
     if ( var == "IN_EFASTA_READS" ) value = dir2 + "/corrected.efasta";
     if ( IsDirectory( dir2 + "/tmp" ) )
     {    TMP = dir2 + "/tmp";
          READS = "#picard";    }

     // Parse LongProto command.  This does not correctly handle backslashes, so
     // hopefully nothing important got slashed.

     fast_ifstream in( dir2 + "/statistics0.txt" );
     String line, command;
     int BREAKS = 0;
     Bool found = False;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( "---------------" ) )
          {    if ( ++BREAKS == 2 ) break;
               continue;    }
          if ( BREAKS < 1 ) continue;
          if ( line.Contains( "LongProto", 0 ) ) found = True;
          if ( BREAKS == 1 && !found ) continue;
          command += line;    }
     vec<String> bits;
     Tokenize( command, bits );
     for ( int j = 0; j < bits.isize( ); j++ )
     {    if ( bits[j].Contains( "SAMPLE=" ) ) SAMPLE = bits[j].After( "SAMPLE=" );
          if ( bits[j].Contains( "X=" ) ) 
          {    X = bits[j].After( "X=" );
               if ( X.nonempty( ) && X.front( ) == '"' && X.back( ) == '"' )
                    X=X.Between( "\"", "\"" );    }
          if ( bits[j].Contains( "GENOME_SUB_PERCENT=" ) )
          {    GENOME_SUB_PERCENT = bits[j].After( 
                    "GENOME_SUB_PERCENT=" ).Double( );    }    }    }

void DexterizeAll( String& IN_SHBV, String& IN_SHBV_FINAL, String& IN_EFASTA_READS,
     const int START_STEP, String& SAMPLE, String& X, double& GENOME_SUB_PERCENT,
     String& TMP, String& READS )
{    if ( IN_SHBV.Contains( "#dexter", 0 ) ) 
     {    Dexterize( "IN_SHBV", IN_SHBV, SAMPLE, X, GENOME_SUB_PERCENT, TMP, READS );
          if ( START_STEP > 1 )
          {    IN_SHBV = IN_SHBV.RevBefore( ".shbv" ).RevBefore( "." )
                    + "." + ToString( START_STEP - 1 ) + ".shbv";    }    }
     if ( IN_SHBV_FINAL.Contains( "#dexter", 0 ) ) 
     {    Dexterize( "IN_SHBV_FINAL", IN_SHBV_FINAL, SAMPLE, X, GENOME_SUB_PERCENT,
               TMP, READS );    }
     if ( IN_EFASTA_READS.Contains( "#dexter", 0 ) ) 
     {    Dexterize( "IN_EFASTA_READS", IN_EFASTA_READS, SAMPLE, 
               X, GENOME_SUB_PERCENT, TMP, READS );    }    }

