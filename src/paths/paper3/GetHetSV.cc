///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// GetHetSV: get list of heterozygous deletions from 1000G site, analyze how they
// appear in the assembly.

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "lookup/LookAlign.h"

int main( )
{    RunTime( );
     vecbasevector ref( 
          "/wga/scr1/ALLPATHS/H.sapiens.NA12878/other_ref/build18.fastb" );
     String dir = "/humgen/1kg/DCC/ftp/pilot_data/release/2010_07/low_coverage/sv";
     fast_pipe_ifstream in( "zcat " + dir 
          + "/CEU.low_coverage.2010_06.deletions.genotypes.vcf.gz "
          "| grep -v IMPRECISE" );
     String line;
     vec<char> seps, seps_semi;
     seps.push_back( '\t' );
     seps_semi.push_back( ';' );
     vecbasevector choice;
     vec<String> CHR;
     vec<int> START, LEN;
     vec<String> REPLACEMENT;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( "#", 0 ) ) continue;
          vec<String> fields;
          TokenizeStrictly( line, seps, fields );
          if ( !fields[51].Contains( "0/1", 0 ) ) continue;
          int score = fields[51].Between( ":", ":" ).Int( );
          if ( score < 40 ) continue;
          vec<String> field7;
          TokenizeStrictly( fields[7], seps_semi, field7 );
          String chr = fields[0];
          int chrint = chr.Int( );
          int start = fields[1].Int( ) - 1;
          int len = fields[3].size( );
          int svlen = 0;
          for ( int j = 0; j < field7.isize( ); j++ )
          {    if ( field7[j].Contains( "SVLEN" ) )
                    svlen = field7[j].After( "=" ).Int( );    }
          String replacement = fields[4];
          CHR.push_back(chr);
          START.push_back(start);
          LEN.push_back(len);
          REPLACEMENT.push_back(replacement);
          basevector c1, c2a, c2b, c2c, c2;
          int flank = 50;
          c1.SetToSubOf( ref[chrint], start - flank, len + 2*flank );
          c2a.SetToSubOf( ref[chrint], start - flank, flank );
          c2b.SetFromString(replacement);
          c2c.SetToSubOf( ref[chrint], start + len, flank );
          c2 = Cat( c2a, c2b, c2c );
          choice.push_back_reserve(c1);
          choice.push_back_reserve(c2);    
          // if ( CHR.size( ) == 100 ) break;
               }

     choice.WriteAll( "ccc.fastb" );

     SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS=ccc.fastb "
          "L=/wga/scr1/ALLPATHS/H.sapiens.NA12878/fullhuman/fullhuman_v1.1/iainm.22jul/ASSEMBLIES/test/linear_scaffolds0.contigs.lookup PARSEABLE=True "
          "VISUAL=True NH=True QUIET=True > ccc.aligns" );

     vec<look_align> aligns;
     vec< vec<int> > aligns_index;
     LoadLookAligns( "ccc.aligns", aligns, aligns_index, choice.size( ) );

     int shorts = 0, longs = 0, boths = 0, neithers = 0, others = 0;
     for ( int i = 0; i < CHR.isize( ); i++ )
     {    String chr = CHR[i];
          int start = START[i];
          int len = LEN[i];
          String replacement = REPLACEMENT[i];

          Bool perfect1 = False, imperfect1 = False;
          Bool perfect2 = False, imperfect2 = False;
          for ( int j = 0; j < aligns_index[2*i].isize( ); j++ )
          {    const look_align& la = aligns[ aligns_index[2*i][j] ];
               if ( la.ErrorRate( ) > 0.1 ) continue;
               if ( la.FullLength( ) && la.ErrorRate( ) <= 0.02 ) perfect1 = True;
               else imperfect1 = True;    }
          for ( int j = 0; j < aligns_index[2*i+1].isize( ); j++ )
          {    const look_align& la = aligns[ aligns_index[2*i+1][j] ];
               if ( la.ErrorRate( ) > 0.1 ) continue;
               if ( la.FullLength( ) && la.ErrorRate( ) <= 0.02 ) perfect2 = True;
               else imperfect2 = True;    }

          if ( perfect1 && !perfect2 && !imperfect2 ) longs++;
          else if ( !perfect1 && !imperfect1 && perfect2 ) shorts++;
          else if ( perfect1 && perfect2 ) boths++;
          else if ( !perfect1 && !imperfect1 && !perfect2 && !imperfect2 ) 
               neithers++;
          else
          {    others++;
               cout << "\n=========================================================="
                    << "==========================\n\n";
               PRINT4( chr, start, len, replacement );    
               for ( int j = 0; j < aligns_index[2*i].isize( ); j++ )
               {    cout << "\nalignment " << j+1 << " of first sequence" << endl;
                    const look_align& la = aligns[ aligns_index[2*i][j] ];
                    la.PrintReadableBrief(cout);    }
               for ( int j = 0; j < aligns_index[2*i+1].isize( ); j++ )
               {    cout << "\nalignment " << j+1 << " of second sequence" << endl;
                    const look_align& la = aligns[ aligns_index[2*i+1][j] ];
                    la.PrintReadableBrief(cout);    }    }    }
     cout << "\n";
     PRINT5( shorts, longs, boths, neithers, others );    }
