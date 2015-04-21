///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FetchExonLocs.  Terrible hack.

#include "FetchReads.h"
#include "MainTools.h"
#include "pairwise_aligners/SmithWatFree.h"

int main(int argc, char *argv[])
{    
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int_Doc(ID, "id of Fosmid");
     EndCommandArguments;

     vecbasevector ref;
     vecString ref_names;
     FetchReads( ref, ref_names, 
          "/wga/dev/jaffe/BroadCRD/fos." + ToString(ID) + ".fasta" );

     String source = "/wga/scr4/blau/share/fosmid_filtering/output.fasta";
     vecbasevector exons_all;
     vecString enames_all;
     FetchReads( exons_all, enames_all, source );
     vec<basevector> exons;
     vec<String> genes;
     vec< pair<int,int> > exon_n;
     vec<int> start, stop;
     for ( int i = 0; i < (int) exons_all.size( ); i++ )
     {    if ( enames_all[i].Before( "." ).Int( ) != ID ) continue;
          int exon_no = enames_all[i].Between( " Exon=", " " ).Int( );
          int total_exons = enames_all[i].Between( " Total_Exon=", " " ).Int( );
          // if ( exon_no == 1 || exon_no == total_exons ) continue;
          exons.push_back( exons_all[i] );
          exon_n.push( exon_no, total_exons );
          String gene = enames_all[i].Between( ".", "." );
          genes.push_back(gene);
          start.push_back( enames_all[i].Between( "Start=", " " ).Int( ) );
          stop.push_back( enames_all[i].After( "End=" ).Int( ) );    } 
     SortSync( start, stop, exons, genes, exon_n );

     for ( int i = 0; i < (int) exons.size( ); i++ )
     {    alignment a;
          int best_loc;
          SmithWatFree( exons[i], ref[0], best_loc, a );
          cout << a.pos2( ) << "-" << a.Pos2( ) << ", gene = " << genes[i] << ", "
               << "exon " << exon_n[i].first << " of " 
               << exon_n[i].second << endl;    }    }
