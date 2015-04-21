///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// GetExons.  Create fasta file of coding exons for a given gene.  Add flanks
// of given size if requested.

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "math/HoInterval.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault(GENE, "");
     CommandArgument_Bool_OrDefault_Doc(PREFIX, False, 
          "if True, gene name has to only start with GENE");
     CommandArgument_String_OrDefault(TRANSCRIPT, "");
     CommandArgument_Int_OrDefault(FLANK, 0);
     CommandArgument_String_OrDefault(OUT, "");
     CommandArgument_Int_OrDefault(START, -1);
     CommandArgument_Int_OrDefault(STOP, -1);
     CommandArgument_String_OrDefault_Doc(MODE, "fasta",
          "fasta or exons or intervals; overlapping intervals are merged");
     CommandArgument_Bool_OrDefault_Doc(CODING, True, "process only coding exons");
     EndCommandArguments;

     fast_ifstream in( "/seq/references/Homo_sapiens_assembly19/v1/"
          "annotation/gencode.v12.annotation.patched_contigs.gff" );

     ostream* out = ( OUT == "" ? &cout : new ofstream( OUT.c_str( ) ) );

     String line;
     vec<char> sep = { '\t' };
     vec<String> fields;
     vec< triple<String,int,int> > exons;
     vec<String> chrs;
     vec<int> ichrs;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( !PREFIX && !line.Contains( "gene_name \"" + GENE + "\"" ) ) 
               continue;
          if ( PREFIX && !line.Contains( "gene_name \"" + GENE ) ) continue;
          TokenizeStrictly( line, sep, fields );
          if ( fields[2] == "gene" || fields[2] == "transcript" ) continue;
          if ( CODING && fields[2] != "CDS" ) continue;
          if ( TRANSCRIPT != "" && !line.Contains( "transcript_name \"" + TRANSCRIPT
               + "\"" ) )
          {    continue;    }
          String chr;
          if ( fields[0] == "MT" ) continue; // IGNORING MITOCHONDRIAL GENES
          else chr = fields[0];
          int start = fields[3].Int( ) - 1, stop = fields[4].Int( );
          if ( START >= 0 && start < START ) continue;
          if ( STOP >= 0 && stop > STOP ) continue;
          chrs.push_back(chr);
          if ( chr == "X" ) ichrs.push_back(22);
          else if ( chr == "Y" ) ichrs.push_back(23);
          else ichrs.push_back( chr.Int( ) - 1 );
          exons.push( chr, start, stop );    }
     UniqueSortSync( chrs, ichrs );

     vecbasevector g;
     if ( MODE != "exons" ) g.Read( "/wga/scr4/bigrefs/human19/genome.fastb", ichrs );

     vec< vec<ho_interval> > intervals(24);
     for ( const auto& x : exons )
     {    int cp = BinPosition( chrs, x.first );
          int start = x.second, stop = x.third;
          String gt = GENE;
          if ( TRANSCRIPT != "" ) gt += "." + TRANSCRIPT;
          String title = gt + ".exon_at." + x.first + ":" + ToString(start) + "-" 
               + ToString(stop) + "_plus_flank_of_" + ToString(FLANK);
          start -= FLANK;
          stop += FLANK;
          if ( MODE == "fasta" )
          {    basevector b( g[cp], start, stop - start );
               b.Print( *out, title );    }
          else if ( MODE == "exons" ) cout << title << "\n";
          else if ( MODE == "intervals" ) 
          {    intervals[cp].push( start, stop );    }    }

     if ( MODE == "intervals" )
     {    int64_t total_intervals = 0, total_bases = 0;
          for ( int c = 0; c < 24; c++ )
          {    if ( intervals[c].empty( ) ) continue;
               int cp = BinPosition( ichrs, c );
               ExtractGivenCoverage( g[cp].size( ), 1, intervals[c], intervals[c] );
               for ( int j = 0; j < intervals[c].isize( ); j++ )
               {    cout << chrs[c] << ":" << intervals[c][j] << "\n";
                    total_intervals++;
                    total_bases += intervals[c][j].Length( );    }    }
          cout << "total intervals = " << ToStringAddCommas(total_intervals) << "\n";
          cout << "total bases = " << ToStringAddCommas(total_bases) 
               << "\n";    }    }
