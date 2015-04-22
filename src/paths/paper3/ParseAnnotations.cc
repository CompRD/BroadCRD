///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"

int main( )
{    
     RunTime( );

     for ( int pass = 1; pass <= 2; pass++ )
     {    
          String dir = ( pass == 1 ? "/wga/scr1/ALLPATHS/M.musculus"
               : "/wga/scr1/ALLPATHS/H.sapiens.NA12878" );
          String genome = ( pass == 1 ? "genome" 
               : "other_ref/build19.mito.female.nonrandom" );
          String build = ( pass == 1 ? "mm9" : "hg19" );
          fast_ifstream gin( dir + "/" + genome + ".fasta" );
          Ofstream( out, dir + "/annotation/annotations" );
          String line;
          vec<String> ids;
          while(1)
          {    getline( gin, line );
               if ( gin.fail( ) ) break;
               if ( !line.Contains( ">", 0 ) ) continue;
               if ( line.Contains( " " ) ) line = line.Before( " " );
               PRINT(line); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               ids.push_back( line.After( ">" ) );    }
          fast_ifstream rin( dir + "/annotation/" + build + ".RepeatMasker" );
          vec<String> tokens;
          getline( rin, line );
          vec<char> tab, comma;
          tab.push_back( '\t' ), comma.push_back( ',' );
          TokenizeStrictly( line, tab, tokens );
          int chrom_pos = Position( tokens, String("genoName") );
          ForceAssert( chrom_pos >= 0 );
          int start_pos = Position( tokens, String("genoStart") );
          int end_pos = Position( tokens, String("genoEnd") );
          ForceAssert( start_pos >= 0 ), ForceAssert( end_pos >= 0 );
          int class_pos = Position( tokens, String("repClass") );
          ForceAssert( class_pos >= 0 );
          while(1)
          {    getline( rin, line );
               if ( rin.fail( ) ) break;
               TokenizeStrictly( line, tab, tokens );
               String chrom = tokens[chrom_pos];
               if ( pass == 2 && chrom.Contains( "chr", 0 ) ) 
                    chrom = chrom.After( "chr" );
               if ( chrom == "M" ) chrom = "MT";
               int chromp = Position( ids, chrom );
               if ( chrom == "Y" && chromp < 0 ) continue;
               if ( chrom.Contains( "_" ) && chromp < 0 ) continue;
               if ( chromp < 0 ) cout << "Can't find " << chrom << "." << endl;
               ForceAssertGe( chromp, 0 );
               int start = tokens[start_pos].Int( ) - 1;
               if ( start == -1 ) start = 0; // not sure how this can happen
               int end = tokens[end_pos].Int( );
               String Class = tokens[class_pos];
               out << chromp << " " << start << " " << end << " " 
                    << Class << "\n";    }
          fast_ifstream in( dir + "/annotation/" + build + ".RefSeqGenes" );
          getline( in, line );
          TokenizeStrictly( line, tab, tokens );
          chrom_pos = Position( tokens, String("chrom") );
          ForceAssert( chrom_pos >= 0 );
          int starts_pos = Position( tokens, String("exonStarts") );
          int ends_pos = Position( tokens, String("exonEnds") );
          ForceAssert( starts_pos >= 0 ), ForceAssert( ends_pos >= 0 );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               TokenizeStrictly( line, tab, tokens );
               String chrom = tokens[chrom_pos];
               if ( pass == 2 && chrom.Contains( "chr", 0 ) ) 
                    chrom = chrom.After( "chr" );
               if ( chrom == "M" ) chrom = "MT";
               int chromp = Position( ids, chrom );
               if ( chrom == "Y" && chromp < 0 ) continue;
               if ( chrom.Contains( "_" ) && chromp < 0 ) continue;
               if ( chromp < 0 ) cout << "Can't find " << chrom << "." << endl;
               ForceAssertGe( chromp, 0 );
               vec<String> starts, ends;
               TokenizeStrictly( tokens[starts_pos], comma, starts );
               TokenizeStrictly( tokens[ends_pos], comma, ends );
               for ( int j = 0; j < starts.isize( ) - 1; j++ )
               {    out << chromp << " " << starts[j].Int( ) - 1 << " " 
                         << ends[j].Int( ) << " EXON" << endl;    }    }    }    }
