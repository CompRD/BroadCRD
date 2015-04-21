///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// VcfFrac.  Messy code to find variants outputed by TiTv that are not called by mom.

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "STLExtensions.h"
#include "TokenizeString.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(VCF,
          "/wga/scr4/human_assemblies/1/v8/combined.filtered.vcf",
          "VCF file to use as input (currently must be produced by DISCOVAR)");
     CommandArgument_Bool_OrDefault_Doc(DETAILS, False, 
          "to print the variants");
     EndCommandArguments;

//     VCF = "/wga/scr4/jaffe/mom.vcf";

     // DISCOVAR constants.

//     const double cutoff1 = 0.995;
//     const double cutoff2 = 0.9;

     // Load DISCOVAR variants.

     fast_ifstream in(VCF);
     vec<String> lines;
     String line, chr, pos, ref, alt, gen, junk, filter;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( "#", 0 ) ) continue;
          if ( line[0] != 'X' && line[0] != 'Y' 
               && ( line[0] < '0' || line[0] > '9' ) )
          {     continue;     }
          lines.push_back(line);    }

     // Filter variants.

     vec< triple<String,String,String> > VARS;
     for ( int z = 0; z < lines.isize( ); z++ )
     {    const String& line = lines[z];
          istrstream iline( line.c_str( ) );
          iline >> chr >> pos >> junk >> ref >> alt >> junk >> filter 
               >> junk >> junk >> gen;
          if(!filter.Contains("PASS")) continue;
          vec<String> fields, q, alts;
          Tokenize( gen, {':'}, fields );
          Tokenize( fields[3], {','}, q );
          q.push_front( fields[2] );
          vec<int> alleles;
          /*
          Bool flaky = False;
          for ( int j = 0; j < q.isize( ); j++ )
          {    int m = Min( (int) q[j].Int( ), 40 );
               double p = 1.0 - pow( 10, -m/10.0 );
               if ( p > 0 && p < cutoff2 ) flaky = True;
               if ( p > cutoff1 ) alleles.push_back(j);    }
          if (flaky) continue;
          */
          
          if ( gen.Contains( ":" ) ) gen = gen.Before( ":" );
          gen.GlobalReplaceBy( "|", "/" );
          ForceAssert(gen.Contains("/"));
          {
              int hi=gen.Before("/").Int();
              int hi2=gen.After("/").Int();
              alleles.push_back(min(hi,hi2));
              if(hi!=hi2) alleles.push_back(max(hi,hi2));
          }
          Tokenize( alt, {','}, alts );
          if ( alleles.empty( ) || alleles.size( ) > 2 ) continue;
          if ( alleles.solo( ) && alleles[0] == 0 ) continue;
          if ( alleles.size( ) > 2 && alleles[0] != 0 ) continue;
          String R = ref, A;
          if ( alleles.solo( ) ) A = alts[ alleles[0] - 1 ];
          else A = alts[ alleles[1] - 1 ];
          if ( R.size( ) > 1 || A.size( ) > 1 ) continue;
          VARS.push( chr + ":" + pos, R, A );    }
     UniqueSort(VARS);

     fast_ifstream in2( "xxx.hom" );
     int total = 0, found = 0;
     while(1)
     {    getline( in2, line );
          if ( in2.fail( ) ) break;
          istringstream iline( line.c_str( ) );
          String junk, loc, ref, alt;
          iline >> junk >> junk >> loc >> junk >> ref >> junk >> alt;
          triple<String,String,String> var( loc, ref, alt );
          total++;
          if ( BinMember( VARS, var ) ) found++;
          else if (DETAILS) cout << line << endl;    }
     PRINT2( found, total );
     cout << PERCENT_RATIO( 3, found, total ) << endl;    }
