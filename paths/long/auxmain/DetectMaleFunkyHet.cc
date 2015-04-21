///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Flag .variant.vcf lines as funky (and print) if any of:
// - non-Mendelian
// - male has two alleles on X or Y
// - neither allele carried by some individual.
//
// Probabilities < 0.01 are rounded down to zero, probabilities > 0.995 are rounded
// up to 1, and for intermediate probabilities, we consider all possibilities.
//
// Problems:
// - hardwired for F3.

#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "math/Functions.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(IN_HEAD, "assembly head");
     EndCommandArguments;

     // Make deductions from X.

     // Start reading the vcf.

     fast_ifstream in( IN_HEAD + ".variant.filtered.vcf" );
     String line;
     int print_count = 0;
     uint64_t count1=0;
     while(1)
     {    
          // Parse a line.

          getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( "#", 0 ) ) continue;
          istringstream iline( line.c_str( ) );
          String chr, junk, ref, alt, Q, defs,filter,gen;
          int pos;
          iline >> chr >> pos >> junk >> ref >> alt >> Q >> filter >> junk >> junk >> gen;
          
          Bool sex = ( chr == "X" || chr == "Y" );
          Bool pseudoautosomal = ( chr == "X" && ( pos < 2709520 || ( pos >= 154584237 && pos < 154913754 ) ) );
          
          if( !filter.Contains("PASS") || !sex || pseudoautosomal ){
              continue;
          }
          
          vec<String> fields, q;
          Tokenize( gen, {':'}, fields );
          Tokenize( fields[3], {','}, q );
          q.push_front( fields[2] );
          
          if ( gen.Contains( ":" ) ) gen = gen.Before( ":" );
          gen.GlobalReplaceBy( "|", "/" );
          ForceAssert(gen.Contains("/"));
          
          int hi=gen.Before("/").Int();
          int hi2=gen.After("/").Int();
          bool good = hi==hi2;
          
          
          vec<vec<String>> rows;
          vec<String> row;
          row.push_back( "", "ref" );
          
          if (good) continue;
          
          vec<String> alts;
          Tokenize( alt, {','}, alts );

          // OK found a bad event, print it.

          if ( print_count == 0 ) cout << "\n";
          else
          {    cout << "----------------------------------------------------------"
                    << "-----------------------\n";    }
          cout << "[" << ++print_count << "] " << chr << ":" << pos << endl;
          vec<String> lefts;
          lefts.push_back( "ref  = " );
          for ( int j = 0; j < alts.isize( ); j++ )
               lefts.push_back( "alt" + ToString(j+1) + " = " );
          vec<String> rights;
          rights.push_back(ref);
          rights.append(alts);
          for ( int i = 0; i < lefts.isize( ); i++ )
          {    cout << lefts[i];
               for ( int j = 0; j < rights[i].isize( ); j++ )
               {    if ( j - 74 % 80 == 0 ) cout << "\n       ";
                    cout << rights[i][j];    }
               cout << "\n";    }
          for ( int j = 1; j <= alts.isize( ); j++ )
               row.push_back( "alt" + ToString(j) );
          rows.push_back(row);
          for ( int i = 0; i < 1; i++ )
          {    vec<String> row;
               row.push_back( "M" );
               for ( int j = 0; j < q.isize( ); j++ )
               {
                   row.push_back(q[j]);
               }
               rows.push_back(row);
          }
          PrintTabular( cout, rows, 2, "l" + String( 1, 'l' ) );    }
}
