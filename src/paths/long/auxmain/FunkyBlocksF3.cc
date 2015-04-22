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
#include "Basevector.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(IN_HEAD, "assembly head");
     CommandArgument_Double_OrDefault_Doc(HI_PROB, 0.995,
          "probability value above this value is treated as exaxtly 1");
     CommandArgument_Double_OrDefault_Doc(LO_PROB, 0.005,
          "probability value above this value is treated as exactly 0");
     CommandArgument_Int_OrDefault_Doc(VERBOSITY, 0, "level of verbosity");
     CommandArgument_Bool_OrDefault_Doc(DISCOVAR_FILTER, False, "ignore sample discrepancy that would've been filtered out by a certain generalization of single-sample filter");
     CommandArgument_Bool_OrDefault_Doc(INVOLVES_MALE_XHET, False, "report only funky blocks with >1 male xhet violation");
     CommandArgument_Bool_OrDefault_Doc(IGNORE_MALE_XHET, False, "ignore funky blocks with >1 male xhet violation");
     EndCommandArguments;
     
     const String X="F3";

     // Make deductions from X.

     int iQCeiling = max(40.0, -10.0*log10(1.-HI_PROB)+1.0) + 0.5;

     int p_total = ( X == "F3" ? 5 : 1 );
     int p_male = ( X == "F3" ? 2 : 1 );
     
     String fReference = "/wga/dev/references/Homo_sapiens/genome.fastb";
     vecbasevector genome;
     genome.ReadAll(fReference);

     // Start reading the vcf.

     fast_ifstream in( IN_HEAD + ".variant.vcf" );
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
          String chr, junk, ref, alt, Q, defs;
          int pos;
          iline >> chr >> pos >> junk >> ref >> alt >> Q >> junk >> junk >> defs;
          int N = ( chr == "Y" ? p_male : p_total );
          vec<String> p;
          if ( X == "F3" )
          {    if ( N == 5 ) p.push_back( "A1", "A2", "M", "D1", "D2" );
               else p.push_back( "D1", "D2" );    }
          else p.push_back( "M" );
          vec<vec<String>> rows;
          vec<String> row;
          row.push_back( "", "ref" );
          vec<String> alts;
          Tokenize( alt, {','}, alts );

          // Extract probabilities.  
          // Values < 0.005 are rounded down and > 0.995 up.

          vec<vec<double>> probs(N);
          
          vec<vec<String>> qs(N);
          
          Bool bMaleXHet = False;
          Bool bUnphysicalAlleleCount = False;
          vec<Bool> vbFlaky(N,false);
          for ( int i = 0; i < N; i++ )
          {    String x;
               iline >> x;
               vec<String> fields, q;
               Tokenize( x, {':'}, fields );
               Tokenize( fields[3], {','}, q );
               q.push_front( fields[2] );
               int allele_count = 0;
               int between_count = 0;
               for ( int j = 0; j < q.isize( ); j++ )
               {    int m = Min( (int) q[j].Int( ), iQCeiling );
                    double p = 1.0 - pow( 10, -m/10.0 );
                    if ( p > LO_PROB && p <= HI_PROB){
                        vbFlaky[i]=True; 
                        ++between_count;
                    }
                    else if ( p <= LO_PROB ){
                        p = 0;
                    }
                    else if ( p > HI_PROB ){
                        p = 1;
                        allele_count++;
                        ++count1;
                    }
                    probs[i].push_back(p);    }
               qs[i]=q;
               if ( allele_count > 2 || allele_count==0 ){
                   vbFlaky[i]=True;
               }
               if ( allele_count > 2 || between_count+allele_count==0 ){
                   bUnphysicalAlleleCount=True;
               }
               if ( i>2 && allele_count >1 ) bMaleXHet=True;
          }
          for( auto&entry : vbFlaky){ entry = entry && DISCOVAR_FILTER; }
          
          if(bUnphysicalAlleleCount){
              continue;
          }

          // Delete variants whose probabilities are all less than one.
/*
          vec<Bool> to_delete( alts.size( ) + 1, False );
          for ( int i = 1; i < alts.isize( ) + 1; i++ )
          {    double s = 0;
               Bool one = False;
               for ( int l = 0; l < N; l++ )
                    if ( probs[l][i] == 1 ) one = True;
               Bool can_be_one = False;
               for ( int l = 0; l < N; l++ )
                    if ( probs[l][i] >0 ) can_be_one = True;
               if ( ! (one || (N>1 &&can_be_one) ) ) to_delete[i] = True;    }
          vec<Bool> alts_del( alts.size( ) );
          for ( int j = 1; j < to_delete.isize( ); j++ )
               alts_del[j-1] = to_delete[j];
          EraseIf( alts, alts_del );
          for ( int i = 0; i < N; i++ ){
               EraseIf( probs[i], to_delete );
               EraseIf( qs[i], to_delete );
          }
          */

//the following line is buggy for the quintet
//          if ( probs[0].solo( ) ) continue;
          
          if ( probs[0].solo( ) ){
              bool bAllOne=true;
              for( size_t j=0;bAllOne&&j<probs.size();++j){
                  bAllOne=(probs[j][0]==1);
              }
              if(  bAllOne         //do not skip reference-only lines
                 ||probs.size()<2) //makes sure does not affect our single-sample analysis
                  continue;
          }

          // Check for X, Y and pseudoautosomal.

          Bool sex = ( chr == "X" || chr == "Y" );
          Bool pseudoautosomal = ( chr == "X" && ( pos < 2709520
               || ( pos >= 154584237 && pos < 154913754 ) ) );

          // Traverse all possibilities for probabilities that are strictly 
          // between 0 and 1.

          vec<pair<int,int>> between;
          for ( int l = 0; l < N; l++ )
          for ( int j = 0; j < probs[l].isize( ); j++ )
               if ( probs[l][j] > 0 && probs[l][j] < 1 ) between.push( l, j );
          Bool good = False;
          int64_t tries = IPow( 2, between.size( ) );
          for ( int64_t t = 0; t < tries; t++ )
          {    int64_t x = t;
               vec<vec<double>> probsx(probs);
               for ( int i = 0; i < between.isize( ); i++ , x/=2)
               {    if ( x % 2 == 0 ) 
                         probsx[ between[i].first ][ between[i].second ] = 0;
                    else probsx[ between[i].first ][ between[i].second ] = 1;    }

               // Test for failure of Mendelian inheritance etc.

               Bool funky = False;
               for ( int l = 0; l < N; l++ )
               {    if (vbFlaky[l]) continue;
                    for ( int j = 0; j < probsx[l].isize( ); j++ )
                         if ( probsx[l][j] > 0 && probsx[l][j] < 1 ) funky = True;
                    int n = Sum( probsx[l] );
                    if ( n == 0 || n > 2 ) funky = True;
                    if ( sex && !pseudoautosomal && n > 1 && ( X == "M" || l > 2 ) )                          funky = True;
                    if ( chr == "Y" && n != 1 ) funky = True;    }
               vec<vec<int>> alleles(N);
               for ( int l = 0; l < N; l++ )
               for ( int j = 0; j < probsx[l].isize( ); j++ )
                    if ( probsx[l][j] == 1 ) alleles[l].push_back(j);
               if ( chr != "Y" && X == "F3" )
               {    if ( !Meet( alleles[0], alleles[2] )  && !vbFlaky[0] && !vbFlaky[2] ) funky = True;
                    if ( !Meet( alleles[0], alleles[3] )  && !vbFlaky[0] && !vbFlaky[3] ) funky = True;
                    if ( !Meet( alleles[1], alleles[2] )  && !vbFlaky[1] && !vbFlaky[2] ) funky = True;
                    if ( !Meet( alleles[1], alleles[4] )  && !vbFlaky[1] && !vbFlaky[4] ) funky = True;
                    vec<int> P1 = alleles[2], P2 = alleles[2];
                    P1.append( alleles[3] ), P2.append( alleles[4] );
                    UniqueSort(P1), UniqueSort(P2);
                    if ( !BinSubset( alleles[0], P1 )&&!vbFlaky[0]&&!vbFlaky[2] &&!vbFlaky[3] ) funky = True;
                    if ( !BinSubset( alleles[1], P2 )&&!vbFlaky[1]&&!vbFlaky[2] &&!vbFlaky[4] ) funky = True;    }
               if ( !funky ) good = True;    }
          if (good) continue;

          if(INVOLVES_MALE_XHET && !bMaleXHet) continue;
          if(IGNORE_MALE_XHET && bMaleXHet) continue;

          // OK found a bad event, print it.
          
          

          if ( print_count == 0 ) cout << "\n";
          else
          {    cout << "----------------------------------------------------------"
                    << "-----------------------\n";    }
          cout << "[" << ++print_count << "] " << chr << ":" << pos << endl;
          
          int64_t rid; 
          if( chr=="X"){ rid = 22; }
          else if(chr=="Y"){ rid = 23; }
          else{ rid = chr.Int()-1; }
          const int64_t loc_id= pos-1;
          auto itr1 = genome[rid].begin() + loc_id-30;
          auto itr2 = genome[rid].begin() + loc_id;
          auto itr3 = genome[rid].begin() + loc_id+1;
          auto itr4 = genome[rid].begin() + loc_id+31;
          String sOut; sOut.reserve(30);
          std::transform(itr1,itr2,std::back_inserter(sOut),BaseToCharMapper());
          std::cout << sOut << " "; sOut.clear();
          std::transform(itr2,itr3,std::back_inserter(sOut),BaseToCharMapper());
          std::cout << sOut << " "; sOut.clear();
          std::transform(itr3,itr4,std::back_inserter(sOut),BaseToCharMapper());
          std::cout << sOut << std::endl;
          
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
          for ( int i = 0; i < N; i++ )
          {    vec<String> row;
               row.push_back( p[i] );
               for ( int j = 0; j < probs[i].isize( ); j++ )
               {    row.push_back(ToString(probs[i][j],5,true)+"("+qs[i][j]+")"); 
/*                    row.push_back(qs[i][j]);
                    if ( probs[i][j] == 1 ) row.push_back( "1" );
                    else if ( probs[i][j] == 0 ) row.push_back( "0" );
                    else row.push_back( ToString( probs[i][j], 2 ) ); */   }
               rows.push_back(row);    }
          PrintTabular( cout, rows, 2, "l" + String( N, 'l' ) );    }
     if(VERBOSITY) std::cout << "number of probability above HI_PROB " << count1 << std::endl;
}
