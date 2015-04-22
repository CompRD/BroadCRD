///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Compute transition/transversion ratios.  This creates a cached file dbsnp.cache.
// Currently hardwired for DISCOVAR and the filtering it uses.
//
// Only uses variants on chr1-22,X,Y.
// Assumes chromosome naming in vcf is 1,...,22,X,Y (without prefix chr).

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "STLExtensions.h"
#include "TokenizeString.h"
#include "util/TextTable.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(VCF,
          "/wga/scr4/human_assemblies/1/v9/v9_combined.filtered.vcf",
          "VCF file to use as input (currently must be produced by DISCOVAR)");
     CommandArgument_String_OrDefault_Doc(CHROM, "",
          "if specified, use only the given chromosome");
     CommandArgument_Int_OrDefault_Doc(MIN_DIST, 0, 
          "exclude variants having another variant closer than this");
     CommandArgument_Int_OrDefault_Doc(MAX_DIST, 1000000000, 
          "exclude variants that do not have another variant closer than this");
     CommandArgument_Int_OrDefault_Doc(MIN_DIST_INDEL, 0, 
          "exclude variants having an indel variant closer than this");
     CommandArgument_String_OrDefault_Doc(TYPE, "all", "hom or het or all");
     CommandArgument_Bool_OrDefault_Doc(DETAILS, False, 
          "to print the variants");
     CommandArgument_Bool_OrDefault_Doc(REMOVE_DBSNP, True, 
          "ignore SNPs in dbSNP");
     CommandArgument_IntSet_OrDefault_Doc(DIST_INDEL_BIN, "{}",
          "upper-bound of bins of ti/tv log");
     EndCommandArguments;

     // DISCOVAR constants.

//     const double cutoff1 = 0.995;
//     const double cutoff2 = 0.9;

     // Load genome.
     
     vecbasevector genome( "/wga/scr4/bigrefs/human19/genome.fastb" );

     // Load dbSNP.

     vec< pair< pair<String,int>, pair<String,String> > > dbsnp;
     String line, junk, chr, ref, alt, gen,filter;
     int pos;
     if ( IsRegularFile( "dbSNP.cache" ) )
          BinaryReader::readFile( "dbSNP.cache", &dbsnp );
     else
     {    fast_ifstream din( "/humgen/gsa-hpprojects/GATK/bundle/current/b37/"
               "dbsnp_138.b37.excluding_sites_after_129.vcf" );
          vec<String> dlines;
          dlines.reserve(14000000);
          while(1)
          {    getline( din, line );
               if ( din.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               if ( line[0] != 'X' && line[0] != 'Y' 
                    && ( line[0] < '0' || line[0] > '9' ) )
               {     continue;     }
               dlines.push_back(line);    }
          dbsnp.reserve( dlines.size( ) );
          for ( int z = 0; z < dlines.isize( ); z++ )
          {    const String& line = dlines[z];
               istrstream iline( line.c_str( ) );
               iline >> chr >> pos >> junk >> ref >> alt;
               if ( ref.size( ) > 1 || alt.size( ) > 1 ) continue;
               dbsnp.push( make_pair(chr,pos), make_pair(ref,alt) );    }
          ParallelSort(dbsnp);
          BinaryWriter::writeFile( "dbSNP.cache", dbsnp );    }

     //binning
     bool bBinning=DIST_INDEL_BIN.size()>0;
     if(bBinning && MIN_DIST_INDEL!=0){
         std::cout << "WARNING: reseting MIN_DIST_INDEL to 0, since we are binning." << std::endl;
         MIN_DIST_INDEL=0;
     }
     UniqueSort(DIST_INDEL_BIN);
     vec<int> ti_dist_bin(DIST_INDEL_BIN.size()+1,0);
     vec<int> tv_dist_bin(DIST_INDEL_BIN.size()+1,0);

     // Load DISCOVAR variants.

     fast_ifstream in(VCF);
     vec<String> lines;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( "#", 0 ) ) continue;
          if ( line[0] != 'X' && line[0] != 'Y' 
               && ( line[0] < '0' || line[0] > '9' ) )
          {     continue;     }
//          istrstream iline( line.c_str( ) );
//          iline >> chr >> pos >> junk >> ref >> alt >> junk >> filter 
//               >> junk >> junk >> gen;
//          if(!filter.Contains("PASS")) continue;
          lines.push_back(line);    }
     vec<String> CHR, REF, ALT;
     vec<int> POS;
     for ( int z = 0; z < lines.isize( ); z++ )
     {    const String& line = lines[z];
          istrstream iline( line.c_str( ) );
          iline >> chr >> pos >> junk >> ref >> alt >> junk >> junk 
               >> junk >> junk >> gen;
          CHR.push_back(chr), POS.push_back(pos);
          REF.push_back(ref), ALT.push_back(alt);    }

     // Find transitions and transversions.

     int transitions = 0, transversions = 0;
     for ( int z = 0; z < lines.isize( ); z++ )
     {    const String& line = lines[z];
          istrstream iline( line.c_str( ) );
          iline >> chr >> pos >> junk >> ref >> alt >> junk >> filter 
               >> junk >> junk >> gen;
          if(!filter.Contains("PASS")) continue;
          if ( CHROM != "" && chr != CHROM ) continue;

          int dist = 1000000000;
          if ( z > 0 && chr == CHR[z-1] ) dist = Min( dist, pos - POS[z-1] );
          if ( z < lines.isize( ) - 1 && chr == CHR[z+1] )
               dist = Min( dist, POS[z+1] - pos );
          if ( dist < MIN_DIST || dist > MAX_DIST ) continue;

          int dist_to_indel=std::numeric_limits<int>::max();

          Bool too_close = False;
          for ( int w = z - 1; w >= 0; w-- )
          {    if ( CHR[w] != CHR[z] ) break;
               if(REF[w].size()!=ALT[w].size()) dist_to_indel = std::min(dist_to_indel,pos-POS[w]);
               if ( pos - POS[w] >= MIN_DIST_INDEL ) break;
               if ( REF[w].size( ) == 1 && ALT[w].size( ) == 1 ) continue;
               too_close = True;    }
          for ( int w = z + 1; w < lines.isize( ); w++ )
          {    if ( CHR[w] != CHR[z] ) break;
               if(REF[w].size()!=ALT[w].size()) dist_to_indel = std::min(dist_to_indel,POS[w]-pos);
               if ( POS[w] - pos >= MIN_DIST_INDEL ) break;
               if ( REF[w].size( ) == 1 && ALT[w].size( ) == 1 ) continue;
               too_close = True;    }
          if (too_close) continue;

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

          String R = ref;
          String A;
          if ( alleles.solo( ) ) A = alts[ alleles[0] - 1 ];
          else A = alts[ alleles[1] - 1 ];

          if ( R.size( ) > 1 || A.size( ) > 1 ) continue;

          Bool hom = alleles.solo( );

          if ( REMOVE_DBSNP && 
               BinMember( dbsnp, make_pair( make_pair(chr,pos), make_pair(R,A) ) ) ) 
          {    continue;    }

          if ( TYPE == "hom" && !hom ) continue;
          if ( TYPE == "het" && hom ) continue;

          Bool transition = False;
          if ( R == "A" && A == "G" ) transition = True;
          if ( R == "G" && A == "A" ) transition = True;
          if ( R == "C" && A == "T" ) transition = True;
          if ( R == "T" && A == "C" ) transition = True;

          size_t bin_idx=0;
          for(;   bin_idx<DIST_INDEL_BIN.size()
               && dist_to_indel >= DIST_INDEL_BIN[bin_idx]
              ;++bin_idx)
          {}

          if (transition){
              transitions++;
              ++ti_dist_bin[bin_idx];
          }
          else{
              transversions++;
              ++tv_dist_bin[bin_idx];
          }

          const int flank = 10;
          String context( 2*flank + 1 );
          int ichr;
          if ( chr == "X" ) ichr = 22;
          else if ( chr == "Y" ) ichr = 23;
          else ichr = chr.Int( ) - 1;
          for ( int j = 0; j < context.isize( ); j++ )
               context[j] = as_base( genome[ichr][ pos - flank + j ] );
          // context[flank] = A[0];

          if (DETAILS)
          {    if (hom)
               {    cout << "hom: " << context << " " << chr << ":" << pos 
                         << " | " << ref << " --> " << A;    }
               else 
               {    cout << "het: " << context << " " << chr << ":" << pos 
                         << " | " << ref << " --> " << A;    }
               cout << " [" << ( transition ? "transition" : "transversion" ) 
                    << "]\n";    }    }

     if (DETAILS) cout << "\n";
     cout << "transitions = " << transitions << ", transversions = "
          << transversions << ", ratio = " 
          << double(transitions)/double(transversions) << endl;
     {
         TextTable table;
         table <<DoubleLine;
         table <<"begin" << Tab << "end" << Tab << "Ti" << Tab << "Tv" << Tab << "Ti/Tv" << EndRow;
         table <<DoubleLine;
         for(size_t bb=0;bb<DIST_INDEL_BIN.size();++bb){
             table << (bb==0?0:DIST_INDEL_BIN[bb-1]) << Tab << DIST_INDEL_BIN[bb] << Tab << ti_dist_bin[bb] << Tab << tv_dist_bin[bb] << Tab << double(ti_dist_bin[bb])/double(tv_dist_bin[bb]) << EndRow;
         }
         table << (DIST_INDEL_BIN.size()==0?0:DIST_INDEL_BIN.back()) << Tab << "inf" << Tab
               << ti_dist_bin.back() << Tab << tv_dist_bin.back() << Tab << double(ti_dist_bin.back())/double(tv_dist_bin.back()) << EndRow;
         table.Print(std::cout);

     }
}
