///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SearchHuman.  Search a bunch of NA12878 read sets for a given sequence S.
// Search only in a given region R, if provided, or else search the entire
// genome (which is very slow).  If S is not provided, return all reads.


// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "paths/long/MakeKmerStuff.h"
#include "system/SortInPlace.h"

void Search( const String& S, const vec<String>& perflets, 
     const vecbasevector& match, const String& R, const String& flowcell, 
     const String& picard_run, const vec<int>& lanes, const vec<String>& libs, 
     const Bool get_partners, const Bool print_headers, const Bool headers_only )
{
     // cout << "\nsearching " << flowcell << endl << endl;
     basevector bfw(S);
     basevector brc(S);
     brc.ReverseComplement( );

     String shead = "samtools view -X /seq/picard/";
     String dbam = ".aligned.duplicates_marked.bam";
     String look;
     if ( S != "" )
          look = " " + R + " | egrep \"" + S + "|" + brc.ToString( ) + "\"";
     else look = " " + R;

     for ( int i = 0; i < lanes.isize( ); i++ )
     for ( int j = 0; j < libs.isize( ); j++ )
     {    if ( print_headers || headers_only )
          {    System( "samtools view -H /seq/picard/" + flowcell + "/" + picard_run 
                    + "/" + ToString(lanes[i]) + "/" + libs[j] + "/" + flowcell
                         + "." + ToString(lanes[i]) + dbam );    }
          if (headers_only) continue;
          fast_pipe_ifstream in1( shead + flowcell + "/" + picard_run + "/"
               + ToString(lanes[i]) + "/" + libs[j] + "/" + flowcell
               + "." + ToString(lanes[i]) + dbam + look );
          String line, x;
          vec<String> ids;
          if ( match.size( ) > 0 )
          {    const int max_pile = 5000000;
               vecbasevector pile;
               vec<String> pile_line;
               pile.reserve(max_pile), pile_line.reserve(max_pile);
               while(1)
               {    while(1)
                    {    getline( in1, line );
                         if ( in1.fail( ) ) break;
                         int tabs = 0, start = 0, stop = 0;
                         for ( int t = 0; t < line.isize( ); t++ )
                         {    if ( line[t] == '\t' )
                              {    tabs++;
                                   if ( tabs == 9 ) start = t + 1;
                                   else if ( tabs == 10 )
                                   {    stop = t;
                                        break;    }    }    }
                         x = line.substr( start, stop - start );
                         if ( !x.Contains( "N" ) )
                         {    pile.push_back( basevector(x) );
                              pile_line.push_back(line);    }
                         if ( (int) pile.size( ) == max_pile ) break;    }
                    const int K = 100;
                    vecbasevector all(match);
                    vec<int> marked( pile.size( ), False );
                    all.Append(pile);
                    vec< triple<kmer<K>,int,int> > kmers_plus;
                    MakeKmerLookup2( all, kmers_plus );
                    for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
                    {    int64_t j;
                         for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                         {    if ( kmers_plus[j].first != kmers_plus[i].first ) 
                                   break;    }
                         Bool valid = False;
                         for ( int64_t k = i; k < j; k++ )
                         {    if ( kmers_plus[k].second < (int) match.size( ) ) 
                                   valid = True;    }
                         if (valid)
                         {    for ( int64_t k = i; k < j; k++ )
                              {    if ( kmers_plus[k].second >= (int) match.size( ) )
                                   {    int id = kmers_plus[k].second 
                                             - (int) match.size( );
                                        marked[id] = True;    }    }    }
                         i = j - 1;    }
                    for ( int i = 0; i < (int) pile.size( ); i++ )
                    {    if ( marked[i] )
                         {    if ( !get_partners ) cout << pile_line[i] << "\n";
                              else ids.push_back( pile_line[i].Before( "\t" ) );
                              }    }
                    pile.clear( );
                    pile_line.clear( );
                    if ( in1.fail( ) ) break;    }    }

          else
          {    while(1)
               {    getline( in1, line );
                    if ( in1.fail( ) ) break;
                    if ( perflets.nonempty( ) )
                    {    int tabs = 0, start = 0, stop = 0;
                         for ( int t = 0; t < line.isize( ); t++ )
                         {    if ( line[t] == '\t' )
                              {    tabs++;
                                   if ( tabs == 9 ) start = t + 1;
                                   else if ( tabs == 10 )
                                   {    stop = t;
                                        break;    }    }    }
                         x = line.substr( start, stop - start );
                         if ( !BinMember( perflets, x ) ) continue;    }
                    if ( !get_partners) cout << line << "\n";
                    else ids.push_back( line.Before( "\t" ) );    }    }

          if (get_partners)
          {    ParallelUniqueSort(ids);
               fast_pipe_ifstream in2( shead + flowcell + "/" + picard_run + "/"
                    + ToString(lanes[i]) + "/" + libs[j] + "/" + flowcell
                    + "." + ToString(lanes[i]) + dbam + " " + R );
               while(1)
               {    getline( in2, line );
                    if ( in2.fail( ) ) break;
                    if ( line.Contains( "#", 0 ) ) continue;
                    if ( BinMember( ids, line.Before( "\t" ) ) )
                         cout << line << "\n";    }    }    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;

     // Screening criteria:

     CommandArgument_String_OrDefault_Doc( S, "", "sequence to search for" );
     CommandArgument_String_OrDefault_Doc( PERF, "", 
          "if provided, require perfect match to this fasta file; only implemented "
          "for 250 base reads" );
     CommandArgument_String_OrDefault_Doc( MATCH, "", 
          "if provided, require that read have 100 base perfect match to something "
          "in this fastb file; optimized for fastb = all reads in Fosmid pools and "
          "running on 512 GB machine" );
     CommandArgument_String_OrDefault_Doc( R, "", "region" );
     CommandArgument_String_OrDefault_Doc( BAMS, "all", "all or free" );
     CommandArgument_String_OrDefault_Doc( PERSON, "child", "child or mom or pop" );
     CommandArgument_Bool_OrDefault_Doc( GET_PARTNERS, False, 
          "return also partners of reads that are found" );
     CommandArgument_Bool_OrDefault_Doc( PRINT_HEADERS, False, 
          "print header lines of bam files" );
     CommandArgument_Bool_OrDefault_Doc( HEADERS_ONLY, False, 
          "print only header lines of bam files" );
     EndCommandArguments;

     for ( int j = 0; j < S.isize( ); j++ )
          S[j] = toupper( S[j] );

     basevector bfw(S);
     basevector brc(S);
     brc.ReverseComplement( );

     vec<String> perflets;
     if ( PERF != "" ) 
     {    vecbasevector G;
          FetchReads( G, 0, PERF );
          const int K = 250;
          int N = 0;
          for ( int g = 0; g < (int) G.size( ); g++ )
               N += 2 * ( G[g].isize( ) - K + 1 );
          perflets.resize( N, String(K) );
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 )
               {    for ( int g = 0; g < (int) G.size( ); g++ )
                         G[g].ReverseComplement( );    }
               vec<pair<int,int>> locs;
               for ( int g = 0; g < (int) G.size( ); g++ )
               for ( int j = 0; j <= G[g].isize( ) - K; j++ )
                    locs.push( g, j );
               int pos = ( pass == 1 ? 0 : N/2 );
               #pragma omp parallel for
               for ( int l = 0; l < locs.isize( ); l++ )
               {    int g = locs[l].first, j = locs[l].second;
                    perflets[pos+l] = basevector( G[g], j, K ).ToString( );    }    }
          sortInPlaceParallel( perflets.begin( ), perflets.end( ) );    }

     vecbasevector match;
     if ( MATCH.size( ) > 0 ) match.ReadAll(MATCH);

     String shead = "samtools view /seq/picard/";
     String dbam = ".aligned.duplicates_marked.bam";
     String look = " " + R + " | egrep \"" + S + "|" + brc.ToString( ) + "\"";
     String flowcell, picard_run;

     if ( PERSON == "mom" )
     {    flowcell = "H06JHADXX";
          picard_run = "C1-508_2013-01-10_2013-01-13";
          {    vec<int> lanes;
               vec<String> libs;
               lanes.push_back( 1, 2 );
               libs.push_back( "Solexa-135853" );
               Search( S, perflets, match, R, flowcell, picard_run, lanes, libs, 
                    GET_PARTNERS, PRINT_HEADERS, HEADERS_ONLY );    }
          flowcell = "H06JUADXX";
          picard_run = "C1-508_2013-01-10_2013-01-13";
          {    vec<int> lanes;
               vec<String> libs;
               lanes.push_back( 2 );
               libs.push_back( "Solexa-135853" );
               Search( S, perflets, match, R, flowcell, picard_run, lanes, libs, 
                    GET_PARTNERS, PRINT_HEADERS, HEADERS_ONLY );    }
          return 0;    }

     if ( PERSON == "pop" )
     {    flowcell = "H03N7ADXX";
          picard_run = "C1-508_2013-01-07_2013-01-10";
          {    vec<int> lanes;
               vec<String> libs;
               lanes.push_back( 1, 2 );
               libs.push_back( "Solexa-135851" );
               Search( S, perflets, match, R, flowcell, picard_run, lanes, libs, 
                    GET_PARTNERS, PRINT_HEADERS, HEADERS_ONLY );    }
          flowcell = "H05F1ADXX";
          picard_run = "C1-508_2013-01-15_2013-01-18";
          {    vec<int> lanes;
               vec<String> libs;
               lanes.push_back( 2 );
               libs.push_back( "Solexa-135851" );
               Search( S, perflets, match, R, flowcell, picard_run, lanes, libs, 
                    GET_PARTNERS, PRINT_HEADERS, HEADERS_ONLY );    }
          return 0;    }

     // PCR-free, HiSeq2500, 250
     flowcell = "H01UJADXX";
     picard_run = "C1-508_2012-11-01_2012-11-04";
     {    vec<int> lanes;
          vec<String> libs;
          lanes.push_back( 1, 2 );
          libs.push_back( "Solexa-125532" );
          Search( S, perflets, match, R, flowcell, picard_run, lanes, libs, 
                    GET_PARTNERS, PRINT_HEADERS, HEADERS_ONLY );    }
     flowcell = "H06HDADXX";
     picard_run = "C1-508_2013-01-10_2013-01-13";
     {    vec<int> lanes;
          vec<String> libs;
          lanes.push_back( 1, 2 );
          libs.push_back( "Solexa-135852" );
          Search( S, perflets, match, R, flowcell, picard_run, lanes, libs, 
               GET_PARTNERS, PRINT_HEADERS, HEADERS_ONLY );    }
     flowcell = "H06JUADXX";
     picard_run = "C1-508_2013-01-10_2013-01-13";
     {    vec<int> lanes;
          vec<String> libs;
          lanes.push_back( 1 );
          libs.push_back( "Solexa-135852" );
          Search( S, perflets, match, R, flowcell, picard_run, lanes, libs, 
               GET_PARTNERS, PRINT_HEADERS, HEADERS_ONLY );    }

     if ( BAMS != "free" )
     {
         // PCR, HiSeq2500, 150
         flowcell = "H01V9ADXX";
         picard_run = "C1-308_2012-10-20_2012-10-22";
         {    vec<int> lanes;
              vec<String> libs;
              lanes.push_back( 1, 2 );
              libs.push_back( "Solexa-123653" );
              Search( S, perflets, match, R, flowcell, picard_run, lanes, libs, 
                    GET_PARTNERS, PRINT_HEADERS, HEADERS_ONLY );    }

         // PCR, HiSeq2500, 250
         flowcell = "H01RHADXX";
         picard_run = "C1-508_2012-10-24_2012-10-31";
         {    vec<int> lanes;
              vec<String> libs;
              lanes.push_back( 1, 2 );
              libs.push_back( "Solexa-126138" );
              Search( S, perflets, match, R, flowcell, picard_run, lanes, libs, 
                    GET_PARTNERS, PRINT_HEADERS, HEADERS_ONLY );    }

         // PCR-free, HiSeq2000, 150
         flowcell = "C0NJNACXX";
         picard_run = "C1-210_2012-06-30_2012-10-05";
         {    vec<int> lanes;
              vec<String> libs;
              lanes.push_back( 1, 2, 3 );
              libs.push_back( "Solexa-97534" );
              Search( S, perflets, match, R, flowcell, picard_run, lanes, libs, 
                    GET_PARTNERS, PRINT_HEADERS, HEADERS_ONLY );    }

         // PCR, HiSeq2500, 150
         flowcell = "H01UTADXX";
         picard_run = "C1-308_2012-10-22_2012-10-24";
         {    vec<int> lanes;
              vec<String> libs;
              lanes.push_back( 1, 2 );
              libs.push_back( "Solexa-123653" );
              Search( S, perflets, match, R, flowcell, picard_run, lanes, libs, 
                    GET_PARTNERS, PRINT_HEADERS, HEADERS_ONLY );    }    }     }
