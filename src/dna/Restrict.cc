///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Restrict.  Given a restriction enzyme specified e.g. by
// SITE=GG^CC
// find restriction fragments that would arise from a complete digest.
//
// The fragment sequence includes single-stranded tails, but the answer is
// presented as the 5'--3' strand of a double-stranded sequence, as if the
// bases complementary to the tail bases had been filled in.
//
// Site may have ambiguous bases in it.  Site must be palindromic.
//
// The substrate is specified by SEQ and CIRCULAR.  Ambiguous bases are allowed
// in the substrate but restriction sites will not match them.
//
// Set SIZES=True to see the fragment sizes and GC contents instead of the
// fragment bases.
//
// Set SHORT=True to get one-line list of sizes.
//
// Set HISTOGRAM={x1,...,xk} to only report a histogram of fragment sizes, using
// x1 < ... < xk as dividers.
//
// Instead of specifying SITE, one may specify TRYALL=1, which causes a builtin
// list of enzymes to be tried in succession.  If instead TRYALL=2,
// then all pairs of digests are done in parallel and the results of each pair
// merged.  Similarly, TRYALL=3 is allowed.  Also a list can be specified in
// ParseIntSet format.
//
// ALL=True: try all restriction enzymes.
// ALL_NEB=True: try all restriction enzymes from NEB.
//
// MIN_RATIO: if set, do not report results for restrictions yielding fragments of
// successive sizes a,b such that b/a < MIN_RATIO.
//
// MAX_GAP: if set, do not report results for restrictions yielding fragments of
// successive sizes a,b such that b-a > MAX_GAP.
//
// MIN_MAX: if set, do not report results for restrictions yielding fragments whose
// maximum size is < MIN_MAX.
//
// MIN_FRAGS: if set, do not report results for restrictions yielding less than
// that number of fragments.
//
// MIN_FRAGS_FLOOR: if set, do not report results for restrictions yielding less
// MIN_FRAGS fragments of length >= MIN_FRAGS_FLOOR.
//
// MIN_FRAG: if set, do not report results for restrictions yielding fragments
// smaller than MIN_FRAG.
//
// ENZYME: if specified and known to this program, use the given enzyme to define
// the site.
//
// RANGE="[a,b)": only return fragments in this size range.

#include <map>

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "math/Functions.h"

vec< pair<String,String> > cutters;

#define CUTTER(N,ENZYME,SITE) cutters.push( #ENZYME, #SITE )

void DefineCutters( )
{
     // Blunt cutters available from NEB and very well-behaved according
     // to their website:

     CUTTER(4,HaeIII,GG^CC);
     CUTTER(4,HpyCH4V,TG^CA);
     CUTTER(4,RsaI,GT^AC);
     CUTTER(4,Cac8I,GCN^NGC);
     CUTTER(4,NlaIV,GGN^NCC);
     CUTTER(5,BsaAI,YAC^GTR);
     CUTTER(5,MspA1I,CMG^CKG);
     CUTTER(5,MslI,CAYNN^NNRTG);
     CUTTER(6,EcoRV,GAT^ATC);
     CUTTER(6,FspI,TGC^GCA);
     CUTTER(6,MscI,TGG^CCA);
     CUTTER(6,PmlI,CAC^GTG);
     CUTTER(6,SmaI,CCC^GGG);
     CUTTER(6,StuI,AGG^CCT);
     CUTTER(6,BsaBI,GATNN^NNATC);
     CUTTER(6,PshAI,GACNN^NNGTC);

     // Other enzymes:

     CUTTER(4,MspI,C^CGG);
     CUTTER(4,Msp1,C^CGG); // alias for MspI
     CUTTER(6,BamHI,G^GATCC);
     CUTTER(6,HindIII,A^AGCTT);
     CUTTER(5,HinFI,G^ANTC);
     CUTTER(6,PstI,CTGCA^G);
     CUTTER(6,NdeI,CA^TATG);
     CUTTER(6,KpnI,GGTAC^C);
     CUTTER(6,XbaI,T^CTAGA);
     CUTTER(6,BglII,A^GATCT);
     CUTTER(6,NcoI,C^CATGG);
     CUTTER(6,PvuII,CAG^CTG);
     CUTTER(4,TaqI,T^CGA);

}

// BaseContains: determine if DNA sequence "target" contains generalized DNA
// sequence "query" at position target_pos.

Bool BaseContains( const String& target, const String& query, int target_pos )
{
    bool result = false;
    if ( target_pos + query.size() <= target.size() )
    {
        result = true;
        String::const_iterator end = query.end();
        String::const_iterator targetItr = target.begin(target_pos);
        for ( String::const_iterator itr = query.begin(); result && itr != end; ++itr, ++targetItr )
            if ( !GeneralizedBase::isGeneralizedBase(*itr) ||
                 !GeneralizedBase::isGeneralizedBase(*targetItr) ||
                 !GeneralizedBase::fromChar(*itr).matches(GeneralizedBase::fromChar(*targetItr)) )
                result = false;
    }
    return result;
}

// Restrict: restrict at a given restriction site.

void Restrict(

     /* input: */ const String& seq, ///< dna sequence to cut with enzyme
                  const Bool circular, ///< is the dna we are cutting circular one?
                  const String& Site, ///< Restriction site (e.g. "C^CGG")
                  const int low, ///< minimum fragment size to accept, [low,high)
                  const int high, ///< maximum fragment size to accept, [low,high)

     /* output (appended): */ vec<String>& fragments,
                              vec<int>& fragment_starts,
                              vec<int>& frag_sizes,
                              vec<double>& frag_gc,
                              int& unbroken_circles )

{
    if ( !Site.Contains( "^" ) ) {
        cout << "Site = " << Site << " does not contain ^.\n";
	exit(1);
    }
    String start = Site.Before( "^" ), stop = Site.After( "^" ), stop_rc;
    StringReverseComplement( stop, stop_rc );
    String site = start + stop, site_rc;
    StringReverseComplement( site, site_rc );

    // sanity checks: restriction site must be an rc-palindrome and
    // site string should contain valid bases
    ForceAssertEq( site, site_rc );
    ForceAssert( GeneralizedBase::areGeneralizedBases(start.begin(),start.end()) );
    ForceAssert( GeneralizedBase::areGeneralizedBases(stop.begin(),stop.end()) );

    String seq2; // sequence to cut

     if (circular) seq2 = seq + seq;
     else seq2 = seq;

     // n=size of the restriction site palindrome; N - size of the sequence to be cut
     int n = site.size( ), N = seq.size( );
     int last = -1;
     int top = ( circular ? 2*N - n : N - n + 1 );

     for ( int i = 0; i <= top; i++ ) {  // loop over dna sequence to be digested and look for restriction sites

         Bool linear_end = ( !circular && i == N - n + 1 ); // true if we reached the end of the non-circ sequence;
	                                                    // clearly, this should be the end of the last fragment

         if ( BaseContains( seq2, site, i ) || linear_end ) {
	     // found restriction site or reached the end of non-circular sequence, need to generate next fragment:

	     if ( last == -1 && circular ) {
	         if ( i >= N ) {
		     ++unbroken_circles;
		     return;
		 }
		 last = i;
	     } else {
	         // we are cutting non-circular dna, or it is circular but it is not the first site found
	         static String frag;
		 frag.resize(0);
		 // if we have not seen restriction site before, then start from the beginning of the dna seq;
		 // otherwise start from the last position where restriction site was found adjusted by the site size:
		 int Start = ( last == -1 ? 0
                         : last + Min( start.isize( ), stop.isize( ) ) );

		 // if we are at the end of non-circular (linear) seq - that's where the fragment ends!
		 // otherwise, the fragment ends at the location of the restriction site we just found -
		 // adjusted by the site prefix/suffix size
		 int Stop = ( linear_end ? N
                         : i + Max( start.isize( ), stop.isize( ) ) );

		 // if the size of the current fragment passes the size filter, record it:
		 if ( Stop - Start >= low && Stop - Start < high ) 
                 {    frag = seq2.substr( Start, Stop - Start );
                      Bool illegal = False;
                      for ( int j = 0; j < frag.isize( ); j++ )
                      {    char c = toupper( frag[j] );
                           if ( c != 'A' && c != 'C' && c != 'G' && c != 'T' )
                                illegal = True;    }
                      if ( !illegal )
		      {    frag_sizes.push_back( frag.size( ) );
		           frag_gc.push_back( GcPercent(frag) );
		           fragments.push_back(frag);
		           fragment_starts.push_back(Start);    }    }
		 last = i;
	     }
	     if ( i >= N ) break;
	 }
     }
     if ( circular && fragments.empty( ) ) ++unbroken_circles;
     SortSync( frag_sizes, frag_gc );
}

// ParallelRestrict: like Restrict, but takes a Site which is of the form s1|...|sn,
// where n >= 1, the si are restriction sites, and the union of the corresponding
// digests is returned.

void ParallelRestrict( const String& seq, const Bool circular, const String& Site,
     const int low, const int high, vec<String>& fragments,
     vec<int>& fragment_starts, vec<int>& frag_sizes, vec<double>& frag_gc,
     int& unbroken_circles )
{    vec<char> sep(1);
     sep[0] = '|';
     vec<String> sites;
     Tokenize( Site, sep, sites );
     for ( int i = 0; i < sites.isize( ); i++ )
     {    Restrict( seq, circular, sites[i], low, high, fragments, fragment_starts,
               frag_sizes, frag_gc, unbroken_circles );    }
     SortSync( frag_sizes, frag_gc );    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault(SITE, "");
     CommandArgument_String_Doc(SEQ, "fasta file to be operated on");
     CommandArgument_Bool(CIRCULAR);
     CommandArgument_Bool_OrDefault(SIZES, False);
     CommandArgument_Bool_OrDefault(MASS, False);
     CommandArgument_Bool_OrDefault(SHORT, False);
     CommandArgument_Bool_OrDefault(ALL, False);
     CommandArgument_Bool_OrDefault(ALL_NEB, False);
     CommandArgument_String_OrDefault(TRYALL, "");
     CommandArgument_Double_OrDefault(MIN_RATIO, 1.0);
     CommandArgument_Int_OrDefault(MAX_GAP, 0);
     CommandArgument_Int_OrDefault(MIN_MAX, 0);
     CommandArgument_Int_OrDefault(MIN_FRAGS, 0);
     CommandArgument_Int_OrDefault(MIN_FRAGS_FLOOR, 0);
     CommandArgument_Int_OrDefault(MIN_FRAG, 0);
     CommandArgument_String_OrDefault(ENZYME, "");
     CommandArgument_String_OrDefault(RANGE, "");
     CommandArgument_String_OrDefault_Doc(OUT_FASTA,"",
                  "if specified, digest fasta will be written into "
		  "this file or into stdout otherwise (for backward compatibility)");
     CommandArgument_String_OrDefault(HISTOGRAM, "");
     EndCommandArguments;

     // Check for RANGE.

     if ( ! IsRegularFile(SEQ) ) 
     {    cout << "ERROR: Sequence file " << SEQ << " does not exist" << endl;
          exit(1);    }
     int low = 0, high = 1000000000;
     if ( RANGE != "" )
     {    low = RANGE.Between( "[", "," ).Int( );
          high = RANGE.Between( ",", ")" ).Int( );    }

     // Set up for histogram.

     vec<int> hist_div, all_sizes;
     if ( HISTOGRAM != "" ) 
     {    ParseIntSet( HISTOGRAM, hist_div );
          hist_div.push_back( 0, 1000000000 );
          UniqueSort(hist_div);    }

     // Deal with OUT_FASTA.

     ofstream fasta_out;
     ostream * out_ptr;
     if ( OUT_FASTA != "" ) 
     {    if ( ! OUT_FASTA.Contains(".fasta",-1) ) OUT_FASTA += ".fasta";
          fasta_out.open(OUT_FASTA.c_str(), ios::out);
          out_ptr = & fasta_out;    } 
     else out_ptr = & cout;

     // Define restriction sites to use.

     DefineCutters( );
     if ( ENZYME != "" )
     {    ForceAssert( SITE == "" );
          for ( int i = 0; i < cutters.isize( ); i++ )
               if ( cutters[i].first == ENZYME ) SITE = cutters[i].second;
          if ( SITE == "" )
          {    cout << "ENZYME unknown.\n";
               exit(1);    }    }
     int opts = 0;
     if ( SITE != "" ) ++opts;
     if ( TRYALL != "" ) ++opts;
     if (ALL) ++opts;
     if (ALL_NEB) ++opts;
     ForceAssertEq( opts, 1 );
     vec<int> tryall;
     if ( TRYALL != "" ) ParseIntSet( TRYALL, tryall );
     for ( int i = 0; i < tryall.isize( ); i++ )
     {    ForceAssertGe( tryall[i], 1 );
          ForceAssertLe( tryall[i], 3 );    }
     vec<String> sites, names;
     if ( TRYALL == "" && !ALL &&!ALL_NEB ) sites.push_back(SITE);
     if ( Member( tryall, 1 ) )
     {    for ( int i = 0; i < cutters.isize( ); i++ )
          {    sites.push_back( cutters[i].second );
               names.push_back( cutters[i].first );    }    }
     if ( Member( tryall, 2 ) )
     {    for ( int i1 = 0; i1 < cutters.isize( ); i1++ )
          {    for ( int i2 = i1 + 1; i2 < cutters.isize( ); i2++ )
               {    sites.push_back( cutters[i1].second + "|" + cutters[i2].second );
                    names.push_back( cutters[i1].first 
                         + "|" + cutters[i2].first );    }    }    }
     if ( Member( tryall, 3 ) )
     {    for ( int i1 = 0; i1 < cutters.isize( ); i1++ )
          {    for ( int i2 = i1 + 1; i2 < cutters.isize( ); i2++ )
               {    for ( int i3 = i2 + 1; i3 < cutters.isize( ); i3++ )
                    {    sites.push_back( cutters[i1].second 
                              + "|" + cutters[i2].second
                              + "|" + cutters[i3].second );
                         names.push_back( cutters[i1].first 
                              + "|" + cutters[i2].first
                              + "|" + cutters[i3].first );    }    }    }    }
     if ( ALL || ALL_NEB )
     {    map<String,String> to_name;
          fast_ifstream in( "/wga/scr4/MolBio/Enzymes/link_withrefm" );
          String line, enzyme, site, suppliers;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( !line.Contains( "<1>", 0 ) ) continue;
               enzyme = line.After( "<1>" );
               getline( in, line );
               getline( in, line );
               ForceAssert( line.Contains( "<3>", 0 ) );
               site = line.After( "<3>" );
               if ( site.Contains( "(" ) || site.Contains( "?" )
                    || site.Contains( "," ) || !site.Contains( "^" ) )
               {    continue;    }
               String start = site.Before( "^" ), stop = site.After( "^" ), stop_rc;
               StringReverseComplement( stop, stop_rc );
               String Site = start + stop, Site_rc;
               StringReverseComplement( Site, Site_rc );
               if ( Site != Site_rc ) continue;
               getline( in, line );
               getline( in, line );
               getline( in, line );
               getline( in, line );
               ForceAssert( line.Contains( "<7>", 0 ) );
               suppliers = line.After( "<7>" );
               if ( ALL_NEB && !suppliers.Contains( "N" ) ) continue;
               if ( to_name[site] == "" ) to_name[site] = enzyme;
               else to_name[site] = to_name[site] + "|" + enzyme;    }
          for ( map<String,String>::iterator i = to_name.begin( );
               i != to_name.end( ); ++i )
          {    sites.push_back( i->first );
               names.push_back( i->second );    }    }

     // Go through the restriction sites.

     for ( int is = 0; is < sites.isize( ); is++ )
     {    int total = 0;
          const String& Site = sites[is];
          vec<String> seq, seqnames;
          fast_ifstream in(SEQ);
          String line;
          String x;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( ">", 0 ) )
               {    if ( seqnames.nonempty( ) ) seq.push_back(x);
                    seqnames.push_back( line.After( ">" ) );
                    x.clear( );    }
               else
               {    for ( int j = 0; j < line.isize( ); j++ )
                         x.push_back( line[j] );    }    }
          seq.push_back(x);
          ForceAssertEq( seq.size( ), seqnames.size( ) );
          ForceAssertGt( seq.size( ), 0u );

          // Go through the sequences.

          for ( int sx = 0; sx < seq.isize( ); sx++ ) 
          {    if ( SHORT || SIZES || MASS ) cout << "Sequence " << sx+1 << "\n";

               // Dice.

               vec<String> fragments;
               vec<int> fragment_starts, frag_sizes;
               vec<double> frag_gc;
               int unbroken_circles = 0;
               ParallelRestrict( seq[sx], CIRCULAR, Site, low, high, fragments,
                    fragment_starts, frag_sizes, frag_gc, unbroken_circles );

               // Check MIN_RATIO, etc. options.

               Bool reject = False;
               for ( int i = 1; i < frag_sizes.isize( ); i++ )
               {    int a = frag_sizes[i-1], b = frag_sizes[i];
                    if ( double(b)/double(a) < MIN_RATIO ) reject = True;
                    if ( MAX_GAP > 0 && b - a > MAX_GAP ) reject = True;    }
               if ( frag_sizes.nonempty( ) && frag_sizes.back( ) < MIN_MAX )
                    reject = True;
               int nfrags = 0;
               for ( int i = 0; i < frag_sizes.isize( ); i++ )
                    if ( frag_sizes[i] >= MIN_FRAGS_FLOOR ) ++nfrags;
               if ( MIN_FRAGS > 0 && nfrags < MIN_FRAGS ) reject = True;
               if ( frag_sizes.nonempty( ) && frag_sizes.front( ) < MIN_FRAG )
                    reject = True;
               if (reject) continue;

               // Print results.

               String SitePlus = Site;
               if ( TRYALL != "" || ALL || ALL_NEB )
                    SitePlus += " (" + names[is] + ")";
               if ( SHORT || SIZES || MASS )
               {    cout << "Cutting at " << SitePlus << ":";
                    cout << "\n";    }
               if ( HISTOGRAM != "" ) all_sizes.append(frag_sizes);
               else if ( unbroken_circles == 1 ) 
               {    cout << "includes unbroken circle\n";    }
               else if ( unbroken_circles > 1 )
               {    cout << "includes " << unbroken_circles
                         << " unbroken circles\n";    }
               else if (SIZES)
               {    for ( int i = 0; i < frag_sizes.isize( ); i++ )
                    {    cout << frag_sizes[i] << " " << setiosflags(ios::fixed)
                              << setprecision(1) << setw(4) << frag_gc[i] << "%"
                              << resetiosflags(ios::fixed) << "\n";    }    }
               else if (MASS)
               {    for ( int i = 0; i < frag_sizes.isize( ); i++ )
                         total += frag_sizes[i];    }
               else if (SHORT)
               {    for ( int i = 0; i < frag_sizes.isize( ); i++ )
                    {    if ( i > 0 ) cout << ",";
                         cout << frag_sizes[i];    }
                    cout << "\n";    }
               else 
               {    for ( int i = 0; i < fragments.isize( ); i++ ) 
                    {    int length = fragments[i].size( );
		         int start = fragment_starts[i];
		         int stop = ( start + length ) % seq[sx].isize( );
		         static basevector b;
		         b.SetFromString( fragments[i] );
		         b.Print( *out_ptr, seqnames[sx] + ", bases " 
                              + ToString(start)
                              + "-" + ToString(stop) + ", " + SitePlus
                              + " res. frag. " + ToString(i+1) );    }    }
          if (MASS) cout << "total bases in fragments = " << total << endl;    }
     if ( OUT_FASTA != "" ) fasta_out.close();    }

     // Print histogram.

     if ( HISTOGRAM != "" )
     {    vec<int> piles( hist_div.isize( ) - 1, 0 );
          for ( int j = 0; j < all_sizes.isize( ); j++ )
          {    for ( int i = 0; i < hist_div.isize( ) - 1; i++ )
               {    if ( all_sizes[j] >= hist_div[i]
                         && all_sizes[j] < hist_div[i+1] )
                    {    piles[i]++;    }    }    }
          cout << "histogram of sizes:\n";
          for ( int j = 0; j < piles.isize( ); j++ )
          {    cout << hist_div[j] << " <= s <= " << hist_div[j+1]
                    << ": " << piles[j] << "\n";    }    }    }
