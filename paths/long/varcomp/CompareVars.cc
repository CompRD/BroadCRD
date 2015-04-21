///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// CompareVars.  Preliminary program to compare variant calls from DISCOVAR, GATK,
// Cortex, and the Fosmid reference sequences.  Output to 
// /wga/scr4/$user/CompareVars/n, where n is the id of a Fosmid.

#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "lookup/PerfectLookup.h"
#include "math/Functions.h"
#include "paths/long/varcomp/CompareVarsTools.h"
#include "paths/long/varcomp/FosmidifyEvents.h"
#include "random/Random.h"

namespace { // open anonymous namespace

const double cutoff1 = 0.995;
const double cutoff2 = 0.9;

// Class varcall.  The alt_nhood and ref_nhood are sequences that are symmetric
// about the locus, on the given allele.

class varcall {

     public:

     varcall( ) { state = "uncalled"; }
     varcall( const String& state ) : state(state) { }

     void clear(){
         state="uncalled";
         alt_nhood.clear();
         ref_nhood.clear();
     };

     String state;
     String alt_nhood;
     String ref_nhood;

     String GetAltNhood( const int left, const int right )
     {    int p = alt_nhood.Position( "|" );
          return alt_nhood.substr( p - left, left ) 
               + alt_nhood.substr( p + 1, right );    }

     String GetRefNhood( const int left, const int right )
     {    int p = ref_nhood.Position( "|" );
          return ref_nhood.substr( p - left, left ) 
               + ref_nhood.substr( p + 1, right );    }

     // Build tiny neighborhoods, at least the given length on each side 
     // (if possible) and long enough to distinguish alt from ref, if possible.

     void BuildMini( const int flank, String& altn, String& refn )
     {    int alt_left_dist = alt_nhood.Position( "|" );
          int alt_right_dist = alt_nhood.isize( ) - alt_left_dist - 1;
          int ref_left_dist = ref_nhood.Position( "|" );
          int ref_right_dist = ref_nhood.isize( ) - ref_left_dist - 1;
          int alt_left = Min( flank, alt_left_dist );
          int alt_right = Min( flank, alt_right_dist );
          int ref_left = Min( flank, ref_left_dist );
          int ref_right = Min( flank, ref_right_dist );
          String alt(alt_nhood), ref(ref_nhood);
          alt.ReplaceBy( "|", "" ), ref.ReplaceBy( "|", "" );
          if ( ref.Contains( GetAltNhood( alt_left, alt_right ) ) )
          {    int alt_left_plus = alt_left;
               Bool fixed = False;
               while( alt_left_plus < alt_left_dist )
               {    alt_left_plus++;
                    if ( !ref.Contains( GetAltNhood( alt_left_plus, alt_right ) ) )
                    {    fixed = True;
                         break;    }    }
               if (fixed) alt_left = alt_left_plus;
               else
               {    int alt_right_plus = alt_right;
                    while( alt_right_plus < alt_right_dist )
                    {    alt_right_plus++;
                         if ( !ref.Contains( 
                              GetAltNhood( alt_left, alt_right_plus ) ) )
                         {    fixed = True;
                              break;    }    }
                    if (fixed) alt_right = alt_right_plus;    }    }
          if ( alt.Contains( GetRefNhood( ref_left, ref_right ) ) )
          {    int ref_left_plus = ref_left;
               Bool fixed = False;
               while( ref_left_plus < ref_left_dist )
               {    ref_left_plus++;
                    if ( !alt.Contains( GetRefNhood( ref_left_plus, ref_right ) ) )
                    {    fixed = True;
                         break;    }    }
               if (fixed) ref_left = ref_left_plus;
               else
               {    int ref_right_plus = ref_right;
                    while( ref_right_plus < ref_right_dist )
                    {    ref_right_plus++;
                         if ( !alt.Contains( 
                              GetRefNhood( ref_left, ref_right_plus ) ) )
                         {    fixed = True;
                              break;    }    }
                    if (fixed) ref_right = ref_right_plus;    }    }
          altn = GetAltNhood( alt_left, alt_right );
          refn = GetRefNhood( ref_left, ref_right );    }

     friend Bool operator<( const varcall& c1, const varcall& c2 )
     {    return c1.state < c2.state;    }
     friend Bool operator>( const varcall& c1, const varcall& c2 )
     {    return c1.state > c2.state;    }
     friend Bool operator==( const varcall& c1, const varcall& c2 )
     {    return c1.state == c2.state;    }
     friend Bool operator!=( const varcall& c1, const varcall& c2 )
     {    return c1.state != c2.state;    }

};

class variant {

     public:

     variant( const int pos, const String& ref, const String& alt,
          const Bool& fosmid, const Bool& confirmed, const String& friends ) 
          : pos(pos), ref(ref), alt(alt), fosmid(fosmid), 
          confirmed(confirmed), friends(friends) { }


     size_t getRefOffset() const { ForceAssert(pos>0); return pos-1; }
     size_t getRefSeqLen() const { return ref.size(); }
     bvec getAltSeq() const { return bvec(alt); }
     bool isHomozygous() const { return    (GATK_100.state=="uncalled" || GATK_100.state=="hom" || GATK_100.state=="unknown")
                                        && (GATK_250.state=="uncalled" || GATK_250.state=="hom" || GATK_250.state=="unknown")
                                        && (DISCOVAR.state=="uncalled" || DISCOVAR.state=="hom" || DISCOVAR.state=="unknown")
                                        && (CORTEX.state=="uncalled" || CORTEX.state=="hom" || CORTEX.state=="unknown") ; }

     int pos;
     String ref;
     String alt;
     Bool fosmid;
     Bool confirmed;
     String het_type;
     varcall GATK_100;
     varcall GATK_250;
     varcall DISCOVAR;
     varcall CORTEX;
     String friends;
     String vartype;

     int Stop( ) const { return pos + ref.isize( ) - 1; }

     friend Bool operator<( const variant& v1, const variant& v2 )
     {    if ( v1.pos < v2.pos ) return True;
          if ( v1.pos > v2.pos ) return False;
          if ( v1.ref < v2.ref ) return True;
          if ( v1.ref > v2.ref ) return False;
          if ( v1.alt < v2.alt ) return True;
          if ( v1.alt > v2.alt ) return False;
          if ( v1.friends < v2.friends ) return True;
          if ( v1.friends > v2.friends ) return False;
          if ( v1.fosmid < v2.fosmid ) return True;
          if ( v1.fosmid > v2.fosmid ) return False;
          if ( v1.confirmed < v2.confirmed ) return True;
          if ( v1.confirmed > v2.confirmed ) return False;
          if ( v1.GATK_100 < v2.GATK_100 ) return True;
          if ( v1.GATK_100 > v2.GATK_100 ) return False;
          if ( v1.GATK_250 < v2.GATK_250 ) return True;
          if ( v1.GATK_250 > v2.GATK_250 ) return False;
          if ( v1.CORTEX < v2.CORTEX ) return True;
          if ( v1.CORTEX > v2.CORTEX ) return False;
          if ( v1.DISCOVAR < v2.DISCOVAR ) return True;
          return False;    }

     friend ostream& operator<<( ostream& out, const variant& v )
     {    out << v.pos << "\t" << v.ref << "\t" << v.alt << "\t";
          Bool printed = False;
          if ( v.GATK_100.state != "uncalled" )
          {    out << "GATK-100=" << v.GATK_100.state;
               printed = True;    }
          if ( v.GATK_250.state != "uncalled" )
          {    if ( printed ) out << ",";
               printed = True;
               out << "GATK-250=" << v.GATK_250.state;    }
          if ( v.CORTEX.state != "uncalled" )
          {    if ( printed ) out << ",";
               printed = True;
               out << "CORTEX=" << v.CORTEX.state;    }
          if ( v.DISCOVAR.state != "uncalled" )
          {    if ( printed ) out << ",";
               printed = True;
               out << "DISCOVAR=" << v.DISCOVAR.state;    }
          if ( v.fosmid )
          {    if ( printed ) out << ",";
               printed = True;
               out << "Fosmid";    }
          if ( v.het_type != "" )
          {    if ( printed ) out << ",";
               printed = True;
               out << v.het_type;    }
          if ( v.confirmed )
          {    if ( printed ) out << ",";
               out << "confirmed";    }
          return out << "," << v.vartype << "\t" << v.friends;    }

     friend Bool VarEq( const variant& v1, const variant& v2 )
     {    return v1.pos == v2.pos && v1.ref == v2.ref && v1.alt == v2.alt;    }

};

Bool Isolated( const vec<variant>& vars, const int i, const int d, 
     const String caller )
{    for ( int j = i - 1; j >= 0; j-- )
     {    if ( vars[j].GATK_100.state == "uncalled"
               && vars[j].GATK_250.state == "uncalled"
               && vars[j].CORTEX.state == "uncalled"
               && vars[j].DISCOVAR.state == "uncalled" )
          {    continue;    }
          if ( vars[j].pos < vars[i].pos - d ) break;
          return False;    }
     for ( int j = i + 1; j < vars.isize( ); j++ )
     {    if ( vars[j].GATK_100.state == "uncalled"
               && vars[j].GATK_250.state == "uncalled"
               && vars[j].CORTEX.state == "uncalled"
               && vars[j].DISCOVAR.state == "uncalled" )
          {    continue;    }
          if ( vars[j].pos > vars[i].pos + d ) break;
          return False;    }
     return True;    }

String RefSeq( const String& s, const variant& v, const int start, 
     const int left_flank, const int right_flank )
{    String ref;
     for ( int j = v.pos - start - left_flank; j < v.pos - start; j++ )
          ref.push_back( s[j] );
     ref += v.ref;
     ForceAssertLe( v.pos - start + v.ref.isize( ) + right_flank, s.isize( ) );
     for ( int j = v.pos - start + v.ref.isize( ); 
          j < v.pos - start + v.ref.isize( ) + right_flank; j++ )
     {    ref.push_back( s[j] );    }
     return ref;    }

String AltSeq( const String& s, const variant& v, const int start, 
     const int left_flank, const int right_flank )
{    String alt;
     for ( int j = v.pos - start - left_flank; j < v.pos - start; j++ )
          alt.push_back( s[j] );
     alt += v.alt;
     for ( int j = v.pos - start + v.ref.isize( ); 
          j < v.pos - start + v.ref.isize( ) + right_flank; j++ )
     {    alt.push_back( s[j] );    }
     return alt;    }

} // close anonymous namespace

void MakeEdits( const basevector& target, const int start, const int stop,
     const vec< triple<int,String,String> >& edits, const int ei, 
     const vec<String>& euse, String& nhood )
{    nhood.clear( );
     for ( int i = start; i < edits[ei].first - 1; i++ )
          nhood.push_back( as_base( target[i] ) );
     nhood.push_back( '|' );
     for ( int i = edits[ei].first - 1; i < stop; i++ )
          nhood.push_back( as_base( target[i] ) );
     for ( int i = edits.isize( ) - 1; i >= 0; i-- )
     {    if ( euse[i] == "-" ) continue;
          int pos = edits[i].first - 1;
          String ref = edits[i].second, alt = edits[i].third;
          int npos = pos - start;
          if ( i >= ei ) npos++;
          String nhoodx(nhood);
          nhoodx.ReplaceBy( "|", "" );
          int px = nhoodx.Position(ref);
          if ( nhood.Contains( ref, npos ) )
          {    nhood = nhood.substr( 0, npos ) + alt + nhood.substr( npos 
                    + ref.isize( ), nhood.isize( ) - npos - ref.isize( ) );    }
          else
          {    int q = nhood.Position( "|" );
               String nhoodx(nhood);
               nhoodx.ReplaceBy( "|", "" );
               if ( !nhoodx.Contains( ref, npos ) )
               {    cout << "\nWell dang, looks like you found a bug!\n"
                         << "As a workaround, you might try raising wing." 
                         << endl;    }
               ForceAssert( nhoodx.Contains( ref, npos ) );
               int nposx(npos);
               if ( i >= ei ) nposx--;
               nhoodx = nhoodx.substr( 0, nposx ) + alt + nhoodx.substr( nposx
                    + ref.isize( ), nhoodx.isize( ) - nposx - ref.isize( ) );    
               q += ( alt.isize( ) - ref.isize( ) ) / 2;
               nhood = nhoodx.substr( 0, q ) + "|" 
                    + nhoodx.substr( q, nhoodx.isize( ) - q );    }    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(ID, "list of Fosmid ids or 'all' or k/n e.g. 1/4 "
          "to use the first quarter of the ids; this can be used to split "
          "computation across multiple servers; also 'parallel' is like 'all' "
          "but runs multiple assemblies at a time (and output is out of order)");
     CommandArgument_Bool_OrDefault_Doc(HET_MAP, False,
          "print map of observed heterozygous sites");
     CommandArgument_Bool_OrDefault_Doc(PAR_MAP, False,
          "print map of parental counts");
     CommandArgument_String_OrDefault_Doc(ROOT, "",
          "if specified, use this as root instead of /wga/scr4/$user/CompareVars");
     CommandArgument_Bool_OrDefault_Doc(COUNT_ONLY, False,
          "print event count instead of events");
     CommandArgument_Bool_OrDefault_Doc(HOMOPOLYMER_VERBOSE, False,
          "verbose output for homopolymer analysis");
     CommandArgument_Bool_OrDefault_Doc(FLANKS_VERBOSE, False,
          "verbose output for flanks of variants");
     CommandArgument_String_OrDefault_Doc(FUN, "",
          "list of cacheable operations, from the following list:\n"
          "     ASSEMBLE -- run DISCOVAR to assemble and call variants\n"
          "     CALL -- run DISCOVAR to call variants\n"
          "     GATK -- gather GATK and CORTEX variants\n"
          "     READS -- gather read sets\n"
          "     ALIGN -- build alignment of Fosmid reference to hg19\n"
          "     ALL -- do everything above\n"
          "     GLOBAL -- rewrite DISCOVAR variants using global vcf\n"
          "     REGIONS -- rewrite DISCOVAR variants after generating merged vcf\n"
          "The default is to do nothing.");
     CommandArgument_String_OrDefault_Doc(R, "",
          "if specified, used the given revision of LongProto");
     CommandArgument_Bool_OrDefault_Doc(PRINT_REGION_INFO, False,
          "print coordinates of Fosmid region on hg19");
     CommandArgument_String_OrDefault_Doc(FILTER_PLUS, "PASS", "filter out variants "
          "that do not contain something in this set; try also "
          "\"{PASS,99.00to,99.90to}\" or \"{PASS,99.90to}\"");
     CommandArgument_String_OrDefault_Doc(FILTER_MINUS, "", "filter out "
          "variants that contain something in this set; try \"Low\" for "
          "GATK-100");
     CommandArgument_Double_OrDefault_Doc(MIN_QUAL_SNP, 180,
          "min QUAL value for SNP VCF line to be accepted; "
          "only applies to GATK-250");
     CommandArgument_Double_OrDefault_Doc(MIN_QUAL_INDEL, 100,
          "min QUAL value for indel VCF line to be accepted; "
          "only applies to GATK-250");
     CommandArgument_Bool_OrDefault_Doc(IGNORE_NEAR_HOMOPOLYMER, False,
          "ignore Fosmid truth variants that are within five bases of a "
          "ten base homopolymer; only makes sense for sensitivity");
     EndCommandArguments;

     // Define root.

     String root = "/wga/scr4/" + Getenv( "USER" ) + "/CompareVars";
     if ( ROOT != "" ) root = ROOT;

     String merged_vcf;
     if (FUN == "GLOBAL") 
       merged_vcf  = "/wga/scr4/human_assemblies/1/v9/v9_combined.filtered.vcf";
     else if (FUN == "REGIONS") {
     // Assemble and call fosmid regions
       String regions_dir = root + "/combined_regions";
       Mkpath(regions_dir);
       merged_vcf = regions_dir + "/combined.filtered.vcf";
       String cmd = "RunRangeByRegion ROOT=" + regions_dir +
	 " TARGET=" + (ID == "all" ? "fosmids" : ID) +
	 " REGION_SIZE=50 REGION_OVERLAP=10 PARALLEL=3"
	 " > " + regions_dir + "/regions.log";
       SystemSucceed( cmd );
       cmd = "combine_region_vcfs.py --output " + merged_vcf +
	 " " + regions_dir + " > " + regions_dir + "/merge.log";
       SystemSucceed( cmd );
     }
     
     // Rewrite DISCOVAR variants using global VCF.
     
     if ( FUN == "GLOBAL" || FUN == "REGIONS")
     {    vec<int> all = AllFosmids( );
          if (FUN == "REGIONS" && ID != "all") 
	    all = vec<int>(1,ID.Int());
          vec<String> chrlist;
          for ( int c = 1; c <= 22; c++ )
               chrlist.push_back( ToString(c) );
          chrlist.push_back( "X" );
          Sort(chrlist);
          vec< vec< triple<int,int,int> > > start_stop( chrlist.size( ) );
          for ( int i = 0; i < all.isize( ); i++ )
          {    String home = root + "/" + ToString(all[i]), loc;
               int g, rstart, rstop;
               String gid;
               GetRegionInfo( ToString(all[i]), g, rstart, rstop, gid, loc );
               align q;
               BinaryReader::readFile( home + "/region.align", &q );
               int start = q.pos2( ) + 1, stop = q.Pos2( );
               int cid = BinPosition( chrlist, gid );
               ForceAssertGe( cid, 0 );
               start_stop[cid].push( start, stop, all[i] );    }
          fast_ifstream in(merged_vcf);
          String line, chr;
          int pos;
          const int flank = 100;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               istrstream iline( line.c_str( ) );
               iline >> chr >> pos;
               int cid = BinPosition( chrlist, chr );
               if ( cid < 0 ) break;
               for ( int j = 0; j < start_stop[cid].isize( ); j++ )
               {    if ( !( pos + flank >= start_stop[cid][j].first
                         && pos - flank <= start_stop[cid][j].second ) )
                    {    continue;    }
                    int id = start_stop[cid][j].third;
                    String home = root + "/" + ToString(id);
                    Ofstream( out, home + "/aaa.final.variant.filtered.vcf" );
                    out << line << "\n";
                    String chr2;
                    int pos2;
                    while(1)
                    {    getline( in, line );
                         // Note that one line is lost at end of each Fosmid,
                         // but shouldn't matter.
                         if ( in.fail( ) ) break;
                         istrstream iline( line.c_str( ) );
                         iline >> chr2 >> pos2;
                         if ( chr2 != chr ) break;
                         if ( !( pos2 + flank >= start_stop[cid][j].first
                              && pos2 - flank <= start_stop[cid][j].second ) )
                         {    break;    }
                         out << line << "\n";    }    }    }
          return 0;    }

     // Parse FILTER.

     vec<String> filter_plus, filter_minus;
     ParseStringSet( FILTER_PLUS, filter_plus );
     ParseStringSet( FILTER_MINUS, filter_minus );

     // Define genotype map.

     vec< pair< triple<String,int,int>, vec<int> > > gmap;
     vec< vec< vec<int> > > groups;
     vec< vec<int> > compat;
     DefineGenotypeMap( root, gmap, groups, compat );

     // Test for non-singleton run.

     if ( !ID.IsInt( ) )
     {    String com = command.TheCommand( );
          Bool parallel = False;
          if ( ID == "parallel" )
          {    parallel = True;
               ID = "all";    }
          RunSome( root, ID, com, COUNT_ONLY, parallel );
          return 0;    }

     // Define paths.

     String fin_dir 
          = "/wga/dev/references/Homo_sapiens/NA12878_Fosmid_Pool.regions.fin";
     String ref_dir = "/wga/scr4/bigrefs/human19";
     String lookup = ref_dir + "/genome.lookup.lookup";
     String fos_fasta = fin_dir + "/fos." + ID + ".fasta";
     String home = root + "/" + ID;
     String logs = home + "/logs";
     Mkdir777(root);
     Mkdir777(home);
     Mkdir777(logs);

     // Define search database.

     disk_map search_db( home + "/search.db" );
     vec<String> reads_id = { "g", "1", "2", "3", "dad", "dadsmom", "dadsdad", "d1", 
          "d2", "d3", "d4", "d5", "s1", "s2", "s3", "s4", "s5", "s6" };
     vec<vecbasevector> reads( reads_id.size( ) );

     // Get region info.  Define Fosmid and genome sequences.

     int g, rstart, rstop;
     String gid, loc;
     GetRegionInfo( ID, g, rstart, rstop, gid, loc );
     vecbasevector fosmid, genome;
     FetchReads( fosmid, 0, fos_fasta );
     const basevector& query = fosmid[0];
     genome.ReadOne( ref_dir + "/genome.fastb", g );
     const basevector& target = genome[0];

     // Run functions.

     RunFunctions( FUN, R, home, ref_dir, loc, ID, query, target, fos_fasta, 
          rstart, rstop );

     // Set up to track variants.

     vec<variant> vars;

     // Load and parse alignment.

     align q;
     BinaryReader::readFile( home + "/region.align", &q );
     int start = q.pos2( ) + 1, stop = q.Pos2( ), p1 = q.pos1( ), p2 = q.pos2( );
     if (PRINT_REGION_INFO)
     {    cout << "FOSMID " << ID << "\t" << gid << "\t" << start << "\t" 
               << stop << endl;    }
     const int iglen = 10;
     const int ignear = 5;
     for ( int j = 0; j < q.Nblocks( ); j++ )
     {    int g = q.Gaps(j);
          if ( g > 0 )
          {    String ref;
               for ( int l = -1; l < g; l++ )
                    ref.push_back( as_base( target[p2+l] ) );
               Bool ignore = False;
               if (IGNORE_NEAR_HOMOPOLYMER)
               {    int start = Max( 0, p1 - iglen - ignear );
                    int stop = Min( query.isize( ), p1 + iglen + ignear );
                    String s = basevector( query, start, stop - start ).ToString( );
                    if ( s.Contains( String( iglen, 'A' ) )
                         || s.Contains( String( iglen, 'C' ) )
                         || s.Contains( String( iglen, 'G' ) )
                         || s.Contains( String( iglen, 'T' ) ) )
                    {    ignore = True;    }    }
               if ( !ignore )
               {    vars.push( p2, ref, String( as_base( query[p1-1] ) ), True, 
                         False, "" );    }
               p2 += g;    }
          if ( g < 0 )
          {    String alt;
               for ( int l = -1; l < -g; l++ )
                    alt.push_back( as_base( query[p1+l] ) );
               Bool ignore = False;
               if (IGNORE_NEAR_HOMOPOLYMER)
               {    int start = Max( 0, p1 - iglen - ignear );
                    int stop = Min( query.isize( ), p1 - g + iglen + ignear );
                    String s = basevector( query, start, stop - start ).ToString( );
                    if ( s.Contains( String( iglen, 'A' ) )
                         || s.Contains( String( iglen, 'C' ) )
                         || s.Contains( String( iglen, 'G' ) )
                         || s.Contains( String( iglen, 'T' ) ) )
                    {    ignore = True;    }    }
               if ( !ignore )
               {    vars.push( p2, String( as_base( target[p2-1] ) ), alt, True, 
                         False, "" );    }
               p1 -= g;    }
          for ( int x = 0; x < q.Lengths(j); x++ )
          {    if ( query[p1] != target[p2] )
               {    Bool ignore = False;
                    if (IGNORE_NEAR_HOMOPOLYMER)
                    {    int start = Max( 0, p1 - iglen - ignear );
                         int stop = Min( query.isize( ), p1 + iglen + ignear );
                         String s 
                              = basevector( query, start, stop - start ).ToString( );
                         if ( s.Contains( String( iglen, 'A' ) )
                              || s.Contains( String( iglen, 'C' ) )
                              || s.Contains( String( iglen, 'G' ) )
                              || s.Contains( String( iglen, 'T' ) ) )
                         {    ignore = True;    }    }
                    if ( !ignore )
                    {    vars.push( p2 + 1, String( as_base( target[p2] ) ),
                              String( as_base( query[p1] ) ), True, 
                              False, "" );    }    }
               ++p1, ++p2;    }    }

     const vec<variant> fosmid_vars = vars;

     // Load and investigate homozygous variants.

     vec< triple<int,String,String> > true_homs;
     vec<VariantCall> fosmidEvents;
     fosmidEvents.reserve(vars.size());
     for ( int i = 0; i < vars.isize( ); i++ ){
          true_homs.push( vars[i].pos, vars[i].ref, vars[i].alt );
          fosmidEvents.emplace_back(vars[i].pos-1,vars[i].ref.size(),bvec(vars[i].alt),true);
     }
     /*
     cout << "Fosmid:\n"; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     for ( int i = 0; i < vars.isize( ); i++ ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     {    PRINT4( i, vars[i].pos, vars[i].ref, vars[i].alt );    } // XXXXXXXXXXXXXX
     */
     Ofstream( hout, home + "/hom_fp" );
     Ofstream( hout2, home + "/hom" );
     for ( int pass = 1; pass <= 4; pass++ )
     {    vec< triple<int,String,String> > all_homs;
          String caller;
          if ( pass == 1 ) caller = "DISCOVAR";
          if ( pass == 2 ) caller = "GATK-100";
          if ( pass == 3 ) caller = "GATK-250";
          if ( pass == 4 ) caller = "CORTEX";
          String vcf;
          if ( pass == 1 ) vcf = "aaa.final.variant.filtered.vcf";
          if ( pass == 2 ) vcf = "vcf.raw.100";
          if ( pass == 3 ) vcf = "vcf.raw.250";
          if ( pass == 4 ) vcf = "vcf.raw.CORTEX";
          Ifstream( gin, home + "/" + vcf );
          vec<String> lines;
          String line, junk;
          while(1)
          {    getline( gin, line );
               if ( gin.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;

               istrstream iline( line.c_str( ) );
               String ref, alt, gen;
               String qual,filter;
               iline >> junk >> junk >> junk >> ref >> alt >> qual >> filter 
                    >> junk >> junk >> gen;
               if ( !gen.Contains( "/" ) && !gen.Contains( "|" ) ) continue;

               Bool filtered = True;
               for ( int f = 0; f < filter_plus.isize( ); f++ )
                    if ( line.Contains( filter_plus[f] ) ) filtered = False;

               for ( int f = 0; f < filter_minus.isize( ); f++ )
                    if ( line.Contains( filter_minus[f] ) ) filtered = True;

               Bool SNP = ( ref.size( ) == 1 && alt.size( ) == 1 );
               if ( pass == 3 && SNP && qual.Double( ) < MIN_QUAL_SNP ) continue;
               if ( pass == 3 && !SNP && qual.Double( ) < MIN_QUAL_INDEL ) continue;

               if ( /*pass >= 2 &&*/ filtered ) continue;
               lines.push_back(line);    }
          int N = lines.size( );
          vec<String> chr(N), pos(N), ref(N), alt(N), gen(N);
          vec< vec<int> > alleles(N);
          vec<Bool> flaky( N, False );
          for ( int z = 0; z < N; z++ )
          {    const String& line = lines[z];
               istrstream iline( line.c_str( ) );
               iline >> chr[z] >> pos[z] >> junk >> ref[z] >> alt[z] >> junk >> junk 
                    >> junk >> junk >> gen[z];
               /*
               if ( pass == 1 )
               {    vec<String> fields, q;
                    Tokenize( gen[z], {':'}, fields );
                    Tokenize( fields[3], {','}, q );
                    q.push_front( fields[2] );
                    for ( int j = 0; j < q.isize( ); j++ )
                    {    int m = Min( (int) q[j].Int( ), 40 );
                         double p = 1.0 - pow( 10, -m/10.0 );
                         if ( p > 0 && p < cutoff2 ) flaky[z] = True;
                         if ( p > cutoff1 ) alleles[z].push_back(j);    }
               }
               */
               if ( gen[z].Contains( ":" ) ) gen[z] = gen[z].Before( ":" );    }
          int count = 0;
          vec<VariantCall> callerEvents;
          for ( int z = 0; z < N; z++ )
          {    if ( alt[z] == "." ) continue;
               if ( pos[z].Int( ) < start || pos[z].Int( ) > stop ) continue;
               
               gen[z].GlobalReplaceBy( "|", "/" );
               String h = gen[z].Before( "/" );
               String h2 = gen[z].After( "/" );
               if ( pass > 1 )
               {    if ( h == "." ){
                        if (h2==".") continue;
                        else h=h2;
                    }
               }
               vec<String> alts;
               Tokenize( alt[z], {','}, alts );
               int hi = h.Int( );
               int h2i = h2.Int( );
               
               /*
               if( pass == 1 ){
                   if( alleles[z].size() == 0 || flaky[z]){
                       continue;
                   }
                   else if( alleles[z].size() == 1){
                       if(alleles[z][0]==0) continue;
                       hi=alleles[z][0];
                       h2i=hi;
                   }
                   else{
                       hi=alleles[z][0];
                       h2i=alleles[z][1];
                   }
               }
               */
               
               if( hi == 0){
                   if(h2i!=0){
                       callerEvents.emplace_back(pos[z].Int()-1,ref[z].size(),bvec(alts[h2i-1]),false);
                   }
               }
               else{
                   if(h2i==0){
                       callerEvents.emplace_back(pos[z].Int()-1,ref[z].size(),bvec(alts[hi-1]),false);
                   }
                   else{
                       if(h2i==hi){
                           callerEvents.emplace_back(pos[z].Int()-1,ref[z].size(),bvec(alts[hi-1]),true);
                       }
                       else{
                           callerEvents.emplace_back(pos[z].Int()-1,ref[z].size(),bvec(alts[hi-1]),false);
                           callerEvents.emplace_back(pos[z].Int()-1,ref[z].size(),bvec(alts[h2i-1]),false);
                       }
                   }
               }
          }
          fosmidifyEvents<VariantCall>( target, fosmidEvents, callerEvents);
          for( const auto& entry: callerEvents){
               if( !entry.isHomozygous()) continue;
               String ALT= entry.getAltSeq().ToString();
               String REF;
               REF.reserve(entry.getRefSeqLen());
               auto begin=target.begin()+entry.getRefOffset();
               std::transform(begin,begin+entry.getRefSeqLen()
                             ,std::back_inserter(REF),BaseToCharMapper());
              
              
               while( ALT.nonempty( ) && REF.nonempty( )
                    && ALT.back( ) == REF.back( ) )
               {    ALT.resize( ALT.size( ) - 1 );
                    REF.resize( REF.size( ) - 1 );    }
               int POS = entry.getRefOffset()+1;
               while( ALT.size( ) > 1 && REF.size( ) > 1
                    && ALT[0] == REF[0] )
               {    ALT = ALT.substr( 1, ALT.size( ) - 1 );
                    REF = REF.substr( 1, REF.size( ) - 1 );
                    POS++;    }
               if ( Member( all_homs, make_triple( POS, REF, ALT ) ) ) continue;
               all_homs.push( POS, REF, ALT );
               count++;
               String h="N/A";
               PRINT6_TO( hout2, ID, caller, POS, REF, ALT, h );
               if ( !Member( true_homs, make_triple( POS, REF, ALT ) ) )
                    PRINT6_TO( hout, ID, caller, POS, REF, ALT, h );    }
          PRINT2_TO( hout, caller, count );    }

     // Load assembly variants.

     String line, junk, pos, ref, alt, prob, rest, friends, genotype;
     const int wing = 600; // manually raised to avoid asserts
     const int max_flank = 60;
     /*
     Ifstream( in, home + "/aaa.final.variant" );
     vec<String> lines;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          lines.push_back(line);    }
     for ( int i = 0; i < lines.isize( ); i++ )
     {    if ( !lines[i].Contains( "BLOCK", 0 ) ) continue;
          int j;
          for ( j = i + 1; j < lines.isize( ); j++ )
               if ( lines[j].size( ) == 0 ) break;
          vec<String> vlines;
          for ( int k = i; k < j; k++ )
               if ( lines[k].Contains( "VAR", 0 ) ) vlines.push_back( lines[k] );
          int nv = vlines.size( );
          vec< vec<String> > gen(nv);
          for ( int k = 0; k < nv; k++ )
          {    istringstream iline( vlines[k].c_str( ) );
               iline >> junk >> junk >> pos >> ref >> alt >> genotype;
               Tokenize( genotype, {'|'}, gen[k] );    }
          int nc = gen[0].size( );
          for ( int l = 0; l < nv; l++ )
          {    int m;
               vec< vec<String> > cols(nc), ucols;
               for ( m = l; m < nv; m++ )
               {    for ( int r = 0; r < nc; r++ )
                         cols[r].push_back( gen[m][r] );
                    vec< vec<String> > colsu(cols);
                    UniqueSort(colsu);
                    if ( colsu.size( ) > 2 ) break;
                    ucols = colsu;    }

               // Now we have a subblock that has either one or two genotypes across
               // it.  Next find the positions of the previous and next variants.

               String last_pos, next_pos;
               if ( l > 0 )
               {    istringstream iline( vlines[l-1].c_str( ) );
                    iline >> junk >> junk >> last_pos;    }
               else
               {    int ip = i;
                    while( ip > 0 )
                    {    ip--;
                         if ( lines[ip].Contains( "VAR", 0 ) )
                         {    istringstream iline( lines[ip].c_str( ) );
                              iline >> junk >> junk >> last_pos;
                              break;    }    }    }
               if ( m < nv )
               {    istringstream iline( vlines[m].c_str( ) );
                    iline >> junk >> junk >> next_pos;    }
               else
               {    int ip = j-1;
                    while( ip < lines.isize( ) - 1 )
                    {    ip++;
                         if ( lines[ip].Contains( "VAR", 0 ) )
                         {    istringstream iline( lines[ip].c_str( ) );
                              iline >> junk >> junk >> next_pos;
                              break;    }    }    }

               // Go through the subblock.

               if (FLANKS_VERBOSE)
               {    cout << "----------------------------------------------------"
                         << "--------------------------------\n";
                    PRINT(last_pos);    }

               vec<String> posx, refx, altx;
               for ( int z = l; z < m; z++ )
               {    istringstream iline( vlines[z].c_str( ) );
                    iline >> junk >> junk >> pos >> ref >> alt;
                    if ( vlines[z].Contains( "REF" ) )
                         ref = vlines[z].Between( "REF=", ";" );
                    if ( vlines[z].Contains( "ALT" ) )
                    {    if ( vlines[z].Contains( " FRIENDS=" ) )
                              alt = vlines[z].Between( "ALT=", " " );
                         else alt = vlines[z].After( "ALT=" );    }
                    posx.push_back(pos);
                    refx.push_back(ref), altx.push_back(alt);    }

               vec<String> ref_flanks( posx.size( ) ), alt_flanks( posx.size( ) );
               if ( ucols.size( ) == 2 )
               {    for ( int y = 0; y < posx.isize( ); y++ )
                    {    if (FLANKS_VERBOSE)
                         {    cout << "\nfinding flanks for variant at " << posx[y]
                                   << endl;    }
                         for ( int w = 0; w < ucols.isize( ); w++ )
                         {    int start = Max( 0, 
                                   (int) posx[0].After( ":" ).Int( ) - wing );
                              int stop = Min( target.isize( ), 
                                   (int) posx.back( ).After( ":" ).Int( ) + wing );
                              vec< triple<int,String,String> > edits;
                              for ( int y = 0; y < posx.isize( ); y++ )
                              {    edits.push( posx[y].After( ":" ).Int( ),
                                        refx[y], altx[y] );    }
                              String nhood;
                              MakeEdits( target, start, stop, edits, y, 
                                   ucols[w], nhood );
                              int l = nhood.Position( "|" );
                              int left = Max( 0, l - max_flank );
                              nhood = nhood.substr( left, nhood.isize( ) - left );
                              int nright 
                                   = nhood.isize( ) - nhood.Position( "|" ) - 1;
                              if ( nright > max_flank )
                              {    nhood.resize( nhood.size( ) 
                                        - ( nright - max_flank ) );    }
                              if (FLANKS_VERBOSE)
                              {    cout << "allele " << ucols[w][y]
                                        << " --> " << nhood << endl;    }
                              if ( ucols[w][y] == "+" ) alt_flanks[y] = nhood;
                              else ref_flanks[y] = nhood;    }    }
                    if (FLANKS_VERBOSE) cout << "\n";    }

               for ( int z = l; z < m; z++ )
               {    istringstream iline( vlines[z].c_str( ) );
                    iline >> junk >> junk >> pos >> ref >> alt;
                    if (FLANKS_VERBOSE) cout << vlines[z] << endl;
                    pos = pos.After( ":" );
                    rest = vlines[z].After( "PROB= " );
                    if ( rest.Contains( " FRIENDS=" ) )
                    {    friends = "FRIENDS=" + rest.After( " FRIENDS=" );
                         rest = rest.Before( " FRIENDS=" );    }
                    if ( rest.Contains( "REF" ) ) prob = rest.Before( " " );
                    else prob = rest;
                    if ( rest.Contains( "REF" ) )
                    {    ref = rest.Between( "REF=", ";" );
                         alt = rest.After( "ALT=" );    }
                    if ( prob.Contains( " ", -1 ) ) prob = prob.Before( " " );
                    const double minp = cutoff1;
                    double p1 = prob.Before( "/" ).Double( );
                    double p2 = prob.After( "/" ).Double( );
                    Bool flaky = False;
                    if ( p1 > 0 && p1 < cutoff2 ) flaky = True;
                    if ( p2 > 0 && p2 < cutoff2 ) flaky = True;
                    if (flaky) continue;
                    if ( p1 <= minp ) p1 = 0;
                    if ( p2 <= minp ) p2 = 0;
                    prob = ToString(p1) + "/" + ToString(p2);

                    if ( prob == "0/1" ) prob = "hom";
                    if ( prob == "1/1" ) prob = "het";
                    if ( prob == "0/0" ) continue; // ??????????????????????????????

                    if ( prob.Contains( "/" ) && prob.After( "/" ) == "0" ) continue;
                    vars.push( pos.Int( ), ref, alt, False, False, friends );    
                    vars.back( ).DISCOVAR.state = prob;
                    vars.back( ).DISCOVAR.alt_nhood = alt_flanks[z-l];
                    vars.back( ).DISCOVAR.ref_nhood = ref_flanks[z-l];    }
               if (FLANKS_VERBOSE) PRINT(next_pos);

               l = m - 1;    }
          i = j - 1;    }
     */

     // Load DISCOVAR variants.
     /*
     {    String caller = "DISCOVAR";
          String vcf = "aaa.final.variant.vcf";
          Ifstream( gin, home + "/" + vcf );
          vec<String> lines;
          while(1)
          {    getline( gin, line );
               if ( gin.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               lines.push_back(line);    }
          int N = lines.size( );
          vec<String> chr(N), pos(N), ref(N), alt(N), gen(N);
          vec< vec<int> > alleles(N);
          vec<Bool> flaky( N, False );
          for ( int z = 0; z < N; z++ )
          {    const String& line = lines[z];
               istrstream iline( line.c_str( ) );
               iline >> chr[z] >> pos[z] >> junk >> ref[z] >> alt[z] >> junk >> junk 
                    >> junk >> junk >> gen[z];
               {    vec<String> fields, q;
                    Tokenize( gen[z], {':'}, fields );
                    Tokenize( fields[3], {','}, q );
                    q.push_front( fields[2] );
                    for ( int j = 0; j < q.isize( ); j++ )
                    {    int m = Min( (int) q[j].Int( ), 40 );
                         double p = 1.0 - pow( 10, -m/10.0 );
                         if ( p > 0 && p < cutoff2 ) flaky[z] = True;
                         if ( p > cutoff1 ) alleles[z].push_back(j);    }    }
               if ( gen[z].Contains( ":" ) ) gen[z] = gen[z].Before( ":" );    }
          int count = 0;
          for ( int z = 0; z < N; z++ )
          {    if ( flaky[z] ) continue;
               // if ( alt[z] == "." ) continue;
               if ( pos[z].Int( ) < start || pos[z].Int( ) > stop ) continue;
               vec<String> alts;
               Tokenize( alt[z], {','}, alts );
               if ( alleles[z].empty( ) ) continue;
               if ( alleles[z].size( ) > 2 ) continue;
               if ( alleles[z].solo( ) && alleles[z][0] == 0 ) continue;
               if ( alleles[z].solo( ) )
               {    vars.push( pos[z].Int( ), ref[z], alts[ alleles[z][0] - 1 ],
                         False, False, "" );
                    vars.back( ).DISCOVAR.state = "hom";    }
               else if ( alleles[z][0] == 0 )
               {    vars.push( pos[z].Int( ), ref[z], alts[ alleles[z][1] - 1 ],
                         False, False, "" );
                    vars.back( ).DISCOVAR.state = "het";    }
               else
               {    vars.push( pos[z].Int( ), ref[z], alts[ alleles[z][0] - 1 ],
                         False, False, "" );
                    vars.back( ).DISCOVAR.state = "hom+";
                    vars.push( pos[z].Int( ), ref[z], alts[ alleles[z][1] - 1 ],
                         False, False, "" );
                    vars.back( ).DISCOVAR.state = "hom+";    }    }    }
     */
     
     {    String caller = "DISCOVAR";
          String vcf = "aaa.final.variant.filtered.vcf";
          Ifstream( gin, home + "/" + vcf );
          vec<String> lines;
          String line,chr,pos,junk,ref,alt,filter,gen;
          
          while(1)
          {    getline( gin, line );
               if ( gin.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               istrstream iline( line.c_str( ) );
               iline >> chr >> pos >> junk >> ref >> alt >> junk >> filter >> junk >> junk >> gen;
               
               if ( gen.Contains( ":" ) ) gen = gen.Before( ":" );
               gen.GlobalReplaceBy( "|", "/" );
               
               if (      !filter.Contains("PASS")
                      || alt=="."
                      || pos.Int( ) < start || pos.Int( ) > stop
                      || !gen.Contains("/")
                  ) continue;
               
               vec<String> alts;
               Tokenize( alt, {','}, alts );
               
               int h1 = gen.Before("/").Int();
               int h2 = gen.After("/").Int();
               if(h1==h2){
                   if( h1!=0) {
                       vars.push( pos.Int( ), ref, alts[ h1 - 1 ], False, False, "" );
                       vars.back( ).DISCOVAR.state = "hom";
                   }
               }
               else{
                   if(h1==0){
                       vars.push( pos.Int( ), ref, alts[ h2 - 1 ], False, False, "" );
                       vars.back( ).DISCOVAR.state = "het";
                   }
                   else if( h2==0){
                       vars.push( pos.Int( ), ref, alts[ h1 - 1 ], False, False, "" );
                       vars.back( ).DISCOVAR.state = "het";
                   }
                   else{
                       vars.push( pos.Int( ), ref, alts[ h1 - 1 ], False, False, "" );
                       vars.back( ).DISCOVAR.state = "hom+";
                       vars.push( pos.Int( ), ref, alts[ h2 - 1 ], False, False, "" );
                       vars.back( ).DISCOVAR.state = "hom+";
                   }
               }
          }
     }

     // Load GATK and CORTEX variants.

     for ( int pass = 1; pass <= 3; pass++ )
     {    int GATK_RL = ( pass == 1 ? 100 : 250 );
          String cs = "GATK-" + ToString(GATK_RL);
          if ( pass == 3 ) cs = "CORTEX";
          String tail = ( pass <= 2 ? ToString(GATK_RL) : "CORTEX" );
          Ifstream( gin, home + "/vcf.raw." + tail );
          vec<String> lines;
          String line, junk;
          while(1)
          {    getline( gin, line );
               if ( gin.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               istrstream iline( line.c_str( ) );
               String ref, alt, g;
               String qual;
               iline >> junk >> junk >> junk >> ref >> alt >> qual >> junk 
                    >> junk >> junk >> g;
               if ( g.Contains( ":" ) ) g = g.Before( ":" );
               g.GlobalReplaceBy( "|", "/" );
               if ( g.Contains( "/" ) && g.Before( "/" ) == "." ) continue;
               if ( g.Contains( "/" ) && g.After( "/" ) == "." ) continue;

               Bool SNP = ( ref.size( ) == 1 && alt.size( ) == 1 );
               if ( pass == 2 && SNP && qual.Double( ) < MIN_QUAL_SNP ) continue;
               if ( pass == 2 && !SNP && qual.Double( ) < MIN_QUAL_INDEL ) continue;

               lines.push_back(line);    }
          int N = lines.size( );
          vec<String> chr(N), pos(N), ref(N), alt(N), gen(N);
          for ( int z = 0; z < N; z++ )
          {    const String& line = lines[z];
               istrstream iline( line.c_str( ) );
               iline >> chr[z] >> pos[z] >> junk >> ref[z] >> alt[z] >> junk >> junk 
                    >> junk >> junk >> gen[z];
               if ( gen[z].Contains( ":" ) ) gen[z] = gen[z].Before( ":" );    }
          vec<String> alt_flanks(N), ref_flanks(N);
          if ( pass == 2 )
          {
          for ( int z = 0; z < N; z++ )
          {    const String& line = lines[z];

               Bool filtered = True;
               for ( int f = 0; f < filter_plus.isize( ); f++ )
                    if ( line.Contains( filter_plus[f] ) ) filtered = False;

               for ( int f = 0; f < filter_minus.isize( ); f++ )
                    if ( line.Contains( filter_minus[f] ) ) filtered = True;
               if (filtered) continue;

               int low = z, high = z;
               while ( low > 0 && chr[low-1] == chr[z] 
                    && gen[low-1].Contains( "|" ) ) 
               {    low--;    }
               while ( high < N - 1 && chr[high+1] == chr[z]
                    && gen[high+1].Contains( "|" ) ) 
               {    high++;    }
               int pos_low = 0, pos_high = pos[high].Int( ) + 10000; // sloppy!
               if ( low > 0 && chr[low-1] == chr[low] ) pos_low = pos[low-1].Int( );
               if ( high < N - 1 && chr[high+1] == chr[high] )
                    pos_high = pos[high+1].Int( );
               vec<vec<int>> h(2);
               int n = high - low + 1;
               for ( int i = low; i <= high; i++ )
               {    if ( gen[i].Contains( "|" ) )
                    {    h[0].push_back( gen[i].Before( "|" ).Int( ) );
                         h[1].push_back( gen[i].After( "|" ).Int( ) );    }
                    else
                    {    h[0].push_back( gen[i].Before( "/" ).Int( ) );
                         h[1].push_back( gen[i].After( "/" ).Int( ) );    }    }
               for ( int c = 0; c < 2; c++ )
               {    vec< triple<int,String,String> > edits;
                    vec<String> euse( high - low + 1 );
                    for ( int i = low; i <= high; i++ )
                    {    int id = h[c][i-low];
                         euse[i-low] = ( id > 0 ? "+" : "-" );
                         if ( id == 0 ) edits.push( pos[i].Int( ), ref[i], ref[i] );
                         else
                         {    vec<String> alts;
                              Tokenize( alt[i], {','}, alts );
                              edits.push( 
                                   pos[i].Int( ), ref[i], alts[id-1] );    }    }
                    int start = Max( 0, (int) pos[low].Int( ) - wing );
                    int stop = Min( target.isize( ), (int) pos[high].Int( ) + wing );
                    String nhood;
                    MakeEdits( target, start, stop, edits, z-low, euse, nhood );
                    int l = nhood.Position( "|" );
                    int left = Max( 0, l - max_flank );
                    nhood = nhood.substr( left, nhood.isize( ) - left );
                    int nright = nhood.isize( ) - nhood.Position( "|" ) - 1;
                    if ( nright > max_flank )
                    {    nhood.resize( nhood.size( ) 
                              - ( nright - max_flank ) );    }
                    if (FLANKS_VERBOSE)
                    {    cout << ( pass <= 2 ? "GATK" : "CORTEX" ) << " allele " 
                              << euse[z-low] << " --> " << nhood << endl;    }
                    if ( euse[z-low] == "+" ) alt_flanks[z] = nhood;
                    else ref_flanks[z] = nhood;    }    }
          }

          for ( int z = 0; z < N; z++ )
          {    const String& line = lines[z];
               Bool filtered = True;
               for ( int f = 0; f < filter_plus.isize( ); f++ )
                    if ( line.Contains( filter_plus[f] ) ) filtered = False;

               for ( int f = 0; f < filter_minus.isize( ); f++ )
                    if ( line.Contains( filter_minus[f] ) ) filtered = True;

               if (filtered) continue;
               if ( alt[z] == "." ) continue;
               int commas = 0;
               for ( int j = 0; j < alt[z].isize( ); j++ )
                    if ( alt[z][j] == ',' ) commas++;
               if ( gen[z].Contains( "1/1", 0 ) || gen[z].Contains( "1|1", 0 ) )
               {    vars.push( pos[z].Int( ), ref[z], alt[z], False, False, "" );
                    if ( pass == 1 ) vars.back( ).GATK_100.state = "hom";
                    else if ( pass == 2 ) vars.back( ).GATK_250.state = "hom";
                    else vars.back( ).CORTEX.state = "hom";    }
               else if (    gen[z].Contains( "0/1", 0 ) || gen[z].Contains( "0|1", 0 )
                         || gen[z].Contains( "1/0", 0 ) || gen[z].Contains( "1|0", 0 ) )
               {    vars.push( pos[z].Int( ), ref[z], alt[z], False, False, "" );
                    if ( pass == 1 ) vars.back( ).GATK_100.state = "het";
                    else if ( pass == 2 ) vars.back( ).GATK_250.state = "het";
                    else vars.back( ).CORTEX.state = "het";    }
               else if (    (    gen[z].Contains("1/2",0) || gen[z].Contains("1|2",0)
                              || gen[z].Contains("2/1",0) || gen[z].Contains("2|1",0)
                            )
                         && commas == 1 )
               {    vars.push( pos[z].Int( ), ref[z], alt[z].Before( "," ), 
                         False, False, "" );
                    if ( pass == 1 ) vars.back( ).GATK_100.state = "hom+";
                    else if ( pass == 2 ) vars.back( ).GATK_250.state = "hom+";
                    else vars.back( ).CORTEX.state = "hom+";
                    vars.push( pos[z].Int( ), ref[z], alt[z].After( "," ), 
                         False, False, "" );
                    if ( pass == 1 ) vars.back( ).GATK_100.state = "hom+";
                    else if ( pass == 2 ) vars.back( ).GATK_250.state = "hom+";
                    else vars.back( ).CORTEX.state = "hom+";    }
               else  if ( ! ( gen[z].Contains("0/0",0) || gen[z].Contains("0|0",0) ) )
               {    vars.push( pos[z].Int( ), ref[z], alt[z], False, False, "" );
                    if ( pass == 1 ) vars.back( ).GATK_100.state = "unknown";
                    else if ( pass == 2 ) vars.back( ).GATK_250.state = "unknown";
                    else vars.back( ).CORTEX.state = "unknown";    }    
               if ( pass == 2 )
               {    vars.back( ).GATK_250.alt_nhood = alt_flanks[z];
                    vars.back( ).GATK_250.ref_nhood = ref_flanks[z];    }    }    }

     // Simplify ref/alt description by trimming excess bases.

     for ( int i = 0; i < vars.isize( ); i++ )
     {    variant& v = vars[i];
          String &r = v.ref, &a = v.alt;
          while ( r.size( ) > 1 && a.size( ) > 1 && r.back( ) == a.back( ) )
          {    r.resize( r.size( ) - 1 );
               a.resize( a.size( ) - 1 );    }    }

     // shift and shrink the ref/alt calls such that there's at most one matching base counting from the left.
     for( auto& entry: vars){
         auto& ref=entry.ref;
         auto& alt=entry.alt;
         int64_t first_diff=0;
         for( const int64_t upper=min(ref.size(),alt.size()) ;first_diff<upper && ref[first_diff]==alt[first_diff] ;++first_diff){ }
         const int64_t real_start = first_diff-1;
         if( real_start > 0){
             ref = ref.substr(real_start,ref.size()-real_start);
             alt = alt.substr(real_start,alt.size()-real_start);
             entry.pos+=real_start;
         }
         if( ref.size()==alt.size() && ref.size()==2 && ref[0]==alt[0]){
             ref = ref.substr(1,1);
             alt = alt.substr(1,1);
             ++entry.pos;
         }
     }
     
     // Apply Ted's transformation per sample
     {
         vec<variant> v100,v250,vC,vD;
         v100.reserve(vars.size()); v250.reserve(vars.size()); vC.reserve(vars.size()); vD.reserve(vars.size());
         for( const auto& entry : vars){
             if( entry.GATK_100.state != "uncalled"){
                 v100.push_back(entry);
                 v100.back().GATK_250.clear();
                 v100.back().CORTEX.clear();
                 v100.back().DISCOVAR.clear();
             }
             if( entry.GATK_250.state != "uncalled"){
                 v250.push_back(entry);
                 v250.back().GATK_100.clear();
                 v250.back().DISCOVAR.clear();
                 v250.back().CORTEX.clear();
             }
             if( entry.CORTEX.state != "uncalled"){
                 vC.push_back(entry);
                 vC.back().GATK_100.clear();
                 vC.back().GATK_250.clear();
                 vC.back().DISCOVAR.clear();
             }
             if( entry.DISCOVAR.state != "uncalled"){
                 vD.push_back(entry);
                 vD.back().GATK_100.clear();
                 vD.back().GATK_250.clear();
                 vD.back().CORTEX.clear();
             }
         }
         auto loc_fosmid_vars=fosmid_vars;
         for( auto& entry: loc_fosmid_vars){
             entry.CORTEX.clear(); entry.DISCOVAR.clear(); entry.GATK_100.clear(); entry.GATK_250.clear();
             entry.GATK_100.state="hom";
         }
         fosmidifyEvents<variant>( target, loc_fosmid_vars, v100);//,100,true);
         for( auto& entry: loc_fosmid_vars){
             entry.CORTEX.clear(); entry.DISCOVAR.clear(); entry.GATK_100.clear(); entry.GATK_250.clear();
             entry.GATK_250.state="hom";
         }
         fosmidifyEvents<variant>( target, loc_fosmid_vars, v250);//,100,true);
         for( auto& entry: loc_fosmid_vars){
             entry.CORTEX.clear(); entry.DISCOVAR.clear(); entry.GATK_100.clear(); entry.GATK_250.clear();
             entry.CORTEX.state="hom";
         }
         fosmidifyEvents<variant>( target, loc_fosmid_vars, vC);//,100,true);
         for( auto& entry: loc_fosmid_vars){
             entry.CORTEX.clear(); entry.DISCOVAR.clear(); entry.GATK_100.clear(); entry.GATK_250.clear();
             entry.DISCOVAR.state="hom";
         }
         fosmidifyEvents<variant>( target, loc_fosmid_vars, vD);//,100,true);
         vars.clear();
         vars.reserve( v100.size()+ v250.size() + vC.size() + vD.size());
         vars.insert(vars.end(),fosmid_vars.begin(),fosmid_vars.end());
         vars.insert(vars.end(),v100.begin(),v100.end());
         vars.insert(vars.end(),v250.begin(),v250.end());
         vars.insert(vars.end(),vC.begin(),vC.end());
         vars.insert(vars.end(),vD.begin(),vD.end());
     }
     for( auto& entry:vars){
         if( entry.DISCOVAR.state == "hom+") entry.DISCOVAR.state="hom";
         if( entry.GATK_100.state == "hom+") entry.GATK_100.state="hom";
         if( entry.GATK_250.state == "hom+") entry.GATK_250.state="hom";
         if( entry.CORTEX.state == "hom+") entry.CORTEX.state="hom";
     }

     // Look for equivalent variants.

     Sort(vars);
     vec<variant> vars0;
     for ( int i = 0; i < vars.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < vars.isize( ); j++ )
          {    if ( vars[j].pos != vars[i].pos ) break;
               if ( vars[j].ref != vars[i].ref ) break;
               if ( vars[j].alt != vars[i].alt ) break;    }
          if ( vars[i].pos < start || vars[i].pos > stop ) continue;
          vars0.push( vars[i].pos, vars[i].ref, vars[i].alt, False, False, "" );
          i = j - 1;    }
     String s =  basevector( target, start - 1, stop - start + 1 ).ToString( );
     vec<String> xs( vars0.size( ) );
     for ( int i = 0; i < vars0.isize( ); i++ )
     {    const variant& v = vars0[i];
          for ( int j = 0; j < v.pos - start; j++ )
               xs[i].push_back( s[j] );
          xs[i] += v.alt;
          for ( int j = v.pos - start + v.ref.isize( ); j < s.isize( ); j++ )
               xs[i].push_back( s[j] );    }
     vec<int> ids( vars0.size( ), vec<int>::IDENTITY );
     SortSync( xs, ids );
     for ( int i = 0; i < xs.isize( ); i++ )
     {    int j = xs.NextDiff(i);
          for ( int k = i + 1; k < j; k++ )
          {    for ( int l = 0; l < vars.isize( ); l++ )
               {    if ( VarEq( vars[l], vars0[ ids[k] ] ) )
                    {    vars[l].pos = vars0[ ids[i] ].pos;
                         vars[l].ref = vars0[ ids[i] ].ref;
                         vars[l].alt = vars0[ ids[i] ].alt;    }    }    }
          i = j - 1;    }

     // Combine variants.

     Sort(vars);
     vec<variant> varsx;
     for ( int i = 0; i < vars.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < vars.isize( ); j++ )
          {    if ( vars[j].pos != vars[i].pos ) break;
               if ( vars[j].ref != vars[i].ref ) break;
               if ( vars[j].alt != vars[i].alt ) break;    }
          vec<String> f;
          vec<Bool> fo;
          for ( int k = i; k < j; k++ )
          {    fo.push_back( vars[k].fosmid );
               f.push_back( vars[k].friends );    }
          UniqueSort(f); // ????????????????????????????????????????????????????????
          String friends;
          for ( int k = 0; k < f.isize( ); k++ )
          {    if ( k > 0 ) friends += ",";
               friends += f[k];    }
          varsx.push( vars[i].pos, vars[i].ref, vars[i].alt, Sum(fo) >= 1, False,
               friends );
          for ( int k = i; k < j; k++ )
          {    if ( vars[k].GATK_100.state != "uncalled" )
                    varsx.back( ).GATK_100 = vars[k].GATK_100;
               if ( vars[k].GATK_250.state != "uncalled" )
                    varsx.back( ).GATK_250 = vars[k].GATK_250;
               if ( vars[k].CORTEX.state != "uncalled" )
                    varsx.back( ).CORTEX = vars[k].CORTEX;
               if ( vars[k].DISCOVAR.state != "uncalled" )
                    varsx.back( ).DISCOVAR = vars[k].DISCOVAR;    }
          i = j - 1;    }

     // Label variants by type.

     for ( int i = 0; i < varsx.isize( ); i++ )
     {    variant& v = varsx[i];
          int r = v.ref.size( ), a = v.alt.size( );
          if ( r == 1 && a == 1 ) v.vartype = "c=sub";
          else if ( r == 1 && a-1 == 1 ) v.vartype = "c=ins-1";
          else if ( r == 1 && a-1 >= 2 && a-1 <= 10 ) v.vartype = "c=ins-2-10";
          else if ( r == 1 && a-1 >= 11 && a-1 <= 100 ) v.vartype = "c=ins-11-100";
          else if ( r == 1 && a-1 > 100 ) v.vartype = "c=ins-gt-100";
          else if ( a == 1 && r-1 == 1 ) v.vartype = "c=del-1";
          else if ( a == 1 && r-1 >= 2 && r-1 <= 10 ) v.vartype = "c=del-2-10";
          else if ( a == 1 && r-1 >= 11 && r-1 <= 100 ) v.vartype = "c=del-11-100";
          else if ( a == 1 && r-1 > 100 ) v.vartype = "c=del-gt-100";
          else v.vartype = "unknown";    }

     // If DISCOVAR=het or GATK=het, and isolated, test for confirmed.

     for ( int i = 0; i < varsx.isize( ); i++ )
     {    variant& v = varsx[i];
          if ( v.pos < start || v.pos > stop ) continue;

          vec<String> callsets;
          callsets.push_back( "GATK-100", "GATK-250", "CORTEX", "DISCOVAR" );
          for ( int c = 0; c < callsets.isize( ); c++ )
          {    const String& cs = callsets[c];
               if ( c == 0 && v.GATK_100.state != "het" ) continue;
               if ( c == 1 && v.GATK_250.state != "het" ) continue;
               if ( c == 2 && v.CORTEX.state != "het" ) continue;
               if ( c == 3 && v.DISCOVAR.state != "het" ) continue;
               if (v.confirmed) continue;
               for ( int pass = 1; pass <= 3; pass++ )
               {    if (v.confirmed) break;
                    int flank = 10;
                    //if ( pass == 1 ) flank = 10;
                    if ( pass == 2 ) flank = 20;
                    if ( pass == 3 ) flank = 50;
                    // int flank = ( pass == 1 ? 10 : 20 );
                    // ---->

               int n = std::max( {v.ref.isize( ), v.alt.isize( ), flank} );
               Bool isolated = Isolated( varsx, i, n, cs );

               Bool have_flanks = False;
               if ( c == 1 && v.GATK_250.state == "het" 
                    && v.GATK_250.alt_nhood.Contains( "|" )
                    && v.GATK_250.ref_nhood.Contains( "|" )
                    && v.GATK_250.alt_nhood.Before( "|" ).isize( ) >= flank
                    && v.GATK_250.alt_nhood.After( "|" ).isize( ) >= flank
                    && v.GATK_250.ref_nhood.Before( "|" ).isize( ) >= flank
                    && v.GATK_250.ref_nhood.After( "|" ).isize( ) >= flank )
               {    have_flanks = True;
                    isolated = True;    }
               if ( c == 2 && v.DISCOVAR.state == "het" 
                    && v.DISCOVAR.alt_nhood.Contains( "|" )
                    && v.DISCOVAR.ref_nhood.Contains( "|" )
                    && v.DISCOVAR.alt_nhood.Before( "|" ).isize( ) >= flank
                    && v.DISCOVAR.alt_nhood.After( "|" ).isize( ) >= flank
                    && v.DISCOVAR.ref_nhood.Before( "|" ).isize( ) >= flank
                    && v.DISCOVAR.ref_nhood.After( "|" ).isize( ) >= flank )
               {    have_flanks = True;
                    isolated = True;    }

               if ( !isolated ) continue;

               if ( v.alt.Contains( "," ) ) continue;
               if ( v.pos - start - flank < 0 ) continue;
               if ( v.pos - start + v.ref.isize( ) + flank > s.isize( ) ) continue;

               String ref = RefSeq( s, v, start, flank, flank );
               String alt = AltSeq( s, v, start, flank, flank );
               if (have_flanks)
               {    if ( c == 1 ) v.GATK_250.BuildMini( flank, alt, ref );
                    else v.DISCOVAR.BuildMini( flank, alt, ref );    }

               int ref1 = Search( search_db, reads, reads_id, home, ref, "1" ); 
               int alt1 = Search( search_db, reads, reads_id, home, alt, "1" );
               int ref2 = Search( search_db, reads, reads_id, home, ref, "2" ); 
               int alt2 = Search( search_db, reads, reads_id, home, alt, "2" );
               int ref3 = Search( search_db, reads, reads_id, home, ref, "3" ); 
               int alt3 = Search( search_db, reads, reads_id, home, alt, "3" );
               const int minc = 8;
               if ( ref1 >= minc && alt1 >= minc )
               {    if ( ( ref2 == 0 && alt2 >= minc )
                         || ( alt3 >= minc && ref3 == 0 ) 
                         ||
                         (ref2 >= minc && alt2 >= minc && ref3 >= minc && alt3 == 0)
                         ||
                         (ref2 >= minc && alt2 == 0 && ref3 >= minc && alt3 >= minc)
                              )
                    {    v.confirmed = True;    }    }    
               if (PAR_MAP)
               {    cout << cs << "\t";
                    cout << v.pos << "\t";
                    cout << "flank=" << flank << "\t";
                    cout << varsx[i] << endl;
                    PRINT6( ref1, alt1, ref2, alt2, ref3, alt3 );    }
               if ( v.het_type == "" )
               {    const int xminc = 10;
                    if ( alt2 >= xminc && alt3 == 0 && ref3 >= xminc )
                         v.het_type = "maternal";
                    if ( alt3 >= xminc && alt2 == 0 && ref2 >= xminc )
                         v.het_type = "paternal";    
                    if ( ref2 == 0 && alt2 >= xminc && ref3 >= xminc )
                         v.het_type = "maternal";
                    if ( ref3 == 0 && alt3 >= xminc && ref2 >= xminc )
                         v.het_type = "paternal";    }    }    }    }

     // Similar test, for the hom case.

     for ( int i = 0; i < varsx.isize( ); i++ )
     {    variant& v = varsx[i];
          if ( v.pos < start || v.pos > stop ) continue;

          vec<String> callsets;
          callsets.push_back( "GATK-100", "GATK-250", "CORTEX", "DISCOVAR" );
          for ( int c = 0; c < callsets.isize( ); c++ )
          {    const String& cs = callsets[c];
               if ( c == 0 && v.GATK_100.state != "hom" ) continue;
               if ( c == 1 && v.GATK_250.state != "hom" ) continue;
               if ( c == 1 && v.CORTEX.state != "hom" ) continue;
               if ( c == 2 && v.DISCOVAR.state != "hom" ) continue;
               if (v.confirmed) continue;
               for ( int pass = 1; pass <= 3; pass++ )
               {    if (v.confirmed) break;
                    int flank = 10;
                    //if ( pass == 1 ) flank = 10;
                    if ( pass == 2 ) flank = 20;
                    if ( pass == 3 ) flank = 50;
                    // ---->

               if ( c != 1 && c != 2 ) continue;

               Bool have_flanks = False;
               if ( c == 1 && v.GATK_250.state == "hom" 
                    && v.GATK_250.alt_nhood.Contains( "|" )
                    && v.GATK_250.ref_nhood.Contains( "|" )
                    && v.GATK_250.alt_nhood.Before( "|" ).isize( ) >= flank
                    && v.GATK_250.alt_nhood.After( "|" ).isize( ) >= flank
                    && v.GATK_250.ref_nhood.Before( "|" ).isize( ) >= flank
                    && v.GATK_250.ref_nhood.After( "|" ).isize( ) >= flank )
               {    have_flanks = True;    }
               if ( c == 2 && v.DISCOVAR.state == "hom" 
                    && v.DISCOVAR.alt_nhood.Contains( "|" )
                    && v.DISCOVAR.ref_nhood.Contains( "|" )
                    && v.DISCOVAR.alt_nhood.Before( "|" ).isize( ) >= flank
                    && v.DISCOVAR.alt_nhood.After( "|" ).isize( ) >= flank
                    && v.DISCOVAR.ref_nhood.Before( "|" ).isize( ) >= flank
                    && v.DISCOVAR.ref_nhood.After( "|" ).isize( ) >= flank )
               {    have_flanks = True;    }
               if ( !have_flanks ) continue;

               String ref, alt;
               if ( c == 1 ) v.GATK_250.BuildMini( flank, alt, ref );
               else v.DISCOVAR.BuildMini( flank, alt, ref );

               int ref1 = Search( search_db, reads, reads_id, home, ref, "1" ); 
               int alt1 = Search( search_db, reads, reads_id, home, alt, "1" );
               int ref2 = Search( search_db, reads, reads_id, home, ref, "2" ); 
               int alt2 = Search( search_db, reads, reads_id, home, alt, "2" );
               int ref3 = Search( search_db, reads, reads_id, home, ref, "3" ); 
               int alt3 = Search( search_db, reads, reads_id, home, alt, "3" );
               const int minc = 8;
               if ( ref1 >= minc && alt1 >= minc )
               {    if ( ( ref2 == 0 && alt2 >= minc )
                         || ( alt3 >= minc && ref3 == 0 ) 
                         ||
                         (ref2 >= minc && alt2 >= minc && ref3 >= minc && alt3 == 0)
                         ||
                         (ref2 >= minc && alt2 == 0 && ref3 >= minc && alt3 >= minc)
                              )
                    {    v.confirmed = True;    }    }    
               if (PAR_MAP)
               {    cout << cs << "\t";
                    cout << v.pos << "\t";
                    cout << "flank=" << flank << "\t";
                    cout << varsx[i] << endl;
                    PRINT6( ref1, alt1, ref2, alt2, ref3, alt3 );    }    }    }    }

     // Create het map.  We accept genotypes either if they match the genotype
     // map, or are very common (seen 5+ times).  There are various reasons why
     // the latter might be needed.  For example, as in #104, it is probably the
     // case that dad's cell line has LOH in the region (see e.g. positions
     // 20:52431468 and 20:52442912, where dad should inherit both alleles).

     if (HET_MAP) cout << "HET MAP\n" << 0 << "\t" << start << "\tstart" << endl;
     vec<int> gts_het, gts_hom;
     for ( int xpass = 1; xpass <= 2; xpass++ )
     {    if ( xpass == 2 )
          {    Sort(gts_het), Sort(gts_hom);
               vec<int> gts_het2, gts_hom2;
               const int min_gt = 5;
               for ( int i = 0; i < gts_het.isize( ); i++ )
               {    int j = gts_het.NextDiff(i);
                    if ( j - i >= min_gt ) gts_het2.push_back( gts_het[i] );
                    i = j - 1;    }
               for ( int i = 0; i < gts_hom.isize( ); i++ )
               {    int j = gts_hom.NextDiff(i);
                    if ( j - i >= min_gt ) gts_hom2.push_back( gts_hom[i] );
                    i = j - 1;    }
               gts_het = gts_het2, gts_hom = gts_hom2;    }
          for ( int i = 0; i < varsx.isize( ); i++ )
          {    variant v = varsx[i];
               if ( v.pos < start || v.pos > stop ) continue;
               if ( v.GATK_100.state != "het" && v.GATK_250.state != "het"
                    && v.DISCOVAR.state != "het" && v.CORTEX.state != "het" )
               {    continue;    }
               v.pos -= start;
               for ( int zpass = 1; zpass <= 3; zpass++ )
               {    int flank;
                    if ( zpass == 1 ) flank = 20;
                    if ( zpass == 2 ) flank = 30;
                    if ( zpass == 3 ) flank = 40;
                    
                    int left = ( i == 0 ? v.pos 
                         : Min( v.pos, varsx[i].pos - varsx[i-1].pos - 1 ) );
                    int n = varsx[i].ref.size( );
                    int right = ( i == varsx.isize( ) - 1 
                         ? s.isize( ) - v.pos - (n-1)
                         : Min( s.isize( ) - v.pos - (n-1), 
                         varsx[i+1].pos - varsx[i].pos - 1 ) );
     
                    // Turned off but seems like a bugfix.
                    // left = Max( 0, left );
                    // right = Max( 0, right );

                    left = Min( left, flank );
                    right = Min( right, flank );
                    Bool isolated = ( left + right >= flank );
                    Bool have_flanks1 = False, have_flanks2 = False;
                    if ( v.GATK_250.state == "het" 
                         && v.GATK_250.alt_nhood.Contains( "|" )
                         && v.GATK_250.ref_nhood.Contains( "|" )
                         && v.GATK_250.alt_nhood.Before( "|" ).isize( ) >= flank
                         && v.GATK_250.alt_nhood.After( "|" ).isize( ) >= flank
                         && v.GATK_250.ref_nhood.Before( "|" ).isize( ) >= flank
                         && v.GATK_250.ref_nhood.After( "|" ).isize( ) >= flank )
                    {    have_flanks1 = True;
                         isolated = True;    }
                    if ( v.DISCOVAR.state == "het" 
                         && v.DISCOVAR.alt_nhood.Contains( "|" )
                         && v.DISCOVAR.ref_nhood.Contains( "|" )
                         && v.DISCOVAR.alt_nhood.Before( "|" ).isize( ) >= flank
                         && v.DISCOVAR.alt_nhood.After( "|" ).isize( ) >= flank
                         && v.DISCOVAR.ref_nhood.Before( "|" ).isize( ) >= flank
                         && v.DISCOVAR.ref_nhood.After( "|" ).isize( ) >= flank )
                    {    have_flanks2 = True;
                         isolated = True;    }
     
                    if ( HET_MAP && xpass == 2 )
                    {    cout << v.pos << "\t";
                         cout << "flank=" << flank << "\t";
                         cout << varsx[i] << ( isolated ? " isolated" : "" ) 
                              << endl;    }
                    if ( !isolated ) continue;
                    if ( !( varsx[i].pos - start + varsx[i].ref.isize( ) + flank 
                         <= s.isize( ) ) )
                    {    continue;    }
     
                    vec<int> hets, homs;
                    for ( int l = 0; l < gmap.isize( ); l++ )
                    {    if ( gmap[l].first.first != "chr" + gid ) continue;
                         if ( varsx[i].pos-1 < gmap[l].first.second ) continue;
                         if ( varsx[i].pos-1 >= gmap[l].first.third ) continue;
                         hets = gmap[l].second;
                         const vec<vec<int>>& x0 = groups[ gmap[l].second[0] ];
                         const vec<vec<int>>& x1 = groups[ gmap[l].second[1] ];
                         vec<int> s1, s2, s3;
                         s1 = x0[2];
                         s1.append( x1[1] );
                         s2 = x0[1];
                         s2.append( x1[2] );
                         Sort(s1), Sort(s2);
                         homs.push_back( GenotypeId( s1, s2, s3 ) );
                         homs.push_back( GenotypeId( s2, s1, s3 ) );
                         s1 = x0[2];
                         s1.append( x1[2] );
                         s2 = x0[1];
                         s2.append( x1[1] );
                         Sort(s1), Sort(s2);
                         homs.push_back( GenotypeId( s1, s2, s3 ) );
                         homs.push_back( GenotypeId( s2, s1, s3 ) );
                         UniqueSort(homs);
                         if ( HET_MAP && xpass == 2 )
                         {    cout << "if dad is het, genotype id should be in "
                                   << "{" << printSeq(hets) << "}\n";
                              cout << "if dad is hom, genotype id should be in "
                                   << "{" << printSeq(homs) << "}\n";    }
                         break;    }

                    // For indel, define period.

                    String indel;
                    if ( v.ref.size( ) > 1 && v.alt.size( ) == 1 )
                         indel = v.ref.substr( 1, v.ref.isize( ) - 1 );
                    if ( v.alt.size( ) > 1 && v.ref.size( ) == 1 )
                         indel = v.alt.substr( 1, v.alt.isize( ) - 1 );
                    int period;
                    for ( period = 1; period < indel.isize( ); period++ )
                    {    if ( indel.isize( ) % period != 0 ) continue;
                         int nc = indel.isize( ) / period;
                         String x;
                         for ( int l = 0; l < nc; l++ )
                              x += indel.substr( 0, period );
                         if ( x == indel ) break;    }
                    String repeat; 
                    if ( indel.size( ) > 0 ) repeat = indel.substr( 0, period );
                    // PRINT(repeat);
                    if ( indel.size( ) > 0 && ( have_flanks1 || have_flanks2 ) )
                    {    String ref_query; 
                         if (have_flanks2) ref_query = v.DISCOVAR.ref_nhood;
                         else ref_query = v.GATK_250.ref_nhood;
                         int p = ref_query.Position( "|" );
                         int rf;
                         for ( rf = 0; ; rf++ )
                         {    if ( p + 2 + rf >= ref_query.isize( ) ) break;
                              if ( ref_query[ p + 2 + rf ] != repeat[ rf % period ] )
                                   break;    }
                         // PRINT(rf);    
                              }

                    // Define query sequences.

                    int left_flank = left, right_flank = right;
                    String ref_query = RefSeq( 
                         s, varsx[i], start, left_flank, right_flank );
                    if (have_flanks1 || have_flanks2) 
                    {    int left_flank = flank, right_flank = flank;
                         if (have_flanks2) ref_query = v.DISCOVAR.ref_nhood;
                         else ref_query = v.GATK_250.ref_nhood;
                         int p = ref_query.Position( "|" );
                         ref_query = ref_query.substr( 
                              p - left_flank, left_flank + right_flank + 1 );
                         ref_query.ReplaceBy( "|", "" );    }
                    String alt_query = 
                         AltSeq( s, varsx[i], start, left_flank, right_flank );
                    if (have_flanks1 || have_flanks2) 
                    {    int left_flank = flank, right_flank = flank;
                         if (have_flanks2) alt_query = v.DISCOVAR.alt_nhood;
                         else alt_query = v.GATK_250.alt_nhood;
                         int p = alt_query.Position( "|" );
                         alt_query = alt_query.substr( 
                              p - left_flank, left_flank + right_flank + 1 );
                         alt_query.ReplaceBy( "|", "" );    }

                    // Search.

                    vec<int> s1, s2, s3, s4;
                    int count = 0;
                    vec< pair<int,int> > counts;
                    for ( auto p : { "dad", "dadsmom", "dadsdad", "d1", "d2", 
                         "d3", "d4", "d5", "s1", "s2", "s3", "s4", "s5", 
                         "s6" } )
                    {    
                         // Search for sequences.
     
                         int ref = Search( 
                              search_db, reads, reads_id, home, ref_query, p );
                         int alt = Search( 
                              search_db, reads, reads_id, home, alt_query, p );
                         if ( HET_MAP && xpass == 2 ) PRINT3( p, ref, alt );    
                         counts.push( ref, alt );
                         if ( count >= 3 )
                         {    if ( ref > 0 && alt > 0 ) s1.push_back(count-3);
                              else if ( ref > 0 ) s2.push_back(count-3);
                              else if ( alt > 0 ) s3.push_back(count-3);
                              else s4.push_back(count-3);    }
                         count++;    }    
                    if ( s4.nonempty( ) ) continue;
          
                    Bool confirmed = False;
                    for ( int pass = 1; pass <= 2; pass++ )
                    {    if ( pass == 2 )
                         {    Bool changed = False;
                              for ( int l = 3; l < counts.isize( ); l++ )
                              {    if ( counts[l].first == 1 )
                                   {    counts[l].first = 0;
                                        changed = True;    }
                                   if ( counts[l].second == 1 )
                                   {    counts[l].second = 0;
                                        changed = True;    }    }
                              if ( !changed ) break;
                              s1.clear( ), s2.clear( ), s3.clear( ), s4.clear( );
                              for ( int l = 3; l < counts.isize( ); l++ )
                              {    if ( counts[l].first > 0 && counts[l].second > 0 )
                                        s1.push_back(l-3);
                                   else if ( counts[l].first > 0 ) s2.push_back(l-3);
                                   else if ( counts[l].second > 0 ) 
                                   {    s3.push_back(l-3);    }
                                   else s4.push_back(l-3);    }
                              if ( s4.nonempty( ) ) break;    }
                         int gt = GenotypeId( s1, s2, s3 );
                         if ( HET_MAP && xpass == 2 ) 
                         {    cout << "see genotype id " << gt 
                                   << " from " << ( v.fosmid ? "Fosmid" : "!Fosmid" )
                                   << endl;    }
                         Bool het = ( counts[0].first > 1 && counts[0].second > 1 );
                         Bool hom = ( counts[0].first > 1 ^ counts[0].second > 1 );
                         Bool z1 = False, nz1 = False, z2 = False, nz2 = False;
                         for ( int l = 3; l < counts.isize( ); l++ )
                         {    if ( counts[l].first == 0 ) z1 = True;
                              if ( counts[l].first > 0 ) nz1 = True;
                              if ( counts[l].second == 0 ) z2 = True;
                              if ( counts[l].second > 0 ) nz2 = True;    }
                         Bool solid1 = ( !z1 || !nz1 ), solid2 = ( !z2 || !nz2 );
                         Bool uninformative = ( solid1 && solid2 );
                         if ( !uninformative && xpass == 1 && het ) 
                              gts_het.push_back(gt);
                         if ( !uninformative && xpass == 1 && hom ) 
                              gts_hom.push_back(gt);
                         if ( het && BinMember(hets, gt) && !uninformative ) 
                              confirmed = True;
                         if ( hom && BinMember(homs, gt) && !uninformative ) 
                              confirmed = True;
                         if ( het && BinMember( gts_het, gt ) ) confirmed = True;
                         if ( hom && BinMember( gts_hom, gt ) ) confirmed = True;
                         if (confirmed) break;    }

                    if ( confirmed && xpass == 2 )
                    {    varsx[i].confirmed = True;
                         if (HET_MAP) cout << "confirmed\n";    }    }    }    }
     if (HET_MAP)
          cout << stop - start << "\t" << stop << "\t" << "stop\n\n" << endl;

     // Examine homopolymer indels.

     vec<String> to_search;
     for ( int pass = 1; pass <= 2; pass++ )
     {    String dataset = "g";
          vecbasevector b;
          vec<String> keys;
          for ( int i = 0; i < to_search.isize( ); i++ )
          {    String query = to_search[i];
               String key = ToString(dataset) + "_" + query;
               if ( search_db.Defined(key) ) continue;
               b.push_back( basevector(query) ), keys.push_back(key);    }
          if ( b.size( ) > 0 )
          {    vec<look_align> aligns;
               PerfectLookup( 12, b,
                    "/wga/scr4/bigrefs/human19/genome.lookup.lookup", 
                    aligns, FW_OR_RC );
               UniqueSort(aligns);
               vec<int> counts( b.size( ) );
               for ( int i = 0; i < aligns.isize( ); i++ )
                    counts[ aligns[i].query_id ]++;
               for ( int i = 0; i < keys.isize( ); i++ )
               {    if ( !search_db.Defined( keys[i] ) )
                         search_db.Set( keys[i], counts[i] );    }    }

          vec< vec<int> > hp( varsx.size( ), vec<int>(4,0) );
          for ( int i = 0; i < varsx.isize( ); i++ )
          {    variant v = varsx[i];
               Bool indel = ( v.ref.size( ) > 1 || v.alt.size( ) > 1 );
               if ( !indel ) continue;
               if ( !( v.ref.size( ) == 1 || v.alt.size( ) == 1 ) ) continue;
               String ind;
               if ( v.ref.size( ) > 1 ) ind = v.ref.substr( 1, v.ref.isize( ) - 1 );
               if ( v.alt.size( ) > 1 ) ind = v.alt.substr( 1, v.alt.isize( ) - 1 );
               vec<char> indx;
               for ( int l = 0; l < ind.isize( ); l++ )
                    indx.push_back( ind[l] );
               UniqueSort(indx);
               if ( indx.size( ) != 1 ) continue;
               hp[i][ as_char( indx[0] ) ]
                    = ( v.ref.size( ) > 1 ? -( v.ref.size( ) - 1 )
                    : v.alt.size( ) - 1 );    }

          for ( int i = 0; i < varsx.isize( ); i++ )
          {    variant v = varsx[i];
               if ( v.pos < start || v.pos > stop ) continue;
               v.pos -= start;
               if (v.confirmed) continue;
               Bool indel = ( v.ref.size( ) > 1 || v.alt.size( ) > 1 );
               if ( !indel ) continue;
               if ( !( v.ref.size( ) == 1 || v.alt.size( ) == 1 ) ) continue;
               String ind;
               if ( v.ref.size( ) > 1 ) ind = v.ref.substr( 1, v.ref.isize( ) - 1 );
               if ( v.alt.size( ) > 1 ) ind = v.alt.substr( 1, v.alt.isize( ) - 1 );
               vec<char> indx;
               for ( int l = 0; l < ind.isize( ); l++ )
                    indx.push_back( ind[l] );
               UniqueSort(indx);
               if ( indx.size( ) != 1 ) continue;
               const int minh = 10;
               const int hflank = 5;
               int pos = v.pos;
               int left = 0, right = 0;
               for ( int j = pos; j >= 0; j-- )
               {    if ( s[j] != indx[0] ) break;
                    left++;    }
               for ( int j = pos+1; j < s.isize( ); j++ )
               {    if ( s[j] != indx[0] ) break;
                    right++;    }
               Bool reflank = False;
               if ( v.alt.size( ) > 1 )
               {    if ( ( left > 0 || right > 0 ) &&
                         left + v.alt.isize( ) - 1 + right >= minh )
                    {    reflank = True;    }    }
               if ( v.ref.size( ) > 1 )
               {    if ( left + right >= minh ) reflank = True;    }
               if ( !reflank ) continue;
               const int uflank = 20;
               String leftq = s.substr( v.pos + 1 - uflank, uflank );
               int rightp = v.pos + 1;
               while( rightp < s.isize( ) - 1 && s[rightp] == indx[0] ) rightp++;
               String rightq = s.substr( rightp, uflank );

               int base = left + right;
               int delta = ( v.ref.size( ) > 1 ? -( v.ref.isize( ) - 1 )
                    : v.alt.isize( ) - 1 );

               int bothc = 0;
               vec<int> allc = { -1, 0, 1, delta-1, delta, delta+1 };
               UniqueSort(allc);
               for ( int l = 0; l < allc.isize( ); l++ )
               {    String x = leftq + String( base + allc[l], ind[0] ) + rightq;
                    if ( pass == 1 ) 
                    {    to_search.push_back(x);
                         continue;    }
                    bothc += Search( search_db, reads, reads_id, home, x, "g" );    }

               if ( pass == 1 ) 
               {    to_search.push_back( leftq, rightq );
                    continue;    }

               int leftc = Search( search_db, reads, reads_id, home, leftq, "g" );
               int rightc = Search( search_db, reads, reads_id, home, rightq, "g" );
               const int eflank = 5;
               vec<int> check;
               check.push_back( delta - 1, delta, delta + 1 );
               UniqueSort(check);
               int ncheck = check.size( );
               vec<int> count1(ncheck), count2(ncheck), count3(ncheck);
               vec<int> count(ncheck);

               for ( int j = 0; j < check.isize( ); j++ )
               {    String query;
                    int d = check[j];
                    if ( d < 0 )
                    {    query = s.substr( v.pos + 1 - eflank, eflank )
                              + s.substr( v.pos + 1 - d,
                              eflank + rightp - v.pos - 1 + d );    }
                    if ( d == 0 )
                    {    query = s.substr( 
                              v.pos + 1 - eflank, 
                              2 * eflank + rightp - v.pos - 1 );    }
                    if ( d > 0 )
                    {    query = s.substr( v.pos + 1 - eflank, eflank )
                              + String( d, indx[0] ) + s.substr( 
                              v.pos + 1, eflank + rightp - v.pos - 1 );    }
                    count1[j] =
                         Search( search_db, reads, reads_id, home, query, "1" );
                    count2[j] = 
                         Search( search_db, reads, reads_id, home, query, "2" );
                    count3[j] = 
                         Search( search_db, reads, reads_id, home, query, "3" );
                    count[j] = count1[j] + count2[j] + count3[j];    }
               if (HOMOPOLYMER_VERBOSE)
               {    cout << "\n";
                    PRINT4( i, v.pos, base, delta );
                    cout << "flanks: " << leftq << "[" << leftc << "], "
                         << rightq << "[" << rightc << "]"
                         << ", bothc = " << bothc << endl;
                    for ( int j = 0; j < check.isize( ); j++ )
                    {    int d = check[j];
                         PRINT4( d, count1[j], count2[j], count3[j] );    }    }
               if ( bothc != 1 ) continue;
               const int cmult = 5;
               const int min_count = 5;
               if ( count[1] < min_count ) continue;
               int alt = 0;
               if ( i > 0 && hp[i-1][ as_char( indx[0] ) ] > 0 
                    && varsx[i-1].pos == varsx[i].pos )
               {    alt = hp[i-1][ as_char( indx[0] ) ];    }
               if ( i < varsx.isize( ) - 1 && hp[i+1][ as_char( indx[0] ) ] > 0
                    && varsx[i+1].pos == varsx[i].pos )
               {    alt = hp[i+1][ as_char( indx[0] ) ];    }
               if ( !( count[1] >= cmult * count[0] || check[0] == alt ) ) continue;
               if ( !( count[1] >= cmult * count[2] || check[2] == alt ) ) continue;
               if (HOMOPOLYMER_VERBOSE) cout << "confirmed" << endl;
               varsx[i].confirmed = True;    }    }

     // Filter and print variants.

     int events = 0;
     Ofstream( aout, home + "/variants.all" );
     Ofstream( aout2, home + "/variants.all.20" );
     int np = 0;
     for ( int i = 0; i < varsx.isize( ); i++ )
     {    const variant& v = varsx[i];
          if ( v.pos < start || v.pos > stop ) continue;

          ostringstream vout, vout2;
          vout << "#" << ID << "  ";
          vout2 << "#" << ID << "  ";
          const int flank = 10;
          const int flank2 = 20;

          for ( int j = 0; j < flank; j++ )
          {    int pos = v.pos - start - flank + j;
               if ( pos >= 0 ) vout << s[pos];    }
          vout << "/";
          for ( int j = 0; j < flank; j++ )
          {    int pos = v.pos - start + v.ref.isize( ) + j;
               if ( pos < s.isize( ) ) vout << s[pos];    }

          for ( int j = 0; j < flank2; j++ )
          {    int pos = v.pos - start - flank2 + j;
               if ( pos >= 0 ) vout2 << s[pos];    }
          vout2 << "/";
          for ( int j = 0; j < flank2; j++ )
          {    int pos = v.pos - start + v.ref.isize( ) + j;
               if ( pos < s.isize( ) ) vout2 << s[pos];    }

          String notes_plus;
          if ( v.GATK_100.state != "uncalled" )
          {    if ( notes_plus != "" ) notes_plus += ",";
               notes_plus += "GATK-100=" + v.GATK_100.state;    }
          if ( v.GATK_250.state != "uncalled" )
          {    if ( notes_plus != "" ) notes_plus += ",";
               notes_plus += "GATK-250=" + v.GATK_250.state;    }
          if ( v.CORTEX.state != "uncalled" )
          {    if ( notes_plus != "" ) notes_plus += ",";
               notes_plus += "CORTEX=" + v.CORTEX.state;    }
          if ( v.DISCOVAR.state != "uncalled" )
          {    if ( notes_plus != "" ) notes_plus += ",";
               notes_plus += "DISCOVAR=" + v.DISCOVAR.state;    }
          if (v.fosmid)
          {    if ( notes_plus != "" ) notes_plus = "Fosmid," + notes_plus;
               else notes_plus = "Fosmid";    }
          if ( v.het_type != "" )
          {    if ( notes_plus != "" ) notes_plus = v.het_type + "," + notes_plus;
               else notes_plus = v.het_type;    }
          if (v.confirmed)
          {    if ( notes_plus != "" ) notes_plus += ",";
               notes_plus += "confirmed";    }
          notes_plus += "," + v.vartype;
          vout << "  \t" << v.pos - start << "\t" << gid << ":" << v.pos << "\t" 
               << v.ref << "\t" << v.alt << "\t" << notes_plus
               << " " << v.friends << endl;
          vout2 << "  \t" << v.pos - start << "\t" << gid << ":" << v.pos << "\t" 
               << v.ref << "\t" << v.alt << "\t" << notes_plus
               << " " << v.friends << endl;

          aout << vout.str( );
          aout2 << vout2.str( );

          // Don't show cases that DISCOVAR calls and which are on the Fosmid
          // or confirmed.

          if ( v.DISCOVAR.state == "het" || v.DISCOVAR.state == "hom" )
          {    if ( v.fosmid || v.confirmed ) continue;    }

          // Don't show cases where DISCOVAR and GATK are in agreement.  The hom
          // cases are generally cases where 'hom' is misinterpreted because TWO 
          // variants are present at the site.  A more correct version would 
          // specifically check for this.

          if ( v.DISCOVAR.state == "het" )
          {    if ( v.GATK_100.state == "het" ) continue;
               if ( v.GATK_250.state == "het" ) continue;    }
          if ( v.DISCOVAR.state == "hom" )
          {    if ( v.GATK_100.state == "hom" ) continue;
               if ( v.GATK_250.state == "hom" ) continue;    }

          // Don't show events if both DISCOVAR and GATK the variant confidently, 
          // even if the genotypes don't agree.  Ditto for Fosmid.

          if ( v.fosmid || v.GATK_100.state == "het" || v.GATK_100.state == "hom"
               || v.GATK_250.state == "het" || v.GATK_250.state == "hom" )
          {    if ( v.DISCOVAR.state == "het" || v.DISCOVAR.state == "hom"
                    || v.DISCOVAR.state.Contains( "/1", -1 ) )
               {    continue;    }    }

          // Don't show isolated homozygous events that are not supported by the
          // Fosmid and which we don't call.

          Bool het = v.GATK_100.state == "het" || v.GATK_250.state == "hom"
               || v.CORTEX.state == "het";
          Bool hom = v.GATK_100.state == "hom" || v.GATK_250.state == "hom"
               || v.CORTEX.state == "hom";
          if ( v.DISCOVAR.state == "uncalled" && !het && hom )
          {    if ( i == 0 || varsx[i-1].Stop( ) < varsx[i].pos )
               {    if ( i == varsx.isize( ) - 1
                         || varsx[i+1].pos > varsx[i].Stop( ) )
                    {    continue;    }    }    }

          // Don't show unconfirmed non-Fosmid events that are not called by
          // DISCOVAR.

          if ( v.DISCOVAR.state == "uncalled" && !v.fosmid && !v.confirmed )
          {    continue;    }

          // Print.

          if (COUNT_ONLY) events++;
          else cout << "[" << ++np << "] " << vout.str( );    }

     if (COUNT_ONLY) 
          cout << "#" << ID << " --> " << events << " events" << endl;

     // Generate variants.all.xxx.

     aout.close( );
     aout2.close( );
     Ofstream( out, home + "/variants.all.xxx" );
     fast_ifstream in1( home + "/variants.all" );
     vec<String> lines1, lines2;
     while(1)
     {    getline( in1, line );
          if ( in1.fail( ) ) break;
          lines1.push_back(line);    }
     fast_ifstream in2( home + "/hom" );
     while(1)
     {    getline( in2, line );
          if ( in2.fail( ) ) break;
          lines2.push_back(line);    }
     for ( int i = 0; i < lines1.isize( ); i++ )
     {    istringstream iline( lines1[i] );
          String id, pos1, ref1, alt1, comments;
          iline >> id >> junk >> junk >> pos1 >> ref1 >> alt1 >> comments;
          for ( int j = 0; j < lines2.isize( ); j++ )
          {    if ( lines2[j].Contains( "count" ) ) continue;
               istringstream jline( lines2[j] );
               String caller, pos2, ref2, alt2;
               jline >> junk >> junk >> junk >> junk >> junk >> caller >> junk
                    >> junk >> pos2 >> junk >> junk >> ref2 >> junk >> junk >> alt2;
               caller = caller.Before( "," );
               pos2 = pos2.Before( "," );
               ref2 = ref2.Before( "," );
               alt2 = alt2.Before( "," );
               if ( pos2 != pos1.After( ":" ) ) continue;
               if ( ref2 != ref1 ) continue;
               if ( alt2 != alt1 ) continue;
               vec<String> com;
               Tokenize( comments, {','}, com );
               for ( int k = 0; k < com.isize( ); k++ )
               {    if ( caller == "DISCOVAR" && com[k] == "DISCOVAR=hom" )
                         com[k] = "DISCOVAR=hom+";
                    if ( caller == "GATK-100" && com[k] == "GATK-100=hom" )
                         com[k] = "GATK-100=hom+";
                    if ( caller == "GATK-250" && com[k] == "GATK-250=hom" )
                         com[k] = "GATK-250=hom+";
                    if ( caller == "CORTEX" && com[k] == "CORTEX=hom" )
                         com[k] = "CORTEX=hom+";    }
               ostringstream bout;
               bout << printSeq(com);
               comments = bout.str( );    }
          vec<String> com;
          Tokenize( comments, {','}, com );
          com.EraseValue( "confirmed" );
          com.EraseValue( "maternal" );
          com.EraseValue( "paternal" );
          for ( int k = 0; k < com.isize( ); k++ )
          {    if ( com[k].Contains( "DISCOVAR=", 0 )
                    && com[k] != "DISCOVAR=hom" && com[k] != "DISCOVAR=het"
                    && com[k] != "DISCOVAR=hom+" )
               {    com.EraseValue( com[k] );
                    break;    }    }
          com.resize( com.size( ) - 1 );
          if ( com.empty( ) ) continue;
          out << id << "\t" << pos1 << "\t" << ref1 << "\t" << alt1 << "\t"
               << printSeq(com) << endl;    }    }
