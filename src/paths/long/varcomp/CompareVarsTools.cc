///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

// MakeDepend: dependency LongProto
// MakeDepend: dependency MakeLookupTable
// MakeDepend: dependency QueryLookupTable

#include <sys/wait.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "ParseSet.h"
#include "PrintAlignment.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "lookup/PerfectLookup.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/long/varcomp/CompareVarsTools.h"
#include "paths/simulation/GenomeVCFToHyperBaseVector.h"
#include "paths/simulation/VCF.h"

// We consider partitions of the 11-member set {d1,...,d5,s1,...,s6}
// (daughters and sons of NA12878) into three subsets.  The first set s1
// corresponds to het sites in the child.  The second set s2 corresponds to 
// homozygous sites (of one type) in the child, and the third set s3 corresponds to
// the other type of homozygous sites.  We assign an integer "genotype id" to
// this partition.

int GenotypeId( const vec<int>& s1, const vec<int>& s2, const vec<int>& s3 )
{    int gt1 = 0, gt2 = 0;
     for ( int j = 0; j < 11; j++ )
     {    if ( j > 0 ) gt1 *= 3;
          if ( BinMember( s1, j ) ) gt1 += 1;
          if ( BinMember( s2, j ) ) gt1 += 2;
          if ( j > 0 ) gt2 *= 3;
          if ( BinMember( s1, j ) ) gt2 += 1;
          if ( BinMember( s3, j ) ) gt2 += 2;    }
     return Min( gt1, gt2 );    }

// Define genotype groups = partititions as above, and compatibility relations.
// We say that (s1,s2,s3) is compatible with (t1,t2,t3) if t1 = s2 union s3
// and s1 = t2 union t3.

void GenotypeGroups( vec<vec<vec<int>>>& groups, vec<vec<int>>& compat )
{    groups.clear( );
     vec<int> gts;
     for ( int i = 0; i < 2048; i++ )
     {    vec<int> s1;
          int x = i;
          for ( int j = 0; j < 11; j++ )
          {    if ( x % 2 == 1 ) s1.push_back(j);
               x /= 2;    }
          for ( int ix = 0; ix < 2048; ix++ )
          {    vec<int> s2, s3;
               int x = ix;
               for ( int j = 0; j < 11; j++ )
               {    if ( x % 2 == 1 ) s2.push_back(j);
                    x /= 2;    }
               if ( Meet( s1, s2 ) ) continue;
               for ( int j = 0; j < 11; j++ )
               {    if ( !BinMember( s1, j ) && !BinMember( s2, j ) )
                         s3.push_back(j);    }
               if ( s3 > s2 ) continue;
               vec< vec<int> > g(3);
               g[0] = s1, g[1] = s2, g[2] = s3;
               int gt = GenotypeId( s1, s2, s3 );
               gts.push_back(gt), groups.push_back(g);    }    }
     vec<vec<int>> groups1( groups.size( ) );
     for ( int i = 0; i < groups.isize( ); i++ )
          groups1[i] = groups[i][0];
     SortSync( groups1, groups, gts );
     compat.clear_and_resize( groups.size( ) );
     #pragma omp parallel for
     for ( int i1 = 0; i1 < groups.isize( ); i1++ )
     {    vec<int> s23 = groups[i1][1];
          s23.append( groups[i1][2] );
          Sort(s23);
          int low = LowerBound( groups1, s23 ), high = UpperBound( groups1, s23 );
          for ( int i2 = low; i2 < high; i2++ )
          {    vec<int> t23 = groups[i2][1];
               t23.append( groups[i2][2] );
               UniqueSort(t23);
               if ( groups[i1][0] == t23 ) 
                    compat[ gts[i1] ].push_back( gts[i2] );    }
          Sort( compat[ gts[i1] ] );    }
     SortSync( gts, groups );    }

void BuildGenotypeMap( const String& root )
{
     // Define groups.

     vec< vec< vec<int> > > groups;
     vec< vec<int> > compat;
     GenotypeGroups( groups, compat );

     // Define samples.
     
     vec<String> samples, aliases;
     samples.push_back( "NA12877", "NA12878", "NA12879", "NA12880",
          "NA12881", "NA12882", "NA12883", "NA12884", "NA12885", "NA12886",
          "NA12887", "NA12888" );
     samples.push_back( "NA12889", "NA12890", "NA12891", "NA12892", "NA12893" );
     aliases.push_back( "dad", "mom", "d1", "d2", "d3", "s1", "s2", "s3", "d4",
          "s4", "d5", "s5" );
     aliases.push_back( "dadsdad", "dadsmom", "momsdad", "momsmom", "s6" );

     // Get raw events.

     vec<event> events;
     #pragma omp parallel for
     for ( int i = 0; i < samples.isize( ); i++ )
     {    const String& sample = samples[i];
          Ifstream( in, "/wga/scr4/human_data/CEPH/" + sample + "_S1.genome.vcf" );
          String line;
          vec<String> fields;
          vec<event> eventsi;
          int count = 0; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( !line.Contains( "chr" , 0 ) || line.Contains( "chrM" , 0 ) ) 
                    continue;
               Tokenize( line, {'\t'}, fields );
               if ( fields.size( ) < 10 ) continue;
               if ( fields[4] == "." || fields[4].Contains( "," ) ) continue;
               if ( fields[3].size( ) > 1 ) continue;
               if ( fields[4].size( ) > 1 ) continue;
               if ( !fields[9].Contains( ":" ) ) continue;
               if ( i == 0 && count++ % 50000 == 0 ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    PRINT2( fields[0], fields[1] ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               eventsi.push( aliases[i], fields[0], fields[1].Int( ), fields[3], 
                    fields[4], fields[9].Before( ":" ), fields[6] );    }
          #pragma omp critical
          {    events.append(eventsi);    }    }

     // Put events in chromosome buckets.  This could be speeded up.

     cout << Date( ) << ": putting events in buckets" << endl; // XXXXXXXXXXXXXXXXXX
     vec<vec<event>> eventsc(24);
     for ( int i = 0; i < events.isize( ); i++ )
     {    String gid = events[i].chr.After( "chr" );
          int g;
          if ( gid == "X" ) g = 22;
          else if ( gid == "Y" ) g = 23;
          else g = gid.Int( ) - 1;
          eventsc[g].push_back( events[i] );    }

     // Loop through chromosomes.

     cout << Date( ) << ": start main for loop" << endl; // XXXXXXXXXXXXXXXXXXXXXXXX
     vec< vec< triple< String, int, vec<int> > > > B(24);
     #pragma omp parallel for
     for ( int c = 0; c < 24; c++ )
     {    
          // Sort.

          Sort( eventsc[c] );
          
          // Group events.

          vec< triple<String,int,int> > G;
          for ( int i = 0; i < eventsc[c].isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < events.isize( ); j++ )
               {    if ( eventsc[c][j].chr != eventsc[c][i].chr ) break;
                    if ( eventsc[c][j].pos != eventsc[c][i].pos ) break;    }
               Bool ok = True;
               for ( int k = i; k < j; k++ )
               {    if ( eventsc[c][k].ref != eventsc[c][i].ref ) ok = False;
                    if ( eventsc[c][k].alt != eventsc[c][i].alt ) ok = False;
                    if ( eventsc[c][k].filter != "PASS" ) ok = False;    }
               if ( !ok )
               {    i = j - 1;
                    continue;    }
               Bool mom_het = False, dad_het = False;
               for ( int k = i; k < j; k++ )
               {    const event& e = eventsc[c][k];
                    if ( e.alias == "mom" && e.genotype == "0/1" ) mom_het = True;
                    if ( e.alias == "dad" && e.genotype == "0/1" ) 
                         dad_het = True;    }
               if ( mom_het && dad_het )
               {    String momsmom = "0/0", momsdad = "0/0"; 
                    String dadsmom = "0/0", dadsdad = "0/0";
                    String mom = "0/1", dad = "0/1";
                    vec<String> d(5, "0/0"), s(6, "0/0");
                    for ( int k = i; k < j; k++ )
                    {    const event& e = eventsc[c][k];
                         if ( e.alias == "momsmom" ) momsmom = e.genotype;
                         if ( e.alias == "momsdad" ) momsdad = e.genotype;
                         if ( e.alias == "dadsmom" ) dadsmom = e.genotype;
                         if ( e.alias == "dadsdad" ) dadsdad = e.genotype;
                         if ( e.alias == "d1" ) d[0] = e.genotype;
                         if ( e.alias == "d2" ) d[1] = e.genotype;
                         if ( e.alias == "d3" ) d[2] = e.genotype;
                         if ( e.alias == "d4" ) d[3] = e.genotype;
                         if ( e.alias == "d5" ) d[4] = e.genotype;
                         if ( e.alias == "s1" ) s[0] = e.genotype;
                         if ( e.alias == "s2" ) s[1] = e.genotype;
                         if ( e.alias == "s3" ) s[2] = e.genotype;
                         if ( e.alias == "s4" ) s[3] = e.genotype;
                         if ( e.alias == "s5" ) s[4] = e.genotype;
                         if ( e.alias == "s6" ) s[5] = e.genotype;    }
                    int gt1 = 0, gt2 = 0;
                    for ( int l = 0; l < 5; l++ )
                    {    if ( l > 0 ) gt1 *= 3;
                         if ( d[l] == "0/1" ) gt1 += 1;
                         if ( d[l] == "0/0" ) gt1 += 2;    }
                    for ( int l = 0; l < 6; l++ )
                    {    gt1 *= 3;
                         if ( s[l] == "0/1" ) gt1 += 1;
                         if ( s[l] == "0/0" ) gt1 += 2;    }
                    for ( int l = 0; l < 5; l++ )
                    {    if ( l > 0 ) gt2 *= 3;
                         if ( d[l] == "0/1" ) gt2 += 1;
                         if ( d[l] == "1/1" ) gt2 += 2;    }
                    for ( int l = 0; l < 6; l++ )
                    {    gt2 *= 3;
                         if ( s[l] == "0/1" ) gt2 += 1;
                         if ( s[l] == "1/1" ) gt2 += 2;    }
                    int gt = Min( gt1, gt2 );
                    Bool verbose = False;
                    if (verbose)
                    {    cout << "\n" << eventsc[c][i].chr << "\t" 
                              << eventsc[c][i].pos << "\t" << eventsc[c][i].ref 
                              << "\t" << eventsc[c][i].alt << endl;
                         PRINT3( momsmom, momsdad, mom );
                         PRINT3( dadsmom, dadsdad, dad );
                         PRINT5( d[0], d[1], d[2], d[3], d[4] );
                         PRINT6( s[0], s[1], s[2], s[3], s[4], s[5] );    
                         PRINT(gt);    }
                    G.push( eventsc[c][i].chr, eventsc[c][i].pos, gt );    }
               i = j - 1;    }

          // Build blocks.

          for ( int i = 0; i < G.isize( ); i++ )
          {    vec<int> friends;
               for ( int j = i - 1; j >= 0; j-- )
               {    if ( G[j].third == G[i].third ) continue;
                    if ( !BinMember( compat[ G[i].third ], G[j].third ) ) break;
                    friends.push_back( G[j].third );    }
               for ( int j = i + 1; j < G.isize( ); j++ )
               {    if ( G[j].third == G[i].third ) continue;
                    if ( !BinMember( compat[ G[i].third ], G[j].third ) ) break;
                    friends.push_back( G[j].third );    }
               UniqueSort(friends);
               for ( int l = 0; l < friends.isize( ); l++ )
               {    vec<int> p;
                    p.push_back( G[i].third );
                    p.push_back( friends[l] );
                    Sort(p);
                    B[c].push( G[i].first, G[i].second, p );    }    }    
          #pragma omp critical
          {    cout << "done with " << c << endl;    }    }

     // Print.

     cout << Date( ) << ": printing" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     Ofstream( out, root + "/genotype.map" );
     vec<vec<vec<String>>> rows(24);
     #pragma omp parallel for
     for ( int c = 0; c < 24; c++ )
     {    for ( int i = 0; i < B[c].isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < B[c].isize( ); j++ )
                    if ( B[c][j].third != B[c][i].third ) break;
               int span = B[c][j-1].second - B[c][i].second + 1;
               const int min_span = 1000;
               if ( span >= min_span )
               {    vec<String> row;
                    row.push_back( B[c][i].first, ToString( B[c][i].second ),
                         ToString( B[c][j-1].second ), ToString(span), 
                         ToString( j - i ), "{" + ToString( B[c][i].third[0] ) + ","
                         + ToString( B[c][i].third[1] ) + "}" );
                    rows[c].push_back(row);    }
               i = j - 1;    }    }
     for ( int c = 0; c < 24; c++ )
          PrintTabular( out, rows[c], 2, "lrrrr" );    }

int Search( disk_map& search_db, vec<vecbasevector>& reads,
     const vec<String>& reads_id, const String& home, const String& query, 
     const String& dataset )
{    String key = ToString(dataset) + "_" + query;
     if ( search_db.Defined(key) ) return search_db.Value(key);
     int pos = Position( reads_id, dataset );
     if ( reads[pos].size( ) == 0 )
     {    String fastb;
          if ( pos == 0 ) fastb = "/wga/scr4/bigrefs/human19/genome.fastb";
          else fastb = home + "/tmp.xxx" + dataset + "/frag_reads_orig.fastb";
          reads[pos].ReadAll(fastb);    }
     int count = 0;
     if ( pos == 0 )
     {    vec<look_align> aligns;
          vecbasevector b;
          b.push_back( basevector(query) );
          PerfectLookup( 12, b,
               "/wga/scr4/bigrefs/human19/genome.lookup.lookup", aligns, FW_OR_RC );
          UniqueSort(aligns);
          count = aligns.size( );    }
     else
     {
          #pragma omp parallel for
          for ( int64_t id = 0; id < (int64_t) reads[pos].size( ); id++ )
          {    String b = reads[pos][id].ToString( );
               if ( b.Contains(query) ) count++;
               else
               {    StringReverseComplement( b, b );
                    if ( b.Contains(query) ) count++;    }    }    }
     search_db.Set( key, count );
     return count;    }

void RunFunctions( const String& FUN, const String& R, const String& home,
     const String& ref_dir, const String& loc, const String& ID, 
     const basevector& query, const basevector& target, const String& fos_fasta, 
     const int rstart, const int rstop )
{
     // Determine functions.

     vec<String> functions;
     vec<String> all = { "ASSEMBLE", "GATK", "READS", "ALIGN" };
     ParseStringSet( FUN, functions );
     if ( Member( functions, String("ALL") ) ) functions = all;
     else
     {    for ( int i = 0; i < functions.isize( ); i++ )
          {    if ( !Member( all, functions[i] ) && functions[i] != "CALL" )
               {    cout << "\nIllegal function: " << functions[i] << endl;
                    Scram(1);    }    }    }

     // Assemble and call variants.

     String longproto = "LongProto";
     if ( R != "" )
     {    ForceAssertEq( R.isize( ), 5 );
          longproto = "/wga/builds/" + R.substr( 0, 3 ) + "/" + R.substr( 3, 2 )
               + "/bin/LongProto";    }
     if ( Member( functions, String("ASSEMBLE") ) )
     {    SystemSucceed( longproto + " SAMPLE=human READS=#picard TMP=" + home 
               + "/tmp.xxx X=" + loc + " OUT_INT_HEAD=" + home + "/aaa "
               + "LOGGING=REFTRACE_VARIANTS=True > " + home + "/LongProto.out" );
          System( "/bin/rm " + home + "/aaa.[0-9].*" );
          System( "/bin/rm " + home + "/aaa.[0-9][0-9].*" );
         }
     if ( Member( functions, String("CALL") ) 
          && !Member( functions, String("ASSEMBLE") ) )
     {    SystemSucceed( longproto + " SAMPLE=human READS=#picard TMP=" + home 
               + "/tmp.xxx X=" + loc + " OUT_INT_HEAD=" + home + "/aaa "
               + "IN_SHBV_FINAL=" + home + "/aaa.final.shbv "
               + "LOGGING=REFTRACE_VARIANTS=True > " + home 
               + "/LongProto.out.vc" );    }

     // Get GATK and CORTEX vcf.raw files.

     if ( Member( functions, String("GATK") ) )
     {    String vcf_dir = "/wga/scr4/NA12878_calls";
          String vcf;
          for ( int pass = 1; pass <= 3; pass++ )
          {    String VCF_CHR_PREFIX;
               int GATK_RL = ( pass == 1 ? 100 : 250 );
               if ( pass <= 2 )
               {    if ( GATK_RL == 250 )
                    {    
                         vcf = "mem4/haplotype-caller-pcr-none/"
                              "reverted.12.aligned.wholegenome.sorted."
                              "indel_cleaned_local.recal.unfiltered."
                              "recal_snp_recal_indel.vcf";

                         // vcf = "mem4/haplotype-caller/reverted.12.aligned."
                         //      "wholegenome.sorted.indel_cleaned_local.recal."
                         //      "unfiltered.recal_snp_recal_indel.vcf";

                         // vcf = "../../../humgen/gsa-hpprojects/dev/depristo/"
                         //      "oneOffProjects/ceuTrioBestPractices/V37.gatk27/"
                         //      "ceuTrio.bwamem.HaplotypeCaller.recalibrated.vcf";
                         // vcf = "heng/gatk-250hc-flt600.vcf";
                         // vcf = "mem4/reverted.12.aligned.wholegenome.sorted."
                         //      "indel_cleaned_local.recal."
                         //      "unfiltered.recal_snp_recal_indel.phased.vcf";    
                              }
                    else 
                    {    vcf = "platinum/NA12878_S1.genome.vcf";
                         VCF_CHR_PREFIX = "chr";    }    }
               if ( pass == 3 )
               {    
                    // BC calls from CORTEX paper.

                    // vcf = "cortex/hg18_to_hg19/vcf/"
                    //      "NA12878_BC_calls.decomp.submission.hg19.vcf";

                    // PD calls from CORTEX paper.

                    // vcf = "cortex/hg18_to_hg19/vcf/"
                    //      "NA12878_PD_calls.decomp.submission.hg19.vcf";

                    // Final calls.

                    //vcf = "cortex/cortex-NA12878-again/vcf/NA12878.bc3161pd31.vcf";
                    vcf = "cortex/cortex-NA12878-again/vcf/"
                         "NA12878.decomp.processed.vcf";

                    // Unclipped corrected k=31.

                    // vcf = "cortex/cortex-NA12878-uncorrected/FINAL_RESULTS/vcfs/"
                    //      "default_vcfname_wk_flow_I_RefCC_FINALcombined_BC_"
                    //      "calls_at_all_k.decomp.vcf";

                    // Unclipped uncorrected.

                    // vcf = "cortex/cortex-NA12878-uncorrected26/FINAL_RESULTS/"
                    //      "vcfs/default_vcfname_wk_flow_I_RefCC_FINALcombined_BC_"
                    //      "calls_at_all_k.decomp.vcf";

                    // Unfiltered k=31.

                    // vcf = "cortex/cortex-NA12878-tmp/vcf-unfiltered/sorted31.vcf";

                    // Filtered k=31.

                    // vcf = "cortex/cortex-NA12878-tmp/vcf/NA12878.K31.vcf";

                    // Filtered k=31 PD.

                    // vcf = "cortex/cortex-NA12878/FINAL_RESULTS/vcfs/"
                    //   "default_vcfname_wk_flow_I_RefCC_FINALcombined_PD_"
                    //   "calls_at_all_k.decomp.vcf"; 

                         }
               String prefix;
               String VCF_FILE = vcf_dir + "/" + vcf;
               vecbasevector genome( ref_dir + "/genome.fastb" );
               vecString chr_names;
               chr_names.ReadAll( ref_dir + "/genome.fastb.names" );
               for_each( chr_names.begin(), chr_names.end(), [] (String& name )
               {    size_t pos = name.find_first_of(' ');
                    if ( pos != string::npos ) name.resize(pos);    }
               );
               Range coords( ID, chr_names );
               String tail = ( pass <= 2 ? ToString(GATK_RL) : "CORTEX" );
               VCF variants( VCF_FILE, VCF_CHR_PREFIX+coords.chr( ), 
                    coords.start1( ), coords.end1( ), home + "/vcf.raw." 
                    + tail +".unsimplified");
               variants.writeSimplifiedVCF(home+"/vcf.raw."+tail);
          }    }

     // Get reads.

     int id = ID.Int( );
     if ( Member( functions, String("READS") ) )
     {    
          // Get trio reads.

          SystemSucceed( longproto + " SAMPLE=human READS=#picard TMP=" + home
               + "/tmp.xxx1 X=" + loc + " EXIT=LOAD DATASET=1,2 > " 
               + home + "/logs/xxx1" );
          SystemSucceed( longproto + " SAMPLE=human READS=#picard TMP=" + home
               + "/tmp.xxx2 X=" + loc + " EXIT=LOAD DATASET=4 > " 
               + home + "/logs/xxx2" );
          SystemSucceed( longproto + " SAMPLE=human READS=#picard TMP=" + home
               + "/tmp.xxx3 X=" + loc + " EXIT=LOAD DATASET=5 > " 
               + home + "/logs/xxx3" );

          // Get Fosmid reads.

          String pool = ( id <= 55 ? "hpool2" : "hpool3" );
          SystemSucceed( longproto + " SAMPLE=" + pool + " READS=#picard X=" + ID
               + " TMP=" + home + "/tmp.fos DATA_SPEC=SELECT_FRAC=1 "
               + "EXIT=LOAD > " + home + "/logs/xxx_fos 2>&1" );

          // Get reads from Platinum datasets.  

          for ( auto p : { "dad", "dadsmom", "dadsdad", "d1", "d2", "d3", "d4", 
               "d5", "s1", "s2", "s3", "s4", "s5", "s6" } )
          {    SystemSucceed( longproto + " SAMPLE=human READS=#picard TMP=" 
                    + home + "/tmp.xxx" + p + " X=chr" + loc 
                    + " EXIT=LOAD DATASET=P_" + p 
                    + " LOGGING=TREAT_AS_UNKNOWN=True > " 
                    + home + "/logs/xxx_" + p + " 2>&1" );    }    }

     // Get alignment.

     if ( Member( functions, String("ALIGN") ) )
     {    basevector region( target, rstart, rstop - rstart );
          {    Ofstream( rout, home + "/region.fasta" );
               region.Print( rout, "region_" + ID );    }
          SystemSucceed( "MakeLookupTable SOURCE=" + home + "/region.fasta "
               "OUT_HEAD=" + home + "/region LO=True > " 
               + home + "/logs/MakeLookupTable.out" );
          SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS=" 
               + fos_fasta + " L=" + home 
               + "/region.lookup PARSEABLE=True SMITH_WAT=True " 
               + " VISUAL=True > " + home + "/region.aligns" );
          vec<look_align> raligns;
          LoadLookAligns( home + "/region.aligns", raligns );
          ForceAssertEq( raligns.isize( ), 1 );
          const look_align& r = raligns[0];
          int start = rstart + r.pos2( ), stop = rstart + r.Pos2( );
          if ( id == 29 ) { start = 7028567, stop = 7066173; }
          if ( id == 34 ) { start = 104039548, stop = 104081406; }
          if ( id == 37 ) { start = 30434721, stop = 30473940; }
          if ( id == 62 ) { start = 10103953, stop = 10143530; }
          if ( id == 99 ) { start = 2110088, stop = 2149370; }
          basevector region2 = basevector( target, start, stop - start );
          alignment a;
          SmithWatAffineParallel( query, region2, a, true, true, 3, 10, 1 );
          Ofstream( out, home + "/region.visual_align" );
          PrintVisualAlignment( True, out, query, region2, a );
          align q(a);
          q.AddToPos2(start);
          BinaryWriter::writeFile( home + "/region.align", q );    }    }

void DefineGenotypeMap( const String& root, 
     vec< pair< triple<String,int,int>, vec<int> > >& gmap,
     vec< vec< vec<int> > >& groups, vec< vec<int> >& compat )
{
     if ( !IsRegularFile( root + "/genotype.map" ) ) 
     {    GenotypeGroups( groups, compat );
          BinaryWriter::writeFile( root + "/genotype.groups", groups );
          BinaryWriter::writeFile( root + "/genotype.compat", compat );
          BuildGenotypeMap( root );    }
     else
     {    BinaryReader::readFile( root + "/genotype.groups", &groups );
          BinaryReader::readFile( root + "/genotype.compat", &compat );    }
     {    fast_ifstream in( root + "/genotype.map" );
          String line;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               istringstream iline( line.c_str( ) );
               String chr, gid;
               int start, stop, len, count;
               iline >> chr >> start >> stop >> len >> count >> gid;
               vec<int> gidp;
               gidp.push_back( gid.Between( "{", "," ).Int( ) );
               gidp.push_back( gid.Between( ",", "}" ).Int( ) );
               gmap.push_back( make_pair( 
                    make_triple( chr, start, stop ), gidp ) );    }    }    }

void RunSome( const String& root, const String& ID, String& com, 
     const Bool COUNT_ONLY, const Bool parallel )
{    vec<int> ids, all = AllFosmids( );
     if ( ID == "all" ) ids = all;
     else if ( ID.Contains( "/" ) )
     {    int k = ID.Before( "/" ).Int( ), n = ID.After( "/" ).Int( );
          ForceAssertGe( k, 1 );
          ForceAssertLe( k, n );
          int N = all.size( );
          for ( int i = ((k-1)*N)/n; i < (k*N)/n; i++ )
               ids.push_back( all[i] );    }
     else ParseIntSet( ID, ids );
     int total_events = 0;
     Bool fail = False;
     String line;
     if ( !parallel )
     {    for ( int i = 0; i < ids.isize( ); i++ )
          {    vec<String> args;
               Tokenize( com, {' '}, args );
               String id_arg;
               for ( int i = 0; i < args.isize( ); i++ )
                    if ( args[i].Contains( "ID=", 0 ) ) id_arg = args[i];
               com.ReplaceBy( id_arg, "ID=" + ToString(ids[i]) );
               if ( !com.Contains( "NH=True" ) ) com += " NH=True";
               Bool print = False;
               int status = System( com + " > CompareVars.one.out" );
               if ( status != 0 ) 
               {    cout << "\nWarning: command failed." << endl;
                    fail = True;
                    print = True;    }
               if ( FileSize( "CompareVars.one.out" ) > 0 ) print = True;
               if (print)
               {    if ( !COUNT_ONLY )
                    {    cout << "-----------------------------------------------"
                         << "-------------------------------------\n";    }
                    fast_ifstream in( "CompareVars.one.out" );
                    while(1)
                    {    getline( in, line );
                    if ( in.fail( ) ) break;
                         total_events++;    }
                    CpAppend( "CompareVars.one.out", cout );    }    }    }
     else
     {    
          // Batch size optimized for FUN=CALL, should be specific to operation.
          // Optimized assuming a single local minimum.
          // Optimized on crd8, a 48-core machine.
          // Optimized on r47305.

          const int batch = 37;
          for ( int i0 = 0; i0 < ids.isize( ); i0 += batch )
          {    vec<int> pids, run_ids, rel_ids;
               vec<Bool> running;
               for ( int i = i0; i < Min( ids.isize( ), i0 + batch ); i++ )
               {    vec<String> args;
                    Tokenize( com, {' '}, args );
                    String id_arg;
                    for ( int j = 0; j < args.isize( ); j++ )
                         if ( args[j].Contains( "ID=", 0 ) ) id_arg = args[j];
                    com.ReplaceBy( id_arg, "ID=" + ToString(ids[i]) );
                    if ( !com.Contains( "NH=True" ) ) com += " NH=True";
                    int rel_id = i - i0 + 1;
                    int pid = Fork( 
                         com + " > CompareVars.one.out." + ToString(rel_id) );
                    pids.push_back(pid);
                    run_ids.push_back( ids[i] );
                    running.push_back(True);
                    rel_ids.push_back(rel_id);    }
               while( Sum(running) > 0 )
               {    int status = -2;
                    int pid = wait( &status );
                    int p = Position( pids, pid );
                    ForceAssert( p >= 0 );
                    running[p] = False;
                    Bool print = False;
                    if ( status != 0 )
                    {    cout << "\nWarning: command for ID " << run_ids[p] 
                              << " failed." << endl;
                         fail = True;
                         print = True;    }
                    if ( FileSize( 
                         "CompareVars.one.out." + ToString(rel_ids[p]) ) > 0 ) 
                    {    print = True;    }
                    if (print)
                    {    if ( !COUNT_ONLY )
                         {    cout << "--------------------------------------------"
                              << "----------------------------------------\n";    }
                         fast_ifstream in( 
                              "CompareVars.one.out." + ToString(rel_ids[p]) );
                         while(1)
                         {    getline( in, line );
                         if ( in.fail( ) ) break;
                              total_events++;    }
                         CpAppend( "CompareVars.one.out." + ToString(rel_ids[p]), 
                              cout );    
                         flush(cout);    }    }    }    }
     if ( ID == "all" )
     {    cout << "\n============================================="
               << "=======================================\n\n";
          cout << "SUMMARY STATS\n";
          if (fail) cout << "*** THERE WAS A FAILURE ***\n";
          cout << "total suspicious DISCOVAR variants = " << total_events 
               << endl;    
          vec<String> callsets;
          callsets.push_back( "GATK-100", "GATK-250", "CORTEX", "DISCOVAR" );
          int A = all.size( ), C = callsets.size( );
          vec< vec<String> > lines(A);
          for ( int i = 0; i < A; i++ )
          {    fast_ifstream in( root + "/" + ToString(all[i]) + "/variants.all" );
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    lines[i].push_back(line);    }    }

          vec< vec<int> > all_het_Fosmid_confirm( C, vec<int>( A, 0 ) );
          vec< vec<int> > all_het_not_Fosmid_confirm( C, vec<int>( A, 0 ) );
          vec< vec<int> > all_het_Fosmid( C, vec<int>( A, 0 ) );
          vec< vec<int> > all_het_not_Fosmid( C, vec<int>( A, 0 ) );
          vec< vec<int> > all_hom_Fosmid( C, vec<int>( A, 0 ) );
          vec< vec<int> > all_Fosmid( C, vec<int>( A, 0 ) );
          for ( int i = 0; i < A; i++ )
          for ( int c = 0; c < C; c++ )
          {    const String& cs = callsets[c];
               for ( int j = 0; j < lines[i].isize( ); j++ )
               {    String line = lines[i][j];
                    Bool confirm = line.Contains( "confirmed" );
                    Bool F = line.Contains( "Fosmid" );
                    if ( line.Contains(cs + "=het") && F && confirm ) 
                         all_het_Fosmid_confirm[c][i]++;
                    if ( line.Contains(cs + "=het") && !F && confirm )
                         all_het_not_Fosmid_confirm[c][i]++;
                    if (F) all_Fosmid[c][i]++;
                    if ( line.Contains( cs + "=het" ) && F ) 
                         all_het_Fosmid[c][i]++;
                    if ( line.Contains( cs + "=het") && !F )
                         all_het_not_Fosmid[c][i]++;
                    if ( line.Contains( cs + "=hom") && F ) 
                         all_hom_Fosmid[c][i]++;    }    }

          int top = 1000;
          for ( int c = 0; c < C; c++ )
          {    vec<double> het_FPs, FNs;
               for ( int pass = 0; pass <= top; pass++ )
               {    const String& cs = callsets[c];
                    if ( pass == top ) cout << "\n" << cs << endl;
                    int het_Fosmid = 0, het_not_Fosmid = 0, hom_Fosmid = 0;
                    int Fosmid = 0, het_Fosmid_confirm = 0;
                    int het_not_Fosmid_confirm = 0;
                    for ( int ip = 0; ip < A; ip++ )
                    {    int i;
                         if ( pass == top ) i = ip;
                         else i = randomx( ) % all.isize( );
                         het_Fosmid_confirm += all_het_Fosmid_confirm[c][i];
                         het_not_Fosmid_confirm += all_het_not_Fosmid_confirm[c][i];
                         het_Fosmid += all_het_Fosmid[c][i];
                         het_not_Fosmid += all_het_not_Fosmid[c][i];
                         hom_Fosmid += all_hom_Fosmid[c][i];
                         Fosmid += all_Fosmid[c][i];    }

                    // if ( pass == top ) PRINT2( het_Fosmid, het_not_Fosmid );
                    int het = het_Fosmid + het_not_Fosmid;
               
                    double het_Fosmid_confirm_rate =
                         double(het_Fosmid_confirm) / double(het_Fosmid);
                    double het_not_Fosmid_confirm_rate =
                         double(het_not_Fosmid_confirm)/double(het_not_Fosmid);
                    if ( pass == top )
                    {    cout << "Fosmid calls = " << het_Fosmid + hom_Fosmid 
                              << endl;
                         cout << "het: Fosmid confirm rate = "
                              << PERCENT_RATIO( 3, het_Fosmid_confirm, het_Fosmid ) 
                              << ", ";
                         cout << "!Fosmid confirm rate = "
                              << PERCENT_RATIO( 3, het_not_Fosmid_confirm,
                              het_not_Fosmid ) << endl;
                         // cout << "hets = " << het << endl;    
                              }

                    double het_error_rate = ( double(het_not_Fosmid) / double(het) )
                         * ( het_Fosmid_confirm_rate - het_not_Fosmid_confirm_rate )
                         / het_not_Fosmid_confirm_rate;
                    if ( pass == top )
                    {    cout << "inferred het FPs = " << fixed 
                              << setprecision(2) << 100.0 * het_error_rate 
                              << "% +/- " << StdDev( het_FPs, Mean(het_FPs) )
                              << "% " << "(" << setprecision(1) 
                              << het * het_error_rate << ")" << endl;
                         cout << "FNs = " << PERCENT_RATIO( 1,
                              Fosmid - het_Fosmid - hom_Fosmid, Fosmid ) 
                              << " +/- " << StdDev( FNs, Mean(FNs) ) << "%" << " (" 
                              << 2 * ( Fosmid - het_Fosmid - hom_Fosmid )
                              << ")" << endl;    }    
                    else 
                    {    het_FPs.push_back( 100.0 * het_error_rate );
                         FNs.push_back( 100.0 * 
                              double(Fosmid - het_Fosmid - hom_Fosmid )
                              / double(Fosmid) );    }    }    }    

          // Analysis of hom FPs.

          cout << "\nhom FPs\n";
          vec<int> count_DISCOVAR( A, 0 ), bad_DISCOVAR( A, 0 );
          vec<int> count_GATK_100( A, 0 ), bad_GATK_100( A, 0 );
          vec<int> count_GATK_250( A, 0 ), bad_GATK_250( A, 0 );
          vec<int> count_CORTEX( A, 0 ), bad_CORTEX( A, 0 );
          for ( int i = 0; i < A; i++ )
          {    fast_ifstream in( root + "/" + ToString(all[i]) + "/hom_fp" );
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    if ( line.Contains( "DISCOVAR" ) && !line.Contains( "count" ) )
                         bad_DISCOVAR[i]++;
                    if ( line.Contains( "DISCOVAR" ) && line.Contains( "count" ) )
                         count_DISCOVAR[i] = line.After( "count = " ).Int( );
                    if ( line.Contains( "GATK-100" ) && !line.Contains( "count" ) )
                         bad_GATK_100[i]++;
                    if ( line.Contains( "GATK-100" ) && line.Contains( "count" ) )
                         count_GATK_100[i] = line.After( "count = " ).Int( );
                    if ( line.Contains( "GATK-250" ) && !line.Contains( "count" ) )
                         bad_GATK_250[i]++;
                    if ( line.Contains( "GATK-250" ) && line.Contains( "count" ) )
                         count_GATK_250[i] = line.After( "count = " ).Int( );
                    if ( line.Contains( "CORTEX" ) && !line.Contains( "count" ) )
                         bad_CORTEX[i]++;
                    if ( line.Contains( "CORTEX" ) && line.Contains( "count" ) )
                         count_CORTEX[i] = line.After( "count = " ).Int( );    }    }
          vec<double> hom_fps_DISCOVAR, hom_fps_GATK_100, hom_fps_GATK_250;
          vec<double> hom_fps_CORTEX;
          for ( int i = 0; i < top; i++ )
          {    int count_DISCOVAR_i = 0, bad_DISCOVAR_i = 0;
               int count_GATK_100_i = 0, bad_GATK_100_i = 0;
               int count_GATK_250_i = 0, bad_GATK_250_i = 0;
               int count_CORTEX_i = 0, bad_CORTEX_i = 0;
               for ( int j = 0; j < A; j++ )
               {    int k = randomx( ) % A;
                    count_DISCOVAR_i += count_DISCOVAR[k];
                    bad_DISCOVAR_i += bad_DISCOVAR[k];
                    count_GATK_100_i += count_GATK_100[k];
                    bad_GATK_100_i += bad_GATK_100[k];
                    count_GATK_250_i += count_GATK_250[k];
                    bad_GATK_250_i += bad_GATK_250[k];
                    count_CORTEX_i += count_CORTEX[k];
                    bad_CORTEX_i += bad_CORTEX[k];    }
               hom_fps_DISCOVAR.push_back( 100.0 * double(bad_DISCOVAR_i)
                    / double(count_DISCOVAR_i) );
               hom_fps_GATK_100.push_back( 100.0 * double(bad_GATK_100_i)
                    / double(count_GATK_100_i) );
               hom_fps_GATK_250.push_back( 250.0 * double(bad_GATK_250_i)
                    / double(count_GATK_250_i) );
               hom_fps_CORTEX.push_back( 100.0 * double(bad_CORTEX_i)
                    / double(count_CORTEX_i) );    }
          cout << "DISCOVAR: " << PERCENT_RATIO( 2, Sum(bad_DISCOVAR),
               Sum(count_DISCOVAR) ) << " +/- " 
               << StdDev( hom_fps_DISCOVAR, Mean(hom_fps_DISCOVAR) ) << "%" 
               << " (" << Sum(bad_DISCOVAR) << " of " << Sum(count_DISCOVAR) 
               << ")" << endl;
          cout << "GATK-100: " << PERCENT_RATIO( 2, Sum(bad_GATK_100),
               Sum(count_GATK_100) ) << " +/- " 
               << StdDev( hom_fps_GATK_100, Mean(hom_fps_GATK_100) ) << "%"
               << " (" << Sum(bad_GATK_100) << " of " << Sum(count_GATK_100) 
               << ")" << endl;
          cout << "GATK-250: " << PERCENT_RATIO( 2, Sum(bad_GATK_250),
               Sum(count_GATK_250) ) << " +/- " 
               << StdDev( hom_fps_GATK_250, Mean(hom_fps_GATK_250) ) << "%"
               << " (" << Sum(bad_GATK_250) << " of " << Sum(count_GATK_250) 
               << ")" << endl;
          cout << "CORTEX:   " << PERCENT_RATIO( 2, Sum(bad_CORTEX),
               Sum(count_CORTEX) ) << " +/- " 
               << StdDev( hom_fps_CORTEX, Mean(hom_fps_CORTEX) ) << "%"
               << " (" << Sum(bad_CORTEX) << " of " << Sum(count_CORTEX) << ")" 
               << endl;
          cout << endl;
          if (fail) cout << "*** THERE WAS A FAILURE ***\n" << endl;

          if (fail) Scram(1);    }    }
