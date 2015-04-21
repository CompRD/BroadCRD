///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeRef.  Use data from parents of NA12878 to phase a simple assembly of a
// region.  By a simple assembly, we mean that its graph is a line, interrupted only
// by two-fold bubbles.  Note that if we can't phase a bubble, we assign it randomly
// Results might be improved by using haplotype knowledge for the human genome.  Or
// by using alignments rather than kmers, as used here.  Jump data might help too.
//
// Then use the original reads to validate the assembly, as compared to a given
// reference sequence.
//
// The input assembly may be presented either as a HyperBasevector or as efasta.

// MakeDepend: dependency SAM2CRDDump

#include "Basevector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "random/Random.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(A_SHBV, "", ".shbv file for assembly");
     CommandArgument_String_OrDefault_Doc(A_EFASTA, "", ".efasta file for assembly");
     CommandArgument_Int_OrDefault_Doc(K, -1, 
          "K value for assembly, needed for A_EFASTA");
     CommandArgument_String_OrDefault_Doc(A_EFASTA_OUT, "", 
          "output A_EFASTA_OUT.{efasta,dot} file for assembly");
     CommandArgument_String_Doc(X, 
          "region of genome to use, following notation of LongProto");
     CommandArgument_String_Doc(TMP, "temporary directory where "
          "frag_reads_orig.{fastb,qualb} are located");
     CommandArgument_Bool_OrDefault_Doc(SHOW_ALL, False, "show alignments of reads "
          "favoring the assembly, not just the reference");
     CommandArgument_String_OrDefault_Doc(OUT, "", 
          "output vhbv file for diploid assembly created by this program");
     CommandArgument_Bool_OrDefault_Doc(SAVE, False,
          "save fastb in canonical location");
     CommandArgument_Bool_OrDefault_Doc(OUR_REF, False,
          "use our reference sequence instead of the Gerstein one");
     CommandArgument_Int_OrDefault_Doc(MAX_MULT, 20,
          "maximum multiplicity for totally indeterminacy");
     EndCommandArguments;

     // Check arguments.

     ForceAssert( A_SHBV != "" ^ A_EFASTA != "" );
     if ( A_EFASTA != "" && K < 0 )
     {    FatalErr("When you specify A_EFASTA, you also need to specify "
                     "the K value for the assembly.\n"
                     "Abort.");    }

     // Load reference.

     vecbasevector R(2);
     vec<HyperBasevector> hb_R;
     int chr = X.Before( ":" ).Int( );
     String CHR = ( chr <= 21 ? ToString(chr+1) : X );
     int START = X.Between( ":", "-" ).Int( ), STOP = X.After( "-" ).Int( );
     const int endk = 100;
     basevector rleft, rright;
     if ( !OUR_REF )
     {    String drefdir = "/wga/scr4/ALLPATHS/H.sapiens.NA12878/other_ref";
          vecbasevector dref, ref;
          dref.ReadOne( drefdir + "/maternal.fastb", chr );
          dref.ReadOne( drefdir + "/paternal.fastb", chr );
          String refdir = "/wga/dev/references/Homo_sapiens";
          ref.ReadOne( refdir + "/genome.fastb", chr );
          rleft = basevector( ref[0], START, endk ); 
          rright = basevector( ref[0], STOP - endk, endk );
          String sm = dref[0].ToString( ), sp = dref[1].ToString( );
          String srleft = rleft.ToString( ), srright = rright.ToString( );
          int left_pos1 = sm.Position(srleft), left_pos2 = sp.Position(srleft);
          int right_pos1 = sm.Position(srright), right_pos2 = sp.Position(srright);
          if ( left_pos1 < 0 || left_pos2 < 0 || right_pos1 < 0 || right_pos2 < 0 )
          {    FatalErr("It would appear that there is a variant near the ends "
                         "of your region.  Try shifting the region ends.\n"
                         "Abort.");    }
          if ( left_pos1 >= right_pos1 || left_pos2 >= right_pos2 )
          {    FatalErr("Failed to map region.\n"
                          "Abort.");    }
          R[0].SetToSubOf( dref[0], left_pos1, right_pos1 - left_pos1 + endk );
          R[1].SetToSubOf( dref[1], left_pos2, right_pos2 - left_pos2 + endk );    }
     else
     {    String rdir = "/wga/dev/references/Homo_sapiens/NA12878_regions";
          vec<String> knowns = AllFiles(rdir);
          int p = BinPosition( knowns, X + ".vhbv" );
          // ***************************** BROKEN! *********************************
          if ( p < 0 )
          {    FatalErr("You've specified OUR_REF but I can't find our "
                         "reference sequence.\n"
                         "Abort.");    }
          BinaryReader::readFile( rdir + "/" + X + ".vhbv", &hb_R );
          if ( hb_R.size( ) != 2 )
          {    FatalErr("\nThe reference needs to consist of two "
                         "HyperBasevectors.\n"
                         "Abort.");    }
          for ( int j = 0; j < 2; j++ )
          {    vec<int> sources, sinks;
               hb_R[j].Sources(sources), hb_R[j].Sinks(sinks);
               if ( !sources.solo( ) || !sinks.solo( ) )
               {    FatalErr("\nEach half of the reference needs to have a "
                              "unique source and a unique sink.\n"
                              "Abort.");    }
               int vstart = sources[0], vstop = sinks[0];
               if ( !hb_R[j].From(vstart).solo( ) || !hb_R[j].To(vstop).solo( ) )
               {    FatalErr("\nI'm confused by the ends of your reference.\n"
                             "Abort.");    }
               const basevector& estart = hb_R[j].EdgeObjectByIndexFrom( vstart, 0 );
               const basevector& estop = hb_R[j].EdgeObjectByIndexTo( vstop, 0 );
               if ( estart.isize( ) < endk || estop.isize( ) < endk )
               {    FatalErr("\nThe ends of your reference sequence are too short.\n"
                              "Abort.");    }
               if ( j == 0 )
               {    rleft = basevector( estart, 0, endk );
                    rright = basevector( estop, estop.isize( ) - endk, endk );    }
               else
               {    if ( rleft != basevector( estart, 0, endk )
                         || rright != 
                         basevector( estop, estop.isize( ) - endk, endk ) )
                    {    FatalErr("It would appear that there is a variant "
                                   "near the ends of your region.\n"
                                   "Abort.");    }    }    }    }

     // Get efasta, either directly, or by converting the HyperBasevector.

     VecEFasta e;
     if ( A_EFASTA != "" ) LoadEfastaIntoStrings( A_EFASTA, e );
     else
     {    SupportedHyperBasevector shb;
          BinaryReader::readFile( A_SHBV, &shb );
          K = shb.K( );
          HyperEfasta he(shb);

          Bool reduce_verbose = True;
          long_logging logc;
          Reduce( he, reduce_verbose, logc );

          /*
          Ofstream( out, "woof.dot" );
          vec<double> lens( he.EdgeObjectCount( ) );
          for ( int i = 0; i < he.EdgeObjectCount( ); i++ )
               lens[i] = he.EdgeObject(i).size( );
          he.PrettyDOT( out, lens, digraphE<efasta>::edge_label_info::DEFAULT, 
               True, True );
          */

          cout << "reduced efasta has " << he.EdgeObjectCount( ) << " edges" << endl;
          if ( A_EFASTA_OUT != "" )
          {    Ofstream( dout, A_EFASTA_OUT + ".dot" );
               const Bool DOT_LABEL_CONTIGS = True;
               const Bool DOT_LABEL_VERTICES = False;
               vec<double> lengths( he.EdgeObjectCount( ) );
               for ( int i = 0; i < he.EdgeObjectCount( ); i++ )
                    lengths[i] = he.EdgeLengthKmers(i);
               vec<String> edge_id_names( he.EdgeObjectCount( ) );
               for ( int i = 0; i < he.EdgeObjectCount( ); i++ )
                    edge_id_names[i] = ToString(i);
               he.PrettyDOT( dout, lengths, HyperEfasta::edge_label_info(
                    HyperEfasta::edge_label_info::DIRECT, &edge_id_names ),
                    DOT_LABEL_CONTIGS, DOT_LABEL_VERTICES );
               Ofstream( out, A_EFASTA_OUT + ".efasta" )
               for ( int v = 0; v < he.N( ); v++ )
               {    for ( size_t j = 0; j < he.From(v).size(); j++ )
                    {    int e = he.EdgeObjectIndexByIndexFrom( v, j );
                         int w = he.From(v)[j];
                         he.EdgeObject(e).Print( 
                              out, ToString(e) + " [vert_" + ToString(v)
                              + "-->vert_" + ToString(w) + "]" );    }    }    }
          e.assign(he.Edges().begin(),he.Edges().end());    }

     // If two bracket expressions are separated by less than K bases, merge them.

     if ( e.size() == 1ul )
     {    while(1)
          {    Bool progress = False;
               vec< vec<String> > blocks;
               e[0].GetBlocks(blocks);
               for ( int i = 2; i < blocks.isize( ) - 2; i += 2 )
               {    if ( blocks[i].solo( ) && blocks[i][0].isize( ) < K )
                    {    vec<String> b3;
                         for ( int j1 = 0; j1 < blocks[i-1].isize( ); j1++ )
                         for ( int j2 = 0; j2 < blocks[i+1].isize( ); j2++ )
                         {    b3.push_back( blocks[i-1][j1] + blocks[i][0] +
                                   blocks[i+1][j2] );    }
                         blocks[i] = b3;
                         vec<Bool> to_delete( blocks.size( ), False );
                         to_delete[i-1] = to_delete[i+1] = True;
                         EraseIf( blocks, to_delete );
                         String n;
                         for ( int l = 0; l < blocks.isize( ); l++ )
                         {    if ( blocks[l].solo( ) ) n += blocks[l][0];
                              else
                              {    n += "{" + blocks[l][0];
                                   for ( int r = 1; r < blocks[l].isize( ); r++ )
                                        n += "," + blocks[l][r];
                                   n += "}";    }    }
                         e[0] = n;
                         progress = True;
                         break;    }    }
               if ( !progress ) break;    }    }

     // Check to see if efasta passes.

     const int max_mult_arb = 15;
     Bool /*simple = True,*/ gen_simple = True;
     vec< vec<basevector> > edges;
     if ( e.size() != 1ul ) gen_simple = False;
     else
     {    e[0].MakeGraph(edges);
          if ( edges.size( ) % 2 != 1 ) gen_simple = False;
          for ( int i = 0; i < edges.isize( ); i++ )
          {    if ( i % 2 == 0 && edges[i].size( ) != 1 ) gen_simple = False;
               if ( i % 2 == 0 && edges[i][0].isize( ) < K ) gen_simple = False;
               if ( i % 2 == 1 && edges[i].size( ) != 2 ) 
               {    if ( edges[i].size( ) == 1 ) gen_simple = False;
                    else if ( edges[i].isize( ) <= max_mult_arb ) gen_simple = True;
                    else
                    {    //simple = False;
                         Bool ok = False;
                         const int max_period = 2;
                         for ( int p = 1; p <= max_period; p++ )
                         {    String motif;
                              for ( int j = 0; j < edges[i].isize( ); j++ )
                              {    if ( edges[i][j].isize( ) >= p )
                                   {    motif = edges[i][j].ToString( );
                                        motif.resize(p);
                                        break;    }    }
                              if ( motif.empty( ) ) continue;
                              Bool eq = True;
                              for ( int j = 0; j < edges[i].isize( ); j++ )
                              {    if ( edges[i][j].isize( ) % p != 0 )
                                   {    eq = False;
                                        break;    }
                                   String s;
                                   for ( int l = 0; 
                                        l < edges[i][j].isize( ) / p; l++ )
                                   {    s += motif;    }
                                        if ( edges[i][j].ToString( ) != s )
                                   {    eq = False;
                                        break;    }    }
                              if (eq)
                              {    ok = True;
                                   break;    }    }
                         if ( !ok ) gen_simple = False;    }    }    }    }
     //if ( !gen_simple ) simple = False;
     int64_t gen_mult = 1;
     if (gen_simple)
     {    cout << "Assembly is generalized simple." << endl;
          Bool exceeded = False;
          for ( int i = 1; i < edges.isize( ); i += 2 )
          {    if ( edges[i].size( ) != 2 ) 
               {    gen_mult *= edges[i].size( );
                    if ( gen_mult > 1000 ) 
                    {    exceeded = True;
                         break;    }    }    }
          if (exceeded) cout << "multiplicity >= " << gen_mult << endl;
          else cout << "multiplicity = " << gen_mult << endl;    }
     else cout << "Assembly is not generalized simple." << endl;
     if ( !gen_simple || gen_mult > MAX_MULT )
     {    FatalErr("Assembly does not meet criteria.\n"
                    "Abort.");    }
     for ( int i = 1; i < edges.isize( ); i += 2 )
     {    for ( int j = 0; j < edges[i].isize( ); j++ )
          {    ForceAssertGe( edges[i-1][0].isize( ), K );
               ForceAssertGe( edges[i+1][0].isize( ), K );
               edges[i][j] = Cat( basevector( edges[i-1][0],
                         edges[i-1][0].isize( ) - (K-1), K-1 ), edges[i][j],
                    basevector( edges[i+1][0], 0, K-1 ) );    }    }

     // Locate reference ends in the assembly, reversing it if need be.

     int left_trim = -1, right_trim = -1;
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 )
          {    edges.ReverseMe( );
               for ( int i = 0; i < edges.isize( ); i++ )
               for ( int j = 0; j < edges[i].isize( ); j++ )
                    edges[i][j].ReverseComplement( );    }
          basevector aleft = edges.front( )[0], aright = edges.back( )[0];
          String saleft = aleft.ToString( ), saright = aright.ToString( );
          String srleft = rleft.ToString( ), srright = rright.ToString( );
          left_trim = -1, right_trim = -1;
          for ( int j = 0; j <= aleft.isize( ) - endk; j++ )
          {    if ( saleft.Contains( srleft, j ) )
               {    left_trim = j;
                    break;    }    }
          for ( int j = aright.isize( ) - endk; j >= 0; j-- )
          {    if ( saright.Contains( srright, j ) )
               {    right_trim = aright.isize( ) - endk - j;
                         break;    }    }
          PRINT3( pass, left_trim, right_trim );
          if ( left_trim >= 0 && right_trim >= 0 ) 
          {    edges.front( )[0].SetToSubOf( edges.front( )[0],
                    left_trim, edges.front( )[0].isize( ) - left_trim );
               edges.back( )[0].resize( edges.back( )[0].isize( ) - right_trim );
               break;    }    }
     if ( left_trim < 0 || right_trim < 0 )
     {    FatalErr("\nUnable to find the reference ends in the assembly.\n"
                    "Abort.");    }

     // Load mom and pop reads for region.

     String G1000 = "/humgen/gsa-hpprojects/NA12878Collection/bams";
     ForceAssert( chr <= 23 );
     for ( int pass = 1; pass <= 2; pass++ ) // mom, then pop
     {    String sample = ( pass == 1 ? "NA12892" : "NA12891" );
          String bam = "CEUTrio.HiSeq.WGS.b37_decoy." + sample
               + ".clean.dedup.recal.bam";
          String head = TMP + "/" + ( pass == 1 ? "mom" : "pop" );
          SystemSucceed( "( cd " + G1000 + "; samtools view " + bam + " " + CHR 
               + ":" + ToString(START) + "-" + ToString(STOP) 
               + " ) | SAM2CRDDump OUT_HEAD="
               + head + " LOG_TO_CERR=False > /dev/null" );    }
     vecbasevector mom( TMP + "/mom.fastb" ), pop( TMP + "/pop.fastb" );

     // Use mom and pop to phase bubbles.

     const int M = 40;
     vecbasevector parents(mom);
     parents.Append(pop);
     vec< triple<kmer<M>,int,int> > kmers_plus;
     MakeKmerLookup1( parents, kmers_plus );
     vec< kmer<M> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;
     for ( int i = 0; i < edges.isize( ); i++ )
     {    if ( edges[i].size( ) != 2 ) continue;
          cout << "\nbubble " << (i+1)/2 << endl;
          vec< vec<int> > count( 2, vec<int>(2,0) );
          for ( int j = 0; j < 2; j++ )
          {    const basevector& e = edges[i][j];
               kmer<M> x;
               for ( int l = 0; l <= e.isize( ) - M; l++ )
               {    int min_count = 1;
                    int momc = 0, popc = 0;
                    x.SetToSubOf( e, l );
                    for ( int xpass = 1; xpass <= 2; xpass++ )
                    {    if ( xpass == 2 ) x.ReverseComplement( );
                         int64_t low = LowerBound( kmers, x );
                         int64_t high = UpperBound( kmers, x );
                         for ( int64_t m = low; m < high; m++ )
                         {    if ( kmers_plus[m].second < (int) mom.size( ) ) momc++;
                              else popc++;    }    }
                    Bool momx = ( momc >= min_count );
                    Bool popx = ( popc >= min_count );
                    if (momx) count[j][0]++;
                    if (popx) count[j][1]++;
                    if ( momx && !popx ) cout << "M";
                    if ( !momx && popx ) cout << "P";
                    if ( momx && popx ) cout << "2";
                    if ( !momx && !popx ) cout << "?";    }
               cout << "\n";    }
          int score = count[0][0] + count[1][1] - count[0][1] - count[1][0];
          if ( score < 0 || ( score == 0 && randomx( ) % 2 == 0 ) )
               swap( edges[i][0], edges[i][1] );    }

     // Build separate chromosomes from the assembly.

     vec<HyperBasevector> hb_A(2);
     for ( int j = 0; j < 2; j++ )
     {    vec< vec<basevector> > edgesx(edges);
          for ( int i = 1; i < edges.isize( ); i += 2 )
          {    if ( edges[i].size( ) != 2 ) continue;
               edgesx[i].erase( edgesx[i].begin( ) + (1-j) );    }
          hb_A[j].SetK(K);
          ( (digraphE<basevector>&) hb_A[j] ).Initialize( edges.size( ) + 1 );
          for ( int i = 0; i < edgesx.isize( ); i++ )
          for ( int k = 0; k < edgesx[i].isize( ); k++ )
               hb_A[j].AddEdge( i, i + 1, edgesx[i][k] );    }
     HyperBasevector hb_A_fw( K, hb_A );

     // Load the reads.

     vecbasevector bases( TMP + "/frag_reads_orig.fastb" );
     vecqualvector quals( TMP + "/frag_reads_orig.qualb" );

     // Set up HyperBasevectors.

     HyperBasevector hb_Ap(hb_A_fw), hb_Rp;
     if ( !OUR_REF )
     {    vec<basevector> Rv(2);
          for ( int j = 0; j < 2; j++ )
               Rv[j] = R[j];
          hb_Rp = HyperBasevector( K, Rv );    }
     else hb_Rp = HyperBasevector( hb_R[0].K( ), hb_R );

     // Evaluate by reads.

     int assembly_count, reference_count;
     EvalByReads( hb_Ap, hb_Rp, bases, quals, assembly_count, reference_count, 
          SHOW_ALL, True );

     // Print summary stats.
     
     cout << "REPORT: favors_assembly = " << assembly_count
          << ", favors reference = " << reference_count << endl;
     if ( reference_count == 0 ) cout << "VALIDATES\n";

     // Output assembly.

     if ( OUT != "" ) BinaryWriter::writeFile( OUT, hb_A );
     if ( SAVE && reference_count == 0 )
     {    cout << "REPORT: saving" << endl;
          String dir = "/wga/dev/references/Homo_sapiens/NA12878_regions";
          vec<String> all = AllFiles(dir);
          Bool conflict = False;
          for ( int i = 0; i < all.isize( ); i++ )
          {    if ( !all[i].Contains( ".vhbv", -1 ) ) continue;
               if ( all[i].Before( ":" ).Int( ) != chr ) continue;
               int start = all[i].Between( ":", "-" ).Int( );
               int stop = all[i].Between( "-", ".vhbv" ).Int( );
               if ( start == START && stop == STOP ) 
               {    cout << "REPORT: duplicate region" << endl;
                    continue;    }
               if ( start >= START && stop <= STOP )
               {    cout << "REPORT: removing " << all[i] << endl;
                    SystemSucceed( "cd " + dir + "; svn remove " + all[i] );    }
               else if ( IntervalOverlap( start, stop, START, STOP ) > 0 )
               {    if ( stop - start < STOP - START )
                    {    cout << "REPORT: removing overlapping " << all[i] << endl;
                         SystemSucceed( 
                              "cd " + dir + "; svn remove " + all[i] );    }
                    else
                    {    cout << "REPORT: conflict, not saving" << endl;
                         conflict = True;
                         break;    }    }    }
          if ( !conflict )
          {    String region = ToString(chr) + ":" 
                    + ToString( double(START)/1000000 ) + "M-" 
                    + ToString( double(STOP)/1000000 ) + "M";
               String outfile = dir + "/" + region + ".vhbv";
               BinaryWriter::writeFile( outfile, hb_A );    }    }    }
