/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// SomaticCallSlave.  See SomaticCall for what this does.

#include "Basevector.h"
#include "Bitvector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "Vec.h"
#include "cancer/AssembleSNP.h"
#include "cancer/SomaticCallTools.h"
#include "lookup/LookAlign.h"
#include "lookup/SAM2CRD.h"
#include "random/Random.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/BinaryStream.h"

class fake_mutation {
     public:

     fake_mutation( ) { }
     fake_mutation( const int tig, const int pos, const char base,
          const double frac ) : tig(tig), pos(pos), base(base), frac(frac) { }

     friend Bool operator<( const fake_mutation& f1, const fake_mutation& f2 )
     {    if ( f1.tig < f2.tig ) return True;
          if ( f1.tig > f2.tig ) return False;
          if ( f1.pos < f2.pos ) return True;
          return False;    }

     int tig, pos;
     char base;
     double frac;

};

int QualDiff( const look_align& la, basevector& rd1, const basevector& rd2,
     qualvector& q1 )
{    int qual_diff = 0;
     if ( la.Rc1( ) ) { rd1.ReverseComplement( ); q1.ReverseMe( ); }
     int p1 = la.pos1( ), p2 = la.pos2( );
     for ( int j = 0; j < la.a.Nblocks( ); j++ )
     {    if ( la.a.Gaps(j) > 0 ) p2 += la.a.Gaps(j);
          if ( la.a.Gaps(j) < 0 ) p1 -= la.a.Gaps(j);
          for ( int x = 0; x < la.a.Lengths(j); x++ )
          {    if ( p2 >= 0 && p2 < rd2.isize( ) )
                    if ( rd1[p1] != rd2[p2] ) qual_diff += q1[p1];
               ++p1; ++p2;    }    }
     if ( la.Rc1( ) ) { rd1.ReverseComplement( ); q1.ReverseMe( ); }
     return qual_diff;    }

int Subs( const look_align& la, basevector& rd1, const basevector& rd2,
     qualvector& q1 )
{    int subs = 0;
     if ( la.Rc1( ) ) { rd1.ReverseComplement( ); q1.ReverseMe( ); }
     int p1 = la.pos1( ), p2 = la.pos2( );
     for ( int j = 0; j < la.a.Nblocks( ); j++ )
     {    if ( la.a.Gaps(j) > 0 ) p2 += la.a.Gaps(j);
          if ( la.a.Gaps(j) < 0 ) p1 -= la.a.Gaps(j);
          for ( int x = 0; x < la.a.Lengths(j); x++ )
          {    if ( p2 >= 0 && p2 < rd2.isize( ) )
                    if ( rd1[p1] != rd2[p2] ) subs++;
               ++p1; ++p2;    }    }
     if ( la.Rc1( ) ) { rd1.ReverseComplement( ); q1.ReverseMe( ); }
     return subs;    }

int main(int argc, char **argv)
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(REF);
     CommandArgument_String(REFLIST);
     CommandArgument_String(TUMOR_BAM);
     CommandArgument_String(NORMAL_BAM);
     CommandArgument_String(DBSNP);
     CommandArgument_Int(C);
     CommandArgument_Int(START);
     CommandArgument_Int(LEN);
     CommandArgument_Int(MAX_ERRS_TUMOR);
     CommandArgument_Int(MAX_ERRS_NORMAL);
     CommandArgument_Int(MAX_QUAL_DIFF);
     CommandArgument_String(OUT_HEAD);
     CommandArgument_Int(MIN_MUTANT_SUM_PRETEST);
     CommandArgument_Int(MIN_MUTANT_SUM);
     CommandArgument_Bool_OrDefault(DUMP_ALIGNS, False);
     CommandArgument_Double(TUMOR_THRESHOLD);
     CommandArgument_Double(NORMAL_THRESHOLD);
     CommandArgument_Bool_OrDefault(ASSEMBLE_SNP, True);
     CommandArgument_Bool_OrDefault(ASSEMBLE_SNP_VERBOSE, False);
     CommandArgument_String_OrDefault(FAKE_MUTATIONS, "");
     EndCommandArguments;

     // Start!

     cout << Date( ) << ": SomaticCallSlave starting" << endl;
     double clock = WallClockTime( );
     if ( OUT_HEAD.Contains( "/", -1 ) ) Mkdir777(OUT_HEAD);

     // Load reference dictionary.

     vec<String> refdict;
     vec<int> refsizes;
     String line;
     fast_ifstream dictin(REFLIST);
     while(1)
     {    getline( dictin, line );
          if ( dictin.fail( ) ) break;
          istrstream iline( line.c_str( ) );
          String chr;
          int n;
          iline >> chr >> n;
          refdict.push_back(chr);
          refsizes.push_back(n);    }

     // Define start and stop.

     int start = START, stop = START + LEN;
     int n = stop - start;
     String RANGES = ToString(C) + ":" + ToString(start) + "-" + ToString(stop);
     vec<int> CHR, POS;
     for ( int i = start; i < stop; i++ )
     {    CHR.push_back(C);
          POS.push_back(i);    }

     // Load fake mutations.

     vec<fake_mutation> fakes;
     if ( FAKE_MUTATIONS != "" )
     {    fast_ifstream fakein(FAKE_MUTATIONS);
          while(1)
          {    getline( fakein, line );
               if ( fakein.fail( ) ) break;
               int tig, pos;
               String base;
               double frac;
               istrstream iline( line.c_str( ) );
               iline >> tig >> pos >> base >> frac;
               if ( tig == C && pos >= start && pos < stop )
                    fakes.push( tig, pos, as_char(base[0]), frac );    }    }
     Sort(fakes);
     vec<int> fake_pos;
     for ( int i = 0; i < fakes.isize( ); i++ )
          fake_pos.push_back( fakes[i].pos );

     // Load reference and dbSNP.

     vecbasevector ref(REF);
     vecbitvector dbsnp(DBSNP);

     // Define region to be extracted.

     const int ext = 100; // WARNING: SHOULD NOT BE HARCODED!
     String region;
     if ( FAKE_MUTATIONS == "" )
     {    int low = Max( start-ext, 0 ), high = Min( stop+ext, refsizes[C] );
          region = refdict[C] + ":" + ToString(low) + "-" + ToString(high);    }
     else
     {    if ( fake_pos.empty( ) )
          {    cout << "No fake mutations in this region, exiting." << endl;
               Remove( OUT_HEAD + "mutations" );
               exit(0);    }
          for ( int j = 0; j < fake_pos.isize( ); j++ )
          {    if ( j > 0 ) region += " ";
               int s = fake_pos[j];
               int low = Max( s-ext, 0 ), high = Min( s+ext, refsizes[C] );
               region += refdict[C] + ":"
                    + ToString(low) + "-" + ToString(high);    }    }

     // Extract region.

     cout << Date( ) << ": extracting region" << endl;
     vecbasevector tbases, nbases;
     vecqualvector tquals, nquals;
     vec<look_align> taligns, naligns;
     vec<int> tmappingscore, nmappingscore;
     vec<look_align_x> talignsz, nalignsz;
     vecString tnames, nnames;
     {    vec<pairinfo> tpairs, npairs;
          {    SAM::BAMFile tumor_bam(TUMOR_BAM,region,false,Logger::nullLogger());
               vecString libNames;
               vec<Bool> first_in_pair;
               SAM2CRD(tumor_bam,tbases,tquals,talignsz,tpairs,tnames,first_in_pair,libNames);
               for ( int i = 0; i < talignsz.isize( ); i++ )
                    talignsz[i].target_id = C;
               for ( int i = 0; i < talignsz.isize( ); i++ )
               {    const look_align_x& la = talignsz[i];
                    int subs = Subs( la, tbases[ la.query_id ],
                         ref[ la.target_id ], tquals[ la.query_id ] );
                    if ( la.mapQ( ) < 10 || subs > MAX_ERRS_TUMOR ) continue;
                    if ( MAX_QUAL_DIFF > 0 )
                    {    int qual_diff = QualDiff( la, tbases[ la.query_id ],
                              ref[ la.target_id ], tquals[ la.query_id ] );
                         if ( qual_diff > MAX_QUAL_DIFF ) continue;    }
                    taligns.push_back(la);
                    tmappingscore.push_back( la.mapQ( ) );    }    }
          {    SAM::BAMFile normal_bam(NORMAL_BAM,region,false,Logger::nullLogger());
               vecString libNames;
               vec<Bool> first_in_pair;
               SAM2CRD(normal_bam,nbases,nquals,nalignsz,npairs,nnames,first_in_pair,libNames);
               for ( int i = 0; i < nalignsz.isize( ); i++ )
                    nalignsz[i].target_id = C;
               for ( int i = 0; i < nalignsz.isize( ); i++ )
               {    const look_align_x& la = nalignsz[i];
                    int subs = Subs( la, nbases[ la.query_id ],
                         ref[ la.target_id ], nquals[ la.query_id ] );
                    if ( subs > MAX_ERRS_NORMAL ) continue;
                    if ( MAX_QUAL_DIFF > 0 )
                    {    int qual_diff = QualDiff( la, nbases[ la.query_id ],
                              ref[ la.target_id ], nquals[ la.query_id ] );
                         if ( qual_diff > MAX_QUAL_DIFF ) continue;    }
                    naligns.push_back(la);
                    nmappingscore.push_back( la.mapQ( ) );    }    }
          cout << Date( ) << ": done" << endl;

          // If there are no reads, exit.

          if ( taligns.empty( ) && naligns.empty( ) )
          {    cout << "No reads, exiting" << endl;
               Remove( OUT_HEAD + "mutations" );
               exit(0);    }

          // Dump alignments.

          if (DUMP_ALIGNS)
          {    for ( int pass = 1; pass <= 2; pass++ )
               {    const vec<look_align>& aligns = (pass == 1 ? naligns : taligns);
                    const vecbasevector& bases = ( pass == 1 ? nbases : tbases );
                    const vec<pairinfo>& pairs
                         = ( pass == 1 ? npairs : tpairs );
                    const vecString& names = ( pass == 1 ? nnames : tnames );
                    vec<int> partner( bases.size( ), -1 );
                    for ( int i = 0; i < pairs.isize( ); i++ )
                    {    partner[ pairs[i].readID1 ] = pairs[i].readID2;
                         partner[ pairs[i].readID2 ] = pairs[i].readID1;    }
                    cout << "\n" << aligns.size( ) << " alignments of "
                         << ( pass == 1 ? "normal" : "tumor" ) << "\n" << endl;
                    for ( int i = 0; i < aligns.isize( ); i++ )
                    {    const look_align& la = aligns[i];
                         int id = la.query_id;
                         ForceAssertLe( la.Pos2( ), ref[ la.target_id ].isize( ) );
                         String title = "read_" + ToString(id) + " = " + names[id];
                         if ( partner[id] >= 0  )
                         {    title += " paired to read_" + ToString( partner[id] )
                                   + " = " + names[ partner[id] ];    }
                         bases[ la.query_id ].Print( cout, title );
                         la.PrintReadableBrief(cout);
                         la.PrintVisual( cout,
                              bases[id], ref[ la.target_id ] );    }    }    }    }

     // Combine data for the tumor and the normal.  This is very confusing: we
     // append normal data to tumor data.

     int ntaligns = taligns.size( ), nnaligns = naligns.size( );
     int ntbases = tbases.size( );
     vecbasevector& allbases = tbases;
     allbases.Append(nbases);
     vecqualvector& allquals = tquals;
     allquals.Append(nquals);
     vec<look_align>& allaligns = taligns;
     for ( int i = 0; i < nnaligns; i++ )
     {    look_align la = naligns[i];
          la.query_id += ntbases;
          allaligns.push_back(la);    }

     // Get maximum read length.

     int maxread = 0;
     for ( int i = 0; i < taligns.isize( ); i++ )
          maxread = Max( maxread, tbases[ taligns[i].query_id ].isize( ) );

     // Introduce fake mutations.

     if ( FAKE_MUTATIONS != "" )
     {    for ( int i = 0; i < ntaligns; i++ )
          {    const look_align& la = taligns[i];
               basevector& rd1 = tbases[ la.query_id ];
               qualvector& q1 = tquals[ la.query_id ];
               if ( la.Rc1( ) ) { rd1.ReverseComplement( ); q1.ReverseMe( ); }
               int p1 = la.pos1( ), p2 = la.pos2( );
               for ( int j = 0; j < la.a.Nblocks( ); j++ )
               {    if ( la.a.Gaps(j) > 0 ) p2 += la.a.Gaps(j);
                    if ( la.a.Gaps(j) < 0 ) p1 -= la.a.Gaps(j);
                    for ( int x = 0; x < la.a.Lengths(j); x++ )
                    {    if ( p2 >= start && p2 < stop )
                         {    int P = BinPosition( fake_pos, p2 );
                              if ( P >= 0 )
                              {    double prob = fakes[P].frac;
                                   int prob_ceil = int( round( prob*1000000.0 ) );
                                   if ( randomx( ) % 1000000 <= prob_ceil )
                                   {    double prob_of_err
                                             = pow( 10.0, -float(q1[p1])/10.0 );
                                        Bool err = ( randomx( ) % 1000000 )
                                             < ( 1000000.0 * prob_of_err );
                                        char base = fakes[P].base;
                                        if (err)
                                        {    base += ( randomx( ) % 3 ) + 1;
                                             base = base % 4;    }
                                        rd1.Set( p1, base );    }    }    }
                         ++p1; ++p2;    }    }
               if ( la.Rc1( ) )
               { rd1.ReverseComplement( ); q1.ReverseMe( ); }    }    }

     // Find the sum of the quality scores for non-reference bases in the tumor.
     // This is used to speed up the subsequent code.

     vec<int> TQsum(n, 0);
     for ( int i = 0; i < ntaligns; i++ )
     {    const look_align& la = taligns[i];
          basevector rd1 = tbases[ la.query_id ];
          qualvector q1 = tquals[ la.query_id ];
          if ( la.Rc1( ) )
          {    rd1.ReverseComplement( );
               q1.ReverseMe( );    }
          int p1 = la.pos1( ), p2 = la.pos2( );
          for ( int j = 0; j < la.a.Nblocks( ); j++ )
          {    if ( la.a.Gaps(j) > 0 ) p2 += la.a.Gaps(j);
               if ( la.a.Gaps(j) < 0 ) p1 -= la.a.Gaps(j);
               for ( int x = 0; x < la.a.Lengths(j); x++ )
               {    if ( p2 >= start && p2 < stop && rd1[p1] != ref[C][p2] )
                         TQsum[p2-start] += q1[p1];
                    ++p1; ++p2;    }    }    }

     // Build quality score vectors.

     vec< vec<unsigned char> > TQ(4*n), NQ(4*n);
     vec< vec<unsigned char> > TQ_rc(4*n), NQ_rc(4*n);
     cout << Date( ) << ": start building Q score vectors" << endl;
     for ( int pass = 1; pass <= 2; pass++ )
     {    vec< vec<unsigned char> >& Q = ( pass == 1 ? TQ : NQ );
          vec< vec<unsigned char> >& Q_rc = ( pass == 1 ? TQ_rc : NQ_rc );
          const vec<look_align>& aligns = ( pass == 1 ? taligns : naligns );
          const vecbasevector& bases = ( pass == 1 ? tbases : nbases );
          const vecqualvector& quals = ( pass == 1 ? tquals : nquals );
          for ( int i = 0; i < ( pass == 1 ? ntaligns : nnaligns ); i++ )
          {    const look_align& la = aligns[i];
               int p1 = la.pos1( ), p2 = la.pos2( );
               int P2 = la.Pos2( );
               Bool have_mutant = False;
               for ( int j = Max( p2, start ); j < Min( P2, stop ); j++ )
               {    if ( TQsum[j-start] >= MIN_MUTANT_SUM_PRETEST )
                    {    if ( FAKE_MUTATIONS != "" )
                         {    if ( BinMember( fake_pos, j ) )
                              {    have_mutant = True;
                                   break;    }    }
                         else
                         {    have_mutant = True;
                              break;    }    }    }
               if ( !have_mutant ) continue;
               basevector rd1 = bases[ la.query_id ];
               qualvector q1 = quals[ la.query_id ];
               if ( la.Rc1( ) )
               {    rd1.ReverseComplement( );
                    q1.ReverseMe( );    }
               for ( int j = 0; j < la.a.Nblocks( ); j++ )
               {    if ( la.a.Gaps(j) > 0 ) p2 += la.a.Gaps(j);
                    if ( la.a.Gaps(j) < 0 ) p1 -= la.a.Gaps(j);
                    for ( int x = 0; x < la.a.Lengths(j); x++ )
                    {    if ( p2 >= start && p2 < stop )
                         {    Q[ 4*(p2-start) + rd1[p1] ].push_back( q1[p1] );
                              Q_rc[ 4*(p2-start) + rd1[p1] ].
                                   push_back( la.Rc1( ) );    }
                         ++p1; ++p2;    }    }    }    }
     cout << Date( ) << ": done" << endl;

     // Look for interesting positions.  All tumor reads covering a putative
     // mutation are flagged (tbases_keep) and those carrying the mutant allele
     // are tagged specially (tbases_carrier).

     vec< pair<int,int> > MUTATIONS;
     vec<String> mutation_reports0, mutation_reports1;
     vec<int> tbases_keep, tbases_carrier, nbases_keep;
     vec<unsigned char> altbases;
     ForceAssertEq( n, CHR.isize( ) );
     for ( int i = 0; i < n; i++ )
     {
          // Pretest for efficiency.

          if ( TQsum[i] < MIN_MUTANT_SUM_PRETEST ) continue;

          // Test for fake run.

          int chr = CHR[i], pos = POS[i];
          if ( FAKE_MUTATIONS != "" && !BinMember( fake_pos, pos ) ) continue;

          // Reject dbSNP positions.

          if ( dbsnp[chr][pos] ) continue;

          // Carry out main test.

          unsigned char refbase = ref[chr][pos];
          String info;
          unsigned char altbase = 4;
          if ( !SomaticMutation( chr, pos, ref, &TQ[4*i], &NQ[4*i], &TQ_rc[4*i],
               &NQ_rc[4*i], refbase, MIN_MUTANT_SUM_PRETEST, MIN_MUTANT_SUM,
               TUMOR_THRESHOLD, NORMAL_THRESHOLD, altbase, info ) )
          {    continue;    }
          mutation_reports0.push_back(info);

          // Check to see if locus is assemblable.

          if ( ASSEMBLE_SNP && !AssembleSNP( maxread, allbases, allquals, allaligns,
               ref, chr, pos, altbase, ASSEMBLE_SNP_VERBOSE ) )
          {    continue;    }

          // Mark the tumor and normal reads covering the mutation.

          for ( int j = 0; j < ntaligns; j++ )
          {    const look_align& la = taligns[j];
               if ( la.target_id == chr && la.pos2( ) <= pos && pos < la.Pos2( ) )
               {    tbases_keep.push_back( la.query_id );
                    basevector rd1 = tbases[ la.query_id ];
                    qualvector q1 = tquals[ la.query_id ];
                    if ( la.Rc1( ) )
                    {    rd1.ReverseComplement( );
                         q1.ReverseMe( );    }
                    int p1 = la.pos1( ), p2 = la.pos2( );
                    for ( int j = 0; j < la.a.Nblocks( ); j++ )
                    {    if ( la.a.Gaps(j) > 0 ) p2 += la.a.Gaps(j);
                         if ( la.a.Gaps(j) < 0 ) p1 -= la.a.Gaps(j);
                         for ( int x = 0; x < la.a.Lengths(j); x++ )
                         {    if ( p2 >= start && p2 < stop
                                   && rd1[p1] == altbase )
                              {    tbases_carrier.push_back( la.query_id );
                                   break;    }
                              ++p1; ++p2;    }    }    }    }
          for ( int j = 0; j < nnaligns; j++ )
          {    const look_align& la = naligns[j];
               if ( la.target_id == chr && la.pos2( ) <= pos && pos < la.Pos2( ) )
                    nbases_keep.push_back( la.query_id );    }

          // Print mutation data.

          MUTATIONS.push( chr, pos );
          altbases.push_back(altbase);
          mutation_reports1.push_back(info);    }

     // Write bases etc.

     BinaryWriter::writeFile( OUT_HEAD + "mutations", MUTATIONS );
     WriteStrings( OUT_HEAD + "mutation_reports0", mutation_reports0 );
     WriteStrings( OUT_HEAD + "mutation_reports1", mutation_reports1 );

     if ( true )
     {
         UniqueSort(tbases_keep);
         String tbasesFilNam(OUT_HEAD+"tbases");
         IncrementalWriter<bvec> tbases_writer( tbasesFilNam.c_str() );
         String tqualsFilNam(OUT_HEAD+"tquals");
         IncrementalWriter<qvec> tquals_writer( tqualsFilNam.c_str() );
         String tnamesFilNam(OUT_HEAD+"tnames");
         ofstream tnames_writer( tnamesFilNam.c_str() );
         vec<int>::iterator end(tbases_keep.end());
         for ( vec<int>::iterator itr(tbases_keep.begin()); itr != end; ++itr )
         {
             int idx = *itr;
             tbases_writer.add(tbases[idx]);
             tquals_writer.add(tquals[idx]);
             tnames_writer << tnames[idx] << '\n';
         }
         tbases_writer.close();
         tquals_writer.close();
         tnames_writer.close();
     }

     if ( true )
     {
         UniqueSort(nbases_keep);
         String nbasesFilNam(OUT_HEAD+"nbases");
         IncrementalWriter<bvec> nbases_writer( nbasesFilNam.c_str() );
         String nqualsFilNam(OUT_HEAD+"nquals");
         IncrementalWriter<qvec> nquals_writer( nqualsFilNam.c_str() );
         String nnamesFilNam(OUT_HEAD+"nnames");
         ofstream nnames_writer( nnamesFilNam.c_str() );
         vec<int>::iterator end(nbases_keep.end());
         for ( vec<int>::iterator itr(nbases_keep.begin()); itr != end; ++itr )
         {
             int idx = *itr;
             nbases_writer.add(nbases[idx]);
             nquals_writer.add(nquals[idx]);
             nnames_writer << nnames[idx] << '\n';
         }
         nbases_writer.close();
         nquals_writer.close();
         nnames_writer.close();
     }

     if ( true )
     {
         UniqueSort(tbases_carrier);
         vec<Bool> thot( tbases_keep.size( ), False );
         for ( int i = 0; i < tbases_carrier.isize( ); i++ )
              thot[ BinPosition( tbases_keep, tbases_carrier[i] ) ] = True;
         BinaryWriter::writeFile( OUT_HEAD + "thot", thot );
     }

     BinaryWriter::writeFile( OUT_HEAD + "altbase", altbases );
     Ofstream( taout, OUT_HEAD + "tqltout" );
     vec<int> tmappingscore_keep, nmappingscore_keep;
     for ( int i = 0; i < ntaligns; i++ )
     {    look_align la = taligns[i];
          int p = BinPosition( tbases_keep, la.query_id );
          if ( p < 0 ) continue;
          la.query_id = p;
          la.PrintParseable(taout);
          tmappingscore_keep.push_back( tmappingscore[i] );    }
     BinaryWriter::writeFile( OUT_HEAD + "tmappingscore", tmappingscore_keep );
     Ofstream( naout, OUT_HEAD + "nqltout" );
     for ( int i = 0; i < nnaligns; i++ )
     {    look_align la = naligns[i];
          int p = BinPosition( nbases_keep, la.query_id );
          if ( p < 0 ) continue;
          la.query_id = p;
          la.PrintParseable(naout);
          nmappingscore_keep.push_back( nmappingscore[i] );    }
     BinaryWriter::writeFile( OUT_HEAD + "nmappingscore", nmappingscore_keep );
     cout << Date( ) << ": SomaticCallSlave finished, elapsed time = "
          << TimeSince(clock) << endl;    }
