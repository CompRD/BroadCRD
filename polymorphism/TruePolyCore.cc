/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// This is the guts of TruePoly.  See TruePoly.cc for what it does.

#include "Alignment.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "ParseSet.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "polymorphism/TruePolyCore.h"

class alignpart {

     public:

     int pos1;
     int pos2;
     int gap;

     alignpart( ) { }
     alignpart( const int pos1, const int pos2, const int gap )
          : pos1(pos1), pos2(pos2), gap(gap) { }

};

void TruePolyCore( vec<poly>& polys, vecbasevector refA, 
     const vecbasevector& refB, const String& TRUSTEDA, const String& TRUSTEDB, 
     const int K, const String& CLASSES_TO_PRINT, vec< vec<poly> >* compound_polys, 
     vecbitvector* untrusted_bases, vecbitvector* nonpolymorphic_bases )
{
     bool doing_compound = false;
     if (compound_polys != NULL) doing_compound = true;
     vec<poly> singleton_polys;

     bool doing_untrusted_bases = false;
     if (untrusted_bases != NULL) 
     {    doing_untrusted_bases = true;
          untrusted_bases->clear();
          Mimic(refB, *untrusted_bases); }

     if ( nonpolymorphic_bases != NULL ) Mimic( refB, *nonpolymorphic_bases );

     // Define size of flanking sequences.

     int flank = Min( K, 20 );

     // Parse print classes.

     Bool print1 = False, print2 = False, print3 = False;
     vec<int> to_print;
     ParseIntSet( CLASSES_TO_PRINT, to_print );
     if ( Member( to_print, 1 ) ) print1 = True;
     if ( Member( to_print, 2 ) ) print2 = True;
     if ( Member( to_print, 3 ) ) print3 = True;

     // Load trusted files.

     int nA = refA.size( ), nB = refB.size( );
     vecbitvector trustedA, trustedB;
     if ( TRUSTEDA != "" ) trustedA.ReadAll(TRUSTEDA);
     if ( TRUSTEDB != "" ) trustedB.ReadAll(TRUSTEDB);

     // Set up counters for event types.

     int n1 = 0, n2 = 0, n3 = 0;

     // To avoid going crazy, we do two passes.  On the first pass, we use refA,
     // and on the second pass its reverse.  We require that alignments are forward
     // on refB.

     polys.clear( );
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 ) ReverseComplement(refA);

          // Build combined reference and paths database for it.

          vecbasevector ref(refA);
          ref.Append(refB);
          vecKmerPath paths, paths_rc;
          vec<big_tagged_rpint> pathsdb;
          ReadsToPathsCoreY( ref, K, paths );
          paths_rc = paths;
          for ( size_t i = 0; i < paths_rc.size( ); i++ )
               paths_rc[i].Reverse( );
          CreateDatabase( paths, paths_rc, pathsdb );

          // Precompute path sums.
     
          vec< vec<int> > path_sum( ref.size( ) );
          for ( size_t i = 0; i < ref.size( ); i++ )
          {    path_sum[i].resize( paths[i].NSegments( ) + 1 );
               path_sum[i][0] = 0;
               for ( int j = 0; j < paths[i].NSegments( ); j++ )
                    path_sum[i][j+1] = path_sum[i][j] + paths[i].Length(j);    }

          // Find polymorphisms and non-polymorphisms.

          for ( int tigA = 0; tigA < nA; tigA++ )
          {    const KmerPath& p1 = paths[tigA];
          
               // We go through every kmer x1 on p1.

               KmerPathLoc l1 = p1.Begin( ), l2;
               for ( int apos1 = 0; apos1 <= refA[tigA].isize( ) - 2*K - 1; apos1++ )
               {    if ( apos1 > 0 ) l1 = l1 + 1;
                    longlong x1 = l1.GetKmer( );

                    // Now we have x1.  Find all occurrences of it in both 
                    // references.  There must be exactly two, one on one reference,
                    // one on the other.  We put them in order.

                    static vec<longlong> con;
                    Contains( pathsdb, x1, con );
                    if ( con.size( ) != 2 ) continue;
                    if ( pathsdb[ con[1] ].ReadId( ) < nA ) swap( con[0], con[1] );
                    if ( pathsdb[ con[1] ].ReadId( ) < nA ) continue;
                    big_tagged_rpint t1 = pathsdb[ con[0] ], u1 = pathsdb[ con[1] ];
                    int tigB = u1.ReadId( ) - nA;
                    if ( u1.PathId( ) < 0 ) continue;

                    // Find non-polymorphisms.

                    if ( nonpolymorphic_bases != 0 )
                    {    KmerPathLoc l1x = l1;
                         int jx;
                         for ( jx = 1; jx < K + 1; jx++ )
                         {    l1x = l1x + 1;
                              longlong x1x = l1x.GetKmer( );
                              static vec<longlong> conx;
                              Contains( pathsdb, x1x, conx );
                              if ( conx.size( ) != 2 ) break;
                              if ( pathsdb[ conx[1] ].ReadId( ) < nA ) 
                                   swap( conx[0], conx[1] );
                              if ( !( pathsdb[ conx[0] ].ReadId( ) < nA ) ) break;
                              if ( pathsdb[ conx[1] ].ReadId( ) < nA ) break;    }
                         if ( jx == K + 1 )
                         {    longlong bpos1 = x1 - u1.Start( )
                                   + path_sum[nA+tigB][ u1.PathPos( ) ];
                              int z;
                              for ( z = K; z < 2*K + 1; z++ )
                              {    if ( refA[tigA][apos1+z] != refB[tigB][bpos1+z] )
                                        break;    }
                              if ( z == 2*K + 1 )
                              {    if ( ( TRUSTEDA == "" || trustedA[tigA][apos1+K] )
                                        && ( TRUSTEDB == "" 
                                        || trustedB[tigB][bpos1+K] ) )
                                   {    (*nonpolymorphic_bases)[tigB].Set( bpos1+K,
                                             True );    }    }    }    }

                    // Now return to finding the polymorphisms.  Find next kmer x2 
                    // on REFA that is in REFB.  It must occur exactly once on each 
                    // reference.

                    Bool found_good_next = False;
                    l2 = l1;
                    big_tagged_rpint t2, u2;
                    longlong x2 = 0;
                    int middleA = -1;
                    Bool overlapping = False;
                    advance: 
                    while(1)
                    {    if ( l2.atEnd( ) ) break;
                         l2 = l2 + 1;
                         ++middleA;
                         x2 = l2.GetKmer( );
                         Contains( pathsdb, x2, con );
                         if ( con.size( ) > 2 ) break;
                         if ( con.size( ) == 1 ) continue;
                         if ( pathsdb[ con[1] ].ReadId( ) < nA ) 
                              swap( con[0], con[1] );
                         t2 = pathsdb[ con[0] ], u2 = pathsdb[ con[1] ];
                         if ( !( t2.ReadId( ) < nA ) ) break;
                         if ( u2.ReadId( ) - nA != tigB ) break;
                         if ( u2.PathId( ) < 0 ) continue;
                         found_good_next = True;
                         break;    }
                    if ( !found_good_next ) continue;
                    
                    // Define positions of x1 and x2 on REFB, verify that they are 
                    // in the correct order.

                    longlong bpos1 = x1 - u1.Start( )
                         + path_sum[nA+tigB][ u1.PathPos( ) ];
                    longlong bpos2 = x2 - u2.Start( )
                         + path_sum[nA+tigB][ u2.PathPos( ) ];
                    if ( !( bpos1 <= bpos2 ) ) continue;
                    if ( middleA == 0 && bpos2 - bpos1 == 1 ) continue;

                    // Verify that kmers between x1 and x2 on REFB do not occur in 
                    // REFA and occur exactly once in REFB.
                    
                    if ( !overlapping )
                    {    KmerPathLoc ib( 
                              paths[ tigB + nA ], u1.PathPos( ), x1 - u1.Start( ) );
                         Bool bad_intermediate = False;
                         while(1)
                         {    ib = ib + 1;
                              longlong z = ib.GetKmer( );
                              if ( z == x2 ) break;
                              Contains( pathsdb, z, con );
                              if ( !con.solo( ) )
                              {    bad_intermediate = True;
                                   break;    }    }
                         if (bad_intermediate) continue;    }

                    // Make sure we're not overlapping.

                    if ( bpos2 - bpos1 < K || middleA < K - 1 )
                    {    found_good_next = False;
                         overlapping = True;
                         goto advance;    }

                    // Make sure event is trusted.
                    
                    // If doing untrusted_snps, check if a SNP, then mark 
                    // if not going to be used because untrusted.

                    if ( doing_untrusted_bases && bpos2 - bpos1 >= K + 1) 
                    {    Bool trusted = True;
                         if ( TRUSTEDA != "" )
                         {    for ( int j = -1; j <= middleA - K + 1; j++ )
                                   if ( !trustedA[tigA][apos1+K+j] ) 
                                        trusted = False;    }
                         if ( TRUSTEDB != "" )
                         {    for ( int j = K - 1; j <= bpos2 - bpos1; j++ )
                                   if ( !trustedB[tigB][bpos1+j] ) 
                                        trusted = False;    }
                         if ( !trusted ) 
                         {     for ( int j = K; j < bpos2 - bpos1; j++ )
                                    (*untrusted_bases)[tigB].Set(bpos1 + j, True);
                                    }    }

                    if ( TRUSTEDA != "" )
                    {    Bool trusted = True;
                         for ( int j = -1; j <= middleA - K + 1; j++ )
                              if ( !trustedA[tigA][apos1+K+j] ) trusted = False;
                         if ( !trusted ) continue;    }
                    if ( TRUSTEDB != "" )
                    {    Bool trusted = True;
                         for ( int j = K - 1; j <= bpos2 - bpos1; j++ )
                              if ( !trustedB[tigB][bpos1+j] ) trusted = False;
                         if ( !trusted ) continue;    }
     
                    // Now we have a polymorphism event.  Output it.
                    
                    if ( middleA - K + 1 == 1 && bpos2 - bpos1 == K + 1 ) 
                    {    ++n1;
                         if ( !print1 ) continue;    }
                    else if ( middleA - K + 1 == 0 ) 
                    {    ++n2;
                         if ( !print2 ) continue;    }
                    else if ( bpos2 - bpos1 == K ) 
                    {    ++n3;
                         if ( !print3 ) continue;    }

                    String left, sample, ref, right;
                    for ( int j = 0; j < flank; j++ )
                         left += as_base( refB[tigB][ (K-flank) + bpos1 + j ] );
                    for ( int j = 0; j < middleA - K + 1; j++ )
                         sample += as_base( refA[tigA][ apos1 + K + j ] );
                    for ( int j = K; j < bpos2 - bpos1; j++ )
                         ref += as_base( refB[tigB][ bpos1 + j ] );
                    for ( int j = 0; j < flank; j++ )
                         right += as_base( refB[tigB][ bpos2 + j ] );
                    polys.push( 
                         tigB, bpos1 + K, left, sample, ref, right, bpos1 + K,
                         apos1 + K, tigA, sample.size( ), 
                         ref.size( ) );    
                    if (doing_compound) 
                    {    singleton_polys.push( tigB, bpos1 + K, left, sample, ref, 
                              right, bpos1 + K, apos1 + K, tigA, sample.size( ), 
                           ref.size( ) );    }    

                    // Mark flanking bases as non-polymorphism bases.

                    if ( nonpolymorphic_bases != 0 )
                    {    for ( int j = 0; j < K; j++ )
                         {    if ( ( TRUSTEDA == "" || trustedA[tigA][ apos1 + j ] )
                                   && ( TRUSTEDB == "" 
                                   || trustedB[tigB][ bpos1 + j ] ) )
                              {    (*nonpolymorphic_bases)[tigB].Set( bpos1+j,
                                        True );    }    }
                         for ( int j = 0; j < K; j++ )
                         {    if ( ( TRUSTEDA == "" ||
                                   trustedA[tigA][ apos1 + middleA + 1 + j ] )
                                   && ( TRUSTEDB == "" 
                                   || trustedB[tigB][ bpos2 + j ] ) )
                              {    (*nonpolymorphic_bases)[tigB].Set( bpos2+j,
                                        True );    }    }    }    }    }    }

     // Split composite polymorphism events.

     int npolys = polys.size( );
     vec<Bool> to_delete( npolys, False );
     vec<Bool> to_delete_singleton( npolys, False );
     Bool split_verbose = False; // TURN ON FOR DEBUGGING
     for ( int i = 0; i < npolys; i++ )
     {    poly p = polys[i];
          if ( p.sample.size( ) == 1 && p.ref.size( ) == 1 ) continue;
          if ( p.sample.size( ) == 0 || p.ref.size( ) == 0 ) continue;
          if (split_verbose) cout << "\nsplit " << polys[i];
          basevector S, T;
          S.SetFromString(p.sample), T.SetFromString(p.ref);
          Bool swapped = False;
          if ( S.size( ) > T.size( ) )
          {    swap( S, T );
               swapped = True;    }
          int best_loc;
          alignment al;
          SmithWatFree( S, T, best_loc, al, True, True );
          align a(al);
          vec<alignpart> parts;
          if ( a.pos1( ) > 0 ) parts.push( 0, 0, a.pos1( ) );
          if ( a.pos2( ) > 0 ) parts.push( 0, 0, -a.pos2( ) );
          int j, p1 = a.pos1( ), p2 = a.pos2( );
          for ( j = 0; j < a.Nblocks( ); j++ )
          {    if ( a.Gaps(j) != 0 ) parts.push( p1, p2, -a.Gaps(j) );
               if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
               if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
               for ( int x = 0; x < a.Lengths(j); x++ )
               {    if ( S[p1] != T[p2] ) parts.push( p1, p2, 0 );
                    ++p1;
                    ++p2;    }    }
          if ( a.Pos1( ) < S.isize( ) )
               parts.push( a.Pos1( ), a.Pos2( ), S.isize( ) - a.Pos1( ) );
          if ( a.Pos2( ) < T.isize( ) )
               parts.push( a.Pos1( ), a.Pos2( ), a.Pos2( ) - T.isize( ) );
          vec<Bool> parts_remove( parts.size( ), False );
          for ( int x = 0; x < parts.isize( ); x++ )
          {    if (swapped)
               {    swap( parts[x].pos1, parts[x].pos2 );
                    parts[x].gap = -parts[x].gap;    }
               if (split_verbose)
                    PRINT3( parts[x].pos1, parts[x].pos2, parts[x].gap );
               if ( !print1 && parts[x].gap == 0 ) parts_remove[x] = True;
               if ( !print2 && parts[x].gap > 0 ) parts_remove[x] = True;
               if ( !print3 && parts[x].gap < 0 ) parts_remove[x] = True;    }
          EraseIf( parts, parts_remove );
          if ( parts.empty( ) )
          {    to_delete[i] = True;
               if (doing_compound) to_delete_singleton[i] = True;
               continue;    }
          vec<poly> temp_compound;  // already cleared.
          poly org_event;
          for ( int x = 0; x < parts.isize( ); x++ )
          {    const alignpart& ap = parts[x];
               int start1 = ap.pos1, start2 = ap.pos2;
               int stop1, stop2;
               if ( ap.gap == 0 )
               {    stop1 = start1 + 1, stop2 = start2 + 1;    }
               else if ( ap.gap > 0 )
               {    stop1 = start1 + ap.gap, stop2 = start2;    }
               else 
               {    stop1 = start1, stop2 = start2 - ap.gap;    }
               String sampleplus = p.left + p.sample + p.right;
               String refplus = p.left + p.ref + p.right;
               poly pd( p.tig, p.ref_start + start2, 
                    refplus.substr( start2, flank ),
                    sampleplus.substr( start1 + flank, stop1 - start1 ),
                    refplus.substr( start2 + flank, stop2 - start2 ),
                    refplus.substr( stop2 + flank, flank ), p.orig_ref_start,
                    p.orig_sample_start, p.sample_tig, p.orig_sample_len,
                    p.orig_ref_len );
               if (split_verbose) cout << "found " << pd;

               if (doing_compound) {
                    if (x == 0) {
                         to_delete_singleton[i] = True;
                         org_event = polys[i];
                    }
                    temp_compound.push_back(pd);
               }
               if ( x == 0 ) {
                    polys[i] = pd;
               } else {    
                    polys.push_back(pd);    
                    to_delete.push_back(False);    
               }    
          } 
          if (doing_compound) {
               Sort(temp_compound);
               // Below done so original event is at the beginning.
               // Sort is not guaranteed to be stable, so cannot assume it
               // will stay first.
               temp_compound.insert(temp_compound.begin(), org_event);
               compound_polys->push_back(temp_compound);
          }

     }
     if (doing_compound) {
        EraseIf( singleton_polys, to_delete_singleton );
        for (unsigned int ii = 0; ii < singleton_polys.size( ); ++ii) {
           vec<poly> temp_compound;
           temp_compound.clear( );
           // Push back twice since original event 
           // is also the compound event
           temp_compound.push_back(singleton_polys[ii]);
           temp_compound.push_back(singleton_polys[ii]);
           compound_polys->push_back(temp_compound);
        }
     }
     EraseIf( polys, to_delete );
     Sort(polys);    }

void TruePolyCore( vec<poly>& polys, const String& REFA, const String& REFB, 
     const String& TRUSTEDA, const String& TRUSTEDB, const int K, 
     const String& CLASSES_TO_PRINT,
     vec< vec<poly> > *compound_polys, vecbitvector *untrusted_bases,
     vecbitvector* nonpolymorphic_bases )
{    vecbasevector refA, refB;
     if ( REFA.Contains( ".fastb", -1 ) ) refA.ReadAll(REFA);
     else FetchReads( refA, 0, REFA );
     if ( REFB.Contains( ".fastb", -1 ) ) refB.ReadAll(REFB);
     else FetchReads( refB, 0, REFB );
     TruePolyCore( polys, refA, refB, TRUSTEDA, TRUSTEDB, K, CLASSES_TO_PRINT, 
          compound_polys, untrusted_bases, nonpolymorphic_bases );    }

void SummarizeClasses( ostream& out, const vec<poly>& polys )
{    vec<int> C(5,0);
     for ( int i = 0; i < polys.isize( ); i++ )
          ++C[ polys[i].Class( ) ];
     out << "#" << C[1] << " SNPs\n";
     out << "#" << C[2] << " deletions from REFA\n";
     out << "#" << C[3] << " deletions from REFB\n";    }

void DumpHeader( ostream& out, const String& REFB )
{    out << "#" << REFB << "\n";
     if ( !REFB.Contains( ".fastb", -1 ) ) 
     {    fast_ifstream inB(REFB);
          String line;
          while(1)
          {    getline( inB, line );
               if ( inB.fail( ) ) break;
               if ( line.Contains( ">", 0 ) ) out << line << "\n";    }    }    }
