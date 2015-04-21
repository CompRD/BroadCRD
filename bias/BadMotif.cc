/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// BadMotif.  This program is a tool to help characterize bases on a genome that
// sequence poorly.  The user specifies a sequence motif, a reference sequence
// for the genome, together with a labeling of certain bases as bad (e.g. from
// BadCoverage).  The program finds all bases on the genome that match the
// motif, then finds those matching bases that are bad.  Two statistics are
// reported:
// (i) the number of matching bases that are bad;
// (ii) the fraction of matching bases that are bad.
//
// A motif consists of a sequence of A, C, G, T, or N, possibly intermixed with
// "fixed length base content restrictions", for example (G+C >= 0.8)^100, which
// specifies a sequence of length 100 having GC content at least 80%.
//
// The motif syntax allows abbreviation using single levels of parentheses and
// exponents.  White space is ignored.  The following are examples of valid 
// motifs:
//
// (GC)^3
// N^100 (GC)^3 N^100
// C^5 G C^4
// (G+C >= 0.8)^100
// (G >= 0.9, G+C = 1)^10
//
// The genome reference sequence is presented to the program by a vecbasevector
// for the bases and a vecbitvector for the ambiguous bases.  Ambiguous bases
// never match.
//
// Composite motifs.  Multiple motifs may be concatenated using "|".  The
// resulting composite motif represents all bases that match one or more of the
// motifs.

#include "Basevector.h"
#include "Bitvector.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "bias/BadMotifCore.h"

int main(int argc, char **argv) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(REF_HEAD,
          "files REF_HEAD.fastb and REF_HEAD.fastamb must exist");
     CommandArgument_String_Doc(BAD, "filename of vecbitvector of bad bases");
     CommandArgument_String_Doc(MOTIF, "motif to be tested");
     CommandArgument_String_OrDefault_Doc(PRINT_INTERVALS, "",
          "print intervals marked by motif to this file");
     CommandArgument_String_OrDefault_Doc(CACHE_DIR, "",
          "to cache vecbitvector of bases marked by motif on genome");
     EndCommandArguments;

     // Load data.

     vecbitvector bad(BAD);
     vecbasevector genome( REF_HEAD + ".fastb" );
     vecbitvector amb( REF_HEAD + ".fastamb" );

     // Parse motif into components.

     vec<String> motifs;
     Tokenize( MOTIF, '|', motifs );

     // Define structure that tracks which bases are in motif.

     vecbitvector in_motif;
     Mimic( genome, in_motif );

     // Check cache.

     String cache_all = CACHE_DIR + "/" + WhiteSpaceFree(MOTIF);
     if ( CACHE_DIR != "" && IsRegularFile(cache_all) )
     {    in_motif.ReadAll(cache_all);
          motifs.clear( );    }

     // Go through the motif components.

     for ( size_t mi = 0; mi < motifs.size( ); mi++ )
     {
          // Check cache.

          String cache = CACHE_DIR + "/" + WhiteSpaceFree( motifs[mi] );
          if ( CACHE_DIR != "" && IsRegularFile(cache) )
          {    vecbitvector im;
               im.ReadAll(cache);
               for ( size_t i = 0; i < in_motif.size( ); i++ )
               {    for ( unsigned int j = 0; j < in_motif[i].size( ); j++ )
                         in_motif[i].Set( j, in_motif[i][j] | im[i][j] );    }
               continue;    }

          // Compute motif.
          VecMatchedLocsVec unused_location_info(genome.size());
          ComputeMotif( motifs[mi], genome, amb, in_motif, unused_location_info );    }

     // Write to cache.

     if ( CACHE_DIR != "" && motifs.solo( ) && !IsRegularFile(cache_all) )
          in_motif.WriteAll(cache_all);

     // Report results.
     
     cout << "marked   " << "\t" << "marked+bad" << "\t" << "(marked+bad)/marked"
          << "\t" << "motif" << endl;
     longlong total = 0, evil = 0;
     for ( size_t i = 0; i < genome.size( ); i++ )
     {    for ( unsigned int j = 0; j < genome[i].size( ); j++ )
          {    if ( in_motif[i][j] )
               {    ++total;
                    if ( bad[i][j] ) ++evil;    }    }    }
     longlong tdigits = 1, t = total, edigits = 1, e = evil;
     while( t >= 10 )
     {    ++tdigits;
          t /= 10;    }
     cout <<      String(std::max( 0L, 9 - tdigits ), ' ') << total << "\t";
     while( e >= 10 )
     {    ++edigits;
          e /= 10;    }
     cout <<      String(std::max( 0L, 9 - edigits ), ' ') << evil << "\t";
     if ( total == 0 ) cout << "  0.0";
     else
     {    double eper = 100.0 * double(evil) / double(total);
          if ( eper == 100.0 ) cout << "100.0";
          else
          {    cout << " " << setiosflags(ios::fixed) << setw(4) << setprecision(1)
                    << eper << resetiosflags(ios::fixed);    }    }
     cout << "%" << "            \t" << MOTIF << endl;

     // Dump intervals.

     if ( PRINT_INTERVALS != "" )
     {    Ofstream( out, PRINT_INTERVALS );
          for ( size_t i = 0; i < in_motif.size( ); i++ )
          {    for ( unsigned int j = 0; j < in_motif[i].size( ); j++ )
               {    if ( !in_motif[i][j] ) continue;
                    unsigned int k = in_motif[i].NextDiff(j);
                    out << i << " " << j << " " << k << "\n";
                    j = k - 1;    }    }    }    }
