// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

// StartPointBias.  For a given 454 run, investigate the bias in start points of
// the reads:
// - Ignore regions of reference which support multiply placed reads.
// - Look only at reads placed on forward strand of reference.
// - Suppose some of these reads are chosen at random so that they are spaced every
//   5 bases.
// - Compute the frequency with which gaps of size 20 bases or 50 bases are
//   encountered, as compared to if the reads fell at random.
//
// For now, we require that there is a single target contig.
//
// Cautionary points:
// - If the DNA source is different from the reference, there could be indels.
// - Theoretically, gaps could arise because of lack of sensitivity in the alignment
//   of reads to the reference.  I've seen no evidence of this.
// - If there is too much duplication within the reference, much or all of the
//   reference could be ignored, leading to results that could be misleading.
// - The code will get stuck an infinite loop if there are not enough reads.

// Three test cases, the last of which exhibits several instances where high
// GC content results in gaps:
//
// StartPointBias
// ALIGNS=/wga/scr2/454work/newbacs/bac1a/610399.1.TCAG.trimA.longer30.orig.qltout
// REF=/wga/454a/xfer/refseq/BAC_1A_AC005865.fasta
//
// StartPointBias
// ALIGNS=/wga/scr2/454work/newbacs/bac2a/610383.1.TCAG.trimA.longer30.orig.qltout
// REF=/wga/454a/xfer/refseq/BAC_2A_AC018698.fasta
//
// StartPointBias
// ALIGNS=/wga/scr2/454work/newbacs/bac5b/610318.2.TCAG.trimA.longer30.orig.qltout
// REF=/wga/454a/xfer/refseq/BAC_5B_AC040977.fasta

#include "FastIfstream.h"
#include "FetchReads.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "random/Random.h"
#include "dna/Bases.h"
#include <sstream>

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(ALIGNS);
     CommandArgument_String(REF);
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     EndCommandArguments;

     // Read in reference.

     vecbasevector ref;
     FetchReads( ref, 0, REF );

     // Get number of reads and number of alignments and number of bases in
     // reference.  Exclude alignments of less than 75% of a read.

     int nreads = 0, naligns = 0, nbases = 0;
     String line;
     {    fast_ifstream in(ALIGNS);
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( !line.Contains( "QUERY", 0 ) ) continue;
               static look_align la;
               la.ReadParseable(line);
               ForceAssertEq( la.target_id, 0 );
               if ( float( la.a.extent1( ) ) / float( la.query_length ) < 0.75 )
                    continue;
               ++naligns;
               nbases = Max( nbases, la.a.Pos2( ) );
               nreads = Max( nreads, la.query_id + 1 );    }    }

     // Read in alignments.

     vec<look_align_plus> aligns(naligns);
     vec< vec<int> > aligns_index(nreads);
     int aligns_count = 0;
     fast_ifstream in(ALIGNS);
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( !line.Contains( "QUERY", 0 ) ) continue;
          aligns[aligns_count].ReadParseable(line);
          const look_align& la = aligns[aligns_count];
          if ( float( la.a.extent1( ) ) / float( la.query_length ) < 0.75 ) continue;
          int id = la.query_id;
          aligns_index[id].push_back(aligns_count);
          ++aligns_count;    }

     // Identify regions which yield ambiguously placed reads.  Exclude those
     // of length < 1000.  Exclude first and last 200 bases.

     vec<ho_interval> mult_cov, clean;
     vec< pair<ho_interval, int> > all_cov;
     for ( int i = 0; i < nreads; i++ )
     {    if ( aligns_index[i].size( ) > 1 )
          {    for ( int j = 0; j < aligns_index[i].isize( ); j++ )
               {    const look_align& la = aligns[ aligns_index[i][j] ];
                    int start = Max( 0, la.a.pos2( ) - 1 );
                    int stop = Min( nbases, la.a.Pos2( ) + 1 );
                    mult_cov.push_back( ho_interval( start, stop ) );    }    }    }
     CondenseIntervals( nbases, mult_cov, all_cov );
     for ( int i = 0; i < all_cov.isize( ); i++ )
     {    if ( all_cov[i].second == 0 && all_cov[i].first.Length( ) >= 1000 )
          {    ho_interval h = all_cov[i].first;
               h.Set( Max( 200, h.Start( ) ), Min( nbases - 200, h.Stop( ) ) );
               clean.push_back(h);    }    }
     int nclean = Sum(clean);
     cout << nclean << " clean bases in " << clean.size( ) << " intervals\n";
     if (VERBOSE)
     {    for ( int j = 0; j < clean.isize( ); j++ )
               cout << clean[j] << "\n";    }

     // Find start points of forward-placed reads that lie in a clean region.

     vec<int> starts;
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    if ( aligns[i].rc1 ) continue;
          int start = aligns[i].a.pos2( );
          if ( Member( clean, start ) ) starts.push_back(start);    }
     Sort(starts);
     cout << "\nforward-placed read starts: " << starts.size( ) << "\n";

     // Generate nclean/5 random reads, and compute the separations between
     // consecutive start points.  Compute number of separations which are at
     // least 20.  Repeat 10 times to assess variability.  Note that the methodology
     // here is slightly different than what we do for the real reads (which is to
     // look at clean intervals), so that some artefactual discrepancy is introduced
     // at this point.

     vec<int> rstarts, rseps;
     vec<float> rtotal1s, rtotal20s, rtotal50s;
     int samples = 100;
     for ( int j = 0; j < samples; j++ )
     {    rstarts.clear( ), rseps.clear( );
          for ( int i = 0; i < nclean/5; i++ )
               rstarts.push_back( randomx( ) % nclean );
          Sort(rstarts);
          int last_rstart = 0, rtotal1 = 0, rtotal20 = 0, rtotal50 = 0;
          for ( int i = 0; i < rstarts.isize( ); i++ )
          {    rseps.push_back( rstarts[i] - last_rstart );
               last_rstart = rstarts[i];    }
          Sort(rseps);
          for ( int i = 0; i < rseps.isize( ); i++ )
          {    if ( rseps[i] >= 1 ) ++rtotal1;
               if ( rseps[i] >= 20 ) ++rtotal20;
               if ( rseps[i] >= 50 ) ++rtotal50;    }
          rtotal1s.push_back(rtotal1);
          rtotal20s.push_back(rtotal20);
          rtotal50s.push_back(rtotal50);    }
     NormalDistribution drtotal1 = SafeMeanStdev(rtotal1s);
     int rtotal1 = int(round(drtotal1.mu_));
     int rtotal1_dev = int(round(drtotal1.sigma_));
     NormalDistribution drtotal20 = SafeMeanStdev(rtotal20s);
     int rtotal20 = int(round(drtotal20.mu_));
     int rtotal20_dev = int(round(drtotal20.sigma_));
     float rtotal20f = 1000.0 * float(rtotal20) / float(nclean);
     float rtotal20f_dev = 1000.0 * float(rtotal20_dev) / float(nclean);
     NormalDistribution drtotal50 = SafeMeanStdev(rtotal50s);
     int rtotal50 = int(round(drtotal50.mu_));
     int rtotal50_dev = int(round(drtotal50.sigma_));
     float rtotal50f = 1000.0 * float(rtotal50) / float(nclean);
     float rtotal50f_dev = 1000.0 * float(rtotal50_dev) / float(nclean);

     // Pick random sample which has one read every 5 bases.  To avoid problems
     // with duplication in the 454 process, we force the number of reads with
     // distinct start points to equal the predicted number (rather than choosing
     // the total number of reads).

     vec<Bool> use( starts.size( ), False ), used_bases( nbases, False );
     int used = 0, used1 = 0;
     while( used1 < rtotal1 )
     {    int r = randomx( ) % starts.size( );
          if ( use[r] ) continue;
          use[r] = True;
          ++used;
          if ( !used_bases[ starts[r] ] ) ++used1;
          used_bases[ starts[r] ] = True;    }
     cout << "using " << used << " forward-placed reads\n";

     // Find separations between consecutive reads.

     vec<int> seps;
     int iclean = 0;
     int last_start = clean[0].Start( );
     int min_high_sep = 50;
     vec<ho_interval> high_seps;
     for ( int i = 0; i < use.isize( ); i++ )
     {    if ( !use[i] ) continue;
          int start = starts[i];
          while( start > clean[iclean].Stop( ) )
          {    int sep = clean[iclean].Stop( ) - last_start;
               seps.push_back(sep);
               if ( sep >= min_high_sep )
               {    high_seps.push_back( ho_interval(
                         last_start, clean[iclean].Stop( ) ) );    }
               ++iclean;
               last_start = clean[iclean].Start( );    }
          int sep = start - last_start;
          seps.push_back(sep);
          if ( sep >= min_high_sep )
               high_seps.push_back( ho_interval( last_start, start ) );
          last_start = start;    }
     Sort(seps);

     // Compare actual to predicted separations.

     int total1 = 0;
     for ( int i = 0; i < seps.isize( ); i++ )
          if ( seps[i] >= 1 ) ++total1;
     if (VERBOSE)
     {    cout << "\nrtotal1: " << rtotal1 << " +/- " << rtotal1_dev << "\n";
          PRINT(total1);    }
     int total20 = 0;
     for ( int i = 0; i < seps.isize( ); i++ )
          if ( seps[i] >= 20 ) ++total20;
     float total20f = 1000.0 * float(total20) / float(nclean);
     cout << "\npredict sep >= 20: " << setprecision(3)
          << rtotal20f << " +/- " << setprecision(2) << rtotal20f_dev << " per kb; ";
     cout << "see sep >= 20: " << setprecision(3) << total20f << " per kb\n\n";
     int total50 = 0;
     for ( int i = 0; i < seps.isize( ); i++ )
          if ( seps[i] >= 50 ) ++total50;
     float total50f = 1000.0 * float(total50) / float(nclean);
     cout << "predict sep >= 50: " << setprecision(3)
          << rtotal50f << " per kb; ";
     cout << "see sep >= 50: " << setprecision(3) << total50f << " per kb\n";
     int max_sep = Max( Max(seps), Max(rseps) );
     vec<int> sep_count( max_sep + 1, 0 );
     vec<int> rsep_count( max_sep + 1, 0 );
     for ( int i = 0; i < seps.isize( ); i++ )
          ++sep_count[ seps[i] ];
     for ( int i = 0; i < rseps.isize( ); i++ )
          ++rsep_count[ rseps[i] ];
     if (VERBOSE)
     {    for ( int i = 0; i <= max_sep; i++ )
          {    if ( sep_count[i] == 0 && rsep_count[i] == 0 ) continue;
               cout << "sep = " << i << ", see " << sep_count[i]
                    << ", predict " << rsep_count[i] << "\n";    }    }

     // Show sequences appearing at high separations.  Show GC content of each
     // sequence, and the rank of sequences of that GC content amongst all
     // reference subsequences of the same size, ordered by GC content.

     cout << "\nsequences appearing at separations >= " << min_high_sep << ":\n";
     for ( int i = 0; i < high_seps.isize( ); i++ )
     {    static basevector b;
          b.SetToSubOf( ref[0], high_seps[i].Start( ), high_seps[i].Length( ) );
          unsigned int GC = b.GcBases();
          int nb = b.size( );

          // Find fraction of sequences of the given length having GC content as
          // high as the given one.

          int total_seq = 0, gc_seq = 0;
          for ( int j = 0; j < clean.isize( ); j++ )
          {    ho_interval h = clean[j];
               for ( int u = h.Start( ); u <= h.Stop( ) - nb; u++ )
               {    ++total_seq;
                    unsigned int gc = ref[0].GcBases(u,u+nb);
                    if ( gc >= GC ) ++gc_seq;    }    }

          ostringstream out;
          out << high_seps[i] << " GC = " << (100 * GC) / b.size( ) << "% "
               << " (" << PERCENT_RATIO( 3, gc_seq, total_seq ) << ")";
          b.Print( cout, out.str( ) );    }    }
