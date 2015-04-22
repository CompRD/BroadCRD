/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <math.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "paths/KmerPath.h"
#include "paths/UnipathCoverageCore.h"
#include "random/Bernoulli.h"


void CopyNumber( int n, int Ilen, double rlen, longlong allreads, 
     longlong gsize, int K, PdfEntryVec& copyno_prob,
     double thresh, double error_rate, double* q_ptr )
{
  ForceAssertGe( n, 0 );
  vec<double> log_nfact( 1, 0 );
     int lnfacts = log_nfact.size( );
     if ( n + 1 > lnfacts )
     {    log_nfact.resize( n + 1 );
          for ( int i = lnfacts; i <= n; i++ )
               log_nfact[i] = log_nfact[i-1] + log( double(i) );    }
     double reads_per_base = double(allreads)/double(gsize);
     double q = reads_per_base * double( Ilen + rlen - 1 );
     if( q_ptr ) *q_ptr=q;
     double logplast = 0.0, sump = 0.0;
     double logthresh = log(thresh);
     double logpmax = -numeric_limits<double>::max();
     static PdfEntryVec cprob;
     cprob.clear( ), copyno_prob.clear( );

     double error_prior = 1 - pow( 1-error_rate, K );
     // naive probability that the a K-mer contains an error;
     // used as a prior against the C=0 hypothesis.

     // First model C=0 (ie error kmer):
     if( error_rate > 0 ) {
          double c = error_rate/4;
          double logp = -q * c + double( n >= 1 ? n-1 : 0 ) * log( q * c ) - log_nfact[n >= 1 ? n-1 : 0];
	  // use n-1 instead of n, since clearly it does appear once
          double p = error_prior * exp(logp);
          cprob.push_back( pdf_entry( 0, logp ) );     }
     // Now model any positive copy number; prior is all are equally likely.
     for ( int C = 1; ; C++ )
     {    double c = C;
          double logp = -q * c + double(n) * log( q * c ) - log_nfact[n];
          if ( (C>1 && logp < logplast) && logp <= logpmax + logthresh ) break;
          cprob.push_back( pdf_entry( C, logp ) );
          logplast = logp;
          logpmax = Max( logpmax, logp );
     }
     // Let prob(C=0) into calculation of pmax, if needed; harmless if not
     logpmax = Max( logpmax, cprob[0].second );

     for ( PdfEntryVec::size_type i = 0; i < cprob.size( ); i++ )
          if ( cprob[i].second >= logpmax + logthresh ) copyno_prob.push_back( cprob[i] );
     // Approximate normalization in log space, using logpmax.
     for ( PdfEntryVec::size_type i = 0; i < copyno_prob.size( ); i++ )
	 copyno_prob[i].second = exp(copyno_prob[i].second - logpmax);
     for ( PdfEntryVec::size_type i = 0; i < copyno_prob.size( ); i++ )
          sump += copyno_prob[i].second;
     for ( PdfEntryVec::size_type i = 0; i < copyno_prob.size( ); i++ )
          copyno_prob[i].second /= sump;    }

void UnipathCoverageCore( 
     // inputs:
     const int K, const vecKmerPath& paths, const vecKmerPath& paths_rc, 
     const vec<tagged_rpint>& pathsdb, const vecKmerPath& unipaths, 
     const vec<tagged_rpint>& unipathsdb, const vec<int>& read_lengths, 
     // output:
     VecPdfEntryVec& p,
     // optional args:
     const int UNIPATH_TO_TRACE, const double THRESH, const double ERROR_RATE, 
     const unsigned int USE_THIS_GENOME_SIZE, const int PLOIDY,
     const vec<double>* p_unipath_bias )

{

     // Compute some numbers.

     longlong readlength_total = 0;
     for ( int j = 1; j < read_lengths.isize( ); j++ )
          readlength_total += read_lengths[j];
     double average_readlength = readlength_total / float(read_lengths.size());
     int nreads = paths.size( ), nuni = unipaths.size( );
     ForceAssertGt( nreads, 0 );
     // Some reads may be empty; probably want to ignore those.
     longlong nonempty_reads=0;
     for ( size_t pi = 0; pi < paths.size( ); pi++ )
          if ( !( paths[pi].IsEmpty( ) ) ) nonempty_reads++;
     ForceAssertGt( nonempty_reads, 0 );

     // Compute coverage of unipaths by reads.  Note that if there are several
     // copies of the unipath, very close together, we tend to underestimate the
     // coverage.  I'm not sure how to compensate for this.

     vec<longlong> coverage_fw( nuni, 0 ), coverage_rc( nuni, 0 );
     vec<longlong> coverage( nuni, 0 );
     vec<longlong> coverage_starts_only( nuni, 0 );
     for ( int pass = 1; pass <= 2; pass++ )
     {    const vecKmerPath& v = ( pass == 1 ? paths : paths_rc );
          for ( size_t i = 0; i < v.size( ); i++ )
          {    static vec< pair<int,int> > uid_offset;
               uid_offset.clear( );
               vec<longlong> locs;
               Contains( unipathsdb, v[i].Start(0), locs );
               if ( locs.empty( ) ) continue;
               ForceAssertEq( locs.size(), 1u );
               coverage_starts_only[ unipathsdb[locs.front()].PathId() ]++;
               for ( int j = 0; j < v[i].NSegments( ); j++ )
               {    const KmerPathInterval& I = v[i].Segment(j);
                    static vec<longlong> locs;
                    Contains( unipathsdb, I, locs );
                    for ( int u = 0; u < locs.isize( ); u++ )
                    {    const tagged_rpint& t = unipathsdb[ locs[u] ];
                         const KmerPath& x = unipaths[ t.PathId( ) ];
                         int xp = t.PathPos( );
                         int pos1 = t.Start( ) - v[i].Start(j);
                         for ( int r = 0; r < j; r++ )
                              pos1 += v[i].Length(r);
                         int pos2 = t.Start( ) - x.Start(xp);
                         for ( int r = 0; r < xp; r++ )
                              pos2 += x.Length(r);
                         int offset = pos2 - pos1;
                         uid_offset.push_back( 
                              make_pair( t.PathId( ), offset ) );    }    }
               UniqueSort(uid_offset);
               for ( int j = 0; j < uid_offset.isize( ); j++ )
               {    if ( uid_offset[j].first == UNIPATH_TO_TRACE )
                    cout << "read " << i << ( pass == 1 ? "fw" : "rc" ) 
                         << " overlaps unipath " << uid_offset[j].first 
                         << " with offset " << uid_offset[j].second << "\n";    }
               for ( int j = 0; j < uid_offset.isize( ); j++ )
               {    if ( pass == 1 ) ++coverage_fw[ uid_offset[j].first ];
                    else ++coverage_rc[ uid_offset[j].first ];    }    }    }

     // Phase 2 of coverage computation.  Try to adjust for the case where we're
     // near the end of a chromosome.  Note that this won't work (for example) if
     // we have exactly two copies of a unipath, in opposite orientation, both near
     // the end.  Note also that a more sophisticated version of this could 
     // explicitly take into account pairing, as close to the end it should be the
     // case that long inserts "reach off the end", but short ones don't.

     for ( int i = 0; i < nuni; i++ )
     {    int c1 = Max( coverage_fw[i], coverage_rc[i] );
          int c2 = Min( coverage_fw[i], coverage_rc[i] );
          if ( i == UNIPATH_TO_TRACE ) {
            PRINT2( coverage_fw[i], coverage_rc[i] );
          }
          if ( c1 + c2 == 0 ) continue;
          double p = pow( 0.5, c1 + c2 - 1 ) * PartialBernoulliSum( c1 + c2, c2 );
          if ( p < 0.001 ) coverage[i] = 2 * Max( coverage_fw[i], coverage_rc[i] );
          else coverage[i] = coverage_fw[i] + coverage_rc[i];     
          if ( i == UNIPATH_TO_TRACE ) {
            PRINT( coverage[i] );
          }
     }

     // Estimate genome size by looking at the top 20 unipaths by length.  It might 
     // be better to use more, but require (by using pairing) that they are not 
     // proximate to each other.

     vec<int> used_for_genome_size_est;

     vec<int> ulen(nuni);
     for ( int i = 0; i < nuni; i++ )
          ulen[i] = unipaths[i].KmerCount( );
     int top = 20, high_length = 0;
     {    vec<int> unipath_lengths = ulen;
          ReverseSort(unipath_lengths);
          if ( nuni >= top ) high_length = unipath_lengths[ top - 1 ];    }
     longlong total_bases = 0, total_reads = 0;
     for ( int i = 0; i < nuni; i++ )
     {    if ( ulen[i] < high_length ) continue;
          used_for_genome_size_est.push_back(i);
          total_bases += ulen[i];
          total_reads += coverage_starts_only[i];    }
     longlong genome_size_estimate = PLOIDY * int( round( 
          double( nonempty_reads ) * double(total_bases) / double(total_reads) ) );

     PRINT(genome_size_estimate);

     if( USE_THIS_GENOME_SIZE ) {
       cout << "Command-line overruling genome size estimate: using "
	    << USE_THIS_GENOME_SIZE << endl;
       genome_size_estimate = USE_THIS_GENOME_SIZE;
     }

     // Compensate for gc bias if necessary.
     if ( p_unipath_bias )
       for ( int i = 0; i < nuni; i++ )
         coverage[i] = (int) round( (double) coverage[i] * 1.0/(*p_unipath_bias)[i] );

     // Estimate number of copies.

     vec<double> q(nuni);  // used for debugging only

     p.clear( );
     p.resize(nuni);
     longlong count = 0;
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 ) p.Reserve( count, nuni );
          for ( int i = 0; i < nuni; i++ )
          {    int n = coverage[i];
               int Ilen = unipaths[i].KmerCount( );
	       double rlen = average_readlength - K + 1;
	       PdfEntryVec copyno_prob;
               CopyNumber( n, Ilen, rlen, nonempty_reads, genome_size_estimate, 
			   K, copyno_prob, THRESH, ERROR_RATE, &q[i] );
               if ( pass == 1 ) count += copyno_prob.size( );
               else
               {    for ( PdfEntryVec::size_type j = 0; j < copyno_prob.size( ); j++ )
                         p[i].push_back( copyno_prob[j] );    }    }    }
     if( UNIPATH_TO_TRACE >= 0 ) {
       int u = UNIPATH_TO_TRACE;
       cout << "Unipath " << u << " had actual coverage " << coverage[u] << ", and "
            << "expected coverage " << q[u] << endl;
       cout << "The resulting pdf was: ";
       for(PdfEntryVec::size_type j=0; j<p[u].size(); j++)
         cout << p[u][j].first << "(" << p[u][j].second << ") ";
       cout << endl;
     }
     if( USE_THIS_GENOME_SIZE ) {
       cout << "Original genome size est was based on:" << endl;
       for( vec<int>::iterator i = used_for_genome_size_est.begin();
	    i != used_for_genome_size_est.end(); i++ ) {
	 cout << "unipath " << *i << ", " << ulen[*i] << " kmers, touches "
	      << coverage[*i] << " reads; q=" << q[*i] << ", pdf ";
	 for(PdfEntryVec::size_type j=0; j<p[*i].size(); j++)
	   cout << p[*i][j].first << "(" << p[*i][j].second << ") ";
	 cout << endl;
       }
     }

}
