/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "pairwise_aligners/PerfectAligner.h"

#include "pairwise_aligners/MatchList.h"
#include "kmers/SortKmers.h"
#include "STLExtensions.h"
#include "kmers/SupportedKmerShapes.h"

// Abstract superclass to define interface for classes that can filter matches.

struct MatchFilter {
  virtual bool Accept( int pos1, int pos2, int len, int len1, int len2 ) const = 0;
  virtual ~MatchFilter() {}
};


// Concrete subclasses to filter matches.

struct ProperFilter : public MatchFilter {
  virtual bool Accept( int pos1, int pos2, int len, int len1, int len2 ) const
  {
    // Two ends must be reached.
    return ( ( pos1 == 0 || pos2 == 0 ) &&
             ( pos1 + len == len1 || pos2 + len == len2 ) );
  }
};

struct SemiproperFilter : public MatchFilter {
  virtual bool Accept( int pos1, int pos2, int len, int len1, int len2 ) const
  {
    // One end must be reached.
    return ( ( pos1 == 0 || pos2 == 0 ) ||
             ( pos1 + len == len1 || pos2 + len == len2 ) );
  }
};

struct ImproperFilter : public MatchFilter {
  virtual bool Accept( int pos1, int pos2, int len, int len1, int len2 ) const
  { 
    return true;
  }
};


// Abstract superclass for implementations of PerfectAligner.

class PerfectAlignerImp {
 public:
  PerfectAlignerImp() {}
  virtual ~PerfectAlignerImp() {};

  virtual void Align( const vecbasevector& sequences, 
                      vec<alignment_plus>& perfectAligns,
                      const int partition = 0,
                      const int maxKmerFreq = 0 ) = 0;
};


// Templatized implementation of PerfectAligner.

template <int K>
class PerfectAlignerImpTemplate : public PerfectAlignerImp {

 public:
  PerfectAlignerImpTemplate( PerfectAligner::Behavior behavior,
                             ostream* pLog )
    : PerfectAlignerImp(),
      m_pFilter( 0 ),
      m_pLog( pLog )
  {
    switch ( behavior ) {
      case PerfectAligner::findProperOnly:
        m_pFilter = new ProperFilter;
        break;
      case PerfectAligner::findSemiproper:
        m_pFilter = new SemiproperFilter;
        break;
      case PerfectAligner::findImproper:
        m_pFilter = new ImproperFilter;
    }
  }

  ~PerfectAlignerImpTemplate()
  {
    delete m_pFilter;
  }

 private:
  MatchFilter* m_pFilter;
  ostream* m_pLog;

  typedef vec< kmer_record<K,2> > kmervector;

 public:
  /// If partition greater than or equal to zero, only create aligns between
  /// sequences on either side of the partition index.
  void Align( const vecbasevector& sequences, vec<alignment_plus>& perfectAligns,
              const int partition = -1, const int maxKmerFreq = 0 ) 
  {
    bool selfCompare = ( partition < 0 );

    vec<int> seqIds( sequences.size() );
    for ( unsigned int i = 0; i < seqIds.size(); ++i ) 
      seqIds[i] = i;

    MatchList perfectMatchList( sequences.size() );
    kmervector kmerRecords;

    bool verbose = ( sequences.sumSizes() > 10 * 1000 * 1000 );
    
    for ( int pass = 0; pass < 100; ++pass )
    {
      if ( m_pLog ) {
        if ( verbose ) {
          *m_pLog << "Pass " << pass+1 << ":" << endl;
          *m_pLog << "Sorting... " << flush;
        }
        else
          Dot( *m_pLog, pass );
      }

      unsigned int numRecords;
      SortKmers( dummy<100>(), sequences, seqIds, pass, kmerRecords, numRecords );
      
      kmerRecords.resize( numRecords );

      if ( verbose && m_pLog ) 
        *m_pLog << "done." << endl;
      
      if ( verbose && m_pLog )
        *m_pLog << "Processing " << numRecords << " kmers in 100 passes." << endl;

      unsigned int passSize = ( numRecords + 99 ) / 100;
      unsigned int nextDotThres = passSize;
      unsigned int nextDotNum = 0;

      typename kmervector::iterator kmerBegin, kmerEnd;
      for ( kmerBegin = kmerRecords.begin(); kmerBegin != kmerRecords.end(); kmerBegin = kmerEnd )
      {
        kmerEnd = upper_bound( kmerBegin, kmerRecords.end(),
                               *kmerBegin );

        if ( maxKmerFreq != 0 &&
             distance( kmerBegin, kmerEnd ) > maxKmerFreq ) {
          if ( verbose && m_pLog ) 
            while ( longlong(distance( kmerRecords.begin(), kmerEnd )) > 
		    longlong(nextDotThres) ) {
              Dot( *m_pLog, nextDotNum++ );
              nextDotThres += passSize;
            }

          continue;
        }
        
        typename kmervector::iterator kmer1, kmer2;
        for ( kmer1 = kmerBegin; kmer1 != kmerEnd-1; ++kmer1 )
        {
          if ( verbose && m_pLog 
	       && longlong(distance( kmerRecords.begin(), kmer1 )) 
	       >= longlong(nextDotThres) ) {
            Dot( *m_pLog, nextDotNum++ );
            nextDotThres += passSize;
          }
          int id1 = kmer1->GetId();

          for ( kmer2 = kmer1+1; kmer2 != kmerEnd; ++kmer2 )
          {
            int id2 = kmer2->GetId();

            if ( selfCompare || ( (id1<partition) != (id2<partition) ) )
              perfectMatchList.ProcessMatchingKmers( kmer1->GetPos(), kmer2->GetPos(), K,
                                                     id1, id2,
                                                     &sequences[id1], &sequences[id2] );
          }
        }
      }

      if ( verbose && m_pLog ) {
          while ( nextDotNum < 100 )
            Dot( *m_pLog, nextDotNum++ );
        *m_pLog << endl;
      }
    }

    if ( m_pLog && ! verbose )
      *m_pLog << endl;

    if ( m_pLog )
      *m_pLog << "Building alignments... " << flush;
    
    vec<Match> matches;
    
    for ( size_t id = 0; id < sequences.size(); ++id )
    {
      perfectMatchList.GetMatches( id, matches );
      
      for ( unsigned int m = 0; m < matches.size(); ++m )
      {
        // TODO: potentially dangerous truncation of index
        int id1 = id;
        int id2 = matches[m].GetId2();
        bool rc1 = false;
        bool rc2 = matches[m].GetRc();
        int pos1 = matches[m].GetPos1();
        int pos2 = matches[m].GetPos2();
        int len  = matches[m].GetLen();
        int len1 = sequences[id1].size();
        int len2 = sequences[id2].size();

        if ( ! m_pFilter->Accept( pos1, pos2, len, len1, len2 ) )
          continue;

        if ( ! selfCompare )
        {
          if ( id1 >= partition )
          {
            id1 -= partition;
            swap( id1, id2 );
            swap( len1, len2 );
            swap( rc1, rc2 );
            swap( pos1, pos2 );
          }
          else
          {
            ForceAssertGe( id2, partition );
            id2 -= partition;
          }
        }

        avector<int> gaps(1), lens(1);
        gaps(0) = 0;
        lens(0) = len;
        int errs = 0;
        alignment theAlignment( pos1, pos2, errs, gaps, lens );
        if ( rc1 )
        {
          theAlignment.ReverseThis( len1, len2 );
          rc1 = !rc1;
          rc2 = !rc2;
        }

        alignment_plus theAlignmentPlus( id1, id2, len1, len2,
                                         rc2, theAlignment, Float( 0.0 ) );

        perfectAligns.push_back( theAlignmentPlus );

        if ( selfCompare && id1 != id2 )
        {
          theAlignmentPlus.Swap( len1, len2 );
          perfectAligns.push_back( theAlignmentPlus );
        }
      }
    }

    if ( m_pLog )
      *m_pLog << "done." << endl;
  }
};


PerfectAligner::PerfectAligner( int kmerSize, Behavior behavior, ostream* pLog )
  : m_maxKmerFreq( 0 ),
    m_pLog( pLog ),
    m_pImp( 0 )
{
  SetKmerSize( kmerSize );
  SetBehavior( behavior );
}


PerfectAligner::~PerfectAligner()
{
  delete m_pImp;
}


void
PerfectAligner::SetKmerSize( const int kmerSize )
{
  if ( kmerSize != m_kmerSize ) {
    if ( m_pImp != 0 ) {
      delete m_pImp;
      m_pImp = 0;
    }
    m_kmerSize = kmerSize;
  }
}


void
PerfectAligner::SetBehavior( const PerfectAligner::Behavior behavior )
{
  if ( behavior != m_behavior ) {
    if ( m_pImp != 0 ) { 
      delete m_pImp;
      m_pImp = 0;
    }
    m_behavior = behavior;
  }
}


void
PerfectAligner::SetMaxKmerFreq( const int max )
{
  m_maxKmerFreq = max;
}


void
PerfectAligner::Align( const vecbasevector& sequences, 
                       vec<alignment_plus>& perfectAligns,
                       const int partition,
		       const bool append )
{
  if ( ! append )
    perfectAligns.clear( );

  ForceAssertLe( partition, (int)sequences.size() );
  
  if ( partition == 0 || partition == (int)sequences.size() ) return;
  
  if ( ! m_pImp ) {
#define CASE(K) m_pImp = new PerfectAlignerImpTemplate<K>( m_behavior, m_pLog )
      DISPATCH_ON_K( m_kmerSize, CASE );
  }

  m_pImp->Align( sequences, perfectAligns, partition, m_maxKmerFreq );
}
