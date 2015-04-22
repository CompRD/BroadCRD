/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef TOTALKMERSET_H
#define TOTALKMERSET_H

#include "Basevector.h"

/// Class: TotalKmerSet
///
/// A class to track what part of the total 4^k space of kmers is
/// actually occupied by a given set of basevectors.  Can take the base sequence
/// of a kmer and say whether that kmer is present.
///
/// *NOTE*: This class used to be called KmerSpace, but was renamed to TotalKmerSet
/// to avoid conflict with paths/KmerSpace.
///
/// See also <KmerIdSet>.
///
template <int K>
class TotalKmerSet {

 private:
  bool Exists( const longlong kmerId ) {
    return ( binary_search( m_exists.begin(), m_exists.begin() + m_sortedEnd, kmerId ) );
  }

 public:
  TotalKmerSet()
    : m_sortedEnd( 0 )
  {}

  TotalKmerSet( const String& filename ) {
    this->Read( filename );
  }

  TotalKmerSet( const vecbasevector& vecBases ) {
    this->Build( vecBases );
  }

  void Build( const vecbasevector& vecBases )
  {
    m_sortedEnd = 0;

    longlong mask = 4;
    for ( int k = 1; k < K; ++k )
      mask <<= 2;
    --mask;
    mask >>= 2;

    longlong numKmers = 0;
    for ( size_t i = 0; i < vecBases.size(); ++i )
      numKmers += vecBases[i].size() - (K-1);
    numKmers *= 2;

    cout << "Filling kmer space in ";
    longlong kmersPerPass = 1000;
    while ( numKmers / kmersPerPass >= 70 )
      kmersPerPass *= 2;
    cout << numKmers / kmersPerPass + 1 << " passes:" << endl;

    longlong kmerCount = 0;
    for ( size_t i = 0; i < vecBases.size(); ++i ) {
      if ( vecBases[i].size() < K ) continue;
      basevector bases = vecBases[i];
      for ( int rc = 0; rc < 2; ++rc ) {
        if ( rc == 1 )
          bases.ReverseComplement();

        longlong kmerId = bases[0];
        for ( int base = 1; base < K; ++base ) {
          kmerId *= 4;
          kmerId += bases[base];
        }
        if ( ! this->Exists( kmerId ) )
          m_exists.push_back( kmerId );
        if ( ++kmerCount % kmersPerPass == 0 )
          Dot( cout, kmerCount / kmersPerPass );

        for ( unsigned int j = K; j < bases.size(); ++j ) {
          kmerId &= mask;
          kmerId *= 4;
          kmerId += bases[j];
          if ( ! this->Exists( kmerId ) ) {
            m_exists.push_back( kmerId );
            if ( m_sortedEnd < m_exists.size() * 9 / 10 )
              this->Pack();
          }
          if ( ++kmerCount % kmersPerPass == 0 )
            Dot( cout, kmerCount / kmersPerPass );
        }
      }
    }
    this->Pack();

    cout << endl;
  }

 public:
  bool Exists( const basevector& kmer ) {
    longlong kmerId = kmer[0];
    for ( int base = 1; base < K; ++base ) {
      kmerId *= 4;
      kmerId += kmer[base];
    }
    return this->Exists( kmerId );
  }

  void Write( const String& filename ) {
    this->Pack();

    ofstream out( filename.c_str() );
    
    size_t size = m_exists.size();
    out.write( (char*)&size, sizeof(size) );
    
    for ( vector<longlong>::const_iterator i = m_exists.begin(); i != m_exists.end(); ++i )
      out.write( (char*)&(*i), sizeof(longlong) );
  }

  void Read( const String& filename ) {
    ifstream in( filename.c_str() );
    
    size_t size;
    in.read( (char*)&size, sizeof(size) );
    
    m_exists.resize( size );
    for ( vector<longlong>::iterator i = m_exists.begin(); i != m_exists.end(); ++i )
      in.read( (char*)&(*i), sizeof(longlong) );
    
    m_sortedEnd = size;
  }

 private:
  void Pack() {
    if ( m_sortedEnd == m_exists.size() ) 
      return;
    vector<longlong>::iterator middle = m_exists.begin() + m_sortedEnd;
    sort( middle, m_exists.end() );
    m_exists.erase( unique( middle, m_exists.end() ),
                    m_exists.end() );
    inplace_merge( m_exists.begin(), middle, m_exists.end() );
    m_sortedEnd = m_exists.size();
  }

  vector<longlong> m_exists;
  size_t m_sortedEnd;
};  // class TotalKmerSet

#endif
// #ifndef TOTALKMERSET_H

