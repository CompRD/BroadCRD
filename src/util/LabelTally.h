/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/** A LabelTally contains an associative array of label objects
    along with their corresponding frequency.
    If a new unique label is added to a LabelTally object then
    a new entry is created with frequency value 1. If a label
    is added that already exists in the LabelTally object then
    the assocated frequency is increased by 1. The frequency of
    a particular label can be retrieved, returning 0 if the label
    does not exist in the LabelTally object.
    Expects Label class to be able to read and write itself.
    \class LabelTally
*/

#ifndef LABELTALLY_H
#define LABELTALLY_H

#include <map>

#include "system/System.h"
#include "system/Types.h"


template <class T> class LabelTally {

private:
  map<T, int> m_lmap;

public:

  typedef pair<int, T>  LabelPair;

  LabelTally()            : m_lmap() {}

  /// Constructor using map of Labels and associated frequencies
  LabelTally(const map<T, int>& m) : m_lmap( m ) {}

  /// Constructor using name of file containing LabelTally object data
  LabelTally(string filename) {
    ifstream in;
    OpenIfstream( in, filename );
    in >> (*this);
  }

  /// Adds a label, increases label frequency by 1
  void AddLabel(T label) {
    ++m_lmap[label];
  }

  /// Adds a label, increases label frequency by tally
  void AddLabel(T label, int tally) {
    m_lmap[label] += tally;
  }

  /// Adds a label, sets label frequency to tally
  void SetTally(T label, int tally) {
    m_lmap[label] = tally;
  }

  /// Returns frequency associated with Label, 0 if label does not exist
  int GetTally(T label) {
    return m_lmap[label];
  }

  /// Returns number of Labels in Object
  int size() const {
    return m_lmap.size();
  }

  /// Returns the sum of Label frequencies
  int Sum() const {
    int tallySum = 0;
  
    typename  map<T, int>::const_iterator pos;
    for (pos = m_lmap.begin(); pos != m_lmap.end(); ++pos) {
      tallySum += pos->second;
    }
    return tallySum;
  }

  /// Generates label frequency table, sorted in order starting most probable.
  void GetFreqTable(vec<pair<int, T> > &freqTable) const {

    // clear and reserve space for the freq table
    freqTable.clear();
    freqTable.reserve(m_lmap.size());

    typename  map<T, int>::const_iterator pos;
    for (pos = m_lmap.begin(); pos != m_lmap.end(); ++pos) {
      freqTable.push_back(LabelPair(pos->second, pos->first));
    }
    sort(freqTable.begin(), freqTable.end(), greater<LabelPair>());

  }    


  /// Generates cumulative frequency tables, sorted in order starting most probable.
  void GenerateCumulativeTables(vec<int> &cumFreqTable, vec<T> &labelTable) const {  
    
    // Create list of labels sorted by frequency (most probable first)
    vec<LabelPair> labelFreqTable;
    GetFreqTable(labelFreqTable);
    
    unsigned int tableSize = labelFreqTable.size();
    
    // clear and reserve space for the cumulative label tables
    cumFreqTable.clear();
    labelTable.clear();
    cumFreqTable.reserve(tableSize);
    labelTable.reserve(tableSize);
    
    // Build cumulative frequency tables
    unsigned int cumFreq = 0;
    for(unsigned int i = 0; i < tableSize; ++i) {
      cumFreq += labelFreqTable[i].first;
      cumFreqTable.push_back(cumFreq);
      labelTable.push_back(labelFreqTable[i].second);
    }
  }


  /// Writes LabelTally object to stream
  friend ostream& operator<< ( ostream& o, const LabelTally<T> &t) {
    typename  map<T, int>::const_iterator pos;
    o << t.size() << "\n";
    for (pos = (t.m_lmap).begin(); pos != (t.m_lmap).end(); ++pos) {
      o << pos->first << " " <<  pos->second << "\n";
    }
    o << flush;
    return o;
  }


  /// Reads LabelTally object from stream
  friend istream& operator>> ( istream& i,  LabelTally<T>& t ) {
    int length;
    i >> length;
    T label;
    int tally;
    for ( int j = 0; j < length; j++ ) {
      i >> label >> tally;
      t.addLabel(label, tally);
    }
    return i;
  }
 

};

#endif
