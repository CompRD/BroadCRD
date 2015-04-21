/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/** The IndexBuilder assists in the creation of indexed lists of unique labels.
    The indexed list is built up one label at a time, with each new and unique
    label assigned the next index number. Duplicate labels are ignored.
    \class IndexBuilder
*/

#ifndef INDEXBUILDER_H
#define INDEXBUILDER_H

#include <map>

#include "VecUtilities.h"
#include "system/System.h"
#include "system/Types.h"


template <class T> class IndexBuilder {

private:
  map<T, int> m_lmap;
  int m_max_index;

public:

  typedef pair<int, T>  LabelPair;

  IndexBuilder()            : m_lmap() { m_max_index = 0; }


  /// Constructor using name of file containing IndexBuilder object data
  IndexBuilder(string filename) {
    ifstream in;
    OpenIfstream( in, filename );
    in >> (*this);
  }

  /// Adds a label, returns the index number for the label
  int AddLabel(T label) {
    int index = m_lmap[label];
    if (!index) {
      m_lmap[label] = ++m_max_index;
      index = m_max_index;
    }
    return index - 1;
  }

  /// Returns the index number associated with the label (-1 if not found)
  int GetIndex(T label) {
    return m_lmap[label] - 1;
  }

  /// Returns number of Labels in Object
  int size() const {
    return m_lmap.size();
  }

  /// Generates sorted indexed table of Labels
  void GetIndexedLabels(vec<T> &labelTable) const {

    // clear and reserve space for the freq table
    labelTable.clear();
    labelTable.reserve(m_lmap.size());

    vec<int> index;
    index.reserve(m_lmap.size());

    typename  map<T, int>::const_iterator pos;
    for (pos = m_lmap.begin(); pos != m_lmap.end(); ++pos) {
      labelTable.push_back(pos->first);
      index.push_back(pos->second);
    }

    // Sort labels by index number
    SortSync(index, labelTable);
  }    

};

#endif
