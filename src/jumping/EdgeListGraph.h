/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <map>
#include <set>
#include <iostream> 

// Graph structure based on multimap.  Edge (a,b) is entered in edge
// lists for vertices a and b.  
struct EdgeListGraph : public
std::map<int, std::set<int, std::greater<int> >, std::greater<int> >
{
  typedef std::set<int, std::greater<int> > Set;
  void add(int a, int b)
  {
    (*this)[a].insert(b);
    (*this)[b].insert(a);
  }
  // Remove vertex a from graph by removing all edges incident on it.
  void remove(int a)
  {
    // The set bs is all the b's for edges (a,b) in graph
    Set &bs = (*this)[a];
    // Now erase the (b,a) edge for all b's
    Set::iterator it=bs.begin(), last=bs.end();
    for ( ; it!= last; ++it) {
      //std::cout << "Erasing " << a << " from edges for " << *it << std::endl;
      Set & bEdges = (*this)[*it];
      bEdges.erase(a);
      if (bEdges.empty()) {
	//std::cout << "Erasing " << *it << " from graph." << std::endl;
	erase(*it);
      }
    }
    // Now we can eliminate the (a,b) edges for all a
    erase(a);
    //std::cout << "Erasing " << a << " from graph." << std::endl;
  }

  bool getEdge(int &a, int &b)
  {
    if (empty())
      return false;
    a = begin()->first;
    b = *(begin()->second.begin());
    return true;
  }
};

std::ostream & operator<< (std::ostream &out, EdgeListGraph & G);
