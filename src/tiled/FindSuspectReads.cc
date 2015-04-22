// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "STLExtensions.h"
#include "tiled/FindSuspectReads.h"

#include <map>
#include <functional>

struct order_by_first
  : public binary_function< pair<int,int>, pair<int,int>, bool>
{
  bool operator() ( const pair<int,int>& lhs, const pair<int,int>& rhs) const
    { return ( lhs.first < rhs.first ); }
};


void FindSnpCycles( const vec< pair<int,int> > &snps, 
                    const vec<read_pairing> &pairs,
                    const vec<int> &pairs_index,
                    vec<snp_cycle> &vecCycles )
{
    vecCycles.clear();

    int numSnp( snps.size() );
    
    pair< vec< pair<int,int> >::const_iterator,  \
        vec< pair<int,int> >::const_iterator > range;
    for ( range.first = snps.begin(); range.first != snps.end(); range.first = range.second )
    {
        int read1 = range.first->first;
        
        range = equal_range( range.first,
                             snps.end(), 
                             *range.first,
                             order_by_first() );
        
        vec<int> snpsWithRead1( distance( range.first, range.second ) );
        transform( range.first, range.second,
                   snpsWithRead1.begin(),
                   select2nd< pair<int,int> >() );
        
        if ( snpsWithRead1.size() > 1 )
        {
            for ( int j=0; j<(int) snpsWithRead1.size(); ++j )
            {
                int read2 = snpsWithRead1[j];
                
                for ( int k=j+1; k<(int) snpsWithRead1.size(); ++k )
                {
                    int read3 = snpsWithRead1[k];
                    
                    pair< vec< pair<int,int> >::const_iterator,  
                        vec< pair<int,int> >::const_iterator > new_range;	
                    
                    // Check for a snp between reads 2 and 3.
                    pair<int,int> target = make_pair( read2, read3 );
                    new_range = equal_range( snps.begin(),
                                             snps.end(),
                                             target );
                    
                    // If one is found, we have a 3-cycle.
                    if ( new_range.first != new_range.second )
                    {
                        set<int> setReadsInCycle;
                        setReadsInCycle.insert( read1 );
                        setReadsInCycle.insert( read2 );
                        setReadsInCycle.insert( read3 );
                        
                        vecCycles.push_back( snp_cycle( setReadsInCycle ) );
                    }
                    
                    // Look for 5-cycles involving the partners of read2 and read3.
                    int pairForRead2 = pairs_index[ read2 ];
                    int pairForRead3 = pairs_index[ read3 ];
                    
                    if ( pairForRead2 < 0 || pairForRead3 < 0 )
                        continue;
                    
                    int read2Partner = pairs[ pairForRead2 ].Partner( read2 );
                    int read3Partner = pairs[ pairForRead3 ].Partner( read3 );
                    
                    // Check for a snp between read2Partner and read3Partner.
                    target = make_pair( read2Partner, read3Partner );
                    new_range = equal_range( snps.begin(),
                                             snps.end(),
                                             target );
                    
                    // If one is found, we have a 5-cycle.
                    if ( new_range.first != new_range.second )
                    {
                        set<int> setReadsInCycle;
                        setReadsInCycle.insert( read1 );
                        setReadsInCycle.insert( read2 );
                        setReadsInCycle.insert( read3 );
                        setReadsInCycle.insert( read2Partner );
                        setReadsInCycle.insert( read3Partner );
                        
                        vecCycles.push_back( snp_cycle( setReadsInCycle ) );
                    }
                }
            }
        }
    }
    
    sort( vecCycles.begin(), vecCycles.end() );
    vecCycles.erase( unique( vecCycles.begin(), vecCycles.end() ), 
                     vecCycles.end() );
}



vec<int> FindSuspectReads( const vec<snp_cycle> &snpCycles, 
                           ostream *logp )
{
  vec<int> vecSuspects;

  // Create a map from the read id to the cycles it participates in.
  multimap<int,const snp_cycle*> mapReadToCyclePtr;

  for ( vec<snp_cycle>::const_iterator cycleIter = snpCycles.begin();
        cycleIter != snpCycles.end(); ++cycleIter )
  {
      for ( set<int>::const_iterator readIter = cycleIter->mSetReadsInCycle.begin();
            readIter != cycleIter->mSetReadsInCycle.end(); ++readIter )
          mapReadToCyclePtr.insert( make_pair( *readIter, &*cycleIter ) );
  }

  int numCyclesRemaining = snpCycles.size();

  // While there are still cycles to break...
  while ( ! mapReadToCyclePtr.empty() )
  {
      pair<multimap<int,const snp_cycle*>::iterator,
          multimap<int,const snp_cycle*>::iterator> range;

      range.first = mapReadToCyclePtr.begin();

      // Find the read participating in the most cycles.
      int worstOffenderId = range.first->first;
      int worstOffenderSize = 1;

      while ( range.first != mapReadToCyclePtr.end() )
      {
          range = mapReadToCyclePtr.equal_range( range.first->first );

          int numCycles = distance( range.first, range.second );
          if ( numCycles > worstOffenderSize )
          {
              worstOffenderId = range.first->first;
              worstOffenderSize = numCycles;
          }

          range.first = range.second;
      }

      if ( logp )
          *logp << "Read " << worstOffenderId
                << " is in " << worstOffenderSize << " cycles"
                << " (out of " << numCyclesRemaining << ")"
                << endl;

      vecSuspects.push_back( worstOffenderId );

      range = mapReadToCyclePtr.equal_range( worstOffenderId );

      // Find the cycles in which the worst offender participates.
      set<const snp_cycle*> setDeadTriPtrs;
      for ( ; range.first != range.second; ++range.first )
          setDeadTriPtrs.insert( range.first->second );

      // Remove them from consideration.
      multimap<int,const snp_cycle*>::iterator mapIter = mapReadToCyclePtr.begin();
      while ( mapIter != mapReadToCyclePtr.end() )
      {
          if ( setDeadTriPtrs.count( mapIter->second ) )
          {
              multimap<int,const snp_cycle*>::iterator deadIter = mapIter++;
              mapReadToCyclePtr.erase( deadIter );
          }
          else
              ++mapIter;
      }
      
      numCyclesRemaining -= setDeadTriPtrs.size();
  }
  
  return vecSuspects;
}   

