// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef ASSEMBLY_DISPLAY_CHANNELER_H
#define ASSEMBLY_DISPLAY_CHANNELER_H

// A Channeler takes a set of objects that responds to the methods
// Begin() and End(), e.g. Intervals, ReadLocations, etc, and finds a
// partitioning of the objects into numbered bins (channels) such that
// no two objects in a bin are within minSpace units of one another.
// It tries to make the minimum number of bins possible, and to put
// objects in the lowest numbered bin possible.  Channels are numbered
// from zero up.

#include <map>
#include <vector>

template <class T>
class Channeler
{
 public:
  // The vector "objects" must be sorted by Begin().
  Channeler( const vec<T> &objects,
             const int minSpace = 1 );

  // Returns the channel for the given object, or -1 if the object was
  // not in the set of objects channeled in the constructor.
  int GetChannelFor( const T &object ) const;
  
  // Returns the channel for the given id (corresponding to an input object),
  // or 1 if the id is out of range.
  int GetChannelFor( const int id ) const;
  
  // Returns the number of channels used, i.e. one more than the
  // largest answer GetChannelFor will ever give.
  int GetNumChannels( ) const { return m_numChannels; }

 private:
  map<T,int> m_channelObjMap;
  vec<int>   m_channelIdMap;
  int m_numChannels;
};

template <class T>
inline
Channeler<T>::Channeler( const vec<T> &objects,
                         const int minSpace )
{
  if ( objects.empty() )
    return;

  vector<int> channelEnds;
  
  for ( unsigned int objIdx = 0; objIdx < objects.size(); ++objIdx )
  {
    const T &object = objects[ objIdx ];

    unsigned int channel = 0;
    while ( channel < channelEnds.size() &&
            channelEnds[ channel ] > object.Begin() )
      ++channel;

    if ( channel == channelEnds.size() )
      channelEnds.push_back( object.End() + minSpace );
    else
      channelEnds[ channel ] = object.End() + minSpace;

    m_channelObjMap.insert( make_pair( object, channel ) );
    m_channelIdMap.push_back( channel );
  }

  m_numChannels = channelEnds.size();
}

template <class T>
inline
int
Channeler<T>::GetChannelFor( const T &object ) const
{
  pair<typename map<T,int>::const_iterator,typename map<T,int>::const_iterator> range;
  
  range = m_channelObjMap.equal_range( object );

  if ( range.first == range.second )
    return -1;
  else
    return range.first->second;
}

template <class T>
inline
int
Channeler<T>::GetChannelFor( const int id ) const
{
  if ( id < 0 || id >= m_channelIdMap.isize( ) )
    return -1;
  return m_channelIdMap[id];
}
      
#endif
