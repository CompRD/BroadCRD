// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

// Finds the minimal extensions (muxes) to the left for all reads.

// An extension is a read that matches the given read perfectly and
// either extends beyond the given read to the left or is coterminal
// with (i.e., extends to the same point as) the given read but is
// shorter.

// The minimal extensions of a given read are the extensions of that
// read that are not an extension of any other extension of that read.

// It is possible for a read to have multiple minimal extensions if
// the genome diverges just to the left of the read.

#include "paths/KmerPathDatabase.h"
#include "paths/OrientedKmerPathId.h"
#include "paths/MuxGraph.h"
#include "STLExtensions.h"
#include "TaskTimer.h"
#include "Histogram.h"


class MuxFinder {

  class LeftExtension;

public:
  MuxFinder( const vecKmerPath& pathsFw,
	     const vecKmerPath& pathsRc,
	     const KmerPathDatabase& pathsDb,
             ostream *p_log = 0 )
    : m_pathsFw( pathsFw ),
      m_pathsRc( pathsRc ),
      m_pathsDb( pathsDb ),
      m_minNumKmers( 1 ),
      m_minOverlap( 1 ),
      m_printHiMuxCount( false ),
      m_equalLengthGaplessReads( false ),
      mp_log( p_log ),
      m_verbose( false )
  { }
  
  void SetMinNumKmers( int num )
  { m_minNumKmers = num; }

  void SetMinOverlap( int over )
  { m_minOverlap = over; }

  void SetVerbose( bool verb )
  { m_verbose = verb; }

  void SetPrintHiMuxCount( bool print )
  { m_printHiMuxCount = print; }

  void SetEqualLengthGaplessReads( bool value )
  { m_equalLengthGaplessReads = value; }

  // This examines the reads and maybe sets m_equalLengthGaplessReads
  void StudyReads( );

  // Process the pathsDb given in the constructor from start to stop,
  // printing out stats every chunkSize entries, putting the results
  // in theMuxGraph.
  longlong FindMuxes( MuxGraph& theMuxGraph,
		      const longlong chunkSize,
		      const longlong start,
		      const longlong stop, 
                      const int partition = 0 );

  // Simplified interface.
  longlong FindMuxes( MuxGraph& theMuxGraph,
                      const int partition );

private:
  void EstablishDominance( LeftExtension& extA, 
			   LeftExtension& extB );

  int EstablishDominance( vec<LeftExtension>& extensions,
			  longlong& dominanceCounter1,
			  TaskTimer& dominanceTimer1, 
			  longlong& dominanceCounter2,
			  TaskTimer& dominanceTimer2 );

private:
  const vecKmerPath& m_pathsFw;
  const vecKmerPath& m_pathsRc;
  const KmerPathDatabase& m_pathsDb;
  int m_minNumKmers;
  int m_minOverlap;
  bool m_printHiMuxCount;
  bool m_equalLengthGaplessReads;
  ostream* mp_log;
  bool m_verbose;




  // helper class: represents an extension to the left of some path
  class LeftExtension {
  public:
    LeftExtension( const OrientedKmerPathId& pathId,
		   const KmerPathLoc& firstMatchingKmer )
      : m_pathId( pathId ),
	m_pathLoc( firstMatchingKmer ),
	m_isDominated( false )
    {}
    
    OrientedKmerPathId GetPathId()  const { return m_pathId; }
    const KmerPath&    GetPath()    const { return m_pathLoc.GetPath(); }
    KmerPathLoc        GetLoc()     const { return m_pathLoc; }
    int                GetSegment() const { return m_pathLoc.GetIndex(); }

    Mux ToMux() const
    {
      int numKmers = DistMin( m_pathLoc.GetPath().Begin(),
			      m_pathLoc );
      return Mux( m_pathId, m_pathLoc.GetIndex(), numKmers );
    }

    Bool IsDominated() const { return m_isDominated; }
    void SetIsDominated( const Bool value ) { m_isDominated = value; }
    
    bool operator< ( const LeftExtension& other ) const 
    {
      return ( m_pathLoc.GetIndex() < other.m_pathLoc.GetIndex() );
    }
    
    friend ostream& operator<< ( ostream& out, const LeftExtension& e )
    {
      out << e.m_pathId << " is " << ( e.m_isDominated ? "submissive" : "dominant  " ) << " ";
      for ( int i = 0; i <= e.m_pathLoc.GetIndex(); ++i )
	out << e.m_pathLoc.GetPath().Segment(i);
      return out;
    }

  private:
    OrientedKmerPathId m_pathId;
    KmerPathLoc m_pathLoc;
    Bool m_isDominated;
  };


  
};
