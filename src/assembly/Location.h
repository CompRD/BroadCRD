///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A Location is the location and orientation of a given entity
// within another.

#ifndef ASSEMBLY_LOCATION
#define ASSEMBLY_LOCATION

#include "assembly/Interval.h"
#include "assembly/Orientation.h"

#include <functional>


using std::binary_function;


template <class ContainerT, class ElementT>
class Location
{
 public:
  Location( ) {}
  Location( const ContainerT &theContainer, const ElementT &theElement, 
            const Interval theInterval, const Orientation theOrientation )
    : m_container( theContainer ), m_element( theElement ), 
      m_interval( theInterval ), m_orientation( theOrientation ) { }
  
  typedef ContainerT ContainerType;
  typedef ElementT ElementType;

  ContainerT  GetContainer( )   const { return m_container; }
  ElementT    GetElement( )     const { return m_element; }
  Interval    GetInterval( )    const { return m_interval; }
  Orientation GetOrientation( ) const { return m_orientation; }
  
  bool IsForward( ) const { return m_orientation == orient_FW; }
  bool IsReverse( ) const { return m_orientation == orient_RC; }
  
  int Begin()  const { return m_interval.Begin(); }
  int End()    const { return m_interval.End(); }
  int Length() const { return m_interval.Length(); }
  int First()  const { return m_interval.First(); }
  int Last()   const { return m_interval.Last(); }
  
  // Output.
  
  friend
  ostream & operator<< ( ostream & out, 
                         const Location<ContainerT,ElementT> &loc )
  {
    return out << loc.m_element << ' '
               << loc.m_orientation << '@'
               << loc.m_container << ':' << loc.Begin() << '-' << loc.End();
  }
  
  // Comparators.
  
  friend 
  bool operator< ( const Location<ContainerT,ElementT> &lhs,
                   const Location<ContainerT,ElementT> &rhs )
  {
    if ( lhs.m_container.GetId() < rhs.m_container.GetId() ) return true;
    if ( lhs.m_container.GetId() > rhs.m_container.GetId() ) return false;

    if ( lhs.m_interval < rhs.m_interval ) return true;
    if ( lhs.m_interval > rhs.m_interval ) return false;

    if ( lhs.m_orientation < rhs.m_orientation ) return true;
    if ( lhs.m_orientation > rhs.m_orientation ) return false;

    if ( lhs.m_element < rhs.m_element ) return true;
    return false;
  }
  
  friend 
  bool operator== ( const Location<ContainerT,ElementT> &lhs,
                    const Location<ContainerT,ElementT> &rhs )
  {
    return ( lhs.m_container == rhs.m_container &&
             lhs.m_element == rhs.m_element &&
             lhs.m_interval == rhs.m_interval &&
             lhs.m_orientation == rhs.m_orientation );
  }
  
  friend 
  bool operator!= ( const Location<ContainerT,ElementT> &lhs,
                    const Location<ContainerT,ElementT> &rhs )
  {
    return ! ( lhs == rhs );
  }
  
  // True if and only if the other read location overlaps this one.
  bool HasOverlapWith( const Location<ContainerT,ElementT> & other ) const
  {
    return( m_interval.HasOverlapWith( other.GetInterval() ) );
  }
  
  // returns 0 if no overlap
  int GetOverlapWith( const Location<ContainerT,ElementT> & other ) const
  {
    return( m_interval.GetOverlapWith( other.GetInterval() ) );
  }
  
  Location<ContainerT,ElementT>
  CopyShiftedBy( const int amount ) const
  {
    return Location<ContainerT,ElementT>( m_container, m_element,
                                          m_interval.CopyShiftedBy( amount ),
                                          m_orientation );
  }
  
  Location<ContainerT,ElementT>
  CopyFlipped() const
  {
    return Location<ContainerT,ElementT>( m_container, m_element,
                                          m_interval,
                                          Flip( m_orientation ) );
  }
  
  protected:
    ContainerT  m_container;
    ElementT    m_element;
    Interval    m_interval;
    Orientation m_orientation;
};

#endif
