#include "ReadLocation.h"
#include "SeqInterval.h"

#include <set>
#include <functional>

#ifndef BAGGED_READS_H
#define BAGGED_READS_H

// This class simplifies the management of collections of read locations.

class bagged_reads 
{
  public:
    bagged_reads();
    bagged_reads( const String &strReadLocsFilename );
    bagged_reads( const vec<read_location> &vecReadLocs );

    void  Load( const String &strReadLocsFilename );
    void  CleanUp( );
    void  SortAndSave( const String &strReadLocsFilename );
    void  Copy( vec<read_location> &vecReadLocsCopy ) const  { vecReadLocsCopy = mVecLocs; }

    size_t  size() const  { return mVecLocs.size(); }

    void  reserve( const size_t size )  { mVecLocs.reserve( size ); }

    // We only provide constant references, otherwise we can't keep
    // track of whether our indexing is accurate.
    const read_location &  operator[] ( const unsigned int index ) const
    {
        return mVecLocs[ index ];
    }

    // Returns number of contigs (some of which may be empty).
    int  GetNumberOfContigs( ) const;
    
    // Fills out vecLocIndices with the indices of locations on the
    // given contig in the order that they start on that contig.
    void  GetIndicesForContig( const int contigId, 
                               vec<int> &vecLocIndices ) const;
    
    // Fills out vecLocIndices with the indices of locations
    // intersecting the given interval on the given contig in the
    // order that they start on that contig.
    void  GetIndicesForIntervalOnContig( const int contigId, const seq_interval &theInterval,
                                         vec<int> &vecLocIndices ) const;
    
    // Fills out vecLocIndices with the indices of locations of the
    // given read in arbitrary order.
    void  GetIndicesForRead( const int readId, 
                             vec<int> &vecLocIndices ) const;

    // Adds the given location to the collection.  Returns the index
    // of the added location.
    int  AddLocation( const read_location &newLoc );

    // Removes the location at the given index from the collection.
    void  RemoveLocationAt( const int locIndex );

    // Removes the locations for a given contig from the collection.
    void  RemoveContig( const int contigId );

    // The JoinContigs() method joins two contigs.  The contig id of
    // all the reads in both contigs are set to firstContigId, and the
    // start on contig of all the reads that were in secondContigId
    // are incremented by offsetOfSecondFromFirst.  The length of the
    // contig is adjusted in each read location to extend from the
    // start of the leftmost read to the stop of the rightmost read.
    void  JoinContigs( const int firstContigId,
                       const int secondContigId,
                       const int offsetOfSecondFromFirst );

    // The SplitContig() method splits the specified contig in two.
    // The locations whose indices are specified in locsToKeep remain
    // in the contig. (Indices specified in locsToKeep that refer to
    // locations not currently in the specified contig are ignored.)
    // All other locations are placed in a new contig.  The start on
    // contigs of all reads in both contigs is then normalized such
    // that the leftmost read of each contig starts at 0.  The lengths
    // of the contigs is adjusted in each read location to extend from
    // the start of the leftmost read to the stop of the rightmost
    // read.  The id of the new contig is returned.
    //
    // Note that no check is done to ensure that either resulting
    // contig is contiguous with respect to its coverage by read
    // locations.
    int  SplitContig( const int contigId, 
                      const vec<int>& locsToKeep );

    // The SplitNonContigs() method splits contigs where the read
    // locations indicate that they are not contiguous, i.e. it walks
    // through the reads in a contig until it finds a read that starts
    // after the maximum stop of all the reads walked over thus far.
    // If it finds such a read, if places that read and all the reads
    // following in a new contig and continues on.
    //
    // The number of new contigs created is returned.
    int SplitNonContigs();

    // The RemoveEmptyContigs() method removes contigs with no reads,
    // renumbering the remaining contigs so that they are contiguous,
    // but preserving their order.
    void RemoveEmptyContigs();

  private:
    // Recreates the by-contig and by-read indices into the vector of read locations.
    void  Index() const;
    
    // Adjusts the starts of the read locations on the contig so that
    // the leftmost read starts at 0.  Also adjusts the length of the
    // contig in all the read locations so that it extends exactly
    // from the start of the leftmost read to the stop of the
    // rightmost read.
    void  NormalizeLocs( const int contigId );

    vec<read_location>  mVecLocs;

    // We may want to switch one or both of these to sets of sets/vecs
    // or a multimap if insertion/deletion/resizing performance
    // becomes an issue.

    mutable bool             mbIndexed;
    mutable vec< vec<int> >  mContigToLocIdx;
    mutable vec< vec<int> >  mReadToLocIdx;
};

#endif
