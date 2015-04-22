// Copyright (c) 2003 Whitehead Institute for Biomedical Research
// 

/// This class manages all the information associated with an assembly.
/// \class Assembly
/// An Assembly can be ReadIn(), which simply attaches it to a given
/// set of files.  As much as possible, it defers actually loading in
/// any information until it is specifically requested.
///
/// Mostly, an Assembly serves as the source for Supers, Contigs, and
/// Reads, read-only entities that answer informational questions.
///
/// Currently, an assembly can be loaded and saved, and entities may be
/// added to it. Entities can also be modified or removed.

#ifndef ASSEMBLY_ASSEMBLY
#define ASSEMBLY_ASSEMBLY

#include "Basevector.h"
#include "CompressedSequence.h"
#include "Qualvector.h"
#include "String.h"

#include "assembly/Super.h"
#include "assembly/AssemblyContig.h"
#include "assembly/AssemblyRead.h"

#include "assembly/ContigLocation.h"
#include "assembly/ReadLocationInContig.h"
#include "assembly/ReadPair.h"
#include "assembly/AddBehavior.h"

class SuperDataManager;
class ContigDataManager;
class ReadDataManager;

class ContigLocationManager;
class ReadLocationManager;
class ReadPairManager;

class PrefetchStrategy;

class Assembly 
{
  public:
    /// Create an empty assembly.
    Assembly( bool verbose = true );

    /// Create an assembly and point it at the given directories.
    Assembly( const String &strFullWorkPath, const String &strSubDir,
              bool verbose = true );

    ~Assembly( );

  private:
    /// To prevent the unintentional copying of assemblies, the
    /// assignment operator is private and undefined.
    Assembly & operator= ( const Assembly & );

  public:

    /// ENTITY CREATION ROUTINES

    Read     NewRead    ( const String &name,
                          const basevector &bases,
                          const qualvector &quals,
                          const pair<int,int> trims,
                          const CompressedSequence &untrimmedBases,
			  AddBehavior behavior = ASSERT_IF_DIFFERENT );

    ReadPair NewReadPair( const ReadToken &read1,
                          const ReadToken &read2,
                          const int sep, 
                          const int stdev );
    
    Contig   NewContig  ( const basevector &bases,
                          const qualvector &quals );
    
    Super    NewSuper   ( );

    /// ENTITY COPYING ROUTINES
    // These can be used to copy entities from one assembly to another.

    Read     AddRead    ( const Read   &theRead,       AddBehavior behavior = ASSERT_IF_DIFFERENT );
    ReadPair AddReadPair( const ReadPair &theReadPair, AddBehavior behavior = ASSERT_IF_DIFFERENT );
    Contig   AddContig  ( const Contig &theContig,     AddBehavior behavior = ASSERT_IF_DIFFERENT );
    Super    AddSuper   ( const Super  &theSuper,      AddBehavior behavior = ASSERT_IF_DIFFERENT );


    /// ENTITY PLACEMENT ROUTINES

    void AddContigLocation( const ContigLocation &theNewLocation );
    void AddReadLocation  ( const ReadLocation   &theNewLocation );


    /// ENTITY OBSCURING ROUTINES
  
    /// Pretend the reads in this pair are unpaired.  Has no effect if
    /// the read pair is already ignored.  The read pair is still
    /// accessible via Assembly::GetReadPair(), but the constituent
    /// Reads do not know about it and will not return it via
    /// Read::GetPair().
    void IgnorePair( const ReadPair &thePair );

    /// Undo the above.  Has no effect if the read is not currently ignored.
    void UsePair( const ReadPair &thePair );

    /// Ignore pairs that appear to be from duplicated inserts.  I.e., if
    /// the positions of the untrimmed ends of some pair of reads are
    /// both within duplicate_template_threshold bases of the untrimmed
    /// ends of some other pair, ignore one of the the two pairs.
    void IgnoreDuplicateTemplates( const int duplicateTemplateThreshold,
                                   ostream *pLog = 0 );

    // LOCATION REMOVAL ROUTINES

    /// Remove the given ReadLocation.  If no such location exists,
    /// nothing is done.
    void RemoveReadLocation( const ReadLocation &theLocation );

    /// Remove all the read locations for a given Contig.  The Contig
    /// is still valid, but contains no reads.
    void Clear( Contig theContig );

    /// Remove the given ContigLocation.  If no such location exists,
    /// nothing is done.
    void RemoveContigLocation( const ContigLocation &theLocation );

    /// Remove all the contig locations for a given Super.  The Super
    /// is still valid, but empty.
    void Clear( Super theSuper );


    // ENTITY DESTRUCTION ROUTINES

    /// Destroy all the Supers in the Assembly.  GetNumSupers() will
    /// return zero after this call, and all Supers in this assembly
    /// will be invalid.
    void DestroyAllSupers();

    /// Destroy all the Supers and Contigs in the Assembly.  The effect
    /// of DestroyAllSupers() is produced, and in addition,
    /// GetNumContigs() will return zero after this call, and all
    /// Contigs in this assembly will be invalid.
    void DestroyAllSupersAndContigs();

    /// Remove all the read locations from the given Contig, set the
    /// length of the consensus sequence of the Contig to zero, and
    /// remove all the ContigLocations of the given Contig.  The Contig
    /// still exists, but contains no reads, has no consensus sequence,
    /// and is not placed in any Super.
    void RemoveContig( Contig theContig );

    /// Remove all the contig locations from the given Super, and also
    /// remove the Contigs (via Assembly::RemoveContig()) in that
    /// Super.  The Super still exists, but contains no contigs, and
    /// the Contigs still exist, but contain no reads, have no
    /// consensus sequence, and are not placed in any Super.
    void RemoveSuper( Super theSuper );


    // ENTITY MODIFICATION ROUTINES

    /// Reverse the ContigLocations in a Super.  (The bases and quals
    /// of the Contigs remain the same, just the orientations and
    /// positions of the ContigLocations are changed.)
    void Reverse( const Super  &theSuper );

    /// Reverse the ReadLocations in a Contig.  (The bases and quals
    /// of the Reads remain the same, just the orientations and
    /// positions of the ReadLocations are changed.)
    void Reverse( const Contig &theContig );

    /// Shift all the ContigLocations in a Super by the given amount,
    /// which may be negative.
    void Shift    ( Super theSuper, const int shiftAmount );

    /// Shift all the ContigLocations in a Super such that the leftmost
    /// starts at zero.
    void Normalize( Super theSuper );

    /// Return a new supercontig that contains the contigs from superA,
    /// the specified gap, then the contigs from superB.  Has the effect
    /// of calling Clear() on superA and superB.
    Super Merge( Super superA, int gap, Super superB );
  
    // ACCESSOR METHODS

    int GetNumReads     ( ) const;
    int GetNumReadPairs ( ) const;
    int GetNumContigs   ( ) const;
    int GetNumSupers    ( ) const;

    /// Read returned will not be valid if no such read is found.
    Read   GetRead   ( const String &name ) const;

    Read     GetRead     ( int id ) const; // Will assert if id is out of range.
    ReadPair GetReadPair ( int id ) const; // Will assert if id is out of range.
    Contig   GetContig   ( int id ) const; // Will assert if id is out of range.
    Super    GetSuper    ( int id ) const; // Will assert if id is out of range.

    void GetAllReads    ( vec<Read>     &vecReads )   const;
    void GetAllReadPairs( vec<ReadPair> &vecReads )   const;
    void GetAllContigs  ( vec<Contig>   &vecContigs ) const;
    void GetAllSupers   ( vec<Super>    &vecSupers )  const;


    // PREFETCHING METHODS

    /// The caller of this method is still responsible for the
    /// management of this pointer.
    void SetReadPrefetchStrategy( PrefetchStrategy *pStrategy ) const;
    void SetContigPrefetchStrategy( PrefetchStrategy *pStrategy ) const;

    void PrefetchAllReads() const;
    void PrefetchReadsByContig() const;
    void PrefetchReadsBySuper() const;
    void PrefetchReadsSingly() const;

    void PrefetchAllContigs() const;
    void PrefetchContigsBySuper() const;
    void PrefetchContigsSingly() const;

    /// This is the default behavior.
    void PrefetchAll() const {
      this->PrefetchAllReads();
      this->PrefetchAllContigs();
    }

    void PrefetchByContig() const {
      this->PrefetchReadsByContig();
      this->PrefetchContigsSingly();
    }

    void PrefetchBySuper() const {
      this->PrefetchReadsBySuper();
      this->PrefetchContigsBySuper();
    }

    void PrefetchSingly() const {
      this->PrefetchReadsSingly();
      this->PrefetchContigsSingly();
    }

    // I/O METHODS

    /// Read in the raw data for an assembly from the filesystem.
    /// Will destroy current contents of assembly.
    void ReadRawData( const String &strFullWorkPath );

    /// Read in a full assembly from the filesystem.
    /// Will destroy current contents of assembly.
    void ReadIn( const String &strFullWorkPath, const String &strSubDir );

    /// By default, the Assembly will purge empty contigs (contigs without
    /// bases and without reads) when saving.  This method controls
    /// whether this purging occurs.
    void SetPurgeBeforeWrite( bool purgeBeforeWrite ) 
    {  mPurgeBeforeWrite = purgeBeforeWrite;  }

    /// By default, the Assembly will compute the deviations of the gap
    /// estimates.  This method controls whether this computation occurs.  If
    /// the deviations are not computed, they are set to the absolute value of
    /// the gap (but at least 1), to indicate their arbitrariness.
    void SetComputeDevs( bool computeDevs ) 
    {  mComputeDevs = computeDevs;  }

    // These are still being implemented.  Use with caution.
    void WriteOutSafely( const String &strFullWorkPath, const String &strSubDir )
    {  this->WriteOut( strFullWorkPath, strSubDir, false, false );  }

    void WriteOutWithSubdirOverwrite( const String &strFullWorkPath, const String &strSubDir ) const
    {  this->WriteOut( strFullWorkPath, strSubDir, true, false );  }

    void WriteOutWithFullOverwrite( const String &strFullWorkPath, const String &strSubDir ) const
    {  this->WriteOut( strFullWorkPath, strSubDir, true, true );  }

  private:
    String GetSourcePath( const String &strFullWorkPath ) const;
    String GetSubdirPath( const String &strFullWorkPath,
			  const String &strSubDir ) const;

    String GetContigLocationsFile( const String &strFullWorkPath,
				   const String &strSubDir ) const;
    String GetReadLocationsFile  ( const String &strFullWorkPath,
				   const String &strSubDir ) const;

    String GetReadPairingsFile( const String &strFullWorkPath,
				const String &strSubDir = "" ) const;
    String GetReadNamesFile   ( const String &strFullWorkPath,
				const String &strSubDir = "" ) const;
    String GetReadLengthsFile ( const String &strFullWorkPath,
			        const String &strSubDir = "" ) const;
    String GetReadFastbFile   ( const String &strFullWorkPath,
			        const String &strSubDir = "" ) const;
    String GetReadQualbFile   ( const String &strFullWorkPath,
				const String &strSubDir = "" ) const;
    String GetReadTrimsFile   ( const String &strFullWorkPath,
			        const String &strSubDir = "" ) const;
    String GetReadFastnFile   ( const String &strFullWorkPath,
				const String &strSubDir = "" ) const;
    String GetReadRepetitiveFlagFile ( const String &strFullWorkPath,
                                       const String &strSubDir = "" ) const;

    String GetContigFastbFile( const String &strFullWorkPath,
			       const String &strSubDir ) const;
    String GetContigQualbFile( const String &strFullWorkPath,
			       const String &strSubDir ) const;

    void WriteOut( const String &strFullWorkPath, const String &strSubDir,
                   const bool overwriteContigsAndSupers, const bool overwriteReads ) const;

    SuperDataManager  *mpSuperDataMgr;
    ContigDataManager *mpContigDataMgr;
    ReadDataManager   *mpReadDataMgr;

    ContigLocationManager *mpContigLocationMgr;
    ReadLocationManager   *mpReadLocationMgr;
    ReadPairManager       *mpReadPairMgr;

    bool  mPurgeBeforeWrite;
    bool  mComputeDevs;
    Bool  mVerbose;
};

#endif
