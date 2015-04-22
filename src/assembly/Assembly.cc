/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "AnnotatedContig.h"
#include "Superb.h"

#include "Vec.h"

#include "TaskTimer.h"

#include "assembly/Assembly.h"

#include "assembly/SuperDataManager.h"
#include "assembly/ContigDataManager.h"
#include "assembly/ReadDataManager.h"

#include "assembly/ContigLocationManager.h"
#include "assembly/ReadLocationManager.h"
#include "assembly/ReadPairManager.h"

#include "assembly/IdManager.h"

#include "assembly/NameManager.h"
#include "assembly/ReadSequenceManager.h"

#include "assembly/ContigSequenceManager.h"

#include "assembly/PrefetchStrategy.h"

#include "assembly/OffsetGraph.h"

Assembly::Assembly( bool verbose )
  : mPurgeBeforeWrite( true ),
    mComputeDevs( true ),
    mVerbose( verbose )
{
    mpContigLocationMgr = new ContigLocationManager( mpSuperDataMgr, mpContigDataMgr );
    mpReadLocationMgr   = new ReadLocationManager( mpContigDataMgr, mpReadDataMgr );
    mpReadPairMgr       = new ReadPairManager( mpReadDataMgr );

    mpSuperDataMgr  = new SuperDataManager( mpContigLocationMgr );
    mpContigDataMgr = new ContigDataManager( mpContigLocationMgr, mpReadLocationMgr );
    mpReadDataMgr   = new ReadDataManager( mpReadLocationMgr, mpReadPairMgr );
}

Assembly::Assembly( const String &strFullWorkPath, const String &strSubDir,
                    bool verbose )
    : mpSuperDataMgr( 0 ),
      mpContigDataMgr( 0 ),
      mpReadDataMgr( 0 ),
      mpContigLocationMgr( 0 ),
      mpReadLocationMgr( 0 ),
      mpReadPairMgr( 0 ),
      mPurgeBeforeWrite( true ),
      mComputeDevs( true ),
      mVerbose( verbose )
{
    this->ReadIn( strFullWorkPath, strSubDir );
}

Assembly::~Assembly( )
{
    delete mpContigLocationMgr;
    delete mpReadLocationMgr;
    delete mpReadPairMgr;
    
    delete mpSuperDataMgr;
    delete mpContigDataMgr;
    delete mpReadDataMgr;
}


///////////////////////////
// ENTITY CREATION ROUTINES
///////////////////////////


Read 
Assembly::NewRead( const String &name,
                   const basevector &bases,
                   const qualvector &quals,
                   const pair<int,int> trims,
                   const CompressedSequence &untrimmedBases,
		   AddBehavior behavior )
{
    return mpReadDataMgr->NewRead( name, bases, quals, trims, untrimmedBases, behavior );
}


ReadPair
Assembly::NewReadPair( const ReadToken &read1,
                       const ReadToken &read2, 
                       const int sep, const int stdev )
{
    // Ensure that the reads in this pair are from this assembly.
    ForceAssert( read1 == this->GetRead( read1.GetId() ) );
    ForceAssert( read2 == this->GetRead( read2.GetId() ) );

    // cout << "Pairing R" << read1.GetId() << " to R" << read2.GetId() << endl;

    return mpReadPairMgr->Add( read1, read2, sep, stdev );
}


Contig
Assembly::NewContig( const basevector &bases,
                     const qualvector &quals )
{
    return mpContigDataMgr->NewContig( bases, quals );
}


Super
Assembly::NewSuper( )
{
    return mpSuperDataMgr->NewSuper( );
}


Read
Assembly::AddRead( const Read &theRead, AddBehavior behavior )
{
  // If we're adding a read that's already in this assembly, we're done.
  Read existingRead = this->GetRead( theRead.GetName() );

  if ( existingRead.IsValid() && behavior != REPLACE_EXISTING )
  {
    if ( behavior == ASSERT_IF_DIFFERENT ) 
    {
      ForceAssert( theRead.GetBases() == existingRead.GetBases() );
      ForceAssert( theRead.GetQuals() == existingRead.GetQuals() );
      ForceAssert( theRead.GetTrims() == existingRead.GetTrims() );
      ForceAssert( theRead.GetUntrimmedBases() == existingRead.GetUntrimmedBases() );
    }    

    return existingRead;
  }

  Read newRead = this->NewRead( theRead.GetName(),
                                theRead.GetBases(),
                                theRead.GetQuals(),
                                theRead.GetTrims(),
                                theRead.GetUntrimmedBases(),
				behavior);
  
  return newRead;
}


ReadPair
Assembly::AddReadPair( const ReadPair &newPair, AddBehavior behavior )
{
  // If we're adding a pair that's already in this assembly, we're done.

  // Get the read in this assembly that corresponds (by name) to the first read
  // in the assembly the pair is in.
  Read thisRead1 = this->GetRead( Read( newPair.GetRead1() ).GetName() );

  // Get the read in this assembly that corresponds (by name) to the second
  // read in the assembly the pair is in.
  Read thisRead2 = this->GetRead( Read( newPair.GetRead2() ).GetName() );
  
  // Get the pair in this assembly for that read.
  ReadPair existingPair = mpReadPairMgr->GetByRead( thisRead1.GetId() );

  if ( existingPair.IsValid() && behavior != REPLACE_EXISTING )
  {
    if ( behavior == ASSERT_IF_DIFFERENT ) 
    {
      // Make sure the pairing information is the same.
      ForceAssert( thisRead2 == existingPair.GetOtherRead( thisRead1 ) );
      ForceAssert( newPair.GetExpectedSep() == existingPair.GetExpectedSep() );
      ForceAssert( newPair.GetExpectedStdDev() == existingPair.GetExpectedStdDev() );
    }

    return existingPair;
  }

  return mpReadPairMgr->Add( thisRead1, thisRead2, 
                             newPair.GetExpectedSep(), 
                             newPair.GetExpectedStdDev() );
}
  
  
Contig
Assembly::AddContig( const Contig &theContig, AddBehavior behavior )
{
    Contig myContig = this->NewContig( theContig.GetBases(),
                                       theContig.GetQuals() );

    // cout << "Copying C" << theContig.GetId() << " to C" << myContig.GetId() << endl;

    vec<ReadLocation> vecReadLocations;
    theContig.GetReadLocations( vecReadLocations );

    for ( vec<ReadLocation>::iterator locIter = vecReadLocations.begin();
          locIter != vecReadLocations.end(); ++locIter )
    {
        Read myRead = this->AddRead( locIter->GetRead(), behavior );

        ReadLocation myLoc( myContig, myRead, 
                            locIter->GetInterval(), locIter->GetOrientation() );

        this->AddReadLocation( myLoc );
    }

    return myContig;
}


Super
Assembly::AddSuper( const Super &theSuper, AddBehavior behavior )
{
    Super mySuper = this->NewSuper();

    // cout << "Copying S" << theSuper.GetId() << " to S" << mySuper.GetId() << endl;

    vec<ContigLocation> vecContigLocations;
    theSuper.GetContigLocations( vecContigLocations );

    for ( vec<ContigLocation>::iterator locIter = vecContigLocations.begin();
          locIter != vecContigLocations.end(); ++locIter )
    {
        Contig myContig = this->AddContig( locIter->GetContig(), behavior );

        ContigLocation myLoc( mySuper, myContig, 
                              locIter->GetInterval(), locIter->GetOrientation() );

        this->AddContigLocation( myLoc );
    }

    return mySuper;
}


////////////////////////////
// ENTITY PLACEMENT ROUTINES
////////////////////////////


void
Assembly::AddContigLocation( const ContigLocation &theNewLocation )
{
    // Ensure that the entities in this location are from this assembly.
    SuperToken theSuper = theNewLocation.GetSuper();
    ForceAssert( theSuper == this->GetSuper( theSuper.GetId() ) );

    ContigToken theContig = theNewLocation.GetContig();
    ForceAssert( theContig == this->GetContig( theContig.GetId() ) );

    return mpContigLocationMgr->Add( theNewLocation );
}


void
Assembly::AddReadLocation( const ReadLocation &theNewLocation )
{
    // Ensure that the entities in this location are from this assembly.
    ContigToken theContig = theNewLocation.GetContig();
    ForceAssert( theContig == this->GetContig( theContig.GetId() ) );

    ReadToken theRead = theNewLocation.GetRead();
    ForceAssert( theRead == this->GetRead( theRead.GetId() ) );

    // cout << "Putting R" << theRead.GetId() << " in C" << theContig.GetId() << endl;

    return mpReadLocationMgr->Add( theNewLocation );
}


////////////////////////////
// ENTITY OBSCURING ROUTINES
////////////////////////////
  

void
Assembly::IgnorePair( const ReadPair &thePair )
{
  ForceAssert( thePair == ReadPair( mpReadPairMgr, thePair.GetId() ) );

  mpReadPairMgr->IgnorePair( thePair.GetId() );
}

void
Assembly::UsePair( const ReadPair &thePair )
{
  ForceAssert( thePair == ReadPair( mpReadPairMgr, thePair.GetId() ) );

  mpReadPairMgr->UsePair( thePair.GetId() );
}


struct order_links_by_contig_orient_start_on_contig : public binary_function<Link,Link,bool>
{
  bool operator() ( const Link &lhs, const Link &rhs ) const 
  {
    ReadLocation lhsFirst = lhs.GetReadLocations().first;
    ReadLocation rhsFirst = rhs.GetReadLocations().first;

    if ( lhsFirst.GetContig() < rhsFirst.GetContig() ) return true;
    if ( lhsFirst.GetContig() > rhsFirst.GetContig() ) return false;
    if ( lhsFirst.GetOrientation() < rhsFirst.GetOrientation() ) return true;
    if ( lhsFirst.GetOrientation() > rhsFirst.GetOrientation() ) return false;
    if ( lhsFirst.Begin() < rhsFirst.Begin() ) return true;
    return false;
  }
};

void
Assembly::IgnoreDuplicateTemplates( const int duplicateTemplateThreshold,
                                    ostream *pLog )
{
  if ( pLog )
    *pLog << "Erasing duplicate templates within " << duplicateTemplateThreshold 
          << " bases of another template." << endl;

  vec<Link> contigLinks;
    
  for ( int contigId = 0; contigId < this->GetNumContigs(); ++contigId )
  {
    Contig theContig = this->GetContig( contigId );
    theContig.GetLinks( contigLinks );
    
    if ( contigLinks.empty() )
      continue;

    sort( contigLinks.begin(), contigLinks.end(),
          order_links_by_contig_orient_start_on_contig() );
    
    vec<Link>::iterator prevLink, thisLink;
    prevLink = thisLink = contigLinks.begin();
    ++thisLink;

    pair<ReadLocation,ReadLocation> prevLocs, thisLocs;

    while ( thisLink != contigLinks.end() )
    {
      prevLocs = prevLink->GetReadLocations();
      thisLocs = thisLink->GetReadLocations();
      
      if ( prevLocs.first.GetContig() != thisLocs.first.GetContig() ||
           prevLocs.first.GetOrientation() != thisLocs.first.GetOrientation() ||
           prevLocs.second.GetContig() != thisLocs.second.GetContig() ||
           prevLocs.second.GetOrientation() != thisLocs.second.GetOrientation() )
      {
        prevLink = thisLink++;
        continue;
      }

      if ( thisLocs.first.GetContig() != theContig ) 
      {
        prevLink = thisLink++;
        continue;
      }

      // We find the first endpoint of the templates represented by prevLink and thisLink.
      Read prevFirstRead( prevLocs.first.GetRead() );
      int prevFirstLeftTrim = prevFirstRead.GetLeftTrim();
      int prevFirstTemplateEnd = ( prevLocs.first.IsForward() ?
                                   prevLocs.first.Begin() - prevFirstLeftTrim :
                                   prevLocs.first.End()   + prevFirstLeftTrim );
      

      Read thisFirstRead( thisLocs.first.GetRead() );
      int thisFirstLeftTrim = thisFirstRead.GetLeftTrim();
      int thisFirstTemplateEnd = ( thisLocs.first.IsForward() ?
                                   thisLocs.first.Begin() - thisFirstLeftTrim :
                                   thisLocs.first.End()   + thisFirstLeftTrim );
      
      if ( abs( prevFirstTemplateEnd - thisFirstTemplateEnd ) > duplicateTemplateThreshold )
      {
        prevLink = thisLink++;
        continue;
      }


      // We find the other endpoint of the templates represented by prevLink and thisLink.
      Read prevSecondRead = prevLocs.second.GetRead();
      int prevSecondLeftTrim = prevSecondRead.GetLeftTrim();
      int prevSecondTemplateEnd = ( prevLocs.second.IsForward() ?
                                    prevLocs.second.Begin() - prevSecondLeftTrim :
                                    prevLocs.second.End()   + prevSecondLeftTrim );
      
      Read thisSecondRead = thisLocs.second.GetRead();
      int thisSecondLeftTrim = thisSecondRead.GetLeftTrim();
      int thisSecondTemplateEnd = ( thisLocs.second.IsForward() ?
                                    thisLocs.second.Begin() - thisSecondLeftTrim :
                                    thisLocs.second.End()   + thisSecondLeftTrim );
      
      // If they're too far apart, we keep both links.
      if ( abs( prevSecondTemplateEnd - thisSecondTemplateEnd ) > duplicateTemplateThreshold )
      {
        prevLink = thisLink++;
        continue;
      }

      this->IgnorePair( thisLink->GetPair() );
      
      if ( pLog )
        *pLog << "\n" 
              << thisLocs.first << "\t" << thisLocs.second << "\t"
              << "(template ends: " 
              << thisLocs.first.GetContig() << ":" << thisFirstTemplateEnd << ", "
              << thisLocs.second.GetContig() << ":" << thisSecondTemplateEnd << ")\n"
              << "appears to duplicate:\n"
              << prevLocs.first << "\t" << prevLocs.second << "\t"
              << "(template ends: " 
              << prevLocs.first.GetContig() << ":" << prevFirstTemplateEnd << ", "
              << prevLocs.second.GetContig() << ":" << prevSecondTemplateEnd << ")\n";

      // Leave prevLink where it was.
      // Increment thisLink.
      ++thisLink;
    }
  }
}


////////////////////////////
// LOCATION REMOVAL ROUTINES
////////////////////////////


void
Assembly::RemoveReadLocation( const ReadLocation &theLocation )
{
    mpReadLocationMgr->Remove( theLocation );
}


void
Assembly::Clear( Contig theContig )
{
    ForceAssert( theContig == this->GetContig( theContig.GetId() ) );

    mpReadLocationMgr->ClearContig( theContig.GetId() );
}


void
Assembly::RemoveContigLocation( const ContigLocation &theLocation )
{
    mpContigLocationMgr->Remove( theLocation );
}


void
Assembly::Clear( Super theSuper )
{
    ForceAssert( theSuper == this->GetSuper( theSuper.GetId() ) );

    mpContigLocationMgr->ClearSuper( theSuper.GetId() );
}


//////////////////////////////
// ENTITY DESTRUCTION ROUTINES
//////////////////////////////


void
Assembly::DestroyAllSupers( )
{
    if ( mVerbose )
      cout << "Resetting ContigLocationManager." << endl;
    
    delete mpContigLocationMgr;
    mpContigLocationMgr = new ContigLocationManager( mpSuperDataMgr, mpContigDataMgr );

    mpContigDataMgr->SetLocationManager( mpContigLocationMgr );

    if ( mVerbose )
      cout << "Resetting SuperDataManager." << endl;
    
    delete mpSuperDataMgr;
    mpSuperDataMgr  = new SuperDataManager( mpContigLocationMgr );
}


void
Assembly::DestroyAllSupersAndContigs( )
{
    this->DestroyAllSupers();

    if ( mVerbose )
      cout << "Resetting ReadLocationManager." << endl;
    
    delete mpReadLocationMgr;
    mpReadLocationMgr = new ReadLocationManager( mpContigDataMgr, mpReadDataMgr );

    mpReadDataMgr->SetLocationManager( mpReadLocationMgr );

    PrefetchReadsByContigStrategy prefetchStrategy( mpReadLocationMgr, this->GetNumReads() );
    mpReadDataMgr->SetPrefetchStrategy( &prefetchStrategy );

    if ( mVerbose )
      cout << "Resetting ContigDataManager." << endl;

    delete mpContigDataMgr;
    mpContigDataMgr  = new ContigDataManager( mpContigLocationMgr,
                                              mpReadLocationMgr );
}



void
Assembly::RemoveContig( Contig theContig )
{
    Clear( theContig );

    vec<ContigLocation> contigLocs;
    theContig.GetSelfLocations( contigLocs );

    for ( vec<ContigLocation>::iterator contigLocIter = contigLocs.begin();
          contigLocIter != contigLocs.end();
          ++contigLocIter )
      this->RemoveContigLocation( *contigLocIter );

    mpContigDataMgr->SetBases( theContig.GetId(), basevector(0) );
    mpContigDataMgr->SetQuals( theContig.GetId(), qualvector(0) );

    ForceAssertEq( mpContigDataMgr->GetLength( theContig.GetId() ), 0 );
}


void
Assembly::RemoveSuper( Super theSuper )
{
    vec<ContigLocation> contigLocs;
    theSuper.GetContigLocations( contigLocs );

    for ( vec<ContigLocation>::iterator contigLocIter = contigLocs.begin();
          contigLocIter != contigLocs.end();
          ++contigLocIter )
      this->RemoveContig( contigLocIter->GetContig() );
}


///////////////////////////////
// ENTITY MODIFICATION ROUTINES
///////////////////////////////


void
Assembly::Reverse( const Contig &theContig )
{
    ForceAssert( theContig == this->GetContig( theContig.GetId() ) );

    mpContigDataMgr->Reverse( theContig.GetId() );
}



void
Assembly::Reverse( const Super &theSuper )
{
    ForceAssert( theSuper == this->GetSuper( theSuper.GetId() ) );

    mpSuperDataMgr->Reverse( theSuper.GetId() );
}


void
Assembly::Shift( Super theSuper, const int shiftAmount )
{
    ForceAssert( theSuper == this->GetSuper( theSuper.GetId() ) );

    mpContigLocationMgr->Shift( theSuper.GetId(), shiftAmount );
}


void
Assembly::Normalize( Super theSuper )
{
    ForceAssert( theSuper == this->GetSuper( theSuper.GetId() ) );

    mpContigLocationMgr->Normalize( theSuper.GetId() );
}


Super
Assembly::Merge( Super superA, int gap, Super superB )
{
  this->Normalize( superA );
  this->Normalize( superB );
  
  this->Shift( superB, superA.GetLength() + gap );

  Super newSuper = this->NewSuper();

  vec<ContigLocation> locsA, locsB;
  superA.GetContigLocations( locsA );
  superB.GetContigLocations( locsB );
  
  vec<ContigLocation> locsToCopy( locsA );
  locsToCopy.append( locsB );

  for ( vec<ContigLocation>::iterator locIter = locsToCopy.begin();
        locIter != locsToCopy.end(); ++locIter )
  {
    this->AddContigLocation( ContigLocation( newSuper,
                                             locIter->GetContig(),
                                             locIter->GetInterval(),
                                             locIter->GetOrientation() ) );
    this->RemoveContigLocation( *locIter );
  }

  return newSuper;
}


///////////////////
// ACCESSOR METHODS
///////////////////


int
Assembly::GetNumReads( ) const
{
    return mpReadDataMgr->GetSize();
}

int
Assembly::GetNumReadPairs( ) const
{
    return mpReadPairMgr->GetSize();
}

int
Assembly::GetNumContigs( ) const
{
    return mpContigDataMgr->GetSize();
}

int
Assembly::GetNumSupers( ) const
{
    return mpSuperDataMgr->GetSize();
}


Read
Assembly::GetRead( const int id ) const
{
    return mpReadDataMgr->GetRead( id );
}

Read
Assembly::GetRead( const String &name ) const
{
    return mpReadDataMgr->GetReadNamed( name );
}

ReadPair
Assembly::GetReadPair( const int id ) const
{
    return mpReadPairMgr->Get( id );
}

Contig
Assembly::GetContig( const int id ) const
{
    return mpContigDataMgr->GetContig( id );
}

Super
Assembly::GetSuper( const int id ) const
{
    return mpSuperDataMgr->GetSuper( id );
}

void
Assembly::GetAllReads( vec<Read> &vecReads ) const
{
    vecReads.clear();
    vecReads.resize( this->GetNumReads() );
    
    for ( unsigned int readId = 0; readId < vecReads.size(); ++readId )
        vecReads[ readId ] = Read( mpReadDataMgr, readId );
}

void
Assembly::GetAllReadPairs( vec<ReadPair> &vecReadPairs ) const
{
    vecReadPairs.clear();
    vecReadPairs.resize( this->GetNumReadPairs() );
    
    for ( unsigned int pairId = 0; pairId < vecReadPairs.size(); ++pairId )
        vecReadPairs[ pairId ] = ReadPair( mpReadPairMgr, pairId );
}

void
Assembly::GetAllContigs( vec<Contig> &vecContigs ) const
{
    vecContigs.clear();
    vecContigs.resize( this->GetNumContigs() );
    
    for ( unsigned int contigId = 0; contigId < vecContigs.size(); ++contigId )
        vecContigs[ contigId ] = Contig( mpContigDataMgr, contigId );
}

void
Assembly::GetAllSupers( vec<Super> &vecSupers ) const
{
    vecSupers.clear();
    vecSupers.resize( this->GetNumSupers() );
    
    for ( unsigned int superId = 0; superId < vecSupers.size(); ++superId )
        vecSupers[ superId ] = Super( mpSuperDataMgr, superId );
}


//////////////////////
// PREFETCHING METHODS
//////////////////////


void 
Assembly::SetReadPrefetchStrategy( PrefetchStrategy *pStrategy ) const
{
    mpReadDataMgr->SetPrefetchStrategy( pStrategy );
}


void 
Assembly::SetContigPrefetchStrategy( PrefetchStrategy *pStrategy ) const
{
    mpContigDataMgr->SetPrefetchStrategy( pStrategy );
}


void
Assembly::PrefetchAllReads() const 
{
  PrefetchAllStrategy strat;
  this->SetReadPrefetchStrategy( &strat );
}

void
Assembly::PrefetchReadsByContig() const
{
  PrefetchReadsByContigStrategy strat( mpReadLocationMgr );
  this->SetReadPrefetchStrategy( &strat );
}

void
Assembly::PrefetchReadsBySuper() const
{
  PrefetchReadsBySuperStrategy strat( mpReadLocationMgr );
  this->SetReadPrefetchStrategy( &strat );
}

void
Assembly::PrefetchReadsSingly() const 
{
  NullPrefetchStrategy strat;
  this->SetReadPrefetchStrategy( &strat );
}

void
Assembly::PrefetchAllContigs() const 
{
  PrefetchAllStrategy strat;
  this->SetContigPrefetchStrategy( &strat );
}

void
Assembly::PrefetchContigsBySuper() const
{
  PrefetchContigsBySuperStrategy strat( mpContigLocationMgr );
  this->SetContigPrefetchStrategy( &strat );
}

void
Assembly::PrefetchContigsSingly() const
{
  NullPrefetchStrategy strat;
  this->SetContigPrefetchStrategy( &strat );
}

//////////////
// I/O METHODS
//////////////


String
Assembly::GetSourcePath( const String &strFullWorkPath ) const
{
    String strPreFindSeedsPath = strFullWorkPath + "/preFindSeeds";

    return ( IsDirectory( strPreFindSeedsPath ) 
             ? strPreFindSeedsPath
             : strFullWorkPath );
}


String
Assembly::GetSubdirPath( const String &strFullWorkPath, const String &strSubDir ) const
{
    return strFullWorkPath + "/" + strSubDir;
}


String
Assembly::GetReadLocationsFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return GetSubdirPath( strFullWorkPath, strSubDir ) + "/mergedcontigs_orig.locs";
}


String
Assembly::GetContigLocationsFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return GetSubdirPath( strFullWorkPath, strSubDir ) + "/mergedcontigs.answer";
}


String
Assembly::GetReadPairingsFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return GetSourcePath( strFullWorkPath ) + "/reads.pairto";
}


String
Assembly::GetReadNamesFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return GetSourcePath( strFullWorkPath ) + "/reads.ids";
}


String
Assembly::GetReadLengthsFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return GetSourcePath( strFullWorkPath ) + "/reads.lengths";
}


String
Assembly::GetReadFastbFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return GetSourcePath( strFullWorkPath ) + "/reads.fastb";
}


String
Assembly::GetReadQualbFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return GetSourcePath( strFullWorkPath ) + "/reads.qualb";
}


String
Assembly::GetReadTrimsFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return strFullWorkPath + "/reads.trim_lr";
}


String
Assembly::GetReadFastnFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return strFullWorkPath + "/reads_orig.fastn";
}


String
Assembly::GetReadRepetitiveFlagFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return strFullWorkPath + "/reads.is_repetitive";
}


String
Assembly::GetContigFastbFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return GetSubdirPath( strFullWorkPath, strSubDir ) + "/mergedcontigs.fastb";
}


String
Assembly::GetContigQualbFile( const String &strFullWorkPath, const String &strSubDir ) const
{
    return GetSubdirPath( strFullWorkPath, strSubDir ) + "/mergedcontigs.qualb";
}


void
Assembly::ReadRawData( const String &strFullWorkPath )
{
    if ( ! IsDirectory( strFullWorkPath ) )
    {
      cout << strFullWorkPath << " is not a directory." << endl;
      TracebackThisProcess();
    }
  
    delete mpContigLocationMgr;
    delete mpReadLocationMgr;
    delete mpReadPairMgr;
    
    delete mpSuperDataMgr;
    delete mpContigDataMgr;
    delete mpReadDataMgr;

    if ( mVerbose )
        cout << "Setting up ContigLocationManager." << endl;

    {
        mpContigLocationMgr = new ContigLocationManager( mpSuperDataMgr, mpContigDataMgr );
    }


    if ( mVerbose )
        cout << "Setting up ReadLocationManager." << endl;

    {
        mpReadLocationMgr = new ReadLocationManager( mpContigDataMgr,
						     mpReadDataMgr );
    }


    if ( mVerbose )
        cout << "Setting up ReadPairManager." << endl;

    {
        String strReadPairFile( GetReadPairingsFile( strFullWorkPath ) );
        
        mpReadPairMgr = new ReadPairManager( mpReadDataMgr,
                                             strReadPairFile );
    }


    if ( mVerbose )
        cout << "Setting up ReadDataManager." << endl;

    {
        IdManager *pReadIdMgr = 
            new IdManager;
        
        String strReadIdsFile( GetReadNamesFile( strFullWorkPath ) );

        NameManager *pReadNameMgr =
            new NameManager( strReadIdsFile );
        
        String strReadLengthsFile( GetReadLengthsFile( strFullWorkPath ) );
        String strReadFastbFile  ( GetReadFastbFile  ( strFullWorkPath ) );
        String strReadQualbFile  ( GetReadQualbFile  ( strFullWorkPath ) );
        String strReadTrimsFile  ( GetReadTrimsFile  ( strFullWorkPath ) );
        String strReadFastnFile  ( GetReadFastnFile  ( strFullWorkPath ) );
 
        ReadSequenceManager *pReadSequenceMgr = 
            new ReadSequenceManager( strReadLengthsFile,
                                     strReadFastbFile,
                                     strReadQualbFile,
                                     strReadTrimsFile,
                                     strReadFastnFile );

        int numReads = MastervecFileObjectCount( strReadFastbFile );

        String strRepetitiveFlagFile( GetReadRepetitiveFlagFile( strFullWorkPath ) );
        if ( ! IsRegularFile( strRepetitiveFlagFile ) )
          strRepetitiveFlagFile.erase();
             
        BinaryVecDataManager<Bool> *pRepetitiveFlagMgr =
          new BinaryVecDataManager<Bool>( strRepetitiveFlagFile );

        pReadIdMgr->UseId( numReads - 1 );

        mpReadDataMgr = new ReadDataManager( mpReadLocationMgr, 
                                             mpReadPairMgr,
                                             pReadNameMgr,
                                             pReadSequenceMgr,
                                             pRepetitiveFlagMgr,
                                             pReadIdMgr );

        this->PrefetchAllReads();
    } 

    
    if ( mVerbose )
        cout << "Setting up ContigDataManager." << endl;

    {
        ContigSequenceManager *pContigSequenceMgr = 
	  new ContigSequenceManager();
        
        IdManager *pContigIdMgr = new IdManager;
        
        mpContigDataMgr = new ContigDataManager( mpContigLocationMgr, 
                                                 mpReadLocationMgr,
                                                 pContigSequenceMgr,
                                                 pContigIdMgr );

        this->PrefetchAllContigs();
    }

    
    if ( mVerbose )
        cout << "Setting up SuperDataManager." << endl;

    {
        mpSuperDataMgr  = new SuperDataManager( mpContigLocationMgr );
    }

}


void
Assembly::ReadIn( const String &strFullWorkPath, const String &strSubDir )
{
    if ( ! IsDirectory( strFullWorkPath ) )
    {
      cout << strFullWorkPath << " is not a directory." << endl;
      TracebackThisProcess();
    }
  
    if ( ! IsDirectory( GetSubdirPath( strFullWorkPath, strSubDir ) ) )
    {
      cout << GetSubdirPath( strFullWorkPath, strSubDir ) << " is not a directory." << endl;
      TracebackThisProcess();
    }

    delete mpContigLocationMgr;
    delete mpReadLocationMgr;
    delete mpReadPairMgr;
    
    delete mpSuperDataMgr;
    delete mpContigDataMgr;
    delete mpReadDataMgr;

    if ( mVerbose )
        cout << "Setting up ContigLocationManager." << endl;

    {
        mpContigLocationMgr = new ContigLocationManager( mpSuperDataMgr, mpContigDataMgr );
    }


    if ( mVerbose )
        cout << "Setting up ReadLocationManager." << endl;

    {
        String strReadLocFile( GetReadLocationsFile( strFullWorkPath, strSubDir ) );
    
	// If the file does not exist, do not try to load it.
	if ( ! IsRegularFile( strReadLocFile ) )
	  strReadLocFile = "";

        mpReadLocationMgr = new ReadLocationManager( mpContigDataMgr,
						     mpReadDataMgr,
						     strReadLocFile );
    }


    if ( mVerbose )
        cout << "Setting up ReadPairManager." << endl;

    {
        String strReadPairFile( GetReadPairingsFile( strFullWorkPath, strSubDir ) );
        
        mpReadPairMgr = new ReadPairManager( mpReadDataMgr,
                                             strReadPairFile );
    }


    if ( mVerbose )
        cout << "Setting up ReadDataManager." << endl;

    {
        IdManager *pReadIdMgr = 
            new IdManager;
        
        String strReadIdsFile( GetReadNamesFile( strFullWorkPath, strSubDir ) );

        NameManager *pReadNameMgr =
            new NameManager( strReadIdsFile );
        
        String strReadLengthsFile( GetReadLengthsFile( strFullWorkPath, strSubDir ) );
        String strReadFastbFile  ( GetReadFastbFile  ( strFullWorkPath, strSubDir ) );
        String strReadQualbFile  ( GetReadQualbFile  ( strFullWorkPath, strSubDir ) );
        String strReadTrimsFile  ( GetReadTrimsFile  ( strFullWorkPath, strSubDir ) );
        String strReadFastnFile  ( GetReadFastnFile  ( strFullWorkPath, strSubDir ) );
 
        ReadSequenceManager *pReadSequenceMgr = 
            new ReadSequenceManager( strReadLengthsFile,
                                     strReadFastbFile,
                                     strReadQualbFile,
                                     strReadTrimsFile,
                                     strReadFastnFile );

        int numReads = MastervecFileObjectCount( strReadFastbFile );

        String strRepetitiveFlagFile( GetReadRepetitiveFlagFile( strFullWorkPath, strSubDir ) );
        if ( ! IsRegularFile( strRepetitiveFlagFile ) )
          strRepetitiveFlagFile.erase();

        BinaryVecDataManager<Bool> *pRepetitiveFlagMgr =
          new BinaryVecDataManager<Bool>( strRepetitiveFlagFile );
        
        pReadIdMgr->UseId( numReads - 1 );

        mpReadDataMgr = new ReadDataManager( mpReadLocationMgr, 
                                             mpReadPairMgr,
                                             pReadNameMgr,
                                             pReadSequenceMgr,
                                             pRepetitiveFlagMgr,
                                             pReadIdMgr );

        this->PrefetchAllReads();
    } 

    
    if ( mVerbose )
        cout << "Setting up ContigDataManager." << endl;

    {
        String strContigFastbFile( GetContigFastbFile( strFullWorkPath, strSubDir ) );
        String strContigQualbFile( GetContigQualbFile( strFullWorkPath, strSubDir ) );

        bool filesExist = ( IsRegularFile( strContigFastbFile ) &&
                            IsRegularFile( strContigQualbFile ) );

        ContigSequenceManager *pContigSequenceMgr = 
          ( filesExist
            ? new ContigSequenceManager( strContigFastbFile,
                                         strContigQualbFile ) 
            : new ContigSequenceManager() );

        // We're about to fetch the length of every contig, so preload the bases.
        if (filesExist)
	  pContigSequenceMgr->LoadAllBases();

        IdManager *pContigIdMgr = new IdManager;
        
        if ( filesExist )
          pContigIdMgr->UseId( MastervecFileObjectCount( strContigFastbFile ) - 1 );
        
        mpContigDataMgr = new ContigDataManager( mpContigLocationMgr, 
                                                 mpReadLocationMgr,
                                                 pContigSequenceMgr,
                                                 pContigIdMgr );

        this->PrefetchAllContigs();
    }

    
    if ( mVerbose )
        cout << "Setting up SuperDataManager." << endl;

    {
        mpSuperDataMgr  = new SuperDataManager( mpContigLocationMgr );
    }


    if ( mVerbose )
        cout << "Filling out ContigLocationManager." << endl;

    {
        String strSupersFile( GetContigLocationsFile( strFullWorkPath, strSubDir ) );

        vec<annotated_supercontig> asupers;

        if ( IsRegularFile( strSupersFile ) ||
             IsRegularFile( strSupersFile + ".gz" ) )
          READX( strSupersFile, asupers ); 
        
        for ( int superId = 0; superId < (int) asupers.size(); superId++)
        {
            annotated_supercontig& as = asupers[superId];
            
            Super theSuper = mpSuperDataMgr->NewSuper( );
            ForceAssertEq( theSuper.GetId(), superId );
            
            int startPos = 0;
            
            for (int contigIdx = 0; contigIdx < as.NumContigs( ); contigIdx++) 
            {
                int contigId = as.Contig(contigIdx).ID();
                
                Contig theContig = mpContigDataMgr->GetContig( contigId );
                
                Interval theInterval( startPos, startPos + theContig.GetLength() );
                
                mpContigLocationMgr->Add( ContigLocation( theSuper, theContig, 
                                                          theInterval, orient_FW ) );
            
                startPos += theContig.GetLength();
                if ( contigIdx < as.NumContigs( ) - 1 )
                    startPos += as.Gap(contigIdx);
            }

	    // Make sure first contig location (in the left-right sense) starts at 0.
	    this->Normalize( theSuper );
        }
    }
}


void 
Assembly::WriteOut( const String &strFullWorkPath, const String &strSubDir,
                    const bool bOverwriteContigsAndSupers, const bool bOverwriteReads ) const
{
    if ( ! IsDirectory( strFullWorkPath ) )
        Mkdir777( strFullWorkPath );

    if ( ! IsDirectory( GetSubdirPath( strFullWorkPath, strSubDir ) ) )
        Mkdir777( GetSubdirPath( strFullWorkPath, strSubDir ) );


    if ( mVerbose )
      cout << "Saving ReadDataManager." << endl;
    {
        String strReadNamesFile  ( GetReadNamesFile  ( strFullWorkPath, strSubDir ) );
        String strReadLengthsFile( GetReadLengthsFile( strFullWorkPath, strSubDir ) );
        String strReadFastbFile  ( GetReadFastbFile  ( strFullWorkPath, strSubDir ) );
        String strReadQualbFile  ( GetReadQualbFile  ( strFullWorkPath, strSubDir ) );
        String strReadTrimsFile  ( GetReadTrimsFile  ( strFullWorkPath, strSubDir ) );
        String strReadFastnFile  ( GetReadFastnFile  ( strFullWorkPath, strSubDir ) );
        
        mpReadDataMgr->Write( bOverwriteReads,
                              strReadNamesFile,
                              strReadLengthsFile,
                              strReadFastbFile,
                              strReadQualbFile,
                              strReadTrimsFile,
                              strReadFastnFile );
    }

    if ( mVerbose )
      cout << "Saving ReadPairManager." << endl;
    {
        String strReadPairFile( GetReadPairingsFile( strFullWorkPath, strSubDir ) );
        
        mpReadPairMgr->Write( bOverwriteReads,
                              strFullWorkPath,
                              strReadPairFile );
    }
    
    if ( mVerbose && mPurgeBeforeWrite )
      cout << "Purging empty contigs (if any)." << endl;

    int contigsPurged = ( mPurgeBeforeWrite ? mpContigDataMgr->PurgeEmptyContigs() : 0 );

    if ( this->GetNumContigs() == 0 )
      cout << "WARNING: there are no contigs in this assembly." << endl;

    if ( mVerbose )
      cout << "Saving ContigLocationManager." << endl;
    {
        String strSupersFile( GetContigLocationsFile( strFullWorkPath, strSubDir ) );

        ForceAssert( bOverwriteContigsAndSupers || ! IsRegularFile( strSupersFile ) );

        String strSummaryFile( GetSubdirPath( strFullWorkPath, strSubDir ) 
                               + "/mergedcontigs.summary" );

        ofstream summaryStrm( strSummaryFile.c_str() );

        summaryStrm <<
            "This file contains a list of all the supercontigs in the final assembly.\n"
            "For each, we give the numbers and lengths of the constituent contigs, as well\n"
            "as the approximate gaps between them (if a positive number is shown in\n"
            "parentheses) or the approximate overlap between them (if a negative number\n"
            "is shown).\n"
            "\n";

        int savedSuperId = 0;
        int emptySuperCount = 0;

        vec<Super> vecSupers;
        this->GetAllSupers( vecSupers );
        
        vec<superb> vecOldSupers;

        vec<ReadPair> allPairs;
        this->GetAllReadPairs( allPairs );
        int maxDev = 0;
        for ( unsigned int i = 0; i < allPairs.size(); ++i )
          if ( maxDev < allPairs[i].GetExpectedStdDev() )
            maxDev = allPairs[i].GetExpectedStdDev();
        maxDev *= 2;
            
        for ( unsigned int superIdx = 0; superIdx < vecSupers.size(); ++superIdx )
        {
            Super &theSuper = vecSupers[ superIdx ];
        
            vec<ContigLocation> vecContigLocs;
            theSuper.GetContigLocations( vecContigLocs );

            if ( vecContigLocs.empty() )
            {
              ++emptySuperCount;
              continue;
            }

            OffsetGraph<Contig> offsetGraph;
            offsetGraph.SetMinLinks(1);
            offsetGraph.SetMaxDev( maxDev );
            offsetGraph.AddSuper( theSuper );

            sort( vecContigLocs.begin(), vecContigLocs.end(), 
                  ContigLocation::OrderBySuperAndLeastGap() );

            vecOldSupers.push_back( superb() );

            ContigLocation firstContigLoc = vecContigLocs.front();
            Contig firstContig = firstContigLoc.GetContig();

	    if ( firstContigLoc.IsReverse() )
	      const_cast<Assembly*>(this)->Reverse( firstContig );

	    vecOldSupers.back().PlaceFirstTig( firstContig.GetId(), firstContig.GetLength() );
            
            ContigLocation prevContigLoc, currContigLoc = firstContigLoc;
            Contig prevContig, currContig = firstContig;
            for ( unsigned int contigIdx = 1; contigIdx < vecContigLocs.size(); ++contigIdx )
            {
                prevContigLoc = currContigLoc;
                prevContig = currContig;

                currContigLoc = vecContigLocs[ contigIdx ];
                currContig    = currContigLoc.GetContig();

                if ( currContigLoc.IsReverse() )
                    const_cast<Assembly*>(this)->Reverse( currContig );

                int gap = currContigLoc.Begin() - prevContigLoc.End();

                // If we do not compute the dev, we set it to be the same as the
                // gap (but at least 1) to indicate its artificiality to anyone
                // who sees it.
                int dev = max( 1, abs(gap) );

                if ( mComputeDevs ) {
                  ContigOffset theOffset;
                  bool found = offsetGraph.FindOffset( prevContig, currContig, theOffset );
                
                  dev = ( found ? theOffset.GetDev() : offsetGraph.GetMaxDev() );
                }

                vecOldSupers.back().AppendTig( currContig.GetId(), currContig.GetLength(),
                                               gap, dev );
            }
        }

        WriteSupercontigFiles( GetSubdirPath( strFullWorkPath, strSubDir ), vecOldSupers );

        if ( emptySuperCount > 0 && mVerbose )
          cout << "  (Ignored " << emptySuperCount << " empty supercontigs.)" << endl;
    }

    if ( mVerbose )
      cout << "Saving ContigDataManager." << endl;
    {
        String strContigFastbFile( GetContigFastbFile( strFullWorkPath, strSubDir ) );
        String strContigQualbFile( GetContigQualbFile( strFullWorkPath, strSubDir ) );

        mpContigDataMgr->Write( bOverwriteContigsAndSupers,
                                strContigFastbFile,
                                strContigQualbFile );
    }

    if ( contigsPurged > 0 && mVerbose )
      cout << "  (Ignored " << contigsPurged << " empty contigs.)" << endl;

    if ( mVerbose )
      cout << "Saving ReadLocationManager." << endl;
    {
        String strReadLocFile( GetReadLocationsFile( strFullWorkPath, strSubDir ) );

        mpReadLocationMgr->Write( bOverwriteContigsAndSupers,
                                  strReadLocFile );
    }
     

}    
