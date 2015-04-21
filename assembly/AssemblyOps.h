// Copyright (c) 2004 The Broad Institute at MIT and Harvard

#ifndef ASSEMBLY_OPS_H
#define ASSEMBLY_OPS_H

#include "assembly/Assembly.h"

// These classes encapsulate operations on assemblies.  It is a
// hierarchy of classes:
//
// AssemblyOp (the abstract base class)
//
// These clases directly descend from AssemblyOp:
//
// ContigSplit (a very flexible splitter of a single contig)
// SuperSplit  (a very flexible splitter of a single super)
//
// Descending from ContigSplit, we have:
//
// ContigSlice    (a more abstract breaker of a single contig)
// ContigTrim     (trims off the ends of a contig)
// ContigSanitize (breaks a given contig into actually contiguous contigs)
// 
// And descending from SuperSplit, we have the analogous:
//
// SuperSlice    (a more abstract breaker of a single super)
// SuperTrim     (trims off the ends of a super)
// SuperSanitize (breaks a given super into actually super-contiguous supers)
//
// Again descending from AssemblyOp:
//
// SuperInsert   (open a supercontig using SuperSlice, insert supercontigs into the
//                opening using SuperAppend, and close)
// SuperAppend   (append one supercontig to another)


// By virtue of their inheritance, you could also store up a set of
// them in a vec<AssemblyOp*> and execute them all at once.  Note that
// they'd have to be non-colliding in that case.

// Conceivably, we could even add an undo operation to the base class,
// allowing one to "try on" a set of modifications, which could then
// be pushed onto a stack to be undone later.  This would be tricky to
// pull off and is left as an exercise to the reader. :)


//------------------------------------------------------------

// The abstract base class.  It takes a single argument in its
// constructor: a pointer to the assembly to be operated upon.

class AssemblyOp
{
 public:
  AssemblyOp( Assembly *p_assembly )
    : mp_assembly( p_assembly )
  {}

  virtual ~AssemblyOp() {}
  
  virtual void Execute( ostream *pLog = 0 ) = 0;


 protected:
  Assembly *mp_assembly;
};

//------------------------------------------------------------

// This class encapsulates the partitioning of a contig's reads into
// one or more new contigs.

class ContigSplit : public AssemblyOp
{
 public:
  // The constructor takes the assembly to be operated upon and a
  // collection of collections of read locations.  All of these read
  // locations must refer to the given contig.
  ContigSplit( Assembly *p_assembly,
               const Contig &theContig,
               const vec< vec<ReadLocation> > &readLocSets );

  virtual ~ContigSplit() {}
  
  // When executed, each sub-collection of read locations will be
  // copied to a new contig, where the bases of the contig are copied
  // from the original contig from the lowest begin on contig to the
  // highest end on contig of all the read locations in the
  // sub-collection.  The original contig is then cleared.
  //
  // An empty sub-collection will create an invalid contig (i.e.,
  // Contig::IsValid() will return false).
  //
  // Note that any read locations in the original contig that are not
  // in any sub-collection are not copied anywhere.  They do not even
  // remain in the original contig.
  //
  // The new contigs are placed in all supers containing the original
  // contig in the corresponding locations.  I.e., if the original
  // contig spanned from 0-1000 in some supercontig and was split into
  // contigs containing its first 300 and last 300 bases, the new
  // contigs would be placed in the same supercontig at 0-300 and
  // 700-1000, respectively.
  virtual void Execute( ostream *pLog = 0 );

  // Copy the contigs created from each sub-collection into the given
  // vector.
  virtual void GetResultingContigs( vec<Contig> &resultingContigs )
  {
    resultingContigs = m_resultingContigs;
  }


 protected:
  Contig m_contig;
  vec< vec<ReadLocation> > m_locSets;
  vec<Contig> m_resultingContigs;
};


//------------------------------------------------------------

// This class encapsulates the partitioning of a super's contigs into
// one or more new supers.

class SuperSplit : public AssemblyOp
{
 public:
  // The arguments and restrictions are similar to those of
  // ContigSplit, except that we are dealing with contig locations
  // rather than read locations.
  SuperSplit( Assembly *p_assembly,
              const Super &theSuper,
              const vec< vec<ContigLocation> > &contigLocSets );

  virtual ~SuperSplit() {}
  
  // When executed, each sub-collection of contig locations will be
  // copied to a new super.  The original super is then cleared.
  //
  // An empty sub-collection will create an invalid super (i.e.,
  // Super::IsValid() will return false).
  //
  // Note that any contig location in the original super that is not
  // in any sub-collection is not copied anywhere.
  virtual void Execute( ostream *pLog = 0 );

  // Copy the supers created from each subcollection into the given
  // vector.
  virtual void GetResultingSupers( vec<Super> &resultingSupers )
  {
    resultingSupers = m_resultingSupers;
  }
    

 protected:
  Super m_super;
  vec< vec<ContigLocation> > m_locSets;
  vec<Super> m_resultingSupers;
};


//------------------------------------------------------------

// This class encapsulates the excision of a portion of a contig,
// creating two new contigs.

class ContigSlice : public ContigSplit
{
 public:
  // The constructor takes the assembly to be operated upon, a contig
  // to be operated upon, and a range of bases to be excised.  These
  // are used to calculate two sets of read locations that, upon
  // execution, will be copied into new contigs as per the rules for
  // ContigSplit above.  Read locations are placed in the two sets as
  // follows:
  //
  // A read location is placed in the first set (the "left" set) if
  // either it ends at or before excisionStart or it begins at or
  // before excisionStart and it is reverse in the contig.
  //
  // A read location is placed in the second set (the "right" set) if
  // either it begins at or after excisionEnd or it ends at or after
  // excisionEnd and it is forward in the contig.
  //
  //           excisionStart  excisionEnd
  //                 |             |
  //  left   <--?--> |             |
  //  left       <---|----         |
  //  left         <-|-------------|---
  //  neither    ----|--->         |
  //  neither        |  <---?--->  |
  //  neither        |         <---|----
  //  right       ---|-------------|->
  //  right          |         ----|--->
  //  right          |             | <--?-->

  ContigSlice( Assembly *p_assembly,
               const Contig &theContig,
               const int excisionStart,
               const int excisionEnd );
  
  virtual ~ContigSlice() {}
};


//------------------------------------------------------------

// This class encapsulates the trimming of the ends of a contig,
// creating a new contig.  The inverse of ContigSlice, essentially.

class ContigTrim : public ContigSplit
{
 public:
  // The constructor takes the assembly to be operated upon, a contig
  // to be operated upon, and the amounts to be trimmed of the left
  // and right ends of the contig.  These define two points in the
  // contig: the "left trim point", which is leftTrim bases after the
  // beginning of the contig, and the "right trim point", which is
  // rightTrim bases before the end of the contig. 
  //
  // If a read location falls within these two points (i.e. the read
  // location begins at or after the left trim point and ends at or
  // before the right trim point), it will be copied into a new contig
  // as per the rules for ContigSplit above.
  //
  // If a read location falls to the left of the left trim point, it
  // will be copied into another contig, and if it falls to the right
  // of the right trim point, a third contig.
  //
  // The contigs will be returned in this order by GetResultingContigs(). 
  ContigTrim( Assembly *p_assembly,
              const Contig &theContig,
              const int leftTrim,
              const int rightTrim );
  
  virtual ~ContigTrim() {}
};


//------------------------------------------------------------

// This class encapsulates the sanitization of a contig;
// the contig is split into new smaller contigs, each of which is
// actually contiguous, i.e. every base is covered by at least one read.

class ContigSanitize : public ContigSplit
{
 public:
  // The constructor takes the assembly to be operated upon and a
  // contig to be operated upon.  The coverage of the contig by its
  // read locations is then computed.  Each zero-coverage region marks
  // a splitting point for the contig.  If a zero-coverage region is
  // at an end of the contig and is no more than endForgiveness bases
  // in size, it is ignored.  If no splitting is necessary, no action
  // is performed when Execute() is called, and GetResultingContigs()
  // will return an empty vector.

  ContigSanitize( Assembly *p_assembly,
		  const Contig &theContig,
                  const int endForgiveness = 0 );

  virtual ~ContigSanitize() {}

  // When executed, no action will be taken if the given contig was
  // covered in its entirety by reads (i.e. if m_needNewContigs is
  // false) otherwise the contig is split as in ContigSplit above.

  virtual void Execute( ostream *pLog = 0 );

 private:
  bool m_needNewContigs;
};


//------------------------------------------------------------

// This class encapsulates the excision of a portion of a supercontig,
// creating two new supercontigs.

class SuperSlice : public SuperSplit
{
 public:
  // The constructor takes the assembly to be operated upon, a super
  // to be operated upon, and a range of bases to be excised.  These
  // are used to calculate two sets of contig locations to be copied
  // into new supers as per the rules for SuperSplit above.  Contig
  // locations are placed in the two sets as follows:
  //
  // A contig location is placed in the first set (the "left" set) if
  // it ends at or before excisionStart.
  //
  // A contig location is placed in the second set (the "right" set) if
  // it begins at or after excisionEnd.
  //
  // If a contig intersects either or both excisionStart and
  // excisionEnd, new contigs will be created (upon execution)
  // containing only the portions of the contig falling to the left of
  // excisionStart (which will be put in the "left" set) and to the
  // right of excisionEnd (which will be placed in the "right" set).
  //
  // Contigs (either original or created by the splitting described in
  // the previous paragraph) that fall entirely between excisionStart
  // and excisionEnd are removed from the assembly.
  //
  //           excisionStart  excisionEnd
  //                 |             |
  //  L      ------  |             |
  //  L/neither  ----|----         |
  //  L/neither/R  --|-------------|---
  //  neither        |  -------    |
  //  neither/R      |          ---|----
  //  R              |             |  ------

  SuperSlice( Assembly *p_assembly,
              const Super &theSuper,
              const int excisionStart,
              const int excisionEnd );

  virtual ~SuperSlice() {}
  

  // When executed, the contigs resulting from excision will be
  // created and the "left" and "right" sets of contig locations are
  // then processed as in SuperSplit above.

  virtual void Execute( ostream *pLog = 0 );


 protected:
  vec<ContigSlice> m_contigSlices;
  vec<Contig> m_contigsToRemove;
};


//------------------------------------------------------------

// This class encapsulates the trimming of the ends of a supercontig,
// creating a new supercontig.

class SuperTrim : public SuperSplit
{
 public:
  // The constructor takes the assembly to be operated upon, a super
  // to be operated upon, and the number of bases to be removed from
  // each end.  These are used to calculate a set of contig locations
  // to be copied into a new super as per the rules for SuperSplit
  // above.  Contig locations are placed in the set if, according to
  // the rules for SuperSlice, they would be placed in neither the
  // left nor the right set.
  //
  // If a contig intersects a trim point, new contigs will be created
  // (upon execution) containing only the portions of the contig
  // falling inside the trimming region.
  // 
  // Three supers will be created (and returned in this order by
  // GetResultingSupers()): the trimmed super, the left trimming, and
  // the right trimming.

  SuperTrim( Assembly *p_assembly,
             const Super &theSuper,
             const int leftTrim,
             const int rightTrim );
  
  virtual ~SuperTrim() {}
  

  // When executed, the contigs resulting from excision will be
  // created and the "to be saved" set of contig locations are then
  // processed as in SuperSplit above.

  virtual void Execute( ostream *pLog = 0 );


 protected:
  vec<ContigTrim> m_contigTrims;
};


//------------------------------------------------------------

// This class encapsulates the sanitization of a super; the super is
// split into new smaller supers, each of which is actually
// connected, i.e. every contig is directly or indirectly
// connected by valid links to every other contig.

class SuperSanitize : public SuperSplit
{
 public:
  // The constructor takes the assembly to be operated upon and a
  // super to be operated upon.  If maxStretch is > 0, then only
  // logical links with stretch <= maxStrech are considered. Two
  // contigs are linked iff there are >= minLinks between them.  If
  // minShortLinks > 0, and if only short links connect the contigs,
  // then there must be >= minShortLinks links for the two
  // contigs to be considered connected.
  
  SuperSanitize( Assembly *p_assembly,
                 const Super &theSuper,
                 const int minLinks = 1,
		 const double maxStretch = -1.0,
		 const int minShortLinks = 0,
		 const int lenShortLink = 10000 );
  
  virtual ~SuperSanitize() {}

  // When executed, no action will be taken if the given super was
  // already connected; otherwise the super is split as in SuperSplit
  // above.

  virtual void Execute( ostream *pLog = 0 );
};

//------------------------------------------------------------

//  Remove ends contigs from a supercontig if the contigs are < minContigSize.  Contigs
//  "popped off" are placed in singleton supercontigs. 
// 

class SuperPop : public SuperSplit
{
 public:
  SuperPop( Assembly *p_assembly,
	    const Super &theSuper,
	    const int minContigSize );

  virtual ~SuperPop() {}

   virtual void Execute( ostream *pLog = 0 );

 private:
   int m_minContigSize;
};
	    


//------------------------------------------------------------

// Append one supercontig to another.  Gap to and orientation of appended supercontig must
// be provided.  Appended supercontig is cleared.

class SuperAppend : public AssemblyOp
{
 public:

  SuperAppend( Assembly *p_assembly,
	       const Super &superToAppendTo,
	       const Super &superToAppend,
	       const int gap,
	       const Orientation theOrientation );
   
  virtual ~SuperAppend() {}

  virtual void Execute( ostream *pLog = 0 );

 private:
  Super m_superToAppendTo, m_superToAppend;
  int m_gap;
  Orientation m_orientation;

};

//------------------------------------------------------------

// Append one supercontig to another.  Gap to and orientation of appended supercontig must
// be provided.  Appended supercontig is cleared.

class SuperSort : public AssemblyOp
{
 public:

  SuperSort( Assembly *p_assembly,
	     const Super &superToSort);
   
  virtual ~SuperSort() {}

  virtual void Execute( ostream *pLog = 0 );

 private:
  Super m_superToSort;

};

//------------------------------------------------------------

// Slice open a supercontig (using SuperSlice), insert supercontigs into the
// opening (using SuperAppend) with specified gaps and orientations, and close.
// 
//                insertion point
// super1 AAAAAAAAAAAA|BBBBBBBBBBBBB
// super_to_insert1 *****@@@ (FW)
// super_to_insert2 ++++--- (RC)
//
// resulting super:AAAAAAAAAAAA(gap)*****@@@(gap)---++++(gap)BBBBBBBBBBBBB
//
//
// Specified insertion point is in supercontig coordinates (could be mid-contig).
// Inserted supercontigs are cleared.
//

class SuperInsert : public AssemblyOp
{
 public:

  SuperInsert( Assembly *p_assembly,
	       const Super &superToInsertInto,
	       const vec<Super> &supersToInsert,
	       const int insertionPoint,
	       const vec<int> gaps,
	       const vec<Orientation> theOrientations );
   
  virtual ~SuperInsert() {}

  virtual void Execute( ostream *pLog = 0 );

 private:
  Super m_superToInsertInto;
  vec<Super> m_supersToInsert;
  int m_insertionPoint;
  vec<int> m_gaps;
  vec<Orientation> m_orientations;

};

#endif
