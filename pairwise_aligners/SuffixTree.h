/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file SuffixTree.h
 * \author tsharpe
 * \date Apr 15, 2009
 *
 * \brief Suffix tree.
 *
 * Data structure for exact string matching.
 */
#ifndef PAIRWISE_ALIGNERS_SUFFIXTREE_H_
#define PAIRWISE_ALIGNERS_SUFFIXTREE_H_

#include "Basevector.h"
#include "system/Assert.h"
#include <list>
#include <string>

/**
 * Implementation of a suffix tree for basevectors.
 * This is a data structure for exact (sub)string matching.
 *
 * This initial version doesn't make any real attempt to be memory efficient (lots of pointers),
 * and it's an implicit Suffix Tree.  The code to make it explicit is there, but right now the
 * nodes aren't storing the string information that would make it useful, so it's commented out.
 */
class SuffixTree
{
public:
    typedef unsigned int uint;

    /// Make an empty one.
    SuffixTree();

    /// Destroy one.
    ~SuffixTree();

    /// Add sequence to an implicit suffix tree.
    void addImplicit( basevector const& seq )
    {
        mSequences.push_back(new bvec(seq));
        buildImplicitTree(*mSequences.back());
    }

    /// Add all sequences to an implicit suffix tree.
    /// Itr is an input-iterator type that derefs to a bvec.
    template <class Itr>
    void addImplicit( Itr start, Itr const& end )
    { while ( start != end ) { addImplicit(*start); ++start; } }

#if 0
    // wait for Node to track it's suffixes

    /// Add a sequence to the tree, and mark its suffixes.
    void add( basevector const& seq )
    {
        addImplicit(seq);
        markTree(*mSequences.back());
    }

    /// Add all sequences to the tree, and mark its suffixes.
    /// Itr is an input-iterator type that derefs to a bvec.
    template <class Itr>
    void add( Itr start, Itr end )
    { while ( start != end ) { add(*start); ++start; } }
#endif

    /// Test for a match.
    /// In detail, tests whether there's an exact match having a specified minimum length,
    /// among any of the sequences in the tree, for any portion of the probe.
    bool matches( basevector const& probe, uint minLen ) const;

    /// check the tree structure (for debugging)
    void validate()
    { validateSuffixLinks(&mRoot,bvec());
      validateSubsequences(); }

    /// print the tree structure (for entertainment)
    void printTree()
    { Node const*const* end = mRoot.endChild();
      for ( Node const*const* itr = mRoot.beginChild(); itr != end; ++itr )
      { printTree(*itr,""); } }

    class SeqItr : public std::iterator<std::input_iterator_tag,bvec,void,bvec const*,bvec const&>
    {
    public:
        SeqItr(list<bvec*>::const_iterator itr)
        : mItr(itr)
        {}

        bvec const& operator*() const
        { return **mItr; }

        bvec const* operator->() const
        { return *mItr; }

        bool operator==( SeqItr const& that ) const
        { return mItr == that.mItr; }

        bool operator!=( SeqItr const& that ) const
        { return mItr != that.mItr; }

        SeqItr& operator++()
        { ++mItr; return *this; }

        SeqItr operator++(int)
        { SeqItr tmp(mItr); ++mItr; return tmp; }

    private:
        list<bvec*>::const_iterator mItr;
    };

    /// iterator over the sequences stored in the tree
    SeqItr beginSequences()
    { return SeqItr(mSequences.begin()); }

    /// endpoint for iteration over the sequences stored in the tree
    SeqItr endSequences()
    { return SeqItr(mSequences.end()); }

    class Node
    {
    public:
        Node()
        : mLabel(&ROOT_BVEC), mLabelStart(0), mLabelEnd(0), mSuffixLink(0)
        { memset(mChildren,0,sizeof(mChildren)); }

        Node( bvec const& label, uint labelStart )
        : mLabel(&label), mLabelStart(labelStart), mLabelEnd(label.size()), mSuffixLink(0)
        { memset(mChildren,0,sizeof(mChildren)); }

        Node( bvec const& label, uint labelStart, uint labelEnd )
        : mLabel(&label), mLabelStart(labelStart), mLabelEnd(labelEnd), mSuffixLink(0)
        { memset(mChildren,0,sizeof(mChildren)); }

        ~Node()
        { delete mChildren[0]; delete mChildren[1]; delete mChildren[2]; delete mChildren[3]; }

        /// is this a leaf node
        bool isLeaf() const
        { return mLabelEnd == mLabel->size(); }

        /// starting-position for iterating through the label
        bvec::const_iterator beginLabel() const
        { return mLabel->Begin(mLabelStart); }

        /// ending-position for iterating through the label
        bvec::const_iterator endLabel() const
        { return mLabel->Begin(mLabelEnd); }

        /// get the base at a given label offset
        base_t getBase( uint posn ) const
        { return mLabel->operator[](mLabelStart+posn); }

        /// length of the label on this node
        uint getLabelSize() const
        { return mLabelEnd - mLabelStart; }

        /// length of the label, capped for phase
        uint getLabelSize( bvec const& seq, uint maxEnd ) const
        { return (&seq == mLabel ? min(mLabelEnd,maxEnd) : mLabelEnd) - mLabelStart; }

        /// get specified child
        Node* getChild( base_t baseVal )
        { return mChildren[baseVal]; }

        /// get specified child
        Node const* getChild( base_t baseVal ) const
        { return mChildren[baseVal]; }

        /// set the specified child
        void setChild( base_t baseVal, Node* pNode )
        { mChildren[baseVal] = pNode; }

        /// beginning-position for iterating through children
        Node const*const* beginChild() const
        { return mChildren; }

        /// ending-position for iterating through children
        Node const*const* endChild() const
        { return mChildren+sizeof(mChildren)/sizeof(mChildren[0]); }

        /// get suffix link
        Node* getSuffixLink()
        { return mSuffixLink; }

        /// get suffix link
        Node const* getSuffixLink() const
        { return mSuffixLink; }

        /// set suffix link
        void setSuffixLink( Node* pNode )
        { mSuffixLink = pNode; }

        /// divide this node at the given label offset
        /// new Node is the head portion of the label, this becomes the tail
        Node* split( uint posn )
        {
            Node* pNode = new Node(*mLabel,mLabelStart,mLabelStart+posn);
            mLabelStart += posn;
            pNode->setChild(mLabel->operator[](mLabelStart),this);
            return pNode;
        }

        /// swizzle label sequence
        void setLeafLabel( bvec const& label, uint labelStart )
        {
            AssertEq(bvec(*mLabel,mLabelStart,mLabelEnd-mLabelStart),bvec(label,labelStart,mLabelEnd-mLabelStart));
            mLabel = &label;
            mLabelStart = labelStart;
            mLabelEnd = label.size();
        }

    private:
        Node( Node const& ); // unimplemented -- no copying
        Node const& operator=( Node const& ); // unimplemented -- no copying

        bvec const* mLabel;
        unsigned int mLabelStart;
        unsigned int mLabelEnd;
        Node* mChildren[4];
        Node* mSuffixLink;
        static bvec const ROOT_BVEC;
    };

private:
    SuffixTree( SuffixTree const& ); // unimplemented -- no copying
    SuffixTree const& operator=( SuffixTree const& ); // unimplemented -- no copying

    void buildImplicitTree( basevector const& seq );
    //void markTree( basevector const& seq );
    void validateSuffixLinks( Node const* node, bvec const& path );
    void validateSubsequences();
    void printTree( Node const* root, std::string prefix );

    list<bvec*> mSequences;
    Node mRoot;
};

#endif /* PAIRWISE_ALIGNERS_SUFFIXTREE_H_ */
