/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file SuffixTree.cc
 * \author tsharpe
 * \date Apr 22, 2009
 *
 * \brief Suffix tree.
 */
#include "pairwise_aligners/SuffixTree.h"
#include <string.h>

using std::string;
using std::cout;
using std::endl;

namespace
{

typedef SuffixTree::uint uint;
typedef SuffixTree::Node Node;

class Posn
{
public:
    Posn( Node const* pRoot )
    : mRoot(pRoot)
    {
        init();
    }

    uint walk( basevector const& seq )
    {
        init();
        bvec::const_iterator end(seq.End());
        for ( bvec::const_iterator itr(seq.Begin()); itr != end; ++itr )
            if ( !step(*itr) )
                break;

        return mPathLen + mIdx;
    }

    bool step( base_t base )
    {
        bool result = false;

        uint labLen = mNode->getLabelSize();
        if ( labLen == mIdx )
        {
            Node const* pNode = mNode->getChild(base);
            if ( pNode )
            {
                mPrevNode = mNode;
                mNode = pNode;
                mPathLen += labLen;
                mIdx = 1;
                result = true;
            }
        }
        else if ( mNode->getBase(mIdx) == base )
        {
            mIdx += 1;
            result = true;
        }

        return result;
    }

    void skip( basevector const& seq, uint start, uint stop )
    {
        if ( stop < start )
        {
            init();
        }
        else
        {
            uint depth = stop - start;
            uint labLen;
            while (  mPathLen+(labLen = mNode->getLabelSize()) < depth )
            {
                mPathLen += labLen;
                mPrevNode = mNode;
                mNode = mNode->getChild(seq[start+mPathLen]);
            }
            mIdx = depth - mPathLen;
        }
    }

    void toSuffix()
    {
        Node const* pNode = mPrevNode->getSuffixLink();
        if ( pNode && pNode != mRoot )
        {
            mPrevNode = mRoot;
            mNode = pNode;
            mPathLen -= pNode->getLabelSize() + 1;
            mIdx = 0;
        }
        else
        {
            init();
        }
    }

private:
    void init()
    {
        mPrevNode = mRoot;
        mNode = mRoot;
        mPathLen = 0;
        mIdx = 0;
    }

    Node const* mRoot;
    Node const* mNode;
    Node const* mPrevNode;
    uint mPathLen;
    uint mIdx;
};

} // end of anonymous namespace

bvec const SuffixTree::Node::ROOT_BVEC(1U);

SuffixTree::SuffixTree()
{}

SuffixTree::~SuffixTree()
{
    list<bvec*>::iterator end(mSequences.end());
    for ( list<bvec*>::iterator itr(mSequences.begin()); itr != end; ++itr )
        delete *itr;
}

bool SuffixTree::matches( basevector const& probe, uint minLen ) const
{
    Posn pos(&mRoot);
    uint len = probe.size();
    uint jjj = 0;
    for ( uint iii = 0; iii < len; ++iii )
    {
        if ( jjj < iii )
            jjj = iii;

        while ( jjj < len && pos.step(probe[jjj]) )
        {
            jjj += 1;
            if ( jjj - iii >= minLen )
                return true;
        }

        pos.toSuffix();
        pos.skip(probe,iii+1,jjj);
    }

    return false;
}

/*
 * This implements Ukkonen's algorithm for building an implicit suffix tree in linear time.
 * See <i>Algorithms on Strings, Trees, and Sequences</i>
 * by Dan Gusfield, Cambridge Univ. Press, reprinted in 1999, pp. 94-107.
 */
void SuffixTree::buildImplicitTree( basevector const& seq )
{
    base_t base1 = 4; // initialized to illegal value.  algorithm always overwrites before using, but compiler doesn't know that.
    base_t base2;
    Node* pNewInternalNode = 0;
    uint len = seq.size();
    uint minExIdx = 0;

    for ( uint phaseIdx = Posn(&mRoot).walk(seq); phaseIdx < len; ++phaseIdx ) // phase loop
    {
        uint seqIdx = minExIdx;
        uint nToMatch = phaseIdx - minExIdx;

        Node* pPrevNode = 0;
        Node* pNode = &mRoot;
        uint labLen = 0;

        for ( uint exIdx = minExIdx; exIdx <= phaseIdx; ++exIdx ) // extension loop
        {
            // skip/count to find seq[exIdx..phaseIdx] (trick 1)
            while ( nToMatch > labLen )
            {
                seqIdx += labLen;
                nToMatch -= labLen;
                pPrevNode = pNode;
                base1 = seq[seqIdx];
                pNode = pNode->getChild(base1);

                // this is our version of trick 3 -- we build a leaf at full length, but adjust
                // its effective length to keep it bounded for the phase we're executing.
                labLen = pNode->getLabelSize(seq,phaseIdx);
            }

            // extend to include seq[phaseIdx]
            if ( nToMatch != labLen ) // we're midway through an edge
            {
                base2 = seq[phaseIdx];
                if ( base2 != pNode->getBase(nToMatch) ) // Rule 2a -- ended partway through an edge
                {                                        // and the next character doesn't match
                    pNode = pNode->split(nToMatch);
                    pPrevNode->setChild(base1,pNode);
                    pNode->setChild(base2,new Node(seq,phaseIdx));
                    if ( pNewInternalNode )
                        pNewInternalNode->setSuffixLink(pNode);
                    pNewInternalNode = pNode;
                }
                else // Rule 3 -- found it in the tree already (within some edge)
                {
                    minExIdx = exIdx;
                    break; // (trick 2)
                }
            }
            else if ( pNode->isLeaf() ) // Rule 1 -- ended at a leaf-end
            {
                // this is a little something extra: since we're building generalized trees
                // we may have to swizzle the current leaf if it's on some other string
                pNode->setLeafLabel(seq,seqIdx);
            }
            else
            {
                if ( pNewInternalNode )
                {
                    pNewInternalNode->setSuffixLink(pNode);
                    pNewInternalNode = 0;
                }

                base2 = seq[phaseIdx];
                if ( !pNode->getChild(base2) ) // Rule 2b -- at end of edge, but there's no
                {                              // path that continues from there
                    pNode->setChild(base2,new Node(seq,phaseIdx));
                }
                else // Rule 3 -- found it in the tree already (at the end of some edge)
                {
                    minExIdx = exIdx;
                    break; // (trick 2)
                }
            }

            if ( exIdx == phaseIdx )
            {
                break;
            }

            // use suffix links to find starting point for next extension
            if ( (pNode = pNode->getSuffixLink()) )
            {
                seqIdx += nToMatch;
            }
            else if ( pPrevNode == &mRoot )
            {
                pNode = pPrevNode;
                seqIdx = exIdx + 1;
            }
            else
            {
                pNode = pPrevNode->getSuffixLink();
            }
            labLen = pNode->getLabelSize();
            seqIdx -= labLen;
            nToMatch = phaseIdx - seqIdx;
        }
    }

    Assert(!pNewInternalNode);
}

#if 0
/*
 * This is like an extra, final phase, extending the tree to encompass a fictitious terminator character.
 * (Gusfield, op. cit., p. 107.)
 * (I.e., it transforms an implicit suffix tree for the string into an explicit one.)
 * While we're doing this, we mark each node that terminates a suffix of this string as such.
 * Note that we're making trees sort of halfway between implicit and explicit trees:  Each node has a list
 * of strings that terminate there (including interior nodes), and we don't actually put the terminator
 * character into the leaves.
 */
void SuffixTree::markTree( basevector const& seq )
{
    uint len = seq.size();
    uint seqIdx = 0;
    uint nToMatch = len;

    Node* pPrevNode = 0;
    Node* pNode = &mRoot;
    uint labLen = 0;

    char base;
    Node* pNewInternalNode = 0;

    for ( uint exIdx = 0; exIdx < len; ++exIdx )
    {
        while ( nToMatch > labLen )
        {
            seqIdx += labLen;
            nToMatch -= labLen;
            pPrevNode = pNode;
            base = seq[seqIdx];
            pNode = pNode->getChild(base);
            labLen = pNode->getLabelSize();
        }

        if ( nToMatch != labLen )
        {
            pNode = pNode->split(nToMatch);
            pPrevNode->setChild(base,pNode);
            if ( pNewInternalNode )
            {
                pNewInternalNode->setSuffixLink(pNode);
            }
            pNewInternalNode = pNode;
        }
        else if ( pNewInternalNode )
        {
            pNewInternalNode->setSuffixLink(pNode);
            pNewInternalNode = 0;
        }

        //pNode->addString(seq);

        pNode = pNode->getSuffixLink();
        if ( pNode )
        {
            seqIdx += nToMatch;
        }
        else if ( pPrevNode == &mRoot )
        {
            pNode = pPrevNode;
            seqIdx = exIdx + 1;
        }
        else
        {
            pNode = pPrevNode->getSuffixLink();
        }
        labLen = pNode->getLabelSize();
        seqIdx -= labLen;
        nToMatch = len - seqIdx;
    }

    if ( pNewInternalNode )
    {
        pNewInternalNode->setSuffixLink(&mRoot);
    }
}
#endif

void SuffixTree::validateSuffixLinks( Node const* pNode, bvec const& path )
{
    if ( !pNode )
        return;

    if ( pNode == &mRoot )
    {
        Node const*const* endChild = pNode->endChild();
        for ( Node const*const* itrChild = pNode->beginChild(); itrChild != endChild; ++itrChild )
            validateSuffixLinks(*itrChild,path);
    }
    else if ( !pNode->isLeaf() )
    {
        bvec myPath(path);
        bvec::const_iterator end(pNode->endLabel());
        for ( bvec::const_iterator itr(pNode->beginLabel()); itr != end; ++itr )
            myPath.AppendBase(*itr);

        Node const* pSuffixNode = pNode->getSuffixLink();
        if ( !pSuffixNode )
        {
            cout << "Internal node has no suffix link." << endl;
        }
        else
        {
            bvec::const_iterator end(myPath.End());
            bvec::const_iterator itr(myPath.Begin(1U)); // Note:  beginning walk at position 1

            if ( itr == end )
            {
                if ( pSuffixNode != &mRoot )
                {
                    cout << "Expected root as suffix link of a length==1 node." << endl;
                }
            }
            else
            {
                Node const* pTest = &mRoot;
                bvec::const_iterator testEnd(pTest->endLabel());
                bvec::const_iterator testItr(pTest->beginLabel());
                for ( ; itr != end; ++itr )
                {
                    if ( testItr == testEnd )
                    {
                        pTest = pTest->getChild(*itr);
                        if ( !pTest )
                        {
                            cout << "Unable to walk to suffix." << endl;
                        }
                        else
                        {
                            testEnd = pTest->endLabel();
                            testItr = pTest->beginLabel();
                            Assert( testItr != testEnd );
                        }
                    }
                    if ( *itr != *testItr )
                    {
                        cout << "Suffix walk doesn't compare." << endl;
                    }
                    ++testItr;
                }
                if ( pTest != pSuffixNode )
                {
                    cout << "Suffix link is wrong." << endl;
                }
            }
        }

        Node const*const* endChild = pNode->endChild();
        for ( Node const*const* itrChild = pNode->beginChild(); itrChild != endChild; ++itrChild )
            validateSuffixLinks(*itrChild,myPath);
    }
}

void SuffixTree::validateSubsequences()
{
    SeqItr endSeq(endSequences());
    for ( SeqItr itrSeq(beginSequences()); itrSeq != endSeq; ++itrSeq )
    {
        bvec seq(*itrSeq);
        while ( seq.size() > 0 )
        {
            if ( !matches(seq,seq.size()) )
            {
                cout << "Failed to match a subsequence." << endl;
            }
            seq = bvec(seq,1,seq.size()-1);
        }
    }
}

void SuffixTree::printTree( Node const* pNode, string prefix )
{
    if ( !pNode )
        return;

    bvec::const_iterator endLabel(pNode->endLabel());
    cout << prefix;
    for ( bvec::const_iterator itrLabel(pNode->beginLabel()); itrLabel != endLabel; ++itrLabel )
        cout << as_base(*itrLabel);
    cout << '\n';

    if ( !pNode->isLeaf() )
    {
        prefix += ' ';
        Node const*const* endChild = pNode->endChild();
        for ( Node const*const* itrChild = pNode->beginChild(); itrChild != endChild; ++itrChild )
            printTree(*itrChild,prefix);
    }
}
