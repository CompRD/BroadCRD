/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// See BadMotif.cc for documentation.

#ifndef BAD_MOTIF_CORE_H
#define BAD_MOTIF_CORE_H

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "TokenizeString.h"
#include "dna/Bases.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

typedef pair<BaseVec::size_type,BaseVec::size_type> SequenceLoc;
typedef SerfVec<SequenceLoc> MatchedLocsVec;
typedef MasterVec<MatchedLocsVec> VecMatchedLocsVec;
extern template class SmallVec< SequenceLoc, MempoolAllocator< SequenceLoc > >;
extern template class OuterVec< MatchedLocsVec >;

class composition_condition
{
    public:
        enum SENSE {EQ, GE, LE, APPROX_EQ, HOMOPOLYMER};
        Bool x[4];
        SENSE sense;
        double bound;
        
        composition_condition( );
        composition_condition( const String& c );
        void ReverseComplement( );
        
        Bool Matches( const int start, const int length,
            const vec< vec<short> >& tally ) const;
};

ostream& operator<<(ostream& s, const composition_condition c);

class FixedLengthMatch
{
    public:
    
        enum MTYPE {NS, BASELIST, COMPOSITION, DISJUNCTION };
        
        int Length( ) const;
        const vec<char>& Bases( ) const;
        char Base( int i ) const;
        MTYPE Mtype( ) const;
        
        void ReverseComplement( );
        
        FixedLengthMatch( const MTYPE mtype, const int length );
        
        FixedLengthMatch( const MTYPE mtype, const vec<char>& bases );
        
        FixedLengthMatch( const MTYPE mtype, const vec<String>& conditions,
          const int length );
        
        Bool Matches( const vec<char>& bases, const int start,
          const vec< vec<short> >& tally ) const;     
        
        friend void Parse( const String& s0, vec<FixedLengthMatch>& m );
        
        
        friend Bool Matches( const vec<char>& bases, int start,
          const vec<FixedLengthMatch>& M, const vec< vec< vec<short> > >& tally );
        
        
        friend ostream& operator<<(ostream& s, const FixedLengthMatch m);
    
    private:
    
        MTYPE mtype_;
        int length_;
        vec<char> bases_;
        vec<composition_condition> conditions_;
};

ostream& operator<<(ostream& s, const vec<FixedLengthMatch> m);

void Parse(const String& s0, vec<FixedLengthMatch>& m);

void ComputeMotif(const String& motif, const vecbasevector& genome,
    const vecbitvector& amb, vecbitvector& in_motif,
    VecMatchedLocsVec& motif_locations);

void ComputeMotif(const String& motif, const BaseVec& chr,
    const BitVec& amb, BitVec& in_motif, MatchedLocsVec& motif_locations);

#endif
