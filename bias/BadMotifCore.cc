/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "bias/BadMotifCore.h"
#include "feudal/OuterVecDefs.h"
#include "feudal/SmallVecDefs.h"

template class SmallVec< SequenceLoc, MempoolAllocator< SequenceLoc > >;
template class OuterVec< MatchedLocsVec >;

composition_condition::composition_condition()
{
    return;
}

composition_condition::composition_condition(const String& c)
{    
    String sum;
    if (c.Contains( ">=" ))
    {
        sum = c.Before( ">=" );
        bound = c.After( ">=" ).Double( );
        sense = GE;
    }
    else if ( c.Contains( "<=" ) )
    {
        sum = c.Before( "<=" );
        bound = c.After( "<=" ).Double( );
        sense = LE;
    }
    else if (c.Contains("="))
    {
        sum = c.Before( "=" );
        bound = c.After( "=" ).Double( );
        sense = EQ;
    }
    else if (c.Contains("%"))
    {
        sum = c.Before("%");
        bound = c.After("%").Double();
        sense = APPROX_EQ;
    }
    else if (c == "HOMOPOLYMER")
    {
        sense = HOMOPOLYMER;
        return;
    }
    else
    {
        std::cerr << "illegal operator in " << c << std::endl;
        CRD::exit(EXIT_FAILURE);
    }

    for ( int i = 0; i < 4; i++ )
    {
        x[i] = False;
    }
    for (size_t i = 0; i < sum.size(); i++ )
    { 
        char c = sum[i];
        ForceAssert(Base::isCanonicalBase(c));
        int j = as_char(c);
        ForceAssert(!x[j]);
        x[j] = True;
        if ( i < sum.size( ) - 1 )
        {
            ForceAssert( sum[++i] == '+' );
        }
    }
}

void composition_condition::ReverseComplement()
{
    swap( x[0], x[3] );
    swap( x[1], x[2] );
}
   
Bool composition_condition::Matches(const int start, const int length,
    const vec<vec<short> >& tally) const
{
    if (sense == HOMOPOLYMER)
    {
        size_t h = 5;
        for (size_t i = 0; i < 4 && h == 5; i++)
        {
            if (tally[i][start] == length)
            {
                h = i;
            }
            else if (tally[i][start] != 0)
            {
                return False;
            }
        }
        
        if ((start == 0 || tally[h][start - 1] < tally[h][start]) &&
            (start + 1 >= tally[h].isize() ||
            tally[h][start + 1] < tally[h][start]))
        {
            return True;
        }
        
        return False;
    }
    else
    {
        int sum = 0;
        for (int i = 0; i < 4; i++)
        {
            if (x[i])
            {
                sum += tally[i][start];
            }
        }
        if (sense == EQ)
        {
            return double(sum) / double(length) == bound;
        }
        else if (sense == GE)
        {
            return double(sum) / double(length) >= bound;
        }
        else if (sense == LE)
        {
            return double(sum) / double(length) <= bound;
        }
        else if (sense == APPROX_EQ)
        {
            return abs(double(sum) / double(length) - bound) <= 0.025;
        }
        else
        {
            CRD::exit(EXIT_FAILURE);
            return False;
        }
    }
}

ostream& operator<<(ostream& s, const composition_condition c)
{
    bool prev_val = false;
    for (int i = 0; i < 4; i++)
    {
        if (c.x[i])
        {
            if (prev_val)
            {
                s << "+";
            }
            s << Base::val2Char(i);
            prev_val = true;
        }
    }
    
    if (c.sense == composition_condition::EQ)
    {
        s << "=";
    }
    else if (c.sense == composition_condition::GE)
    {
        s << ">=";
    }
    else if (c.sense == composition_condition::LE)
    {
        s << "<=";
    }
    else if (c.sense == composition_condition::APPROX_EQ)
    {
        s << "%";
    }
    
    s << c.bound;
    
    return s;
}

int FixedLengthMatch::Length( ) const
{
    return length_;
}

const vec<char>& FixedLengthMatch::Bases( ) const
{
    return bases_;
}

char FixedLengthMatch::Base( int i ) const
{
    return bases_[i];
}

FixedLengthMatch::MTYPE FixedLengthMatch::Mtype( ) const
{
    return mtype_;
}

void FixedLengthMatch::ReverseComplement()
{    
    if (mtype_ == BASELIST)
    {
        GeneralizedBase::reverseComplement(bases_.begin(), bases_.end());
    }
    else if (mtype_ == COMPOSITION || mtype_ == DISJUNCTION)
    {    
        for (size_t i = 0; i < conditions_.size(); i++)
        {
            conditions_[i].ReverseComplement();
        }
    }
}

FixedLengthMatch::FixedLengthMatch( const MTYPE mtype, const int length )
{    
    ForceAssert( mtype == NS );
    mtype_ = mtype;
    length_ = length;
}

FixedLengthMatch::FixedLengthMatch( const MTYPE mtype, const vec<char>& bases )
{    
    ForceAssert( mtype == BASELIST );
    mtype_ = mtype;
    bases_ = bases;
    length_ = bases.size( );
}

FixedLengthMatch::FixedLengthMatch(const MTYPE mtype,
                                   const vec<String>& conditions,
                                   const int length) : mtype_(mtype), 
                                                       length_(length)
{    
    ForceAssert(mtype == COMPOSITION || mtype == DISJUNCTION);
    for (size_t i = 0; i < conditions.size(); i++)
    {
        conditions_.push(conditions[i]);
    }
}

Bool FixedLengthMatch::Matches(const vec<char>& bases, const int start,
    const vec<vec<short> >& tally) const
{    
    if (mtype_ == NS)
    {
        return True;
    }
    else if (mtype_ == BASELIST)
    {    
        for (size_t i = 0; i < bases_.size(); i++)
        {
            if (Base(i) != 'N' && Base(i) != bases[i+start])
            {
                return False;
            }
        }
        return True;
    }
    else if (mtype_ == COMPOSITION)
    {
        for (size_t i = 0; i < conditions_.size(); i++)
        {
            if (!conditions_[i].Matches(start, length_, tally))
            {
                return False;
            }
        }
        return True;
    }
    else if (mtype_ == DISJUNCTION)
    {
        for (size_t i = 0; i < conditions_.size(); i++)
        {
            if (conditions_[i].Matches(start, length_, tally))
            {
                return True;
            }
        }
        return False;
    }
    
    ForceAssert(0 == 1); 
    return False;
}

int read_exponent(String& s, size_t& i)
{
    int exp = 1;                                        
    if ( i < s.size( ) - 1 && s[i+1] == '^' )
    {
        String digits;                                
        size_t j;                                           
        for ( j = i+2; j < s.size() && isdigit(s[j]); j++ )             
        {    
            digits += s[j];    
        }                        
        i = j-1;                                         
        ForceAssert( digits.size( ) >= 1 );              
        exp = digits.Int( );
    }
    return exp;
}

void Parse( const String& s0, vec<FixedLengthMatch>& m )
{    
    String s = WhiteSpaceFree(s0);
    m.clear( );
    vec<char> bases;
    for (size_t i = 0; i < s.size(); i++)
    {
        char c = s[i];
        if ( c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' )
        {
            int exp = read_exponent(s, i);
            if ( c == 'N' )
            {
                if ( bases.nonempty( ) ) m.push( FixedLengthMatch::BASELIST,
                    bases );
                bases.clear( );
                m.push( FixedLengthMatch::NS, exp );
            }
            else bases.push_back_copies( c, exp );
        }
        else
        {
            ForceAssert( c == '(' );
            size_t j = i + 1;
            while (j < s.size() && s[j] != ')')
            {
                j++;
            }

            if (j == s.size())
            {
                cout << "Unmatched parentheses.\nAbort." << endl;
                exit(1);
            }
            String inside = s.substr( i+1, j-i-1 );
            i = j;
            int exp = read_exponent(s, i);
            if (!inside.Contains("=") && !inside.Contains("%") &&
                inside != "HOMOPOLYMER")
            {
                for ( int u = 0; u < exp; u++ )
                {   
                    for (size_t k = 0; k < inside.size(); k++)
                    {
                        char c = inside[k];
                        ForceAssert( Base::isCanonicalBase(c) || c == 'N' );
                        bases.push_back(c);   
                    }
                }
            }
            else
            {
                if (bases.nonempty())
                {
                    m.push(FixedLengthMatch::BASELIST, bases);
                }
                
                bases.clear();

                if (exp > 30000)
                {
                    std::cerr << "Maximum exponent for fixed length base "
                         << "content restrictions is 30000." << endl;
                    std::cerr << "Abort." << endl;
                    CRD::exit(1);
                }
                
                if (inside.Contains("|") && inside.Contains("&"))
                {
                    std::cerr << inside << " cannot be both a conjunction "
                        "and a disjunction" << std::endl;
                    CRD::exit(1);
                }
                else if (inside.Contains("|"))
                {
                    vec<String> parts;
                    Tokenize(inside, '|', parts);
                    m.push(FixedLengthMatch::DISJUNCTION, parts, exp);
                }
                else
                {
                    vec<String> parts;
                    Tokenize(inside, '&', parts);
                    m.push(FixedLengthMatch::COMPOSITION, parts, exp);
                }
            }
        }
    }
    if ( bases.nonempty( ) ) m.push( FixedLengthMatch::BASELIST, bases );    
}

Bool Matches(const vec<char>& bases, int start,
          const vec<FixedLengthMatch>& M, const vec< vec< vec<short> > >& tally)
{
    for (size_t i = 0; i < M.size(); i++)
    {
        if ( !M[i].Matches( bases, start, tally[i] ) ) return False;
               start += M[i].Length( );
    }
    return True;
}

ostream& operator<<(ostream& s, const FixedLengthMatch m)
{   
    s << "(";
    if (m.mtype_ == FixedLengthMatch::BASELIST)
    {
        for (size_t i = 0; i < m.bases_.size(); i++)
        {
            s << m.bases_[i];
        }
    }
    else if (m.mtype_ == FixedLengthMatch::COMPOSITION
        || m.mtype_ == FixedLengthMatch::DISJUNCTION)
    {
        for (size_t i = 0; i < m.conditions_.size(); i++)
        {
            s << "[" << m.conditions_[i] << "]" << "^" << m.length_;
            if (i < m.conditions_.size() - 1)
            {
                if (m.mtype_ == FixedLengthMatch::COMPOSITION)
                {
                    s << "&";
                }
                else if (m.mtype_ == FixedLengthMatch::DISJUNCTION)
                {
                    s << "|";
                }
            }
        }
    }
    else 
    {
        s << "N^" << m.length_;
    }
    s << ")";
    
    return s;
}
    

ostream& operator<<(ostream& s, const vec<FixedLengthMatch> v)
{
    for (size_t i = 0; i < v.size(); i++)
    {
        s << v[i];
    }
    
    return s;
}

void ComputeMotif(const String& motif, const vecbasevector& genome,
    const vecbitvector& amb, vecbitvector& in_motif,
    VecMatchedLocsVec& motif_locations)
{
    for (size_t i = 0; i < genome.size(); i++)
    {
        ComputeMotif(motif, genome[i], amb[i], in_motif[i], motif_locations[i]);
    }
    return;
}

void ComputeMotif(const String& motif, const BaseVec& chr,
    const BitVec& amb, BitVec& in_motif, MatchedLocsVec& motif_locations)
{
    // Build reverse complement of motif.
    
    vec<vec<FixedLengthMatch> > M(2);
    Parse(motif, M[0]);
    M[1] = M[0];
    M[1].ReverseMe();
    for (size_t i = 0; i < M[1].size(); i++)
    {
        M[1][i].ReverseComplement();
    }
    
    // Compute bases that are in motif.     
    int len = 0;
    for (size_t k = 0; k < M[0].size(); k++)
    {
        len += M[0][k].Length();
    }
    
    // return if the contig/chromosome is smaller than the pattern,
    // failure to do so will cause terrible crash or data corruption
    if (static_cast<unsigned int>(len) > chr.size())
    {
        return;
    }
    
    vec<char> G(chr.size());
    for (size_t j = 0; j < chr.size(); j++)
    {
        G[j] = as_base(chr[j]);
    }      

    // Precompute running tallies of base counts.
    vec<vec<vec<short> > > tally0(M[0].size());
    vec<vec<vec<short> > > tally1(M[0].size());
    
    for (size_t m = 0; m < M[0].size(); m++)
    {
        if (M[0][m].Mtype() == FixedLengthMatch::COMPOSITION || 
            M[0][m].Mtype() == FixedLengthMatch::DISJUNCTION)
        {
            tally0[m].resize(4);
            size_t len = M[0][m].Length();
            if (len <= chr.size())
            {
                for (size_t v = 0; v < 4; v++)
                {
                    tally0[m][v].resize(chr.size() - len + 1, 0);
                }

                for (size_t x = 0; x < len; x++)
                {
                    tally0[m][chr[x]][0]++;
                }
                for (size_t y = 1; y <= chr.size() - len; y++)
                { 
                    for (size_t v = 0; v < 4; v++)
                    {
                        tally0[m][v][y] = tally0[m][v][y - 1];
                    }
                    tally0[m][chr[y - 1]][y]--;
                    tally0[m][chr[y + len - 1]][y]++;
                }
                tally1[M[0].size() - m - 1] = tally0[m];   
            }
        }
    }
    
    // Find first ambiguous base.
    size_t next_ambig = 0;
    while (next_ambig < chr.size() && !amb[next_ambig])
    {
        next_ambig++;
    }
    
    // Go through the bases.
    for (size_t j = 0; j <= chr.size() - len; j++)
    {
        // Avoid ambiguous bases.
        if (next_ambig < j + len)
        {   
            if (next_ambig == j)
            {
                next_ambig++;
                while (next_ambig < chr.size() && !amb[next_ambig])
                {
                    next_ambig++;
                }
            }
        }
        else if (Matches(G, j, M[0], tally0) || Matches(G, j, M[1], tally1))
        {   
            motif_locations.push_back(SequenceLoc(j, j + len - 1));
            for (int k = 0; k < len; k++ )
            {    
                 in_motif.Set(j + k, True);
            }
        }
    }
    
    return;
}
