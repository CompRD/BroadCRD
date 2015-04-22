///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Nov 19, 2014 - <crdhelp@broadinstitute.org>
//

#ifndef NANODATA_H_
#define NANODATA_H_

#include <array>
#include "paths/long/magic/SimpleDataFrame.h"
#include "Basevector.h"
#include "Qualvector.h"

class NanoData {
public:
    NanoData() {
        temp_indicator.fill(false);
    };

    typedef enum { TEMP, COMP, CONS, SIZEOF } ReadType;
    void AddRead(
            basevector const& tbase, qualvector const& tqual,
            basevector const& cbase, qualvector const& cqual,
            basevector const& base, qualvector const& qual ) {
        bases[TEMP].push_back(tbase);
        quals[TEMP].push_back(tqual);
        bases[COMP].push_back(cbase);
        quals[COMP].push_back(cqual);
        bases[CONS].push_back(base);
        quals[CONS].push_back(qual);
    }

    size_t GetSize() { return bases[CONS].size(); }

    void AddRead() {
        ForceAssertEq(static_cast<size_t>(SIZEOF), temp_indicator.size());
        for ( size_t i = 0; i < SIZEOF; ++i )
            if ( !temp_indicator[i] ) FatalErr( "didn't set all components of a read");
        for ( size_t i = 0; i < SIZEOF; ++i ) {
            bases[i].push_back(temp_bases[i]);
            quals[i].push_back(temp_quals[i]);
            temp_indicator[i] = false;
        }
    }

    void AddPartial( ReadType const type,
            basevector const& base,
            qualvector const& qual ) {
        temp_bases[type] = base;
        temp_quals[type] = qual;
        temp_indicator[type] = true;
    }

    void AddModel( SimpleDataFrame const& c, ReadType const type ) {
        if ( type == TEMP ) tmodels.push_back(c);
        else if ( type == COMP ) cmodels.push_back(c);
        else FatalErr("unknown model ReadType");
    }

    vec<SimpleDataFrame>& GetModels(ReadType const type) {
        if ( type == TEMP ) return tmodels;
        else if ( type != COMP ) FatalErr("unknown model ReadType");
        return cmodels;
    }

    void AddEvents( SimpleDataFrame const& e, ReadType const type ) {
        if ( type == TEMP ) tevents.push_back(e);
        else if ( type == COMP ) cevents.push_back(e);
        else FatalErr("unknown event ReadType");
    }

    vec<SimpleDataFrame>& GetEvents(ReadType const type) {
        if ( type == TEMP ) return tevents;
        else if ( type != COMP ) FatalErr("unknown model ReadType");
        return cevents;
    }

    void AddAligns( SimpleDataFrame const& a ) { aligns.push_back(a); }

    void AddHairpin( SimpleDataFrame const& a ) { hairpin.push_back(a); }

    void writeBinary( BinaryWriter& writer ) const {
        for ( size_t i = 0; i < SIZEOF; ++i ) {
            writer.write(bases[i]);
            writer.write(quals[i]);
        }
        writer.write(tmodels);
        writer.write(cmodels);
        writer.write(tevents);
        writer.write(cevents);
        writer.write(aligns);
        writer.write(hairpin);
    }
    void readBinary( BinaryReader& reader ) {
        for ( size_t i = 0; i < SIZEOF; ++i ) {
            reader.read(&bases[i]);
            reader.read(&quals[i]);
        }
        reader.read(&tmodels);
        reader.read(&cmodels);
        reader.read(&tevents);
        reader.read(&cevents);
        reader.read(&aligns);
        reader.read(&hairpin);
    }
    static size_t externalSizeof() { return 0u; }

    void CopyFrom(NanoData const& src, size_t idx) {
        for ( size_t i = 0; i < SIZEOF; ++i ) {
            bases[i].push_back( src.bases[i][idx] );
            quals[i].push_back( src.quals[i][idx] );
        }
        tmodels.push_back(src.tmodels[idx]);
        cmodels.push_back(src.cmodels[idx]);
        tevents.push_back(src.tevents[idx]);
        cevents.push_back(src.cevents[idx]);
        aligns.push_back(src.aligns[idx]);
        hairpin.push_back(src.hairpin[idx]);
    }

    vecbasevector& GetBases(ReadType const type) { return bases[type]; }
    vecqualvector& GetQuals(ReadType const type) { return quals[type]; }
    vec<SimpleDataFrame>& GetTempModels() { return tmodels; }
    vec<SimpleDataFrame>& GetCompModels() { return cmodels; }
    vec<SimpleDataFrame>& GetTempEvents() { return tevents; }
    vec<SimpleDataFrame>& GetCompEvents() { return cevents; }
    vec<SimpleDataFrame>& GetAligns() { return aligns; }
    vec<SimpleDataFrame>& GetHairpin() { return hairpin; }


private:
    std::array<vecbasevector, SIZEOF> bases;
    std::array<vecqualvector, SIZEOF> quals;

    vec<SimpleDataFrame> tmodels;
    vec<SimpleDataFrame> cmodels;

    vec<SimpleDataFrame> tevents;
    vec<SimpleDataFrame> cevents;

    vec<SimpleDataFrame> aligns;

    vec<SimpleDataFrame> hairpin;

    std::array<basevector, SIZEOF> temp_bases;
    std::array<qualvector, SIZEOF> temp_quals;
    std::array<bool, SIZEOF>       temp_indicator;
};
SELF_SERIALIZABLE(NanoData);



#endif /* NANODATA_H_ */
