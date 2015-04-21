///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Jan 5, 2015 - <crdhelp@broadinstitute.org>
//

#ifndef NANORCMODEL_H_
#define NANORCMODEL_H_
#include "paths/long/magic/SimpleDataFrame.h"
#include "system/Assert.h"

class NanoRCModel {
public:
    NanoRCModel(SimpleDataFrame const& fw) : rc(fw) {
        ForceAssertEq(fw.size(),1024u);
        vec<Bool> done(1024, false);
        for ( size_t row = 0; row < 1024u; ++row ) {
            if ( done[row] ) continue;

            unsigned int irc = 0u;
            unsigned int ifw = row;
            for ( size_t i = 0; i < 5; ++i ) {
                irc = ( irc << 2 ) | (ifw & 0x3);
                ifw >>= 2;
            }
            irc = (~irc & 0x3ff);
            ifw = row;

            // swap the fw and rc row and mark done
            rc.SwapRows(irc, ifw);
            done[irc] = done[ifw] = True;
        }
    };

    static void testSwap(vec<double>& a) {
        ForceAssertEq(a.size(),1024u);
        vec<Bool> done(1024, false);
        for ( size_t row = 0; row < 1024u; ++row ) {
            if ( done[row] ) continue;

            unsigned int irc = 0u;
            unsigned int ifw = row;
            for ( size_t i = 0; i < 5; ++i ) {
                irc = ( irc << 2 ) | (ifw & 0x3);
                ifw >>= 2;
            }
            irc = (~irc & 0x3ff);
            ifw = row;

            // swap the fw and rc row and mark done
            std::swap(a[irc], a[ifw]);
            done[irc] = done[ifw] = True;
        }
    }

    SimpleDataFrame Get() { return rc; }


private:
    SimpleDataFrame rc;
};




#endif /* NANORCMODEL_H_ */
