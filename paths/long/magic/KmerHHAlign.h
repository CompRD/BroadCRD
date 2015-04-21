///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Dec 18, 2014 - <crdhelp@broadinstitute.org>
//

// implementation of Soding's Profile HMM <-> Profile HMM alignment

#ifndef KMERHHALIGN_H_
#define KMERHHALIGN_H_
#include "paths/long/magic/KmerProfileHMM.h"
#include "PackAlign.h"


class KmerHHAlign {
public:
    typedef array<double,1024> Profile;
    KmerHHAlign() : bkgf(uniform()) {};
    KmerHHAlign(Profile const& back) : bkgf(back) {};

    double Align(KmerProfileHMM const& q, KmerProfileHMM const& p, align& a, bool debug = false);

    double Saa( KmerProfileHMM::Module const& pmod,
            KmerProfileHMM::Module const& qmod);

    // make a uniform distribution for initialization
    static Profile uniform() {
        Profile back;
        // uniform background distribution
        double len = back.size();
        for ( auto& bkg : back ) bkg = 1. / len;
        return back;
    }

private:
    Profile bkgf;       // background frequency

};



#endif /* KMERHHALIGN_H_ */
