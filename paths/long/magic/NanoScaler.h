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

#ifndef NANOSCALER_H_
#define NANOSCALER_H_

class NanoScaler {
public:
    NanoScaler(double shift, double scale, double drift) :
        _shift(shift), _scale(scale), _drift(drift) {
    };

    vec<double> Scale(vec<double> const& raw, vec<double> const& times,
            Bool const reverse = False ) {
        ForceAssertEq(raw.size(), times.size());
        vec<double> out;
        out.reserve(raw.size());
        for ( size_t i = 0; i < raw.size(); ++i )
            out.push_back( _shift + _scale*raw[i] +
                    _drift*( times[i] - times[0]));

        if ( reverse ) std::reverse(out.begin(), out.end() );

        return out;
    }

    vec<double> ReverseAndScale( vec<double> const& raw, vec<double> const& times ) {
        return Scale(raw, times, True);
    }

private:
    double _shift;
    double _scale;
    double _drift;
};




#endif /* NANOSCALER_H_ */
