///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Dec 15, 2014 - <crdhelp@broadinstitute.org>
//

#ifndef PROFILEHMM_H_
#define PROFILEHMM_H_
#include <array>
#include <cmath>
#include "Vec.h"
#include "Basevector.h"
#include "paths/long/magic/Nanomer.h"
#include "VecUtilities.h"

class KLaplace {
public:
    double operator()( double x, double mean, double sd ) {
        return exp(-1.*std::abs(x-mean) / sd);
    }
};


class KmerProfileHMM {
public:
    struct Module {
        std::array<double,1024> emit;
        std::array<double,3> trans;
        enum { DEL, MATCH, INS };

        void dump( String const& head = "" ) const {
            cout << head << endl;
            for ( size_t i = 0; i < emit.size(); ++i )
                cout << i << ": " << emit[i] << endl;
        }

        String top_kmer() const {
            unsigned imax = std::max_element(emit.begin(), emit.end()) - emit.begin();
            return Nanomer<5>(imax).ToString();
        }

        String top_kmers(size_t n) const {
            auto sorted = emit;
            vec<int> index(1024,vec<int>::IDENTITY);
            ReverseSortSync(sorted, index);
            std::ostringstream s;
            for ( size_t i = 0; i < n; ++i ) {
                s << Nanomer<5>(index[i]).ToString();
                if ( i < n-1 ) s << ", ";
            }
            return s.str();
        }

        void el_mult(Module const& that) {
            for ( size_t i = 0; i < emit.size(); ++i )
                emit[i] *= that.emit[i];
        }

    };

    explicit KmerProfileHMM(int nmod) : modules(nmod) {};

    KmerProfileHMM( vec<double> const& levels,
            vec<double> const& model_mean,
            vec<double> const& model_stdev) {

        ForceAssertGt(levels.size(), 0u);
        modules.resize(levels.size());

        ForceAssertEq(model_mean.size(), modules[0].emit.size());
        ForceAssertEq(model_stdev.size(), modules[0].emit.size());

        KLaplace kernel;
        for ( size_t i = 0; i < levels.size(); ++i ) {
            long double sum = 0.;
            auto& module = modules[i];
            auto level = levels[i];
            for ( size_t j = 0; j < module.emit.size(); j++ ) {
                module.emit[j] = kernel( level, model_mean[j], model_stdev[j] );
                sum += module.emit[j];
            }
#if 1
            for ( size_t j = 0; j < module.emit.size(); j++ )
                module.emit[j] += ( sum / module.emit.size() * 0.000005 );
            for ( size_t j = 0; j < module.emit.size(); j++ )
                module.emit[j] /= 1.000005*sum;
#else
            for ( size_t j = 0; j < module.emit.size(); j++ )
                module.emit[j] /= sum;
#endif
            module.trans = {{0.008, 0.88, 0.112}};
        }
    }

    Module const& Mod(size_t i) const {
        ForceAssertLt(i, modules.size());
        return modules[i];
    }

    size_t size() const { return modules.size(); }

private:
    vec<Module> modules;

};



#endif /* PROFILEHMM_H_ */
