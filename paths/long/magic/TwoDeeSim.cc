///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "paths/long/magic/NanoData.h"
#include "IntPairVec.h"
#include "lookup/LookAlign.h"
#include "paths/long/CreateGenome.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "kmers/KmerRecord.h"
#include "paths/RemodelGapTools.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "PrintAlignment.h"
#include "paths/long/ultra/MultipleAligner.h"
#include "paths/long/ultra/FounderAlignment.h"
#include "paths/long/magic/KmerProfileHMM.h"
#include "paths/long/magic/NanoScaler.h"
#include "paths/long/magic/KmerHHAlign.h"
#include "paths/long/magic/NanoRCModel.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "PackAlign.h"
#include "feudal/BaseVec.h"
#include "random/NormalRandom.h"

namespace {

void MakeModel( double lo, double hi, double std_mean, double std_std,
        size_t nelem, vec<double>& means, vec<double>& stds )
{
    double segm = (hi - lo) / nelem;
    NormalRandom gaussian(std_mean, std_std);

    for ( size_t i = 0; i < nelem; ++i ) {
        double val = lo + segm*i + segm/2.0;
        means.push_back(val);
        stds.push_back(gaussian.value());
    }
}

void SimulateNano( BaseVec const& genome, double irate, double drate,
        vec<double> const& level_mean, vec<double> const& level_std,
        vec<double>& signal, vec<String>& truth )
{
    unsigned int kindex = 0u;
    for (size_t k = 0; k < 5; ++k ) { // first K-1 bases
        kindex <<= 2;
        kindex &= 0x3ff;
        kindex |= genome[k];
    }
    for ( size_t k=4; k < genome.size(); ) {
        BaseVec kmer(genome.begin() + k - 4, genome.begin() + k + 1);
        unsigned toss = randomx();
        ForceAssertLt(kindex, 1024u);
        double mean = level_mean[kindex];
        double std = level_std[kindex];
        NormalRandom gaussian(mean, std);
        if ( toss < irate*RNGen::RNGEN_RAND_MAX ) {
            // insertion
            signal.push_back( gaussian.value() );
            truth.push_back("I"+ToString(mean));
            continue;   // NON-STRUCTURED CODE!
        } else if ( toss < (drate+irate)*RNGen::RNGEN_RAND_MAX ) {
            // deletion
            ++k;
            truth.push_back("D");
        } else {
            // match/sub
            signal.push_back( gaussian.value() );
            ++k;
            truth.push_back("M"+ToString(mean));
        }
        kindex <<= 2;
        kindex &= 0x3ff;
        kindex |= genome[k];
    }
}

};

int main(int argc, char* argv[]) {
    RunTime();

    BeginCommandArguments;
    CommandArgument_String_Doc(NANO, ".nano dataset (if READ specified)");
    CommandArgument_String_OrDefault_Doc(INGENOME, "", "read genome");
    CommandArgument_String_OrDefault_Doc(OUTGENOME, "", "read genome");
    CommandArgument_Bool_OrDefault_Doc(RC, False, "simulate RC genome");
    CommandArgument_Int_OrDefault_Doc(READ, -1,
            "read to use for MODEL (or simulate if not specified)");
    EndCommandArguments;

    //
    vec<double> m0_mean, m0_std, m1_mean, m1_std;
    if ( READ != -1 ) {
        cout << Date() << ": reading nano file " << NANO << endl;
        NanoData nano;
        BinaryReader::readFile(NANO, &nano);
        ForceAssertGt((int)nano.GetSize(), READ);
        cout << Date() << ": read " << nano.GetSize() << " reads" << endl;

        auto& m0 = nano.GetModels(NanoData::TEMP)[READ];
        auto m1_fw = nano.GetModels(NanoData::COMP)[READ];
        auto m1 = NanoRCModel( m1_fw ).Get();

        m0_mean = m0["level_mean"]->get_vec_double();
        m0_std = m0["level_stdv"]->get_vec_double();
        m1_mean = m1["level_mean"]->get_vec_double();
        m1_std = m1["level_stdv"]->get_vec_double();
    } else {
        MakeModel( 0., 1024., 0.005, 0.001, 1024, m0_mean, m0_std );
#if 0
        for ( size_t i = 0; i < m0_mean.size(); ++i ) {
            cout << "m0_mean[" << i << "]=" << m0_mean[i] <<
                    ", m0_std[" << i << "]=" << m0_std[i] << endl;
        }
#endif
        m1_mean = m0_mean;
        m1_std = m0_std;
    }

    ///
    /// manage the GENOME
    ///

    size_t genomelen = 100u;
    vecbasevector genome;
    if ( INGENOME != "" ) {
        genome.ReadAll(INGENOME);
        if ( genome.size() > 1 )
            cout << "warning: only using first component as genome" << endl;
        else if ( genome.size() == 0 )
            FatalErr("nothing read from " + INGENOME);
        genomelen = genome[0].size();
    } else {
        genome.push_back( basevector(genomelen) );
        for ( size_t i = 0; i < genome[0].size(); ++i )
            genome[0].Set(i, randomx() % 4 );
    }

    if ( OUTGENOME != "" ) {
        genome.WriteAll(OUTGENOME);
    }



    vec<double> e0, e1;
    vec<String> true0, true1;
    double drate = 0.112;
    double irate = 0.008;

    size_t start1 = randomx() % static_cast<long int>(0.05 * genomelen);
//    size_t start1 = 0u;
    size_t end_offset1 = randomx() % static_cast<long int>(0.05 * genomelen);
//    size_t end_offset1 = 0u;
    size_t end1 = genomelen - end_offset1 - 1;

    cout << "start1=" << start1 << ", end1=" << end1 << ", end_offset1="
            << end_offset1 << endl;

    // e0 is perfect, e1 has whatever
    auto g = genome[0];
    if ( RC ) g.ReverseComplement();

    SimulateNano( g, 0.0, 0.0,
            m0_mean,
            m0_std,
            e0, true0);

    SimulateNano( BaseVec(genome[0].begin() + start1,
            genome[0].begin() + end1), irate, drate,
            m1_mean,
            m1_std,
            e1, true1 );

    cout << "True0: " << endl;
    Bool done = False;
    for ( size_t j = 0; !done; ++j ) {
        done = True;
        for ( size_t i = 0; i < true0.size(); ++i ) {
            if ( true0[i].size() > j ) { cout << true0[i][j]; done = False; }
            else cout << " ";
        }
        cout << endl;
    }

    cout << "True1: " << endl;
    done = False;
    for ( size_t j = 0; !done; ++j ) {
        done = True;
        for ( size_t i = 0; i < true1.size(); ++i ) {
            if ( true1[i].size() > j ) { cout << true1[i][j]; done = False; }
            else cout << " ";
        }
        cout << endl;
    }

    if ( RC ) {
        cout << "REVERSE COMPLEMENTING MODEL & EVENTS" << endl;
        cout << "test reverse e0; forward=" << endl;
        for ( size_t i = 0; i < 5; ++i )
            cout << e0[i] << " ";
        cout << endl;
        std::reverse(e0.begin(), e0.end());
        cout << "reverse=" << endl;
        for ( size_t i = 0; i < 5; ++i )
            cout << e0[e0.size()-i-1] << " ";
        cout << endl;
        auto temp = m0_mean;
        NanoRCModel::testSwap(m0_mean);
        NanoRCModel::testSwap(m0_std);
        cout << "RC model:" << endl;
        for ( size_t i = 0; i < 1024; ++i )
            cout << i << ": " << temp[i] << "<->" << m0_mean[i] << endl;
    }

    KmerProfileHMM p( e0,
                      m0_mean,
                      m0_std);

#if 0
    for ( size_t i = 0; i < e0.size(); ++i ) {
        cout << "e0[" << i << "]=" << e0[i] << endl;
        for ( size_t j = 0; j < p.Mod(i).emit.size(); ++j )
            cout << "       j=" << j << ", =>" <<
            p.Mod(i).emit[j] << endl;
    }
#endif


    KmerProfileHMM q( e1,
                      m1_mean,
                      m1_std);

    KmerHHAlign aligner;
    align a;
    auto score = aligner.Align(q, p, a);
}
