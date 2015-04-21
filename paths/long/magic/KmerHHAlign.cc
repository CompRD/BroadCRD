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
#include "paths/long/magic/KmerHHAlign.h"
#include "math/Array.h"
#include "system/System.h"


double KmerHHAlign::Saa( KmerProfileHMM::Module const& pmod,
        KmerProfileHMM::Module const& qmod )
{
    long double sum = 0.;

    for ( size_t j = 0; j < 1024; ++j ) {
//        cout << "qmod.emit[" << j << "]=" << qmod.emit[j];
//        cout << ", pmod.emit[" << j << "]=" << pmod.emit[j];
//        cout << ", bkgf[" << j << "]=" << bkgf[j] << endl;
        sum += qmod.emit[j]*pmod.emit[j]/bkgf[j];
    }

//    cout << "sum=" << sum << endl;

    return log(sum);
}

double KmerHHAlign::Align( KmerProfileHMM const& q, KmerProfileHMM const& p, align& a, bool debug )
{
    cout << "debug=" << ((debug)?"true":"false") << endl;

    // annoying q/p order to be consistent with paper
    size_t nq = q.size();
    size_t np = p.size();

    RecArray<double> Smm(nq,np);
    RecArray<double> Smi(nq,np);
    RecArray<double> Sim(nq,np);
    RecArray<double> Sdg(nq,np);
    RecArray<double> Sgd(nq,np);
    RecArray<uchar>  Trc(nq,np);

    // for now we're going to use these hard-coded limits AND ignoring
    // those in the module
    double qMM, pMM, qMI, pMI, qIM, pIM, qII, pII, \
                qMD, pMD, qDM, pDM, qDD, pDD;

//    double drate = 0.112;
//    double irate = 0.008;
    double drate = 0.333;
    double irate = 0.112;

    pMM = 1.0 - drate - irate; pMI = irate; pMD = drate;
    pDM = 1.0 - drate; pDD = drate;
    pIM = 1.0 - irate; pII = irate;

    qMM = pMM; qMI = pMI; qMD = pMD;
    qDM = pDM; qDD = pDD;
    qIM = pIM; qII = pII;

    auto neginf = std::numeric_limits<double>::lowest();

    for ( size_t i = 0; i < nq; ++i )  {
        Smm[i][0] = 0;
        Smi[i][0] = neginf; Sim[i][0] = neginf;
        Sgd[i][0] = neginf; Sdg[i][0] = neginf;
        Trc[i][0] = ' ';
    }

    for ( size_t j = 0; j < np; ++j ) {
        Smm[0][j] = 0;
        Sim[0][j] = neginf; Smi[0][j] = neginf;
        Sdg[0][j] = neginf; Sgd[0][j] = neginf;
        Trc[0][j] = ' ';
    }

    vec<double> tmp5(5);
    vec<double> tmp2(2);

    vec<double> log5(5);
    log5[0] = log( qMM * pMM );
    log5[1] = log( qMM * pIM );
    log5[2] = log( qIM * pMM );
    log5[3] = log( qDM * pMM );
    log5[4] = log( qMM * pDM );

    vec<double> log2A(2);
    log2A[0] = log( qMM * pMI );
    log2A[1] = log( qMM * pII );

    vec<double> log2B(2);
    log2B[0] = log( qMI * pMM );
    log2B[1] = log( qII * pMM );

    vec<double> log2C(2);
    log2C[0] = log(qMD);
    log2C[1] = log(qDD);

    vec<double> log2D(2);
    log2D[0] = log(pMD);
    log2D[1] = log(pDD);

    cout << Date() << ": starting profile HMM alignment" << endl;
    for ( size_t j = 1; j < np; ++j ) {
//        cout << "j = " << j << endl;
        for ( size_t i = 1; i < nq; ++i ) {
            if ( debug ) cout << "i,j=" << i << "," << j << endl;

            tmp5[0] = Smm[i-1][j-1] + log5[0]; // log( qMM + pMM );
            tmp5[1] = Smi[i-1][j-1] + log5[1]; // log( qMM + pIM );
            tmp5[2] = Sim[i-1][j-1] + log5[2]; // log( qIM + pMM );
            tmp5[3] = Sdg[i-1][j-1] + log5[3]; // log( qDM + pMM );
            tmp5[4] = Sgd[i-1][j-1] + log5[4]; // log( qMM + pDM );

            auto maxp = std::max_element(tmp5.begin(), tmp5.end());
            Trc[i][j] = maxp - tmp5.begin();

            if ( debug ) cout << "tmp5: " << printSeq(tmp5) << " [" << (int)Trc[i][j] << "]" << endl;

            double saa = Saa(q.Mod(i), p.Mod(j));
            Smm[i][j] = saa + *maxp;

            if ( debug ) {
                cout << "saa: " << saa << endl;
                cout << "Smm=" << Smm[i][j] << endl;
            }

            tmp2[0] = Smm[i-1][j] + log2A[0]; // log( qMM + pMI );
            tmp2[1] = Smi[i-1][j] + log2A[1]; // log( qMM + pII );
            Smi[i][j] = Max(tmp2);

            if (debug) cout << "Smi tmp2: " << printSeq(tmp2) << endl;

            tmp2[0] = Smm[i][j-1] + log2B[0]; // log( qMI + pMM );
            tmp2[1] = Sim[i][j-1] + log2B[1]; // log( qII + pMM );
            Sim[i][j] = Max(tmp2);

            if (debug) cout << "Sim tmp2: " << printSeq(tmp2) << endl;

            tmp2[0] = Smm[i-1][j] + log2C[0]; // log( qMD );
            tmp2[1] = Sdg[i-1][j] + log2C[1]; // log( qDD );
            Sdg[i][j] = Max(tmp2);

            if (debug) cout << "Sdg tmp2: " << printSeq(tmp2) << endl;

            tmp2[0] = Smm[i][j-1] + log2D[0]; // log( pMD );
            tmp2[1] = Sgd[i][j-1] + log2D[1]; // log( pDD );
            Sgd[i][j] = Max(tmp2);

            if (debug) cout << "Sgd tmp2: " << printSeq(tmp2) << endl;
        }
    }

    cout << "debugging Smm: " << endl;
    for ( size_t j = 0; j < 11; ++j ) {
        for ( size_t i = 0; i < 11; ++i ) {
            cout << Smm[i][j] << " ";
        }
        cout << endl;
    }

    cout << "debugging Sgd: " << endl;
    for ( size_t j = 0; j < 11; ++j ) {
        for ( size_t i = 0; i < 11; ++i ) {
            cout << Sgd[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Smm[X][j]=";
    for ( size_t j = 0; j < np; ++j )
        cout << Smm[nq-1][j] << " ";
    cout << endl << endl;

    cout << "Smm[i][X]=";
    for ( size_t i = 0; i < nq; ++i )
        cout << Smm[i][np-1] << " ";
    cout << endl;

    // decode the winning joint-path

    double max_over_j = Smm[nq-1][0];
    size_t max_index_j = 0;
    for (size_t j = 1; j < np; ++j )
        if ( Smm[nq-1][j] > max_over_j ) {
            max_index_j = j;
            max_over_j = Smm[nq-1][j];
        }

    double max_over_i = Smm[0][np-1];
    size_t max_index_i = 0;
    for ( size_t i = 1; i < nq; ++i )
        if ( Smm[i][np-1] > max_over_i ) {
            max_index_i = i;
            max_over_i = Smm[i][np-1];
        }

    PRINT2(max_over_i, max_over_j);

    size_t j,i;
    if ( max_over_i > max_over_j ) {
        i = max_index_i;
        j = np-1;
    } else {
        i = nq-1;
        j = max_index_j;
    }

    vec<pair<size_t,size_t>> path;
    vec<String> path_states;

    cout << "start i=" << i << ", j=" << j << endl;

    path.push_back(make_pair(i,j));
    path_states.push_back("MM");

    while ( i > 0 && j > 0 ) {
        // Smm Smi Sim Sdg Sgd
        if ( Trc[i][j] == 0 ) {
            path_states.push_back("MM");
            i--; j--;
        } else if ( Trc[i][j] == 1 ) {
            path_states.push_back("MI");
            i--;
        } else if ( Trc[i][j] == 2 ) {
            path_states.push_back("IM");
            j--;
        } else if ( Trc[i][j] == 3 ) {
            path_states.push_back("DG");
            i--;
        } else if ( Trc[i][j] == 4 ) {
            path_states.push_back("GD");
            j--;
        } else
            FatalErr("bad Trc="+ToString((int)Trc[i][j])
                    + " at " + ToString(i) + "," + ToString(j) );
        path.push_back(make_pair(i,j));

        PRINT3(path_states.back(), path.back().first, path.back().second);
    }

    cout << "final i=" << i << ", j=" << j << endl;
    for ( size_t k = 0; k < path.size(); ++k ) {
        int this_i = path[k].first;
        int this_j = path[k].second;
        double dot = Saa(q.Mod(this_i), p.Mod(this_j));
//        q.Mod(this_i).dump("Q_i:");
//        p.Mod(this_j).dump("P_j:");
        auto tmp = q.Mod(this_i);
        tmp.el_mult(p.Mod(this_j));
        cout << path_states[k] << " (" << path[k].first << "," <<
                path[k].second << ")"  << ", dot=" << dot
                << ", top_q_i=" << q.Mod(this_i).top_kmers(3)
                << ", top_p_j=" << p.Mod(this_j).top_kmers(3)
                << ", top_joint=" << tmp.top_kmers(3) << endl;

    }

#if 0
    size_t nstates =
    RecArray<double> dscore;
    RecArray<size_t> dindex;
#endif

    return 0.;
}


