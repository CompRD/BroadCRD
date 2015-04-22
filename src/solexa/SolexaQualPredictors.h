/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/** Calculate predictors for quality for each base

 \class SolexaQualPredictors

 Unified predictor calculation here for use by both QualsFromSolexa
 and ConstructPhredTableSolexa.

 Currently hardcoded at 5 parameters, could be made more flexible.

 */

#ifndef SOLEXA_QUAL_PREDICTORS_H
#define SOLEXA_QUAL_PREDICTORS_H

#include "Basevector.h"
#include "Qualvector.h"
#include "solexa/FourBase.h"
#include "solexa/SolexaTools.h"
#include "feudal/VirtualMasterVec.h"

class SolexaQualPredictors
{
private:
    const int POSSTART;
    const int RANGE;
    const int MINQ10I;
    const int LOCALI;
    const int O1; ///first other, local by default.
    const int O2; ///second other, none by default.
    int currentRead;
    VirtualMasterVec<qvec> oldq;
    qvec curOldQ;
    VirtualMasterVec<bvec> reads;
    bvec curRead;
    VirtualMasterVec<FourBaseVec> I;
    FourBaseVec curI;
    VirtualMasterVec<FourBaseVec> I2;
    FourBaseVec curI2;
    vec<float> hScore;
    vec<float> values;

public:
    /// \param HEAD: the root for all files
    /// \param predParamHandler: to provide POS_BASE and RANGE.
    SolexaQualPredictors( const String & HEAD,
                          const SolexaPredictorParameters * handler )
    : POSSTART(handler->Get("posstart")),
      RANGE(handler->Get("localrange")),
      MINQ10I(handler->Get("minq10i")),
      LOCALI(handler->Get("locali")),
      O1(handler->Get("o1")),
      O2(handler->Get("o2")),
      currentRead(-1),
      oldq( (HEAD+".qualb").c_str() ),
      reads( (HEAD+".fastb").c_str() ),
      I( (HEAD+".intensities").c_str() ),
      I2(IsRegularFile(HEAD+".intensities2")?
              (HEAD+".intensities2").c_str() :
              (HEAD+".intensities").c_str()),
      hScore(),
      values(N())
    { LoadRead(0); }

    /// Number of reads
    int size() const { return reads.size(); }

    bool empty() const { return reads.empty(); }

    /// length of current read.
    int readSize() const { return curRead.size(); }

    int N() const { return (0 == O2) ? 5 : 6; }

    /// Make this the current read, and cache relevant information
    void LoadRead( int i )
    {
        if ( currentRead == i )
            return;

        currentRead = i;
        curOldQ = oldq[i];
        curRead = reads[i];
        curI = I[i];
        curI2 = I2[i];

        // Calculate the max homopolymer of me and my adjacent bases.
        vec<int> h(curRead.size());
        for ( unsigned int b = 0; b != curRead.size(); ++b )
            h[b] = curRead.Homopol(b);
        hScore.resize(h.size());
        MaxWindow(h.begin(), h.end(), hScore.begin(), 1);

        //set the read quality value from _sig or _sig2.
        FourBaseVec & v = (1 == MINQ10I) ? curI : curI2;
        //PRINT3(v.size(), N(), values.size());
        values[1] = -MinQuality(v, 10);
    }

    /// Get the vector of predictor values.
    /// Note that this is a mutable member vector: this is needed because it
    /// is eventually passed to CapValues() in PhredTableReader, which may need
    /// to modify it. It's OK if that happens, because we reset all but
    /// values[1] each time anyway, and the cap for values[1] would be valid
    /// for the whole read.
    vec<float> & Values( int read, int b )
    {
        LoadRead(read);
        values[0] = abs(POSSTART - b);
        values[2] = -curOldQ[b];
        values[4] = hScore[b];
        switch ( O1 )
        {
        case 1:
            switch ( LOCALI )
            {
            case 1:
                values[3] = -MinQuality(curI, b + RANGE + 1, b - RANGE);
                break;
            case 2:
                values[3] = -MinQuality(curI2, b + RANGE + 1, b - RANGE);
                break;
            default:
                FatalErr("Bad LOCALI " << LOCALI)
                ;
            }
            break;
        case 2:
            values[3] = Height(b);
            break;
        default:
            FatalErr("Bad O1 " << O1)
            ;
            break;
        }
        switch ( O2 )
        {
        case 0:
            break;
        case 1:
            values[5] = Height(b);
            break;
        case 2:
            values[5] = PhasingCorrection(b);
            break;
        default:
            FatalErr("Bad O2 " << O2)
            ;
            break;
        }

        return values;
    }

private:
    float Height( int b )
    {
        float ret;
        switch ( LOCALI )
        {
        case 1:
            ret = curI[b].MaxInt();
            break;
        case 2:
            ret = curI2[b].MaxInt();
            break;
        default:
            FatalErr("Bad LOCALI " << LOCALI)
            ;
            break;
        }
        return ret;
    }

    float PhasingCorrection( int b )
    {
        int call = curI2[b].Call();
        float diff = abs(curI[b].base(call) - curI2[b].base(call));
        float sum = curI[b].base(call) + curI2[b].base(call);
        return diff / sum;
    }

};

#endif // SOLEXA_QUAL_PREDICTORS_H
