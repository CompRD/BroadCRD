///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// TheProgram. A new program.
//

// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "math/Functions.h"
#include "paths/long/ultra/ThreadedBlocks.h"
#include "paths/long/ultra/ConsensusScoreModel.h"

double Combination( int N, int m ) {
    double result = 1;
    for( int k = 0; k < m; k++ ) {
        result *= N - k;
        result /=  k+1;
    }
    return result;
}

int main(int argc, char *argv[])
{
    RunTime( );
    // essential arguments
    BeginCommandArguments;
    CommandArgument_Int_OrDefault_Doc(NREADS, 50, "Number of reads");
    CommandArgument_Int_OrDefault_Doc(TLEN, 24, "Len of the true read");
    CommandArgument_Bool_OrDefault_Doc(SIM_READ_LEN, False, "Simulate the read length");
    CommandArgument_Double_OrDefault_Doc(DEL, 0.03, "Deletion rate");
    CommandArgument_Double_OrDefault_Doc(INS, 0.002, "Insertion rate");
    CommandArgument_Double_OrDefault_Doc(SUB, 0.008, "Substitution rate");
    // logging options
    CommandArgument_Int_OrDefault(VERBOSITY, 0);
    EndCommandArguments;
    // Check command line arguments
    // Thread control, etc.
    double clock = WallClockTime( );
    
    double del_rate = DEL, ins_rate = INS, sub_rate = SUB;
    ConsensusScoreModel model( del_rate, ins_rate, sub_rate, false, false );

    vec<BaseVec> reads( NREADS );
    int len = TLEN;
    if ( SIM_READ_LEN ) {
        int sum  = 0;
        srand(-1);
        int rand_size = 100000000;
        for ( size_t i = 0; i < reads.size(); ++i ) {
            int ndel = 0;
            for( int j = 0; j < len; j++ ) {
                double r = double( rand() % rand_size ) / rand_size;
                if ( r < del_rate ) ndel++;
                else if ( r > 1 - ins_rate ) ndel--;
            }
            reads[i] = BaseVec( String(len - ndel, 'A') );
            sum += len - ndel;
        }
        double mean_len = double( sum ) / double( reads.size() );
        cout << "Simulated mean read length " << mean_len << endl;
    }
    else {
        int rlen = len * ( 1 - del_rate );
        for ( size_t i = 0; i < reads.size(); ++i ) 
            reads[i] = BaseVec( String(rlen, 'A') );
        cout << "Use read length " << rlen << " for all reads" << endl;
    }

    if ( VERBOSITY >= 1 )
        for ( size_t i = 0; i < reads.size(); ++i ) 
            cout << reads[i].ToString() << " " << ToString( reads[i].size() ) << endl;

    vec< int > lens;
    vec< double > probs0;
    vec< double > probs1;
    for( int x = len - 5; x <= len + 5; x++ ) {
        //cout << "x= " << x << endl;
        BaseVec concensus( String(x, 'A') );
        double score0 = 0;
        double score1 = 0;
        for ( size_t i = 0; i < reads.size(); ++i ) {
            int N = concensus.size();
            int M = reads[i].size();
            int score = model.Score( concensus, reads[i] );
            score0 += score;
            //// DEBUG
            //double s = Min(M,N) *log(1 - ins_rate - sub_rate - del_rate) ;
            //if ( M < N ) s += log( del_rate ) * (N-M);
            //if ( N < M ) s += log( ins_rate ) * (M-N);

            // corrections
            // There are C( N, abs(M-N) ) ways of generating the length M from a length N homopolymer
            // hence p(N->M) = p0(N->M) * Correction
            // log p(N-M) = log( p0(N->M) ) + log(Correction)
            double correction = Combination( N, abs(N-M) );
            score1 += score - log(correction) * 100;
            //cout << score << " " << -s * 100 << " " << correction << "(" <<  score - log(correction) * 100  << ")" << endl;
        }
        lens.push( x );
        probs0.push( exp( -score0/100 ) );
        probs1.push( exp( -score1/100 ) );
    }
    double prob_max0 = Max( probs0 );
    double prob_max1 = Max( probs1 );
    cout << "len\tprob0   \tprob1(re-weight):" << endl;
    cout << setprecision( 3 );
    for ( size_t i = 0; i < lens.size(); ++i ) {
        String mark = ( lens[i] == TLEN ? "*": "" );
        cout << scientific << lens[i] << mark << "\t" << probs0[i] / prob_max0 << "\t(" << probs1[i]/prob_max1 << ")" << endl;
    }
    cout << "\n" << Date( ) << ": time used = " << TimeSince(clock) << endl;
    cout << Date() << ": TheProgram done!" << endl;    
}
