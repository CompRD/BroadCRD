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
#include "paths/long/magic/KmerHHAlign.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "Alignment.h"

namespace {

class NanoScaler {
public:
    NanoScaler(double shift, double scale, double drift) :
        _shift(shift), _scale(scale), _drift(drift) {};

    vec<double> Scale(vec<double> const& raw, vec<double> const& times) {
        cout << "Scale() called" << endl;
        PRINT3(_shift,_scale,_drift);
        PRINT2(raw.size(), times.size());
        ForceAssertEq(raw.size(), times.size());
        vec<double> out;
        out.reserve(raw.size());
        for ( size_t i = 0; i < raw.size(); ++i )
            out.push_back( _shift + _scale*raw[i] +
                    _drift*( times[i] - times[0]));

        return out;
    }

private:
    double _shift;
    double _scale;
    double _drift;
};

template<int K> void MakeKmerLookup0SingleX(const vecbasevector& unibases,
        vec<triple<kmer<K>, int, int> >& kmers_plus) {
    vec<int64_t> starts;
    starts.push_back(0);
    for (size_t i = 0; i < unibases.size(); i++) {
        const basevector& u = unibases[i];
        starts.push_back(starts.back() + Max(0, u.isize() - K + 1));
    }
    kmers_plus.resize(starts.back());
    for (size_t i = 0; i < unibases.size(); i++) {
        const basevector& u = unibases[i];
        for (int j = 0; j <= u.isize() - K; j++) {
            int64_t r = starts[i] + j;
            kmers_plus[r].first.SetToSubOf(u, j);
            kmers_plus[r].second = i;
            kmers_plus[r].third = j;
        }
    }
    Sort(kmers_plus);
}


matrix<double>  rates_to_confusion_matrix( double p_d, double p_i, double p_s )
{
     // rates to confusion matrix p(D|T):
     // p_good = 1.0 - p_ind - p_del - p_sub
     // p( base_1 | base_2 ) = p_good      \forall base_1 == base_2
     // p( _ | _ ) = p_good
     // p( base_1 | base_2 ) = p_sub/3.    \forall base_1 != base_2
     // p( base | _ ) = p_ins
     // p( _ | base ) = p_del
     matrix<double> cprob(5,5,0.);

     double p_good = 1.0 - p_i - p_s - p_d;

     double bias = 0.001;       // ensure that nothing has zero probability

     for ( size_t i = 0; i < 5; ++i )
         for ( size_t j = 0; j < 5; ++j ) {
             if ( i == j ) cprob[i][j]=p_good+bias;
             else if ( i == 4 ) cprob[i][j] = p_d+bias;
             else if ( j == 4 ) cprob[i][j] = p_i+bias;
             else cprob[i][j] = p_s / 3.0 + bias;
         }

     // ensure that columns sum to one
     for (size_t j = 0; j < 5; ++j ) {
          double sum = 0.;
          for ( size_t i = 0; i < 5; ++i )
               sum += cprob[i][j];
          for ( size_t i = 0; i < 5; ++i )
               cprob[i][j] /= sum;
     }


     return cprob;
}

void calc_error_rates( VecUCharVec const& multi, UCharVec const& cons,
                        vec<triple<double,double,double>>& rates  )
{
     size_t Nj = multi.size();
     size_t Ni = multi[0].size();

     rates.resize(Nj);
     for ( size_t j = 0; j < Nj; ++j ) {
          // del-ins-sub
          int del = 0;
          int ins = 0;
          int sub = 0;
          double count = 0.;
          for ( size_t i = 0; i < Ni; ++i ) {
               auto call = multi[j][i];
               if ( call != 5 && cons[i] != 5 ) {
                    count+=1.;
                    if ( call < 4 && cons[i] < 4 && call != cons[i] )
                         sub++;
                    else if ( call < 4 && cons[i] == 4 )
                         ins++;
                    else if ( cons[i] < 4 && call == 4 )
                         del++;
               }
          }

          if ( count > 0 ) rates[j] = make_triple( del / count, ins / count, sub / count );
          else rates[j] = make_triple(0.,0.,0.);
     }
}

void bayes_consensus(vec<matrix<double>> const& p_d_t, VecUCharVec const& multi,
          UCharVec& cons)
// note cons is both input and output -- it's expected that an esimtated
// consensus is here simply for calculation of the prior
{
     // consensus:
     //     Pr{ T_i = t | D; \theta } = p(t) \prod_j p(D_j | T=t )
     // ----------------------------------------------------------
     // \sum_t' Pr{ T_i = t' | D; \theta } = p(t') \prod_j p(D_j | T=t' )
     size_t Nj = multi.size();
     size_t Ni = multi[0].size();

     vec<double> priors(5,0.);
     // global prior based on the current estimated consensus
     for ( size_t i = 0; i < Ni; ++i )
          if (cons[i] < 5) priors[cons[i]] += 1.;
     double sum = Sum(priors);
     ForceAssertGe( sum, 1.0 );
     for ( auto& prior : priors ) prior = log(prior/sum);


     for ( size_t i = 0; i < Ni; ++i ) {

          vec<double> calls = priors;
          for ( size_t j = 0; j < Nj; ++j ) {
               for ( size_t k = 0; k < 5; ++k )
                    if ( multi[j][i] < 5 )
                         calls[k] += log( p_d_t[j][ multi[j][i] ][k] );
          }
          vec<int> ident(calls.size(), vec<int>::IDENTITY);
          if ( i == 0 ) {
              /*
              cout << "COL0: ";
              for ( size_t k = 0; k < calls.size(); ++k )
                  cout << exp(calls[k]) <<  " ";
              cout << endl;
              */
          }
          ReverseSortSync(calls, ident );
          if ( i == 0 ) {
              /*
              cout << "sorted COL0: ";
              for ( size_t k = 0; k < calls.size(); ++k )
                  cout << exp(calls[k]) <<  " ";
              cout << endl;
              cout << "sorted log COL0: ";
              for ( size_t k = 0; k < calls.size(); ++k )
                  cout << calls[k] <<  " ";
              cout << endl;
              cout << "sorted ident COL0: ";
              for ( size_t k = 0; k < calls.size(); ++k )
                  cout << ident[k] <<  " ";
              cout << endl;
              */
          }
          cons[i] = ident[0];
     }
}

// begs to be in a class...
typedef enum { CONSENSUS_NONE, CONSENSUS_MAJORITY,
     CONSENSUS_STAT, CONSENSUS_ITER } ConsensusType;

ConsensusType parseConsensusType( String const& s )
{
     if ( s == "MAJORITY" ) return CONSENSUS_MAJORITY;
     else if ( s == "STAT" ) return CONSENSUS_STAT;
     else if ( s == "ITER" ) return CONSENSUS_ITER;
     else if ( s != "NONE" )
          FatalErr("unrecognized consensus type: " + s);

     return CONSENSUS_NONE;

}

String consensusTypeStr( ConsensusType const& ctype )  {
     if ( ctype == CONSENSUS_MAJORITY) return "MAJORITY";
     else if ( ctype == CONSENSUS_STAT ) return "STAT";
     else if ( ctype == CONSENSUS_ITER ) return "ITER";
     else if ( ctype != CONSENSUS_NONE)
          FatalErr("unrecognized consensus code: " + ctype);

     return "NONE";
}

void Consensus( double p_d, double p_i, double p_s, VecUCharVec const& multi,
          size_t founder, vec<triple<double,double,double>>& rates, UCharVec& cons,
           ConsensusType mode = CONSENSUS_MAJORITY, vec<double>* pWeights = 0)
{
     if ( multi.size() == 0 ) {
          cons.clear();
          rates.clear();
          return;
     }
     ForceAssertLt(founder, multi.size() );

     size_t Nj = multi.size();
     size_t Ni = multi[0].size();
     rates.resize( Nj );
     cons.resize( Ni );

     // initialization value for consensus
     // majority vote
     for ( size_t i = 0; i < Ni; ++i )  {
          vec<double> score(6,0.);
          vec<int> ident(6, vec<int>::IDENTITY);
          for ( size_t j = 0; j < Nj; ++j ) {
               ForceAssertLt( multi[j][i], 6 );
               if ( pWeights == 0 ) score[ multi[j][i] ] += 1.0;
               else score[ multi[j][i] ] += (*pWeights)[j];
          }
          ReverseSortSync(score, ident );
          // don't score missing bases, unless they're all missing bases
          // 5 == missing, 4 == gap, <4 == base
          double top_score = score[0];
          int top_base = ident[0];
          size_t top_idx = 0;
          if ( top_base == 5 && score[1] != 0 ) {
               top_score = score[1];
               top_base = ident[1];
               top_idx = 1;
          }
          cons[i] = top_base;           // set provisionally to first high scorer...
          // ... but we make any ties go to the founder base
          for ( size_t k = top_idx; k < score.size() && score[k] == top_score; ++k )
               if ( ident[k] == multi[founder][i] )
                    cons[i] = ident[k];
     }

     // initialization values for rates
     calc_error_rates( multi, cons, rates );

     if ( mode == CONSENSUS_MAJORITY ) return;  // EARLY RETURN

     for ( size_t iter = 0; iter < 10; ++iter ) {
     vec<matrix<double>> p_d_t( Nj );
     for ( size_t j = 0; j < Nj; ++j ) {
         p_d_t[j] = rates_to_confusion_matrix(rates[j].first, rates[j].second, rates[j].third);
         // debugging below
         /*
         PRINT4(j, rates[j].first, rates[j].second, rates[j].third );
         cout << "matrix " << j << ":" << endl;
         for ( int k1 = 0; k1 < p_d_t[j].Nrows(); ++k1 ) {
             for ( int k2 = 0; k2 < p_d_t[j].Ncols(); ++k2 )
                 cout << p_d_t[j][k1][k2] << " ";
             cout << endl;
         }
         */

     }

     bayes_consensus( p_d_t, multi, cons );

     if ( mode == CONSENSUS_STAT ) return; //EARLY RETURN
     calc_error_rates( multi, cons, rates );
     }
}

void write_uchar_reads( ostream& out, String const& label, UCharVec read )
{
    size_t count = 0;
    out << ">" << label << endl;
    for (size_t i = 0; i < read.size(); ++i ) {
        if ( read[i] < 4 ) {
            out << as_base( read[i] );
            if ( ++count % 80 == 0 ) out << endl;
        }
    }
    if ( count % 80 != 0 ) out << endl;
}


void PrintMultipleAlignment(vecbasevector const& gang, vec<align> const& aligns)
{
    cout << "\nmultiple alignment\n\n";
    const double del_rate = 0.05;
    const double ins_rate = 0.02;
    const double sub_rate = 0.05;
    Scorer scorer(sub_rate, del_rate, ins_rate);
    VecUCharVec multi;
    const int bandwidth = 200;
    AlignFriendsToFounder(gang, 0, aligns, scorer, &multi);

    // Delete columns having only one base, and at least five entries.

    PRINT3(gang.size(), gang[0].size(), aligns.size());

    vec<Bool> to_delete(multi[0].size(), False);
    for (int j = 0; j < (int) multi[0].size(); j++) {
        int non45 = 0;
        for (int i = 0; i < (int) multi.size(); i++)
            if (multi[i][j] != 4 && multi[i][j] != 5)
                non45++;
        if (non45 == 1 && multi.size() >= 5)
            to_delete[j] = True;
    }
    for (int i = 0; i < (int) multi.size(); i++) {
        SerfVec<uchar> m = multi[i];
        vec<uchar> mx;
        for (int l = 0; l < (int) m.size(); l++)
            mx.push_back(m[l]);
        EraseIf(mx, to_delete);
        multi[i].resize(0);
        for (int l = 0; l < mx.isize(); l++)
            multi[i].push_back(mx[l]);
    }

    // weights are for majority voting (CONSENSUS_MAJORITY) and
     // for initialization of other methods (which are initialized
     // with majority vote).  A weight of 2, will score as if the
     // read appeared twice in the stack.

    UCharVec cons1,cons2;
    vec<triple<double,double,double>> rates;
    vec<double> weights( multi.size(), 1.0);
    Consensus( del_rate, ins_rate, sub_rate, multi, 0, rates, cons1,
            CONSENSUS_MAJORITY, &weights );
    Consensus( del_rate, ins_rate, sub_rate, multi, 0, rates, cons2,
            CONSENSUS_STAT, &weights );
    multi.push_back( cons1 );
    multi.push_back( cons2 );

    Ofstream( debug_out, "neil.fasta" );
    write_uchar_reads(debug_out, "founder", multi[0] );
    write_uchar_reads(debug_out, "majority", cons1);
    write_uchar_reads(debug_out, "stat", cons2);

    const int width = 80;
    for (int js = 0; js < (int) multi[0].size(); js += width) {
        for (int i = 0; i < (int) multi.size(); i++) {
            Bool nothing = True;
            for (int j = js; j < Min((int) multi[0].size(), js + width); j++) {
                uchar c = multi[i][j];
                if (c != 5)
                    nothing = False;
            }
            if (!nothing) {
                for (int j = js; j < Min((int) multi[0].size(), js + width);
                        j++) {
                    uchar c = multi[i][j];
                    if (c < 4)
                        cout << as_base(c);
                    if (c == 4)
                        cout << "-";
                    if (c == 5)
                        cout << "=";
                }
                cout << "\n";
            }
        }
        cout << "\n";
    }

}



};

int main(int argc, char* argv[]) {
    RunTime();

    BeginCommandArguments;
    CommandArgument_String_Doc(NANO, ".nano dataset");
    CommandArgument_Int_Doc(READ, "read to correct");
    CommandArgument_Int_OrDefault_Doc(MINHITS, 5,
            "minimum number of 12-mer matches between reads");
    CommandArgument_Double_OrDefault_Doc(ERRTOL,0.80,
            "dump reads with score > ERRTOL*overlap_length");
    CommandArgument_Int_OrDefault_Doc(MIN_OVERLAP_BASES,300,
            "don't include aligns with overlap < MIN_OVERLAP_BASES")
    CommandArgument_Bool_OrDefault_Doc(VISUAL, False,
            "print SmithWatBandedA visual alignment");
    CommandArgument_Bool_OrDefault_Doc(MULTIPLE, False,
            "print multiple visual alignment");
    CommandArgument_Bool_OrDefault_Doc(DEBUG, False, "debug alignment");
    CommandArgument_String_OrDefault_Doc(CHECK_HEAD,"", "HEAD to checkpoint to (for later RESTART)");
    CommandArgument_String_OrDefault_Doc(RESTART_HEAD,"", "start by loading a file from CHECKPOINT");

    EndCommandArguments;

    ForceAssertGe(READ,0);

    if ( CHECK_HEAD != "" && CHECK_HEAD == RESTART_HEAD ) {
        FatalErr("CHECK_HEAD should not equal RESTART_HEAD (and a weird thing to do, anyway");
    }

    cout << Date() << ": reading nano file " << NANO << endl;
    NanoData nano;
    vec<int> read_ids;
    vec<align> aligns;
    vecbasevector gang;

    if ( RESTART_HEAD == "" ) {
        BinaryReader::readFile(NANO, &nano);
        // TODO: DEBUG
//        cout << "DUMPING one" << endl;
//        nano.GetModels(NanoData::TEMP)[0].Dump();
//        for ( size_t i = 0; i < nano.GetSize(); ++i )
//            PRINT(nano.GetModels(NanoData::TEMP)[i]("shift")->get_double());


        ForceAssertGt((int)nano.GetSize(), READ);
        cout << Date() << ": read " << nano.GetSize() << " reads" << endl;

#if 0
        if (READ != 0)
            std::swap(nano.GetBases(NanoData::CONS)[0], nano.GetBases(NanoData::CONS)[READ]);
#endif

        int const K1 = 16;
        int const N = nano.GetSize();

        cout << Date( ) << ": building kmers" << endl;
        vec< triple<kmer<K1>,int,int> > kmers_plus;
        auto bases = nano.GetBases(NanoData::CONS);
#if 0
        bases.reserve( 2*bases.size() );
        size_t count = bases.size();
        cout << Date() << ": raw reads: " << bases.size() << endl;
        for ( size_t i = 0; i < count; ++i ) {
            bases.push_back(bases[i]);
            bases.back().ReverseComplement();
        }
        cout << Date() << ": reads+rc: " << bases.size() << endl;
#else
        size_t count = bases.size();
        cout << Date() << ": raw reads: " << bases.size() << endl;
#endif
        MakeKmerLookup0( bases, kmers_plus );

        // do a one-to-many alignment using consensuses
        cout << Date() << ": building matches" << endl;
        vec<triple<int, int, int> > match;
        vec<kmer<K1> > km;
        vec<int> kmult;
        vec<int> pos1s;
        for (int i = 0; i < kmers_plus.isize(); i++) {
            int j;
            for (j = i + 1; j < kmers_plus.isize(); j++)
                if (kmers_plus[j].first != kmers_plus[i].first)
                    break;
            for (int k1 = i; k1 < j; k1++)
                for (int k2 = i; k2 < j; k2++) {
                    if (k2 == k1)
                        continue;
                    int id1 = kmers_plus[k1].second, id2 = kmers_plus[k2].second;
                    if (id1 >= N)
                        continue;
                    int pos1 = kmers_plus[k1].third, pos2 = kmers_plus[k2].third;
                    int offset = pos1 - pos2;
                    match.push(id1, id2, offset);
                    km.push_back(kmers_plus[i].first);
                    kmult.push_back(j - i);
                    pos1s.push_back(pos1);
                }
            i = j - 1;
        }
        // Sort(match);
        SortSync(match, km, kmult, pos1s);

        cout << Date() << ": matches for read " << READ << endl;

        read_ids.push_back(READ);
        gang.push_back( bases[READ]);
        aligns.push_back( align() );

        for ( int i = 0; i < match.isize();  i++) {
            if ( match[i].first != READ ) continue;
            int j;
            for ( j = i+1; j < match.isize(); ++j )
                if ( match[j].second != match[i].second ) break;
            if ( j-i > MINHITS ) {
                int id1 = match[i].first;
                int id2 = match[i].second;
                vec<int> offsets;
                for ( int k = i; k < j; ++k ) offsets.push_back(match[k].third);

                // okay, align them!
                int offset_low = Min(offsets);
                int offset_high = Max(offsets);
                int offset = (offset_low + offset_high) / 2;
                int bandwidth = (offset_high - offset_low) / 2;
                int errsx;
                int BW_ADD = 300;
                const int MISMATCH_PENALTY = 2;
                const int GAP_PENALTY = 3;
                align a;
                SmithWatBandedA( bases[id2], bases[id1], -offset, bandwidth + BW_ADD, a, errsx,
                        0, MISMATCH_PENALTY, GAP_PENALTY );
                auto overlap = a.EndOnTarget() - a.StartOnTarget();
                auto errfrac = errsx/static_cast<float>(overlap);
                if ( overlap >= MIN_OVERLAP_BASES && errfrac <= ERRTOL ) {
                    auto sz1 = bases[id1].size();
                    auto sz2 = bases[id2].size();
                    cout << "id1=" << match[i].first << ", id2=" << match[i].second <<
                            ", nhits=" << j-i << ": ";
                    cout << printSeq(offsets) << endl;
                    cout << "size1=" << sz1 << ", size2=" << sz2 <<
                            ", overlap(on 1)=" << overlap <<
                            ", errors=" << errsx << ", frac_overlap=" <<
                            overlap/static_cast<float>(sz1) <<
                            ", frac_errors=" << errfrac << endl;
                    if (VISUAL) PrintVisualAlignment( True, cout, bases[id2], bases[id1], a );
                    read_ids.push_back(id2);
                    cout << "before push_back gang.size()=" << gang.size() << endl;
                    gang.push_back(bases[id2]);
                    cout << "after push_back gang.size()=" << gang.size() << endl;
                    a.Flip();       // necessary for multiple alignment
                    aligns.push_back(a);
                }
            }
            i=j-1;          // i=j after loop increment
        }
    } else {
        BinaryReader::readFile(RESTART_HEAD+".nano", &nano);
        BinaryReader::readFile(RESTART_HEAD+".read_ids", &read_ids);
        BinaryReader::readFile(RESTART_HEAD+".gang", &gang);
        BinaryReader::readFile(RESTART_HEAD+".aligns", &aligns);
    }

    if (MULTIPLE) { PrintMultipleAlignment(gang, aligns); }

    if ( CHECK_HEAD != "" ) {
        NanoData tmp_nano;
        tmp_nano.CopyFrom(nano, read_ids[0]);
        tmp_nano.CopyFrom(nano, read_ids[1]);
        PRINT2(nano.GetModels(NanoData::TEMP)[read_ids[0]]("shift")->get_double(),
               tmp_nano.GetModels(NanoData::TEMP)[0]("shift")->get_double() );
        BinaryWriter::writeFile(CHECK_HEAD+".nano",tmp_nano);
        BinaryWriter::writeFile(CHECK_HEAD+".read_ids", read_ids);
        BinaryWriter::writeFile(CHECK_HEAD+".gang", gang);
        BinaryWriter::writeFile(CHECK_HEAD+".aligns", aligns);
    }

    // Instantiate profile HMM
    int id0 = RESTART_HEAD == "" ? read_ids[0] : 0;
    cout << "on restore or continue..." << endl;
    PRINT(nano.GetModels(NanoData::TEMP)[id0]("shift")->get_double());
    int id1 = RESTART_HEAD == "" ? read_ids[1] : 1;
    ForceAssert(RESTART_HEAD != "" || id0 == READ);
    auto& e0_raw = nano.GetEvents(NanoData::TEMP)[id0];
    auto& m0 = nano.GetModels(NanoData::TEMP)[id0];
    auto e0_bases = nano.GetBases(NanoData::TEMP)[id0];

    PRINT2(id0,id1);

    cout << "SHIFT: " << m0("shift")->get_double() << endl;

    auto e0 = NanoScaler(
            m0("shift")->get_double(),
            m0("scale")->get_double(),
            m0("drift")->get_double() ).Scale(
                    e0_raw["mean"]->get_vec_double(),
                    e0_raw["start"]->get_vec_double());

    auto& e1_raw = nano.GetEvents(NanoData::TEMP)[id1];
    auto& m1 = nano.GetModels(NanoData::TEMP)[id1];
    auto e1_bases = nano.GetBases(NanoData::TEMP)[id1];

    alignment al;
    auto errs = SmithWatAffine(e0_bases, e1_bases, al, true, true);
    cout << "Errors=" << errs << endl;
    cout << "alignment=" << al << endl;
    PrintVisualAlignment(False, cout, e0_bases, e1_bases, al );

    Ofstream(debug_out, "neil2.fasta");
    e0_bases.Print(debug_out, "e0");
    e1_bases.Print(debug_out, "e1");

    auto e1 = NanoScaler(
            m1("shift")->get_double(),
            m1("scale")->get_double(),
            m1("drift")->get_double() ).Scale(
                    e1_raw["mean"]->get_vec_double(),
                    e1_raw["start"]->get_vec_double() );

    KmerProfileHMM p( e0,
                      m0["level_mean"]->get_vec_double(),
                      m0["level_stdv"]->get_vec_double() );
    KmerProfileHMM q( e1,
                      m1["level_mean"]->get_vec_double(),
                      m1["level_stdv"]->get_vec_double() );

    KmerHHAlign aligner;
    align a;
    auto score = aligner.Align(p,q,a,DEBUG);
}
