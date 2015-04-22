///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Inspect Kmer Parcels.

#include <climits>
#include <strstream>
#include <math.h>

#include "Basevector.h"

#include "MainTools.h"

#include "feudal/QualNibbleVec.h"  // QualNibbleVec

#include "kmers/KmerParcels.h" // KmerParcelReader, KmerBatch
#include "kmers/KmerParcelsTools.h"
#include "system/WorklistN.h"

static inline 
String Tag(String S = "KPI") { return Date() + " (" + S + "): "; } 


template<class T>
class LockingVec : public vec<T>
{
private:
  LockedData _locks[256];

public:
  LockingVec(size_t i, size_t v) : vec<T>(i, v) {}

  LockedData & Lock(size_t i) 
  {
    return _locks[i % 256];
  }
};


template<class T>
class Bins : public vec<T>
{
public:
  T Sum() const 
  {
    T sum = 0;
    for (typename vec<T>::const_iterator it = (*this).begin();
         it != (*this).end(); ++it)
      sum += (*it);
    return sum;
  }

  void Print(std::ostream& out = cout) const
  {
    out << "# 1:y 2:y/total 3:cum(y) 4:cum(y)/total" << endl;
    
    T total = Sum();
    T cum = 0;
    for (typename vec<T>::const_iterator it = (*this).begin();
         it != (*this).end(); ++it) {
      cum += (*it);
      out << " " << (*it)
          << " " << setprecision(10) << static_cast<double>(*it) / total
          << " " << cum
          << " " << setprecision(10) << static_cast<double>(cum) / total
          << endl;
    }
  }

};






class BinnedHits
{
private:
  static const size_t _n_super_bins = 19;

  vec< Bins<size_t> > _bins;
  LockedData _locks[_n_super_bins];

public:
  BinnedHits() 
  { 
    _bins.resize(_n_super_bins);
    for (size_t i = 0; i != _n_super_bins; i++) 
      _bins[i].resize(NumBins(i), 0);
  }
  
  ~BinnedHits() {}

  size_t SuperBinID(const size_t i) const 
  {
    if (i ==         1) return  0;
    if (i ==         2) return  1;
    if (i ==         3) return  2;
    if (i < (1L <<  4)) return  3;
    if (i < (1L <<  6)) return  4;
    if (i < (1L <<  8)) return  5;
    if (i < (1L << 10)) return  6;
    if (i < (1L << 12)) return  7;
    if (i < (1L << 14)) return  8;
    if (i < (1L << 16)) return  9;
    if (i < (1L << 18)) return 10;
    if (i < (1L << 20)) return 11;
    if (i < (1L << 22)) return 12;
    if (i < (1L << 24)) return 13;
    if (i < (1L << 26)) return 14;
    if (i < (1L << 28)) return 15;
    if (i < (1L << 30)) return 16;
    if (i < (1L << 32)) return 17;
    return 18;
  }

  size_t NumBins(const size_t s_bin_ID) const 
  { 
    if (s_bin_ID == 0) return 1;
    if (s_bin_ID == 1) return 2;
    if (s_bin_ID == 2) return 2;
    if (s_bin_ID == 3) return 8;
    if (s_bin_ID == 4) return 32;
    return 50;
  }
    

  void AddHit(const size_t n, const size_t m) 
  {
    size_t s_bin_ID = SuperBinID(n);
    double mu = 0.5 * n;           // binomial: n.p
    double sig = sqrt(0.5 * mu);   // binomial: sqrt(n.p.(1-p))
    double im = static_cast<size_t>(fabs(m - mu) / sig);

    size_t binID = (im < _bins[s_bin_ID].size()) ? im : _bins[s_bin_ID].size() - 1;
 
    Locker lock(_locks[s_bin_ID]);
    _bins[s_bin_ID][binID]++;
  }

  void Write(const String & filename) const
  {
    ofstream ofs;
    ofs.open(filename.c_str());
    
    for (size_t i = 0; i != _n_super_bins; i++) {
      _bins[i].Print(ofs);
      ofs << endl << endl;
    }
    ofs.close();
  }


};







class BatchTest
{
private:
  const size_t _min_kf;
  const size_t _max_kf;
  const size_t _min_kf_FWRC;
  const size_t _max_kf_FWRC;
  const double _min_n_sig_FWRC;
  const double _max_n_sig_FWRC;
  const size_t _read_len;
  const double _min_n_sig_BEG;
  const double _max_n_sig_BEG;

public:
  BatchTest(const size_t min_kf,
            const size_t max_kf,
            const size_t min_kf_FWRC,
            const size_t max_kf_FWRC,
            const double min_n_sig_FWRC,
            const double max_n_sig_FWRC,
            const size_t read_len,
            const double min_n_sig_BEG,
            const double max_n_sig_BEG)
    :_min_kf(min_kf),
     _max_kf(max_kf),
     _min_kf_FWRC(min_kf_FWRC),
     _max_kf_FWRC(max_kf_FWRC),
     _min_n_sig_FWRC(min_n_sig_FWRC),
     _max_n_sig_FWRC(max_n_sig_FWRC),
     _read_len(read_len),
     _min_n_sig_BEG(min_n_sig_BEG),
     _max_n_sig_BEG(max_n_sig_BEG) {}


  bool operator () (const KmerBatch & batch) const
  {
    const double kf   = batch.GetKmerFreq();
    
    if (kf < _min_kf) return false;
    if (kf > _max_kf) return false;

    const double kfRC = batch.GetKmerFreqRC();

    const size_t low_kf = (kfRC < 0.5 * kf) ? kfRC : kf - kfRC; 
    if (low_kf < _min_kf_FWRC) return false;
    if (low_kf > _max_kf_FWRC) return false;


    // test asymmetry of FW and RC frequencies
    if (_min_n_sig_FWRC > 0 || _max_n_sig_FWRC > 0) {

      const double n_sig_FWRC = NumSigmaBinomial(kf, kfRC);
    
      //cout << "_min_n_sig_FWRC = " << _min_n_sig_FWRC << endl;
      //cout << "_max_n_sig_FWRC = " << _max_n_sig_FWRC << endl;
      //cout << "n_sig_FWRC      = " << n_sig_FWRC << endl;
    
      if (_min_n_sig_FWRC > 0 && n_sig_FWRC < _min_n_sig_FWRC) return false;
      if (_max_n_sig_FWRC > 0 && n_sig_FWRC < _max_n_sig_FWRC) return false;
    }

    // test asymmetry of positions in the first and second half of the read
    if (_read_len >= batch.GetK() && 
        (_min_n_sig_BEG > 0 || _max_n_sig_BEG > 0)) {

      const size_t nk = _read_len - batch.GetK() + 1;
      const double prob = static_cast<double>(nk / 2) / static_cast<double>(nk);

      const double kf_beg = batch.GetKmerFreqPosLT(nk / 2);
      const double n_sig_BEG = NumSigmaBinomial(kf, kf_beg, prob); 

      /*
      cout << "kf_beg = " << kf_beg << endl;
      cout << "kf_beg_mu = " << kf_beg_mu << endl;
      cout << "kf_beg_sig = " << kf_beg_sig << endl;
      cout << "n_sig_BEG = " << n_sig_BEG << endl;
      */

      if (_min_n_sig_BEG > 0 && n_sig_BEG < _min_n_sig_BEG) return false;
      if (_max_n_sig_BEG > 0 && n_sig_BEG < _max_n_sig_BEG) return false;
    }

    return true;
  }


};




void QualsRead(const String & quals_fn, 
               VecQualNibbleVec * quals)
{
  quals->clear();

  typedef VirtualMasterVec<qualvector> FvvQv;

  cout << Tag() << "Loading quality scores from " << quals_fn << endl;
  FvvQv quals_full(quals_fn.c_str());
  for (FvvQv::const_iterator qItr = quals_full.begin();
       qItr != quals_full.end(); qItr++)
    quals->push_back(QualNibbleVec(*qItr));
}





template<class TEST>
class ParcelProc
{
private:
  const KmerParcelsDiskStore & _parcels;
  const KmerParcelsDiskStore & _parcels_cmp;

  const bool _do_and;
  const bool _do_not;


  const TEST & _batch_test;
  
  const bool                _do_stats;
  KmerFrequencyStatistics * _p_stats_yes;
  KmerFrequencyStatistics * _p_stats_no;

  BinnedHits              * _p_bhits_yes;
  BinnedHits              * _p_bhits_no;

  const bool             _do_reads;
  LockingVec<uint32_t> * _reads_n_hits;

  const bool             _do_quals;

  const size_t _read_len;
  const bool   _print_all;
  NaiveThreadPool & _thread_pool;

public:
  
  ParcelProc(const KmerParcelsDiskStore & parcels, 
             const KmerParcelsDiskStore & parcels_cmp, 
             const bool do_and,
             const bool do_not,
             const TEST & batch_test,
             const bool do_stats,
             KmerFrequencyStatistics * p_stats_yes,
             KmerFrequencyStatistics * p_stats_no,
             BinnedHits * p_bhits_yes,
             BinnedHits * p_bhits_no,
             const bool do_reads,
             LockingVec<uint32_t> * reads_n_hits,
             const bool do_quals,
             const size_t read_len,
             const bool print_all,
             NaiveThreadPool & thread_pool)
    : _parcels(parcels),
      _parcels_cmp(parcels_cmp),
      _do_and(do_and),
      _do_not(do_not),
      _batch_test(batch_test),
      _do_stats(do_stats),
      _p_stats_yes(p_stats_yes),
      _p_stats_no(p_stats_no),
      _p_bhits_yes(p_bhits_yes),
      _p_bhits_no(p_bhits_no),
      _do_reads(do_reads),
      _reads_n_hits(reads_n_hits),
      _do_quals(do_quals),
      _read_len(read_len),
      _print_all(print_all),
      _thread_pool(thread_pool)
  {}


  void operator() (const size_t parcel_ID) 
  {
    const size_t thread_ID = _thread_pool.AssignThread(parcel_ID);

    String tag = "FPP";
    //if (_verbose)
    cout << Tag(tag) 
         << "Processing parcel " << parcel_ID 
         << " on thread " << thread_ID << endl;

    KmerParcelReader parcel_reader(_parcels, parcel_ID);

    KmerParcelReader parcel_cmp_reader;

    if (_do_not || _do_and) 
      parcel_cmp_reader.Init(_parcels_cmp, parcel_ID);

    VecQualNibbleVec kmer_quals;
    if (_do_quals)
      QualsRead(_parcels.GetDirectoryName() + "/" + ToString(parcel_ID) + ".qualb", &kmer_quals);


    typedef vec<KmerLoc>::const_iterator KLIter;

    size_t kmer_ID = 0;
    while (parcel_reader.GetNextKmerBatch()) {
      const KmerBatch & batch = parcel_reader.CurrentKmerBatch();

      if (true) {
      }
      
      // ---- kmer matches all kmer criteria
      if (_batch_test(batch)) {

        const bool found_match = (_do_and || _do_not) ? 
          FindMatchingKmerBatch(parcel_cmp_reader, batch.Kmer()) : false;
        
        
        if ((!_do_and ||  found_match) &&
            (!_do_not || !found_match)) { 
          
          if (_do_stats) {
            _p_stats_yes->UpdateFreqs(thread_ID, batch.GetKmerFreq());
            
            //_p_bhits_yes->AddHit(batch.GetKmerFreq(), batch.GetKmerFreqRC());
          }
          
          for (KLIter it = batch.begin(); it != batch.end(); it++) {
            size_t i = (*it).GetReadID();
            Locker lock((*_reads_n_hits).Lock(i));
            (*_reads_n_hits)[i]++;
          }
          
          
          if (!_do_stats && !_do_reads) {
            cout << "   ";

            if (_do_quals) batch.Print(cout, kmer_quals[kmer_ID], 
                                       (_read_len) ? _read_len : 101);
            else           batch.Print(cout, (_read_len) ? _read_len : 101);


            cout << " Parcel " << parcel_ID
                 << " Batch " << parcel_reader.GetKmerBatchID() 
                 << endl;
          }
          
        }
        else { // kmer doesn't match kmer criteria from other parcels
              
          if (_print_all && !_do_stats && !_do_reads) {
            cout << "** ";

            if (_do_quals) batch.Print(cout, kmer_quals[kmer_ID], 
                                       (_read_len) ? _read_len : 101);
            else           batch.Print(cout, (_read_len) ? _read_len : 101);

            cout << " Parcel " << parcel_ID
                 << " Batch " << parcel_reader.GetKmerBatchID() 
                 << endl;
          }
        }
      }
      // ---- kmer does NOT match any criteria
      else {
        
        if (_do_stats) {
          _p_stats_no->UpdateFreqs(thread_ID, batch.GetKmerFreq());
          
          //_p_bhits_no->AddHit(batch.GetKmerFreq(), batch.GetKmerFreqRC());
        }
      }

      kmer_ID++;
    } 
    
    _thread_pool.UnassignThread(parcel_ID);
  }


};












int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;

  CommandArgument_UnsignedInt_Doc(K, "kmer size");
  CommandArgument_String_Doc(HEAD, "looks for <HEAD>.<K>merParcels directory");
  CommandArgument_Int_OrDefault_Doc(PARCEL_ID, -1, "the parcel index (-1: all)");
  CommandArgument_Int_OrDefault_Doc(BATCH_ID, -1, "the batch index (-1: all)");
  
  CommandArgument_Bool_OrDefault_Doc(DO_READS, false, "write all matching and non-matching reads.");
  CommandArgument_Bool_OrDefault_Doc(DO_STATS, false, "write statistics on matching and non-matching kmers.");
  CommandArgument_Bool_OrDefault_Doc(DO_QUALS, false, "display kmers with bases colored by average qual.");

  CommandArgument_String_OrDefault_Doc(HEAD_NOT, "", "exclude kmers from <HEAD_NOT>.<K>merParcels");
  CommandArgument_String_OrDefault_Doc(HEAD_AND, "", "exclude kmers not in <HEAD_AND>.<K>merParcels");
  
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_FREQ, 0, "minimum kmer frequency.");
  CommandArgument_UnsignedInt_OrDefault_Doc(MAX_FREQ, (INT_MAX), "maximum kmer frequency.");

  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_FREQ_FWRC, 0, "minimum FW or RC kmer frequency.");
  CommandArgument_UnsignedInt_OrDefault_Doc(MAX_FREQ_FWRC, (INT_MAX), "maximum FW or RC kmer frequency.");

  CommandArgument_Double_OrDefault_Doc(MIN_N_SIG_FWRC, -1.0, "number of minimum standard deviations in FW-RC.");
  CommandArgument_Double_OrDefault_Doc(MAX_N_SIG_FWRC, -1.0, "number of maximum standard deviations in FW-RC.");

  CommandArgument_Double_OrDefault_Doc(MIN_N_SIG_BEG, -1.0, "number of minimum standard deviations in BEG-END.");
  CommandArgument_Double_OrDefault_Doc(MAX_N_SIG_BEG, -1.0, "number of maximum standard deviations in BEG-END.");
  CommandArgument_UnsignedInt_OrDefault_Doc(READ_LEN, 101, "read_size.");

  CommandArgument_Bool_OrDefault_Doc(PRINT_ALL, false, "print all kmers that match the kmer based criteria");

  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 4, "number of parallel threads to use.");
  EndCommandArguments;

  
  cout << Tag() << "BEGIN" << endl;

  bool do_and = (HEAD_AND != "");
  bool do_not = (HEAD_NOT != "");

  if (do_and && do_not) {
    cout << "**** You can't specify both HEAD_AND and HEAD_NOT." << endl;
    ::exit(1);
  }

  KmerParcelsDiskStore parcels(K, HEAD);
  String parcels_dn = parcels.GetDirectoryName();

  string head_cmp = (do_and) ? HEAD_AND : HEAD_NOT;

  KmerParcelsDiskStore parcels_cmp(K, head_cmp);

  size_t n_parcels = parcels.GetNumParcels();
  size_t parcel0_ID = (PARCEL_ID >= 0) ? PARCEL_ID : 0;
  size_t parcel1_ID = (PARCEL_ID >= 0) ? PARCEL_ID+1 : n_parcels;

  cout << Tag() << "n_parcels = " << n_parcels << endl;
  cout << Tag() << "parcel0_ID = " << parcel0_ID << endl;
  cout << Tag() << "parcel1_ID = " << parcel1_ID << endl;


  if (!DO_READS && !DO_STATS) NUM_THREADS = 1; // stdout only 

  String label = "";
  if (DO_READS || DO_STATS) {
    if (do_and) label += "and_" + head_cmp;
    if (do_not) label += "not_" + head_cmp;
    
    if (MIN_N_SIG_FWRC > 0) {
      if (label != "") label += ".";
      char buf[10];
      sprintf(buf, "minnsigRC%.1f", MIN_N_SIG_FWRC);
      label += buf;
    }
  }
  else {
    label = "match";
  }
 
  size_t n_reads = parcels.GetNumReads();
  LockingVec<uint32_t> reads_n_hits(n_reads, 0);

  BaseVecVec reads;
  if (DO_READS) reads.ReadAll(HEAD + ".fastb");

  KmerFrequencyStatistics kmer_stats_yes(n_reads, NUM_THREADS);
  KmerFrequencyStatistics kmer_stats_no(n_reads, NUM_THREADS);

  BinnedHits bhits_yes, bhits_no;


  cout << Tag() << "n_reads = " << n_reads << endl;


  NaiveThreadPool thread_pool(NUM_THREADS, n_parcels);

  
  BatchTest batch_test(MIN_FREQ, MAX_FREQ,
                       MIN_FREQ_FWRC, MAX_FREQ_FWRC,
                       MIN_N_SIG_FWRC, MAX_N_SIG_FWRC,
                       READ_LEN, MIN_N_SIG_BEG, MAX_N_SIG_BEG);


  ParcelProc<BatchTest> proc(parcels, parcels_cmp,
                             do_and, do_not,
                             batch_test,
                             DO_STATS,
                             &kmer_stats_yes,
                             &kmer_stats_no,
                             &bhits_yes,
                             &bhits_no,
                             DO_READS,
                             &reads_n_hits,
                             DO_QUALS,
                             READ_LEN,
                             PRINT_ALL,
                             thread_pool);
  {    
    parallelFor(parcel0_ID,parcel1_ID,proc,NUM_THREADS);
  }
  
  if (DO_STATS) {
    cout << Tag() << "writing kmer stats." << endl;

    kmer_stats_yes.AllKmerFrequencyCounters().Write
      (parcels_dn + "/kmer_frequencies." + label + ".count.dat");
    kmer_stats_no.AllKmerFrequencyCounters().Write
      (parcels_dn + "/kmer_frequencies.no." + label + ".count.dat");

    bhits_yes.Write(parcels_dn + "/binned_hits." + label + ".dat");
    bhits_no.Write(parcels_dn + "/binned_hits.no." + label + ".dat");
  }

  if (DO_READS) {
    cout << Tag() << "writing matching and non matching reads." << endl;

    BaseVecVec yes_reads, no_reads;
    for (size_t i = 0; i != n_reads; i++) {
      if (reads_n_hits[i] > 0)
        yes_reads.push_back(reads[i]);
      else
        no_reads.push_back(reads[i]);
    }

    yes_reads.WriteAll(HEAD + "." + label + ".fastb");
    no_reads.WriteAll(HEAD + ".no_" + label + ".fastb");
  }
  
  cout << Tag() << "writing number of reads with a specific number of hits." << endl;
  MapOfCounters counts(reads_n_hits);
  counts.Write(parcels_dn + "/reads_n_hits." + label + ".count.dat");


  size_t n_yes = 0;
  size_t n_no = 0;
  for (size_t i = 0; i != n_reads; i++) {
    if (reads_n_hits[i] > 0)
      n_yes++;
    else
      n_no++;
  }

  cout << Tag() << "n_reads                     = " << setw(12) << n_reads 
       << " [" << setw(5) << fixed << setprecision(1) << 100.0  << " %]" << endl;
  cout << Tag() << "n_reads with matching kmers = " << setw(12) << n_yes
       << " [" << setw(5) << fixed << setprecision(1) 
       << 100.0 * static_cast<double>(n_yes)/ static_cast<double>(n_reads) << " %]" << endl;
  cout << Tag() << "n_reads w/o  matching kmers = " << setw(12) << n_no 
       << " [" << setw(5) << fixed << setprecision(1) 
       << 100.0 * static_cast<double>(n_no)/ static_cast<double>(n_reads) << " %]" << endl;
                
  cout << Tag() << "END" << endl;

}

