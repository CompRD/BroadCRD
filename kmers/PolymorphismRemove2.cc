///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

//
//    Finds SNPs and Indels and removes them.  work in progress.
//    Also, correlate them!
//
//  DEPRECATED after correlation implementation


#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "Map.h"
#include "graph/GraphAlgorithms.h"

#include "kmers/naif_kmer/KmerFreqAffixesMap.h"

#include "paths/PolymorphismRemoveCore.h"


static inline 
String Tag(String S = "PR") { return Date() + " (" + S + "): "; } 








class NPhases_t 
{
  uint64_t _nk_in     : 16;
  uint64_t _nk_out    : 16;
  uint64_t _nr_in     : 16;
  uint64_t _nr_out    : 16;
  
public:
  NPhases_t() : _nk_in(0), _nk_out(0), _nr_in(0), _nr_out(0) {}

  bool in_phase() const { return _nr_in > _nr_out; }

  void add_nk_in (const size_t nk) {  _nk_in  += nk;  _nr_in++; }
  void add_nk_out(const size_t nk) {  _nk_out += nk;  _nr_out++; }

  size_t nk_in (const bool in_phase = true) const { return in_phase ? _nk_in : _nk_out; };
  size_t nk_out(const bool in_phase = true) const { return in_phase ? _nk_out : _nk_in; };
  size_t nr_in (const bool in_phase = true) const { return in_phase ? _nr_in : _nr_out; };
  size_t nr_out(const bool in_phase = true) const { return in_phase ? _nr_out : _nr_in; };

  bool good() const 
  {
    return true;
    const unsigned nr_max = (_nr_out > _nr_in) ? _nr_out : _nr_in;
    const unsigned nr_min = (_nr_out < _nr_in) ? _nr_out : _nr_in;
    return (nr_max > 2 && nr_max > 10 * nr_min);
  }
  
  
  float score(const bool in_phase = true) const 
  {
    return (in_phase) ?
      float(_nr_out) / float(_nr_in  * _nr_in  + 1) :
      float(_nr_in)  / float(_nr_out * _nr_out + 1);
    
    /*
      if (nr_out * nr_in != 0) {
      if (nr_out > nr_in) {
      score += float(nr_out) / float(nr_in);
      }
      if (nr_out < nr_in) {
      score -= float(nr_in) / float(nr_out);
      }
      else {
      score += nr_out - nr_in;
      }
    */
  }
  
  float strength() const 
  {
    const unsigned nr_max = (_nr_out > _nr_in) ? _nr_out : _nr_in;
    const unsigned nr_min = (_nr_out < _nr_in) ? _nr_out : _nr_in;
    
    return float(nr_max) / float(nr_min*nr_min + 1);
  }

};










class PolymorphismCorrelations
{
  typedef StdMap<unsigned, NPhases_t>    PolysNPhases_t;
  typedef vec<PolysNPhases_t>         PolyCorr_t;

  PolyCorr_t  _poly_corr;  // _poly_corr[i1_poly][i2_poly].{nk,nr}_{in,out}
  vec<size_t> _poly_conflicts;

  struct NKmersAB_t 
  {
    uint16_t nk_a;
    uint16_t nk_b;
    NKmersAB_t() : nk_a(0), nk_b(0) {}
  };


public:
  PolymorphismCorrelations(const Polymorphisms & polys,
			   const PairsManager & pairs, 
			   const BaseVecVec & bvv) 
    : _poly_corr(polys.size()),
      _poly_conflicts(polys.size(), 0)
  {
    typedef map<unsigned, NKmersAB_t>  PolysNKmersAB_t;

    const unsigned K = polys.K();
    const size_t n_pairs = pairs.nPairs();

    size_t n_pairs_conflict = 0;
    size_t n_pairs_corr = 0;
    

    for (size_t i_pair = 0; i_pair < n_pairs; i_pair++) {
      vec<size_t> i_bv(2);
      i_bv[0] = pairs.ID1(i_pair);
      i_bv[1] = pairs.ID2(i_pair);

      PolysNKmersAB_t nk_polys_ab[3];
      
      // ---- parse both reads in pair to count poly kmers hits for each poly version (a or b)

      for (unsigned i = 0; i < 2; i++) {  
        
        // the kmers in the read
	SubKmers<BaseVec, Kmer_t> kmer_bv(K, bvv[i_bv[i]]);

	const size_t nk = kmer_bv.n_kmers();
	while (kmer_bv.not_done()) {

          // check if read kmer exists in the 'a' or 'b' poly databases 
	  const KmerIPoly<Kmer_t> & poly_a_rec = polys.kmer_poly_a(kmer_bv.canonical());
	  const KmerIPoly<Kmer_t> & poly_b_rec = polys.kmer_poly_b(kmer_bv.canonical());

	  if (poly_a_rec.is_valid_kmer() && poly_b_rec.is_valid_kmer()) {
            // something went wrong in the poly computation. 
            // poly_a and poly_b should be mutually exclusive.
	    cout << i_pair << " read pair kmer present in both poly versions..." << endl;
	    exit(1);
	  }
	  else if (poly_a_rec.is_valid_kmer()) {
            const size_t i_poly = poly_a_rec.i_poly();
            nk_polys_ab[i][i_poly].nk_a++;
            nk_polys_ab[2][i_poly].nk_a++;
          }
          else if (poly_b_rec.is_valid_kmer()) {
            const size_t i_poly = poly_b_rec.i_poly();
            nk_polys_ab[i][i_poly].nk_b++;
            nk_polys_ab[2][i_poly].nk_b++;
	  }

          kmer_bv.next();
	}

        for (PolysNKmersAB_t::const_iterator it = nk_polys_ab[i].begin(); 
             it != nk_polys_ab[i].end(); it++) {
          const NKmersAB_t & n_ab = it->second;
          if (n_ab.nk_a > 0 && n_ab.nk_b > 0) {
            const size_t i_poly = it->first;
            cout << "conflicting correlation:" 
                 << " i_pair= " << i_pair
                 << "," << i
                 << " i_poly= " << i_poly 
                 << " nk_{a,b}= " << n_ab.nk_a << "  " << n_ab.nk_b 
                 << endl;
            polys.print(i_poly);
          }
        }


        
      }
      
      
      // ---- make sure that there aren't conflicting kmer hits from the read pair 
      
      bool pair_conflict = false;
      for (PolysNKmersAB_t::const_iterator it = nk_polys_ab[2].begin(); 
           it != nk_polys_ab[2].end(); it++) {
	const NKmersAB_t & n_ab = it->second;
	if (n_ab.nk_a > 0 && n_ab.nk_b > 0) {
	  const size_t i_poly = it->first;
          _poly_conflicts[i_poly]++;
          /*
	  cout << "conflicting correlation:" 
	       << " i_pair= " << i_pair 
	       << " i_poly= " << i_poly 
	       << " nk_{a,b}= " << n_ab.nk_a << "  " << n_ab.nk_b 
	       << endl;
	  polys.print(i_poly);
          */
	  pair_conflict = true;
	}
      }


      // ---- for each pair of different polys, add strength of phase correlation

      if (pair_conflict) {
        n_pairs_conflict++;
      }
      else {

        bool pair_correlates = false;
	for (PolysNKmersAB_t::const_iterator it1 = nk_polys_ab[2].begin(); 
             it1 != nk_polys_ab[2].end(); /* -- empty on purpose -- */ ) {

	  const size_t &   i_poly1 = it1->first;
	  const NKmersAB_t & n_ab1 = it1->second;
	  
	  it1++; // the starting point for it2
	  for (PolysNKmersAB_t::const_iterator it2 = it1; it2 != nk_polys_ab[2].end(); it2++) {

            pair_correlates = true;
	    const size_t &   i_poly2 = it2->first;
	    const NKmersAB_t & n_ab2 = it2->second;
            
            // for each poly the nk_a and nk_b should be mutually exclusive
            // i.e., no conflicts at this point

	    const bool in_phase  = ((n_ab1.nk_a && n_ab2.nk_a) || (n_ab1.nk_b && n_ab2.nk_b));
	    const bool out_phase = ((n_ab1.nk_a && n_ab2.nk_b) || (n_ab1.nk_b && n_ab2.nk_a));
	    const size_t nk_hits = (n_ab1.nk_a + n_ab1.nk_b) + (n_ab2.nk_a + n_ab2.nk_b);
	  
	    if (in_phase && !out_phase) {
	      _poly_corr[i_poly1][i_poly2].add_nk_in(nk_hits);
	      _poly_corr[i_poly2][i_poly1].add_nk_in(nk_hits);
	    } 
	    if (out_phase && !in_phase) {
	      _poly_corr[i_poly1][i_poly2].add_nk_out(nk_hits);
	      _poly_corr[i_poly2][i_poly1].add_nk_out(nk_hits);
	    }
	  }
	}

        if (pair_correlates) n_pairs_corr++;
      }


    } // for (size_t i_pair = 0; i_pair < n_pairs; i_pair++)

    cout << Tag() << setw(12) << n_pairs << " read pairs total." << endl;
    cout << Tag() << setw(12) << n_pairs_corr << " read pairs with good correlation info." << endl;
    cout << Tag() << setw(12) << n_pairs_conflict << " read pairs with bad correlation info." << endl;
    cout << Tag() << setw(12) << (n_pairs - n_pairs_conflict - n_pairs_corr)
         << " read pairs with no correlation info." << endl;
  } 


  size_t size() const { return _poly_corr.size(); }

  size_t n_read_correlations(const size_t i_poly) const 
  { 
    size_t nr_corr = 0;
    
    const PolysNPhases_t & phs = _poly_corr[i_poly];
    for (PolysNPhases_t::const_iterator it2 = phs.begin(); it2 != phs.end(); it2++) {
      const NPhases_t & nph = it2->second;
      nr_corr += nph.nr_in() + nph.nr_out();
    }
    return nr_corr;
  }







  void optimize(vec<bool> * a_is_strong_p) const 
  {
    const size_t n_poly = _poly_corr.size();
    
    // ---- build weighted edges 

    vec<EdgeWeight<float> > edges;
      
    for (size_t i_poly1 = 0; i_poly1 < n_poly; i_poly1++) {
      const PolysNPhases_t & phs = _poly_corr[i_poly1];
      for (PolysNPhases_t::const_iterator it2 = phs.begin(); it2 != phs.end(); it2++) {
        const size_t i_poly2 = it2->first;
        const NPhases_t & nph = it2->second;
        if (i_poly2 > i_poly1 && nph.good()) {
          edges.push_back(EdgeWeight<float>(i_poly1, i_poly2, 1.0 / nph.strength()));
        }
      }
    }

        
    // ---- find spanning tree/forest of polymorphisms
    
    UnionFind poly_tree(n_poly); // ---- initial forest; each node is a tree

    vec<size_t> i_edges;         // ---- indexes of selected edges

    minimum_spanning_tree_kruskal(edges, & poly_tree, & i_edges);

    poly_tree.report_print();

    // ---- label strong edges as strong

    const size_t n_edges_strong = i_edges.size();
    set<pair<unsigned, unsigned> > edge_is_strong;

    for (size_t ii = 0; ii < n_edges_strong; ii++) {
      const EdgeWeight<float> & edge = edges[i_edges[ii]];
      const size_t i_poly1 = edge.iv0;
      const size_t i_poly2 = edge.iv1;
      edge_is_strong.insert(pair<unsigned, unsigned>(i_poly1, i_poly2));
      edge_is_strong.insert(pair<unsigned, unsigned>(i_poly2, i_poly1));
    }

    // ---- propagate phasing info, with a depth-first traversal

    vec<bool> poly_visited(n_poly, false);
    
    for (size_t i_poly1 = 0; i_poly1 < n_poly; i_poly1++) {
      if (!poly_visited[i_poly1]) {
        poly_visited[i_poly1] = true;  // locks phase as 'a is strong', the default

        vec<size_t> i_poly_stack(1, i_poly1);
        
        while (i_poly_stack.size() > 0) {
          const size_t i_poly_cur = i_poly_stack.back();
          
          bool go_deeper = false;
          const PolysNPhases_t & phs = _poly_corr[i_poly_cur];
          for (PolysNPhases_t::const_iterator it2 = phs.begin(); 
               it2 != phs.end() && !go_deeper; it2++) {

            const size_t  i_poly2 = it2->first; 
            const NPhases_t & nph = it2->second;

            // consider only strong links
            // skip already visited
            if (edge_is_strong.count(pair<unsigned, unsigned>(i_poly_cur, i_poly2)) && 
                !poly_visited[i_poly2]) {  

              (*a_is_strong_p)[i_poly2] = 
                (nph.in_phase() == (*a_is_strong_p)[i_poly_cur]); // lock phase

              i_poly_stack.push_back(i_poly2); // put in the stack
              poly_visited[i_poly2] = true;
              go_deeper = true;
            }            
          }
          if (!go_deeper)
            i_poly_stack.pop_back();
        }

      }
    }
    



  }






  float score(const vec<bool> & a_is_strong) const
  {
    float score_tot = 0;
    
    const size_t n_poly = _poly_corr.size();
    for (size_t i_poly1 = 0; i_poly1 < n_poly; i_poly1++) {

      float score = 0;
      
      const PolysNPhases_t & phs = _poly_corr[i_poly1];
      
      for (PolysNPhases_t::const_iterator it2 = phs.begin(); it2 != phs.end(); it2++) {

	const size_t i_poly2 = it2->first;
	const NPhases_t & nph = it2->second;
        const bool in_phase = (a_is_strong[i_poly1] == a_is_strong[i_poly2]);
        score += nph.score(in_phase);
      }

      score_tot += score;
    }
    return score_tot;
  }







  void print(const Polymorphisms & polys, const vec<bool> & a_is_strong) const
  {
    float score_tot = 0;
    const size_t n_poly = _poly_corr.size();
    for (size_t i_poly1 = 0; i_poly1 < n_poly; i_poly1++) {
      cout << "i_poly1= " << setw(6) << i_poly1 
           << " nb_ab= " << setw(6) << polys.base_vec_a(i_poly1).size() 
           << " "        << setw(6) << polys.base_vec_b(i_poly1).size() 
           << fixed << setprecision(1)
           << " kf_ab= " << setw(6) << polys.kmer_freq_a(i_poly1)
           << " "        << setw(6) << polys.kmer_freq_b(i_poly1)
           << " n_conflicts= " << _poly_conflicts[i_poly1] 
           << endl;

      float score = 0;
      const PolysNPhases_t & phs = _poly_corr[i_poly1];
      
      for (PolysNPhases_t::const_iterator it2 = phs.begin(); it2 != phs.end(); it2++) {

	const size_t i_poly2 = it2->first;
	const NPhases_t & nph = it2->second;

        const bool in_phase = (a_is_strong[i_poly1] == a_is_strong[i_poly2]);

        score += nph.score(in_phase);
	cout << "i_poly2= " << setw(6) << i_poly2
             << " nb_ab= " << setw(6) << polys.base_vec_a(i_poly2).size() 
             << " "        << setw(6) << polys.base_vec_b(i_poly2).size() 
             << fixed << setprecision(1)
             << " kf_ab= " << setw(6) << polys.kmer_freq_a(i_poly2)
             << " "        << setw(6) << polys.kmer_freq_b(i_poly2)
	     << " nk_in_out= " << setw(6) << nph.nk_in(in_phase)
	     << " "            << setw(6) << nph.nk_out(in_phase)
	     << " nr_in_out= " << setw(6) << nph.nr_in(in_phase)
	     << " "            << setw(6) << nph.nr_out(in_phase)
             << " n_conflicts= " << _poly_conflicts[i_poly2] 
	     << endl;
      }
      cout << "score= " << score << endl;
      score_tot += score;
    }
    cout << endl << "score_tot= " << score_tot << endl;
    
  }
  

};







class PolymorphismPhases
{
  vec<bool> _phase_a;
  vec<bool> _is_free;
  
public:
  PolymorphismPhases(const Polymorphisms & polys,
                     const PolymorphismCorrelations & poly_corrs) :
    _phase_a(polys.size()),
    _is_free(polys.size(), true)
  {
    const size_t n_polys = polys.size();

    // ---- sort polys by number and quality of correlation links

    vec< pair<unsigned, unsigned> > score_i_poly;
    for (size_t i_poly = 0; i_poly != n_polys; i_poly++) {
      score_i_poly.push_back(make_pair(poly_corrs.n_read_correlations(i_poly), 
                                       i_poly));
    }
    sort(score_i_poly.rbegin(), score_i_poly.rend());

    
    // ---- starting from strongest poly, apply phases and recursively propagate

    for (size_t j_poly = 0; j_poly != n_polys; j_poly++) {
      const size_t i_poly = score_i_poly[j_poly].second;

      if (_is_free[i_poly]) {
	_is_free[i_poly] = false;
	_phase_a[i_poly] = true;     // arbitrary

	vec<size_t> i_poly_stack;
      }
      else {




      }

    }
  }
};














int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_UnsignedInt_Doc(K, "Kmer size.");
  CommandArgument_String_Doc(HEAD_IN, "Looks for reads in <HEAD_IN>.fastb");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT, 
                                       (HEAD_IN != "" ? HEAD_IN + ".no_poly.k" + ToString(K) :
                                        "<HEAD_IN>.no_poly.k<K>"), 
                                       "Output to <HEAD_OUT>.fastb");
  CommandArgument_UnsignedInt_OrDefault_Doc(VERBOSITY, 1, "Verbosity level.");
  CommandArgument_Bool_OrDefault_Doc(REUSE_MAP, False, "Whether to reuse a previously saved map.");



  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;

  NUM_THREADS = configNumThreads(NUM_THREADS);
  
  if (!(K & 1u)) {
    cout << Tag() << "K = " << K << " is not odd. Changing K to " << (K + 1) << "." << endl;
    K++;
  }
 
  // ---- Loading bases and pairing information
  
  cout << Tag() << "Loading read bases." << endl;
  BaseVecVec bvv_frags(HEAD_IN + ".fastb");

  cout << Tag() << "Loading pairing info." << endl;
  PairsManager pairs_frags(HEAD_IN + ".pairs");



  // ---- Build hash map of all read kmers and their affixes

  const String HEAD_IN_K = HEAD_IN + ".k" + ToString(K);
  

  Polymorphisms polys(K);
  KmerSpectrum kspec(K);
  if (REUSE_MAP) {
    polys.read(HEAD_IN_K);
  }
  else {
    polymorphisms_find_parallel(K, bvv_frags, &polys, &kspec, VERBOSITY, NUM_THREADS);
    polys.print_stats();
    
    // ---- Saving map and ambiguities
    {
      cout << Tag() << "Saving ambiguities statistics." << endl;
      polys.write_stats(HEAD_IN_K + ".alt_stats"); 

      cout << Tag() << "Saving polymorphisms." << endl;
      polys.write(HEAD_IN_K);
    }
  }

  

  

  // ---- Correlate polys with paired reads

  cout << Tag() << "Finding polymorphism correlations." << endl;
  const PolymorphismCorrelations poly_corr_frags(polys, pairs_frags, bvv_frags);

  vec<bool> a_is_strong(polys.size(), false);
  //poly_corr_frags.print(polys, a_is_strong);
  

  cout << Tag() << "Optimizing phases." << endl;
  poly_corr_frags.optimize(&a_is_strong);

  poly_corr_frags.print(polys, a_is_strong);

  



  const String head_poly_in  = HEAD_IN + ".poly.k" + ToString(K) + ".in";
  const String head_poly_out = HEAD_IN + ".poly.k" + ToString(K) + ".out";
  
  cout << Tag() << "Saving kept ambiguity branches to '" << head_poly_in << ".fasta'." << endl;
  cout << Tag() << "Saving discarded ambiguity branches to '" << head_poly_out << ".fasta'." << endl;
  polys.write_fastas(head_poly_in, head_poly_out);
  


  /*



  // ---- Fixing reads

  const size_t n_fixed = 
  polymorphisms_remove_parallel(polys, &bases, &quals, VERBOSITY, NUM_THREADS);


  // ---- Saving reads

  bvv_frags.WriteAll(HEAD_OUT + ".fastb");
  */

  cout << Tag() << "Done." << endl;


}
