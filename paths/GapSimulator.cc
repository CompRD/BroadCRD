///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// author: Filipe Ribeiro 04/2011
//
//

#include <stdlib.h> // rand(), RAND_MAX

#include "MainTools.h"
#include "Vec.h" 
#include "feudal/BinaryStream.h"

#include "math/IntDistribution.h"
#include "paths/OffsetDistribution.h"

static inline 
String Tag(String S = "GS") { return Date() + " (" + S + "): "; } 


vec<IntDistribution> distributions_read(const String & HEAD) 
{
  const String fn = HEAD + ".distribs";
  
  cout << Tag() << "Reading binary distributions from '" << fn << "'." << endl;
  BinaryReader reader(fn.c_str());
  size_t nl;
  reader.read(&nl);
  vec<IntDistribution> dists(nl);
  reader.readItr(dists.begin(), dists.end());
  return dists;
}






int rand_max(const int max) { return rand() % max; }
double rand_01() { return double(rand()) / double(RAND_MAX); }



vec<GapBridge> bridges_set(const int S1, 
                           const int S2, 
                           const int G, 
                           const int COV, 
                           const vec<IntDistribution> & dists)
{

   vec<GapBridge> bridges;
   size_t i_dist = 0;
   int len = 101;
   int fw0 = 500;
   int fw1 = fw0 + len;
   int bw0 = 500;
   int bw1 = bw0 - len;


   bridges.push_back(GapBridge(i_dist, S1, S2,
                               ContigReadIndex::contig2(S2, fw0, fw1),
                               ContigReadIndex::contig1(S1, bw0, bw1)));


   return bridges;

   bridges.push_back(GapBridge(i_dist, S1, S2,
                               ContigReadIndex::contig1(S1, fw0, fw1),
                               ContigReadIndex::contig2(S2, bw0, bw1)));

   bridges.push_back(GapBridge(i_dist, S1, S2,
                               ContigReadIndex::contig1(S1, fw0, fw1),
                               ContigReadIndex::not_aligned(len)));


   bridges.push_back(GapBridge(i_dist, S1, S2,
                               ContigReadIndex::contig2(S2, fw0, fw1),
                               ContigReadIndex::not_aligned(len)));
   return bridges;
}



void bridges_simulate(const int S1, 
                      const int S2, 
                      const int G, 
                      const int COV, 
                      const vec<IntDistribution> & dists,
                      vec<GapBridge> * p_bridges)
{
  const int i0_1 = 0;
  const int i1_1 = S1;
  const int i0_2 = i1_1 + G;
  const int i1_2 = i0_2 + S2;

  const int i0 = (i0_1 < i0_2) ? i0_1 : i0_2;
  const int i1 = (i1_1 > i1_2) ? i1_1 : i1_2;

  const int S_tot = i1 - i0;
  const int len = 101;

  const int sz_min = dists[0].x_min();
  const int sz_max = dists[0].x_max();
  
  const int fw0_min = (sz_max > 0) ? i0 - sz_max : i0;
  const int fw0_max = (sz_min < 0) ? i1 - sz_min : i1;
  const int fw0_range = fw0_max - fw0_min;

  const size_t n = COV * fw0_range / len;
  cout << Tag() << "Simulating " << n << " bridges." << endl;
  cout << Tag() << "S_tot = " << S_tot << endl;
  cout << Tag() << "fw0_min = " << fw0_min << endl;
  cout << Tag() << "fw0_max = " << fw0_max << endl;
  cout << Tag() << "fw0_range = " << fw0_range << endl;

  size_t n10 = 0;
  size_t n12 = 0;
  size_t n20 = 0;
  size_t n21 = 0;

  for (size_t i = 0; i != n; i++) {
    const size_t i_dist = 0;
    const IntDistribution & dist = dists[i_dist];
    
    
    // uniformly distributed fw read in (fw0_min, fw0_max)
    const int fw0_tot = fw0_min + rand_max(fw0_range);
    const int fw1_tot = fw0_tot + len - 1;

    // sample invariant size from the empirical distribution
    const int inv_sz = dist.x_prob_le(rand_01()); 

    const int bw0_tot = fw0_tot + inv_sz - 1;
    const int bw1_tot = bw0_tot - len + 1;

    if (fw1_tot < i1_1 && fw0_tot >= i0_1) { // fw read in contig 1
      const int fw0 = fw0_tot;
      const int fw1 = fw1_tot;
      if (bw0_tot < i1_2 && bw1_tot >= i0_2) { // bw read in contig 2
        const int bw0 = bw0_tot - i0_2;
        const int bw1 = bw1_tot - i0_2;
        
        p_bridges->push_back(GapBridge(i_dist, S1, S2,
                                       ContigReadIndex::contig1(S1, fw0, fw1),
                                       ContigReadIndex::contig2(S2, bw0, bw1)));
        n12++;
      }
      else {
        p_bridges->push_back(GapBridge(i_dist, S1, S2,
                                       ContigReadIndex::contig1(S1, fw0, fw1),
                                       ContigReadIndex::not_aligned(len)));
        n10++;
      }
    } 
    else if (fw1_tot < i1_2 && fw0_tot >= i0_2) { // fw read in contig 2
      const int fw0 = fw0_tot - i0_2;
      const int fw1 = fw1_tot - i0_2;

      if (bw0_tot < i1_1 && bw1_tot >= i0_1) { // bw read in contig 1
        const int bw0 = bw0_tot;
        const int bw1 = bw1_tot;
        
        p_bridges->push_back(GapBridge(i_dist, S1, S2,
                                       ContigReadIndex::contig2(S2, fw0, fw1),
                                       ContigReadIndex::contig1(S1, bw0, bw1)));
        n21++;
      }
      else {
        p_bridges->push_back(GapBridge(i_dist, S1, S2,
                                       ContigReadIndex::contig2(S2, fw0, fw1),
                                       ContigReadIndex::not_aligned(len)));
        n20++;
      }
    }

  }

  cout << Tag() << "n12 [1] = " << n12 << endl;
  cout << Tag() << "n10 [2] = " << n10 << endl;
  cout << Tag() << "n21 [3] = " << n21 << endl;
  cout << Tag() << "n20 [4] = " << n20 << endl;
}









int main(int argc, char *argv[]) 
{
  RunTime();
  BeginCommandArguments;
  CommandArgument_String_Doc        (HEAD,      "Looks for <HEAD>.distribs.");
  CommandArgument_Int_Doc           (S1,        "Size of contig 1.");
  CommandArgument_Int_Doc           (S2,        "Size of contig 1.");
  CommandArgument_Int_Doc           (OFFSET,    "Offset of S2 realtive to S1.");
  CommandArgument_Int_OrDefault_Doc (COV, 30,   "Coverage.");
  CommandArgument_Int_OrDefault_Doc (CYCLES, 1, "Number of cycles.");
  EndCommandArguments;

  const int G = OFFSET - S1;

  const vec<IntDistribution> dists = distributions_read(HEAD);
  
  //vec<IntDistribution> dists(1, IntDistribution::gaussian(2000, 20));
  //dists[0].to_text_file("gaussian");

  vec<double> qs(CYCLES);
  vec<double> qs_short(CYCLES);

  if (1) {
  for (int cycle = 0; cycle != CYCLES; cycle++) {
    cout << "cycle= " << cycle;

    vec<GapBridge> bridges;
    bridges_simulate(S1, S2, G, COV, dists, &bridges);

    IntDistribution dist_offset = 
      offset_distribution_compute(dists, bridges);
    
    const int os_min = dist_offset.x_min();
    const int os_max = dist_offset.x_max();
    
    qs[cycle] = dist_offset.prob_le(OFFSET);
    //qs_short[cycle] = dist_offset_short.prob_le(G);
    
    dist_offset.to_text_file("offset");
    //dist_offset_short.to_text_file("gap_short_test");
    

    if (cycle == 0) {
      {
        //IntDistribution dist_offset[4];
        for (size_t j = 0; j != 4; j++) {
          IntDistribution dist =
            offset_distribution_compute(dists, bridges, NULL, j+1);
          dist.to_text_file("offset_test" + ToString(j+1));
        }
      }

      if (1) {
        {
          IntDistribution dist_bridge1_2 = 
            distribution_bridge_fw_given_offset(dists[0],
                                                S1, S2, 
                                                101, 101,
                                                os_min, os_max, 1);
        
          dist_bridge1_2.to_text_file("bridge1_test");
        }
        {
          IntDistribution dist_bridge1_01 =
            distribution_bridge_fw_given_offset(dists[0],
                                                S1, S2, 
                                                101, 101,
                                                os_min, os_max, 2);
        
          dist_bridge1_01.to_text_file("bridge2_test");
        }
        {
          IntDistribution dist_bridge2_1 =
            distribution_bridge_fw_given_offset(dists[0],
                                                S1, S2, 
                                                101, 101,
                                                os_min, os_max, 3);
        
          dist_bridge2_1.to_text_file("bridge3_test");
        }
        {
          IntDistribution dist_bridge2_02 = 
            distribution_bridge_fw_given_offset(dists[0],
                                                S1, S2, 
                                                101, 101,
                                                os_min, os_max, 4);
          
          dist_bridge2_02.to_text_file("bridge4_test");
        }
      }
    }
    

  }
  cout << endl;
  

  sort(qs.begin(), qs.end());
  {
    ofstream os;
    os.open("qs.txt");
    for (size_t c = 0; c != qs.size(); c++)
      os << qs[c] << endl;
    os.close();
  }
  /*
  sort(qs_short.begin(), qs_short.end());
  {
    ofstream os;
    os.open("qs_short.txt");
    for (size_t c = 0; c != qs_short.size(); c++)
      os << qs_short[c] << endl;
    os.close();
  }
  */

  }
  cout << Tag() << "Done!" << endl;

}
