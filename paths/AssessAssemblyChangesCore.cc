///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AssessAssemblyChangesCore.
// Experimental.

// MakeDepend: dependency QueryLookupTable

#include "Basevector.h"
#include "Fastavector.h"
#include "efasta/EfastaTools.h"
#include "FastIfstream.h"
#include "lookup/LookAlign.h"

void assess_assembly_changes(const String & TARGETS, 
                             const String & scaffolds_tigs_file,
                             const String & efasta_file,
                             const String & fasta_file,
                             const String & data_dir,
                             const vec<int> & tigs)
{
  const size_t n_used_tigs = tigs.size();
  

  // ---- load modified contigs

  cout << Date() << ": loading modified contigs." << endl;
  vec<FastaVec> combo;
  vec<int> exp_start(n_used_tigs + 1);
  exp_start[0] = 0;

  {
    // ---- load the efasta for the changed contigs
    
    VecEFasta econtigs_new;
    LoadEfastaIntoStrings(efasta_file, econtigs_new);
    
    for (size_t it = 0; it < n_used_tigs; it++) {
      int tig = tigs[it];
      FastaVecVec contigs;
      contigs.ReadOne(scaffolds_tigs_file, tig);
      combo.push_back(contigs[0]);
    }

    for (size_t it = 0; it < n_used_tigs; it++) {
      // Expand the contig.  Note that this is obviously unsafe, as there could
      // be a very large number of elements in the expansion.
      vec<FastaVec> contig_expanded;
      econtigs_new[it].ExpandTo(contig_expanded);
      combo.append(contig_expanded);
      exp_start[it + 1] = exp_start[it] + contig_expanded.size();
    }
  }


  // ---- call QueryLookupTable for alignments

  cout << Date() << ": aligning to reference" << endl;
  const String TARGETS_arg = (TARGETS == "" ? ""
                              : " TARGETS_TO_PROCESS=\"" + TARGETS + "\"");
  fast_pipe_ifstream in("QueryLookupTable K=12 MM=12 MC=0.15 "
                        "SEQS=" + fasta_file + " L=" + data_dir + "/genome.lookup "
                        "PARSEABLE=True" + TARGETS_arg);
  vec<look_align> aligns;
  while(1) {
    String line;
    getline(in, line);
    if (in.fail()) break;
    if (!line.Contains("QUERY", 0)) continue;
    look_align la;
    la.ReadParseable(line);
    aligns.push_back(la);
  }




  cout << Date() << ": separate new and old alignments" << endl;
  vec< vec<look_align> > old_aligns(n_used_tigs);
  vec< vec<look_align> > new_aligns(n_used_tigs);
  for (size_t i = 0; i < aligns.size(); i++) {
    int id = aligns[i].query_id;
    if (id < (signed)n_used_tigs) {
      old_aligns[id].push_back(aligns[i]);
    }
    else {
      for (size_t it = 0; it < n_used_tigs; it++) {
        if (id - (signed)n_used_tigs < exp_start[it+1]) {
          new_aligns[it].push_back(aligns[i]);
          break;
        }
      }
    }
  }


  cout << Date() << ": pick best" << endl;
  vec< vec<look_align> > old_aligns_best(n_used_tigs);
  vec< vec<look_align> > new_aligns_best(n_used_tigs);
  for (size_t it = 0; it < n_used_tigs; it++) {
    if (old_aligns[it].nonempty()) 
      old_aligns_best[it] = old_aligns[it];
    look_align best;
    int min_new_errors = 1000000000;
    for (size_t ia = 0; ia < new_aligns[it].size(); ia++) {
      if (new_aligns[it][ia].Errors() < min_new_errors) {
        best = new_aligns[it][ia];
        min_new_errors = new_aligns[it][ia].Errors();
      }
    }
    if (min_new_errors < 1000000000) 
      new_aligns_best[it].push_back(best);
  }
  vec<int> tids;
  for (size_t it = 0; it < n_used_tigs; it++) {
    for (size_t ia = 0; ia < old_aligns_best[it].size(); ia++)
      tids.push_back(old_aligns_best[it][ia].target_id);
    for (size_t ia = 0; ia < new_aligns_best[it].size(); ia++)
      tids.push_back(new_aligns_best[it][ia].target_id);
  }

  UniqueSort(tids);

  cout << Date() << ": get " << data_dir << "/genome.fastb" << endl;
  BaseVecVec targets;
  targets.Read(data_dir + "/genome.fastb", tids);
  cout << "\n=============================================================="
       << "======================\n";
  cout << "\nALIGNMENTS OF OLD AND NEW CONTIGS\n" << endl;
  int old_errs = 0;
  int new_errs = 0;
  int old_fails = 0;
  int new_fails = 0;
  for (size_t it = 0; it < n_used_tigs; it++) {
    int tig = tigs[it];
    cout << "==========================================================="
         << "=========================\n" << "\nCONTIG " << tig << "\n" 
         << endl;
    Bool same = False;
    if (old_aligns_best[it].empty() && new_aligns_best[it].empty())
      same = True;
    if (old_aligns_best[it].nonempty() && new_aligns_best[it].nonempty()) {
      const look_align& la1 = old_aligns_best[it][0];
      const look_align& la2 = new_aligns_best[it][0];
      if (combo[la1.query_id] == combo[la2.query_id]) 
        same = True;
    }
    cout << "OLD" << (same ? " = NEW" : "") << "\n" << endl;
    if (old_aligns_best[it].empty()) {
      cout << "(no alignment found)" << endl;
      if (!same) old_fails++;
    }
    else {
      const look_align& la = old_aligns_best[it][0];
      la.PrintReadableBrief(cout);
      la.PrintVisual(cout, combo[la.query_id],
                     targets[Position(tids, la.target_id)]);
      if (!same) old_errs += la.Errors();
    }
    if (!same) {
      cout << "\nNEW\n" << endl;
      if (new_aligns_best[it].empty()) {
        cout << "(no alignment found)" << endl;
        new_fails++;
      }
      else { 
        const look_align& la = new_aligns_best[it][0];
        la.PrintReadableBrief(cout);
        la.PrintVisual(cout, combo[la.query_id],
                       targets[Position(tids, la.target_id)]);
        new_errs += la.Errors();
      }
    }
  }
  cout << "\nimprovement = " << (old_errs > new_errs ? "+" : "");
  cout << old_errs - new_errs << ", ";
  PRINT2(old_fails, new_fails);
}

