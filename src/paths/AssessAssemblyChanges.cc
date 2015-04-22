///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AssessAssemblyChanges
// Wrapper module for call to assess_assembly_changes()


#include "MainTools.h"
#include "Basevector.h"
#include "Fastavector.h"
#include "efasta/EfastaTools.h"

#include "paths/AssessAssemblyChangesCore.h"






int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN_HEAD);

  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String_OrDefault(SUBDIR, "test");
  CommandArgument_String_OrDefault(HEAD, "extended40.shaved");
  CommandArgument_String_OrDefault(SCAFFOLDS_IN, "linear_scaffolds0.patched");
  CommandArgument_String_OrDefault(PATCHDIR, "patch");
  CommandArgument_String_OrDefault_Doc(TIGS, "", "if unspecified, process all "
                                       "contigs; otherwise it is one of the following: \n(a) a list of contig "
                                       "ids (in ParseIntSet format) or \n(b) the letter s followed by a list of "
                                       "scaffolds or \n(c) s<scaffold id>.<list of indices of contigs in the "
                                       "scaffold");
  CommandArgument_String_OrDefault_Doc(TARGETS, "", "for assessment, list of "
                                       "reference targets to use, in ParseIntSet format, to speed up alignment");
  CommandArgument_Bool_OrDefault(VERBOSE, True);
  EndCommandArguments;

  // Begin.

  double clock = WallClockTime();

  // Define directories, etc.

  String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  cout << Date() << ": " << run_dir << endl;
  String pdir = run_dir + "/" + PATCHDIR;
  String ch_head = pdir + "/PostPatcher." + SCAFFOLDS_IN + ".";

  // Note some files that are needed.

  String TIGS_file = ch_head + "TIGS";
  String tigsa_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";

  String efasta_file = IN_HEAD + ".efasta";
  String fasta_file = IN_HEAD + ".fasta";

  String scaffolds_tigs_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.vecfasta";

  // Get total number of contigs.

  const size_t n_tigs = MastervecFileObjectCount(scaffolds_tigs_file);

  
  // Parse arguments.

  vec<int> tigs_in;

  {
    if (TIGS == "") {
      for (size_t j = 0; j < n_tigs; j++)
        tigs_in.push_back(j);
    }
    else if (TIGS.Contains("s", 0)) {
      TIGS = TIGS.After("s");
      vec<superb> scaffolds;
      ReadSuperbs(sub_dir + "/" + SCAFFOLDS_IN + ".superb", scaffolds);
      if (TIGS.Contains(".")) {
        int scaffold = TIGS.Before(".").Int();
        ForceAssertLt(scaffold, scaffolds.isize());
        vec<int> spos;
        ParseIntSet(TIGS.After("."), spos);
        for (int j = 0; j < spos.isize(); j++)
          tigs_in.push_back(scaffolds[scaffold].Tig(spos[j]));
      }
      else {
        vec<int> s;
        ParseIntSet(TIGS, s);
        for (int i = 0; i < s.isize(); i++) {
          int scaffold = s[i];
          ForceAssertLt(scaffold, scaffolds.isize());
          for (int j = 0; j < scaffolds[scaffold].Ntigs(); j++)
            tigs_in.push_back(scaffolds[scaffold].Tig(j));
        }
      }
    }
    else 
      ParseIntSet(TIGS, tigs_in);
  }

  // const reference garantees no change
  const vec<int>& tigs = tigs_in;


  assess_assembly_changes(TARGETS, scaffolds_tigs_file, efasta_file, fasta_file, data_dir, 
                          tigs);

  cout << "time used = " << TimeSince(clock) << endl;

}
