///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "efasta/EfastaTools.h"
#include "math/Functions.h"
#include "reporting/PerfStat.h"


int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( HEAD );
  EndCommandArguments;

  String assembly_file = HEAD + ".assembly.efasta";
  String contigs_file = HEAD + ".contigs.efasta";

  cout << Date() << " Loading " << assembly_file << endl;
  VecEFasta assembly;
  LoadEfastaIntoStrings(assembly_file, assembly);
  u_int n_scaffolds = assembly.size();

  cout << Date() << " Loading " << contigs_file << endl;
  VecEFasta contigs;
  LoadEfastaIntoStrings(contigs_file, contigs);
  u_int n_contigs = contigs.size();

  cout << Date() << " Computing assembly stats:" << endl;
  size_t gapped_length = 0;
  vec<size_t> ungapped_scaffold_lengths(n_scaffolds);
  vec<size_t> gapped_scaffold_lengths(n_scaffolds);
  for (size_t i = 0; i < assembly.size(); ++i) {
    size_t ulen = assembly[i].Length1(False);
    size_t glen = assembly[i].Length1(True);
    ungapped_scaffold_lengths[i] = ulen;
    gapped_scaffold_lengths[i] = glen;
    gapped_length += glen;
  }
  cout << "Length including gaps: " << gapped_length << endl;
  PerfStat::log( ) << PerfStat( "gapped_length", "assembly length including gaps", gapped_length );

  size_t ungapped_length = 0;
  vec<size_t> contig_lengths(n_contigs);
  for (size_t i = 0; i < contigs.size(); ++i) {
    size_t len = contigs[i].Length1();
    contig_lengths[i] = len;
    ungapped_length += len;
  }

  ReverseSort(contig_lengths);
  ReverseSort(ungapped_scaffold_lengths);
  ReverseSort(gapped_scaffold_lengths);

  vec<size_t> contig_nstats = NStatistics(contig_lengths);
  vec<size_t> ungapped_nstats = NStatistics(ungapped_scaffold_lengths);
  vec<size_t> gapped_nstats = NStatistics(gapped_scaffold_lengths);

  cout << "Length excluding gaps: " << ungapped_length << endl;
  PerfStat::log( ) << PerfStat( "ungapped_length", "assembly length excluding gaps", ungapped_length );

  cout << "Number of scaffolds: " << n_scaffolds << endl;
  PerfStat::log( ) << PerfStat( "n_scaffolds", "number of scaffolds in assembly", n_scaffolds );

  cout << "Scaffold N50 including gaps: " << gapped_nstats[4] << endl;
  PerfStat::log( ) << PerfStat( "scaffold_n50_gapped", "scaffold N50 including gaps", gapped_nstats[4]);
  cout << "Scaffold N90 including gaps: " << gapped_nstats[8] << endl;
  PerfStat::log( ) << PerfStat( "scaffold_n90_gapped", "scaffold N90 including gaps", gapped_nstats[8]);

  cout << "Scaffold N50 excluding gaps: " << ungapped_nstats[4] << endl;
  PerfStat::log( ) << PerfStat( "scaffold_n50", "scaffold N50 excluding gaps", ungapped_nstats[4]);
  cout << "Scaffold N90 excluding gaps: " << ungapped_nstats[8] << endl;
  PerfStat::log( ) << PerfStat( "scaffold_n90", "scaffold N90 excluding gaps", ungapped_nstats[8]);

  cout << "Number of contigs: " << n_contigs << endl;
  PerfStat::log( ) << PerfStat( "n_contigs", "number of contigs in assembly", n_contigs);

  cout << "Contig N50: " << contig_nstats[4] << endl;
  PerfStat::log( ) << PerfStat( "contig_n50", "contig N50", contig_nstats[4]);
  cout << "Contig N90: " << contig_nstats[8] << endl;
  PerfStat::log( ) << PerfStat( "contig_n90", "contig N90", contig_nstats[8]);

  cout << "Longest contig: " << contig_lengths[0] << endl;
  PerfStat::log( ) << PerfStat( "longest_contig", "longest contig in assembly", contig_lengths[0]);

  cout << "Longest ungapped scaffold: " << ungapped_scaffold_lengths[0] << endl;
  PerfStat::log( ) << PerfStat( "longest_scaffold", "longest scaffold excluding gaps", ungapped_scaffold_lengths[0]);

  int snp_count = 0,indel_count = 0, ambiguous_bases = 0;
  for (size_t i = 0; i < contigs.size(); ++i) {
    int snps = 0, indels = 0;
    int amb =  contigs[i].AmbCount(snps, indels);
    ambiguous_bases += amb;
    snp_count += snps;
    indel_count += indels;
  }
  if (ambiguous_bases) {
    double rate = (double) ambiguous_bases / (double) ungapped_length;
    int rate1 = round(1.0/rate);
    cout << "Ambiguous base rate: " << setprecision(1) << scientific << rate
	 << " (1/" << rate1 << ")" << endl;
    PerfStat::log( ) << PerfStat( "ambig_base_rate", "fraction of ambiguous bases", rate);
  }
  if (snp_count + indel_count) {
    double rate = (double) (snp_count + indel_count) / (double) ungapped_length;
    int rate1 = round(1/rate);
    cout << "Ambiguity event rate: " << setprecision(1) << scientific << rate
	 << " (1/" << rate1 << ")" << endl;
    PerfStat::log( ) << PerfStat( "ambig_event_rate", "rate of ambiguity events", rate);
  }
  if (snp_count) {
    double rate = (double) snp_count / (double) ungapped_length;
    int rate1 = round(1/rate);
    cout << "SNP event rate: " << setprecision(1) << scientific << rate
	 << " (1/" << rate1 << ")" << endl;
    PerfStat::log( ) << PerfStat( "ambig_snp_rate", "rate of snp-like ambiguity", rate);
  }
  if (indel_count) {
    double rate = (double) indel_count / (double) ungapped_length;
    int rate1 = round(1.0/rate);
    cout << "Indel event rate: " << setprecision(1) << scientific << rate
	 << " (1/" << rate1 << ")" << endl;
    PerfStat::log( ) << PerfStat( "ambig_indel_rate", "rate of indel-like ambiguity", rate);
  }
}


