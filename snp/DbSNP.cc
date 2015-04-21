/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// DbSNP.cc
/// Method implementations. See DbSNP.h
///

#include "snp/DbSNP.h"
#include "FastIfstream.h"
#include "String.h"
#include "TokenizeString.h"

/// Load dbSNP data from a flat file. Currently, the only recognized
/// file format is: space or tab separated lines with the following fields
/// contig pos id strand alleles.
/// contig, pos are integers (position on the genome), SNP id is unused, strand
/// must be "+" or "-" sign, and alleles is a "/" separated (no spaces!) list
/// of alternative alleles, e.g. A/C. NOTE: position is always assumed to be given
/// on the positive strand, the "strand" field affects only how alleles are interpreted
/// (complement is taken for each allele base if strand=="-").
void DbSNP::Load(const String & fname) {
     vec<char> allele_separator;
     allele_separator.push_back('/');

     fast_ifstream snp_in(fname);
     String line;
     unsigned int line_number = 0;
     unsigned int SNPs_read = 0;

     while ( 1 ) {
       line_number++;
       //       if ( line_number == 10001 ) break;
       getline( snp_in, line);
       if ( snp_in.fail() ) break;
       line = line.Trim(" \t");
       vec<String> tokens;
       Tokenize(line,tokens);

       unsigned int contig = tokens[0].Int();
       unsigned int pos = tokens[1].Int() - 1;
       bool rc = false;
       if ( tokens[3] == "-" ) {
	 rc=true;
       } else {
	 if ( tokens[3] != "+" ) {
	   cout << fname << " : line " << line_number << ": unrecognized strand symbol" << endl;
	   exit(1);
	 }
       }

       String allele_set = tokens[4];

       tokens.clear();
       Tokenize(allele_set,allele_separator,tokens);

       if ( tokens.size() > 4 ) {
	 cout << fname << " : line " << line_number << ": " << tokens.size() << " alleles found" << endl;
	 exit(1);
       }
       switch (tokens.size()) {
       case 0: cout << fname << " : line " << line_number << ": no alleles found" << endl; exit(1);
       case 1: cout << fname << " : line " << line_number << ": no alternative alleles found" << endl; exit(1);
       }

       for ( unsigned int i = 0 ; i < tokens.size() ; i++ ) {
	 if ( tokens[i].size() > 1 || !Base::isCanonicalBase(tokens[i][0]) ) {
	   cout << fname << " : line " << line_number << ": unrecognized allele " << tokens[i] << endl;
	   exit(1);
	 }
	 if ( rc ) tokens[i][0] = GetComplementaryBaseChar(tokens[i][0]);
       }

       SNPAlleles SNP;
       if ( ! HasFeature(contig,pos) ) {
	 SNP.SetAlleles(tokens); // initialize with reported alleles
	 SNPs_read++;
       } else {
	 SNP = Get(contig, pos);  // pull out what we got so far
	 for ( unsigned int i = 0 ; i< tokens.size(); i++ ) SNP.AddAllele(tokens[i][0]);
       }
       Add(contig,pos,SNP);
       //       cout << tokens[0]<<"-"<<tokens[1]<<":"<<tokens[3]<<":"<<tokens[4]<<endl;
     }

     cout << "done." << endl;
     cout << "lines read: " << line_number << endl;
     cout << "unique SNP positions read: " << SNPs_read << endl;

}

/// Check SNP collection against the provided reference: verify if the reference base
/// is reported among the alternative alleles. All the SNP positions stored in this map
/// are looked up directly in the <ref>, so the reference must be consistent with
/// genomic positions used in the map. The second argument specifies what course
/// of action should be taken if the reference base is missing (this are bitflags,
/// that can be joined with OR when it makes sense). ABORT: abort immediately when SNP
/// missing the reference allele is found (cannot be OR'd - any other flag overrides this);
/// WARN: issue a warning into stdout for each SNP missing the reference base;
/// UPDATE: add reference base if it is missing (silently; use with WARN if needed);
/// DISCARD: discard all SNPs missing the reference base (silently; use with WARN if needed).
/// UPDATE and DISCARD can not be used together.
///
void DbSNP::CheckAgainstRef(const vecbasevector & ref, int action) {
  if ( ( (action & UPDATE) != 0 ) && ( (action & DISCARD) != 0 ) ) {
    cout << "DbSNP::CheckAgainstRef: Can not use UPDATE and DISCARD flags simultaneously." << endl;
    exit(1);
  }

  unsigned int total_SNPs = 0;
  unsigned int missing_ref = 0;
  for ( unsigned int i = 0 ; i < storage.size() ; i++ ) {
    storage_map::iterator it_end = storage[i].end();
    vec<unsigned int> to_remove;
    const basevector & contig = ref[i];

    for ( storage_map::iterator it = storage[i].begin() ; it != it_end; ++it ) {
      total_SNPs++;
      base_t ref_base = contig[ ( it->first ) ];
      if ( (it->second).Contains( ref_base ) ) continue;
      // ooops, ref base is missing:
      missing_ref++;
      if ( action == ABORT || (action & WARN) != 0 ) {
	cout << "DbSNP: no reference base among the alternative alleles" << endl;
	cout << "       SNP:" << i <<":" << (it->first) << " " << (it->second)
	     << " REF: " << as_char(ref_base) << endl;
      }
      if ( action == ABORT ) exit(1);
      if ( (action & UPDATE) != 0 ) (it->second).AddAllele(ref_base);
      if ( (action & DISCARD) != 0 ) to_remove.push_back(it->first);
    }
    for ( unsigned int i = 0 ; i < to_remove.size() ; i++ ) {
      storage_map::iterator it = storage[i].find(to_remove[i]);
      if ( it == storage[i].end() ) {
	cout << "Fatal internal error: can not locate previously found discordant SNP" << endl;
	exit(1);
      }
      storage[i].erase(it);
    }
  }

  cout << "Reference check:" << endl;
  cout << "  Total SNPs:                  " << setw(8) << total_SNPs << endl;
  cout << "  SNPs missing reference base: " << setw(8) << missing_ref;
  if ( (action & UPDATE)  != 0 ) cout << " (updated)" << endl ;
  if ( (action & DISCARD)  != 0 ) {
    cout << " (discarded)" << endl ;
    cout << "  SNPs remaining:              " << setw(8) << (total_SNPs-missing_ref) << endl;
  }
}
