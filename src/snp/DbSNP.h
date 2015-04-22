/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// DbSNP.  Simple wrapper class that can load and parse dbsnp data from flat
/// table and keep SNP positions alleles. Memory efficient and can be queried.
///

#include "Basevector.h"
#include "feudal/TrackingAllocator.h"
#include <ext/hash_map>

/// A wrapper class that keeps collection of up to 4 bases ("alternative alleles")
/// and thus represents a SNP. This class does not store allele frequency information,
/// but it's instance takes only one byte of memory. Methods allow adding alternative
/// alleles, checking the total number of recorder alternative alleles, and querying
/// if a specific base is kept by a given instance. NOTE: for the sake of efficiency,
/// this class is unsafe in a sense that there is no way to have "empty" SNP with 0
/// alleles. Default constructor creates a "SNPAllele" that already contains one base,
/// but it is not the actual "reference base", but rather it's always an A. Make sure
/// you use setter methods to initialize content correctly! The intended use is to keep
/// the reference base in the list among all the alternative alleles (all the animals are equal),
/// but this contract is never actually enforced. Bases (alleles) are unsorted and stored in the
/// order they were placed into the list.

struct SNPAlleles {

  /// ctor; in order to use memory more efficiently, a concession is made:
  /// default constructor creates not an "empty" SNP record, but one that
  /// already contains one base - A. Make sure that Set... method is used to
  /// corectly initialize alleles in a SNP created with default constructor!!
  SNPAlleles() : data(0) {}

  /// returns base reported as i-th allele at this SNP;
  /// 0-based.
  base_t operator[] (unsigned int i) const {
    // implementation of this class tries to save memory, so all
    // alternative alleles (up to 4, apparently), are stored in one byte
    // of data; the actual number of stored alleles is not stored separately.
    // Instead, it is guaranteed by the contract that if K alleles are reported
    // for the given SNP (K<=4), then bit pairs at offsets 2*[0,...,K-1] store
    // distinct allele bases, and all the unused bit pairs (if any) store the
    // same base as found at bit offset 0. This contract is used below:
    ForceAssertLt(i,4u);
    base_t alt1 = data & 0x3;
    base_t alt2 = data;
    switch(i) {
    case 0: return alt1;
    case 1: alt2 >>= 2 ;
    case 2: alt2 >>= 2 ;
    case 3: alt2 >>= 2 ;
    }
    alt2 &= 0x3;
    if ( alt2 == alt1 ) {
      cout << "ERROR: allele variant " << (i+1) << " does not exist in the SNP" << endl;
      exit(1);
    }
    return alt2;
  }

  /// Returns the total number of alleles at this SNP; if there are no alternative
  /// alleles (this is not actually a SNP), then returns 1, i.e. the reference allele
  /// (or "default" allele, see constructor docs) is counted.
  unsigned int NAlleles() const {
    base_t alt1 = data & 0x3;
    unsigned int n = 1;
    unsigned char data_copy = data;
    // below is nothing by a loop, unrolled:
    data_copy >>= 2;
    if ( ( data_copy & 0x3 ) == alt1 ) return n; // [1]==[0]
    n++; // n=2
    data_copy >>= 2;
    if ( ( data_copy & 0x3 ) == alt1 ) return n; // [2]==[0]
    n++; // n=3
    data_copy >>= 2;
    if ( ( data_copy & 0x3 ) == alt1 ) return n; // [3]==[0]
    return 4;
  }

  /// Returns <true> if alternative allele <base> is already stored for this SNP
  bool Contains(base_t base) const {
    base_t alt1 = data & 0x3;
    if ( base == alt1 ) return true;

    base_t alt2 = ( data >> 2 ) & 0x3;
    if ( alt2 == alt1 ) return false; // reached the end
    if ( alt2 == base ) return true; //  allele found
    alt2 = ( data >> 4 ) & 0x3;
    if ( alt2 == alt1 ) return false; // reached the end
    if ( alt2 == base ) return true; //  allele found
    alt2 = ( data >> 6 ) & 0x3;
    if ( alt2 == base ) return true; //  allele found
    return false;
  }


  /// adds allele <alt> to the list of this SNP's alternative alleles. If
  /// <alt> is already among this SNP's alleles, nothing will be done.
  /// NOTE: this method is unsafe as it can only add alternatives to a
  /// non-empty SNP (i.e. at least one allele should have been set up previously
  /// with Set...). An uninitialized SNP (default conctructor) does look like a
  /// "non-empty" one, so using this method can lead to unexpected results.
  void AddAllele(base_t alt) {
    base_t alt1 = data & 0x3;
    if ( alt1 == alt ) return; // we already have it
    base_t alt2 = (data >> 2) & 0x3;
    // if past end, add <alt>
    if ( alt2 == alt1 ) data |= (( alt & 0x3  ) << 2);
    if ( alt2 == alt ) return;

    alt2 = (data >> 4) & 0x3;
    if ( alt2 == alt1 ) data |= (( alt & 0x3  ) << 4);
    if ( alt2 == alt ) return;

    alt2 = (data >> 6) & 0x3;
    if ( alt2 == alt1 ) data |= (( alt & 0x3  ) << 6);
    // no need for the check below: if all 4 alleles
    // are already recorded, one of them (this one since
    // we got to this point) must be <alt>!
    // if ( alt2 == alt ) return;
  }

  /// sets a single allele for this SNP (so strictly
  /// speaking this is not a SNP then!). Aborts if <alt> is not a valid base.
  /// This method is safe as it overrides previous content, if any.
  void SetAlleles(base_t alt) {
    ForceAssertLt(alt,4);
    data = alt;
    data <<= 2;
    data |= alt;
    data <<= 2;
    data |= alt;
    data <<= 2;
    data |= alt;
  }

  // sets two alleles for this SNP;
  // aborts if any of the alternative alleles is not a valid base or
  // if the two allele bases are the same.
  /// This method is safe as it overrides previous content, if any.
  void SetAlleles(base_t alt1, base_t alt2) {
    ForceAssertLt(alt1,4);
    ForceAssertLt(alt2,4);
    ForceAssertNe(alt1,alt2);
    data = alt1;
    data <<= 2;
    data |= alt1;
    data <<= 2;
    data |= alt2;
    data <<= 2;
    data |= alt1;
  }

  /// sets three alternative alleles for this SNP;
  /// aborts if any of the alternative alleles is not a valid base or
  /// if any pair of allele bases are the same.
  /// This method is safe as it overrides previous content, if any.
  void SetAlleles(base_t alt1, base_t alt2, base_t alt3) {
    ForceAssertLt(alt1,4);
    ForceAssertLt(alt2,4);
    ForceAssertLt(alt3,4);
    ForceAssertNe(alt2, alt1);
    ForceAssertNe(alt3, alt1);
    ForceAssertNe(alt2, alt3);
    data = alt1;
    data <<= 2;
    data |= alt3;
    data <<= 2;
    data |= alt2;
    data <<= 2;
    data |= alt1;
  }


  /// sets four alternative alleles for this SNP;
  /// aborts if any of the alternative alleles is not a valid base or
  /// if any pair of passed allele bases are the same.
  /// This method is safe as it overrides previous content, if any.
  void SetAlleles(base_t alt1, base_t alt2, base_t alt3, base_t alt4) {
    ForceAssertLt(alt1,4);
    ForceAssertLt(alt2,4);
    ForceAssertLt(alt3,4);
    ForceAssertLt(alt4,4);
    ForceAssertNe(alt2, alt1);
    ForceAssertNe(alt3, alt1);
    ForceAssertNe(alt4, alt1);
    ForceAssertNe(alt2, alt3);
    ForceAssertNe(alt2, alt4);
    ForceAssertNe(alt3, alt4);
    data = alt4;
    data <<= 2;
    data |= alt3;
    data <<= 2;
    data |= alt2;
    data <<= 2;
    data |= alt1;
  }

  /// sets alternative alleles for this SNP from a container holding 1-letter
  /// (character) representations of the allele bases. Container must be at most
  /// 4-elements long (4 alternative alleles max) and each string in the container
  /// must be 1 character long (single base). Like other Set... methods, if passed
  /// symbol is not a base or if same base is submitted twice, the method will abort.
  /// This method is safe as it overrides previous content, if any.
  void SetAlleles(vec<String> alleles) {
    ForceAssertLe(alleles.size(), 4u);
    ForceAssertGt(alleles.size(),0u);
    for ( unsigned int i = 0 ; i < alleles.size() ; i++ ) {
      ForceAssertEq(alleles[i].size(),1u);
      ForceAssert(Base::isCanonicalBase(alleles[i][0]));
    }

    switch (alleles.size() ) {
    case 1: SetAlleles( as_char(alleles[0][0]) ); break;
    case 2: SetAlleles( as_char(alleles[0][0]), as_char(alleles[1][0]) ); break;
    case 3: SetAlleles( as_char(alleles[0][0]), as_char(alleles[1][0]), as_char(alleles[2][0]) ); break;
    case 4: SetAlleles( as_char(alleles[0][0]), as_char(alleles[1][0]), as_char(alleles[2][0]), as_char(alleles[3][0])); break;
    default: cout << "Vector of alleles is too long" << endl; exit(1); // should never get here
    }
  }

  friend ostream & operator << (ostream & s, const SNPAlleles & a) {

    base_t alt1 = a.data & 0x3;
    s << as_base( alt1 ) ;
    base_t alt = ( a.data >> 2 ) & 0x3;
    if ( alt == alt1 ) return s;
    s << '/' << as_base( alt );

    alt = ( a.data >> 4 ) & 0x3;
    if ( alt == alt1 ) return s;
    s << '/' << as_base( alt );

    alt = ( a.data >> 6 ) & 0x3;
    if ( alt == alt1 ) return s;
    s << '/' << as_base( alt );
    return s;
  }

private:
  unsigned char data;
};

/// Wrapper class: collection of features mapped onto a genome. Main purpose is hiding implementation and
/// providing slightly more convenient interface; logically, it IS a map. The <Feature> must have default
/// constructor, and a single (always the same!!) instance of default-constructed "null feature" will be
/// returned when the map is queried at genome position where no feature was previously recorded: in contrast
/// to STL hash_map implementation, querying for non-existing key does not create a new entry. If default
/// constructor does not create an "empty" feature that is generically distinguishable from "real" feature,
/// then it will always look like Get() returns something, for any position on the genome,
/// so that HasFeature() query should be used first!

template <class Feature>
class FeatureMap {
public:
  /// ctor; takes the numbre of contigs in the genome.
  FeatureMap(unsigned int ncontigs) : null_feature() { Resize(ncontigs); }

  /// sets the number of contigs in the genome
  void Resize(unsigned int ncontigs) { storage.resize(ncontigs); }

  /// sets the number of contigs and number of hash buckets per contig
  void Resize(unsigned int ncontigs, unsigned int features_per_contig) {
    storage.resize(ncontigs);
    for ( unsigned int i = 0 ; i < storage.size() ; i++ ) storage[i].resize(features_per_contig);
  }

  /// sets the number of hash buckets per each contig
  void ResizeContigs(unsigned int features_per_contig) {
    for ( unsigned int i = 0 ; i < storage.size() ; i++ ) storage[i].resize(features_per_contig);
  }

  /// adds feature at the specified position
  void Add(unsigned int contig, unsigned int pos, const Feature & f) {
    storage[contig][pos]=f;
  }

  /// returns true if a feature exists at the specified position in the genome
  bool HasFeature(unsigned int contig, unsigned int pos) {
    return ( storage[contig].find(pos) != storage[contig].end() );
  }

  /// returns the feature at the specified position in the genome. If
  /// no feature exists at this position returns an instance (this is always
  /// the same single instance for the whole lifetime of the FeatureMap object!!)
  /// of default-constructed <Feature>. If there is no way to tell "empty"
  /// default-constructed <Feature> from a "real" one, use HasFeature() first
  /// to make sure that the feature at the specified position does exist.
  const Feature & Get(unsigned int contig, unsigned int pos) {
    typename storage_map::iterator it = storage[contig].find(pos);
    if ( it == storage[contig].end() ) return null_feature;
    else return it->second;
  }

protected:
  typedef __gnu_cxx::hash_map<unsigned,Feature,std::hash<unsigned>,std::equal_to<unsigned>,typename DefaultAllocator<Feature>::type> storage_map;
  vec< storage_map > storage;
  Feature null_feature;
};

/// Utility wrapper class: this is a FeatureMap of SNPAlleles with
/// a few utility methods added (loading from file, comparing to the
/// reference, etc)
class DbSNP : public FeatureMap<SNPAlleles> {
 public:

  DbSNP(unsigned int ncontigs) : FeatureMap<SNPAlleles>(ncontigs) {}

  enum MissingRefAction {
    ABORT = 0,  ///< abort if reference base is missing
    WARN = 0x1, ///< issue a warning if reference base is missing, don't abort, keep the SNP
    UPDATE = 0x2, ///< if reference base is missing, add it to the list of alleles, keep the SNP
    DISCARD=0x4 ///< discard SNP if reference base is missing from the list of alleles
  };

  /// Load dbSNP data from a flat file. Currently, the only recognized
  /// file format is: space or tab separated lines with the following fields
  /// contig pos id strand alleles.
  /// contig, pos are integers (position on the genome), SNP id is unused, strand
  /// must be "+" or "-" sign, and alleles is a "/" separated (no spaces!) list
  /// of alternative alleles, e.g. A/C. NOTE: position is always assumed to be given
  /// on the positive strand, the "strand" field affects only how alleles are interpreted
  /// (complement is taken for each allele base if strand=="-").
  void Load(const String & fname);

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
  void CheckAgainstRef(const vecbasevector & ref, int action = (WARN | UPDATE) );


};
