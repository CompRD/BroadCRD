/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#ifndef ERRORGENERATOR_H
#define ERRORGENERATOR_H

#include "system/Types.h"
#include "random/Random.h"
#include "AlignmentProfile.h"


/**
   Class: ErrorGenerator

   Uses an error template to produce random read alignment profiles
   and quality scores. These can then be used to generate simulated
   reads. Error generator tables come in 3 parts:
   tablename.error_index, tablename.error_profiles,
   tablename.error_quals These tables are produced by
   <MakeErrorGenerator>.

*/
class ErrorGenerator {

private:
  vecqualvector m_quals;
  vec<int> m_profile_index;
  vec<AlignmentProfile> m_profiles;
  longlong m_qual_table_size;

public:

  ErrorGenerator() {
    m_qual_table_size = 0;
  }

  /// General Constructor
  ErrorGenerator(const vecqualvector& quals, const vec<int>& profile_index,
		 const vec<AlignmentProfile>& profiles) : m_quals(quals),
							  m_profile_index(profile_index),
							  m_profiles(profiles) {
    m_qual_table_size = m_quals.size();
  }


  /// Constructor for special case error generator where all reads have no errors and Q=40
  ErrorGenerator(int num_bases) : m_quals(1),
				  m_profile_index(1,0),
				  m_profiles(1, AlignmentProfile(string(num_bases, '.'))) {
    m_quals[0] = qualvector(num_bases, 40);
    m_qual_table_size = m_quals.size();
  }


  /// Get random alignment profile and quality scores
  void GetErrorProfile(AlignmentProfile& profile, qualvector& quals) const {
    int qual_index = (randomx() % m_qual_table_size);
    quals = m_quals[qual_index];
    profile = m_profiles[m_profile_index[qual_index]];
  }

  /// Get random alignment profile only - no quality scores
  void GetErrorProfile(AlignmentProfile& profile) const {
    int qual_index = (randomx() % m_qual_table_size);
    profile = m_profiles[m_profile_index[qual_index]];
  }

  /// Get random alignment profile and quality scores and trim to required read size
  void GetErrorProfile(AlignmentProfile& profile, qualvector& quals, const int rsize) const {
    GetErrorProfile(profile, quals);
    quals.resize(rsize);
    profile.TrimToReadSize(rsize);
  }

  /// Import Error Generator Tables
  ErrorGenerator(string filename_head) {
    Read(filename_head);
  }

  int QualSize() const {
    return m_qual_table_size;
  }

  // TODO: Potentially dangerous truncation of quals size
  int ProfileSize() const {
    return m_quals.size();
  }

  /// What is the maximum reads size we can generate
  int MaxReadSize() const {
    if (m_profiles.size() != 0)
      return m_profiles[0].GetReadSize();
    else
      return 0;
  }

  // What is the minimum and maximum number of bases from the reference
  // we need to generate reads of the specified size using this generator
  void GetMinMaxRefBaseCount(int& min, int& max, int size = 0) const;

  void Write(const String& filename_head) const;

  void Read(const String& filename_head);

};

#endif
