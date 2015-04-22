/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "paths/simulation/ErrorGenerator.h"
#include "feudal/BinaryStream.h"

void ErrorGenerator::GetMinMaxRefBaseCount(int& min, int& max, int size) const {
  ForceAssertLe(size, MaxReadSize());
  if (size == 0)
    size = MaxReadSize();

  min = max = size;
  if (m_profiles.size() != 0)
    for (int i = 0; i < m_profiles.isize(); ++i) {
      AlignmentProfile ap = m_profiles[i];
      ap.TrimToReadSize(size);
      int ref_bases = ap.GetRefBaseCount();
      min = MIN(min, ref_bases);
      max = MAX(max, ref_bases);
    }
}

void ErrorGenerator::Write(const String& filename_head) const {
  String profile_file = filename_head + ".error_profiles";
  String quals_file = filename_head + ".error_quals";
  String index_file = filename_head + ".error_index";
  WriteAll(profile_file, m_profiles);
  m_quals.WriteAll(quals_file);
  BinaryWriter::writeFile(index_file, m_profile_index);
}

void ErrorGenerator::Read(const String& filename_head) {
  String profile_file = filename_head + ".error_profiles";
  String quals_file = filename_head + ".error_quals";
  String index_file = filename_head + ".error_index";
  ReadAll(profile_file, m_profiles);
  m_quals.ReadAll(quals_file);
  BinaryReader::readFile(index_file, &m_profile_index);
  m_qual_table_size = m_quals.size();
}
