/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/** An AlignmentProfile is a record of where all the mutations,
    insertions, and deletions occur in an alignment.  Storage is
    very inefficient; you should never, for example, have one of
    these for every read sitting around.  The intended use is in
    creating simulated data which has the same error profiles as
    some experimentally-obtained data.

    Main Functionality
    1. Given an align of a read to a reference, create an AP
    2. Given an AP and a new reference, create a simulated read
       whose alignment to that ref matches the original align.

    \class AlignmentProfile
*/

#ifndef ALIGNMENTPROFILE_H
#define ALIGNMENTPROFILE_H

#include <string>
#include "PackAlign.h"

class AlignmentProfile {

  string m_s;

public:

  enum Status { OK = '.', MUTATION = '*', INSERTION = 'I', DELETION = 'D' };

  AlignmentProfile() {}
  
  /// Constructor using an align
  AlignmentProfile( const align& a, 
		    const basevector& read, 
		    const basevector& ref ) 
  { SetFromAlign(a, read, ref); }

  /// Constructor using a string
  explicit AlignmentProfile( const String& profile_string )
  { SetFromString( profile_string ); }


  /// Create AlignmentProfile from an align of a read to a reference
  void SetFromAlign( const align& a, 
		     const basevector& read, 
		     const basevector& ref );

  /// Create a simulated read using the AlignmentProfile and a new reference
  void CreateRead( const basevector& new_ref, int pos, 
		   basevector& new_read ) const;

  /// Create a simulated rc read using the AlignmentProfile and a new reference
  void CreateRcRead( const basevector& new_ref, int end_pos, 
		     basevector& new_read ) const;

  /// Trims the AlignmentProfile so that it will create reads of this size, if possible.
  /// Returns the new size of simulated reads produced by this AP. May not be requested
  /// size if original AP was too short or contained too many deletions.
  int TrimToReadSize(int size);

  /// Returns the size of the simulated read that will be created using this Alignment Profile
  int GetReadSize() const;

  /// Returns the number of bases on the reference required by this Alignment Profile
  int GetRefBaseCount() const;

  /// Returns number of base deletions, substitutions, insertions and matches.
  void GetStatistics(int& match, int& mut, int& ins, int& del) const;

  /// Returns total number of base deletions, substitutions and insertions.
  int GetErrorCount() const;

  /// Returns base positions of deletions, substitutions and insertions.
  void GetErrorPostions(vec<int>& mut, vec<int>& ins, vec<int>& del) const;

  /// Does the AlignmentProfile contain base insertions?
  bool HasInsertions() const {
    return m_s.find(INSERTION) != string::npos;
  }

  /// Does the AlignmentProfile contain base deletions?
  bool HasDeletions() const {
    return m_s.find(DELETION) != string::npos;
  }

  /// Does the AlignmentProfile contain base mutations?
  bool HasMutations () const {
    return m_s.find(MUTATION) != string::npos;
  }

  /// Does the AlignmentProfile contain base insertions, deletions and/or mutations?
  bool HasErrors () const {
    return m_s.find_first_not_of(OK) != string::npos;
  }


  /// Returns AlignmentProfile as a string
  const string& ProfileString() const { return m_s; }

  /// Creates AlignmentProfile from a string
  void SetFromString( const string& profile_string );


  /// Writes AlignmentProfile as string
  friend ostream& operator<<(ostream& out, const AlignmentProfile& prof)
  { return( out << prof.m_s ); }

  /// Reads AlignmentProfile from string
  friend istream& operator>>(istream& in, AlignmentProfile& prof)
  { 
    string s;
    in >> s;
    prof.SetFromString(s);
    return in;
  }


  friend bool operator==( const AlignmentProfile& lhs, 
			  const AlignmentProfile& rhs )
  { return( lhs.m_s == rhs.m_s ); }

  friend bool operator<( const AlignmentProfile& lhs, 
			 const AlignmentProfile& rhs )
  { return( lhs.m_s < rhs.m_s ); }

};

void WriteAll(String filename, const vec<AlignmentProfile>& profiles);

void ReadAll(String filename, vec<AlignmentProfile>& profiles);


#endif
