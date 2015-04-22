/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "AlignmentProfile.h"
#include "random/Random.h"

void AlignmentProfile::SetFromAlign( const align& a, 
				     const basevector& read,
				     const basevector& ref ) {
  int len=0;
  for(int j=0; j<a.Nblocks(); j++)
    len += a.Lengths(j) + abs(a.Gaps(j));
  m_s.reserve(len);
  
  int p1 = a.pos1(), p2 = a.pos2();
  for( int j = 0; j < a.Nblocks(); j++ ) {
    
    if ( a.Gaps(j) > 0 ) {
      m_s.resize( m_s.size() + a.Gaps(j), DELETION );
      p2 += a.Gaps(j);
    }
    else if ( a.Gaps(j) < 0 ) {
      m_s.resize( m_s.size() - a.Gaps(j), INSERTION );
      p1 -= a.Gaps(j);
    }

    for( int x = 0; x < a.Lengths(j); x++ )
      m_s += ( (read[p1++]==ref[p2++]) ? OK : MUTATION );
  }
}


void AlignmentProfile::CreateRead( const basevector& new_ref, int pos, 
				   basevector& new_read ) const {
  new_read.Setsize( m_s.size() );
  int readpos=0;
  
  for( unsigned int i=0; i<m_s.size(); i++)
    switch( m_s[i] ) {
    case OK:
      new_read.Set( readpos++, new_ref[pos++] );
      break;
    case MUTATION:
      new_read.Set( readpos++, (new_ref[pos++] + 1 + (randomx()%3)) % 4);
      break;
    case INSERTION:
      new_read.Set( readpos++, randomx() % 4);
      break;
    case DELETION:
      pos++;
      break;
    default:
      FatalErr("Corrupt value " << m_s[i] << " in AlignmentProfile!");
    }
  new_read.resize( readpos );
}

void AlignmentProfile::CreateRcRead( const basevector& new_ref, int end_pos, 
				   basevector& new_read ) const {
  new_read.Setsize( m_s.size() );
  int readpos=0;
  int pos = end_pos;
  
  for( unsigned int i=0; i<m_s.size(); i++)
    switch( m_s[i] ) {
    case OK:
      new_read.Set( readpos++, new_ref[pos--]^3);
      break;
    case MUTATION:
      new_read.Set( readpos++, ((new_ref[pos--]^3) + 1 + (randomx()%3)) % 4);
      break;
    case INSERTION:
      new_read.Set( readpos++, randomx() % 4);
      break;
    case DELETION:
      pos--;
      break;
    default:
      FatalErr("Corrupt value " << m_s[i] << " in AlignmentProfile!");
    }
  new_read.resize( readpos );
}


void AlignmentProfile::GetStatistics(int& match, int& mut, int& ins, int& del) const{
  match = mut = ins = del = 0;
  for( unsigned int i=0; i<m_s.size(); ++i) {
    switch( m_s[i] ) {
    case OK:
      ++match;
      break;
    case MUTATION:
      ++mut;
      break;
    case INSERTION:
      ++ins;
      break;
    case DELETION:
      ++del;
      break;
    }
  }
}


int AlignmentProfile::GetErrorCount() const{
  int errors = m_s.size();
  for( unsigned int i=0; i<m_s.size(); ++i)
    if ( m_s[i] == OK )
      --errors;
  return errors;
}


void AlignmentProfile::GetErrorPostions(vec<int>& mut, vec<int>& ins, vec<int>& del) const{
  mut.clear();  ins.clear();  del.clear();
  for( unsigned int i=0; i<m_s.size(); ++i) {
    switch( m_s[i] ) {
    case MUTATION:
      mut.push_back(i);
      break;
    case INSERTION:
      ins.push_back(i);
      break;
    case DELETION:
      del.push_back(i);
      break;
    }
  }
}


void AlignmentProfile::SetFromString( const string& profile_string ) {
  for(unsigned int i=0; i<profile_string.size(); i++)
    switch( profile_string[i] ) {
    case OK: case MUTATION: case INSERTION: case DELETION: break;
    default:
      FatalErr( "Bad character '" << profile_string[i]
		<< "' in AlignmentProfile string " 
		<< profile_string );
    }
  m_s = profile_string;
}


int AlignmentProfile::TrimToReadSize(int size) {
  int length = 0;
  for( unsigned int i=0; i<m_s.size(); ++i) {
    if (m_s[i] != DELETION)
      ++length;
    if (length == size) {
      m_s = m_s.substr(0, i+1);
      break;
    }
  }
  return length;
}


int AlignmentProfile::GetReadSize() const {
  int length = 0;
  for( unsigned int i=0; i<m_s.size(); i++) {
    if (m_s[i] != DELETION)
      ++length;
  }
  return length;
}

int AlignmentProfile::GetRefBaseCount() const {
  int length = 0;
  for( unsigned int i=0; i<m_s.size(); i++) {
    if (m_s[i] != INSERTION)
      ++length;
  }
  return length;
}

void WriteAll(String filename, const vec<AlignmentProfile>& profiles)
{
  Ofstream( out, filename );
  out << profiles.size() << "\n";
  for (int i = 0; i < profiles.isize(); ++i)
    out << profiles[i] << "\n";
  out.flush();
}

void ReadAll(String filename, vec<AlignmentProfile>& profiles)
{
  Ifstream( in, filename );
  int table_size;
  in >> table_size;
  profiles.clear();
  profiles.resize(table_size);
  for (int i = 0; i < table_size; ++i)
    in >> profiles[i];
}
