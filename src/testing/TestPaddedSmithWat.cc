// Copyright (c) 2003 Whitehead Institute for Biomedical Research
//

#include "MainTools.h"
#include "PaddedSmithWat.h"

int main( int argc, char** argv )
{
  BeginCommandArguments;
  CommandArgument_String( SEQ1 );
  CommandArgument_String( SEQ2 );
  CommandArgument_Bool_OrDefault( LEFT, False );
  CommandArgument_Bool_OrDefault( RIGHT, False );
  CommandArgument_UnsignedInt_OrDefault( BANDWIDTH, 10 );
  CommandArgument_UnsignedInt_OrDefault( OFFSET, 0 );
  CommandArgument_UnsignedInt_OrDefault( GAP_PENALTY, 3 );
  CommandArgument_UnsignedInt_OrDefault( MISMATCH_PENALTY, 2 );
  CommandArgument_Bool_OrDefault( NEG_OFFSET, False );
  EndCommandArguments;

  vec<char> bases1;
  for ( unsigned int b=0; b<SEQ1.size(); ++b )
    bases1.push_back( SEQ1[b] );

  vec<char> bases2;
  for ( unsigned int b=0; b<SEQ2.size(); ++b )
    bases2.push_back( SEQ2[b] );
  
  alignment a;

  int bandwidth = BANDWIDTH;
  int offset = OFFSET;
  if ( NEG_OFFSET )
    offset = -offset;

  float score = PaddedSmithWat( bases1, bases2, 
				offset, bandwidth, a,
                                MISMATCH_PENALTY, GAP_PENALTY );
  
  int pos1(0), pos2(0), errors(0);
  avector<int> gaps(0), lengths(0);

  a.Unpack(pos1, pos2, errors, gaps, lengths);

  cout << "seq1: ";
  copy( bases1.begin(), bases1.end(), ostream_iterator<char>( cout, "" ) );
  cout << endl;

  cout << "seq2: ";
  copy( bases2.begin(), bases2.end(), ostream_iterator<char>( cout, "" ) );
  cout << endl;

  cout << "pos1: " << pos1 << endl;
  cout << "pos2: " << pos2 << endl;

  cout << "gaps: ";
  for ( int i = 0; i < (int) gaps.length; ++i )
    cout << gaps(i) << " ";
  cout << endl;

  cout << "lengths: ";
  for ( int i = 0; i < (int) lengths.length; ++i )
    cout << lengths(i) << " ";
  cout << endl;
    
  for ( int j = 0; j < pos2 - pos1; ++j )
    cout << " ";
  for ( int j = 0; j < pos1; ++j )
    cout << bases1[j];
  int pos = pos1;
  for ( int i = 0; i < (int) gaps.length; ++i )
  {
    if ( gaps(i) > 0 )
      for ( int j = 0; j < gaps(i); ++j )
	cout << "_";
    else if ( gaps(i) < 0 )
      for ( int j = 0; j < -gaps(i); ++j )

	cout << bases1[pos++];
    for ( int j = 0; j < lengths(i); ++j )
      cout << bases1[pos++];
  }
  while ( pos < (int) bases1.size() )
    cout << bases1[pos++];
  cout << endl;

  for ( int j = 0; j < pos1 - pos2; ++j )
    cout << " ";
  for ( int j = 0; j < pos2; ++j )
    cout << bases2[j];
  pos = pos2;
  for ( int i = 0; i < (int) gaps.length; ++i )
  {
    if ( gaps(i) < 0 )
      for ( int j = 0; j < -gaps(i); ++j )
	cout << "_";
    else if ( gaps(i) > 0 )
      for ( int j = 0; j < gaps(i); ++j )
	cout << bases2[pos++];
    for ( int j = 0; j < (int) lengths(i); ++j )
      cout << bases2[pos++];
  }
  while ( pos < (int) bases2.size() )
    cout << bases2[pos++];
  cout << endl;

  cout << score << endl;
}
  
