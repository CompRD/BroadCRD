// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


#include "CompressedSequence.h"
#include "system/RunTime.h"

bool test( const vec<char> test_vec )
{
  CompressedSequence test_cs(test_vec);
  vec<char> result_vec = test_cs.asVecChar();
  basevector result_basevec = test_cs.asBasevector();

  if ( test_vec.size() > 0 )
    cout << "Compressed sequence is 50% the size of the original." << endl;

  if ( test_vec.size() != result_vec.size() ) {
    cout << "Original: " << test_vec << endl;
    cout << "Reproduced: " << result_vec << endl;
    return false;
  }

  for ( unsigned int i = 0; i < test_vec.size(); i++ ) {
    if ( toupper(test_vec[i]) != toupper(result_vec[i]) ) {
      cout << "Original: " << test_vec << endl;
      cout << "Reproduced: " << result_vec << endl;
      return false;
    }
  }

  String result_basevec_string = result_basevec.ToString();

  if ( test_vec.size() != result_basevec.size() ) {
    cout << "Original: " << test_vec << endl;
    cout << "Basevector produced: "
	 << result_basevec_string.size() << endl
	 << result_basevec_string << endl;
    return false;
  }

  for ( unsigned int i = 0; i < test_vec.size(); i++ ) {
    if ( toupper(test_vec[i]) != 'N' &&
	 toupper(test_vec[i]) != toupper(result_basevec_string[i]) ) {
      cout << "Original: " << test_vec << endl;
      cout << "Basevector produced: "
	   << result_basevec_string.size() << endl
	   << result_basevec_string << endl;
      return false;
    }
  }

  return true;
}

int main(int argc, char** argv)
{
  RunTime(0);

  cout << "Testing zero-length sequence... " << endl;
  vec<char> test_vec_1;

  if ( test(test_vec_1) )
    cout << "passed." << endl;
  else
    cout << "FAILED!" << endl;

  char upper_chars[6] = "ACTGN";
  for ( int test_char = 0; test_char < 5; test_char++ ) {

    cout << "Testing sequence '" << upper_chars[test_char] << "'... " << flush;
    vec<char> test_vec_2;
    test_vec_2.push_back(upper_chars[test_char]);

    if ( test(test_vec_2) )
      cout << "passed." << endl;
    else
      cout << "FAILED!" << endl;
  }

  char lower_chars[6] = "actgn";
  for ( int test_char = 0; test_char < 5; test_char++ ) {

    cout << "Testing sequence '" << lower_chars[test_char] << "'... " << flush;
    vec<char> test_vec_2;
    test_vec_2.push_back(lower_chars[test_char]);

    if ( test(test_vec_2) )
      cout << "passed." << endl;
    else
      cout << "FAILED!" << endl;
  }

  cout << "Testing length 2 sequence... " << flush;
  vec<char> test_vec_3;
  test_vec_3.push_back('A');
  test_vec_3.push_back('C');

  if ( test(test_vec_3) )
    cout << "passed." << endl;
  else
    cout << "FAILED!" << endl;

  cout << "Testing length 3 sequence... " << flush;
  test_vec_3.push_back('N');

  if ( test(test_vec_3) )
    cout << "passed." << endl;
  else
    cout << "FAILED!" << endl;

  cout << "Testing length 4 sequence with lower... " << flush;
  test_vec_3.push_back('t');

  if ( test(test_vec_3) )
    cout << "passed." << endl;
  else
    cout << "FAILED!" << endl;

  cout << "Testing length 5 sequence... " << flush;
  test_vec_3.push_back('C');

  if ( test(test_vec_3) )
    cout << "passed." << endl;
  else
    cout << "FAILED!" << endl;

  cout << "Testing length 6 sequence... " << flush;
  test_vec_3.push_back('G');

  if ( test(test_vec_3) )
    cout << "passed." << endl;
  else
    cout << "FAILED!" << endl;

  char const* large_sequence = "NGCATAGCGTTGTGCCTGCGGTCGACTCTAGAGGATCCCCTGAGTTTGCCAGCCAACTATTCTGGATACCTGTGTCTCCAACTGAAGAAACCAAAAGCAGGATTACTACACTACCTAGGAATTTTGACTCTGTTGACCAGGAGCAAATCAAATAGTTACTTACACTAGAGGGAAATAATAAACTTATGTTTGGAATAAGTTTGTGTTAAGCTGTCCCACTTGGTACTACAGTAGCTTATAACTAGCATGGCTGGGAATTACCACATTTAGCATAAGTAGAAATGACTCTTAATTAAGCCTTCCTTAATTAAGAAGAAAACATACTGTAGTCAGTATTTTGAGTTCAGGACCAAATTGAAATTTGAACTTCTGAAGGATTGCTATAAATTTTGAACTGGACTGGCTAGCCTTACCCTTATGTAAACAAATCTCTTGTAATACATTTCTCAATATCTATTCCTGACTATTTTTGTTTCTCAGGACTTTTATTGACACATGATACATGATTCAATATATTTAAGGAGAGATGAATTTTTAAGAATCATGAATAAGGTGAAAATTATCCTGTGAGAATTCAAT";
  int large_seq_size = strlen(large_sequence);

  test_vec_3.clear();
  for ( int i = 0; i < large_seq_size; i++ )
    test_vec_3.push_back(large_sequence[i]);

  cout << "Testing large sequence..." << flush;
  if ( test(test_vec_3) )
    cout << "passed." << endl;
  else
    cout << "FAILED!" << endl;

  cout << "Testing copy constructor... " << flush;
  CompressedSequence test_cs_1( test_vec_3 );
  CompressedSequence test_cs_2( test_cs_1 );

  if ( test_cs_1 == test_cs_2 )
    cout << "passed." << endl;
  else
    cout << "FAILED!" << endl;

  cout << "Testing operator=... " << flush;
  CompressedSequence test_cs_3 = test_cs_1;

  if ( test_cs_1 == test_cs_3 )
    cout << "passed." << endl;
  else
    cout << "FAILED!" << endl;

  READ( "cs_test.data", vec<String>, bases );

  cout << "Testing const char* constructor... " << flush;
  CompressedSequence test_cs_4( bases[0].c_str() );

  vec<char> test_vec_4 = test_cs_4.asVecChar();
  unsigned int i;
  for ( i = 0; i < bases[0].size(); ++i )
  {
    if ( bases[0][i] != test_vec_4[i] )
    {
      cout << "FAILED! at element " << i << "." << endl;
      break;
    }
  }
  if ( i < bases[0].size() )
  {
    cout << "String:  " << bases[0] << endl;
    cout << "VecChar: " << test_vec_4 << endl;
  }
  else
  {
    cout << "passed." << endl;
  }

  cout << "Testing const char* constructor again... " << flush;

  veccompseq cs_mv_orig;
  int total_bases = 0;
  for ( unsigned int i = 0; i < bases.size(); ++i )
  {
    cs_mv_orig.push_back( CompressedSequence( bases[i].c_str() ) );
    total_bases += bases[i].size();
  }

  bool passed = true;
  for ( unsigned int i = 0; i < bases.size(); ++i )
  {
    if ( cs_mv_orig[i].asString() != bases[i] )
    {
      passed = false;
      cout << "FAILED! at element " << i << "." << endl;
      cout << "String:   " << bases[i] << endl;
      cout << "Original: " << cs_mv_orig[i].asString() << endl;
    }
  }

  if ( passed )
    cout << "passed." << endl;

  cout << "Testing mastervec read and write... " << flush;

  cs_mv_orig.WriteAll( "cs_test.fastn" );

  veccompseq cs_mv_in;
  cs_mv_in.Reserve( total_bases/5+bases.size(), bases.size() );
  cs_mv_in.ReadAll( "cs_test.fastn" );

  passed = true;
  for ( unsigned int i = 0; i < bases.size(); ++i )
  {
    if ( cs_mv_in[i] != cs_mv_orig[i] ||
	 cs_mv_in[i].asString() != bases[i] )
    {
      passed = false;
      cout << "FAILED! at element " << i << "." << endl;
      cout << "String:   " << bases[i] << endl;
      cout << "Original: " << cs_mv_orig[i].asString() << endl;
      cout << "Read in:  " << cs_mv_in[i].asString() << endl;
    }
  }

  if ( passed )
    cout << "passed." << endl;

  cout << "Testing mastervec sparse read... " << flush;

  cs_mv_in.clear();
  vec<int> subset;
  subset.push_back(1);
  subset.push_back(4);
  subset.push_back(5);
  subset.push_back(9);

  cs_mv_in.SparseRead( "cs_test.fastn", subset, 0 );

  passed = true;
  for ( unsigned int i = 0; i < subset.size(); ++i )
  {
    if ( cs_mv_in[subset[i]] != cs_mv_orig[subset[i]] ||
	 cs_mv_in[subset[i]].asString() != bases[subset[i]] )
    {
      passed = false;
      cout << "FAILED! at element " << i << "." << endl;
      cout << "String:   " << bases[subset[i]] << endl;
      cout << "Original: " << cs_mv_orig[subset[i]].asString() << endl;
      cout << "Read in:  " << cs_mv_in[subset[i]].asString() << endl;
    }
  }

  if ( passed )
    cout << "passed." << endl;

  cout << "Testing asBasevector()... " << flush;
  basevector b4 = test_cs_4.asBasevector();
  String test_cs_4_string = test_cs_4.asString();
  String b4_string = b4.ToString();

  for ( unsigned int i = 0; i < b4_string.size(); ++i )
    if ( test_cs_4_string[i] == 'N' )
      b4_string[i] = 'N';

  if ( b4_string != test_cs_4_string )
  {
    cout << "FAILED!" << endl;
    cout << "basevector: " << b4_string << endl;
    cout << "comp seq:   " << test_cs_4_string << endl;
  }
  else
    cout << "passed." << endl;


  char short_sequence[6];

  for ( int testCase = 0; testCase < 5; ++testCase )
  {
    switch ( testCase )
    {
      case 0:
        short_sequence[0] = 'N';
        break;
      case 1:
        short_sequence[1] = 'A';
        break;
      case 2:
        short_sequence[2] = 'C';
        break;
      case 3:
        short_sequence[3] = 'G';
        break;
      case 4:
        short_sequence[4] = 'T';
        break;
    }
    short_sequence[ testCase+1 ] = 0;

    cout << "Testing reverse complementation... " << flush;

    CompressedSequence test_cs_5( short_sequence );

    basevector b5 = test_cs_5.asBasevector();
    b5.ReverseComplement();
    test_cs_5.ReverseComplement();
    String test_cs_5_string = test_cs_5.asString();
    String b5_string = b5.ToString();

    for ( unsigned int i = 0; i < b5_string.size(); ++i )
      if ( test_cs_5_string[i] == 'N' )
        b5_string[i] = 'N';

    if ( b5_string != test_cs_5_string )
    {
      cout << "FAILED!" << endl;
      cout << "original:   " << short_sequence << endl;
      cout << "basevector: " << b5_string << endl;
      cout << "comp seq:   " << test_cs_5_string << endl;
    }
    else
      cout << "passed." << endl;
  }

  cout << "Testing reverse complementation... " << flush;

  char const *long_sequence=
"TTCTGTTTCTATTTTGTGGTTACTTTGAGGAGAGTTGGAATTAGGTCTTCTTTGAAGGTCTGGTAGAACTCTGCATTAAA"
"CCCATCTGGTCCTGGGCTTTTTTTTTTTTTTTTTTTTTTTTTTGGGTGGGAGACTATTGATGACTGCCTCTATTTCTTTA"
"GGGGAAATGGGACTTTTAGTCCATGAATCTGATCCTGATTTAGCTTTGGTACCTGGTATCTGTCTAGGAAGTTGTCCATT"
"TCATCCAGGTTTTCCTGGTTTTTTTTTAGTATAGCCTTTCATAGTAAAATCTGATGATGTTTTTGATATCCTCATGTTCT"
"GTTGGTATGTCTCCTTTTTCATTTCTGATTTTGTTAATTATAGTACAGTCCCTATGCCCTCTAGTTAGTCTGGCTAAGGG"
"TTTATCTATCTTGTTGACTTTCTCAAAGAACCAGCTACTATTTTGGTTGATTCTTTGAATATTTCTTTTTGTTTCCACTT"
"GGTTGATTTCAGCTCTGAGTTTGATTATTTCCTGCTGTCTACTCATCTTGGGTGAATTTGCTTCCTTTTGTTCTAGAGCT"
"TCTAGATTTGCTGTCAGGCTGCTAGTGTATACTCTAGTTTCCTTTTGGAGGCACACAGGCCTGTGAGTTTTACTCTTAGG"
"ACTGCCTCATTGTGCCCCATATGTTTGGCTATGTTGTGGATTTATTTTCATTAAACTTTAAAACATCTTTAATTTTTTTC"
"TTTATTTCATCATTGACCAAGCTATCATTAAGTAGAGTATTGTTCAGTTTCCAGGTGAATGTTGGCTTTCTATTATTTAT"
"GCTGTTATTGAAGATCAGCCTTAGTCCGTAGTTATCTGAAAGGATGCATGGGAAAATTTCAATATTTTTGTATCTGTTGA"
"GGACTTTTTGTGAGTGACTATATGGTCAATTTTGGAGGATTTGGTACTGAGAAGAAGGTATATATCCTTTTGTCTTATGA"
"TAAAATGTTCTGTAGATATCTATTAAATTCATTTGTTTCATAACTTCGGTTAGTGTCCGTGTGTCTCTGTTCAGTTTCTG"
"CTTCCAGGATCTGTCTCTTGGTGAGAGTGTGGTCTTGAAGTCTCCCAGTATTATTTTATGAGGTGCAATGTGTGCTTTGA"
"TCTTTAGCAAAGTGTATTTAATGAATGTGGCTGCTCTTGCATTTAGAGCATAGACATTCAGAATTGAGAGGTCATCTTGG"
"TAGATTTTGCCTTTGATGAGTATGAAGTGCCCCTCATTTTTTTTTTGATAACTTTGAGCTGAAAGTAAATTTGGTTCGAT"
"ATTAGAATTGCTACTCCAGCTTAATTCTTCATCCCATTTGCTTGGAAATTTGTTTTCCAACCTTTAACTCTGAGGTAGTG"
"TCTGTCTTTTTCCCTGAGGTGGGTTTCCTGTAAGCAACAAAATGTTGGGTCCTGTTTGTGTAGCCAGTCTATTAGTCTAT"
"GTCTTTTTACTGGGGAATTGAGTCCATTGATGTTAAGAGAAATTAAAGACAAGTAATTGTTGCTTCCTGTTATTTTTGTT"
"GTTAAATTTGGGATTCTGTTCTTGAGGCTGTCTTCTTTTAGGTTTGTTAAAGGATTGCATTCTTGCTTTTTCTAAGGTGT"
"AGTTTCTGTCCTTGTGTTTGTGTTATCCCGTTATTATCCTTTGAAAGGCTAGATTCATGGAAAGATATTGTGTGAATTTT"
"GTTTGGTCATGGAATACTTTGGTTTCTCCATCTATGGTAATTGAGAGTTTGGCTGGGTATAGCAACCTGGGCTGGCACTT"
"CTGCTTTCTTAGGGTCTGTATAACATCTGTCCAGGATCTTCTGGCTTTCATAGTCTCTGGTGAGAAGTCTGGTGTAATTC"
"TAATAGGCCTGCATTTATATGTTACTTGACCTTTTTCCCTTACTGCTAAAGATATTCTATCTTTATTTAGTGAATTTGTT"
"GTTCTGATTATTATGTGTCAGGAGGAATTTCTTTTCTGGTCCAATCTATTTGGAATTCTGTAGGCTTCTTTTATGTTCAT"
"GGATATCTCTTTCTTTAAGTTTGGGAAGTTTTCTTCTCTAATTTTGTTAAAGATATTTGCTGGTCCTTTAAGTTGAAAAT"
"CTTCATTCTCACCTACTCTTGTTGTCCATATGTTTGGGCTTTTCATTGC"
;
  CompressedSequence test_cs_6 = CompressedSequence( long_sequence );
  basevector b6 = test_cs_6.asBasevector();

  b6.ReverseComplement();
  test_cs_6.ReverseComplement();

  String test_cs_6_string = test_cs_6.asString();
  String b6_string = b6.ToString();

  for ( unsigned int i = 0; i < b6_string.size(); ++i )
    if ( test_cs_6_string[i] == 'N' )
      b6_string[i] = 'N';

  if ( b6_string != test_cs_6_string )
  {
    cout << "FAILED!" << endl;
    cout << "original:   " << long_sequence << endl;
    cout << "basevector: " << b6_string << endl;
    cout << "comp seq:   " << test_cs_6_string << endl;
  }
  else
    cout << "passed." << endl;

  cout << "Testing basevector ctor... " << flush;

  long_sequence = "ACGTACGTACTTTACGCCCTAACGTCCACGTACGT";
  basevector b7;
  b7.SetFromString( long_sequence );
  CompressedSequence test_cs_7( b7 );

  if ( b7 != test_cs_7.asBasevector() )
  {
    cout << "FAILED!" << endl;
    cout << "original:   " << long_sequence << endl;
    cout << "basevector: " << b7.ToString() << endl;
    cout << "comp seq:   " << test_cs_7.asString() << endl;
  }
  else
    cout << "passed." << endl;


  cout << "Done." << endl;
}
