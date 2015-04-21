// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

// Basic tools on twelvmeres (as convert twelvmere to integer).

#ifndef TWELVMERE_H
#define TWELVMERE_H



/*
 * TmerToInt
 *
 * Convert the twelvmere of contig beginning at start into a unique integer
 * identifier. The vector fourPow stores the first twelve powers of 4.
 */
int TmerToInt( const basevector &contig,
	       const vec<int> &fourPow,
	       int start )
{
  int result = 0;

  for (int ii=start; ii<start+12; ii++) {
      result += contig[ii] * fourPow[ii-start];
  }

  return result;
}



/*
 * IntToTmer
 *
 * Inverse to TmerToInt: it converts the given integer into the proper
 * twelvmere.
 */
void IntToTmer( int id,
		basevector &bases )
{
  bases.resize(12);

  for (int ii=0; ii<12; ii++) {
    bases.Set(ii,id&3);
    id >>= 2;
  }
}



#endif
