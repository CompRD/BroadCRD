#include "system/ParsedArgs.h"
#include "system/RunTime.h"
#include "String.h"
#include "system/System.h"
#include "Vec.h"
#include "lookup/LookAlign.h"
#include <ext/hash_map>
using __gnu_cxx::hash_map;



// This program takes an alignment in maq format (produced by the mapview tool) and converts it to qltout.
// We use the comment (non-parseable) section to record extra information that the maq aligner produces, namely
// mapping quality scores. 
// This program produces two qltout files, one for each end.  In the case that one end has a placement and the other does nto, we write this
// alignment into the qltout for end 1. The problem here is that we don't know which lane the read came from. 
// We assume the the maq aligner was run in paired mode. Some columns are interrupted differently between single end and paired end modes. 

struct StringHashFunction
{
  size_t operator()(const String& s) const
  {
    __gnu_cxx::hash<const char*> h;
    return h(s.c_str());
  }
};

// Need to break this out into a separate module. 
// Make sure that the hash function returns the same value for fw and rc kmers.
typedef hash_map<String,int,StringHashFunction> chromNameMap;
// Mapping of maq names to Broad names (chrX, chr1, .. chr22, chrX, chrY, chr1_random, ..., chr22_random, chrX_random) excluding chr12_random and chr14_random
void InitializeChromNames(chromNameMap& names) {
     int i = 0;
     String chrPrefix = String("chr");
     String chrSuffix = String("_random");
     names["chrM"] = 0;
     for (int j = 1; j <= 22; j ++) {
       names[chrPrefix + ToString(j)] = i++;
     }
     names[chrPrefix + ToString("X")] = i++;
     names[chrPrefix + ToString("Y")] = i++;
     for (int j = 1; j <= 22; j ++) {
         if (j == 12 || j == 14 || j == 20) continue;
         names[chrPrefix + ToString(j) + chrSuffix] = i++;
     }
     names[chrPrefix + ToString("X") + chrSuffix] = i++;
}

#define MAQ_SMITH_WAT_CODE 130
int main( int argc, char *argv[] )
{
     RunTime();

     BeginCommandArguments;
     CommandArgument_String(REFERENCE);
     CommandArgument_String(MAQ_ALIGNS);
     CommandArgument_String(QLTOUT_END1);
     CommandArgument_String(QLTOUT_END2);
     EndCommandArguments;

     vecbasevector ref(REFERENCE);
     chromNameMap chrNames;

     InitializeChromNames(chrNames);

     Ifstream(fin, MAQ_ALIGNS);
     Ofstream(fout1, QLTOUT_END1);
     Ofstream(fout2, QLTOUT_END2);
     int prevIndexToRead = -1;
     while(!fin.eof()) {
          String line; 
          vec<String> tokens;
          getline(fin, line);
          Tokenize(line, tokens);
          // Tolerate blank lines.  All others are expected to include all columns.   
          if (tokens.size() == 0) continue;
            // A single line includes the following:
            //   [0] read name
            //   [1] chromosome
            //   [2] position
            //   [3] strand
            //   [4] insert size from the outer coordinates of a pair
            //   [5] paired flag
            //   [6] mapping quality
            //   [7] single-end mapping quality
            //   [8] alternative mapping quality
            //   [9] number of mismatches of the best hit
            //   [10] sum of qualities of mismatched bases of the best hit
            //   [11] number of 0-mismatch hits of the first 24bp
            //   [12] number of 1-mismatch hits of the first 24bp on the reference
            //   [13] length of the read
            //   [14] read sequence
            //   [15] read qualities
            // See the Maq manual for more detail
  	    int indexToRead = tokens[0].Int();
            unsigned int readLength = tokens[13].Int();
            unsigned int alignedFromReadBase = 0;               // Qltout coordinate system is 0-based...
            unsigned int alignedToReadBase = readLength;         // ... and open-ended.

            String strand = tokens[3];
            Bool orientationRc = strand.StartsWith("-");

            String chrom = tokens[1];
            int indexOfReferenceSequence = chrNames[chrom];
            int referenceSequenceStart = tokens[2].Int() - 1; // Maq *appears* to be 1-based
                                                               // (based on blatting representative sequences
                                                               // on the UCSC Genome browser)
            unsigned int referenceSequenceEnd = referenceSequenceStart + readLength;
            unsigned int lengthOfReferenceSeq = ref[indexOfReferenceSequence].isize();

            int numAlignmentBlocks = 1;  // Because Maq reports at most on gap
            
            int pairCode = tokens[5].Int();
            int gapStart = tokens[6].Int();
            int signedGapSize = tokens[7].Int();
            avector<int> gapSteps(0);
   	    avector<int> gapSizes(0);
            // Make block representing entire read. 
            // Add block for an indel, if any 
            if (pairCode == MAQ_SMITH_WAT_CODE  && gapStart) {
	        gapSteps.Append(0);
	        gapSteps.Append(-signedGapSize);
	        gapSizes.Append(gapStart - 1);
	        gapSizes.Append(readLength - (abs(signedGapSize) + gapStart));
                numAlignmentBlocks += 1;
            } else {
                gapSteps.Append(0);
                gapSizes.Append(readLength);
            }

            long lengthOfMatchingPortion = readLength;
            nmuts_t numMismatches = tokens[9].Int();
#if 0 
 look_align( int query_id_arg,
	      int target_id_arg,
	      unsigned int query_length_arg,
	      unsigned int target_length_arg,
	      Bool rc1_arg,
	      const align& a_arg,
	      int nhits_arg,
	      nmuts_t mutations_arg,
	      int indels_arg )
#endif
            look_align la(indexToRead,
                          indexOfReferenceSequence,
	                  readLength,
                          lengthOfReferenceSeq,
                          orientationRc,
            	          align(alignedFromReadBase, referenceSequenceStart, gapSteps, gapSizes),
                          0,
                          numMismatches,
                          numAlignmentBlocks);
            ostream& endStream = (prevIndexToRead == indexToRead ? fout2 : fout1);
            prevIndexToRead = indexToRead;
            la.PrintParseable(endStream);
            endStream << "PAIR CODE=" <<  pairCode << "\t" <<
	      "MAPPING_QUALITY=" << tokens[6].Int() << "\t" <<
	      "SINGLE_END_MAPPING_QUALITY=" << tokens[7] <<  "\t" <<
	      "ALTERNATIVE_MAPPING_QUALITY=" << tokens[8] << endl;

    }

    fin.close();
    fout1.close();
    fout2.close();

}
