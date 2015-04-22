// Copyright (c) 2003 Broad Institute / Massachusetts Institute of Technology

// FastnToFast: Input a file_path.fastn file and output a file_path_untrimmed.fastb
//   and file_path_untrimmed.qualb (with ambiguous bases assigned a quality of 0 ).
//   Process BATCH_SIZE objects at time so don't run out of  memory.
//
// NOTE: make sure there is enough disc space since since writing out the new qualb and
//       fastb files probably will take up a lot of space!

#include "CompressedSequence.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "feudal/IncrementalWriter.h"

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(FILE);
  CommandArgument_String(OUTPUT_DIR);
  CommandArgument_UnsignedInt_OrDefault(BATCH_SIZE, 10000);
  EndCommandArguments;

  // see if output dir exists
  if ( !IsDirectory(OUTPUT_DIR) ) FatalErr("bad OUTPUT_DIR");

  String fastn_filename = FILE;
  ForceAssert( fastn_filename.Contains( "_orig.fastn", -1 ) );
  String qualb_filename = fastn_filename.Before ( ".fastn")  + "_fastn.qualb";
  ForceAssert( qualb_filename.Contains( "_orig_fastn.qualb", -1 ) );
  String filename = fastn_filename.Before ( ".fastn" );
  while( filename.Contains( "/" ) ) filename = filename.After( "/" );
  String new_qualb_filename =  OUTPUT_DIR + "/" + filename  + "_untrimmed.qualb";
  String fastb_filename = OUTPUT_DIR + "/" + filename + "_untrimmed.fastb";

  cout << "new qualb = " << new_qualb_filename << endl
       << "new fastb = " << fastb_filename << endl;

  IncrementalWriter<bvec> bWriter(fastb_filename.c_str());
  IncrementalWriter<qvec> qWriter(new_qualb_filename.c_str());
  veccompseq sequences;
  vecqualvector quals;
  vecbasevector bases;

  int fastn_objs = MastervecFileObjectCount( fastn_filename );
  int qualb_objs = MastervecFileObjectCount( qualb_filename );
  ForceAssertEq( fastn_objs, qualb_objs );
  cout << fastn_objs << " items" << endl;

  vec<char> temp_chars;
  for ( unsigned int start=0; start < static_cast<unsigned int>(fastn_objs); start += BATCH_SIZE )
  {
      unsigned int stop = start + BATCH_SIZE;
      if ( stop > static_cast<unsigned int>(fastn_objs) ) stop = fastn_objs;

      // clear since ReadRange appends
      sequences.clear();

      // read in [start, stop) objects
      sequences.ReadRange( fastn_filename, start, stop );
      quals.ReadRange( qualb_filename, start, stop );

      // reserve space for bases
      bases.reserve( quals.size() );

      for ( unsigned int seq_idx = 0; seq_idx < static_cast<unsigned int>(sequences.size()); ++seq_idx )
      {
          // bases, ambiguous ones are N
          temp_chars.clear();
          temp_chars = sequences[seq_idx].asVecChar();

          for ( unsigned int base_idx = 0; base_idx < temp_chars.size(); ++base_idx )
          {
              if ( temp_chars[base_idx] == 'N' )
	      {
	          //cout << seq_idx << " " << base_idx << endl;
                  quals[seq_idx][base_idx] = 0;
		}
          }

          // bases, ambiguous ones are randomly ACTG
          bases.push_back(sequences[seq_idx].asBasevector());
      }

      // write out bases and quals
      bWriter.add(bases.begin(),bases.end());
      qWriter.add(quals.begin(),quals.end());
      bases.clear();
      quals.clear();
      cout << "." << flush;
      if ( (start+1)%50 == 0) cout << endl;
    }

  bWriter.close();
  qWriter.close();
}
