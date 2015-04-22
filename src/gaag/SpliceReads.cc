///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "util/RunCommand.h"
// MakeDepend: dependency FastqToFastbQualb

static inline 
String Tag(String S = "SR") { return Date() + " (" + S + "): "; } 

int min_overlap, min_score, min_score_delta, min_score_per_base;

/**
 * SpliceReads
 * 
 * Original documentation from Bruce
 *
 * SpliceReads does try to find the best overlapping alignment.  It
 * slides the reads over each other, scoring each overlap by summing
 * the quality of the matching bases (using the best of the qualities
 * between the reads) minus the quality of the mismatching bases
 * (using the lower of the qualities between the reads).
 *
 * It picks the alignment producing the best score, with a few
 * additional heuristic restrictions specified as arguments:
 * 
 * 1) MIN_OVERLAP: only alignments with at least this many overlapping
 * bases are considered (default 5)
 *
 * 2) MIN_SCORE: only alignments scoring at least this value are
 * considered (default 100, only comes into play for small-overlap
 * candidates).
 *
 * 3) MIN_SCORE_DELTA: the best alignment score must be at least this
 * much higher than the 2nd best score (default 100).  This is to
 * protect against over-extending alignments with simple sequence at
 * the ends (e.g., mono- or di-nuc repeats), when many overlaps
 * produce high scores.
 *
 * 4) MIN_SCORE_PER_BASE: the best score divided by the number of
 * overlapping bases must be at least this value (default 10).  Note
 * that the MIN_SCORE effectively raises this for small overlaps.
 *
 * Once the best (if any) alignment is selected, it merges the
 * overlapping portion of the reads using the base with the highest
 * quality between the two.  If the bases match, it assigns a base
 * quality which is the max quality between the two reads, and if
 * there was a mismatch, it uses the difference in quality for that
 * base.  That's pretty conservative while being mathematically
 * defensible; not sure if you are using the qualities or just the
 * bases.
 *
 * I'll be the first to admit that the heuristics above might not be
 * optimal; they are just what seemed to work reasonably well
 * initially.  I originally wrote this as an experiment for a
 * different purpose (trying to enhance error correction of
 * overlapping fragments for assembly), but the data I was using often
 * didn't have sufficient overlap, so I dropped it.  When Sarah told
 * me about this new datatype, it seemed to fit the need and work
 * reasonably well on the initial data she provided.  However, if you
 * have feedback, I'm open to new suggestions for better default
 * values or additional restrictions!
 *
 * It does not slide alignments beyond the read lengths; this is why
 * it will probably not find alignments for fragments which are
 * smaller than the read lengths.  As we told Dirk, this could be
 * handled with a bit of extra development work.
*/
void 
splice_pair( const bvec & r1, const qvec & q1,
	     const bvec & r2, const qvec & q2,
	     bvec & splice, qvec & qual,
	     u_int verbosity )
{
  splice.resize(0);
  
  size_t r1size = r1.size();
  size_t r2size = r2.size();
  size_t max_overlap = min(r1size, r2size);
  vec<int> scores(max_overlap, -1000000000);
  vec<int> indicies(max_overlap, -1);

  for (size_t overlap = min_overlap; overlap < max_overlap; ++overlap) {
    int score = 0;
    for (size_t i2 = 0; i2 < overlap; ++i2) {
      size_t i1 = r1size - overlap + i2;
      if (r1[i1] == r2[i2]) score += max(q1[i1], q2[i2]);
      else score -= min(q1[i1], q2[i2]);
    }
    indicies[overlap] = overlap;
    scores[overlap] = score;
    if (verbosity >= 2) {
      cout << "overlap " << overlap << " score " << score << endl;
    }
  }
  ReverseSortSync(scores, indicies);
  int dscore = scores[0] - scores[1];
  if (scores[0] >= min_score && scores[0] - scores[1] >= min_score_delta 
      && scores[0] / indicies[0] >= min_score_per_base) {
    size_t overlap = indicies[0];
    size_t newsize = r1size + r2size - overlap;
    splice.resize(newsize);
    qual.resize(newsize);
    for (u_int i = 0; i < r1size - overlap; ++i) {
      splice.set(i, r1[i]);
      qual.set(i, q1[i]);
    }
    for (size_t i2 = 0; i2 < overlap; ++i2) {
      size_t i1 = r1size - overlap + i2;
      int qi1 = q1[i1], qi2 = q2[i2];
      
      if (r1[i1] == r2[i2]) {
	splice.set(i1, r2[i2]);
	qual.set(i1, max(qi1, qi2));
      } else if (qi2 > qi1) {
	splice.set(i1, r2[i2]);
	qual.set(i1, qi2 - qi1);
      } else {
	splice.set(i1, r1[i1]);
	qual.set(i1, qi1 - qi2);
      }
    }
    for (u_int i2 = overlap; i2 < r2size; ++i2) {
      size_t i1 = r1size - overlap + i2;
      splice.set(i1, r2[i2]);
      qual.set(i1, q2[i2]);
    }

    if (verbosity >= 1) {
      PRINT4(indicies[0], scores[0], indicies[1], scores[1]);

      r1.PrintCol(cout, r1size);
      for (u_int i = 0; i < r1size - overlap; ++i) cout << " ";
      r2.PrintCol(cout, r2size);
      splice.PrintCol(cout, newsize);
    }
  }
}

int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;

  // File locations:
  CommandArgument_String_Doc(IN_HEAD, "looks for either IN_HEAD.{fastb,qualb} (if FASTQ=False, or IN_HEAD.{fastq.1,fastq.2} (if FASTQ=True)");
  CommandArgument_String_Doc(OUT_HEAD, "writes OUT_HEAD.{names,fastb,qualb.fastq}");

  // Fastq options:
  CommandArgument_Bool_OrDefault_Doc(FASTQ, False, "load and save .fastq files, rather than fastb/qualb files" );
  CommandArgument_UnsignedInt_OrDefault(PHRED_OFFSET, 33);
  
  // Heuristics:
  CommandArgument_UnsignedInt_OrDefault(MIN_OVERLAP, 5);
  CommandArgument_UnsignedInt_OrDefault(MIN_SCORE, 100);
  CommandArgument_UnsignedInt_OrDefault(MIN_SCORE_DELTA, 100);
  CommandArgument_UnsignedInt_OrDefault(MIN_SCORE_PER_BASE, 10);

  // Computational performance:
  CommandArgument_UnsignedInt_OrDefault(VERBOSITY, 0);

  EndCommandArguments;

  // Dir and file names.
  String fastq1_fn = IN_HEAD + ".fastq.1";
  String fastq2_fn = IN_HEAD + ".fastq.2";
  String fastb_in_fn = IN_HEAD + ".fastb";
  String qualb_in_fn = IN_HEAD + ".qualb";

  String tmp_head1 = OUT_HEAD + ".tmp1";
  String tmp_fastb1_fn = tmp_head1 + ".fastb";
  String tmp_qualb1_fn = tmp_head1 + ".qualb";

  String tmp_head2 = OUT_HEAD + ".tmp2";
  String tmp_fastb2_fn = tmp_head2 + ".fastb";
  String tmp_qualb2_fn = tmp_head2 + ".qualb";

  String names_out_fn = OUT_HEAD + ".names";
  String fastb_out_fn = OUT_HEAD + ".fastb";
  String qualb_out_fn = OUT_HEAD + ".qualb";
  String fastq_out_fn = OUT_HEAD + ".fastq";

  // Check args, and if files exist.
  if ( FASTQ ) {
    if ( ! ( IsRegularFile( fastq1_fn ) && IsRegularFile( fastq2_fn ) ) ) {
      cout << "Fatal error: .fastq.1 and/or .fastq.2 are missing.\n"
	   << "Maybe you want to run with FASTQ=False?\n" << endl;
      return 1;
    }
  }
  else {
    if ( ! ( IsRegularFile( fastb_in_fn ) && IsRegularFile( qualb_in_fn ) ) ) {
      cout << "Fatal error: .fastb and/or .qualb are missing.\n"
	   << "Maybe you want to run with FASTQ=True?\n" << endl;
      return 1;
    }
  }
  if ( PHRED_OFFSET != 33 && PHRED_OFFSET != 64 ) {
    cout << "Fatal error: PHRED_OFFSET may only be 33 or 64.\n" << endl;
    return 1;
  }
  
  // Global variables.
  min_overlap = MIN_OVERLAP;
  min_score = MIN_SCORE;
  min_score_delta = MIN_SCORE_DELTA;
  min_score_per_base = MIN_SCORE_PER_BASE;

  // Load.
  vecString names;   // names of pairs (with size = reads.size( ) / 2 )
  vecbvec reads;     // reads[2*i] and reads[1+2*1] are mates
  vecqvec quals;     // in sync with quals
  if ( ! FASTQ ) {
    cout << Tag() << "loading reads and quals" << endl;
    reads.ReadAll( fastb_in_fn );
    quals.ReadAll( qualb_in_fn );
    const size_t n_pairs = reads.size( ) / 2;
    names.reserve( n_pairs );
    for (size_t ii=0; ii<n_pairs; ii++)
      names.push_back( "read_" + ToString( ii ) );
  }
  else {
    cout << Tag( ) << "parsing fastq.1 and fastq.2" << endl;
    ifstream in1( fastq1_fn.c_str( ) );
    while ( in1 ) {
      String name;
      getline( in1, name );
      if ( ! in1 ) break;
      names.push_back( name );
      for (int ii=1; ii<4; ii++)
	getline( in1, name );
    }
    in1.close( );

    for (size_t ii=0; ii<names.size( ); ii++) {
      ForceAssert( names[ii].Contains( "/1", -1 ) );
      ForceAssert( names[ii].Contains( "@" ) );
      names[ii] = names[ii].Before( "/1" ).After( "@" );
    }
    
    cout << Tag( ) << "converting to fastb/qualb and reloading" << endl;
    String theCommand
      = "FastqToFastbQualb FASTQ=" + fastq1_fn
      + " OUT_HEAD=" + tmp_head1;
    RunCommandWithLog( theCommand, "/dev/null" );

    theCommand
      = "FastqToFastbQualb FASTQ=" + fastq2_fn
      + " OUT_HEAD=" + tmp_head2;
    RunCommandWithLog( theCommand, "/dev/null" );
    
    vecbvec bases1( tmp_fastb1_fn );
    vecbvec bases2( tmp_fastb2_fn );
    vecqvec quals1( tmp_qualb1_fn );
    vecqvec quals2( tmp_qualb2_fn );

    const size_t n_pairs = names.size( );
    const size_t n_reads = 2 * n_pairs;
    reads.reserve( n_reads );
    quals.reserve( n_reads );

    for (size_t ii=0; ii<n_pairs; ii++) {
      reads.push_back( bases1[ii] );
      reads.push_back( bases2[ii] );
      quals.push_back( quals1[ii] );
      quals.push_back( quals2[ii] );
    }
  }

  // Sanity check.
  ForceAssertEq( reads.size( ), quals.size( ) );
  ForceAssertEq( reads.size( ), 2 * names.size( ) );

  // Splice pairs.
  size_t dotter = 100000;
  cout << Tag() << "splicing "
       << ToStringAddCommas( names.size( ) ) << " pairs (.="
       << ToStringAddCommas( dotter ) << " pairs)"
       << endl;

  vecString spliced_names;
  vecbvec spliced_reads;
  vecqvec spliced_quals;
  for (size_t i = 0; i < reads.size(); i += 2) {
    if ( (i/2) % dotter == 0 ) Dot( cout, (i/2) / dotter );

    bvec splice;
    qvec qual;
    reads[i+1].ReverseComplement();
    quals[i+1].ReverseMe();
    splice_pair(reads[i], quals[i], reads[i+1], quals[i+1], splice, qual, VERBOSITY);
    if (splice.size()) {
      spliced_names.push_back( names[i/2] );
      spliced_reads.push_back(splice);
      spliced_quals.push_back(qual);
    }
  }
  cout << endl;

  // Report.
  size_t npairs = reads.size() / 2;
  size_t count = spliced_reads.size();
  cout << "\n"
       << ToStringAddCommas( count ) << " pairs spliced, out of "
       << ToStringAddCommas( npairs ) << " pairs in input ("
       << ToString( (100.0 * count) / npairs, 1 ) << "%)\n"
       << endl;

  // Save.
  cout << Tag( ) << "saving" << endl;
  spliced_names.WriteAll( names_out_fn );
  spliced_reads.WriteAll( fastb_out_fn );
  spliced_quals.WriteAll( qualb_out_fn );
  if ( FASTQ ) {
    ofstream outfq( fastq_out_fn.c_str( ) );

    for (size_t ii=0; ii<count; ii++) {
      String str_quals;
      size_t n_bases = spliced_quals[ii].size( );
      for (size_t jj=0; jj<n_bases; jj++) 
	str_quals += spliced_quals[ii][jj] + PHRED_OFFSET;
      
      outfq << "@" << spliced_names[ii] << "\n"
	    << spliced_reads[ii].ToString( ) << "\n"
	    << "+\n"
	    << str_quals << "\n";
    }
    
    outfq.close( );
  }
  
  cout << Tag( ) << "SpliceReads done" << endl;
  
}
