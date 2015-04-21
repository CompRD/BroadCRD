///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* 
   --- ExciseTransposons ---

   This program trims off those parts of reads which properly match a library of
   transposons (or other sequences).

   (The following information needs to be updated.)

   --- Motivation and problems ---

   DNA from bacterial (E. coli) transposons may appear in reads, in any of several 
   possible ways:

   [1]: The read may be an accidental all-bacterial read, which happens to
   include part of a transposon.

   [2]: There may be a bacterial transposon in the target genome.

   [3]: All copies of a BAC may have an inserted transposon.

   [4]: Some copies of a BAC may have an inserted transposon.

   It is not at all clear how this rather complicated situation should be dealt
   with.  ExciseTransposons is a first (and not very good) attempt.

   Here are two questions which should probably be answered before proceeding
   further:

   [a]: Can [2] really happen?  The official sequence for human genome 22 includes
   an E. coli transposon, but perhaps this was an error.  Are we confident that
   E. coli transposons could never occur in mammalian DNA?

   [b]: Life would be simpler if the E. coli strains in use had been modified so
   that their transposons were not viable.  It is presumably far too late to
   carry out such a modification.

   Note.  This program can generate zero-length reads.  It is not yet clear that
   they are properly handled.

*/

#include <stdlib.h>

#include "Alignment.h"
#include "system/Assert.h"
#include "Badness.h"
#include "Basevector.h"
#include "ExciseTransposons.h"
#include "FetchReads.h"
#include "math/Functions.h"
#include "pairwise_aligners/MakeAligns.h"
#include "Overlap.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "Quality.h"
#include "Qualvector.h"
#include "ShortVector.h"
#include "String.h"
#include "system/System.h"
#include "TaskTimer.h"
#include "TrimAlignmentEnds.h"
#include "Vec.h"

class nobbit {    

public:

  align a;
  int RC;
  int length1, length2;
  int id1, id2;
  int pos1, pos2;
  int Pos1, Pos2;

  nobbit() { }
  nobbit( align a_arg,
	  int RC_arg,
	  int length1_arg,
	  int length2_arg,
	  int id1_arg, 
          int id2_arg,
	  int pos1_arg,
	  int pos2_arg,
	  int Pos1_arg,
	  int Pos2_arg )
    : a(a_arg), 
      RC(RC_arg),
      length1(length1_arg),
      length2(length2_arg),
      id1(id1_arg),
      id2(id2_arg),
      pos1(pos1_arg),
      pos2(pos2_arg),
      Pos1(Pos1_arg),
      Pos2(Pos2_arg)
  { }
};

void TrimReadByMatchToContig( const basevector& read,
			      const String& readname,
                              const basevector& contig,
			      const String& contigname,
			      nobbit& the_nobbit,
			      const bool kill_all,
			      const int ml,
			      bool* trimmed_to_nothing,
			      short& left_trim,
			      short& right_trim,
			      ostream& log );

void ExciseTransposons( vecbasevector& reads,
			vecqualvector& Q,
                        const vecString& ids,
			String transposon_file,
			vec< vec<int>* >& killed_reads,
                        vec<int>& start,
			vec<int>& stop,
			ostream& log,
                        String run_dir,
			int ml,
			const vec<int>& min_overlap,
                        const vec<int>& min_perfect_overlap,
                        const vec<Bool>& kill_all, const Bool quiet )
{
  if ( !quiet ) cout << Date( ) << ": Excising transposons..." << endl;

  String contigs_file = transposon_file;

  // Read in the "transposon" contigs.

  const vecbasevector & EE = reads;

  // Read in the contigs and their names.

  vecbasevector contigs;
  FetchReads( contigs, 0, contigs_file );
  // Read in known contig names
  vec<String> contig_names;
  Ifstream( contigs_f, contigs_file );
  while(1)
  {    String s;
       getline( contigs_f, s );
       if ( !contigs_f ) break;
       if ( s.Contains( ">", 0 ) ) 
       {    s = s.After( ">" );
	    DeleteLeadingWhiteSpace(s);
	    DeleteTrailingWhiteSpace(s);
	    contig_names.push_back(s);    }    }
     
  int num_reads = EE.size(), num_contigs = contigs.size();
  ForceAssert( (int) min_overlap.size( ) == num_contigs );
  ForceAssert( (int) min_perfect_overlap.size( ) == num_contigs );
  ForceAssert( (int) kill_all.size( ) == num_contigs );
  ForceAssert( (int) killed_reads.size( ) == num_contigs );

  // Based on limited tests, the following value of block_size seems to be
  // about optimal.  In the tests, it was better than 18000, and about the
  // same as 13000.

  int block_size = 15000;

  if ( !quiet )
  {    cout << "Working on a set of " << num_reads << " reads " 
            << block_size << " at a time." << endl;    }

  temp_file aligns_file( run_dir + "/excision_aligns_file_temp_XXXXXXX" );

  vector<Bool> trimmed((size_type) num_reads, false);
  vector<short> left_trim((size_type) num_reads, 0);
  vector<short> right_trim((size_type) num_reads, 0);
  vector<int> who_trimmed((size_type) num_reads, 0);  // default 0 is lazy - should never be used

  int passes = (num_reads+block_size-1)/block_size;

  // Define chunk and precompute space to be reserved.

  vecbasevector chunk;
  longlong maxtotalbases = 0;
  int maxseqs = block_size + num_contigs;

  for ( int pass = 0; pass < passes; pass++ )
  {
    longlong totalbases = 0;
    int chunk_start = block_size * pass;
    int chunk_size = Min( block_size, num_reads - chunk_start );
    for ( int i = 0; i < chunk_size; i++ )
      totalbases += EE[ chunk_start + i ].size( );
    for ( int i = 0; i < num_contigs; i++ )
      totalbases += contigs[i].size( );
    maxtotalbases = Max( totalbases, maxtotalbases );    
  }
  chunk.Reserve( maxtotalbases/16 + maxseqs, maxseqs );
  chunk.resize(maxseqs);

  if ( !quiet ) cout << "Running " << passes << " passes:" << endl;
  TaskTimer timer;

  for ( int pass = 0; pass < passes; pass++ )
  {
    timer.Start();

    int chunk_start = block_size * pass;
    int chunk_size = Min( block_size, num_reads - chunk_start );
    chunk.resize( chunk_size + num_contigs );
    for ( int i = 0; i < chunk_size; i++ )
      chunk[i] = EE[ chunk_start + i ];
    for ( int i = 0; i < num_contigs; i++ )
      chunk[ i + chunk_size ] = contigs[ i ];

    // cout << "Pre MakeAligns(): " ;
    // PrintMemUsage();
	 
    // Note.  If you see vector which is not being trimmed, it may be
    // because the parameters in the following MakeAligns call are not
    // lenient enough.
	 
    MakeAligns<2, 24, 50>( 10, 10, chunk, 
          to_compare(FIRST_VS_SECOND, chunk_size), 1000, 1000, aligns_file, log, 
          200, 10000, 0, 100 );
	
    // cout << "Post MakeAligns(): " ;
    // PrintMemUsage();
	 
    int aligns_length;
    Ifstream( aligns_in, aligns_file );
    int id1, id2;
    Bool rc = False;

    int aligns_in_this_pass = 0;
	 
    align an_align;
    while(1)
    {
      BinRead( aligns_in, aligns_length );

      if ( !aligns_in ) 
	break;

      for ( int ll = 0; ll < aligns_length; ll++ )
      {
	aligns_in_this_pass++;

        int errors;
	an_align.Read( aligns_in, errors, id1, id2, rc );
	int length1 = chunk[id1].size( ), length2 = chunk[id2].size( );
	bool swapped = false;
		 
	if ( id1 >= id2 )
	{
	  swap( length1, length2 );
	  swap( id1, id2 );
	  swapped = true;    
	}

	if ( trimmed[id1+chunk_start] || EE[id1+chunk_start].size( ) == 0 ) 
             continue;
		 
	// We require an overlap of at least min_overlap bases.
		 
	if ( an_align.Pos2( ) - an_align.pos2( ) < min_overlap[id2 - chunk_size] ) 
	     continue;

        int pos1 = an_align.pos1( ), pos2 = an_align.pos2( );
        const avector<int> &gaps = an_align.Gaps( ), &lengths = an_align.Lengths( );
        int nblocks = an_align.Nblocks( );

	int Pos1 = pos1, Pos2 = pos2;
	for ( int jj = 0; jj < nblocks; jj++ )
	{
	  if ( gaps(jj) > 0 ) Pos2 += gaps(jj);
	  else if ( gaps(jj) < 0 ) Pos1 -= gaps(jj);
	  Pos1 += lengths(jj); Pos2 += lengths(jj);    
	}

	if (swapped)
	{
	  swap( pos1, pos2 );
	  swap( Pos1, Pos2 );    
	}

	nobbit the_nobbit;
	
	if ( !swapped )
	{
	  if ( !rc ) 
	  {
	    the_nobbit = nobbit( an_align, rc, length1, 
				 length2, id1+chunk_start, id2-chunk_size, 
				 pos1, pos2, Pos1, Pos2 );    
	  }
	  else 
	  {
	    align r;
            r = an_align;
	    r.ReverseThis( chunk[id1].size( ), 
			   chunk[id2].size( ) );
            int xpos1 = r.pos1( ), xpos2 = r.pos2( );
            int xPos1 = r.Pos1( ), xPos2 = r.Pos2( );
	    the_nobbit = nobbit( r, rc, length1, length2, 
				 id1+chunk_start, id2-chunk_size, xpos1, xpos2, 
				 xPos1, xPos2 );    
	  }
	}
      	else
	{
	  align a;
          a = an_align;
	  a.Flip( );
	  the_nobbit = nobbit( a, rc, length1, length2, 
			       id1+chunk_start, id2-chunk_size, 
			       pos1, pos2, Pos1, Pos2 );   
	}
	
	int read_id = the_nobbit.id1, contig_id = the_nobbit.id2;

	TrimAlignmentEnds(the_nobbit.a, the_nobbit.RC, EE[read_id], 
             contigs[contig_id]);
        the_nobbit.pos1 = the_nobbit.a.pos1( );
        the_nobbit.pos2 = the_nobbit.a.pos2( );
        the_nobbit.Pos1 = the_nobbit.a.Pos1( );
        the_nobbit.Pos2 = the_nobbit.a.Pos2( );
  
	if ( the_nobbit.Pos2 - the_nobbit.pos2 < min_overlap[contig_id] ) continue;
	if ( the_nobbit.Pos2 - the_nobbit.pos2 < min_perfect_overlap[contig_id] ) 
             continue;

        // We require a perfect match of a size specified by min_perfect_overlap.

        if ( min_perfect_overlap[contig_id] > 0 )
        {    int maxp = MaxPerfectMatch( the_nobbit.RC, the_nobbit.a, 
                  EE[read_id], contigs[contig_id] );
             if ( maxp < min_perfect_overlap[contig_id] ) continue;    }
	
	bool trimmed_to_nothing;
	TrimReadByMatchToContig( EE[read_id], ids[read_id], contigs[contig_id], 
                                 contig_names[contig_id],
				 the_nobbit, kill_all[contig_id], ml,
				 &trimmed_to_nothing,
				 left_trim[read_id], right_trim[read_id], log );
        who_trimmed[read_id] = contig_id;
	
	if ( trimmed_to_nothing )
	{
	  reads[read_id].Setsize(0);
	  Destroy(Q[read_id]);    
	  if ( !trimmed[read_id] ) killed_reads[contig_id]->push_back(read_id);
	  trimmed[read_id] = True;
	}
      }
    }
    
    //    Dot( cout, pass );
    if ( !quiet ) cout << Date() << ": pass " << pass << endl;
    timer.Stop();
    if ( !quiet ) cout << timer << endl;
    timer.Reset();
  }
     
  if ( !quiet ) cout << endl;

  if ( !quiet ) cout << Date() << ": tidying up all reads..." << endl;

  for ( int id1 = 0; id1 < num_reads; id1++ ) 
  {
    if ( trimmed[id1] == True ) 
      continue;
    
    int new_length = reads[id1].size( );
    new_length -= left_trim[id1] + right_trim[id1];
    if ( new_length < ml ) // if it's so trimmed, it's useless
    {    
      reads[id1].Setsize(0);
      Destroy(Q[id1]);

      // Note that the following way of determining the "killer" of id1 is
      // inaccurate, because it does not distinguish between the left-trimmer(s)
      // of id1, and the right-trimmer(s) of id1.  It would be best to declare that
      // the largest trimmer is the killer.  Instead, one of the trimmers is chosen
      // (essentially at random) and declared the killer.

      killed_reads[ who_trimmed[id1] ]->push_back(id1);

      continue; 
    }
    basevector b = reads[id1];
    if ( left_trim[id1] > 0 || right_trim[id1] > 0 )  // if we're trimming
    {
      reads[id1].SetToSubOf( b, left_trim[id1], new_length);
      start[id1] += left_trim[id1];
      stop[id1] -= right_trim[id1];
      if ( left_trim[id1] > 0 )
        for ( int j = 0; j < new_length; j++ )
  	  Q[id1][j] = Q[id1][j+left_trim[id1]];
      Q[id1].resize(new_length);
    }
  }
  if ( !quiet ) cout << Date() << ": tidying up all reads...done." << endl;

  if ( !quiet ) cout << Date( ) << ": Excising transposons...done." << endl;
  if ( !quiet ) PrintMemUsage();
}


void TrimReadByMatchToContig(const basevector& read,
			     const String& readname,
                             const basevector& contig,
			     const String& contigname,
			     nobbit& the_nobbit,
			     const bool kill_all,
			     const int ml,
			     bool* trimmed_to_nothing,
			     short& left_trim,
			     short& right_trim,
			     ostream& log)
{
  *trimmed_to_nothing = false;

  int RC = the_nobbit.RC;
  int id1 = the_nobbit.id1,   id2 = the_nobbit.id2;
  int pos1 = the_nobbit.pos1, Pos1 = the_nobbit.Pos1;
  int pos2 = the_nobbit.pos2, Pos2 = the_nobbit.Pos2;
  int length1 = the_nobbit.length1;

  // int pos1r = RC ? (Pos1 - 1) : pos1;
  // int Pos1r = RC ? pos1 : (Pos1 - 1);
  int pos1r = pos1;
  int Pos1r = Pos1;

  if (RC) 
  {    
    pos1r = length1 - Pos1;
    Pos1r = length1 - pos1;    
  }

  log << "\n***** " << (RC ? "rc to " : "") << "read " << id1
      << " (" << readname << ")" << " vs foreign sequence " << id2 
      << " (" << contigname << ")" << " *****\n\n";
  log << "start on read = " << pos1r << "\n";
  log << "length of read = " << read.size( ) << "\n";
  log << "overlap_length = " << Pos2 - pos2 << "\n";

  // The following logging is turned off because it generates monster log
  // files.  But if you are experimenting with the ExciseTransposons code,
  // you may wish to turn it back on.

  /*
  if ( !RC ) 
    PrintVisualAlignment( True, log, read, contig, a );
  else
  { 
    basevector b = read;
    b.ReverseComplement( );
    PrintVisualAlignment( True, log, b, contig, a );    
  }
  */

  if (kill_all)
  {    
    log << "kill_all set, so trimming entire read\n";
    *trimmed_to_nothing = true;
    return;
  }

  else if ( (pos1r <= 1 && Pos1r >= length1-2) ||
	    (pos1r <= 1 && length1 - Pos1r < ml) ||
	    (Pos1r >= length1-2 && pos1r < ml) )
  {
    log << "would trim entire read\n";
    *trimmed_to_nothing = true;
    return;
  }

  else if ( pos1r <= 1 )
  { 
    log << "would trim positions 0 through " << Pos1r - 1 << " on read\n";
    left_trim = Max( left_trim, (short) Pos1r );    
  }

  else if ( Pos1r >= length1-2 )
  {
    log << "would trim positions " << pos1r << " through end of read\n";
    right_trim = Max( right_trim, (short) (length1 - pos1r) );    
  }

  else 
  {
    log << "CONFUSING OVERLAP\n";
    if ( Pos1r <= left_trim )
      log << "would trim positions 0 through " << Pos1r - 1 
	  << " on read\n";
    else if ( length1 - pos1r <= right_trim )
      log << "would trim positions " << pos1r 
	  << " through end of read\n";
    else
    { 
      int if_left = Pos1r - left_trim;
      int if_right = length1 - pos1r - right_trim;
      if ( if_left <= if_right )
      {
	log << "would trim positions 0 through " << Pos1r - 1 
	    << " on read\n";
	left_trim = Pos1r;    
      }
      else
      {
	log << "would trim positions " << pos1r 
	    << " through end of read\n";
	right_trim = length1 - pos1r;   
      }
    }
  }    
}
