///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// EndSeq.  For NA12878 Fosmid pool data, find the sam records that bound ends of 
// Fosmid inserts and put them in a file endseqs.sam.  Because these sequences align
// poorly to the genome, they cannot all be obtained by directly searching the picard 
// bam file at the appropriate location on the genome.
//
// This program takes ~20 minutes to run.

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"

int main( )
{
     RunTime( );

     // The following two sequences should bound the ends of the Fosmid insert
     // as s1-genomic-s2'.

     basevector s1( "ACACGACGCTCTTCCGATCTAGTTGCTT" );
     basevector s2( "GACGTGTGCTCTTCCGATCTAGTTGCTT" );

     // Get their reverse complements.

     basevector s1rc(s1), s2rc(s2);
     s1rc.ReverseComplement( ), s2rc.ReverseComplement( );

     // Search the entire read set for these sequences.  Note copying of flowcell
     // etc. from LoadAndCorrect.cc.  On this first pass we get the pair ids.

     String picard = "/wga/scr4/picard";
     String flowcell = "A2925";
     String picard_run = "C1-508_2012-11-12_2012-11-14";
     String lane = "1";
     String lib = "Solexa-127359";
     String bam = picard + "/" + flowcell + "/" + picard_run + "/" + lane
          + "/" + lib + "/" + flowcell + "." + lane 
          + ".aligned.duplicates_marked.bam";
     String pattern = s1.ToString( ) + "|" + s1rc.ToString( ) + "|"
          + s2.ToString( ) + "|" + s2rc.ToString( );
     fast_pipe_ifstream in1( 
          "samtools view " + bam + " | egrep '" + pattern + "' | Col 1" );
     vec<String> ids;
     String line;
     while(1)
     {    getline( in1, line );
          if ( in1.fail( ) ) break;
          ids.push_back(line);    }
     UniqueSort(ids);

     // Now get all the reads that are in these pairs.

     vec<String> sams;
     fast_pipe_ifstream in2( "samtools view " + bam );
     while(1)
     {    getline( in2, line );
          if ( in2.fail( ) ) break;
          String id = line.Before( "\t" );
          if ( BinMember( ids, id ) ) sams.push_back(line);    }
     Sort(sams);
     
     // Output the lines.

     Ofstream( out, "/wga/scr4/human_data/NA12878_A2925/endseqs.sam" );
     for ( int j = 0; j < sams.isize( ); j++ )
          out << sams[j] << "\n";    }
