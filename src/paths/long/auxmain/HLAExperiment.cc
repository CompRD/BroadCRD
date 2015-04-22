///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// HLAExperiment.  Try to find the 250 base reads in the HLA region
// chr6: 32.43 Mb - 32.63 Mb.  Use as 'bait' the corresponding two regions in the
// PNAS assembly, plus the reference sequence.

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "paths/long/MakeKmerStuff.h"

int main( )
{
     RunTime( );

     // Setup output.

     Ofstream( out, "bait.sam" );

     // From pnas assembly, get all the contigs in scaffold 603, plus all the contigs 
     // in scaffold 94 that are between contigs 8653 and 321866.

     String pnas = "/wga/scr4/ALLPATHS/H.sapiens.NA12878/fullhuman/fullhuman_v1/"
          "iainm.22jul/ASSEMBLIES/test/linear_scaffolds0.patched";
     vec<int> tigs = {1247,3810,5302,7483,8653,18368,18607,19903,24719,25495,28320,
          29710,30094,32073,34687,35214,41879,46603,48607,49876,50108,50275,58940,
          59838,59839,59840,61753,63945,75402,75938,78422,80655,83735,84193,85449,
          88869,93965,120230,126086,127132,143727,203490,205017,210648,222105,225347,
          248381,281464,289729,290847,317429,321866,329280};
     vecbasevector bait;
     bait.Read( pnas + ".contigs.fastb", tigs );

     // Append the region from the reference genome.

     vecbasevector chr;
     chr.ReadOne( "/wga/scr4/bigrefs/human19/genome.fastb", 5 );
     bait.push_back( chr[0], 32430000, 200000 );

     // Create list of kmers.

     vecbasevector bait2(bait);
     vecbasevector baitrc(bait);
     for ( int i = 0; i < (int) baitrc.size( ); i++ )
          baitrc[i].ReverseComplement( );
     bait2.Append(baitrc);
     const int K = 100;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup1( bait2, kmers_plus );
     vec< kmer<K> > kmers( kmers_plus.size( ) );
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          kmers[i] = kmers_plus[i].first;

     // Go through the reads.

     vec<String> hits;
     vec< kmer<K> > source;
     String dir = "/wga/scr4/picard/H01UJADXX/C1-508_2012-11-01_2012-11-04";
     for ( int lane = 1; lane <= 2; lane++ )
     {    String bam = dir + "/" + ToString(lane) + "/Solexa-125532/H01UJADXX."
               + ToString(lane) + ".aligned.duplicates_marked.bam";
          fast_pipe_ifstream in( "samtools view " + bam );
          String line;
          vec<char> sep = {'\t'};
          vec<String> tokens;
          kmer<K> x;
          int count = 0;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               TokenizeStrictly( line, sep, tokens );
               for ( int i = 0; i < tokens[9].isize( ); i++ )
                    if ( tokens[9][i] == 'N' ) tokens[9][i] = 'A';
               basevector b( tokens[9] );
               if ( ( tokens[1].Int( )/16 ) % 2 == 0 ) x.SetToSubOf( b, 0 );
               else x.SetToSubOf( b, b.isize( ) - K );
               if ( BinMember( kmers, x ) ) 
               {    hits.push_back( tokens[0] );    
                    source.push_back(x);    }
               count++;
               if ( count % 100000 == 0 ) DPRINT2( count, hits.size( ) );
                    }    }

     // Summarize what was found.

     cout << Date( ) << ": sorting" << endl;
     SortSync( source, hits );
     for ( int i = 0; i < source.isize( ); i++ )
     {    int j = source.NextDiff(i);
          PRINT2( i, j-i );
          i = j - 1;    }

     // Now go back and fetch all the relevant pairs.

     Sort(hits);
     for ( int lane = 1; lane <= 2; lane++ )
     {    String bam = dir + "/" + ToString(lane) + "/Solexa-125532/H01UJADXX."
               + ToString(lane) + ".aligned.duplicates_marked.bam";
          fast_pipe_ifstream in( "samtools view " + bam );
          String line;
          vec<char> sep = {'\t'};
          vec<String> tokens;
          kmer<K> x;
          int count = 0;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               TokenizeStrictly( line, sep, tokens );
               if ( BinMember( hits, tokens[0] ) ) out << line << "\n";    }    }

     // Make a proper bam file.

     SystemSucceed( "samtools view -H /wga/scr4/picard/H01UJADXX/C1-508_2012-11-01_2012-11-04/1/Solexa-125532/H01UJADXX.1.aligned.duplicates_marked.bam > head1.sam" );
     SystemSucceed( "samtools view -H /wga/scr4/picard/H01UJADXX/C1-508_2012-11-01_2012-11-04/2/Solexa-125532/H01UJADXX.2.aligned.duplicates_marked.bam > head2.sam" );
     SystemSucceed( "cat head1.sam head2.sam bait.sam > all.sam" );
     SystemSucceed( "samtools view -Sb all.sam > bait.bam" );
     SystemSucceed( "samtools sort bait.bam bait.sorted" );
     SystemSucceed( "samtools index bait.sorted.bam" );    }
