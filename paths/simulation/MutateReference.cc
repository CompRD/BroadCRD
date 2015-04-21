/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: MutateReference

   Simulate a diploid genome from a haploid genome, 
   by introducing "polymorphic" mutations into given sequences.
   A diploid genome has two copies of each chromosome which are mostly the same
   except for some differences; this program takes a single copy of each chromosome
   (or of a subset of chromosome regions) and creates a slightly mutated copy of each
   chromosome or chromosome region, which we pretend is the second copy of that
   chromosome/chromosome-region in the now-diploid genome.
   
   Create as output a fastb file containing both the original sequences and
   the mutated ones.  Simplistic - to be done better later.

   See also <SimulateReads>, <SimulatePolymorphicReads>.

   @file
*/

#include "Basevector.h"
#include "MainTools.h"
#include "Vec.h"
#include "random/Random.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandDoc("Creates an artifical diploid genome fastb from a haploid one.");

     CommandArgument_String_Doc(FASTB_IN, "Haploid genome fastb file.");
     CommandArgument_String_Doc(FASTB_OUT, "Diploid genome fastb file.");
     CommandArgument_Double_Doc(SUB_RATE, "Per base mutation rate.");
     EndCommandArguments;

     cout << "Loading haploid genome file:\n" << FASTB_IN << endl;

     vecbasevector b1(FASTB_IN);
     vecbasevector b2(b1);

     cout << "Haploid genome size: " << b1.sumSizes() << endl;

     cout << "Creating SNPs (" << SUB_RATE * 100.0 
	  << "% polymorphic = 1 SNP per " << 1/SUB_RATE << " bases)" << endl;
     // TODO: potentially dangerous truncation of index by snp_locs
     vec< pair<int,int> > snp_locs;
     
     int M = 1000000;
     int m = int( round( double(M) * SUB_RATE ) );
     for ( size_t i = 0; i < b2.size( ); i++ )
       for ( int j = 0; j < b2[i].isize( ); j++ )
         if ( randomx( ) % M <= m ) {
           b2[i].Set( j, ( b2[i][j] + randomx( ) % 3 + 1 ) % 4 );
           snp_locs.push_back( make_pair( i, j ) );
           snp_locs.push_back( make_pair( b1.size() + i, j ) );
         }
     
     b1.Append(b2);

     cout << "Created " << snp_locs.isize()/2 << " SNPs." << endl;

     cout << "Saving diploid genome file:\n" << FASTB_OUT << endl;

     b1.WriteAll(FASTB_OUT);    

     cout << "Saving SNP location file file:\n" << FASTB_OUT + ".snp_locs" << endl;
    
     sort( snp_locs.begin(), snp_locs.end() );
     BinaryWriter::writeFile( FASTB_OUT + ".snp_locs", snp_locs );
}
