/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// FastbSizes: print out the number of records, sizes (in bases) of the 
//	       records and average record size in a fastb file.
//
// SHOW_INDEX: print also the index (in the vector) of the item
// SUMMARY_ONLY: show only overall statistics (number of items, etc.)
// COUNT_SIZES: outputs a histogram of sizes.
// MAX_CONTIG_PRINT: limits the histogram to the first MAX_CONTIG_PRINT
//                   contigs.  Adds zero boxes padding if this is larger
//                   than the largest contig.

#include "Basevector.h"
#include "MainTools.h"
#include "math/Functions.h"
#include <fstream>

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(FILE);
     CommandArgument_String(REPORT);
     CommandArgument_Bool_OrDefault( SHOW_INDEX, False );
     CommandArgument_Bool_OrDefault( TOTAL_ONLY, False );
     CommandArgument_Bool_OrDefault( SUMMARY_ONLY, False );
     CommandArgument_Bool_OrDefault( TABULAR, False );
     CommandArgument_Bool_OrDefault( COUNT_SIZES, False );
     CommandArgument_UnsignedInt_OrDefault( MAX_CONTIG_PRINT, 0 );
     EndCommandArguments;

     if ( !FILE.Contains(".fastb") )
          cout << "Warning: " << FILE << " does not end in fastb." << endl;
     if ( !IsRegularFile(FILE) )
	FatalErr( "File " << FILE.c_str() << " doesn't exist" << endl );

     std::ofstream os(REPORT.c_str());

     vecbasevector b;
     b.ReadAll(FILE);

     os.setf(ios::internal);

     vec<int> sizes;
     vec<longlong> counts(1);
     longlong maxlength = 0;
     longlong sum = 0;
     for ( size_t i = 0; i < b.size( ); i++ )
     {    sum += b[i].size();
          if ( COUNT_SIZES ) {
             if ( b[i].size( ) > maxlength ) {
                maxlength = b[i].size( );
                counts.resize( maxlength+1 );
             }
             ++counts[ b[i].size( ) ];
          }
          if ( b[i].size( ) > 0 ) sizes.push_back( b[i].size( ) );    
     }
     Sort(sizes);

     if (TOTAL_ONLY)

       os << sum << endl;

     else { 

       if ( ! SUMMARY_ONLY ) {
	 for ( size_t i = 0; i < b.size( ); i++ ) {
	   if (SHOW_INDEX)
	     os << setw(8) << i << "   ";
	   
	   os << setw(10) << b[i].size( ) << "\n";
	 }
       }

       long double ave =  static_cast<long double>(sum) / b.size();
       
       os.setf(ios::fixed);
       if (TABULAR)
	 {
	   os << "number_of_items total_bases average" << endl;
	   os << b.size() << " " << sum << " " << setprecision(2) << ave << endl;
	 }
       else
	 {
	   os << "Number of items = " << b.size() << endl
		<< "Total Bases = " << sum << endl
		<< "Average = " << setprecision(2) << ave << "\n" 
		<< "N50 = " << N50(sizes) << "\n" << endl;
	 }
       if (COUNT_SIZES)
	 {
	   os << "Contig length      Count\n";
	   os << "-------------      -----\n";
	   if ( MAX_CONTIG_PRINT == 0 ) 
             MAX_CONTIG_PRINT = maxlength;
	   counts.resize( MAX_CONTIG_PRINT + 1);
	   for (unsigned int i = 0; i <= MAX_CONTIG_PRINT; ++i)
             os << "     " << i << "               " << counts[i] << endl;
	 }
     }
     os.close();
}
