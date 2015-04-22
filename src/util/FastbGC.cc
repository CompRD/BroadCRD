/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// FastbGC: print out the records-by-record and average GC content in a fastb file.

#include "Basevector.h"
#include "MainTools.h"

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(FILE);
  CommandArgument_Bool_OrDefault( SHOW_INDEX, False );
  CommandArgument_Bool_OrDefault( COMPUTE_AVERAGE_PER_POSITION, False );
  // If positive, remove this many bases from each end of each sequence
  CommandArgument_UnsignedInt_OrDefault( TRIM, 0 );
  EndCommandArguments;

  if ( !FILE.Contains(".fastb") )
    FatalErr( "File " << FILE.c_str() << " must end in .fastb"  << endl);
  if ( !IsRegularFile(FILE) )
    FatalErr( "File " << FILE.c_str() << " doesn't exist" << endl );

  vecbasevector b;
  b.ReadAll(FILE);

  cout.setf(ios::internal);

  size_t longest_read = b[0].size();
  for (size_t i = 0; i < b.size(); i++)
  {
      if (b[i].size() > longest_read) { longest_read = b[i].size(); }
  }

  vec<int> GC_counts(longest_read);
  vec<int> GC_totals(longest_read);

  ulonglong totalLength = 0, totalGC = 0;
  int start, end;
  unsigned int gc, size;
  for ( size_t i = 0; i < b.size( ); i++ ) {
    if (SHOW_INDEX)
      cout << setw(8) << i << "   ";
    size = b[i].size();
    start = min(TRIM, size);
    end = max(start, static_cast<int>(size) - static_cast<int>(TRIM));
    size = end-start;
    gc = GcBases(b[i], start, end);
    cout << PERCENT_RATIO(3, gc, size) << "\n";
    totalLength += size;
    totalGC += gc;

    if (COMPUTE_AVERAGE_PER_POSITION)
    {
        bvec const& bv = b[i];
        for (int j = 0; j < bv.isize(); j++)
        {
            if  ( IsGC(bv[j]) )
                GC_counts[j] += 1;
            GC_totals[j] += 1;
        }
    }
  }

  cout << "Number of items = " << b.size() << endl
       << "Total Bases = " << totalLength << endl
       << "Average GC content = "
       << PERCENT_RATIO(3, totalGC, totalLength) << endl;

    if (COMPUTE_AVERAGE_PER_POSITION)
    {
        cout << endl;
        for (size_t i = 0; i < longest_read; i++)
        {
            cout << i << " " << GC_totals[i] << " " << GC_counts[i] << " " << PERCENT_RATIO(3, GC_counts[i], GC_totals[i]) << endl;
        }
        cout.flush();
    }

  return 0;
}
