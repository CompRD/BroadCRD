///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BreakAtNs: read in a fasta file.  For each instance in a record where at least
// MIN_TO_BREAK Ns or ns occur in a row, break the record.

#include "MainTools.h"
#include "FastIfstream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_String(OUT);
     CommandArgument_UnsignedInt(MIN_TO_BREAK);
     CommandArgument_String_OrDefault(N_HISTOGRAM, "");
     EndCommandArguments;

     // Generate histogram of gap sizes (Optional)
     const size_t max_gap = 65536;
     size_t last_bin = 0;
     vec<size_t> counts;
     if (N_HISTOGRAM != "")
	 counts.resize(max_gap, 0);

     fast_ifstream in(IN);
     Ofstream( out, OUT );
     vec<char> record;
     String line, header;
     while(1)
     {    getline( in, line );
          if ( !in.fail( ) && !line.Contains( ">", 0 ) )
          {    for ( int i = 0; i < (int) line.size( ); i++ )
                    record.push_back( line[i] );    }
          else
          {    int part = 1, last = 0;
               for ( int i = 0; i < record.isize( ); i++ )
               {    if ( record[i] == 'N' || record[i] == 'n' )
                    {    int j;
                         for ( j = i + 1; j < record.isize( ); j++ )
                              if ( record[j] != 'N' && record[j] != 'n' ) break;
			 if (N_HISTOGRAM != "") {
			     size_t bin_id = Min(max_gap - 1, (size_t) j - i );
			     counts[bin_id]++;
			     last_bin = Max(last_bin, bin_id);
			 }
                         if ( j - i >= (int) MIN_TO_BREAK )
                         {    if ( i > last )
                              {    if ( header.Contains( " " ) )
                                   {    out << header.Before( " " ) << ".part"
                                             << part << " " << header.After( " " )
                                             << "\n";     }
                                   else out << header << ".part" << part << "\n";
                                   for ( int k = last; k < i; k++ )
                                   {    if ( k - last > 0 && (k - last) % 80 == 0 ) 
                                             out << "\n";
                                        out << record[k];    }
                                   out << "\n";
                                   ++part;    }
                              last = j;    }
                         i = j - 1;    }    }
               if ( last < record.isize( ) )
               {    if ( header.Contains( " " ) )
                    {    out << header.Before( " " ) << ".part"
                              << part << " " << header.After( " " )
                              << "\n";     }
                    else out << header << ".part" << part << "\n";
                    for ( int k = last; k < record.isize( ); k++ )
                    {    if ( k - last > 0 && (k - last) % 80 == 0 ) out << "\n";
                         out << record[k];    }
                    out << "\n";    }
               if ( in.fail( ) ) break;
               header = line;
               record.clear( );    }    }    

     // Write histogram of gap sizes (Optional)    
     if (N_HISTOGRAM != "") {
	 Ofstream( histo_out, N_HISTOGRAM );
	 for (size_t i = 0; i <= last_bin; i++) {
	     histo_out << i << " " << counts[i] << endl;
	 }
     }
}
