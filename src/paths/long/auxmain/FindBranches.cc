///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FindBranches.  Given a region of the human genome, see if it has 500-600
// base perfect repeats somewhere else in the genome.  Hardcoded and writes
// temporary files.  Bad.

// MakeDepend: dependency QueryLookupTable

#include "Basevector.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"

int main(int argc, char *argv[])
{
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_String_Doc(IN, "fastb file for chunk of genome");
     EndCommandArguments;

     const int K = 500;
     const int step = 100;
     int count = 0;
     vecbasevector genome(IN);
     {    Ofstream( out, "woof.fasta" );
          for ( int g = 0; g < (int) genome.size( ); g++ )
          for ( int i = 0; i <= genome[g].isize( ) - K; i += step )
          {    basevector b( genome[g], i, K );
               b.Print( out, ToString(i) + "-" + ToString(i+K) );
               count++;    }    }
     SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS=woof.fasta "
          "L=/wga/scr4/bigrefs/human19/genome.lookup.lookup PARSEABLE=True "
          " MF=5000 > woof.aligns" );
     vec<look_align> aligns;
     vec< vec<int> > aligns_index;
     LoadLookAligns( "woof.aligns", aligns, aligns_index, count );
     for ( int i = 0; i < count; i++ )
     {    int n = aligns_index[i].size( );
          if ( n != 1 ) PRINT2( i, n );    }    }
