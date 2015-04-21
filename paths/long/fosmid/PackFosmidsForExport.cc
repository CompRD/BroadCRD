///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// PackFosmidsForExport.  Create a fasta file for the 'finished' NA12878 Fosmids.
// Writes to standard output.

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"

int main( )
{    const String version = "1.0";
     String fdir = "/wga/dev/references/Homo_sapiens/NA12878_Fosmid_Pool.regions.fin";
     vec<String> all = AllFiles(fdir);
     for ( int i = 0; i < all.isize( ); i++ )
     {    if ( !all[i].Contains( ".fasta", -1 ) ) continue;
          int id = all[i].Between( ".", "." ).Int( );
          vecbasevector b;
          FetchReads( b, 0, fdir + "/" + all[i] );
          ForceAssertEq( (int) b.size( ), 1 );
          b[0].Print( cout, "NA12878_Fosmid_" + ToString(id) 
               + "_version_" + version );    }    }
