/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// FindCutters.  Find N-cutters having a K-base overhang (where K positive is 3',
// K zero is blunt, K negative is 5'), that are available from New England Biolabs.
// Ambiguous bases are not allowed unless AMB=True is specified.
//
// This assumes that the restriction enzyme list
// http://rebase.neb.com/rebase/link_withrefm
// has been downloaded into /wga/scr4/MolBio/Enzymes.

#include "MainTools.h"
#include "FastIfstream.h"
#include "dna/Bases.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int(N);
     CommandArgument_Int_OrDefault(K, 0);
     CommandArgument_Bool_OrDefault(AMB, False);
     EndCommandArguments;

     if ( (N-K) % 2 != 0 )
     {    cout << "N - K must be even.\n";
          exit(1);    }

     fast_ifstream in ( "/wga/scr4/MolBio/Enzymes/link_withrefm" );
     String line, enzyme, site, suppliers;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( !line.Contains( "<1>", 0 ) ) continue;
          enzyme = line.After( "<1>" );
          getline( in, line );
          getline( in, line );
          ForceAssert( line.Contains( "<3>", 0 ) );
          site = line.After( "<3>" );
          if ( site.isize( ) != N + 1 ) continue;
          if ( site[ (N+K)/2 ] != '^' ) continue;
          Bool ambiguous = False;
          for ( int i = 0; i <= N; i++ )
               if ( i != (N+K)/2 && !Base::isCanonicalBase(site[i]) ) ambiguous = True;
          if ( !AMB && ambiguous ) continue;
          getline( in, line );
          getline( in, line );
          getline( in, line );
          getline( in, line );
          ForceAssert( line.Contains( "<7>", 0 ) );
          suppliers = line.After( "<7>" );
          if ( !suppliers.Contains( "N" ) ) continue;
          cout << enzyme << " " << site << "\n";    }    }
