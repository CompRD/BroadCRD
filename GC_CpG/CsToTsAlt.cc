/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// CsToTsAlt.  Given a fasta file HEAD.fasta, generate a fasta file
// HEAD.fw.CtoT.fasta by changing the Cs to Ts in HEAD.fasta, and a
// fasta file HEAD.rc.CtoT.fasta, by changing the Cs to Ts in the reverse
// complement of HEAD.fasta.

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"

int main( int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(HEAD);
     EndCommandArguments;

     fast_ifstream in( HEAD + ".fasta" );
     vec<String> seq, seqnames;
     String line, x;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( ">", 0 ) )
          {    if ( seqnames.nonempty( ) ) seq.push_back(x);
               seqnames.push_back( line.After( ">" ) );
               x = "";    }
          else x += line;    }
     seq.push_back(x);
     ForceAssertEq( seq.size( ), seqnames.size( ) );
     ForceAssertGt( seq.size( ), 0u );
     Ofstream ( fw, HEAD + ".fw.CtoT.fasta" );
     Ofstream ( rc, HEAD + ".rc.CtoT.fasta" );
     for ( int i = 0; i < seq.isize( ); i++ )
     {    for ( int pass = 1; pass <= 2; pass++ )
          {    String S = seq[i];
               if ( pass == 2 ) StringReverseComplement( S, S );
               for ( int j = 0; j < S.isize( ); j++ )
               {    if ( S[j] == 'c' ) S[j] = 't';
                    if ( S[j] == 'C' ) S[j] = 'T';    }
               ( pass == 1 ? fw : rc ) << ">" << seqnames[i] 
                    << ( pass == 2 ? " [rc]" : "" ) << " [C-->T]" << "\n";
               for ( int i = 0; i < S.isize( ); i++ )
               {    if ( i > 0 && i % 80 == 0 ) ( pass == 1 ? fw : rc ) << "\n";
                    ( pass == 1 ? fw : rc ) << S[i];    }
               ( pass == 1 ? fw : rc ) << "\n";    }    }    }
