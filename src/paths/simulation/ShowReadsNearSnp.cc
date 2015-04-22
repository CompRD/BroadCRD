/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// ShowReadsNearSnp.  Given a simulated and error-corrected diploid data set, show 
// the reads that overlap a given SNP position.  This position is specified by the
// POS argument, which should be of the form a.b, where a is a genome contig id
// (< genome.fastb.size/2), and b is a base position on it.

#include "Basevector.h"
#include "math/Functions.h"
#include "MainTools.h"
#include "ReadLocation.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(POS);
     EndCommandArguments;

     int tig1 = POS.Before( "." ).Int( ), pos = POS.After( "." ).Int( );

     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;

     vecbasevector genome( data_dir + "/genome.fastb" );
     ForceAssertLt( static_cast<size_t>(tig1), genome.size( )/2 );
     ForceAssertGe( pos, 0 );
     ForceAssertLt( pos, genome[tig1].isize( ) );

     int tig2 = tig1 + genome.size( )/2;
     const basevector &g1 = genome[tig1], &g2 = genome[tig2];

     if ( g1[pos] == g2[pos] )
     {    cout << "That's not a SNP position.\n";
          exit(1);    }

     vecbasevector reads( run_dir + "/reads.fastb" );
     vecbasevector reads_orig( run_dir + "/reads.fastb.orig" );
     READ( run_dir + "/reads.before_removing_suspicious.ref.locs", 
          vec<read_location>, locs );
     BREAD2( run_dir + "/reads.id_map", vec<int>, id_map );
     vec<int> del_lowfreq;
     BinaryReader::readFile(run_dir + "/reads.deleted_by_lowfreq", &del_lowfreq );

     int maxr = 0;
     for ( size_t i = 0; i < reads_orig.size( ); i++ )
          maxr = Max( reads_orig[i].isize( ), maxr );

     cout << "Base codes:\n";
     cout << "1 means that the base is at a SNP position, and is the first allele,"
          << "\n  and was correct and stayed correct\n";
     cout << "2 means that the base is at a SNP position, and is the second allele,"
          << "\n  and was correct and stayed correct\n";
     cout << "Otherwise:\n";
     cout << ". means that the base on the read was correct and stayed correct\n";
     cout << "C means that the base on the read was corrected\n";
     cout << "E means that the base on the read was wrong and uncorrected\n";
     cout << "M means that the base was correct, then miscorrected\n";
     cout << "X means that the base was incorrect, then miscorrected\n\n";

     int maxidsize = ToString( reads.size( ) - 1 ).size( );
     for ( int i = 0; i < locs.isize( ); i++ )
     {    int start = locs[i].Start( ), stop = locs[i].Stop( ) + 1;
          if ( locs[i].Contig( ) != tig1 && locs[i].Contig( ) != tig2 ) continue;
          if ( pos < start || pos >= stop ) continue;
          const basevector& g = ( locs[i].Contig( ) == tig1 ? g1 : g2 );
          int id = locs[i].ReadId( );
          String idstring = ToString(id);
          for ( int j = 0; j < maxidsize - idstring.isize( ); j++ )
               cout << " ";
          cout << idstring << " ";
          basevector ro = reads_orig[id];
          for ( int j = 0; j < maxr + start - pos; j++ )
               cout << " ";
          if ( locs[i].Contig( ) == tig1 ) cout << "[1] ";
          else cout << "[2] ";
          if ( locs[i].Rc( ) ) ro.ReverseComplement( );
          if ( id_map[id] < 0 )
          {    for ( int j = 0; j < ro.isize( ); j++ )
               {    int p = j + start;
                    if ( ro[j] == g[p] ) 
                    {    if ( g1[p] == g2[p] ) cout << ".";
                         else if ( locs[i].Contig( ) == tig1 ) cout << "1";
                         else cout << "2";    }
                    else cout << "E";    }
               cout << " [deleted";
               if ( BinMember( del_lowfreq, id ) ) cout << " - low freq kmer";
               cout << "]";    }
          else
          {    basevector r = reads[ id_map[id] ];
               if ( locs[i].Rc( ) ) r.ReverseComplement( ); 
               for ( int j = 0; j < r.isize( ); j++ )
               {    int p = j + start;
                    if ( r[j] == g[p] && ro[j] == g[p] ) 
                    {    if ( g1[p] == g2[p] ) cout << ".";
                         else if ( locs[i].Contig( ) == tig1 ) cout << "1";
                         else cout << "2";    }
                    if ( r[j] == g[p] && ro[j] != g[p] ) cout << "C";
                    if ( r[j] != g[p] && ro[j] == g[p] ) cout << "M";
                    if ( r[j] != g[p] && ro[j] != g[p] && r[j] == ro[j] ) 
                         cout << "E";
                    if ( r[j] != g[p] && ro[j] != g[p] && r[j] != ro[j] ) 
                         cout << "X";    }    }
          cout << "\n";    }    }
