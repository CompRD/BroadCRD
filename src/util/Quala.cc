/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Quala: convert a file reads.qualb into a human-readable reads.qual.

// If run with CLEAN=True, generate a blastable file (no blank lines,
// no empty reads).  If NAMES=True, use original read names.

// If MAXREADS > 0, print no more than MAXREADS reads.

// If SUBSET is not empty, it is parsed as an IntSet (c.f. ParseSet.h)
// and only those entries are printed.

#include "MainTools.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "math/Functions.h"

int main( int argc, char *argv[] )
{
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(HEAD, "reads");
     CommandArgument_Bool_OrDefault(CLEAN, False);
     CommandArgument_Bool_OrDefault(NAMES, False);
     CommandArgument_String_OrDefault(SUBSET, "");
     CommandArgument_LongLong_OrDefault(LL_MAXREADS, 0);
     EndCommandArguments;

     vecqvec::size_type MAXREADS = LL_MAXREADS;
     String run_dir = PRE + "/" + DATA + "/" + RUN;

     String qualb_file = run_dir + "/" + HEAD + ".qualb";
     String quala_file = run_dir + "/" + HEAD + ".qual";
     String ids_file = run_dir + "/" + HEAD + ".ids";

     int n_reads = MastervecFileObjectCount( qualb_file );

     vec<int> select_ids;
     if ( SUBSET == "" ) {
       select_ids.resize( n_reads );
       iota( select_ids.begin( ), select_ids.end( ), 0 );
     }
     else
       ParseIntSet( SUBSET, select_ids );

     vec<Bool> select;
     select.resize( n_reads, False );
     for (int ii=0; ii<select_ids.isize( ); ii++)
       select[ select_ids[ii] ] = True;

     vecString ids;
     if (NAMES) {
       if ( IsAsciiVec( ids_file ) ) {
	 READ( ids_file, vec<String>, oldseqnames );
	 ids.assign(oldseqnames.begin(),oldseqnames.end());
       } else {
	 ids.ReadAll( ids_file );
       }
     }

     vecqualvector Q;
     Q.ReadAll( qualb_file );

     Ofstream( out, quala_file );
     if ( MAXREADS == 0 ) MAXREADS = Q.size( );
     MAXREADS = Min(MAXREADS,Q.size());
     for ( unsigned long id = 0; id < MAXREADS; id++ ) {
       if ( ! select[id] ) continue;
       if ( CLEAN && Q[id].size( ) == 0 ) continue;
       if ( !NAMES ) out << ">sequence_" << id << "\n";
       else out << ">" << ids[id] << "\n";
       for ( qvec::size_type j = 0; j < Q[id].size( ); ++j ) {
	 if ( j % 25 == 0 && j > 0 ) out << "\n";
	 else if ( j > 0 ) out << " ";
	 out << (int) Q[id][j];
       }
       out << ( CLEAN ? "\n" : "\n\n" );
     }

}
