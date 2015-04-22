// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// Fasta: convert a file reads.fastb into a human-readable reads.fasta.
// If run with CLEAN=True, generate a blastable file (no blank lines,
// no empty reads).  

// If NAMES=True, use original read names.

// If MAXREADS > 0, print no more than MAXREADS reads.

// If SUBSET is not empty, it is parsed as an IntSet (c.f. ParseSet.h) and only 
// those entries are printed.

// If FOLD_SIZE specified, fold fasta sequences into chunks of that size.  This
// is needed for some external programs.

#include "Basevector.h"
#include "MainTools.h"
#include "math/Functions.h"
#include "ParseSet.h"
#include "STLExtensions.h"

int main( int argc, char *argv[] )
{
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBSET, "");
     CommandArgument_String_OrDefault(HEAD, "reads");
     CommandArgument_Bool_OrDefault(CLEAN, False);
     CommandArgument_Bool_OrDefault(NAMES, False);
     CommandArgument_UnsignedInt_OrDefault(MAXREADS, 0);
     CommandArgument_Bool_OrDefault(GZIPPED, False);
     CommandArgument_UnsignedInt_OrDefault(FOLD_SIZE, 0);
     CommandArgument_UnsignedInt_OrDefault(MAX_BASES, 0);
     CommandArgument_String_OrDefault(SUFFIX, "");
     EndCommandArguments;

     String run_dir = PRE + "/" + DATA + "/" + RUN;

     String infile = run_dir + "/" + HEAD + ".fastb";
     String outfile = run_dir + "/" + HEAD + ".fasta";
     
     ForceAssert( IsRegularFile( infile ) );

     // Prep output filestream.
     ostream *pOut = 0;
     procbuf *pBuf = 0;
     if ( GZIPPED )
     { 
       outfile += ".gz";
       pBuf = new procbuf( outfile.c_str(), ios::out );
       pOut = new ostream( pBuf );
     }
     else
       pOut = new ofstream( outfile.c_str() );
     ostream &out = *pOut;

     int num_seqs = MastervecFileObjectCount( infile );

     // Gather ids.
     // TODO: potentially dangerous truncation of index by ids var
     vec<int> ids;
     if ( ! SUBSET.empty() )
       ParseIntSet( SUBSET, ids );
     else
     {
       ids.resize( num_seqs );
       iota( ids.begin(), ids.end(), 0 );
     }

     // Truncate ids if needed.
     if ( MAXREADS > 0 && 
          MAXREADS < ids.size() )
       ids.resize( MAXREADS );

     // Load data.
     vecbasevector bases;
     bases.SparseRead( infile, ids, 0 );

     vecString seqnames;
     if (NAMES) {
       String namesfile = run_dir + "/" + HEAD + ".ids";
       seqnames.SparseRead( namesfile, ids, 0 );
     }

     // Print sequences.
     for ( size_t id = 0; id < bases.size( ); id++ )
     { 
       if ( CLEAN && bases[id].size( ) == 0 ) continue;
       if ( MAX_BASES > 0 && bases[id].size( ) > MAX_BASES ) continue;
       if ( SUBSET != "" && !BinMember( ids, id ) ) continue;
       if ( (int) FOLD_SIZE > 0 && bases[id].size( ) > FOLD_SIZE )
       {    int part = 1;
            for ( unsigned int i = 0; i < bases[id].size( ); i += FOLD_SIZE )
            {    static basevector b;
                 b.SetToSubOf( bases[id], i, Min( FOLD_SIZE, bases[id].size( ) - i ) );
                 static String bn;
                 if (NAMES) bn = seqnames[id];
                 else bn = "sequence_" + ToString(id);
                 bn += ".part" + ToString(part++);
                 b.Print( out, bn );    }    }
       else
       {    if (NAMES) bases[id].Print( out, seqnames[id] );
            else bases[id].Print( out, ToString(id) + SUFFIX );    }
       if ( ! CLEAN ) out << "\n";
     }

     // Deleting these pointers will flush the ostream and close the procbuf.
     delete pOut;
     delete pBuf;
}
