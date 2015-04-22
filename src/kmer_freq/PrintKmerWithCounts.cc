/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "feudal/BinaryStream.h"
#include "kmers/KmerRecord.h"

template <int K>
void PrintKWC( const String &kmers_file, ostream &out );

/**
 * PrintKmerWithCounts
 *
 * Load the binary vec< kmer_with_count<K> > and print it in human
 * readable format.
 *
 * KSIZE: kmer size
 * KMERS: kmers, loaded with BinaryReader
 * ARCHIVE: save to a file parallel to KMERS (use cout if False)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( KMERS );
  CommandArgument_Bool_OrDefault( ARCHIVE, True );
  EndCommandArguments;
  
  // Out stream.
  ofstream archive_out;
  if ( ARCHIVE ) {
    String archive_file =  KMERS + ".human";
    archive_out.open( archive_file.c_str( ) );
    PrintCommandPretty( archive_out );
  }
  ostream &out = ARCHIVE ? archive_out : * (ostream *) &cout;

  // Call core function.
#define CORE( K_ ) PrintKWC<K_> ( KMERS, out );
  DISPATCH_ON_K( K, CORE );
  
}

/**
 * PrintKWC
 *
 * kmers_file: full path name of kmers
 * out: output stream
 */
template <int K>
void PrintKWC( const String &kmers_file, ostream &out )
{
  // Load and print.
  out << Date( ) << ": loading" << endl;
  vec< kmer_with_count<K> > kmers;
  BinaryReader::readFile( kmers_file, &kmers );
  out << "\n";

  bvec bases;
  for (int ii=0; ii<(int)kmers.size( ); ii++) {
    kmers[ii].GetBasevector( bases );
    out << ii << "/" << kmers.size( ) - 1 << "   "
	 << bases.ToString( ) << "   "
	 << kmers[ii].Count( ) << "\n";
  }
  out << "\n";

  // Done.
  out << Date( ) << ": done" << endl;
}

/**
 * Instantiate templates.
 */
#define INSTANTIATE_PRINT_KWC( K_, dummy )  \
template void PrintKWC<K_>( const String &, ostream & )

FOR_ALL_K( INSTANTIATE_PRINT_KWC, unused );

