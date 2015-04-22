///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Tool designed for Assemblathon competition to clean up ambigiuity clumps "
  "introduced by FlattenHKP.";

#include "Fastavector.h"
#include "MainTools.h"

bool is_ambiguous(char base) {
  return !(base == 'A' || base == 'T' || base == 'C' || base == 'G') ;
}

size_t tidy_up(fastavector &contig, size_t window, double min_N) {

  size_t events = 0;
  size_t last_clump = 0;
  for ( size_t j = 0; j < contig.size( ) - window; j++ ) {
    size_t N = 0, total = 0;
    if ( !is_ambiguous(contig[j]) ) continue;
    size_t k;
    for ( k = j; k < j + window; k++ ) {
      total++;
      if ( is_ambiguous(contig[k]) ) N++;
    }
    if ( k == contig.size( ) )
	k--;
    while( !is_ambiguous(contig[k]) ) 
      k--;
    if (  double (N) / double (total) >= min_N) {
      if (j - last_clump != 1)
	events++;
      last_clump = j;
      for ( size_t u = j; u <= k; u++ )
	contig[u] = 'N';
    }
  }

  return events;
} 


int main( int argc, char *argv[] )
{
  RunTime( );
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(HEAD_IN,
    "Assembly to tidy");
  CommandArgument_String_Doc(HEAD_OUT,
    "Tidied assembly");
  CommandArgument_Int_OrDefault_Doc(WINDOW, 20,
    "Size of window in which to look for clumps");
  CommandArgument_Double_OrDefault_Doc(MIN_N, 0.6,
    "Minimum fraction of ambiguities in window to identify a clump ");

  EndCommandArguments;

  String fasta_in = HEAD_IN + ".fasta";
  String fasta_out = HEAD_OUT + ".fasta";

  cout << "Processing: " << fasta_in << endl;
  
  vec<String> names;
  vec<fastavector> tigsa;
  LoadFromFastaFile( fasta_in, tigsa, names );

  size_t initial_amb_count = 0;
  size_t final_amb_count = 0;
  size_t event_count = 0;

  Ofstream( out, fasta_out );
  for ( size_t i = 0; i < tigsa.size( ); i++ )     { 
    fastavector S = tigsa[i];

    initial_amb_count += S.AmbCount();
    event_count += tidy_up(S, WINDOW, MIN_N);
    final_amb_count += S.AmbCount();

    S.Print( out, names[i] );    
  }    

  cout << "  Ambiguity Clump Count: " << event_count << endl;
  cout << "Initial Ambiguity Count: " << initial_amb_count << endl;
  cout << "  Final Ambiguity Count: " << final_amb_count << endl;

  cout << "Wrote: " << fasta_out << endl;

}


