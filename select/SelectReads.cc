// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


// SelectReads: create reads for a new "genome" by taking a subset of the
// reads for a given genome.  If EXPAND >= 1, also include all reads which
// overlap the given reads, and force closure under pairing.  Iterate the expansion 
// process, EXPAND times.

// Note: the code now reads gzipped input files.  It can easily be switched
// back, or (better), allow for gzipped or not.

#include <algorithm>

#include "Alignment.h"
#include "system/Assert.h"
#include "Basevector.h"
#include "FastIfstream.h"
#include "system/ParsedArgs.h"
#include "Quality.h"
#include "ReadPairing.h"
#include "system/RunTime.h"
#include "Set.h"
#include "String.h"
#include "system/System.h"
#include "system/Types.h"
#include "Vec.h"

int main( int argc, char *argv[] )
{
     RunTime();

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(FILE_OF_IDS);
     CommandArgument_String(NEW_DATA);
     CommandArgument_UnsignedInt_OrDefault(EXPAND, 0);
     CommandArgument_Double_OrDefault(MAX_SCORE_FOR_EXPAND, 0);
     CommandArgument_UnsignedInt_OrDefault(MIN_OVERLAP_FOR_EXPAND, 0);
     CommandArgument_Bool_OrDefault(PAIR_EXPAND, True);
     EndCommandArguments;

     // Read in the read names.

     String run_dir = PRE + "/" + DATA + "/" + RUN;
     vecString ids( run_dir + "/orig/reads.ids" );
     int N = ids.size( );

     // Expand the list of ids.

     vec<int> source;
     set<int> source2, source3;
     Ifstream( nids, FILE_OF_IDS );
     while(1)
     {    int id;
          nids >> id;
          if ( !nids ) break;
          source.push_back(id);
          source2.insert(id);    }

     if ( EXPAND > 0 )
     {    
          String run_dir = PRE + "/" + DATA + "/" + RUN;
          READ( run_dir + "/aligns.total2", vec<alignment_plus>, all_aligns );
          vec<int> all_aligns_index(N, -1);
          for ( int i = (int) all_aligns.size( ) - 1; i >= 0; i-- )
               all_aligns_index[ all_aligns[i].Id1( ) ] = i;

          Ifstream( pairs_f, run_dir + "/reads.pairto" );
          int pairs_length;
          pairs_f >> pairs_length;
          vec<read_pairing> pairs(pairs_length);
          for ( int i = 0; i < pairs_length; i++ ) 
               pairs_f >> pairs[i];
          vec<int> pairs_index(N, -1);
          for ( unsigned int i = 0; i < pairs.size( ); i++ )
               if ( pairs[i].Alive( ) )
                    pairs_index[ pairs[i].id1 ] = pairs_index[ pairs[i].id2 ] = i;
          
          for ( unsigned int j = 0; j < EXPAND; j++ )
          {    cout << "starting expansion\n";
               for ( unsigned int i = 0; i < source.size( ); i++ )
               {    int id = source[i];
                    int ind = all_aligns_index[id];
                    for (unsigned int l = ind; ind >= 0 && l < all_aligns.size( ); 
                         l++)
                    {    const alignment_plus& ap = all_aligns[l];
                         if ( ap.Id1( ) > id ) break;
                         if ( MAX_SCORE_FOR_EXPAND > 0 
                              && ap.score > MAX_SCORE_FOR_EXPAND ) continue;
                         if ( ap.a.Pos1( ) - ap.a.pos1( ) 
                              < (int) MIN_OVERLAP_FOR_EXPAND )
                              continue;
                         if ( source2.insert( ap.Id2( ) ).second )
                              cout << "adding read " << ap.Id2( ) << "\n";    }    }
               for ( set<int>::iterator i = source2.begin( ); i != source2.end( ); 
                    ++i )
               {    source3.insert(*i);
                    int p = pairs_index[*i];
                    if ( PAIR_EXPAND && p >= 0 )
                    {    source3.insert( pairs[p].id1 );
                         source3.insert( pairs[p].id2 );    }    }
               source.clear( );
               for ( set<int>::iterator i = source3.begin( ); i != source3.end( ); 
                    ++i )
                    source.push_back(*i);

               cout << "iteration " << j+1 
                    << ",  initial list of reads expanded to " 
                    << source.size( ) << endl;    }    }
     
     // Read in the ids from FILE_OF_IDS.  Determine the corresponding read names.

     vec<String> read_names;
     for ( unsigned int i = 0; i < source.size( ); i++ )
          read_names.push_back( ids[ source[i] ] );
     sort( read_names.begin( ), read_names.end( ) );

     // Create the new files of reads and quality scores.

     vec<Bool> found( read_names.size( ), False );
     Mkdir777( PRE + "/" + NEW_DATA );
     System( "echo " + DATA + " > " + PRE + "/" + NEW_DATA + "/source.data" );
     Ofstream( bases_out, PRE + "/" + NEW_DATA + "/reads.fasta" );
     Ofstream( quals_out, PRE + "/" + NEW_DATA + "/reads.qual" );

     // fast_ifstream bases_in( PRE + "/" + DATA + "/reads.fasta" );
     // fast_ifstream quals_in( PRE + "/" + DATA + "/reads.qual" );

     fast_pipe_ifstream bases_in( "zcat " + PRE + "/" + DATA + "/reads.fasta.gz" );
     fast_pipe_ifstream quals_in( "zcat " + PRE + "/" + DATA + "/reads.qual.gz" );

     String line, s;
     getline( bases_in, line );
     Assert( line.Contains( ">", 0 ) );
     while(1)
     {    if ( bases_in.fail( ) ) break;
          s = line.After( ">" );
          if ( s.Contains( ".exp" ) ) s = s.Before( ".exp" );
          if ( s.Contains( ".scf" ) ) s = s.Before( ".scf" );
          int p = BinPosition( read_names, String(s) );
          if ( p >= 0 ) 
          {    found[p] = True;
               bases_out << line << "\n";    }
          while(1)
          {    getline( bases_in, line );
               if ( bases_in.fail( ) ) break;
               if ( line.Contains( ">", 0 ) ) break;
               if ( p >= 0 ) bases_out << line << "\n";    }    }
     Ofstream( ids_out, PRE + "/" + NEW_DATA + "/source.ids" );
     for ( unsigned int i = 0; i < found.size( ); i++ )
     {    if ( !found[i] ) cout << "read " << i << " not found\n";
          else ids_out << i << "\n";    }
     getline( quals_in, line );
     Assert( line.Contains( ">", 0 ) );
     while(1)
     {    if ( quals_in.fail( ) ) break;
          s = line.After( ">" );
          if ( s.Contains( ".exp" ) ) s = s.Before( ".exp" );
          if ( s.Contains( ".scf" ) ) s = s.Before( ".scf" );
          int p = BinPosition( read_names, String(s) );
          if ( p >= 0 ) 
          {    found[p] = True;
               quals_out << line << "\n";    }
          while(1)
          {    getline( quals_in, line );
               if ( quals_in.fail( ) ) break;
               if ( line.Contains( ">", 0 ) ) break;
               if ( p >= 0 ) quals_out << line << "\n";    }    }

     // Copy other source files.

     if ( IsRegularFile( PRE + "/" + DATA + "/contigs.fasta" ) )
          Cp( PRE + "/" + DATA + "/contigs.fasta", PRE + "/" + NEW_DATA );
     if ( IsRegularFile( PRE + "/" + DATA + "/rel" ) )
          Cp( PRE + "/" + DATA + "/rel", PRE + "/" + NEW_DATA );
     if ( IsRegularFile( PRE + "/" + DATA + "/bad_forwards" ) )
          Cp( PRE + "/" + DATA + "/bad_forwards", PRE + "/" + NEW_DATA );
     if ( IsRegularFile( PRE + "/" + DATA + "/bad_reverses" ) )
          Cp( PRE + "/" + DATA + "/bad_reverses", PRE + "/" + NEW_DATA );    }
