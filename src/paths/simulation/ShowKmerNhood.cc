// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

// ShowKmerNhood: pick count k-mers at random from the collection of sequences G.  
// For each, show all the occurrences of the k-mer, together with w bases on each 
// side.
//
// This is intended as a exploratory tool to guide the development of 
// error-correction algorithms.
//
// Both G.fastb and G.lookup should exist.

#include "Basevector.h"
#include "MainTools.h"
#include "lookup/LookupTable.h"
#include "random/Random.h"
     
#define ABORT(MSG)                                  \
{    cout << MSG << "  Abort." << endl << endl;     \
     exit(1);    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(G);
     CommandArgument_Int(w);
     CommandArgument_Int(K);
     CommandArgument_Int(count);
     EndCommandArguments;

     // Select random K-mers from the genome.

     vec<basevector> kmers(count);
     vecbasevector genome( G + ".fastb" );
     longlong gbases = 0;
     for ( size_t i = 0; i < genome.size( ); i++ )
     {    if ( genome[i].isize( ) < K ) continue;
          gbases += genome[i].isize( ) - K + 1;    }
     for ( int i = 0; i < count; i++ )
     {    longlong pos = big_random( ) % gbases;
          for ( size_t j = 0; j < genome.size( ); j++ )
          {    const basevector& g = genome[j];
               if ( g.isize( ) < K ) continue;
               if ( pos > g.isize( ) - K )
               {    pos -= ( g.isize( ) - K + 1 );
                    continue;    }
               kmers[i].SetToSubOf( g, pos, K );
               break;    }    }

     // Get header info from lookup table.

     if ( !IsRegularFile( G + ".lookup" ) ) 
          ABORT( "I can't find the file " << G << ".lookup." );
     lookup_table look(G + ".lookup");
     int Klook = look.K( );
     ForceAssertLe( Klook, K );

     // Find hits.

     vec< vec<basevector> > hits(count);
     for ( unsigned int i = 0; i < look.NChunks( ); i++ )
     {    look.ReadChunk(i);
          for ( int j = 0; j < kmers.isize( ); j++ )
          {    basevector s = kmers[j];
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 ) s.ReverseComplement( );
                    unsigned int index = Index( s, 0, Klook );
                    unsigned int start = look.StartLocs(index);
                    unsigned int stop = look.StopLocs(index);
                    for ( unsigned int l = start; l < stop; l++ )
                    {    unsigned int c, cpos;
                         look.GetContigPos( look.Locs(l), c, cpos );
                         if ( (int) cpos + K + w > genome[c].isize( ) ) continue;
                         if ( (int) cpos < w ) continue;
                         Bool match = True;
                         for ( int u = Klook; u < K; u++ )
                         {    if ( s[u] != genome[c][cpos+u] )
                              {    match = False;
                                   break;    }    }
                         if ( !match ) continue;
                         static basevector b;
                         b.SetToSubOf( genome[c], cpos - w, K + 2 * w );
                         if ( pass == 2 ) b.ReverseComplement( );
                         hits[j].push_back(b);    }    }    }    }

     // Print hits.

     for ( int i = 0; i < count; i++ )
     {    cout << "\n" << i+1 << " = " << kmers[i].ToString( ) << ":\n\n";
          Sort( hits[i] );
          for ( int j = 0; j < hits[i].isize( ); j++ )
          {    int j2;
               for ( j2 = j + 1; j2 < hits[i].isize( ); j2++ )
                    if ( hits[i][j] != hits[i][j2] ) break;
               for ( int u = 0; u < w; u++ )
                    cout << as_base( hits[i][j][u] );
               cout << "-";
               for ( int u = w + K; u < 2 * w + K; u++ )
                    cout << as_base( hits[i][j][u] );
               cout << " [" << j2-j << "]\n";
               j = j2 - 1;    }    }    }
