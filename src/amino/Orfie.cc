///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Find open reading frames of minimum length, exploring the graph up to a fixed
// depth.  An open reading frame is defined as starting at a start codon ATG
// and stopping at a stop code TAA or TAG or TGA.
//
// Note that ORFs are reported on an edge, and also separately on its reverse
// complement.
//
// Very inefficiently implemented.
//
// Hardcoded for now.
//
// Will also find Stop Stop Frames (SSFs).

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "amino/Amino.h"
#include "paths/HyperBasevector.h"

int main( )
{    RunTime( );

     // Directory.

     String dir = "/wga/scr4/jaffe/GapToy/mouse.1456abcde/a.final";

     // Load data.

     HyperBasevector hb;
     BinaryReader::readFile( dir + "/a.hbv", &hb );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     int K = hb.K( );

     // Heuristics.

     const int max_path = 8;
     const int min_orf = 100;
     Bool stop_stop = True;

     // Find ORFs.

     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( hb.Bases(e) == 0 ) continue;
          vec<vec<int>> paths;
          paths.push_back( {e} );
          for ( int d = 0; d < max_path - 1; d++ )
          {    vec<vec<int>> paths2;
               for ( int j = 0; j < paths.isize( ); j++ )
               {    int v = to_right[ paths[j].back( ) ];
                    int exts = 0;
                    for ( int l = 0; l < hb.From(v).isize( ); l++ )
                    {    if ( hb.Bases( hb.IFrom(v,l) ) > 0 )
                         {    vec<int> p = paths[j];
                              p.push_back( hb.IFrom(v,l) );
                              paths2.push_back(p);    
                              exts++;    }    }
                    if ( exts == 0 ) paths2.push_back( paths[j] );    }
               paths = paths2;    }
          vec< triple<int,vec<int>,int> > orfs;
          for ( int u = 0; u < paths.isize( ); u++ )
          {    basevector E = hb.Cat( paths[u] );
               for ( int j = 0; j <= E.isize( ) - 3; j++ )
               {    
                    if ( !stop_stop )
                    {    if ( as_base( E[j] ) != 'A' || as_base( E[j+1] ) != 'T'
                              || as_base( E[j+2] ) != 'G' )
                         {    continue;    }    }

                    else
                    {    Bool stop = False;
                         if ( as_base( E[j] ) == 'T' )
                         {    if ( as_base( E[j+1] ) == 'A' 
                                   && as_base( E[j+2] ) == 'A' )
                              {    stop = True;    }
                              if ( as_base( E[j+1] ) == 'A' 
                                   && as_base( E[j+2] ) == 'G' )
                              {    stop = True;    }
                              if ( as_base( E[j+1] ) == 'G' 
                                   && as_base( E[j+2] ) == 'A' )
                              {    stop = True;    }    }
                         if (stop) continue;    }

                    int n = -1;
                    for ( int k = j + 3; k <= E.isize( ) - 3; k++ )
                    {    Bool stop = False;
                         if ( as_base( E[k] ) == 'T' )
                         {    if ( as_base( E[k+1] ) == 'A' 
                                   && as_base( E[k+2] ) == 'A' )
                              {    stop = True;    }
                              if ( as_base( E[k+1] ) == 'A' 
                                   && as_base( E[k+2] ) == 'G' )
                              {    stop = True;    }
                              if ( as_base( E[k+1] ) == 'G' 
                                   && as_base( E[k+2] ) == 'A' )
                              {    stop = True;    }    }
                         n++;
                         if (stop)
                         {    if ( n >= min_orf )
                              {    vec<int> p = paths[u];
                                   int l;
                                   for ( l = 1; l < p.isize( ); l++ )
                                   {    int len = hb.Bases( p[0] );
                                        for ( int m = 1; m < l; m++ )
                                             len += hb.Kmers( p[m] );
                                        if ( j + 3*n <= len ) break;    }
                                   p.resize(l);
                                   orfs.push( j, p, n );    }
                              break;    }    }    }    }    
          UniqueSort(orfs);
          if ( orfs.nonempty( ) )
          {    
               #pragma omp critical
               {    for ( int i = 0; i < orfs.isize( ); i++ )
                    {    int j = orfs[i].first, n = orfs[i].third;
                         const vec<int>& p = orfs[i].second;
                         if ( j >= hb.Kmers( p[0] )
                              && hb.From( to_right[ p[0] ] ).nonempty( ) )
                         {    continue;    }
                         basevector E = hb.Cat(p);
                         String aa;
                         for ( int l = 0; l < n; l++ )
                         {    int s = j + 3 * (l+1);
                              uint codon = 4*E[s] + 2*E[s+1] + E[s+2];
                              aa.push_back(
                                   AminoAcid::forCodon(codon).getSymbol( ) );    }
                         cout << "ORF of length " << n << " on " << j << " : "
                              << printSeq(p) << ", peptide = " << aa
                              << endl;    }    }    }    }    }
