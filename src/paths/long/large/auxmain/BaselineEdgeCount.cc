///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BaselineEdgeCount.  Find the averge number of edges in a 50 kb assembly of a 
// human region from the NA12878 assisted regional assemblies that we ran.
// Failed assemblies are ignored.
// Answer: median = 79, mean = 77.36.
// Inferred baseline goal: 1.55 edges/kb.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "FastIfstream.h"
#include "MainTools.h"

int main( )
{    RunTime( );
     String dir = "/wga/scr4/human_assemblies/1/v9";
     vec<String> chrs;
     for ( int c = 1; c <= 22; c++ )
          chrs.push_back( ToString(c) );
     chrs.push_back( "X" );
     vec<int> NEDGES;
     #pragma omp parallel for
     for ( int i = 0; i < chrs.isize( ); i++ )
     {    String cdir = dir + "/" + chrs[i];
          vec<String> reg = AllFiles(cdir);
          vec<int> ne;
          for ( int j = 0; j < reg.isize( ); j++ )
          {    String rdir = cdir + "/" + reg[j];
               if ( !IsDirectory(rdir) ) continue;
               String fasta = rdir + "/out.final.fasta";
               if ( !IsRegularFile(fasta) ) continue;
               fast_ifstream in(fasta);
               int nedges = 0;
               String line;
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    if ( line.Contains( ">", 0 ) ) nedges++;    }
               ne.push_back(nedges);    }
          #pragma omp critical
          {    NEDGES.append(ne);    }    }
     Sort(NEDGES);
     PRINT2( Median(NEDGES), Mean(NEDGES) );    }
