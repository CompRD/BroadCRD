///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// DisplayPlus: display a neighborhood of a SupportedHyperBasevector, seeded on 
// given edges (to be rendered green) and extended to given depth in the graph.  
// Boundary vertices are labeled red.  Compare DisplayNhood.
// Looks for IN_HEAD.{hbv,efx}.

#include "Basevector.h"
#include "MainTools.h"
#include "paths/long/DisplayTools.h"

int main(int argc, char *argv[])
{
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_String_Doc(IN_HEAD, "looks for IN_HEAD.{hbv,efx}");
     CommandArgument_String_Doc(OUT_HEAD, "generate OUT_HEAD.{dot,fasta}");
     CommandArgument_String_Doc(SEEDS, "use these edge ids as seeds, or 'all'; "
          "the notation x..y is interpreted as all edges between x and y; "
          "the notation random:n is also allowed, to pick n random edges, "
          "and random:n:l, to pick n random edges of size in kmers at least l");
     CommandArgument_Int_OrDefault_Doc(RANDOM_SEED, -1,
          "set set for random number generator");
     CommandArgument_String_OrDefault_Doc(SEEDS_MINUS, "", "exclude these seeds");
     CommandArgument_Int_Doc(DEPTH, "expand to this depth");
     CommandArgument_Bool_OrDefault_Doc(GREEN_SEEDS, True, "color seeds green");
     EndCommandArguments;

     // Load data.

     cout << Date( ) << ": loading" << endl;
     HyperBasevector hb;
     if ( isReadable(IN_HEAD+".hbv") ) BinaryReader::readFile( IN_HEAD+".hbv", &hb );
     else FatalErr("Can't find your input file.");
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);

     // Find seeds, build neighborhood, and color edges.

     vec<int> seeds;
     ParseSeeds( hb, to_right, SEEDS, RANDOM_SEED, SEEDS_MINUS, seeds );
     vec<Bool> invisible;
     BuildNhood( hb, to_left, to_right, seeds, DEPTH, invisible );
     vec<String> edge_color( hb.EdgeObjectCount( ) );
     for ( int i = 0; i < seeds.isize( ); i++ )
          if (GREEN_SEEDS) edge_color[ seeds[i] ] = "green";
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    if ( !invisible[e] ) 
          {    if ( hb.EdgeLengthBases(e) == 0 ) edge_color[e] = "brown";    }    }

     // Name edges.

     cout << Date( ) << ": defining names" << endl;
     vec<int> used;
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          if ( !invisible[e] ) used.push_back(e);
     vec<String> edge_names( hb.EdgeObjectCount( ) );
     for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
     {    if ( !BinMember( used, i ) ) continue;
          if ( hb.EdgeLengthBases(i) == 0 ) continue;
          edge_names[i] = ToString(i);    }

     // Generate output.

     {    Ofstream( fout, OUT_HEAD + ".fasta" );
          for ( int i = 0; i < used.isize( ); i++ )
          {    int e = used[i];
               hb.EdgeObject(e).Print( fout, 
                    ToString(e) + ": " + ToString(to_left[e]) + " --> " 
                    + ToString(to_right[e]) );    }    }
      vec<double> lengths( hb.EdgeObjectCount( ) );
      for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
           lengths[i] = hb.EdgeLengthKmers(i);
      Ofstream( dout, OUT_HEAD + ".dot" );
      cout << Date( ) << ": calling prettyDOT" << endl;
      hb.PrettyDOT( dout, lengths, HyperBasevector::edge_label_info(
          HyperBasevector::edge_label_info::DIRECT, &edge_names ), True, False, 
          NULL, NULL, NULL, NULL, &invisible, &edge_color, NULL, "", hb.K( ) + 50 );
     Scram(0);    }
