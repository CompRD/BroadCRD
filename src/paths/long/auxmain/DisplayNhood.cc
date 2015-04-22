///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// DisplayNhood: display a neighborhood of a SupportedHyperBasevector, seeded on 
// given edges (to be rendered green) and extended to given depth in the graph.  
// Boundary vertices are labeled red.

// MakeDepend: dependency QueryLookupTable

#include "Basevector.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/DisplayTools.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector.h"
#include <fstream>

int main(int argc, char *argv[])
{
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_String_Doc(IN_HEAD, "looks for IN_HEAD.[s]hbv");
     CommandArgument_String_Doc(OUT_HEAD, "generate OUT_HEAD.{dot,fasta}");
     CommandArgument_String_Doc(SEEDS, "use these edge ids as seeds, or 'all'; "
          "the notation x..y is interpreted as all edges between x and y; "
          "the notation random:n is also allowed, to pick n random edges, "
          "and random:n:l, to pick n random edges of size in kmers at least l");
     CommandArgument_Int_OrDefault_Doc(RANDOM_SEED, -1,
          "set set for random number generator");
     CommandArgument_String_OrDefault_Doc(SEEDS_MINUS, "", "exclude these seeds");
     CommandArgument_Int_Doc(DEPTH, "expand to this depth");
     CommandArgument_Bool_OrDefault_Doc(HIDE_INV, True, "hide inv edges");
     CommandArgument_String_OrDefault_Doc(LOOKUP, "",
          "lookup file for reference; if provided, align edges to reference and "
          "show summary in the dot file");
     CommandArgument_String_OrDefault_Doc(TMP, "", 
          "temporary directory, needed for lookup option");
     CommandArgument_Bool_OrDefault_Doc(GREEN_SEEDS, True, "color seeds green");
     EndCommandArguments;

     SupportedHyperBasevector shb;
     bool isSupported = false;

     cout << Date( ) << ": loading" << endl;
     if ( isReadable(IN_HEAD+".shbv") )
     {
         BinaryReader::readFile( IN_HEAD+".shbv", &shb );
         isSupported = true;
     }
     else if ( isReadable(IN_HEAD+".hbv") )
         BinaryReader::readFile( IN_HEAD+".hbv",
                                 static_cast<HyperBasevector*>(&shb) );
     else
         FatalErr("Can't find your input file.");

     HyperBasevector const& hbv = shb;
     cout << Date( ) << ": setting up" << endl;
     vec<int> to_left, to_right;
     hbv.ToLeft(to_left);
     hbv.ToRight(to_right);

     // Find seeds, build neighborhood, and color edges.

     vec<int> seeds;
     ParseSeeds( hbv, to_right, SEEDS, RANDOM_SEED, SEEDS_MINUS, seeds );
     vec<Bool> invisible;
     BuildNhood( hbv, to_left, to_right, seeds, DEPTH, invisible );
     vec<String> edge_color( hbv.EdgeObjectCount( ) );
     for ( int i = 0; i < seeds.isize( ); i++ )
          if (GREEN_SEEDS) edge_color[ seeds[i] ] = "green";
     for ( int e = 0; e < hbv.EdgeObjectCount( ); e++ )
     {    if ( !invisible[e] ) 
          {    if ( hbv.EdgeLengthBases(e) == 0 ) edge_color[e] = "brown";    }    }

     long_logging logc( "" );

     cout << Date( ) << ": defining names" << endl;
     vec<int> used;
     for ( int e = 0; e < hbv.EdgeObjectCount( ); e++ )
          if ( !invisible[e] ) used.push_back(e);

     {    Ofstream( fout, OUT_HEAD + ".fasta" );
          for ( int i = 0; i < used.isize( ); i++ )
          {    int e = used[i];
               hbv.EdgeObject(e).Print( fout, 
                    ToString(e) + ": " + ToString(to_left[e]) + " --> " 
                    + ToString(to_right[e]) );    }    }

     vec<String> edge_names( hbv.EdgeObjectCount( ) );
     for ( int i = 0; i < hbv.EdgeObjectCount( ); i++ )
     {    if ( !BinMember( used, i ) ) continue;
          if ( hbv.EdgeLengthBases(i) == 0 ) continue;
          if ( !isSupported || !shb.InvDef(i) || !HIDE_INV )
              edge_names[i] = ToString(i);
          else
              edge_names[i] = ToString(i) + "=" + ToString( shb.Inv(i) ) + "'";    }

     if ( LOOKUP != "" )
     {    ForceAssert( TMP != "" );
          vecbasevector edges;
          for ( int e = 0; e < hbv.EdgeObjectCount( ); e++ )
               if ( !invisible[e] ) edges.push_back( hbv.EdgeObject(e) );
          edges.WriteAll( TMP + "/edges.fastb" );
          SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS=" + TMP
               + "/edges.fastb L=" + LOOKUP + " PARSEABLE=True > " + TMP 
               + "/edges.aligns" );
          vec<look_align> aligns;
          vec< vec<int> > aligns_index;
          LoadLookAligns( 
               TMP + "/edges.aligns", aligns, aligns_index, edges.size( ) );
          vec<String> labels( edges.size( ) );
          for ( int i = 0; i < (int) edges.size( ); i++ )
          for ( int j = 0; j < aligns_index[i].isize( ); j++ )
          {    if ( j > 0 ) labels[i] += ",";
               const look_align& la = aligns[ aligns_index[i][j] ];
               labels[i] += ToString(la.target_id) + ( la.Fw1( ) ? "fw" :"rc" )
                    + ":" + ToString( la.pos2( ) )
                    + "-" + ToString( la.Pos2( ) );    }

          vec<String> edge_names( hbv.EdgeObjectCount( ) );
          for ( int i = 0; i < hbv.EdgeObjectCount( ); i++ )
          {    if ( !BinMember( used, i ) ) continue;
               int p = BinPosition( used, i );
               if ( labels[p] != "" ) edge_names[i] += " = " + labels[p];    }    }

      vec<Bool> hide;
      if (isSupported && HIDE_INV)
          FlagEdgesForHiding( shb, shb.Inv( ), hide, logc );
      else hide.resize( hbv.EdgeObjectCount( ), False );
      const Bool DOT_LABEL_CONTIGS = True;
      const Bool DOT_LABEL_VERTICES = False;
      vec<double> lengths( hbv.EdgeObjectCount( ) );
      for ( int i = 0; i < hbv.EdgeObjectCount( ); i++ )
           lengths[i] = hbv.EdgeLengthKmers(i);
      Ofstream( dout, OUT_HEAD + ".dot" );
      cout << Date( ) << ": calling prettyDOT" << endl;
      hbv.PrettyDOT( dout, lengths, HyperBasevector::edge_label_info(
           HyperBasevector::edge_label_info::DIRECT, &edge_names ),
           DOT_LABEL_CONTIGS, DOT_LABEL_VERTICES, NULL, NULL, NULL, &hide,
           &invisible, &edge_color, NULL, logc.LAYOUT, hbv.K( ) + 50 );
     Scram(0);    }
