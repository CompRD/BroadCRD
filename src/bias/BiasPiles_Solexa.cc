/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// BiasPiles_Solexa.  Show piles of read placements that drive the bias metric,
/// showing the coordinates of each cluster in the pile and the first cycle
/// intensity sums.  Use the first N bases from each read.
///
/// Also compute bias, after discarding proximate clusters, for various
/// definitions of proximate.
///
/// Three metrics are defined:
/// * opt_dup_27b
/// * meta_dup_27b
/// * bias_dup_27b (no longer computed by default).

#include "Basevector.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "feudal/BinaryStream.h"
#include "feudal/VirtualMasterVec.h"
#include "bias/StartBias.h"
#include "lookup/PerfectCount.h"
#include "math/Functions.h"
#include "random/FindRandomReads.h"
#include "random/Random.h"
#include "solexa/FourBase.h"
#include "solexa/SolexaTools.h"
#include "solexa/SolexaMetrics.h"
#include "solexa/SolexaPosition.h"
#include "solexa/SolexaTools.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(HEAD, "flowcell.lane");
     CommandArgument_String_OrDefault_Doc(GENOME,"",
          "Use reference GENOME.reference.lookup; by default GENOME=HEAD.");
     CommandArgument_Int_OrDefault_Doc(N, 27, "length to truncate reads to");
     CommandArgument_String_OrDefault_Doc(REF, "", 
          "Use alternative reference REF.lookup "
          "(legacy/use with non-conforming reference file names)");
     CommandArgument_Bool_OrDefault_Doc(WRITE, True, "to write metrics");
     CommandArgument_Bool_OrDefault_Doc(VERBOSE, True, "to print interesting stuff");
     CommandArgument_Bool_OrDefault_Doc(BIAS_DUP, False, 
          "to compute bias_dup metrics");
     EndCommandArguments;

     // Load metrics.  Quit if there is no reference or genome is big.

     solexa_metric_db db( HEAD + ".metrics" );
     if ( GENOME != "" && REF != "" ) {
       cout << "Illegal arguments: setting both REF and GENOME is not allowed" 
            << endl;
     }
     if ( GENOME != "" ) REF=GENOME+".reference";
     else {
       if ( REF == "" ) REF=HEAD+".reference";
     }
 
     Bool have_ref = IsRegularFile( REF +".lookup" );
     if ( !have_ref ) {
       cout << "Reference file " << REF << ".lookup does not exist" << endl;
       exit(0);
     }
     ForceAssert( db.Defined( "genome_size" ) );
     if ( db.Value( "genome_size" ).Int( ) > 100 * 1000 * 1000 ) exit(0);

     // Load bases.

     vecbasevector bases( HEAD + ".fastb" );
     vecbasevector bases_orig;
     if (VERBOSE) bases_orig = bases;

     // Shrink to desired size. 
     for ( size_t i = 0; i < bases.size( ); i++ )
          bases[i].resize(N);

     // Kill noise reads.
     /*
     {
       int HomopolReads = 0;
       VEC_FROM_PARAMS_FILE(HEAD+".paramsByRead", zScoreDecayHeight);
       for ( int i = 0; i < bases.size( ); i++ ) {
	 if (bases[i].HomopolPercent().first>90  && zScoreDecayHeight[i] >= 2 ) {
	   ++HomopolReads;
	   bases[i].resize(0);
	 }
       }
       PRINT(HomopolReads);
     }
     */

     // Define Intensity1.

     vec<float> Intensity1( bases.size( ) );
     if ( true )
     {    
       VirtualMasterVec<FourBaseVec> I( (HEAD + ".intensities").c_str() );
       for ( size_t i = 0; i < bases.size( ); i++ )
	 Intensity1[i] = I[i][0].Sum( );    
     }

     // Load more data.

     vec<SolexaPosition> pos;
     BinaryReader::readFile( HEAD + ".position", &pos );
     String refname = REF ;
     vecbasevector ref( refname + ".fastb" );
     vecbasevector rref = ref;
     ReverseComplement(rref);

     // Align reads.

     vec< pair<placement_mark,int> > places;
     vec<Bool> qPerfect;
     PerfectPick( bases, refname + ".lookup", FW_OR_RC, places, qPerfect );
     Sort(places);

     // Report first base stats.

     if (VERBOSE) {    
         vec<int> base( 4, 0 ), rbase( 4, 0 );
	 vec<double> basefrac(4), rbasefrac(4);
	 for ( size_t i = 0; i < ref.size( ); i++ ) {
	    for ( unsigned int j = 0; j < ref[i].size( ); j++ ) {
	      ++rbase[ ref[i][j] ];
	      ++rbase[ 3 - ref[i][j] ];    
	    }    
	  }
          cout << "bases of reference: {";
          for ( int i = 0; i < 4; i++ )  {    
	      if ( i > 0 ) cout << ",";
	      rbasefrac[i] = double( rbase[i] ) / double( Sum(rbase) );
	      cout << ToString( 100.0 * rbasefrac[i], 2 ) + "%";    
	  }
          cout << "}\n";
          for ( int i = 0; i < places.isize( ); i++ )  {    
	    int tig = places[i].first.Tig( );
	    Bool fw1 = places[i].first.Fw1( );
	    int pos = places[i].first.Pos( );
	    int id = places[i].second;
	    if (fw1) ++base[ ref[tig][pos] ];
	    else ++base[ 3 - ref[tig][ pos + bases[id].size( ) - 1 ] ];
	  }
          cout << "first base of perfect reads: {";
          for ( int i = 0; i < 4; i++ )  {    
	    if ( i > 0 ) cout << ",";
	    basefrac[i] = double( base[i] ) / double( Sum(base) );
	    cout << ToString( 100.0 * basefrac[i], 2 ) + "%";    
	  }
          cout << "}\n";
          vec<double> fracdiff(4);
          double first_base_bias = 0.0;
          for ( int i = 0; i < 4; i++ )  {    
	    fracdiff[i] = basefrac[i] - rbasefrac[i];
	    first_base_bias += fracdiff[i] * fracdiff[i];    
	  }
          first_base_bias = 100.0 * sqrt(first_base_bias);
          cout << "first base bias: " << ToString( first_base_bias, 2 ) 
               << "\n\n";   
     }

     // Find piles.

     vec<int> pilesizes;
     if (VERBOSE)
     {    cout << "id = read id\n";
          cout << "int1 = first cycle intensity\n";
          cout << "tile = tile from which read came\n";
          cout << "xy = coordinates on tile\n";
          cout << "tail = read bases after the first " << N << "\n";
          for ( int i = 0; i < places.isize( ); i++ ) {    
	      int j;
	      for ( j = i + 1; j < places.isize( ); j++ ) {
		  if ( places[j].first != places[i].first ) break;
	      }
	      if ( j - i >= 2 )  {    
		  cout << "\n[" << j-i << "] " << places[i].first;
		  int id0 = places[i].second;
		  for ( int u = 0; u < N; u++ ) {
		      cout << as_base( bases[id0][u] );
		  }
		  cout << "\n";
		  pilesizes.push_back(j-i);
		  for ( int k = i; k < j; k++ ) {    
		      int id = places[k].second;
		      float int1 = Intensity1[id];
		      int tile = pos[id].tile, x = pos[id].x, y = pos[id].y;
		      String tail;
		      String xy = "(" + ToString(x) + "," + ToString(y) + ")";
		      for ( unsigned int u = N; u < bases_orig[id].size( ); u++ ) {
			  tail += as_base( bases_orig[id][u] );
		      }
		      PRINT5( id, int1, tile, xy, tail );    
		  }    
	      }
	      i = j - 1;    
	  }    
     }

     // Compute bias.

     vec<placement_mark> places0;
     for ( int i = 0; i < places.isize( ); i++ )
          places0.push_back( places[i].first );
     if (VERBOSE)
     {    cout << "\nbias:\n" << StartBias( bases, ref, rref, places0, 
               True, "", True, True ) << endl;    }

     // Compute bias after removing proximate duplicates, for various definitions of
     // proximate.

     longlong opt_dups = 0, nplaces = places.size( );
     if ( nplaces == 0 )
     {    if (VERBOSE) cout << "No reads are placed perfectly.\n";
          exit(0);    }
     vec<double> tilesep( db.ValueInt( "tiles" ), 0 );
     vec<double> tilesep_predicted( tilesep.size( ), 0 );
     vec<int> tilesize( tilesep.size( ) + 1, 0 );
     for ( int i = 0; i < places.isize( ); i++ )
          ++tilesize[ pos[ places[i].second ].tile ];
     for ( int i = 0; i < 1000000; i++ )
     {    int j1 = randomx( ) % places.size( ), j2 = randomx( ) % places.size( );
          if ( j1 == j2 ) continue;
          int id1 = places[j1].second, id2 = places[j2].second;
          ++tilesep_predicted[ Abs( pos[id1].tile - pos[id2].tile ) ];    }
     tilesep_predicted /= Sum(tilesep_predicted);
     for ( int pass = 1; pass <= 3; pass++ )
     {    vec<Bool> to_remove( places.size( ), False );
          for ( int i = 0; i < places.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < places.isize( ); j++ )
                    if ( places[j].first != places[i].first ) break;
               for ( int k1 = i; k1 < j; k1++ )
               {    for ( int k2 = k1+1; k2 < j; k2++ )
                    {    int id1 = places[k1].second, id2 = places[k2].second;
                         int tilediff = Abs( pos[id1].tile - pos[id2].tile );
                         int xdiff = Abs( pos[id1].x - pos[id2].x );
                         int ydiff = Abs( pos[id1].y - pos[id2].y );
                         if ( pass == 1 )
                         {    if ( tilediff == 0 
                                   && sqrt( xdiff*xdiff + ydiff*ydiff ) <= 5.0 )
                              {    to_remove[k2] = True;
                                   ++opt_dups;    }
                              else ++tilesep[tilediff];    }
                         if ( pass == 2 && tilediff <= 5 ) to_remove[k2] = True;
                         if ( pass == 3 && tilediff <= 30 )
                              to_remove[k2] = True;    }    }
               i = j - 1;    }
          EraseIf( places, to_remove );
          places0.clear( );
          for ( int i = 0; i < places.isize( ); i++ )
               places0.push_back( places[i].first );
          if (VERBOSE)
          {    if ( pass == 1 )
               {    cout << "\nbias after removing optical duplicates "
                         << "(<= 5 pixels apart):\n";    }
               else if ( pass == 2 )
               {    cout << "\nbias after removing medium-range duplicates "
                         << "(<= 5 tiles apart):\n";    }
               else 
               {    cout << "\nbias after removing longer-range duplicates "
                         << "(<= 30 tiles apart):\n";    }
               cout << StartBias( bases, ref, rref, places0, True, "", True, True ) 
                    << endl;    }    }

     // Estimate probability that read has an optical duplicate.

     double P_opt_dup = double(opt_dups) / double(nplaces);
     String metric1 = "opt_dup_" + ToString(N) + "b";
     String metric1_val = ToString( P_opt_dup * 100.0, 2 ) + "%";
     if (VERBOSE)
     {    cout << "\n" << metric1 << " = chance that read has optical duplicate = "
               << metric1_val << endl;    }

     // Estimate probability that read has a metastasized duplicate.

     vec< pair<double,int> > tsd;
     for ( int j = 0; j < tilesep.isize( ); j++ )
          tsd.push( tilesep[j] / tilesep_predicted[j], j );
     Sort(tsd);
     double high_tile_dist_dups = 0.0, tcount = 0.0;
     int nt = tsd.size( );
     for ( int j = nt/2; j < Min( nt/2 + 10, nt ); j++ )
     {    int k = tsd[j].second;
          high_tile_dist_dups += tilesep[k];
          tcount += tilesep_predicted[k];    }
     vec<double> mtilesep( tilesep.size( ) );
     for ( int i = 0; i < tilesep.isize( ); i++ )
     {    mtilesep[i] = Max( 0.0, 
               tilesep[i] - high_tile_dist_dups * tilesep_predicted[i]/tcount );    }
     String metric2 = "meta_dup_" + ToString(N) + "b";
     String metric2_val 
          = ToString( 100.0 * Sum(mtilesep) / double(nplaces), 2 ) + "%";
     if (VERBOSE)
     {    cout << "\n" << metric2 << " = chance that read has metastasized "
               << "duplicate = " << metric2_val << endl;    }

     // Estimate bias due to optical and metasized duplicates.  This code is now
     // turned off by default, because DUP can be very large, and hence the 
     // subsequent computation can be very expensive.

     if (BIAS_DUP)
     {    vecbasevector reads;
          FindRandomReads F;
          F.Reads( nplaces, ref, reads, db.ValueInt( "read_length" ) );
          double DUP = metric1_val.Double( ) + metric2_val.Double( );
          if (VERBOSE) PRINT(DUP);
          int extra = int( round( DUP * double(nplaces) ) );
          if (VERBOSE) PRINT2( nplaces, extra );
          for ( int i = 0; i < extra; i++ )
               reads.push_back_reserve( reads[ randomx( ) % reads.size( ) ] );
          vec<placement_mark> xplaces;
          PerfectPick( reads, refname +".lookup", FW_OR_RC, xplaces );
          String metric3 = "bias_dup_" + ToString(N) + "b";
          String metric3_val 
               = StartBias( reads, ref, rref, xplaces, True, "", True, True );
          if ( N == 27 ) db.SetValue( metric3, metric3_val );    
          if (VERBOSE)
          {    cout << "\nEstimated bias due to optical and metasized duplicates: ";
               cout << metric3_val << "\n";    }    }

     // Print tile separation histogram, etc.

     if (VERBOSE)
     {    cout << "\nRate of non-optical duplicates (including all other sources):"
               << "\n";
          double rate_sum = 0;
          for ( int i = 0; i < tilesep.isize( ); i++ )
          {    double denom = 0;
               for ( int j = 0; j < tilesize.isize( ) - i; j++ )
                    denom += tilesize[j] * tilesize[i+j];
               double rate = 100.0 * double(tilesep[i]) / double(nplaces);
               rate_sum += rate;
               cout << i << " " << rate << "%\n";    }    
          cout << "Sum = " << rate_sum << "%\n";    }
     tilesep /= Sum(tilesep);
     if (VERBOSE) cout << "\ndistance in tiles between non-optical-duplicates:\n";
     vec< vec<String> > rows;
     vec<String> row0;
     row0.push_back( "distance", "observed", "random" );
     rows.push_back(row0);
     double ts_sum = 0;
     for ( int i = 0; i < tilesep.isize( ); i++ )
     {    vec<String> row;
          row.push_back( ToString(i) );
          row.push_back( ToString( 100.0 * tilesep[i], 2 ) + "%" );
          row.push_back( ToString( 100.0 * tilesep_predicted[i], 2 ) + "%" );
          rows.push_back(row);
          double diff = tilesep[i] - tilesep_predicted[i];
          ts_sum += diff*diff;    }
     if (VERBOSE) PrintTabular( cout, rows, 2, "rrr" );
     if (VERBOSE)
     {    cout << "\n" << metric2 << " = sqrt(sum of squares of differences) = "
               << metric2_val << endl;    }

     // Describe pile sizes.

     if (VERBOSE)
     {    Sort(pilesizes);
          cout << "\npile-size count-of-piles:\n";
          for ( int i = 0; i < pilesizes.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < pilesizes.isize( ); j++ )
                    if ( pilesizes[j] != pilesizes[i] ) break;
               cout << pilesizes[i] << " " << j - i << "\n";
               i = j - 1;    }    }

     // Save metrics.

     if ( N == 27 )
     {    db.SetValue( metric1, metric1_val );
          db.SetValue( metric2, metric2_val );     }
     if (WRITE) db.WriteMetrics( HEAD + ".metrics" );    
}
