/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: SNPRateVsUnipaths

   A code for simulating the effect of a SNP rate on lengths of unipaths with copy
   numbers 1 or 2.


   @file
*/

#include "Basevector.h"
#include "MainTools.h"
#include "Vec.h"
#include "random/Random.h"
#include "paths/KmerPath.h"
#include "paths/HyperKmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include <time.h>

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandDoc("Simulates SNPs and resultant unipaths.");

     CommandArgument_String_OrDefault_Doc(FASTB_IN, "",
            "Haploid genome fastb file.");
     CommandArgument_String_Doc(FASTB_OUT, "Diploid genome fastb file.");
     CommandArgument_Double_Doc(SUB_RATE, "Per base mutation rate.");
     CommandArgument_Int_Doc(K, "Kmer size.");
     CommandArgument_Int_OrDefault_Doc(N_SIMUL, 10,
	    "Number of independent simulations.");
     CommandArgument_Int_OrDefault_Doc(FASTB_SIZE, 10000000,
            "Simulated haploid genome fastb size.");
    
     EndCommandArguments;

    

     vecbasevector b1;
     if ( FASTB_IN.nonempty() ){
        cout << "Loading haploid genome file:\n" << FASTB_IN << endl;
       b1.ReadAll(FASTB_IN);       
     }else{
       cout << "Simulating haploid genome" << endl;
       String bs( FASTB_SIZE );
       for ( size_t i = 0; i < bs.size(); i++ ){
	 int val = randomx() % 4;
	 if ( val == 0 ) bs[i] = 'A';
	 else if ( val == 1 ) bs[i] = 'T';
	 else if ( val == 2 ) bs[i] = 'C';
	 else if ( val == 3 ) bs[i] = 'G';
	 else{ FatalErr("ERROR"); }
       }
       basevector b( bs );
       b1.push_back(b);
     }

     // structures for holding results
     vec<float> ratios, first_cn1s, first_high_cns, first_cn2s;

     // run simulations
     for ( int l = 1; l <= N_SIMUL; l++ ){  
       srandomx( l );
       vecbasevector b2;
       
       for ( size_t i = 0; i < b1.size(); i++ )
	 b2.push_back( b1[i] );
      
       
       cout << "Haploid genome size: " << b1.sumSizes() << endl;
       
       cout << "Creating SNPs (" << SUB_RATE * 100.0 
	    << "% polymorphic = 1 SNP per " << 1/SUB_RATE << " bases)" << endl;

       
       vec< pair<int,int> > snp_locs;
       
       int genome_size = 0;
       int M = 1000000;
       int m = int( round( double(M) * SUB_RATE ) );
       for ( size_t i = 0; i < b2.size( ); i++ ){
	 genome_size += b2[i].isize();
	 for ( int j = 0; j < b2[i].isize( ); j++ )
	   if ( randomx( ) % M <= m ) {
	     b2[i].Set( j, ( b2[i][j] + randomx( ) % 3 + 1 ) % 4 );
	     snp_locs.push_back( make_pair( i, j ) );
	     snp_locs.push_back( make_pair( b1.size() + i, j ) );
	   }
       }
       
       vecbasevector b = b1;
       b.Append(b2);
       genome_size = b.sumSizes();
       
       // Generate unipaths.
       
       vecKmerPath paths, paths_rc, unipaths;
       vec<big_tagged_rpint> pathsdb, unipathsdb;
       cout << Date( ) << ": calling RTP" << endl; 
       ReadsToPathsCoreY( b, K, paths, paths_rc, pathsdb );
       cout << Date( ) << ": calling Unipath" << endl; // XXX
       Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb );
       digraph AG;
       BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths, unipathsdb, AG );
       

       vec<int> ulen( unipaths.size( ) );
       for ( size_t id = 0; id < unipaths.size( ); id++ )
	 ulen[id] = K - 1 + unipaths[id].KmerCount( );
       
       vec<unipath_id_t> to_rc( unipaths.size() );
       UnipathInvolution( unipaths, unipathsdb, to_rc );
       vec<short> u2use( unipaths.size(), 1 );
       PRINT( Sum( u2use ) );
       for ( size_t id = 0; id < unipaths.size(); id++ ){
	 if ( u2use[id] == 0 || (int)id == to_rc[id] ) continue;
	 u2use[ to_rc[id] ] = 0;
       }
       PRINT( Sum( u2use ) );

       // The thing we're going to generate:
       vec<int> ucns( unipaths.size(), 0 );
       vec<unipath_interval_id_t> hits;
       for( size_t i=0; i < unipaths.size(); i++ ) {
	 if( unipaths[i].IsEmpty() )
	   hits.clear();
	 else
	   Contains( pathsdb, unipaths[i].Start(0), hits );
	 ucns[i] = hits.size();
       }
       
       vec<int> uord( unipaths.size(), vec<int>::IDENTITY );
       vec<int> sulen( ulen );
       ReverseSortSync( sulen, uord );
       //ulen.Print(cout); cout << endl;

       

       //ucns.Print(cout); cout << endl;
       double genome_uni_size = 0;
       double est_genome_size = 0;
       
       for ( size_t i = 0; i < uord.size(); i++ ){
	 int uid = uord[i];
	 ForceAssertGt( ucns[uid], 0 );
	 if ( ! u2use[uid] ) continue;
	 genome_uni_size += ulen[uid] -K + 1;
	 est_genome_size += (ulen[uid] - K +1) * ucns[uid];
	 if ( AG.ToMutable(uid).size() == 0 )
	   est_genome_size += K -1;
       }

       double uni_cn1_size = 0, uni_cumul_size = 0;
       double first_high_cn = -1, first_cn1 = -1, first_cn2 = -1;
       
       for ( size_t i = 0; i < uord.size(); i++ ){
	 int uid = uord[i];
	 if ( ! u2use[uid] ) continue;
	 if ( ucns[uid] == 1 ){
	   uni_cn1_size += ulen[uid] -K +1;
	   if ( first_cn1 == -1 )
	     first_cn1 = (double) uni_cumul_size / (double) genome_uni_size * 100.0;
	 }
	 if ( ucns[uid] > 2 ){
	   PRINT3( "HIGH CN", ucns[uid], ulen[uid] );
	   if ( first_high_cn == -1 )
	     first_high_cn = (double) uni_cumul_size / (double) genome_uni_size * 100.0;
	 }
	 if ( ucns[uid] == 2 ){
	   if ( first_cn2 == -1 )
	     first_cn2 = (double) uni_cumul_size / (double) genome_uni_size * 100.0;
	 }
	 uni_cumul_size += ulen[uid] - K +1;
       }
       
       double cnratio = (double)uni_cn1_size / (double) genome_uni_size * 100.0;
       PRINT5( genome_size, est_genome_size, genome_uni_size, uni_cn1_size, cnratio );
       ratios.push_back( cnratio );
       PRINT3( first_cn1, first_high_cn, first_cn2 );
       first_cn1s.push_back( first_cn1 );
       first_high_cns.push_back( first_high_cn );
       first_cn2s.push_back( first_cn2 );
     } 
     
     NormalDistribution ratio = SafeMeanStdev( ratios );
     PRINT3( "ratios", ratio.mu_, ratio.sigma_ );
     NormalDistribution cn1s = SafeMeanStdev( first_cn1s );
     PRINT4( "first_cn1", cn1s.mu_, cn1s.sigma_, Min( first_cn1s) );
     NormalDistribution high_cns = SafeMeanStdev( first_high_cns );
     PRINT4( "firt_high_cn", high_cns.mu_, high_cns.sigma_, Min( first_high_cns) );
     NormalDistribution cn2s = SafeMeanStdev( first_cn2s );
     PRINT4( "first_cn2", cn2s.mu_, cn2s.sigma_, Min( first_cn2s ) );
     
     cout << Date() << ": done with SNPRateVsUnipaths!" << endl;
}
