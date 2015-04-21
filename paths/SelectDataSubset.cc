/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


/**
   Program: SelectDataSubset

   
   Extracts a set of specified reads from fastb, qualb, and pairto(b,_index)
   files. If EXPAND_PAIRED option is set True then partner paired reads are
   are extracted. Otherwise, only those specified reads that also have their 
   partner (pair) reads specified will be extracted.

   Required Inputs:
   - reads file in fastb format (default: reads.fastb)
   - pairto (or pairtob ) file 
   - .lanes, .lane_infos, .lane_index, .lengths files
   - flat file containing read ids to be selected seperated by white space
  
   Assumed input:

   - pairing informations (reads.pairto)

   Output:
   (reads).fastb, .qualb, .qualb, .pairto, .pairtob, .pairto_index
   .lanes, .lane_infos, .lane_index in the DIR_OUT directory
*/

#include "MainTools.h"
#include "FeudalMimic.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "ReadPairing.h"
#include "CommonSemanticTypes.h"
#include "solexa/SolexaMetrics.h"
#include "math/Functions.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
     RunTime( );
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(IDS);
     CommandArgument_String(DATA);
     CommandArgument_String_Doc(DIR_OUT, "output directory (PRE will be added to it)");
     CommandArgument_Bool_OrDefault(EXPAND_PAIRED, True);
     CommandArgument_String_OrDefault_Doc( READS, "reads",
	  "prefix of the file containing the reads. normally the file is reads.fastb" );
          
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     EndCommandArguments;

     /// Load data.
     String ReadsFile = PRE + "/" + DATA + "/" + READS + ".fastb";
     String QualsFile = PRE + "/" + DATA + "/" + READS + ".qualb";
     String PairsFile = PRE + "/" + DATA + "/" + READS + ".pairto";
     String IdsFile;
     if ( IsRegularFile(PRE + "/" + DATA + "/" + IDS ) )
       IdsFile = PRE + "/" + DATA + "/" + IDS;
     else if ( IsRegularFile( IDS ) )
       IdsFile = IDS;
     String LaneInfosFile  = PRE + "/" + DATA + "/" + READS + ".lane_infos";
     String LanesFile      = PRE + "/" + DATA + "/" + READS + ".lanes";
     String LaneIndexFile  = PRE + "/" + DATA + "/" + READS + ".lane_index";
     String PropertiesFile = PRE + "/" + DATA + "/" + READS + ".props";
     
     String sub_dir = PRE + "/" + DIR_OUT;
     Mkdir777( sub_dir );

     String ReadsFileOut      = PRE + "/" + DIR_OUT + "/" + READS + ".fastb";
     String QualsFileOut      = PRE + "/" + DIR_OUT + "/" + READS + ".qualb";
     String LaneInfosFileOut  = PRE + "/" + DIR_OUT + "/" + READS + ".lane_infos";
     String LanesFileOut      = PRE + "/" + DIR_OUT + "/" + READS + ".lanes";
     String LaneIndexFileOut  = PRE + "/" + DIR_OUT + "/" + READS + ".lane_index";
     String LengthsFileOut    = PRE + "/" + DIR_OUT + "/" + READS + ".lengths";
     String PropertiesFileOut = PRE + "/" + DIR_OUT + "/" + READS + ".props";

     ForceAssertNe( ReadsFile, ReadsFileOut );


     
     int nreads    = MastervecFileObjectCount( ReadsFile );
     cout << "number of reads in the original file = " << nreads << endl;
     

     vec<Bool> ids2keep( nreads, False );
     Ifstream( idin, IdsFile );
     while(1){    
       int n;
       idin >> n;
       if ( !idin ) break;
       ids2keep.at( n ) = True;    
     }

     
     long npairs;
     vec<read_pairing> pairs;
     cout << Date() << " Loading pair information..." << endl;
     ReadPairsFile( PairsFile, pairs );
     npairs = pairs.isize();
     cout << Date() << " Loaded " << pairs.isize() << " pairs." << endl;

     cout << Date() << " processing pairs" << endl;
     for ( int i = 0; i < pairs.isize( ); i++ ){   
       int id1 = pairs[i].id1, id2 = pairs[i].id2;
       if ( ids2keep[ id1 ] && ! ids2keep[ id2 ] ){
	 if ( EXPAND_PAIRED )
	   ids2keep[ id2 ] = True;
	 else 
	   ids2keep[ id1 ] = False;
       }
       if ( ! ids2keep[ id1 ] && ids2keep[ id2 ] ){
	 if ( EXPAND_PAIRED )
	   ids2keep[ id1 ] = True;
	 else
	   ids2keep[ id2 ] = False;
       }
     }


     vec<int>  old2new_ids( nreads, -1 );
     int nid = -1;
     for ( int i = 0; i < nreads; i++ ){
       if ( ids2keep[ i ] ){
	 nid++;
	 old2new_ids.at( i ) = nid; // to make sure that IDS are not out of range
       }
     }
     int nreads_select = nid +1;
     cout << "number of reads selected = " << nreads_select << endl;

     vec<read_pairing> select_pairs;
     select_pairs.reserve( pairs.isize() );
     for ( int i = 0; i < pairs.isize( ); i++ ){   
       int oid1 = pairs[i].id1, oid2 = pairs[i].id2;
       if ( ids2keep[ oid1 ] && ids2keep[ oid2 ] ){
	 int nid1 = old2new_ids[ oid1 ];
	 int nid2 = old2new_ids[ oid2 ];
	 select_pairs.push_back( read_pairing( nid1, nid2, pairs[i].sep, pairs[i].sd ) );
       }
     }
     WritePairs( PRE + "/" + DIR_OUT, select_pairs, nreads_select, False, READS );
     pairs.resize( 0 );
     select_pairs.resize( 0 );

     int nlanes = FirstLineOfFile( LanesFile ).Int( );
     vec<int> lane_size( nlanes, 0 );

     vec<int> index_in, index_out( nreads_select );
     BinaryReader::readFile ( LaneIndexFile, &index_in );
     ForceAssertEq( index_in.isize(), nreads );
     for ( int oid = 0; oid < nreads; oid++){
       int nid = old2new_ids[ oid ];
       if ( nid == -1 ) continue;
       index_out[ nid ] = index_in[ oid ];
       lane_size[ index_in[ oid ] ]++;
     }
     BinaryWriter::writeFile( LaneIndexFileOut, index_out );


     vec<solexa_metric_db> lane_infos;
     solexa_metric_db::ReadMetricsMult( LaneInfosFile, lane_infos );
     ForceAssert( lane_infos.isize( ) == nlanes );
     for ( int i = 0; i < nlanes; i++ )
       lane_infos[i].SetValue( "NREADS", lane_size[i] );
     
     solexa_metric_db::WriteMetricsMult( LaneInfosFileOut, lane_infos );

     Cp( LanesFile, LanesFileOut );
     Cp( PropertiesFile, PropertiesFileOut );

     cout << Date() << " selecting reads and quals" << endl;
     vecbasevector select_reads( nreads_select );
     vecqualvector select_quals( nreads_select );
     vec<int> read_lengths( nreads_select );
     int batchsize = 10000000;
     int nbatches = nreads / batchsize;
     if ( nreads % batchsize > 0 ) nbatches++;
     int si = -1; 
     int bct = 0;
     for ( int b = 0; b < nreads; b += batchsize ){
       bct++;
       if ( VERBOSE )
	 cout << Date() << " processing batch " << bct  << " out of " << nbatches << endl;
       vecbasevector readsTmp;
       vecqualvector qualsTmp;
       int start = b, stop = Min( nreads, b + batchsize );            
       readsTmp.ReadRange( ReadsFile, start, stop );
       qualsTmp.ReadRange( QualsFile, start, stop );
       for ( size_t j = 0; j < readsTmp.size( ); j++ )
	 if ( ids2keep[ b + j ] ){
	   si++;
	   select_reads[ si ] = readsTmp[ j ];
	   select_quals[ si ] = qualsTmp[ j ];
	   read_lengths[ si ] = readsTmp[ j ].isize();
	 }
     }
     
     cout << Date() << " writing reads and quals files" << endl;
     select_reads.WriteAll( ReadsFileOut );
     select_quals.WriteAll( QualsFileOut );
     BinaryWriter::writeFile( LengthsFileOut, read_lengths );
}
