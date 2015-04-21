///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Take in input a set of reads (bases, quals, and pairing info), and trim "
  "them to TRIM_SIZE bases, starting from the beginning of the read. If TAIL "
  "is set to true, the bases are selected at the end of the read (the reads "
  "will not be rc-ed either way). Can optionally update pairing and alignment "
  "information.";

#include "Set.h"
#include "Alignment.h"
#include "Basevector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "PairsManager.h"
#include "ReadPairing.h"
#include "PairsManager.h"
#include "lookup/LookAlign.h"

/**
 * TrimReads
 *
 * HEAD_IN: it needs <HEAD_IN>.{fastb,qualb,pairto}
 * HEAD_OUT: it saves <HEAD_OUT>.{fastb,qualb,pairto}
 * TRIM_SIZE: how many bases to keep
 * TRIM_START: first base to keep
 * TAIL: select bases at the end of the read, rather than the beginning 
 *       (overrides TRIM_START)
 * UPDATE_QLTS: whether to update aligns file (.qltout)
 * UPDATE_PAIR_SEPS: whether to update pairto separations
 * PRESERVE_SHORTER: preserve reads that are shorter than TRIM_SIZE instead 
 *                   of replacing them with empty reads
 */

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String( HEAD_IN );
  CommandArgument_String( HEAD_OUT );
  CommandArgument_Int_OrDefault( TRIM_SIZE, -1 );
  CommandArgument_Int_OrDefault( TRIM_START, 0 );
  CommandArgument_Bool_OrDefault( TAIL, False );
  CommandArgument_Bool_OrDefault( UPDATE_QLTS, False );
  CommandArgument_Bool_OrDefault( UPDATE_PAIR_SEPS, False );
  CommandArgument_Bool_OrDefault( PRESERVE_SHORTER, False );
  EndCommandArguments;
  
  // Dir and file names.
  String bases_in = HEAD_IN + ".fastb";
  if ( ! IsRegularFile( bases_in ) ) {
    cout << bases_in << " not found. Abort.\n" << endl;
    return -1;
  }

  String quals_in = HEAD_IN + ".qualb";
  Bool FoundQuals = True;
  if ( ! IsRegularFile( quals_in ) ) {
    FoundQuals = False;
    cout << quals_in << " not found. Trimming bases only" << endl;
  }
  
  String qlt_in = HEAD_IN + ".qltout";
  Bool FoundQlts = UPDATE_QLTS;
  if ( UPDATE_QLTS && ! IsRegularFile( qlt_in ) ) {
    FoundQlts = False;
    cout << qlt_in << " not found. Skipping qlts." << endl;
  }
  
  String out_dir = ".";
  if ( HEAD_OUT.Contains( "/" ) ) out_dir = HEAD_OUT.RevBefore( "/" );
  Mkpath( out_dir );
  
  String info_out = HEAD_OUT + ".log";
  String bases_out = HEAD_OUT + ".fastb";
  String quals_out = HEAD_OUT + ".qualb";
  String qlt_out = HEAD_OUT + ".qltout";
  
  ofstream out( info_out.c_str( ) );
  PrintCommandPretty( out );
  out.close( );

  // Original read lengths will be created (and used) only if TAIL is false.
  vec<int> orig_lens;
  longlong nbases = MastervecFileObjectCount( bases_in );
  if ( ! TAIL ) orig_lens.reserve( nbases );

  // Bases.
  {
    cout << Date( ) << ": trimming bases" << endl;
    vecbvec bases( bases_in );
    bvec empty;
    IncrementalWriter<bvec> bwriter( bases_out.c_str( ) );
    for (ulonglong ii=0; ii<bases.size( ); ii++) {
      const bvec &orig = bases[ii];
      if ( ! TAIL ) orig_lens.push_back( (int)orig.size( ) );
      if ( TRIM_SIZE < 0 ) TRIM_SIZE = orig.isize( ) - TRIM_START;
      if ( orig.isize( ) < TRIM_SIZE ) {
	if ( ! PRESERVE_SHORTER )
	  bwriter.add( empty );
	else
	  bwriter.add( bases[ii] );
	continue;
      }
      size_type start = TAIL ? orig.size() - TRIM_SIZE : TRIM_START;
      bwriter.add( bvec( orig, start, TRIM_SIZE ) );
    }
    bwriter.close( );
  }

  // Quals.
  if ( FoundQuals) {
    cout << Date( ) << ": trimming quals" << endl;
    vecqvec quals( quals_in );
    qvec empty;
    qvec newq;
    IncrementalWriter<qvec> qwriter( quals_out.c_str( ) );
    for (ulonglong ii=0; ii<quals.size( ); ii++) {
      const qvec &orig = quals[ii];
      if ( TRIM_SIZE < 0 ) TRIM_SIZE = (int) orig.size( ) - TRIM_START;
      if ( (int)orig.size( ) < TRIM_SIZE ) {
	if ( ! PRESERVE_SHORTER )
	  qwriter.add( empty );
	else
	  qwriter.add( quals[ii] );
	continue;
      }
      size_type start = TAIL ? orig.size( ) - TRIM_SIZE : 0;
      newq.resize(TRIM_SIZE);
      CopyQuals( orig, start, newq, 0, TRIM_SIZE );
      qwriter.add( newq );
    }
    qwriter.close( );
  }  
  
  // Pairing info.
  if (UPDATE_PAIR_SEPS) {
    cout << Date() << ": updating pair separations (NOTE: this assumes that reads are pointing inwards" << endl;
    PairsManager pairs;
    String pairs_in  = HEAD_IN + ".pairs";
    String pairto_in = HEAD_IN + ".pairto";
    Bool IsPairto = False;
    if ( ! IsRegularFile( pairs_in ) ) {
      if ( ! IsRegularFile( pairto_in ) ) {
	cout << "Neither " << pairs_in << " nor " << pairto_in << " were found. Abort.\n" << endl;
	return -1;
      }
      pairs.ReadFromPairtoFile( pairto_in, nbases );
      IsPairto = True;
    }else{
      pairs.Read( pairs_in );
    }
    
    if ( ! TAIL ) {
      
      vec< StdSet<int> > newLibSeps( pairs.nLibraries() );
      for(ulonglong ii=0; ii<pairs.nPairs( ); ii++) {
	size_t id1    = pairs.ID1( ii );
	size_t id2    = pairs.ID2( ii );
	size_t libID  = pairs.libraryID( ii );
	int origSep   = pairs.getLibrarySep( libID );
        int TRIM_SIZE1 = TRIM_SIZE, TRIM_SIZE2 = TRIM_SIZE;
        if ( TRIM_SIZE < 0 ) TRIM_SIZE1 = orig_lens[id1] - TRIM_START;
        if ( TRIM_SIZE < 0 ) TRIM_SIZE2 = orig_lens[id2] - TRIM_START;
	int amt1      = Max( 0, orig_lens[id1] - TRIM_SIZE1 );
	int amt2      = Max( 0, orig_lens[id2] - TRIM_SIZE2 );
        int newSep    = origSep + amt1 + amt2;
	if ( orig_lens[id1] > 0 && orig_lens[id2] > 0 )
	  newLibSeps[ libID ].insert( newSep );
      }
      
      for ( size_t libID = 0; libID < pairs.nLibraries(); libID++ ){
	if ( newLibSeps[libID].size() != 1 ){
	  cout << "separations for libID=" << libID << ": ";
	  for ( set<int>::iterator i = newLibSeps[libID].begin(); i != newLibSeps[libID].end(); i++ )
	    cout << *i << " ";
	  cout << "\n";
	  FatalErr("Inconsistent separations in read library id " + libID );
	}
	int origSep = pairs.getLibrarySep( libID );
	int origSd  = pairs.getLibrarySD( libID );
	int newSep  = *(newLibSeps[libID].begin());
	PRINT4( libID, origSep, origSd, newSep );
	pairs.changeLibrarySepSd( libID, newSep, origSd );
      }
      
    }
    if ( IsPairto ){
      vec<read_pairing> pairto = pairs.convert_to_read_pairings();
      WritePairs( pairto, pairs.nReads(), HEAD_OUT );
    }
    else
      pairs.Write( HEAD_OUT + ".pairs" );
  }
  
  // Update aligns.
  if ( FoundQlts ) {
    cout << Date( ) << ": updating aligns" << endl;
    vec<look_align_plus> qlts;
    LoadLookAlignPlus( qlt_in, qlts );

    ofstream qout( qlt_out.c_str( ) );
    for (size_t ii=0; ii<qlts.size( ); ii++) {
      look_align_plus &la = qlts[ii];
      int qlen = la.QueryLength( );
      int trim_size = TRIM_SIZE;
      if ( trim_size < 0 ) trim_size = qlen - TRIM_START;
      if ( qlen < trim_size ) {
	if ( PRESERVE_SHORTER )
	  la.WriteParseable( qout );
	continue;
      }
      int begin = TAIL ? qlen - trim_size : TRIM_START;
      if ( la.Rc1( ) ) begin = la.query_length - ( begin + trim_size );
      la.a = la.a.TrimmedTo1( begin, trim_size );
      la.WriteParseable( qout );
    }
    qout.close( );
  }

  // Done.
  cout << Date( ) << ": done" << endl;
  
}

