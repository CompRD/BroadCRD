///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SearchFastb2: Given fastb files F1 and F2, find all sequences in F1 that align 
// perfectly from end to end to a sequence in F2.  Ignores sequences of length < K.
// Parallelized.
//
// The allowed values of K are hardwired, but other values could easily be added.
//
// Ways to make the code run faster:
// - If all sequences in F1 have the same size, set K to that value.
// - Use the MAX_PLACEMENTS argument.
//
// Output:
// ALIGN id1 id2 orientation start stop
// where orientation is + (forward) or - (reverse) and [start,stop) gives the
// position of F1[id1] or rc(F1[id1]) on F2[id2].

#include "Basevector.h"
#include "MainTools.h"
#include "util/SearchFastb2Core.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(F1);
     CommandArgument_String(F2);
     CommandArgument_String_OrDefault_Doc(Q2, "", "qualb associated with F2");
     CommandArgument_Int(K);
     CommandArgument_Int_OrDefault_Doc(MAX_PLACEMENTS, -1,
          "If specified, and if some F1[i] has more than MAX_PLACEMENTS placements, "
          "return no placements of it.");
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Double_OrDefault_Doc(MEM_FRAC_TO_USE, 0.5, 
          "Try to use at most this fraction of total memory.");
     CommandArgument_Bool_OrDefault(VERBOSE, True);
     EndCommandArguments;

     vec< triple<int64_t,int64_t,int> > ALIGNS;
     SearchFastb2( F1, F2, K, &ALIGNS, 0, MAX_PLACEMENTS, MEM_FRAC_TO_USE, VERBOSE );
     VirtualMasterVec<QualVec>* q2p = 0;
     if ( Q2 != "") q2p = new VirtualMasterVec<QualVec>(Q2);

     if (WRITE)
     {    if (VERBOSE) cout << Date( ) << ": reloading F1" << endl;
          vecbasevector f1(F1);
          if (VERBOSE) cout << Date( ) << ": printing alignments" << endl;
          for ( int i = 0; i < ALIGNS.isize( ); i++ ) 
          {    bool fw = ( ALIGNS[i].third >= 0 );
               int begin = fw ? ALIGNS[i].third : - ALIGNS[i].third - 1;
               int end = begin + f1[ ALIGNS[i].first ].size( );
               size_t target = ALIGNS[i].second;
               cout << "ALIGN " << ALIGNS[i].first << " " 
                    << target << " " << ( fw ? "+" : "-" )
                    << " " << begin << " " << end;

               if ( q2p ) {
		   long qsum = 0;
		   for ( int j = begin; j < end; ++j )
		       qsum += (*q2p)[target][j];

		   cout << " QSUM " << qsum << " QMEAN " 
                        << qsum / (end-begin) 
                        << " QMID " << int((*q2p)[target][ begin + (end-begin)/2 ]);
                         ;
               }

               cout << "\n";
          }
     }
     if (VERBOSE) cout << Date( ) << ": done" << endl;


     delete q2p;
}
