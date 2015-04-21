///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"

#include "Basevector.h"

int main( int argc, char* argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(FASTB);
     CommandArgument_String_OrDefault(LCONTEXT,"");
     CommandArgument_String_OrDefault(RCONTEXT,"");
     CommandArgument_UnsignedInt_OrDefault(MIN_LEN,1);
     CommandArgument_String(REPEAT);

     EndCommandArguments;


     vecbasevector reads(FASTB);
     std::map<size_t,size_t> counts;

//     LCONTEXT += REPEAT;        // must find at least one copy and this ensures that
                                // the match is "right-aligned"

     for ( size_t i = 0; i < reads.size(); ++i ) {
         auto& r1 = reads[i];
//         if ( i == 417 ) cout << "READ " << i << ": "<< r1.ToString() << endl;

         for ( auto& r : { r1, ReverseComplement(r1) } ) {
             std::string s( r.ToString() );

             auto itr = s.begin();
             auto end = s.end();

	     if ( LCONTEXT != "" ) {
		 // use LCONTEXT to define the beginning of the repeat
		 // locus
		 itr = std::search(itr, end, LCONTEXT.begin(), LCONTEXT.end() );
		 if ( itr == end ) continue;
                 advance(itr,LCONTEXT.size());
	     }

	     if ( RCONTEXT != "" ) {
		 // use RCONTEXT to define the end of the repeat locus
		 // and move to next read if not found
		 end = std::search( itr, end, RCONTEXT.begin(), RCONTEXT.end() );
		 if ( end == s.end() ) continue;
	     }


	     vec<size_t> these_counts;
             while ( itr != end ) {

		 itr = std::search( itr, end, REPEAT.begin(), REPEAT.end() );
		 if ( itr == end ) break;
		 advance(itr, REPEAT.size());

                 // count repeats
                 size_t count = 1;
                 while ( itr < end ) {    // ick!
                     if ( !std::equal(REPEAT.begin(), REPEAT.end(), itr ) )
                         break;
                     count++;
                     advance(itr, REPEAT.size() );
                 }

//                 if ( count > MIN_LEN ) counts[count]++;
                 these_counts.push_back(count);
             }
	     // now just counts the *longest* repeat in the read
	     size_t count = 0;
	     if ( these_counts.size() ) count = Max(these_counts);
	     if ( count > MIN_LEN ) counts[count]++;

         }

     }

     for ( auto& c : counts )
         cout << REPEAT << "(" << c.first << ") happens " << c.second << endl;

     return 0;
}
