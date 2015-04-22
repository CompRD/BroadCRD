/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// MakeFake.  Make fake mutations, suitable for SomaticCall.  This will make 
// duplicates at low frequency.

#include "Basevector.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "random/Random.h"

int main(int argc, char **argv) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(REF_FASTA);
     CommandArgument_Int(COUNT);
     CommandArgument_Double(FRAC);
     CommandArgument_String_OrDefault_Doc(CHRS,
          "chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}",
          "list of chromosomes to search for somatic mutations; the default list "
          "is appropriate for a human male; note exclusion of mitochondrial genome "
          "and 'random' chromosomes");
     CommandArgument_String_OrDefault_Doc(REF_MASK, "", "vecbitvector, if "
          "provided, only make mutations at bases that are turned on");
     EndCommandArguments;

     vecbasevector ref;
     vecString refnames;
     vecbitvector ref_amb, ref_mask;
     if ( REF_MASK != "" ) ref_mask.ReadAll(REF_MASK);
     FetchReads( ref, refnames, REF_FASTA);
     FetchReadsAmb( ref_amb, REF_FASTA );
     size_t N = ref.sumSizes();
     vec<String> chrs;
     ParseStringSet( CHRS, chrs );
     for ( int i = 0; i < COUNT; i++ )
     {    size_t P = big_random( ) % N;
          size_t t;
          for ( t = 0; t < ref.size( ); t++ )
          {    if ( P < ref[t].size( ) ) break;
               P -= ref[t].size( );    }
          if ( !Member( chrs, refnames[t] ) || ref_amb[t][P] )
          {    --i;
               continue;    }
          if ( REF_MASK != "" && !ref_mask[t][P] )
          {    --i;
               continue;    }
          cout << t << " " << P << " "
               << as_base( ( ref[t][P] + ( randomx( ) % 3 ) + 1 ) % 4 )
               << " " << FRAC << "\n";    }    }
