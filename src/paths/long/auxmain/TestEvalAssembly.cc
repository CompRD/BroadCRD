///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * TestEvalAssembly.cc
 *
 *  Created on: Jan 23, 2014
 *      Author: blau
 */

// MakeDepend: private

#include "MainTools.h"

#include "paths/long/EvalAssembly.h"
#include "paths/long/SupportedHyperBasevector.h"
namespace{
    void FixInversion( const HyperBasevector& hb, vec<int>& inv2 )
    {    inv2.resize( hb.EdgeObjectCount( ), -1 );
         vec< triple<basevector,int,int> > uni2;
         vec<Bool> used;
         hb.Used(used);
         for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
         {    if ( !used[e] ) continue;
              uni2.push( hb.EdgeObject(e), 0, e );
              basevector b = hb.EdgeObject(e);
              b.ReverseComplement( );
              uni2.push( b, 1, e );    }
         Sort(uni2);
         for ( int i = 0; i < uni2.isize( ); i++ )
         {    int j;
              for ( j = i + 1; j < uni2.isize( ); j++ )
                   if ( uni2[j].first != uni2[i].first ) break;
              if ( j - i == 2 && uni2[i].second != uni2[i+1].second )
              {    inv2[ uni2[i].third ] = uni2[i+1].third;
                   inv2[ uni2[i+1].third ] = uni2[i].third;    }
              else if ( j - i == 1 ) inv2[ uni2[i].third ] = -1;
              // else PRINT(j-i);
              i = j - 1;    }    };
}

int main(int argc, char *argv[])
{
    RunTime( );
    BeginCommandArguments;
    CommandArgument_String_Doc(IN, "hbv or shbv file for assembly");
    CommandArgument_String_Doc(REF_FASTB, "fastb file containing reference sequences");
     CommandArgument_Int_OrDefault_Doc(ALIGNER_K, 501, "if set to !=0, use "
          "SWA-super with K=ALIGNER_K, note that SWA-super has default value K=501");
    CommandArgument_Int_OrDefault_Doc(VERBOSITY, 0, "verbosity level");
    CommandArgument_Bool_OrDefault_Doc(CIRCULAR, False, "Set references to be circular");
    EndCommandArguments;

    HyperBasevector hb;
    if ( IN.Contains( ".shbv", -1 ) )
    {    SupportedHyperBasevector shb;
         BinaryReader::readFile( IN, &shb );
         hb = shb;    }
    else if ( IN.Contains( ".hbv", -1 ) )
    {    BinaryReader::readFile( IN, &hb );    }
    else
    {    cout << "Illegal suffix for IN." << endl;
         Scram(1);    }

    vec<int> inv2;
    FixInversion(hb,inv2);
    vecbasevector ref_seqs(REF_FASTB);

    for(const auto& ref_seq: ref_seqs){
        ref_data ref;
        ref.G.push_back( ref_seq );
        ref.G3 = ref.G;
        ref.is_circular.resize( 1, CIRCULAR );
        CreateGlocs( ref.G, ref.LG, ref.Glocs );
        CreateGlocs( ref.G3, ref.LG, ref.G3locs );
        CreateGlocs( ref.G3plus, ref.LG, ref.G3pluslocs );
        vec<basevector> p;
        p.push_back( ref_seq );
        const int KX = 80; // no justification for this
        ref.GH.push_back( HyperBasevector( KX, p ) );
        EvalAssembly(ref,hb,inv2,cout,ALIGNER_K,VERBOSITY);
    }
}
