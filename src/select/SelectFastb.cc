// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

/** 
\file SelectFastb.cc Select some of the entries from a fastb file.
*/

const char *DOC =
"Select some of the entries from a fastb file.  There are three ways "
"to use this program: "
" 1. Specify MIN_LENGTH.  In that case, sequences of length >= MIN_LENGTH are selected."
" 2. Specify IDS.  In that case, numbers read from the file IDS determine the selected sequences."
" 3. Specify TOP.  In that case, approximately the TOP largest contigs are selected."
" 4. Specify N.  In that case, N sequences are selected at random. "
"If this option is given, option SEED is used to initialize the "
"random number generator so that many subsets of the same file "
"can be generated.";


#include "MainTools.h"
#include "Basevector.h"
#include "random/Shuffle.h"
#include <fstream>

int main( int argc, char *argv[] )
{
     RunTime();

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandDoc(DOC);
     CommandArgument_String_OrDefault_Doc(INPUT, "", "Filename relative to PRE.");
     CommandArgument_String_Abbr_OrDefault_Doc(LINPUT, IN, "",
					       "Filename relative to current directory; alternative to INPUT.");
     CommandArgument_String_OrDefault(MIN_LENGTH, "");
     CommandArgument_String_OrDefault_Doc(IDS, "", "A file of whitespace-separated read ids.");
     CommandArgument_String_OrDefault_Doc(OUTPUT_IDS, "/dev/null", "A file of whitespace-separated read ids, listing which reads were chosen.");
     CommandArgument_UnsignedInt_OrDefault(TOP, 0);
     CommandArgument_UnsignedInt_OrDefault(N, 0);
     CommandArgument_UnsignedInt_OrDefault(SEED, 0);

     CommandArgument_String_OrDefault(FILTER_VECTOR, "");

     EndCommandArguments;

     ForceAssert( INPUT == "" ^ LINPUT == "" );
     String input = PRE + "/" + INPUT;
     if ( LINPUT != "" ) input = LINPUT;

     int options = 0;
     if ( MIN_LENGTH != "" ) ++options;
     if ( IDS != "" ) ++options;
     if ( TOP > 0 ) ++options;
     if ( N > 0 ) ++options;
     ForceAssert( options == 1 );

     vec<int> ids;
     int min_length = -1;
     size_t top = TOP;

     vecbasevector EE;

     if ( IDS != "" )
     {    Ifstream( idin, IDS );
          while(1)
          {    int n;
               idin >> n;
               if ( !idin ) break;
               ids.push_back(n);    }
          sort( ids.begin( ), ids.end( ) );
          EE.Read( input, ids, 0 );
          for ( int i = 0; i < (int) EE.size( ); i++ )
               EE[i].Print( cout, "sequence_" + ToString( ids[i] ) );
          EXIT_MAIN_NORMALLY;    }

     if ( MIN_LENGTH != "" ) min_length = MIN_LENGTH.Int( );

     EE.ReadAll( input );

    vec<int> filter_vector(EE.size());
    vecbasevector pfReads;
    if (FILTER_VECTOR != "")
    {
        std::ifstream fvs(FILTER_VECTOR.c_str());
        for (size_t i = 0; i < EE.size(); i++)
        {
            int flag;
            if ( !(fvs >> flag) )
                FatalErr("Can't read FILTER_VECTOR file: " << FILTER_VECTOR);
            filter_vector[i] = flag;
            if (flag == 1) pfReads.push_back(EE[i]);
        }
    }
    else
    {
        for (size_t i = 0; i < EE.size(); i++)
        {
            filter_vector[i] = 1;
        }
    }


     if ( N > 0 ) 
     {
       vecbasevector& ee = (pfReads.size() > 0 ? pfReads : EE);
       std::ofstream outIDs(OUTPUT_IDS.c_str());
       ForceAssertLe(N, EE.size());
       Shuffle(ee.size(), ids, SEED); // Shuffle ids to ensure no duplicates
       ids.resize(N);
       sort(ids.begin(), ids.end()); // put them in sorted order just for fun
       for (unsigned int i=0; i<N; ++i) 
       {
	        ee[ids[i]].Print( cout, "sequence_" + ToString( ids[i] ) );
	        outIDs << ids[i] << '\n';
       }
       outIDs.close();
       ForceAssert(outIDs);
       EXIT_MAIN_NORMALLY;
      }
     
     if ( top > 0 )
     {    min_length = 0;
          if ( top <= EE.size( ) )
          {    vec<int> EE_sizes( EE.size( ) );
               for ( int i = 0; i < (int) EE.size( ); i++ )
                    EE_sizes[i] = EE[i].size( );
               sort( EE_sizes.rbegin( ), EE_sizes.rend( ) );
               min_length = EE_sizes[top-1];    
               cout << "outputting contigs of length >= " << min_length
                    << endl;    }    }

     for ( int i = 0; i < (int) EE.size( ); i++ )
     {    if ( min_length >= 0 )
          {    if ( (int) EE[i].size( ) >= min_length )
                    EE[i].Print( cout, "sequence_" + ToString(i) );    }
          else
          {    if ( BinPosition( ids, i ) >= 0 )
                    EE[i].Print( cout, "sequence_" + ToString(i) );    }    }    }
