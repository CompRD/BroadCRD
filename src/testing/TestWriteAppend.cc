#include "MainTools.h"
#include "Alignment.h"
#include "Vec.h"

int main( int argc, char** argv )
{
    RunTime();

    BeginCommandArguments;
    CommandArgument_String(FILE1);
    CommandArgument_String(FILE2);
    CommandArgument_String(RESULT);
    EndCommandArguments;

    cout << "Reading..." << endl;
    READ( FILE1, vec<alignment_plus>, vec1 );
    cout << "Writing..." << endl;
    WriteAppend( RESULT, vec1 );

    cout << "Reading..." << endl;
    READ( FILE2, vec<alignment_plus>, vec2 );
    cout << "Writing..." << endl;
    WriteAppend( RESULT, vec2 );

    cout << "Reading..." << endl;
    READ( RESULT, vec<alignment_plus>, result );
    
    cout << "Copying..." << endl;
    vec<alignment_plus> concat( vec1 );
    concat.insert( concat.end(), vec2.begin(), vec2.end() );

    cout << "Sorting..." << endl;
    sort( concat.begin(), concat.end() );
    sort( result.begin(), result.end() );
    
    cout << "Comparing..." << endl;
    vec<alignment_plus> not_in_result;

    set_difference( concat.begin(), concat.end(),
                    result.begin(), result.end(),
                    back_inserter( not_in_result ) );

    cout << not_in_result.size() << " elements missing." << endl;
}
