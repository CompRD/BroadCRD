// Copyright (c) 2004 Broad Institute / Massachusetts Institute of Technology

#include "lookup/MinimalAlign.h"

void AssemblyAligns::Read(const String & fileName)
{
     if ( !IsRegularFile(fileName) )
        FatalErr( fileName << " is not a regular file" );

     String line;
     unsigned int alignsCount = 0;
     FlatMinimalAlign align; 

     for ( unsigned int pass = 1; pass <= 2; ++pass )
     {
        if ( pass == 2 ) {
            cout << alignsCount << " alignments counted in " 
                 << fileName.c_str() << endl;  
            m_minAligns.clear();
            m_minAligns.resize( alignsCount );
            alignsCount = 0;
        }
        fast_ifstream in( fileName );
        while(1)
        {    
            getline( in, line );
            if ( in.fail( ) ) break;
            if ( !line.Contains( "QUERY", 0 ) ) continue;
            if ( pass == 2 )
            {
                align.Read(line); 
                //align.Print(cout);
                m_minAligns[alignsCount] = align;                    
            }
            ++alignsCount;    
        }  
     }
}
