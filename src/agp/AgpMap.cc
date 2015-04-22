// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#include "agp/AgpMap.h"
#include "TokenizeString.h"


void
AgpMapBuilder::BuildFrom( const String &strAgpFile, 
                          map< String, set<AgpMapContig> > &mapChromToAgpMap )
{
    ifstream agpStrm( strAgpFile.c_str() );

    String line;
    vec<String> tokens;

    while ( agpStrm )
    {
        getline( agpStrm, line );

        Tokenize( line, tokens );
        
        const int chromosomeTokenIdx = 0;
        const int startCoordTokenIdx = 1;
        const int stopCoordTokenIdx = 2;
        const int contigTokenIdx = 5;
        const int contigOrientIdx = 8;
        const int maxNeededTokenIdx = 8;

        if ( tokens.size() < maxNeededTokenIdx+1 ) continue;
        if ( ! tokens[contigTokenIdx].Contains( "contig_", 0 ) ) continue;
        if ( ! agpStrm ) continue;

        int contigId = tokens[contigTokenIdx].After( "contig_" ).Int();
        int startCoord = tokens[startCoordTokenIdx].Int();
        int stopCoord = tokens[stopCoordTokenIdx].Int();
        bool contigIsRc = ( tokens[contigOrientIdx] == "-" );

        String strChromosome = tokens[chromosomeTokenIdx].After( "chr" );
        
        map<String, set<AgpMapContig> >::iterator mapIter;

        mapIter = mapChromToAgpMap.find( strChromosome );

        if ( mapIter == mapChromToAgpMap.end() )
        {
            mapChromToAgpMap.insert( make_pair( strChromosome, set<AgpMapContig>() ) );
            mapIter = mapChromToAgpMap.find( strChromosome );
        }

        mapIter->second.insert( AgpMapContig( contigId, startCoord, stopCoord, contigIsRc ) );
    }
}
