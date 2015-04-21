/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// Translate binary file containging look_aligns or alignlets to a text file.
/// 
/// \file Qltb2Qlta.cc
///
/// Default for OUT is IN's prefix + .qltout or .qltoutlet

#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"
#include "paths/Alignlet.h"
#include "MainTools.h"

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN);
  CommandArgument_String_OrDefault(OUT, "");
  CommandArgument_Bool_OrDefault_Doc( ALIGNLET, False,
       "does the input file contain alignlet format?");
  EndCommandArguments;

  if (OUT.empty()){
    if ( ! ALIGNLET)
      OUT = IN.SafeBefore(".qlt") + ".qltout";
    else{
      OUT = IN.SafeBefore(".qlt") + ".qltoutlet";
      if ( ! IsRegularFile( IN + ".index" ) )
	FatalErr("Index file " + IN + ".index not found." );
    }
  }

  if ( IsRegularFile( OUT ) )
    FatalErr("Output file " + OUT + " already exists. Please specify a different name." );
  

  ofstream os( OUT.c_str() );

  if ( ALIGNLET ){
    String index_file = IN + ".index";
    if ( ! IsRegularFile( index_file ) )
    FatalErr("Index file " + index_file + " not found." );
    vec<alignlet> alignlets;
    BinaryReader::readFile( IN, &alignlets );
    vec<int> index;
    BinaryReader::readFile( index_file, &index );
    for ( size_t id = 0; id < index.size(); id++ )
      if ( index[id] >= 0 )
	alignlets.at( index[id] ).PrintReadableBrief( os, ToString(id) );
    
  }else{
    vec<look_align> aligns;
    LoadLookAlignBinary( IN, aligns );
    for ( size_t i = 0; i < aligns.size(); i++ ){
      aligns[i].PrintParseable( os );
      aligns[i].PrintReadableBrief( os );
    }
  }
  
}
