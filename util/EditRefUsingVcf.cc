///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// EditRefUsingVcf.  Parse vcf files to find *homozygous* changes, then edit
// reference accordingly.

#include "Basevector.h"
#include "Fastavector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "TokenizeString.h"

class edit {

     public:

     edit( ) { }

     edit( const String& chrom, const int pos, const String& ref,
          const String& alt ) : chrom(chrom), pos(pos), ref(ref), alt(alt) { }

     String chrom;
     int pos;
     String ref, alt;

     friend Bool operator<( const edit& e1, const edit& e2 )
     {    if ( e1.chrom < e2.chrom ) return True;
          if ( e1.chrom > e2.chrom ) return False;
          return e1.pos < e2.pos;     }

};

class editx {

     public: 

     editx( ) { }

     editx( const int start, const int stop, const String& alt )
          : start(start), stop(stop), alt(alt) { }

     int start, stop;
     String alt;

};

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(SAMPLE);
     CommandArgument_String(VCFS);
     CommandArgument_String(REF_IN);
     CommandArgument_String(REF_OUT);
     EndCommandArguments;

     // Set up.

     vec<fastavector> ref_in;
     LoadFromFastaFile( REF_IN, ref_in );
     // for ( size_t t = 0; t < ref_in.size( ); t++ ) XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     // {    for ( size_t j = 0; j < ref_in[t].size( ); j++ ) XXXXXXXXXXXXXXXXXXXXXX
     //           ref_in[t][j] = toupper( ref_in[t][j] );    } XXXXXXXXXXXXXXXXXXXXX
     vecbasevector ref_in_dummy;
     vecString refnames0;
     FetchReads( ref_in_dummy, refnames0, REF_IN );
     vec<String> refnames;
     for ( size_t i = 0; i < refnames0.size( ); i++ )
     {    String r = refnames0[i];
          if ( r.Contains( " " ) ) r = r.Before( " " );
          if ( !r.Contains( "chr", 0 ) ) r = "chr" + r;
          refnames.push_back(r);    }
     vec<edit> edits;
     vec<String> vcfs;
     ParseStringSet( VCFS, vcfs );
     Ofstream( out, REF_OUT );

     // Parse vcf files.

     for ( int v = 0; v < vcfs.isize( ); v++ )
     {    fast_ifstream in( vcfs[v] );
          String line;
          vec<char> seps_tab, seps_comma;
          seps_tab.push_back( '\t' ), seps_comma.push_back( ',' );
          int pos_CHROM = -1, pos_POS = -1, pos_REF = -1, pos_ALT = -1; 
          int pos_SAMPLE = -1;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) )
               {    cout << "Failed to find CHROM line." << endl << "Abort.";
                    exit(1);    }
               if ( line.Contains( "#CHROM", 0 ) )
               {    line = line.After( "#" );
                    vec<String> fields;
                    TokenizeStrictly( line, seps_tab, fields );
                    for ( int j = 0; j < fields.isize( ); j++ )
                    {    if ( fields[j] == "CHROM" ) pos_CHROM = j;
                         else if ( fields[j] == "POS" ) pos_POS = j;
                         else if ( fields[j] == "REF" ) pos_REF = j;
                         else if ( fields[j] == "ALT" ) pos_ALT = j;
                         else if ( fields[j] == SAMPLE ) pos_SAMPLE = j;    }
                    break;    }    }
          ForceAssert( pos_CHROM >= 0 ); ForceAssert( pos_POS >= 0 );
          ForceAssert( pos_REF >= 0 ); ForceAssert( pos_ALT >= 0 );
          ForceAssert( pos_SAMPLE >= 0 );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               vec<String> fields, alt_fields;
               TokenizeStrictly( line, seps_tab, fields );
               String chrom = fields[ pos_CHROM ];
               if ( !chrom.Contains( "chr", 0 ) ) chrom = "chr" + chrom;
               int pos = fields[ pos_POS ].Int( ) - 1;
               String ref = fields[ pos_REF ], alt = fields[ pos_ALT ];
               String genotype = fields[ pos_SAMPLE ];
               if ( genotype.Contains( ":" ) ) genotype = genotype.Before( ":" );
               genotype.GlobalReplaceBy( "|", "/" );
               if ( genotype.Before( "/" ) != genotype.After( "/" ) ) continue;
               int id = genotype.Before( "/" ).Int( );
               if ( id == 0 ) continue;
               TokenizeStrictly( alt, seps_comma, alt_fields );
               String this_alt = alt_fields[id-1];
               edits.push( chrom, pos, ref, this_alt );    }    }
     Sort(edits);

     // Review edits.

     vec< vec<editx> > editxs( ref_in.size( ) );
     for ( int j = 0; j < edits.isize( ); j++ )
     {    String chrom = edits[j].chrom;
          int t = Position( refnames, chrom );
          if ( t < 0 )
          {    cout << "Couldn't find " << chrom << " in reference." << endl;
               cout << "Abort." << endl;
               exit(1);    }
          ForceAssertLt( t, (int) ref_in.size( ) );
          int pos = edits[j].pos;
          ForceAssert( pos >= 0 );
          String ref = edits[j].ref, alt = edits[j].alt;
          ForceAssertLe( pos + ref.isize( ), (int) ref_in[t].size( ) );
          for ( int z = 0; z < ref.isize( ); z++ )
               ForceAssertEq( toupper( ref[z] ), toupper( ref_in[t][ pos + z ] ) );
          int start = pos, stop = pos + ref.size( );
          // PRINT4( chrom, start, stop, alt );

          // Punt in the rare cases where there are conflicting edits.
          // Some of these could be resolved.

          Bool conflict = False;
          for ( int k = j + 1; k < edits.isize( ); k++ )
          {    if ( edits[k].chrom != chrom ) break;
               int pos_k = edits[k].pos;
               if ( pos_k >= stop ) break;
               cout << "Warning: may conflict with next!" << endl;
               conflict = True;
               j = k;    }
          if (conflict) continue;

          // Save edit.

          editxs[t].push( start, stop, alt );    }

     // Edit and write reference.

     for ( size_t t = 0; t < ref_in.size( ); t++ )
     {    fastavector R;
          int ei = 0;
          for ( int j = 0; j < (int) ref_in[t].size( ); j++ )
          {    Bool edited = False;
               while ( ei < editxs[t].isize( ) ) 
               {    if ( editxs[t][ei].start > j ) break;
                    if ( editxs[t][ei].start < j ) ei++;
                    else
                    {    for ( int z = 0; z < editxs[t][ei].alt.isize( ); z++ )
                              R.push_back( editxs[t][ei].alt[z] );
                         j += editxs[t][ei].stop - editxs[t][ei].start - 1;
                         edited = True;
                         break;    }    }
               if ( !edited ) R.push_back( ref_in[t][j] );    }
          R.Print( out,  refnames[t] );    }    }
