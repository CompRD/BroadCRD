///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Build a human reference sequence by combining the maternal reference from
// Gerstein's group with the homozygous SNP calls from DePristo's group.

#include "Basevector.h"
#include "FastIfstream.h"
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

int main( )
{    
     RunTime( );

     String human_dir = "/wga/scr1/ALLPATHS/H.sapiens.NA12878";
     String other = human_dir + "/other_ref";
     String SNP_file = human_dir + "/variants/snps.v9.vcf";

     cout << Date( ) << ": loading genome" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     vecbasevector ref( other + "/build18.fastb" );

     vec< vec<char> > maternal(23);
     String line;
     cout << Date( ) << ": reading maternal reference" << endl; // XXXXXXXXXXXXXXXXX
     for ( int ch = 1; ch <= 23; ch++ )
     {    String chr = ( ch <= 22 ? ToString(ch) : "X" );
          fast_ifstream in( 
               other + "/mom_and_pop/chr" + chr + "_NA12878_maternal.fa" );
          getline( in, line );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               for ( int j = 0; j < line.isize( ); j++ )
                    maternal[ch-1].push_back( line[j] );    }    }

     fast_ifstream in(SNP_file);
     vec<char> seps_tab, seps_comma;
     seps_tab.push_back( '\t' ), seps_comma.push_back( ',' );
     int pos_CHROM = -1, pos_POS = -1, pos_REF = -1, pos_ALT = -1; 
     int pos_SAMPLE = -1;
     cout << Date( ) << ": reading SNPs" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     String SAMPLE = "NA12878";
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

     vec<edit> edits;
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
          String genotype = fields[ pos_SAMPLE ].Before( ":" );
          genotype.GlobalReplaceBy( "|", "/" );
          if ( genotype.Before( "/" ) != genotype.After( "/" ) ) continue;
          int id = genotype.Before( "/" ).Int( );
          if ( id == 0 ) continue;
          TokenizeStrictly( alt, seps_comma, alt_fields );
          String this_alt = alt_fields[id-1];
          edits.push( chrom, pos, ref, this_alt );    }

     cout << Date( ) << ": inserting SNPs" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     for ( int ch = 1; ch <= 23; ch++ )
     {    String chr = ( ch <= 22 ? ToString(ch) : "X" );
          cout << Date( ) << ": processing chr" << chr << endl; // XXXXXXXXXXXXXXXXX

          String map_file = other + "/mom_and_pop/chr" + chr + "_NA12878.map";
          fast_ifstream xin(map_file);
          vec< triple<int,int,int> > maplist;
          while(1)
          {    getline( xin, line );
               if ( xin.fail( ) ) break;
               int p, m, r;
               istringstream iline( line.c_str( ) );
               iline >> p >> m >> r;
               maplist.push( p-1, m-1, r-1 );    }

          vec<int> ref_to_maternal( ref[ch].size( ), -1 );
          int maternal_pos = 0, ref_pos = 0;
          for ( int i = 0; i < maplist.isize( ); i++ )
          {    int m = maplist[i].second, r = maplist[i].third;
               if ( r < 0 )
               {    if ( m >= 0 )
                    {    while ( m > maternal_pos )
                         {    
                              /*
                              PRINT3( 1, ref_pos, maternal_pos ); // XXXXXXXXXXXXXXX
                              */
                              ref_to_maternal[ref_pos++] = maternal_pos++;    }
                         maternal_pos++;    }    }
               else if ( m < 0 )
               {    while( r > ref_pos ) 
                    {    
                         /*
                         PRINT3( 2, ref_pos, maternal_pos ); // XXXXXXXXXXXXXXXXXXXX
                         */
                         ref_to_maternal[ref_pos++] = maternal_pos++;    }
                    ref_pos++;    }
               else
               {    while( r > ref_pos ) 
                    {    
                         /*
                         PRINT3( 3, ref_pos, maternal_pos ); // XXXXXXXXXXXXXXXXXXXX
                         */
                         ref_to_maternal[ref_pos++] = maternal_pos++;    }    }    }
          while( ref_pos < ref[ch].isize( ) )
          {    
               /*
               PRINT3( 4, ref_pos, maternal_pos ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               */
               ref_to_maternal[ref_pos++] = maternal_pos++;    }

          int total = 0, diff = 0; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          for ( int r = 0; r < ref[ch].isize( ); r++ ) // XXXXXXXXXXXXXXXXXXXXXXXXXX
          {    int m = ref_to_maternal[r]; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               if ( m < 0 ) continue; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               if ( maternal[ch-1][m] == 'n' || maternal[ch-1][m] == 'N' ) // XXXXXX
                    continue; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               total++; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               if ( toupper( as_base( ref[ch][r] ) ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    != toupper( maternal[ch-1][m] ) ) // XXXXXXXXXXXXXXXXXXXXXXXXXXX
               {    diff++; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    // PRINT2( r, m ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                         }    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          cout << "SNP rate = " << PERCENT_RATIO( 4, diff, total ) << endl; // XXXXX

          for ( int i = 0; i < edits.isize( ); i++ )
          {    const edit& e = edits[i];
               if ( e.chrom != "chr" + ToString(ch) ) continue;
               if ( ref_to_maternal[e.pos] >= 0 )
               {    maternal[ch-1][ ref_to_maternal[e.pos] ] 
                         = e.alt[0];    }    }    }

     Ofstream( out, human_dir + "/genome.fasta" );
     out << ">chrM\n";
     for ( int j = 0; j < ref[0].isize( ); j++ )
     {    if ( j > 0 && j % 80 == 0 ) out << "\n";
          out << as_base( ref[0][j] );    }
     out << "\n";
     for ( int ch = 1; ch <= 23; ch++ )
     {    String chr = ( ch <= 22 ? ToString(ch) : "X" );
          out << ">" << chr << "\n";
          for ( int j = 0; j < maternal[ch-1].isize( ); j++ )
          {    if ( j > 0 && j % 80 == 0 ) out << "\n";
               out << maternal[ch-1][j];    }
          out << "\n";    }    }
