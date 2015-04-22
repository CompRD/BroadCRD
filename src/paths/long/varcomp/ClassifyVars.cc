///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "FastIfstream.h"
#include "MainTools.h"

int main( )
{
     String dir = "/wga/scr4/vartests/CompareVars";
     String line;
     
     vec< vec<int> > irows;
     for ( int pass = 1; pass <= 16; pass++ )
     {    fast_pipe_ifstream in( "cat " + dir + "/*/variants.all" );
          int c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, c6 = 0, c7 = 0, c8 = 0, c9 = 0;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( !line.Contains( "Fosmid" ) ) continue;
               Bool D = line.Contains( "DISCOVAR=het" ) 
                    || line.Contains( "DISCOVAR=hom" );
               Bool G1 = line.Contains( "GATK-100=het" ) 
                    || line.Contains( "GATK-100=hom" );
               Bool G2 = line.Contains( "GATK-250=het" ) 
                    || line.Contains( "GATK-250=hom" );
               Bool C = line.Contains( "CORTEX=het" ) 
                    || line.Contains( "CORTEX=hom" );

               if ( pass == 1 && !( D && G1 && G2 && C ) ) continue;
               if ( pass == 2 && !( !D && !G1 && !G2 && !C ) ) continue;
               if ( pass == 3 && !( D && !G1 && !G2 && !C ) ) continue;
               if ( pass == 4 && !( !D && G1 && !G2 && !C ) ) continue;
               if ( pass == 5 && !( !D && !G1 && G2 && !C ) ) continue;
               if ( pass == 6 && !( !D && !G1 && !G2 && C ) ) continue;
               if ( pass == 7 && !( D && G1 && !G2 && !C ) ) continue;
               if ( pass == 8 && !( D && !G1 && G2 && !C ) ) continue;
               if ( pass == 9 && !( D && !G1 && !G2 && C ) ) continue;
               if ( pass == 10 && !( !D && G1 && G2 && !C ) ) continue;
               if ( pass == 11 && !( !D && G1 && !G2 && C ) ) continue;
               if ( pass == 12 && !( !D && !G1 && G2 && C ) ) continue;
               if ( pass == 13 && !( D && G1 && G2 && !C ) ) continue;
               if ( pass == 14 && !( D && G1 && !G2 && C ) ) continue;
               if ( pass == 15 && !( D && !G1 && G2 && C ) ) continue;
               if ( pass == 16 && !( !D && G1 && G2 && C ) ) continue;

               String cat = line.RevAfter( "," );
               cat = cat.Before( " " );
               Bool cat1 = ( cat == "c=sub" );
               Bool cat2 = ( cat == "c=ins-1" );
               Bool cat3 = ( cat == "c=ins-2-10" );
               Bool cat4 = ( cat == "c=ins-11-100" );
               Bool cat5 = ( cat == "c=ins-gt-100" );
               Bool cat6 = ( cat == "c=del-1" );
               Bool cat7 = ( cat == "c=del-2-10" );
               Bool cat8 = ( cat == "c=del-11-100" );
               Bool cat9 = ( cat == "c=del-gt-100" );

               if (cat1) c1++;
               if (cat2) c2++;
               if (cat3) c3++;
               if (cat4) c4++;
               if (cat5) c5++;
               if (cat6) c6++;
               if (cat7) c7++;
               if (cat8) c8++;
               if (cat9) c9++;    }

          vec<int> irow;
          irow.push_back( c1, c2, c3, c4, c5, c6, c7, c8, c9 );
          irows.push_back(irow);
          /*
          vec<String> row;
          row.push_back( ToString(c1), ToString(c2), ToString(c3), ToString(c4),
               ToString(c5), ToString(c6), ToString(c7), ToString(c8),
               ToString(c9) );
          rows.push_back(row);
          */
               }

     for ( int i = 0; i < irows.isize( ); i++ )
          irows[i].push_front( Sum( irows[i] ) );
     vec<int> irow;
     for ( int j = 0; j < irows[0].isize( ); j++ )
     {    int sum = 0;
          for ( int i = 0; i < irows.isize( ); i++ )
               sum += irows[i][j];
          irow.push_back(sum);    }
     irows.push_back(irow);

     vec< vec<String> > rows;
     for ( int i = 0; i < irows.isize( ); i++ )
     {    vec<String> row;
          for ( int j = 0; j < irows[i].isize( ) ; j++ )
               row.push_back( ToString( irows[i][j] ) );
          rows.push_back(row);    }
     PrintTabular( cout, rows, 2, "rrrrrrrrrr" );    }
