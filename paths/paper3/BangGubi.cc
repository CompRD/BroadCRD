#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "VecUtilities.h"

int main( )
{    
     RunTime( );

     String dir = "/wga/scr1/ALLPATHS/H.sapiens.NA12878/other_ref/bang_gubi";

     fast_ifstream in1( dir + "/BangGubi.fasta" );
     fast_ifstream in2( dir + "/agp04d.supercontigs" );
     Ofstream( out, dir + "/BangGubi.assembly.fasta" );

     vec<String> contig_names, contigs, scaffolds;

     String line, record;
     int line_count = 0;
     while(1)
     {    getline( in1, line );
          if ( in1.fail( ) ) break;
          if ( line_count++ % 10000 == 0 )
               cout << "line " << line_count << " of 1424780" << endl;
          if ( line.Contains( ">", 0 ) )
          {    contig_names.push_back( line.After( ">" ) );
               if ( record.size( ) > 0 ) contigs.push_back(record);
               record.clear( );    }
          else record.append(line);    }
     contigs.push_back(record);
     record.clear( );
     cout << Date( ) << ": sorting" << endl;
     SortSync( contig_names, contigs );
     cout << Date( ) << ": done" << endl;

     while(1)
     {    getline( in2, line );
          if ( in2.fail( ) ) break;
          if ( line.Contains( "super", 0 ) )
          {    if ( record.size( ) > 0 ) scaffolds.push_back(record);
               if ( scaffolds.size( ) % 10000 == 0 )
                    cout << scaffolds.size( ) << " scaffolds of 375768" << endl;
               record.clear( );    }
          else if ( line.Contains( "contig", 0 ) )
          {    int tig = BinPosition( contig_names, line.After( "contig " ) );
               ForceAssertGe( tig, 0 );
               record.append( contigs[tig] );    }
          else if ( line.Contains( "gap", 0 ) )
          {    int gap = line.Between( "gap ", " " ).Int( );
               record.append( String( gap, 'N' ) );    }    }

     cout << Date( ) << ": writing scaffolds" << endl;
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    out << ">" << i << "\n";
          int base_count = 0;
          for ( int j = 0; j < scaffolds[i].isize( ); j++ )
          {    if ( scaffolds[i][j] == 'n' ) continue;
               if ( scaffolds[i][j] == 'x' ) scaffolds[i][j] = 'N';
               if ( base_count > 0 && base_count % 80 == 0 ) out << "\n";
               base_count++;
               out << scaffolds[i][j];    }
          out << "\n";    }    }
