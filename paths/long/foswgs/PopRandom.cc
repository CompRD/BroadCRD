#include "Basevector.h"
#include "MainTools.h"
#include "FastIfstream.h"
#include "ParallelVecUtilities.h"
#include "paths/long/fosmid/FosmidPool.h"
#include "random/Random.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PERSON);
     CommandArgument_Int_OrDefault(N, -1);
     CommandArgument_String_OrDefault(REGION, "");
     CommandArgument_String_OrDefault(FOS, "");
     CommandArgument_String(INSTANCE);
     EndCommandArguments;

     int opts = 0;
     if ( N >= 0 ) opts++;
     if ( REGION != "" ) opts++;
     if ( FOS != "" ) opts++;
     ForceAssertLe( opts, 1 );

     vec<String> person;
     ParseStringSet( PERSON, person );

     if ( FOS != "" )
     {    vec<int> fos;
          ParseIntSet( FOS, fos );
          vec<String> regions;
          vec< vec< pair<String,String> > > junctions, breaks, edits;
          ParseFosmidPoolMetainfo( regions, junctions, breaks, edits );
          for ( int j = 0; j < fos.isize( ); j++ )
              REGION += " " + ToString( regions[fos[j]] );    }

     if ( N == -1 )
     {    if ( person.solo( ) )
          {    SystemSucceed( "samtools view -b /wga/scr4/jaffe/fos_filter/" 
                    + person[0] + ".sorted.bam " + REGION 
                    + " > /wga/scr4/jaffe/fos_filter/sub." + INSTANCE + ".bam" );
                         }
          else
          {    Remove( "/wga/scr4/jaffe/fos_filter/sub." + INSTANCE + ".sam" );
               for ( int p = 0; p < person.isize( ); p++ )
               {    SystemSucceed( "samtools view -H /wga/scr4/jaffe/fos_filter/" 
                         + person[p] + ".sorted.bam " + REGION 
                         + " >> /wga/scr4/jaffe/fos_filter/sub." 
                         + INSTANCE + ".sam" );    }
               for ( int p = 0; p < person.isize( ); p++ )
               {    SystemSucceed( "samtools view /wga/scr4/jaffe/fos_filter/" 
                         + person[p] + ".sorted.bam " + REGION 
                         + " >> /wga/scr4/jaffe/fos_filter/sub." 
                         + INSTANCE + ".sam" );    }
               SystemSucceed( "samtools view -b -S "
                    "/wga/scr4/jaffe/fos_filter/sub." + INSTANCE + ".sam "
                    + " > /wga/scr4/jaffe/fos_filter/sub." 
                    + INSTANCE + ".bam" );    }
          return 0;    }

     String line;
     vec<String> header;
     cout << Date( ) << ": extracting header" << endl;
     for ( int p = 0; p < person.isize( ); p++ )
     {    fast_ifstream in( "/wga/scr4/jaffe/fos_filter/" + person[p] + ".sam" );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( !line.Contains( "@", 0 ) ) break;
               header.push_back(line);    }    }
     Ofstream( out, "/wga/scr4/jaffe/fos_filter/sub." + INSTANCE + "sam" );
     for ( int i = 0; i < header.isize( ); i++ )
          out << header[i] << "\n";
     for ( int p = 0; p < person.isize( ); p++ )
     {    vec<String> ids;
          cout << Date( ) << ": getting ids" << endl;
          {    fast_ifstream in( 
                    "/wga/scr4/jaffe/fos_filter/" + person[p] + ".sam" );
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    ids.push_back( line.Before( "\t" ) );    }    }
          cout << Date( ) << ": sorting ids" << endl;
          ParallelUniqueSort(ids);
          cout << Date( ) << ": picking ids" << endl;
          vec<String> keep;
          for ( int i = 0; i < N; i++ )
               keep.push_back( ids[ randomx( ) % ids.isize( ) ] );
          UniqueSort(keep);
          cout << Date( ) << ": finding reads" << endl;
          fast_ifstream in( "/wga/scr4/jaffe/fos_filter/" + person[p] + ".sam" );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( BinMember( keep, line.Before( "\t" ) ) )
                    out << line << endl;     }    }
     SystemSucceed( "cat /wga/scr4/jaffe/fos_filter/sub." + INSTANCE 
          + ".sam | wc --lines" );
     cout << Date( ) << ": bamming" << endl;
     SystemSucceed( "samtools view -b -S /wga/scr4/jaffe/fos_filter/sub." 
          + INSTANCE + ".sam " "> /wga/scr4/jaffe/fos_filter/sub." 
          + INSTANCE + ".bam" );
     cout << Date( ) << ": done" << endl;    }
