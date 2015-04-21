// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// FetchGenbankSequence: get fasta for a genbank entry.  This has some
// paths hardcoded for the Genome Center which will probably break at 
// some point. 

// If you set DIRECT=True, get directly from Genbank rather than from the
// local copy.

// Based partially on a script given to me by Seth Purcell,
// /seq/software/bin/Annotation/ImportExport/gb2fasta.pl

#include <ctype.h>
#include <csignal>
#include <stdlib.h>

#include <strstream>


#include "system/CommandLoader.h"
#include "MainTools.h"
#include "FastIfstream.h"
#include "TokenizeString.h"

const String yank_bin_dir = "/seq/util/bin/alpha";
const String yank_index_path = "/seq/blastdb/genbank/yank";
const String yank_index = "yank_index_gb";

int main( int argc, char *argv[] )
{
     RunTime();

     BeginCommandArguments;
     CommandArgument_String_OrDefault(NAME, "");
     CommandArgument_Bool_OrDefault(DIRECT, False);
     CommandArgument_Bool_OrDefault(PAGEDUMP, False);
     
     CommandArgument_String_OrDefault(BATCH, "");
     CommandArgument_UnsignedInt_OrDefault(CHUNK_SIZE, 100);
     
     //  full path to accs when done (only with DIRECT=True)
     CommandArgument_String_OrDefault(OUT_FILE, "");
     
     //  num secs to wait before trying to fetch a chunk of accessions
     //  (only with DIRECT=True)
     CommandArgument_UnsignedInt_OrDefault(WAIT_TIME,30);

     EndCommandArguments;

     if ( NAME.empty() && BATCH.empty() ||
          ! NAME.empty() && ! BATCH.empty() )
         InputErr( "Either NAME or BATCH must be specified." );

     if ( ! OUT_FILE.empty() && ! DIRECT )
       InputErr( "If OUT_FILE is used, DIRECT must be True." );

     if ( OUT_FILE.empty() && DIRECT )
       InputErr( "If DIRECT is True, you must specify OUT_FILE." );

     vec<String> names;
     if ( ! NAME.empty() )
         names.push_back( NAME );
     else
     {
         ifstream batch_stream( BATCH.c_str() );
         String name;
         while ( batch_stream >> name )
             names.push_back( name );
     }

     if ( ! DIRECT )
     {
         // From man putenv: "The string pointed to by [the argument
         // to putenv] becomes part of the environment," so really
         // weird things may happen to the environment if we were to
         // use env_command.c_str_mutable() when env_command goes out
         // of scope.
         String env_command = "YANK_INDEX_PATH=" + yank_index_path;
	 char *p_command = new char[env_command.size()+1];
	 memcpy(p_command,env_command.c_str(),env_command.size()+1);
         putenv( p_command );

         for ( vec<String>::iterator name_iter = names.begin(); 
               name_iter != names.end(); 
               ++name_iter )
         {
             String name = *name_iter;

             String load_command =
                 yank_bin_dir + "/yank -i " + yank_index + " -v -a " + name;

             fast_pipe_ifstream in(load_command);

             Bool start(False);
             String line;
             while(1)
             {    
                 getline( in, line );
                 if ( line.Contains( "not found" ) )
                 {    
                     cout << "Not found." << endl;
                     EXIT_MAIN_NORMALLY;
                 }
                 if ( in.fail( ) ) break;

                 if ( PAGEDUMP )
                 {
                     cout << line << endl;
                     continue;
                 }

                 if ( !start )
                 {   
                     if ( line.Contains( "Cannot find" ) )
                     {
                         cout << "Not found." << endl;
                         EXIT_MAIN_NORMALLY;
                     }
                     if ( !line.Contains( "LOCUS", 0 ) ) 
                         continue;
                     else 
                         start = True;
                 }
                 else // code to extract fasta from genbank entry
                 {   
                     if ( line.Contains( "VERSION", 0 ) )
                     {   
                         NAME = line.After( "VERSION" );
                         DeleteLeadingWhiteSpace(NAME);    
                     }
                     else if ( line.Contains( "ORIGIN", 0 ) )
                     {    
                         cout << ">" << NAME << "\n";
                         while(1)
                         {    
                             getline( in, line );
                             if ( in.fail( ) ) break;
                             if ( line.Contains( "//", 0 ) ) break;
                             for ( int i = 0; i < (int) line.size( ); i++ )
                                 if ( isalpha( line[i] ) ) 
                                     cout << line[i];
                             cout << "\n";
                         }
                         break;
                     }
                 }
             }
         }
     }

     else // DIRECT == True
     {   
         vec<String> tokens;
         vec<char> delimiters;
         delimiters.push_back('|');

	 String acc_file(OUT_FILE);
	 ofstream acc_strm( acc_file.c_str() );

	 CommandLoader *loader = CommandLoader::GetInstance( ); 

         vec<String>::iterator chunk_begin = names.begin();
         vec<String>::iterator chunk_end;
         for ( ; chunk_begin != names.end(); chunk_begin = chunk_end )
         {
             if ( distance( chunk_begin, names.end() ) > int(CHUNK_SIZE) )
                 chunk_end = chunk_begin + CHUNK_SIZE;
             else
                 chunk_end = names.end();
               
             String load_command( "lynx -dump -nolist -dont_wrap_pre "
                                  "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?"
                                  "extrafeat=-1&SendTo=on&txt=on" );
             
             load_command += "&dispmax=" + ToString( (int)distance( chunk_begin, chunk_end ) );

             if ( PAGEDUMP )
                 load_command += "&view=def";
             else
                 load_command += "&view=fasta";

             load_command += "&list_uids=";
             while ( chunk_begin != chunk_end )
             {
                 load_command += *chunk_begin;
                 ++chunk_begin;
                 if ( chunk_begin != chunk_end )
                     load_command += ',';
             }

	     //             fast_pipe_ifstream in(load_command);
	    
	     String tmp_file(OUT_FILE+".tmp");
	     Bool timed_out(True);  
	     while ( timed_out )
	     {
	       pid_t l_pid = loader->Run( load_command, tmp_file, false );
	       CommandLoader::WaitInfo w_info( l_pid, (int) WAIT_TIME );
	       if ( loader->Wait( w_info ) )
		 timed_out = False;
	       else
		 {
		 cout << "\nTimed out " << load_command << ".\n" << endl;
		 kill( l_pid, SIGTERM );
		 }
	     }
	     
	     ifstream in( tmp_file.c_str() );
             String line;
             while(1)
             {    
                 getline( in, line );
                 if ( line.Contains( "Cannot find" ) )
                 {    
                     cout << line << endl;
                     EXIT_MAIN_NORMALLY;
                 }
                 if ( in.fail( ) ) break;

                 if ( PAGEDUMP )
                 {
                     if ( ! line.empty() )
                         acc_strm << line << "\n";
                 }
                 else if ( line[0] == '>' )
                 {
                     TokenizeStrictly( line, delimiters, tokens );
                     acc_strm << ">" << tokens[3] << "\n";
                 }
                 else if ( ! line.empty () && ! isspace( line[0] ) )
                     acc_strm << line << "\n";
             }
	     in.close();
         }     
	 acc_strm.close();
	 System( "rm " + OUT_FILE + ".tmp" );
     }
}
