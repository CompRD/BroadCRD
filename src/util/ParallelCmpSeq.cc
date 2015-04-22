// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// ParallelCmpSeq: thread CmpSeq over multiple processors.
//
// Usage: ParallelCmpSeq FILE1=... QUAL1=... FILE2=... QUAL2=... PROCESSORS=n
//                       NO_HEADER=...
// (followed by other arguments which are passed to CmpSeq).
// FILE1, FILE2, and PROCESSORS are required.
// 
// We require that FILE2 is in fastb format.
//
// The following arguments are always passed to CmpSeq: SILENT=True NO_HEADER=True.
//
// If BINARY_ALIGNMENTS_FILE is not specified, SUMMARY_ALIGNMENTS_BRIEF=True is
// passed to CmpSeq.
//
// ParallelCmpSeq in effect breaks up FILE2 (and QUAL2 if present) into n piles, 
// and in parallel calls CmpSeq on each of the n piles.  Then ParallelCmpSeq 
// reintegrates the output, so that it is essentially the same as what would
// have been obtained through a single call to CmpSeq.  (The mechanism for this
// uses the options START2 and COUNT2 for CmpSeq, thereby trivializing the process.)
//
// If LOG_FILE is set, ParallelCmpSeq will pass this argument to CmpSeq, 
// appending the run number to the file.

#include <map>
#include <strstream>

#include <sys/wait.h>

#include "Alignment.h"
#include "Basevector.h"
#include "MainTools.h"
#include "FastIfstream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(FILE1);
     CommandArgument_String_OrDefault(QUAL1, "");
     CommandArgument_String(FILE2);
     CommandArgument_String_OrDefault(QUAL2, "");
     CommandArgument_UnsignedInt(PROCESSORS);
     // EndCommandArguments; // commented out to allow for extra arguments

     if ( command.GetHelpOnly( ) )
     {    command.PrintArgHelp( );
          cout << "plus arguments to CmpSeq." << endl << endl;
          exit(0);    }

     ForceAssert( FILE2.Contains( ".fastb", -1 ) );
     ForceAssertGt( PROCESSORS, 0u );

     // Determine residual arguments.  Note that quoted strings will not be
     // handled correctly.

     String passed_args = "SILENT=True NO_HEADER=True";
     String logfile, binary_alignments_file;
     Bool bin_append = False;
     {    String command_string = command.TheCommand( );
          istrstream in( command_string.c_str( ) );
          String stuff;
          in >> stuff;
          while(1)
          {    in >> stuff;
               if ( stuff.Contains( "PROCESSORS=", 0 ) ) continue;
               if ( stuff.Contains( "NO_HEADER=", 0 ) ) continue;
               if ( stuff.Contains( "BINARY_ALIGNMENTS_FILE=", 0 ) )
               {    binary_alignments_file = stuff.After( "=" );
                    continue;    }
               if ( stuff.Contains( "BINARY_ALIGNMENTS_FILE_APPEND=True", 0 ) )
               {    bin_append = True;
                    continue;    }
               ForceAssert( !stuff.Contains( "START2=", 0 ) );
               ForceAssert( !stuff.Contains( "COUNT2=", 0 ) );
               ForceAssert( !stuff.Contains( "SUMMARY_ALIGNMENTS_BRIEF=", 0 ) );
               ForceAssert( !stuff.Contains( "SILENT=", 0 ) );
               if ( stuff.Contains( "LOG_FILE=", 0 ) )
               {    logfile = stuff.After( "=" );
                    continue;    }
               passed_args += " " + stuff;    
               if ( in.fail( ) ) break;    }    }
     if ( binary_alignments_file == "" ) 
          passed_args += " SUMMARY_ALIGNMENTS_BRIEF=True";

     // Determine division points for FILE2 and QUAL2.  In rare circumstances we
     // lower the number of processors.

     vec<int> start(PROCESSORS), stop(PROCESSORS);
     {    vecbasevector b2;
          b2.ReadAll(FILE2);
          ForceAssert( b2.size( ) > 0 );
          int current = 0;
          for ( int i = 0; i < (int) PROCESSORS; i++ )
          {    start[i] = current;
               size_t bases_loaded = 0;
               while(1)
               {    bases_loaded += b2[current++].size( );
                    if ( current == (int) b2.size( )
                         || ( i < (int) PROCESSORS - 1 
                              && bases_loaded 
                                   >= (b2.sumSizes() * 16) / PROCESSORS ) )
                    {    stop[i] = current;
                         break;    }    }    
               if ( current == (int) b2.size( ) ) 
               {    PROCESSORS = i + 1;
                    break;    }    }    }

     // Run CmpSeq.

     map<int, int> pid_to_index;
     vec<temp_file*> output_file( PROCESSORS );
     String CmpSeq = String( argv[0] ).Before( "ParallelCmpSeq" ) + "CmpSeq";
     for ( int i = 0; i < (int) PROCESSORS; i++ )
     {    output_file[i] = new temp_file( "/tmp/ParallelCmpSeq.log_XXXXXX" );
          String logfile_arg;
          if ( logfile != "" ) logfile_arg = " LOG_FILE=" + logfile + ToString(i+1);
          String binary_arg;
          if ( binary_alignments_file != "" )
               binary_arg = " BINARY_ALIGNMENTS_FILE=" + binary_alignments_file
                    + ToString(i+1);
          int pid = Fork( CmpSeq + " START2=" + ToString( start[i] )
               + " COUNT2=" + ToString( stop[i] - start[i] ) + " " + passed_args 
               + logfile_arg + binary_arg + " > " + *output_file[i] );
          pid_to_index.insert( make_pair( pid, i ) );    }

     while ( pid_to_index.size( ) > 0 )
     {    int pid, status;
          pid = wait( &status );
          map<int, int>::iterator found = pid_to_index.find(pid);
          if ( found != pid_to_index.end( ) ) pid_to_index.erase(found);    }

     for ( int i = 0; i < (int) PROCESSORS; i++ )
     {    CpAppend( *output_file[i], cout );
          delete output_file[i];    }

     if ( binary_alignments_file != "" )
     {    if ( !bin_append ) Remove(binary_alignments_file);
          for ( int i = 0; i < (int) PROCESSORS; i++ )
          {    String bfn = binary_alignments_file + ToString(i+1);
               READ( bfn, vec<alignment_plus>, aligns );
               WriteAppend( binary_alignments_file, aligns );
               Remove(bfn);    }    }    }
