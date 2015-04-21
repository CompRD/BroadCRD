///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "CmpAssemblies.  Compare two assemblies in directories DIR1 and DIR2 to "
  "find files that are different, or absent from DIR2.  Report the file or "
  "files with the earliest timestamp.  This only tests files that are in a "
  "given list of 'cmp-able' suffixes.";

// MakeDepend: dependency DiffFastb
// MakeDepend: dependency DiffVecKmerPath

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParseSet.h"

int main(int argc, char *argv[])
{
     RunTime();

     BeginCommandArguments;
     CommandDoc(DOC);
     CommandArgument_String_Doc(DIR1, "path of assembly 1 directory");
     CommandArgument_String_Doc(DIR2, "path of assembly 2 directory");
     CommandArgument_String_OrDefault_Doc(DIR_PRE, "", 
          "prefix DIR1 and DIR2 by this");
     CommandArgument_String_OrDefault_Doc(DIR_POST, "", 
          "suffix DIR1 and DIR2 by this");
     CommandArgument_String_OrDefault_Doc(SUFFIXES, 
          "{distribs,efasta,fasta,fastb,paths.k40,paths.k96,paths_rc.k40,"
          "paths_rc.k96,qualb,superb,unibases.k40,unibases.k96,unipaths.k40,"
          "unipaths.k96}", "list of suffixes to use");
     CommandArgument_Bool_OrDefault_Doc(SYM, False,
          "if True, report if files in DIR2 are not present in DIR1");
     CommandArgument_Bool_OrDefault(QUIET, False);
     EndCommandArguments;

     if ( !QUIET )
          cout << "Looking at files with suffixes in\n" << SUFFIXES << ".\n" << endl;
     if ( DIR_PRE != "" ) 
     {    DIR1 = DIR_PRE + "/" + DIR1, DIR2 = DIR_PRE + "/" + DIR2;    }
     if ( DIR_POST != "" ) 
     {    DIR1 = DIR1 + "/" + DIR_POST, DIR2 = DIR2 + "/" + DIR_POST;    }
     vec<String> suffixes;
     ParseStringSet( SUFFIXES, suffixes );
     vec< pair<int,String> > diffs;
     for ( int i = 0; i < suffixes.isize( ); i++ )
     {    String command;
          if ( !SYM )
          {    command = "cd " + DIR1 + "; find . -name \"*." + suffixes[i]
                    + "\" -print";    }
          else
          {    command = "( cd " + DIR1 + "; find . -name \"*." + suffixes[i]
                    + "\" -print; cd " + DIR2 + "; find . -name \"*." + suffixes[i]
                    + "\" -print ) | sort -u";    }
          fast_pipe_ifstream in(command);
          while(1)
          {    String path;
               getline( in, path );
               if ( in.fail( ) ) break;
               int status;
               String fn1 = DIR1 + "/" + path, fn2 = DIR2 + "/" + path;
               if ( IsSomeSortOfFile(fn2) )
               {    
                    // Special handling for fastb files.  I wish this wasn't
                    // needed.

                    if ( suffixes[i] == "fastb" 
                         || suffixes[i] == "unibases.k40"
                         || suffixes[i] == "unibases.k96" )
                    {    status = System( "DiffFastb IN1=" + fn1 + " IN2=" + fn2 
                              + " FAIL=True > /dev/null" );    }

                    // Special handling for kmerpath files.  Even worse.

                    else if ( suffixes[i] == "paths.k40" 
                         || suffixes[i] == "paths.k96" 
                         || suffixes[i] == "paths_rc.k40"
                         || suffixes[i] == "paths_rc.k96"
                         || suffixes[i] == "unipaths.k40"
                         || suffixes[i] == "unipaths.k40"
                         || suffixes[i] == "unipaths.k96" )
                    {    status = System( "DiffVecKmerPath FILE1=" + fn1 + " FILE2=" 
                              + fn2 + " FAIL=True > /dev/null" );    }

                    // General case.

                    else
                    {    status = System( "cmp " + fn1 + " " + fn2 
                              + " > /dev/null" );    }    }
               else status = 1;
               if ( status != 0 ) 
                    diffs.push( LastModified( DIR1 + "/" + path ), path );    }    }
     Sort(diffs);
     if ( diffs.empty( ) && !QUIET ) cout << "No differences found." << endl;
     else
     {    int i;
          for ( i = 1; i < diffs.isize( ); i++ )
               if ( diffs[i].first != diffs[0].first ) break;
          cout << "The following file(s) are the first to have changed:\n";
          for ( int j = 0; j < i; j++ )
               cout << diffs[j].second << endl;    
          return 1;    }    }
