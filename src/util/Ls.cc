///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Ls: print files, as if a make clean had been done, and truncate long file names.
// The point is to get the entire list to fit in a single window.

// If a single argument is supplied, cd to it first.

// Directories are shown in orange.
// Files which do not appear as SVN entries are shown in magenta.

// If woof.cc and woof.h both have SVN entries (or both don't), 
// just woof.cc+h is shown.

// #include <fstream>

#include "MainTools.h"

// The appropriateness of these numbers depends on your preferences and window 
// width:

const int MAX_FILENAME_WIDTH = 27;
const int COLUMNS = 3;
const int BETWEEN_COLUMN_SPACING = 1;

int main( int argc, char *argv[] )
{    if ( argc == 2 ) chdir( argv[1] );

     String START_DIRECTORY = START_GREEN;

     vector<String> allfiles = AllFiles("");
     int nfiles = allfiles.size( );
     
     // Remove files that would be deleted by make clean.

     /*
     temp_file tempfile( "/tmp/Ls_tmp_XXXXXXX" );
     vector<String> deletes;
     int status = System( "make -n clean > " + String(tempfile) + " 2>&1" );
     if ( status == 0 )
     {    Ifstream( tm, tempfile );
          String line, clean;
          while(1)
          {    getline(tm, line);
               clean += line; 
               if ( clean[ clean.size( ) - 1 ] == '\\' ) 
                    clean[ clean.size( ) - 1 ] = ' ';
               else break;    }
          tm.close( );
          if ( !clean.Contains( "rm -f", 0 ) )
               FatalErr( "I can't understand the clean target." );
          clean.erase( 0, 5 );
          System( "csh -c \"echo " + clean + " > " + tempfile + "\"" );
          String destroy;
          ifstream tmp( tempfile.c_str( ) );
          while(tmp)
          {    tmp >> destroy;
               if ( destroy != "" ) deletes.push_back(destroy);    }
          tmp.close( );    }
     vector<String> somefiles0(nfiles);
     nfiles = 0;
     for ( unsigned int i = 0; i < somefiles0.size( ); i++ )
     {    String s = allfiles[i];
          unsigned int j;
          for ( j = 0; j < deletes.size( ); j++ )
               if ( s == deletes[j] ) break;
          if ( j < deletes.size( ) ) continue;

          // I just don't want to see this:
          if ( s == "cxx_repository" ) continue;

          if ( s.size( ) >= 3 && s[s.size( ) - 1] == 'o' 
               && s[s.size( ) - 2] == '.' ) continue;
          somefiles0[nfiles++] = s;    }
     */
     vec<String> somefiles0(allfiles);

     // Identify files for which there is a CVS entry.

     vector<bool> entered0(nfiles, false);
     for ( int i = 0; i < nfiles; i++ )
     {    if ( IsSomeSortOfFile( ".svn/text-base/" + somefiles0[i] + ".svn-base" ) )
               entered0[i] = true;    }

     // Combine .cc and .h where possible.

     vector<String> somefiles(nfiles);
     vector<bool> entered(nfiles);

     nfiles = 0;
     for ( unsigned int i = 0; i < somefiles.size( ); i++ )
     {    String s = somefiles0[i];
          int sz = s.size( );
          if ( i < somefiles0.size( ) - 1 && sz >= 4 && s[sz-1] == 'c' &&
               s[sz-2] == 'c' && s[sz-3] == '.' )
          {    String t = somefiles0[i+1];
               int tz = t.size( );
               if ( tz == sz - 1 && t[tz-1] == 'h' && t[tz-2] == '.' )
               {    int j;
                    for ( j = 0; j < tz-2; j++ )
                         if ( s[j] != t[j] ) break;
                    if ( j == tz-2 && entered0[i] == entered0[i+1] )
                    {    somefiles[nfiles] = somefiles0[i] + "+h";
                         entered[nfiles++] = entered0[i];
                         ++i;
                         continue;    }    }    }
          somefiles[nfiles] = somefiles0[i];
          entered[nfiles++] = entered0[i];    }

     // Truncate long file names.

     int M = MAX_FILENAME_WIDTH;
     vector<String> shortfiles(nfiles);
     for ( int i = 0; i < nfiles; i++ )
     {    int n = somefiles[i].size( );
          if ( n > M )
          {    int Mm = M-4;
               String snew;
               snew.append( somefiles[i], 0, n/2 - (n-Mm)/2 );
               snew += "...";
               snew.append( somefiles[i], n/2 + (n-Mm)/2, n-1 );
               shortfiles[i] = snew;    }
          else shortfiles[i] = somefiles[i];    }

     // Columnate and print.

     int ncols = COLUMNS;
     int colwidth = MAX_FILENAME_WIDTH + BETWEEN_COLUMN_SPACING;
     int nrows = (nfiles+ncols-1)/ncols;
     String command;
     if ( nrows > 25 ) command = "less -r -F";
     else command = "cat";
     procbuf outp( command.c_str( ), ios::out );
     ostream out( &outp );
     out << START_WHITE;
     int total_width = COLUMNS * (MAX_FILENAME_WIDTH + BETWEEN_COLUMN_SPACING);
     for ( int i = 0; i < total_width; i++ )
          out << "-";
     out << END_ESCAPE << "\n";
     for ( int i = 0; i < nrows; i++ )
     {    for ( int j = 0; j < ncols; j++ )
          {    int where = j*nrows + i;
               if ( where < nfiles )
               {    int isdir = IsDirectory( somefiles[where] );
                    if (isdir) out << START_DIRECTORY;
                    else if ( !entered[where] ) out << START_MAGENTA;
                    String s = shortfiles[where];
                    out << s;
                    if (isdir) out << END_ESCAPE;
                    else if ( !entered[where] ) out << END_ESCAPE;
                    for ( unsigned int k = 0; k < colwidth-s.size( ); k++ )
                         out << " ";    }    }
          out << "\n";    }
     out << START_WHITE;
     for ( int i = 0; i < total_width; i++ )
          out << "-";
     out << END_ESCAPE << "\n";    }
