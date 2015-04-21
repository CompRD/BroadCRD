///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BigBuilder.  Compile successive revisions of BroadCRD codebase on directories
// dirs[j] ( j = 0, ... ) local to some machine.  Search for "hardcoded stuff" below 
// to see what these and other parameters are set to now.  
//
// Revision xxxyy is stored in dirs[0] (or linked to from dirs[0]) as xxx/yy, with
// subdirectories:
// - BroadCRD   = source code
// - bin        = executables
// - bin_debug  = debug executables.
// Within these directories, soft links are used were possible to reduce space 
// usage.
//
// Upon starting up, BigBuilder will attempt to pick up where it left off before.
// If started from scratch, it will begin from 500 revisions behind the current
// revision.
//
// If BigBuilder runs out of space, it will delete the oldest directory xxx,
// typically containing 100 revisions.  Note that under pathological circumstances,
// this might not be the right thing to do.
//
// The directory /local/crdscratch/builds on crd12 is exported read-only as
// /wga/builds.  This is a soft mount, so occasionally you may observe that the
// mount fails.  This could be associated with high load on crd12.
//
// Known problems:
// (a) does not recognize a badly broken build
//
// To do:
// (1) What to do about dead executables.
// (2) Figure out where the running copying of BigBuilder lives, and who runs it.
// (3) Start from cron file.
// (4) Better error handling.

#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "system/AllOfOutput.h"

void SymlinkForceRelative( String existing_file, String name_of_symbolic_link )
{    existing_file.GlobalReplaceBy( "/./", "/" );
     vec<char> sep;
     sep.push_back( '/' );
     vec<String> x1, x2;
     TokenizeStrictly( existing_file, sep, x1 );
     TokenizeStrictly( name_of_symbolic_link, sep, x2 );
     int nagree;
     for ( nagree = 0; nagree < Min( x1.isize( ), x2.isize( ) ); nagree++ )
          if ( x1[nagree] != x2[nagree] ) break;
     String new_existing_file;
     for ( int j = nagree; j < x1.isize( ) - 1; j++ )
     {    new_existing_file += "..";
          if ( j < x1.isize( ) - 2 ) new_existing_file += "/";    }
     for ( int j = nagree; j < x1.isize( ); j++ )
          new_existing_file += "/" + x1[j];
     SymlinkForce( new_existing_file, name_of_symbolic_link );    }

// Define core utilities.

int64_t FileSystemFreeSpaceKB( const String& fs )
{    fast_pipe_ifstream in( "df -P " + fs );
     String line, junk;
     getline( in, line ), getline( in, line );
     istringstream iline( line.c_str( ) );
     int64_t free_space_KB;
     iline >> junk >> junk >> junk >> free_space_KB;
     return free_space_KB;    }

String FileSystemFracUsed( const String& fs )
{    fast_pipe_ifstream in( "df -P " + fs );
     String line, junk;
     getline( in, line ), getline( in, line );
     istringstream iline( line.c_str( ) );
     String frac_used;
     iline >> junk >> junk >> junk >> junk >> frac_used;
     return frac_used;    }

Bool EqualFiles( const String& f1, const String& f2 )
{    return System( "cmp " + f1 + " " + f2 + " > /dev/null" ) == 0;    }

// Global variables.

String daily_logfile, email_list, last_day_emailed, svn_rep, build_dir, log_dir;
ofstream* outp;
vec<String> dirs;
int mult, logmail;

// Logger.

void Log( const String& msg )
{    static vec<String> logs;
     String day = LineOfOutput( "date +%D" );
     *outp << Date( ) << ": " << msg << endl;
     logs.push_back( Date( ) + ": " + msg );
     if ( day != last_day_emailed )
     {    {    Ofstream( out, daily_logfile );
               if ( logs.isize( ) <= logmail )
               {    for ( int i = 0; i < logs.isize( ); i++ )
                         out << logs[i] << "\n";    }
               else
               {    for ( int i = 0; i < logmail/2; i++ )
                         out << logs[i] << "\n";
                    out << "\n... (" << logs.isize( ) - logmail << " lines) ...\n\n";
                    for ( int i = logs.isize( ) - logmail/2; i < logs.isize( ); i++ )
                         out << logs[i] << "\n";    }    }
          SystemSucceed( "mail -s \"daily BigBuilder report\" " + email_list 
               + " < " + daily_logfile );
          logs.clear( );
          last_day_emailed = day;    }    }

// Update build directory.  This take as input the revision id to be started with.
// If id <= 0 then we set id = current revision + id.  If the revision id is in the
// future, wait.

void Update( int& id )
{    String line;
     while ( id <= 0 )
     {    int status = System(
               "svn info " + svn_rep + " > " + log_dir + "/info.log" );
          if ( status != 0 )
          {    Log( "svn info FAILED" ); sleep(300); continue;    }
          fast_ifstream in( log_dir + "/info.log" );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( "Revision: ", 0 ) )
               {    id = Max( 1, int( line.After( "Revision: " ).Int( ) ) + id );
                    Log( "to start at revision " + ToString(id) );
                    break;    }    }
          if ( id <= 0 )
          {    Log( "svn info FAILED" ); sleep(300); continue;    }    }
     String at_rev = "At revision ", up_to_rev = "Updated to revision ";
     if ( !IsDirectory(build_dir) )
     {    while(1)
          {    Log( "attempting svn checkout" );
               int status = System( "cd " + dirs[0] + "; svn co -r " + ToString(id) 
                    + " " + svn_rep + " > " + log_dir + "/checkout.log" );
               if ( status != 0 || !IsDirectory(build_dir) ) 
                    Log( "svn checkout FAILED" );
               else
               {    if ( id == 0 )
                    {    String svn_out =
                              LineOfOutput( "cd " + build_dir + "; svn update" );
                         if ( !svn_out.Contains( at_rev, 0 ) )
                         {    Log( "svn update FAILED 1, output = " + svn_out );
                              status = 1;
                              if ( svn_out.Contains( "svn cleanup" ) )
                              {    Log( "attempting svn cleanup" );
                                   String svn_cleanup_out 
                                        = LineOfOutput( "cd " + build_dir + "; svn cleanup" );
                                   Log( "cleanup output = " + svn_cleanup_out );    }    }
                         else id = svn_out.Between( at_rev, "." ).Int( );    }
                    if ( status == 0 ) break;    }
               sleep(300);    }    }
     else
     {    String svn_update = log_dir + "/svn_update";
          Log( "waiting for update to revision " + ToString(id) );
          int sleeptime = 1, failcount = 0;
          while(1)
          {    sleep(sleeptime);
               int status = System( "cd " + build_dir + "; svn update -r" 
                    + ToString(id) + " > " + svn_update + " 2>&1" );
               fast_ifstream in(svn_update);
               int lines = 0;
               String xline;
               while(1)
               {    getline( in, xline );
                    if ( in.fail( ) ) break;
                    line = xline;
                    lines++;    }
               Bool OK1 = ( lines >= 1 && line.Contains( at_rev, 0 ) );
               Bool OK2 = ( lines >= 1 && line.Contains( up_to_rev, 0 ) );
               Bool OK = OK1 || OK2;
               int rev = -1;
               if (OK)
               {    if (OK1) rev = line.Between( at_rev, "." ).Int( );
                    if (OK2) rev = line.Between( up_to_rev, "." ).Int( );    }
               Bool not_there = ( lines >= 1 && line.Contains("No such revision") );
               if ( status != 0 || !OK )
               {    if ( !not_there )
                    {    failcount++;
                         String svn_out = AllOfOutput1( "cat " + svn_update );
                         if ( svn_out.size( ) > 0 && svn_out.back( ) == '\n' )
                              svn_out.resize( svn_out.isize( ) - 1 );
                         Log( "svn update FAILED 2: " + svn_out );    
                         if ( svn_out.Contains( "svn cleanup" ) )
                         {    Log( "attempting svn cleanup" );
                              String svn_cleanup_out 
                                   = LineOfOutput( "cd " + build_dir + "; svn cleanup" );
                              Log( "cleanup output = " + svn_cleanup_out );    }    }
                    sleeptime = ( failcount <= 1 ? 60 : 300 );
                    continue;    }
               id = rev;
               break;    }    }    }

// Report free space in all build directories.

void ReportFreeSpace( )
{    for ( int j = 0; j < dirs.isize( ); j++ )
     {    int64_t free_space_GB = FileSystemFreeSpaceKB( dirs[j] ) / 1000000;
          Log( "free space in " + dirs[j] 
               + " = " + ToString(free_space_GB) + " GB (" 
               + FileSystemFracUsed( dirs[j] ) + " used)" );    }    }

// Extract build id from path.  A kind of flaky thing to do.

int GetId( const String& path )
{    vec<char> separators;
     separators.push_back( '/' );
     vec<String> x;
     TokenizeStrictly( path, separators, x );
     for ( int j = 0; j < x.isize( ) - 1; j++ )
     {    if ( x[j].IsInt( ) && x[j+1].IsInt( ) )
               return x[j].Int( ) * mult + x[j+1].Int( );    }
     PRINT(path);
     ForceAssert( 0 == 1 );
     return -1;    }

// Main program.

int main( )
{
     RunTime( );

     // Define hardcoded stuff.

     String machine = "crd12";           /* machine this runs on */
     const int back_builds = 1500;       /* builds to go back if clean start */
     mult = 100;                         /* builds per directory */
     const int nproc = 50;               /* parallelization level (TO REDUCE) */
     const int64_t min_free_GB = 200;    /* min free space for big directory */
     const int max_back = 500;           /* max revisions to link back through */
     email_list = "jaffe";               /* recipients for daily email */
     logmail = 40;                       /* max log entries to mail */
     dirs;                               /* directories that can be used */
     dirs.push_back( "/local/crdscratch/builds" );
     svn_rep = "https://svn.broadinstitute.org/comprd/trunk/BroadCRD";

     // Check free space.

     for ( int j = 0; j < dirs.isize( ); j++ )
     {    int64_t free_space_GB = FileSystemFreeSpaceKB( dirs[j] ) / 1000000;
          if ( free_space_GB < 10 )
          {    cout << "There is very little free space.  You may need to "
                    << "make space manually.  Giving up." << endl;
               exit(1);    }    }

     // Create directories.

     for ( int i = 0; i < dirs.isize( ); i++ )
          Mkdir777( dirs[i] );

     // Start logging.

     log_dir = dirs[0] + "/logs";
     Mkdir777(log_dir);
     String logfile = dirs[0] + "/logs/BigBuilder.log";
     daily_logfile = dirs[0] + "/logs/BigBuilder.daily.log";
     outp = new ofstream;
     ofstream& out = *outp;
     OpenOfstream( out, logfile );
     Log( "begin BigBuilder" );
     ReportFreeSpace( );

     // Determine which build directories exist.

     int min_id = 1000000000, max_id = -1;
     vec<int> build_ids;
     vec<String> build_locs;
     {    vec<String> all = AllFiles( dirs[0] );
          for ( int i = 0; i < all.isize( ); i++ )
          {    if ( !all[i].IsInt( ) ) continue;
               int n = all[i].Int( );
               String diri = dirs[0] + "/" + all[i];
               vec<String> alln = AllFiles(diri);
               for ( int j = 0; j < alln.isize( ); j++ )
               {    if ( alln[j].Contains( ".active" , -1 ) )
                    {    Log( "deleting active directory" );
                         SystemSucceed( "/bin/rm -rf " + diri + "/" + alln[j] );    }
                    if ( !alln[j].IsInt( ) ) continue;
                    int k = alln[j].Int( );
                    if ( k < 0 || k >= mult ) continue;
                    int id = ( n * mult ) + k;
                    build_ids.push_back(id);
                    build_locs.push_back( diri + "/" + alln[j] );
                    min_id = Min( id, min_id ), max_id = Max( id, max_id );    }    }
          SortSync( build_ids, build_locs );
          if ( max_id >= 0 )
          {    Log( "found builds " + ToString(min_id) 
                    + " ... " + ToString(max_id) );    }
          else Log( "found no builds" );    }

     // Get code.

     int id = ( max_id >= 1 ? max_id + 1 : -back_builds );
     build_dir = dirs[0] + "/BroadCRD";
     Update(id);

     // Loop forever.

     while(1)
     {    
          // Suppose no directory for id1 exists.  Find space for it.

          int id1 = id / mult, id2 = id % mult;
          String dirid = dirs[0] + "/" + ToString(id1);
          Log( "making directory" );

          // Delete stuff to make space if needed, then make top-level directory.

          int64_t free_space_GB = FileSystemFreeSpaceKB( dirs[0] ) / 1000000;
          int home = 0;
          if ( free_space_GB < min_free_GB )
          {    home = -1;
               while ( home < 0 )
               {    if ( build_ids.empty( ) )
                    {    Log( "completely out of space, giving up" );
                         exit(1);    }
                    int first_id = build_ids[0];
                    int first_id1 = first_id / mult;
                    int last_id = (first_id1+1) * mult - 1;
                    Log( "out of space, deleting revisions "
                         + ToString(first_id) + "-" + ToString(last_id) );
                    String d = dirs[0] + "/" + ToString(first_id1);
                    if ( IsSymbolicLink(d) )
                    {    String links_to = ReadSymbolicLink(d);
                         Remove(d);
                         SystemSucceed( "/bin/rm -rf " + links_to );    }
                    else SystemSucceed( "/bin/rm -rf " + d );
                    vec<Bool> to_delete( build_ids.size( ), False );
                    for ( int i = 0; i < build_ids.isize( ); i++ )
                         if ( build_ids[i]/mult == first_id1 ) to_delete[i] = True;
                    EraseIf(build_ids, to_delete), EraseIf(build_locs, to_delete);
                    for ( int j = 0; j < dirs.isize( ); j++ )
                    {    int64_t free_space_KB = FileSystemFreeSpaceKB( dirs[j] );
                         if ( free_space_KB / 1000000 >= min_free_GB )
                         {    home = j;
                              break;    }    }    }
               Log( "creating " + ToString(id1) + "/ in " + dirs[home] );
               Mkdir777( dirs[home] + "/" + ToString(id1) );
               if ( home > 0 ) 
               {    SymlinkForceRelative( 
                         dirs[home] + "/" + ToString(id1), dirid );    }    }
          
          // Make directory for build.

          Mkdir777( dirs[home] + "/" + ToString(id1) );
          String sid2 = ToString(id2);
          if ( sid2.size( ) == 1 ) sid2 = "0" + sid2;
          String this_build_dir = dirid + "/" + sid2;
          Mkdir777(this_build_dir); 
          Mkdir777( this_build_dir + "/BroadCRD" );

          // Copy or link source code.

          Log( "copying or linking source code for revision " + ToString(id) );
          SystemSucceed( "cd " + build_dir 
               + "; find . -name \"*\" -print > " + log_dir + "/source_finds 2>&1" );
          fast_ifstream fin( log_dir + "/source_finds" );
          while(1)
          {    String f;
               getline( fin, f );
               if ( fin.fail( ) ) break;
               if ( IsDirectory( build_dir + "/" + f ) ) continue;
               if ( f.Contains( "/.svn/" ) ) continue;
               if ( f.Contains( "./bin/", 0 ) 
                    || f.Contains( "./bin_debug/", 0 ) 
                    || f.Contains( "./obj/", 0 )
                    || f.Contains( "./obj_debug/", 0 ) ) 
               {    continue;     }
               String dirhead = f.RevBefore( "/" );
               Mkpath( this_build_dir + "/BroadCRD/" + dirhead );
               Bool linked = False;
               if ( build_ids.nonempty( ) )
               {    int idx = build_ids.back( );
                    String bf = build_locs.back( ) + "/BroadCRD/" + f;
                    if ( IsSymbolicLink(bf) ) 
                    {    String links_to = RealPath(bf);
                         if ( id - GetId(links_to) <= max_back
                              && LastModified( build_dir + "/" + f )
                              == LastModified( links_to ) )
                         {    SymlinkForceRelative( links_to, 
                                   this_build_dir + "/BroadCRD/" + f );
                              linked = True;    }    }
                    else if ( IsRegularFile(bf) 
                         && LastModified( build_dir + "/" + f ) == LastModified(bf) )
                    {    SymlinkForceRelative( 
                              bf, this_build_dir + "/BroadCRD/" + f );    
                         linked = True;    }    }
               if ( !linked )
               {    String command = "cp -p " + build_dir + "/" + f + " " 
                         + this_build_dir + "/BroadCRD/" + dirhead 
                         + " > " + log_dir + "/cp_report 2>&1";    
                    int status = System(command);
                    if ( status != 0 )
                    {    Log( "FAILED to copy " + f );
                         Log( "full command = " + command );
                         Log( "giving up" );
                         exit(1);    }    }    }

          // There are two passes, one for regular compile, and one for a debug
          // compile.

          for ( int pass = 1; pass <= 2; pass++ )
          {    String bin_name = "bin", obj_name = "obj";
               if ( pass == 2 ) 
               {    bin_name += "_debug", obj_name += "_debug";    }
               String this_bin_dir = this_build_dir + "/" + bin_name + ".active";
               Mkdir777(this_bin_dir);

               // Run make.

               Log( String("starting ") + ( pass == 1 ? "regular" : "debug" ) 
                    + " make" );
               vec<int> lastmod;
               String bin_dir = build_dir + "/" + bin_name;
               vec<String> bin;
               if ( IsDirectory(bin_dir) ) bin = AllFiles(bin_dir);
               for ( int i = 0; i < bin.isize( ); i++ )
                    lastmod.push_back( LastModified( bin_dir + "/" + bin[i] ) );
               SortSync( bin, lastmod );
               int status = System( "cd " + build_dir + "; make -k -j" 
                    + ToString(nproc) + ( pass == 1 ? "" : " DEBUGGING=yes" )
                    + " > " + log_dir + "/make.log 2>&1" );
               Log( String("make ") + ( status == 0 ? "succeeded" : "FAILED" ) );
               Log( "copying or linking executables" );
               vec<String> bin2 = AllFiles(bin_dir);
               for ( int i = 0; i < bin2.isize( ); i++ )
               {    String ex = bin2[i];
                    int p = BinPosition( bin, ex );
                    Bool linked = False;
                    if ( build_ids.nonempty( ) && p >= 0 
                         && lastmod[p] == LastModified( bin_dir + "/" + ex ) )
                    {
                         // The executable 'ex' is unchanged.  Try to find it in the
                         // last build, and in that case link to it.  Otherwise
                         // copy it.

                         int idx = build_ids.back( );
                         if ( id - idx > max_back ) break;
                         String bex = build_locs.back( ) + "/" + bin_name + "/" + ex;
                         if ( IsSymbolicLink(bex) ) 
                         {    String links_to = RealPath(bex);
                              if ( id - GetId(links_to) <= max_back )
                              {    SymlinkForceRelative( 
                                        links_to, this_bin_dir + "/" + ex );
                                   linked = True;    }    }
                         if ( !linked && IsRegularFile(bex) )
                         {    SymlinkForceRelative( bex, this_bin_dir + "/" + ex );
                              linked = True;    }    }
                    if ( !linked )
                    {    SystemSucceed( "cp -p " + bin_dir + "/" + ex + " "
                              + this_bin_dir );    }    }
               Rename( this_bin_dir, this_build_dir + "/" + bin_name );    }

          // Advance to next build.

          ReportFreeSpace( );
          build_ids.push_back(id), build_locs.push_back(this_build_dir);
          Update(++id);    }    }
