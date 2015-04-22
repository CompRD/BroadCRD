/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "TaskTimer.h"
#include "system/System.h"
#include "system/file/FileReader.h"
#include "system/file/FileWriter.h"

/**
 * TestUnitHD
 *
 * Run a simple filesystem test. All temporary files will be deleted
 * at the end of the process, and output will be sent to a log file in
 * OUT_DIR.
 * 
 * How the test works:
 *
 *  (1) - start timer
 *      - write N_FILES (each of size SIZE_GB / N_FILES, using the given
 *        BUFFER_SIZE as buffer size)
 *      - stop timer and report
 *
 *  (2) - create and write a "garbage" SIZE_GB file, to confuse the OS
 *        caching system
 *
 *  (3) - start timer
 *      - read the N_FILES
 *      - stop timer and report
 *
 *  (4) - start timer
 *      - read the N_FILES (this time faster, because of the OS' caching)
 *      - stop timer and report
 *
 *  (5) - start timer
 *      - delete N_FILES
 *      - stop timer and report
 *
 * OUT_DIR: where tests are run
 * SIZE_GB : total size of files (in Gb)
 * GARBAGE_SIZE_GB : size of garbage file (in Gb)
 * N_FILES: how many files
 * BUFFER_SIZE: size of buffer (in bytes) when writing/reading 
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( OUT_DIR );
  CommandArgument_UnsignedInt_OrDefault( SIZE_GB, 40 );
  CommandArgument_UnsignedInt_OrDefault( GARBAGE_SIZE_GB, 50 );
  CommandArgument_UnsignedInt_OrDefault( N_FILES, 1000 );
  CommandArgument_UnsignedInt_OrDefault( BUFFER_SIZE, 8192 );
  EndCommandArguments;
  
  // Dir and file names.
  String tmp_dir = OUT_DIR + "/tmp";

  String log_file = OUT_DIR + "/TestUnitHD.log";

  Mkdir777( OUT_DIR );
  Mkdir777( tmp_dir );

  // Useful const.
  const String node_name = StringOfOutput( "uname -n" );

  String fs_comm = "df -T " + OUT_DIR + " | grep -v Filesystem";
  const String fs_name = StringOfOutput( fs_comm );
  
  const longlong giga = 1000000000;
  const longlong mega = 1000000;
  const longlong total_size = giga * SIZE_GB;
  const longlong garbage_size = giga * GARBAGE_SIZE_GB;
  const longlong file_size = total_size / N_FILES;
  
  // Check initial data.
  if ( total_size % N_FILES != 0 ) {
    cout << "Fatal error: SIZE_GB Gb is not divisible by N_FILES.\n" << endl;
    return 0;
  }

  if ( BUFFER_SIZE % 4 != 0 ) {
    cout << "Fatal error: BUFFER_SIZE must be a mulitple of 4.\n" << endl;
    return 0;
  }

  // Buffers.
  int n_corebuffers = file_size / BUFFER_SIZE;
  int rem_size = file_size % BUFFER_SIZE;
  char *corebuffer = new char[ BUFFER_SIZE ];
  char *rembuffer = new char[ rem_size ];
  
  // Log file.
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );

  cout << "Sending log onto " << log_file << "\n" << endl;

  log << "Testing " << fs_name << " from " << node_name << "\n"
      << "\n"
      << "number of files: " << N_FILES << "\n"
      << "size of each file: " << file_size / mega << " Mb\n"
      << "total size: " << SIZE_GB << " Gb\n"
      << "size of garbage file: " << GARBAGE_SIZE_GB << " Gb\n"
      << "buffer size: " << BUFFER_SIZE << "\n"
      << endl;
  
  // The timer.
  TaskTimer timer;

  // Create N_FILES.
  log << Date( ) << ": writing... " << flush;

  timer.Start( );
  for (int ii=0; ii<(int)N_FILES; ii++) {
    String outfile = tmp_dir + "/file_" + ToString( ii );
    FileWriter fw(outfile);
    for (int jj=0; jj<n_corebuffers; jj++)
      fw.write( corebuffer, BUFFER_SIZE );
    fw.write( rembuffer, rem_size );
    fw.close();
  }
  timer.Stop( );
  
  log << timer.Elapsed( ) << " seconds" << endl;
  timer.Reset( );

  // IO garbage to confuse the OS caching tool.
  log << Date( ) << ": write garbage (confuse OS caching)... " << flush;

  int garbage_bufsize = 8192;
  char *garbagebuffer = new char[garbage_bufsize];
  String garbage_file = tmp_dir + "/garbage";

  timer.Start( );
  FileWriter fw(garbage_file);
  longlong current_size = 0;
  while ( current_size < garbage_size ) {
    fw.write( garbagebuffer, garbage_bufsize );
    current_size += garbage_bufsize;
  }
  fw.close();
  timer.Stop( );

  log << timer.Elapsed( ) << " seconds" << endl;
  timer.Reset( );

  log << Date( ) << ": read garbage (confuse OS caching)... " << flush;
  
  timer.Start( );
  if ( true )
  {
      FileReader fr(garbage_file.c_str());
      current_size = 0;
      while ( current_size < garbage_size ) {
        fr.read( garbagebuffer, garbage_bufsize );
        current_size += garbage_bufsize;
      }
  }
  unlink( garbage_file.c_str( ) );
  timer.Stop( );
  
  if ( garbagebuffer ) delete ( garbagebuffer );

  log << timer.Elapsed( ) << " seconds" << endl;
  timer.Reset( );

  // Read twice.
  for (int pass=0; pass<2; pass++) {
    String str_pass = ( pass == 0 ? "first" : "second" );
    log << Date( ) << ": reading (" << str_pass << " pass)... " << flush;
    
    timer.Start( );
    for (int ii=0; ii<(int)N_FILES; ii++) {
      String infile = tmp_dir + "/file_" + ToString( ii );
      FileReader fr(infile.c_str());
      for (int jj=0; jj<n_corebuffers; jj++)
	fr.read( corebuffer, BUFFER_SIZE );
      fr.read( rembuffer, rem_size );
    }
    timer.Stop( );

    log << timer.Elapsed( ) << endl;
    timer.Reset( );
  }
  
  // Clean up.
  log << Date( ) << ": deleting files... " << flush;

  timer.Start( );
  for (int ii=0; ii<(int)N_FILES; ii++) {
    String outfile = tmp_dir + "/file_" + ToString( ii );
    unlink( outfile.c_str( ) );
  }
  timer.Stop( );
  
  log << timer.Elapsed( ) << " seconds\n" << endl;
  timer.Reset( );

  // Done.
  cout << Date( ) << ": done" << endl;
  log << Date( ) << ": done" << endl;
  if ( corebuffer ) delete corebuffer;
  if ( rembuffer ) delete rembuffer;
  log.close( );

}
