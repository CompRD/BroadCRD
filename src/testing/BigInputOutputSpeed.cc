///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BigInputOutputSpeed.  Fill memory (for 0.5 TB box) to prevent caching, write
// lots of random stuff to disk, then read it back in.
//
// This code has stuff copied from old versions of System.{cc,h}.
//
// Results:
//
// box        filer               minutes
// crd5       /local/crdscratch   39.4
// crd5       /wga/scr4           24.9
// crd13      /wga/scr4           25.7
// crd13      /local/tmp          18.5
// crd9+10+11 /wga/scr4           29.9, 28.6, 28.3

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "random/Random.h"

String ErrnoExplanation( int e )
{    String why = "with errno = " + ToString(e);
     why += " [";
     why += strerror(e);
     why += "]";
     #ifdef ENOSPC
          if ( e == ENOSPC ) why = "because no space was left in the filesystem";
     #endif
     #ifdef EDQUOT
          if ( e == EDQUOT ) why = "because a disk quota was exceeded\n"
               "(which may actually mean that a filesystem overflowed)";
     #endif
     return why;    }

int Open( const String& filename, int flags, int mode = 0664 )
{    int ntries = 0, fd;
     retry: fd = open( filename.c_str( ), flags, mode );
     ++ntries;
     if ( fd < 0 )
     {    String msg = "Attempt to open file " + filename + " using flags = "
               + ToString(flags) + " failed " + ErrnoExplanation(errno) + ".";

          // Handle case of interrupted function call.

          if ( errno == EINTR && ntries <= 25 )
          {    cout << msg << endl << "Retrying." << endl;
               if ( ntries <= 5 )
               {    cout << "retrying in 10 seconds" << endl;
                    sleep(10);    }
               else if ( ntries <= 15 )
               {    cout << "retrying in 5 minutes" << endl;
                    sleep(300);    }
               else if ( ntries <= 25 )
               {    cout << "retrying in 10 minutes" << endl;
                    sleep(600);    }
               goto retry;    }

          // Handle case where there are too many open files.

          if ( errno == EMFILE )
          {    cout << msg << endl;
               close(50); // Have to close *some* file or rest will fail!
               int pid = getpid( );
               cout << "\nHere is a list of all but one open file:\n" << endl;
               System( "lsof -p " + ToString(pid) );
               cout << "\nThis error is so bad that we have to abort." << endl;
               cout << "Bye." << endl;
               CRD::exit(1);    }

          // Bail.

          FatalErr( "Attempt to open file " << filename
               << " using flags = " << flags << " failed "
               << ErrnoExplanation(errno) << "." );    }
     return fd;    }

inline int OpenForRead( const String& filename )
{    return Open( filename, O_RDONLY );    }

inline int OpenForWrite( const String& filename, int mode = 0664)
{    return Open( filename, O_WRONLY | O_CREAT, mode );    }

void Close( int fd )
{    int status = close(fd);
     if ( status != 0 )
     {    FatalErr( "Attempt to close file descriptor " << fd
               << " failed\n" << ErrnoExplanation(errno) << "." );    }    }

const longlong max_read_or_write = 2147479552; // 2^31 - 2^12

void ReadBytes( int filedes, const void* buffer0, longlong nbytes )
{    char* buffer = (char*) buffer0;
     longlong orig_nbytes = nbytes;
     ssize_t answer;
     while( nbytes > max_read_or_write )
     {    answer = read( filedes, buffer, max_read_or_write );
          if ( answer != max_read_or_write ) goto fail;
          buffer += max_read_or_write;
          nbytes -= max_read_or_write;    }
     answer = read( filedes, buffer, nbytes );
     if ( answer >= 0 ) return;
     fail: FatalErr( "Read of " << orig_nbytes << " bytes using file descriptor "
          << filedes << " failed\n" << ErrnoExplanation(errno) << "." );    }

void WriteBytes( int filedes, const void* buffer0, longlong nbytes )
{    char* buffer = (char*) buffer0;
     longlong orig_nbytes = nbytes;
     ssize_t answer;
     while( nbytes > max_read_or_write )
     {    answer = write( filedes, buffer, max_read_or_write );
          if ( answer < 0 ) goto fail;
          buffer += max_read_or_write;
          nbytes -= max_read_or_write;    }
     answer = write( filedes, buffer, nbytes );
     if ( answer >= 0 ) return;
     fail: FatalErr( "Write of " << orig_nbytes << " bytes using file descriptor "
          << filedes << " failed\n" << ErrnoExplanation(errno) << "." );    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(XFILE);
     EndCommandArguments;

     // Define parameters.

     int64_t n1 = 1000000, n2 = 1000, n3 = 20;

     // Fill up memory.

     cout << Date( ) << ": filling up memory" << endl;
     int64_t G = 1000000000;
     int64_t c = 490;
     vec<char> buf( c * G );
     #pragma omp parallel for
     for ( int64_t i = 0; i < 1000 * c; i++ )
          memset( &buf[1000000*i], 1, 1000000 );
     cout << Date( ) << ": done filling memory" << endl;

     // Generate four random megabytes.

     unsigned int seed;
     timeval t;
     gettimeofday( &t, NULL );
     seed = t.tv_sec + t.tv_usec;
     srandomx(seed);
     vec<int> r(n1);
     for ( int i = 0; i < n1; i++ )
          r[i] = randomx( );

     // Concatenate 1000 copies of the four random megabytes.

     vec<int> x(n1*n2);
     #pragma omp parallel for
     for ( int64_t i2 = 0; i2 < n2; i2++ )
     for ( int64_t i1 = 0; i1 < n1; i1++ )
          x[ i2*n1 + i1 ] = r[i1];

     // Write 20 times, then read back in.
     
     cout << Date( ) << ": start write" << endl;
     double clock = WallClockTime( );
     int fd1 = OpenForWrite(XFILE);
     for ( int64_t i = 0; i < n3; i++ )
          WriteBytes( fd1, &x[0], n1 * n2 * sizeof(int) );
     Close(fd1);
     cout << Date( ) << ": start read" << endl;
     int fd2 = OpenForRead(XFILE);
     vec<int> y(n1*n2);
     for ( int64_t i = 0; i < n3; i++ )
          ReadBytes( fd2, &x[0], n1 * n2 * sizeof(int) );
     Close(fd2);
     cout << TimeSince(clock) << " used" << endl;

     // Clean up.

     Remove(XFILE);    }
