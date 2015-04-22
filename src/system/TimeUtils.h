/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef TIME_UTILS_H
#define TIME_UTILS_H

#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <signal.h>
#include <stdlib.h>
#include "system/Types.h"
#include "String.h"
#include <time.h>
#include <string.h>


/**
 * PrintableFileDate
 *
 * Returns file date as in "2007/3/29".
 */
inline String PrintableFileDate( const String &in_file )
{
  struct tm* clock;
  struct stat attrib;
  stat( in_file.c_str( ), &attrib );
  clock = localtime( &( attrib.st_mtime ) );

  String str_time
    = ToString( clock->tm_year + 1900 ) + "/"
    + ToString( 1 + clock->tm_mon ) + "/"
    + ToString( clock->tm_mday );

  return str_time;
}

// Logical type: intSeconds_t
// Number of seconds (for timeout alarms), as an integer.
typedef int intSeconds_t;

/**
   Class: TimeoutTimer

   Handles timeouts.  Provides a single global timer,
*/
class TimeoutTimer {
 public:

  static void SetAlarm( intSeconds_t timeoutInSeconds ) {
    timeoutOccurred_ = False;
    signal( SIGALRM, handleTimeoutSignal );
    alarm( timeoutInSeconds );
  }

  static Bool TimeoutOccurred() { return timeoutOccurred_; }

 private:
  TimeoutTimer() { }

  static void handleTimeoutSignal(int) {
    timeoutOccurred_ = True;
  }

  static volatile Bool timeoutOccurred_;

};  // class Timer

// our suupported date formats
// the following sizes are assumed, and the code will drop an
// error if they're different:
//
// year = 4 or 2, it tries to figure out which you mean
// day = 2
// month = 2
// hour = 2
// minute = 2
//
enum DATE_FORMAT {
    YEAR_MONTH_DAY_HYPHENED,
    YEAR_MONTH_DAY,
    DAY_MONTH_YEAR,
    MONTH_DAY_YEAR, // who doesn't think this format is annoying?
    YEAR_MONTH_DAY_HOUR_MINUTE, // wow this one is a big one
    // this next one looks like: year-mm-dd-hh:mn
    YEAR_MONTH_DAY_HYPHENED_HOUR_MINUTE_COLONED
};


/**
 * \class Date
 * \brief the Date class, which transforms dates between formats
 *
 * The date class provides a bunch of date transforms that let you
 * go between different date formats, and retrieve a string representing
 * the current system date.  This class is useful when interacting with
 * the file system, to decode timestamps.
 */
class DateString {
 private:
  // the raw time we currently have - all time should be
  // converted to this time
  time_t myRawtime;
  struct tm myTimeInfo;

  // a function that does the common initialization
  void init();

  // the default format type
  DATE_FORMAT defaultFormat;

 public:
  // convert the date to a specified format
  String to_s(DATE_FORMAT format) const;

  // default constructor, which sets it to the current date
  DateString();

  // the constructor given the month, day, year
  DateString(int year, int day, int month);

  // the constructor given a file pointer
  DateString(const String &date, DATE_FORMAT format);

  // take a time_t object in the constructor
  DateString(time_t time);

  // sets the date
  bool SetDateString(const String &date, DATE_FORMAT format);

  // the dreaded copy assignment operator
  DateString operator= (const DateString &d) {
    // this is the easist way (and safest way)
    SetDateString(d.to_s(YEAR_MONTH_DAY),YEAR_MONTH_DAY);
    return *this;
  }
  // the less than operator
  bool operator< (const DateString &d) const;

  // the greater than operator
  // you would think this could simply be !<
  // (not less than), but the equals case should be false
  // for both
  bool operator> (const DateString &d) const;

  // the equals than operator
  // we can make this operator out of the other two
  // (as long as they're both false for the equals case)
  bool operator== (const DateString &d) const  {
    if (!(*this < d) &&
	!(*this > d)) {
      return true;
    }
    return false;
  }
  // a way to serialize this data out
  friend ostream& operator<< (ostream& os, DateString & obj );

  // a way to serialize data in
  friend istream& operator>> (istream& os, DateString & obj );

};



#endif
