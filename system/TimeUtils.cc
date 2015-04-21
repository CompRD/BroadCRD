/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "system/TimeUtils.h"
#include "system/System.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <istream>
volatile Bool TimeoutTimer::timeoutOccurred_ = False;



// takes a date format, and returns the string represention
// of the DateString object given that format.  The available formats
// are listed in the DATE_FORMAT enumeration.
//
String DateString::to_s(DATE_FORMAT format) const {
  // we're all out of magic, so we'll need verbose code to do this

  ostringstream year;  // our year
  ostringstream month;  // our month
  ostringstream day;  // our day
  ostringstream hour;  // our hour
  ostringstream minute;  // our minute

  // set the default widths and default fills, and
  // then put the current values into each
  year << setw(4) << setfill('0') << myTimeInfo.tm_year+1900;
  month << setw(2) << setfill('0') << myTimeInfo.tm_mon;
  day << setw(2) << setfill('0') << myTimeInfo.tm_mday;
  hour << setw(2) << setfill('0')<< myTimeInfo.tm_hour;
  minute << setw(2) << setfill('0') << myTimeInfo.tm_min;

  // now switch and put them together in the target order
  switch ( format )
    {
    case YEAR_MONTH_DAY_HYPHENED :
	return year.str() + "-" + month.str() + "-" + day.str();
    case  YEAR_MONTH_DAY :
      return year.str() + month.str() + day.str();
    case  DAY_MONTH_YEAR :
      return day.str() + month.str() + year.str();
    case  MONTH_DAY_YEAR :
      return month.str() + day.str() + year.str();
    case YEAR_MONTH_DAY_HOUR_MINUTE :
      return year.str() + month.str() + day.str() + hour.str() + minute.str();
    case YEAR_MONTH_DAY_HYPHENED_HOUR_MINUTE_COLONED :
      return year.str() + "-" + month.str() + "-" + day.str() + "-" + hour.str() + ":" + minute.str();
    default:
      FatalErr("DateString:to_s() something bad happened!");
    }

  // this line will never be reached, but it prevents warnings about
  // no return value
  return "";
}
/// default constructor, which sets it to the current date
DateString::DateString() {
  init();
}

/// the constructor given the month, day, year
DateString::DateString(int year, int day, int month) {
  init();
  myTimeInfo.tm_year = year - 1900;
  myTimeInfo.tm_mon = month - 1;
  myTimeInfo.tm_mday = day;
}

/// the constructor
DateString::DateString(const String &date, DATE_FORMAT format) {
  init();
  SetDateString(date, format);
}

/// the constructor given a time_t object
DateString::DateString(time_t time) {
  defaultFormat = YEAR_MONTH_DAY;
  memcpy(&myRawtime, &time, sizeof(time_t));
  struct tm *pt = localtime ( &myRawtime );
  memcpy(&myTimeInfo,pt,sizeof(myTimeInfo));
  ++myTimeInfo.tm_mon; // to make month 1-12 based
}


// set the date, given a string like "080601" or "2008_08_01" etc
// and format type it represents
//
bool DateString::SetDateString(const String &date, DATE_FORMAT format) {
  // have we found an appropriate conversion?
  bool foundConverstion = false;
  try {
    // we're all out of magic, so we'll need verbose code to do this
    switch ( format )
      {
      case YEAR_MONTH_DAY_HYPHENED :
	if (date.isize() == 10) {
	  foundConverstion = true;
	  myTimeInfo.tm_year = date.Before("-").Int() - 1900;
	  myTimeInfo.tm_mon = date.After("-").Before("-").Int() - 1;
	  myTimeInfo.tm_mday = date.After("-").After("-").Int();
	}
	break;
      case  YEAR_MONTH_DAY :
	if (date.isize() == 8) {
	  foundConverstion = true;
	  myTimeInfo.tm_year = date.substr(0,4).Int() - 1900;
	  myTimeInfo.tm_mon = date.substr(4,2).Int();
	  myTimeInfo.tm_mday = date.substr(6,2).Int();
	}
	else if (date.isize() == 6) {
	  foundConverstion = true;
	  myTimeInfo.tm_year = date.substr(0,2).Int();
	  // correct for small date format, ie 08 means 2008 not 1908
	  if (date.substr(0,2).Int() < 70) {myTimeInfo.tm_year += 100;}
	  myTimeInfo.tm_mon = date.substr(2,2).Int();
	  myTimeInfo.tm_mday = date.substr(4,2).Int();
	}
	break;
      case  DAY_MONTH_YEAR :
	if (date.isize() == 8) {
	  foundConverstion = true;
	  myTimeInfo.tm_year = date.substr(4,4).Int() - 1900;
	  myTimeInfo.tm_mon = date.substr(2,2).Int();
	  myTimeInfo.tm_mday = date.substr(0,2).Int();
	}
	else if (date.isize() == 6) {
	  foundConverstion = true;
	  myTimeInfo.tm_year = date.substr(4,2).Int();
	  // correct for small date format, ie 08 means 2008 not 1908
	  if (date.substr(0,2).Int() < 70) {myTimeInfo.tm_year += 100;}
	  myTimeInfo.tm_mon = date.substr(2,2).Int();
	  myTimeInfo.tm_mday = date.substr(0,2).Int();
	}
	break;
      case  MONTH_DAY_YEAR :
	if (date.isize() == 8) {
	  foundConverstion = true;
	  myTimeInfo.tm_year = date.substr(4,4).Int() - 1900;
	  myTimeInfo.tm_mon = date.substr(0,2).Int();
	  myTimeInfo.tm_mday = date.substr(2,2).Int();
	}
	else if (date.isize() == 6) {
	  foundConverstion = true;
	  myTimeInfo.tm_year = date.substr(4,2).Int();
	  // correct for small date format, ie 08 means 2008 not 1908
	  if (date.substr(0,2).Int() < 70) {myTimeInfo.tm_year += 100;}
	  myTimeInfo.tm_mon = date.substr(0,2).Int();
	  myTimeInfo.tm_mday = date.substr(2,2).Int();
	}
	break;
      case YEAR_MONTH_DAY_HOUR_MINUTE :
	if (date.isize() == 12) {
	  foundConverstion = true;
	  myTimeInfo.tm_year = date.substr(0,4).Int();
	  myTimeInfo.tm_mon = date.substr(4,2).Int();
	  myTimeInfo.tm_mday = date.substr(6,2).Int();
	  myTimeInfo.tm_hour = date.substr(8,2).Int();
	  myTimeInfo.tm_min = date.substr(8,2).Int();
	}
      case YEAR_MONTH_DAY_HYPHENED_HOUR_MINUTE_COLONED :
	if (date.isize() == 16) {
	  foundConverstion = true;
	  myTimeInfo.tm_year = date.substr(0,4).Int();
	  myTimeInfo.tm_mon = date.substr(5,2).Int();
	  myTimeInfo.tm_mday = date.substr(8,2).Int();
	  myTimeInfo.tm_hour = date.substr(11,2).Int();
	  myTimeInfo.tm_min = date.substr(15,2).Int();
	}
      }

    // make sure we found a valid conversion
    if (foundConverstion != true) {
      cout << "no conversion" << endl;
      FatalErr("No valid DateString conversion found. (DateString.cc)");
    }
  } // try
  catch (...) { // catch all possible exceptions
    cout << "Failed parsing" << endl;
  }
  return foundConverstion;
}


/// a function that does the common initialization
/// which is to simply set the variables to the current time
void DateString::init() {
  defaultFormat = YEAR_MONTH_DAY;
  time ( &myRawtime );
  struct tm *pt = localtime ( &myRawtime );
  memcpy(&myTimeInfo,pt,sizeof(myTimeInfo));
  ++myTimeInfo.tm_mon; // to make month 1-12 based
}
// the less than operator
bool DateString::operator< (const DateString &d) const {
  if (myTimeInfo.tm_year < d.myTimeInfo.tm_year) {
    return true;
  } else if (myTimeInfo.tm_year > d.myTimeInfo.tm_year) {
    return false;
  }
  if (myTimeInfo.tm_mon < d.myTimeInfo.tm_mon) {
    return true;
  } else if (myTimeInfo.tm_mon > d.myTimeInfo.tm_mon) {
    return false;
  }
  if (myTimeInfo.tm_mday < d.myTimeInfo.tm_mday) {
    return true;
  } else if (myTimeInfo.tm_mday > d.myTimeInfo.tm_mday) {
    return false;
  }
  if (myTimeInfo.tm_hour < d.myTimeInfo.tm_hour) {
    return true;
  } else if (myTimeInfo.tm_hour > d.myTimeInfo.tm_hour) {
    return false;
  }
  if (myTimeInfo.tm_min < d.myTimeInfo.tm_min) {
    return true;
  }
  // well they appear complete equal, or the min of other is greater,
  // either way return false
  return false;
}


// the greater than operator
bool DateString::operator> (const DateString &d) const {
  if (myTimeInfo.tm_year < d.myTimeInfo.tm_year) {
    return false;
  } else if (myTimeInfo.tm_year > d.myTimeInfo.tm_year) {
    return true;
  }
  if (myTimeInfo.tm_mon < d.myTimeInfo.tm_mon) {
    return false;
  } else if (myTimeInfo.tm_mon > d.myTimeInfo.tm_mon) {
    return true;
  }
  if (myTimeInfo.tm_mday < d.myTimeInfo.tm_mday) {
    return false;
  } else if (myTimeInfo.tm_mday > d.myTimeInfo.tm_mday) {
    return true;
  }
  if (myTimeInfo.tm_hour < d.myTimeInfo.tm_hour) {
    return false;
  } else if (myTimeInfo.tm_hour > d.myTimeInfo.tm_hour) {
    return true;
  }
  if (myTimeInfo.tm_min > d.myTimeInfo.tm_min) {
    return true;
  }

  // well they appear complete equal, or the minute of other is greater,
  // either way return false
  return false;
}

// a way to serialize this data out
ostream& operator<< (ostream& os, DateString & obj ) {
  os << obj.to_s(obj.defaultFormat);
  return os;
}

// a way to serialize data in
istream& operator>> (istream& os, DateString & obj ) {
  string tt;
  os >> tt;
  obj.SetDateString(tt,obj.defaultFormat);
  return os;
}
