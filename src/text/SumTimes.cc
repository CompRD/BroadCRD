/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2012) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// SumTimes: read from standard input, expecting a number in column 1 and a time
// unit in column 2 (seconds, minutes, hours, or days).  Add.
//
// If F is specified, read from file or files F rather than standard input.  Faster.
//
// If called with CLASSIFY=True, find lines in standard input that look like time
// lines (xxx [seconds,minutes,hours,days] yyy), and sum them by their yyy field,
// ignoring anything that comes after "=".
//
// If instead called with CLASSIFY=n, where n is a positive integer, report only 
// the top n entries.

#include <map>
#include <sstream>

#include "FastIfstream.h"
#include "math/Functions.h"
#include "MainTools.h"

void PrintTime( double T )
{    double t;
     static String units, s;
     if ( T < 60.0 ) { t = T, units = "seconds"; }
     else if ( T < 3600.0 ) { t = T/60.0, units = "minutes"; }
     else if ( T < 86400.0 ) { t = T/3600.0, units = "hours"; }
     else { t = T/86400.0, units = "days"; }
     ostringstream out;
     out << setprecision(3) << t;
     s = out.str( );
     if ( !s.Contains( "." ) && s.isize( ) == 2 ) s += ".0";
     if ( s.Contains( "." ) && s.isize( ) == 3 ) s += "0";
     cout << s << " " << units;    }

int main( int argc, char *argv[] )
{    
     RunTime( );

     BeginCommandArgumentsAcceptEmptyArgList;
     CommandArgument_String_OrDefault(CLASSIFY, "False");
     CommandArgument_String_OrDefault(F, "");
     EndCommandArguments;

     int classify;
     if ( CLASSIFY == "False" ) classify = 0;
     else if ( CLASSIFY == "True" ) classify = 1000000000;
     else if ( CLASSIFY.IsInt( ) ) classify = CLASSIFY.Int( );
     else FatalErr( "Illegal value for CLASSIFY." );
     long double x, sum = 0.0;
     String line, time_unit;
     map<String,double> used;
     vec<String> f;
     if ( F != "" ) 
     {    f = AllFilesInSource(F);
          if ( f.empty( ) )
          {    cout << "F doesn't match any files" << endl;
               exit(1);    }    }
     for ( int i = 0; i < Max( 1, f.isize( ) ); i++ )
     {    fast_ifstream* pIN = 0;
          if ( F != "" ) pIN = new fast_ifstream(f[i]);
          while(1)
          {    if ( F == "" )
               {    getline( cin, line );
                    if ( !cin ) break;    }
               else
               {    getline( *pIN, line );
                    if ( pIN->fail( ) ) break;    }
               istrstream iline( line.c_str( ) );
               iline >> x >> time_unit;
               double T;
               if ( time_unit == "seconds" ) T = x;
               else if ( time_unit == "minutes" ) T = 60.0 * x;
               else if ( time_unit == "hours" ) T = 3600.0 * x;
               else if ( time_unit == "days" ) T = 24.0 * 3600.0 * x;
               else continue;
               sum += T;
               if ( classify > 0 )
               {    String tag = line.After(time_unit);
                    if ( tag.Contains( "=" ) ) tag = tag.Before( "=" ) + "= ...";
                    used[tag] += T;    }    }
           delete pIN; }
     if ( classify > 0 )
     {    vec< pair<double,String> > times;
          for ( map<String,double>::iterator i = used.begin( ); 
               i != used.end( ); ++i )
          {    times.push_back( make_pair( i->second, i->first ) );    }
          ReverseSort(times);
          double post = 0.0;
          for ( int i = 0; i < times.isize( ); i++ )
          {    if ( i < classify )
               {    cout << setiosflags(ios::fixed) << setprecision(1) << setw(5) 
                         << 100.0 * double( times[i].first ) / double(sum) 
                         << resetiosflags(ios::fixed) << "% = ";
                    PrintTime( times[i].first );
                    cout << times[i].second << "\n";    }
               else post += times[i].first;    }
          if ( classify < times.isize( ) ) 
          {    cout << setiosflags(ios::fixed) << setprecision(1) << setw(5) 
                    << 100.0 * double( post ) / double(sum) 
                    << resetiosflags(ios::fixed) << "% = ";
               PrintTime( post );
               cout << " OTHER" << "\n";    }
          cout << "--------------------------------------------------------\n";
          cout << "100.0% = ";
          PrintTime(sum);
          cout << " total\n\n";    }
     else 
     {    PrintTime(sum);
          cout << "\n";    }    }
