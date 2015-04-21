///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Trav.  This, along with files trav* comprise a system for testing and maintaining
// results of LongProto on a set of regions.  The file trav.data defines the regions 
// and the results of their assembly.  This is to be run from a checkout of 
// BroadCRD.
//
// Required environment variables (see usage in trav and travc and travr):
// - TRAV_CRD: directory for putting some stuff
// - TRAV_LCRD: directory for putting some other stuff.
// Note that in parallel mode (the default), TRAV_LCRD needs to point to a separate
// directory on each machine.

// MakeDepend: dependency LongProto

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParseSet.h"

vec<int> running;
vec<String> mach;
String head;

void trav_signal_handler( int signal_number, siginfo_t* info, void* context )
{
     if ( signal_number == SIGCHLD && info->si_code != CLD_KILLED &&
         info->si_code != CLD_DUMPED )
     {    return;    }

     cout << "\nInterrupt received, cleaning up, please wait." << endl;
     for ( int m = 0; m < running.isize( ); m++ )
     {    if ( running[m] >= 0 )
               System( "ssh " + mach[m] + " killall -9 -q LongProto." + head );    }
     _exit(1);     }
 
int main(int argc, char *argv[])
{
     RunTime( 1, &trav_signal_handler );

     BeginCommandArgumentsAcceptEmptyArgListNoHeader;
     CommandArgument_String_OrDefault(HEAD, "aaa");
     CommandArgument_String_OrDefault_Doc(EXTRA, "",
          "extra arguments for LongProto; if first is 'control', "
          "ignore trav built-in extra arguments");
     CommandArgument_Int_OrDefault(VERBOSITY, 0);
     CommandArgument_String_OrDefault_Doc(MACH, 
          "{crd5,crd7,crd8,crd9,crd10,crd11,crd4}",
          "optional list of machines to parallelize over; copies bin/LongProto and "
          "paths/long/trav/{trav,travc,travr} for use on all machines");
     CommandArgument_String_OrDefault_Doc(VERBOSE_OUT, "", 
          "send verbose output to this file");
     CommandArgument_String_OrDefault_Doc(DATA, "paths/long/trav/trav.data",
          "trav.data file");
     CommandArgument_String_OrDefault_Doc(MASTER, "", "if specified, copy all "
          "LongProto output to this file; single-machine mode only");
     CommandArgument_String_OrDefault_Doc(TEST, "",
          "can be set to trav, travc or travr to restrict tests");
     EndCommandArguments;

     // Set up logging.

     ostream* out 
          = ( VERBOSE_OUT == "" ? &cout : new ofstream( VERBOSE_OUT.c_str( ) ) );

     // Get list of machines.

     ParseStringSet( MACH, mach );
     head = HEAD;

     // Parse trav.data.

     String line;
     vec<String> CHR, START, STOP, EXPECT, CONTROLS;
     fast_ifstream in(DATA);
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( "//" ) ) line = line.Before( "//" );
          while ( line.Contains( " ", -1 ) ) line.resize( line.isize( ) - 1 );
          if ( line.empty( ) || line[0] == '#' ) continue;
          if ( !line.Contains( " " ) )
          {    FatalErr("\nCan't parse \"" << line << "\", giving up");    }
          String chr = line.Before( " " ), start, stop, controls, expect;
          if ( chr.Contains( "r" ) )
          {    expect = line.After( " " );
               while ( expect[0] == ' ' ) expect = expect.After( " " );    }
          else
          {    start = line.Between( " ", " " );
               stop = line.After( " " ).Between( " ", " " );
               if ( !line.Contains( "and" ) )
               {    expect = line.After( " " ).After( " " ).After( " " );
                    while ( expect[0] == ' ' ) expect = expect.After( " " );    }
               else
               {    line = line.After( " " ).After( " " ).After( " " );
                    controls = line.Before( " " );
                    expect = line.After( " " );
                    while ( expect[0] == ' ' ) expect = expect.After( " " );    }    }
          Bool trav = controls.empty( ), travr = chr.Contains( "r" );
          String exname = ( travr ? "travr" : ( trav ? "trav" : "travc" ) );
          if ( TEST != "" && exname != TEST ) continue;
          CHR.push_back(chr);
          START.push_back(start);
          STOP.push_back(stop);
          CONTROLS.push_back(controls);
          EXPECT.push_back(expect);    }

     // Single machine version.

     if ( mach.empty( ) )
     {    if ( MASTER != "" ) Remove(MASTER);
          for ( int i = 0; i < CHR.isize( ); i++ )
          {    String chr = CHR[i], start = START[i], stop = STOP[i]; 
               String controls = CONTROLS[i], expect = EXPECT[i];
               Bool trav = controls.empty( ), travr = chr.Contains( "r" );
               String exname = ( travr ? "travr" : ( trav ? "trav" : "travc" ) );
               if ( chr.Contains( "r" ) ) chr = chr.After( "r" );
               if ( VERBOSITY >= 1 ) 
               {    *out << Date( ) << ": trying " << exname << " " << chr << " " 
                         << start << " " << stop 
                         << ( trav ? "" : " " + controls ) << endl;    }
               SystemSucceed( "paths/long/trav/" + exname + " " + chr + " " + start 
                    + " " + stop + ( trav ? "" : " " + controls )
                    + " " + HEAD + " LongProto " + EXTRA + " > trav.out" );
               if ( MASTER != "" )
                    CpAppend( Getenv( "TRAV_LCRD" ) + "/xxx", MASTER );
               String result_line = LineOfOutput( "cat trav.out | grep STATUS" );
               if ( !result_line.Contains( "STATUS" ) )
               {    FatalErr("Can't find STATUS, giving up.");    }
               String result = result_line.After( "STATUS: " );

               if ( result != expect )
               {    cout << exname << " " << chr << " " << start << " " 
                         << stop << " changed from \"" << expect << "\" to \"" 
                         << result << "\"" << endl;    }    }    }

     // Multiple machine version.  

     else
     {    
          // Initialize.

          running.resize( mach.size( ), -1 );
          vec<double> load( mach.size( ) ), load_sampled( mach.size( ), 0 );
          vec<double> load_one( mach.size( ) );
          vec<String> RESULT( CHR.size( ) );
          vec<Bool> STARTED( CHR.size( ), False ), PRINTED( CHR.size( ), False );
          for ( int id = 0; id < CHR.isize( ); id++ )
               Remove( "trav." + ToString(id) + ".out" );
          String executable = "~/crd/bin/LongProto." + HEAD;
          String travex = "~/crd/bin/trav." + HEAD;
          String travcex = "~/crd/bin/travc." + HEAD;
          String travrex = "~/crd/bin/travr." + HEAD;
          SystemSucceed( "cp bin/LongProto " + executable );
          SystemSucceed( "cp paths/long/trav/trav " + travex );
          SystemSucceed( "cp paths/long/trav/travc " + travcex );
          SystemSucceed( "cp paths/long/trav/travr " + travrex );

          // Loop.

          while(1)
          {    
               // Collect output.

               if ( VERBOSITY >= 4 ) *out << Date( ) << ": at top of loop" << endl;
               for ( int m = 0; m < mach.isize( ); m++ )
               {    int id = running[m];
                    if ( id < 0 ) continue;
                    if ( !IsRegularFile( "trav." + ToString(id) + ".out" ) ) continue;
                    String result_line 
                         = LineOfOutput( "cat trav." + ToString(id) + ".out" );
                    if ( result_line.Contains( "STATUS: " ) )
                    {    RESULT[id] = result_line.After( "STATUS: " );
                         running[m] = -1;    }    }

               // Print output.

               for ( int id = 0; id < CHR.isize( ); id++ )
               {    if ( PRINTED[id] ) continue;
                    if ( RESULT[id].empty( ) ) break;
                    String chr = CHR[id], start = START[id], stop = STOP[id]; 
                    String controls = CONTROLS[id], expect = EXPECT[id]; 
                    String result = RESULT[id];
                    Bool trav = controls.empty( ), travr = chr.Contains( "r" );
                    String exname = ( travr ? "travr" : ( trav ? "trav" : "travc" ) );
                    if ( chr.Contains( "r" ) ) chr = chr.After( "r" );

                    String val = " ";
                    if ( exname == "travc" )
                    {    int n1 = expect.Before( " errors" ).Int( );
                         int g1 = expect.Between( "and ", " gaps" ).Int( );
                         int n2 = result.Before( " errors" ).Int( );
                         int g2 = result.Between( "and ", " gaps" ).Int( );
                         if ( g2 < g1 ) val = "+";
                         else if ( g2 > g1 ) val = "-";
                         else if ( n2 < n1 ) val = "+";
                         else if ( n2 > n1 ) val = "-";    }
                    if ( exname == "trav" )
                    {    Bool valid1 = expect.Contains( "validates" );
                         Bool valid2 = result.Contains( "validates" );
                         Bool favor1 = expect.Contains( "favor" );
                         Bool favor2 = result.Contains( "favor" );
                         Bool extra1 = expect.Contains( "extra" );
                         Bool extra2 = result.Contains( "extra" );
                         Bool inv1 = expect.Contains( "inversion" );
                         Bool inv2 = result.Contains( "inversion" );
                         Bool gap1 = expect.Contains( "gap" );
                         Bool gap2 = result.Contains( "gap" );
                         int e1 = -1, e2 = -1, f1 = -1, f2 = -1;
                         if ( valid1 && !valid2 ) val = "-";
                         if ( !valid1 && valid2 ) val = "+";
                         if (extra1) e1 = expect.Before( " extra" ).Int( );
                         if (extra2) e2 = result.Before( " extra" ).Int( );
                         if (favor1) f1 = expect.Before( " favor" ).Int( );
                         if (favor2) f2 = result.Before( " favor" ).Int( );
                         if ( favor1 && extra2 ) val = "-";
                         if ( favor2 && extra1 ) val = "+";
                         if ( extra1 && extra2 && e1 < e2 ) val = "-";
                         if ( extra1 && extra2 && e1 > e2 ) val = "+";
                         if ( inv1 && extra2 ) val = "+";
                         if ( inv2 && extra1 ) val = "-";
                         if ( favor1 && gap2 ) val = "-";
                         if ( favor2 && gap1 ) val = "+";
                         if ( favor1 && favor2 && f1 < f2 ) val = "-";
                         if ( favor1 && favor2 && f2 < f1 ) val = "+";
                         if ( favor1 && inv2 ) val = "-";
                         if ( favor2 && inv1 ) val = "+";
                         if ( gap1 && inv2 ) val = "+";
                         if ( gap2 && inv1 ) val = "-";    }
                    if ( exname == "travr" )
                    {    Bool perf1 = expect.Contains( "perfect" );
                         Bool perf2 = result.Contains( "perfect" );
                         int n1 = -1, n2 = -1;
                         if (perf1)
                         {    n1 = expect.Before( "-" ).Int( )
                                   + expect.Between( "plus ", " other" ).Int( );    }
                         if (perf2)
                         {    n2 = result.Before( "-" ).Int( )
                                   + result.Between( "plus ", " other" ).Int( );    }
                         if ( perf1 && !perf2 ) val = "-";
                         if ( !perf1 && perf2 ) val = "+";
                         if ( perf1 && perf2 && n1 < n2 ) val = "-";
                         if ( perf1 && perf2 && n1 > n2 ) val = "+";    }

                    if ( result != expect )
                    {    cout << val << exname << " " << chr << " " << start << " " 
                              << stop << ( trav ? "" : " " + controls )
                              << " changed from \"" << expect << "\" to \"" 
                              << result << "\"" << endl;    }
                    else if ( VERBOSITY >= 3 )
                    {    *out << exname << " " << chr << " " << start << " " << stop
                              << ( trav ? "" : " " + controls )
                              << ": no change" << endl;    }
                    PRINTED[id] = True;    }
               if ( PRINTED.back( ) == True ) break;

               // Launch processes.

               if ( VERBOSITY >= 4 ) *out << Date( ) << ": launching" << endl;
               for ( int id = 0; id < CHR.isize( ); id++ )
               {    if ( STARTED[id] || RESULT[id].nonempty( ) ) continue;
                    String chr = CHR[id], start = START[id], stop = STOP[id]; 
                    String controls = CONTROLS[id];
                    Bool trav = controls.empty( ), travr = chr.Contains( "r" );
                    String exname 
                         = ( travr ? "travr" : ( trav ? "trav" : "travc" ) );
                    if ( chr.Contains( "r" ) ) chr = chr.After( "r" );
                    for ( int m = 0; m < mach.isize( ); m++ )
                    {    if ( running[m] >= 0 ) continue;

                         // Ignore hosts whose load was recently checked and found
                         // to be too high.

                         // const double load_max_sum = 2.0;
                         const double load_max_sum = 1.0;
                         const double load_max_one = 0.5;
                         const double load_interval = 30;
                         double clock = WallClockTime( );
                         if ( load_sampled[m] > 0 
                              && clock - load_sampled[m] <= load_interval
                              && ( load[m] > load_max_sum 
                                   || load_one[m] > load_max_one ) )

                         {    continue;    }

                         // Don't launch to a reserved host.

                         if ( VERBOSITY >= 4 ) 
                         {    *out << Date( ) << ": checking to see if " << mach[m]
                                   << " is reserved" << endl;    }
                         String res = "/wga/dev/crds-stats/" + mach[m] + "/reserve";
                         ifstream rin;
                         rin.open( res.c_str( ) );
                         if (rin)
                         {    getline( rin, line );
                              String user = line.After( " " );
                              int hour = LineOfOutput( "date +%H" ).Int( );
                              if ( user != "dexter" ) continue;
                              if ( hour >= 19 || hour <= 7 ) continue;     }

                         // Don't launch to an overloaded host.

                         if ( VERBOSITY >= 4 ) 
                              *out << Date( ) << ": getting load" << endl;
                         load_sampled[m] = clock;
                         load[m] = 0;
                         load_one[m] = 0;

                         /*
                         fast_pipe_ifstream in( "ssh " + mach[m] + " ps aux "
                              + "| sort -n -r -k3 | grep -v root" );
                         */

                         double lm = 0;
                         fast_pipe_ifstream tin( "ssh " + mach[m] + " top -b -n1 "
                              + "| grep -v root" );
                         if ( VERBOSITY >= 4 ) 
                              *out << Date( ) << ": stream open" << endl;
                         while(1)
                         {    getline( tin, line );
                              if ( tin.fail( ) )
                              {    cout << "Unexpected top output." << endl;
                                   // need to change behavior if this happens!
                                   cout << "Abort." << endl;
                                   Scram(1);    }
                              if ( line.Contains( "PID" ) ) break;    }
                         if ( VERBOSITY >= 4 ) 
                              *out << Date( ) << ": found PID line" << endl;
                         double lo = 0;
                         while(1)
                         {    getline( tin, line );
                              if ( tin.fail( ) ) break;
                              istringstream iline( line.c_str( ) );
                              String junk;
                              double x;
                              for ( int i = 0; i < 8; i++ )
                                   iline >> junk;
                              iline >> x;
                              lo = Max( lo, x / 100.0 );
                              lm += x / 100.0;    }
                         load[m] = lm;
                         load_one[m] = lo;
                         if ( VERBOSITY >= 4 ) 
                              *out << Date( ) << ": found load" << endl;

                         /*
                         while(1)
                         {    getline( in, line );
                              if ( in.fail( ) ) break;
                              istringstream iline( line.c_str( ) );
                              String junk;
                              double x;
                              iline >> junk >> junk >> x;
                              load[m] += x / 100.0;    }
                         */

                         if ( load[m] > load_max_sum ) 
                         {    if ( VERBOSITY >= 2 )
                              {    *out << "not using " << mach[m] 
                                        << " because load = " << load[m] 
                                        << endl;    }
                              continue;    }
                         if ( load_one[m] > load_max_one ) 
                         {    if ( VERBOSITY >= 2 )
                              {    *out << "not using " << mach[m] 
                                        << " because load_one = " << load_one[m] 
                                        << endl;    }
                              continue;    }
                         if ( VERBOSITY >= 3 )
                         {    *out << TimeSince(clock) << " used getting load on "
                                   << mach[m] << endl;    }

                         int pid = Fork( "ssh " + mach[m] 
                              + " \"cd ~/crd; nice +20 " 
                              + ( travr ? travrex : ( trav ? travex : travcex ) )
                              + " " + chr + " " + start + " " + stop 
                              + ( trav ? "" : " " + controls ) + " " + HEAD + "." 
                              + ToString(id) + " " + executable + " " + EXTRA 
                              + "\" > trav." + ToString(id) + ".out" );
                         if ( VERBOSITY >= 1 ) 
                         {    *out << Date( ) << ": trying " << exname << " " << chr
                                   << " " << start << " " << stop 
                                   << ( trav ? "" : " " + controls )
                                   << " (id=" << id << ") on " << mach[m] 
                                   << ", load = " << load[m] 
                                   << ", pid = " << pid << endl;    }
                         running[m] = id;
                         STARTED[id] = True;
                         break;    }
                    break;    }

               // Sleep.

               if ( VERBOSITY >= 4 ) *out << Date( ) << ": sleeping" << endl;
               sleep(1);    }    }    }
