/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Agilent.  Highly experimental code for playing with Agilent chip data.

#include <strstream>

#include "FastIfstream.h"
#include "MainTools.h"
#include "math/Matrix.h"
#include "math/PolynomialFunction.h"

void CallPeaks( const vec< pair<double,double> >& xy,
     vec<double>& peak, vec<double>& peak_at, const double min_peak )
{    peak.clear( ), peak_at.clear( );
     for ( int j = 0; j < xy.isize( ); j++ )
     {    if ( xy[j].second < min_peak ) continue;
          int k;
          for ( k = j + 1; k < xy.isize( ); k++ )
               if ( xy[k].second < xy[k-1].second ) break;
          if ( k - 3 >= 0 && k + 1 < xy.isize( ) && xy[k-3].second < xy[k-2].second 
               && xy[k-2].second < xy[k-1].second && xy[k].second < xy[k-1].second 
               && xy[k+1].second < xy[k].second )
          {    double a, b, c;
               vec<double> x, y;
               for ( int u = k-3; u <= k+1; u++ )
               {    x.push_back( xy[u].first ), y.push_back( xy[u].second );    }
               BestQuadratic( x, y, a, b, c );
               double t = -b / ( 2.0 * a );
               double p = a*t*t + b*t + c;
               peak.push_back(p), peak_at.push_back(t);    }
          j = k;    }    }

int main( int argc, char *argv[] )
{    
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_Int(RUN);
     CommandArgument_Int(LANE);
     CommandArgument_Double_OrDefault(MIN_PEAK, 20.0);
     EndCommandArguments;

     String agilent_dir = "/wga/dev1/jaffe/Arachne/agilent";

     String run = "run" + ToString(RUN);

     String run_dir = agilent_dir + "/" + run;

     const double tmin = 20.0, tdelta = 0.05, tmax = 94.0;
     const double eps = 0.000001;
     const int N( int(round( (tmax-tmin)/tdelta + 1.0 )) );

     // Define the ladder.  Units are bp, nM.  Lane 0 is the ladder provided with
     // the instrument.  Lane 10 is an alternate ladder for run1, consisting of an 
     // RsaI digest of phiX174, plus the control fragments (first and last).  The 
     // concentrations for lane 10 are fabricated.

     int lad_lane = 0;
     vec<int> ladder;
     vec<double> lconc;
     if ( lad_lane == 0 )
     {    ladder.push_back(50);    lconc.push_back(251.5);
          ladder.push_back(100);   lconc.push_back(60.6);
          ladder.push_back(300);   lconc.push_back(20.2);
          ladder.push_back(500);   lconc.push_back(12.1);
          ladder.push_back(700);   lconc.push_back(8.7);
          ladder.push_back(1000);  lconc.push_back(6.1);
          ladder.push_back(1500);  lconc.push_back(4.0);
          ladder.push_back(2000);  lconc.push_back(3.0);
          ladder.push_back(3000);  lconc.push_back(2.0);
          ladder.push_back(5000);  lconc.push_back(1.2);
          ladder.push_back(7000);  lconc.push_back(0.9);
          ladder.push_back(10380); lconc.push_back(0.6);    }
     if ( lad_lane == 10 )
     {    ladder.push_back(50);    lconc.push_back(251.5);
          ladder.push_back(89);    lconc.push_back(10.0);
          ladder.push_back(138);   lconc.push_back(10.0);
          ladder.push_back(157);   lconc.push_back(10.0);
          ladder.push_back(197);   lconc.push_back(10.0);
          ladder.push_back(247);   lconc.push_back(10.0);
          ladder.push_back(392);   lconc.push_back(10.0);
          ladder.push_back(472);   lconc.push_back(10.0);
          ladder.push_back(525);   lconc.push_back(10.0);
          ladder.push_back(645);   lconc.push_back(10.0);
          ladder.push_back(964);   lconc.push_back(10.0);
          ladder.push_back(1560);  lconc.push_back(10.0);
          ladder.push_back(10380); lconc.push_back(0.6);    }

     // Determine number of samples, filename header, and chip type.

     vec<String> files = AllFiles(run_dir), headers;
     vec<int> sample_ids;
     for ( int i = 0; i < files.isize( ); i++ )
     {    if ( files[i].Contains( "_Sample" ) )
          {    sample_ids.push_back( 
                    files[i].Between( "_Sample", "." ).Int( ) );
               headers.push_back( files[i].Before( "Sample" ) );    }    }
     int nsamples = sample_ids.size( );
     UniqueSort(sample_ids), UniqueSort(headers);
     if ( !headers.solo( ) || sample_ids.isize( ) != nsamples
          || sample_ids[0] != 1 || sample_ids.back( ) != nsamples )
     {    cout << "The labeling of the sample files doesn't make sense.\n";
          exit(1);    }
     String head = run_dir + "/" + headers[0];
     int chip_type = headers[0].Between( "2100 expert_DNA ", "_" ).Int( );

     // Load the runs, thus creating fu = {fu} and tfu = {(t,fu)}.

     vec< vec<double> > fu(nsamples+1);
     vec< vec< pair<double,double> > > tfu(nsamples+1);
     for ( int i = 0; i <= nsamples; i++ )
     {    String filename = head;
          if ( i == 0 ) filename += "Ladder";
          else filename += "Sample" + ToString(i);
          filename += ".csv";
          fast_ifstream in(filename);
          String line;
          while(1)
          {    getline( in, line );
               ForceAssert( !in.fail( ) );
               if ( line.Contains( "Time,Value", 0 ) ) break;    }
          double T = tmin;
          while(1)
          {    getline( in, line );
               ForceAssert( line.Contains( "," ) );
               line.ReplaceBy( ",", " " );
               istrstream iline( line.c_str( ) );
               double t, y;
               iline >> t >> y;
               ForceAssert( Abs( t - T ) < eps );
               fu[i].push_back(y);
               tfu[i].push_back( make_pair( t, y ) );
               T += tdelta;
               if ( T - tmax > eps ) break;    }    }
    
     // Call peaks.  Note that this will only call "unambiguous" peaks: if two
     // peaks are sufficiently close, no peak at all will be called.

     vec< vec<double> > peak(nsamples+1), peak_at(nsamples+1);
     for ( int i = 0; i <= nsamples; i++ )
     {    CallPeaks( tfu[i], peak[i], peak_at[i], MIN_PEAK );
          ForceAssertGe( peak[i].size( ), 2u );    }

     // Do linear transformation on time axis to get the sample lanes to align to
     // the ladder, yielding stfu.

     vec< vec< pair<double,double> > > stfu(tfu);
     for ( int i = 1; i <= nsamples; i++ )
     {    double r1 = peak_at[i].front( ), r2 = peak_at[i].back( );
          double s1 = peak_at[lad_lane].front( ), s2 = peak_at[lad_lane].back( );
          double a = ( s1 - s2 ) / ( r1 - r2 );
          double b = s1 - a * r1;
          for ( int j = 0; j < N; j++ )
               stfu[i][j].first = a * tfu[i][j].first + b;    }

     // Define a function that maps ln(bp) to time.

     vec<double> logladder( ladder.size( ) ), time( ladder.size( ) );
     for ( int i = 0; i < ladder.isize( ); i++ )
     {    logladder[i] = log( ladder[i] );
          time[i] = peak_at[lad_lane][i];    }
     piecewise_polynomial_function mobility_inv;
     NaturalCubicSpline( logladder, time, mobility_inv );

     // Generate lnbpfu = ( ln(bp), fu ).

     vec< vec< pair<double,double> > > lnbpfu(nsamples+1), stfu0(nsamples+1);
     for ( int i = 1; i <= nsamples; i++ )
     {    for ( int j = 0; j < N; j++ )
          {    double x = stfu[i][j].first;
               if ( x >= mobility_inv( logladder.front( ) )
                    && x <= mobility_inv( logladder.back( ) ) )
               {    lnbpfu[i].push_back( make_pair( 
                         mobility_inv.Inverse( stfu[i][j].first ), 
                              stfu[i][j].second ) );    
                    stfu0[i].push_back( stfu[i][j] );    }    }    }

     // Compute bpfu, and interpolate.

     vec< vec< pair<double,double> > > bpfu(nsamples+1);
     vec<piecewise_polynomial_function> Fbpfu(nsamples+1), Ftfu(nsamples+1);
     for ( int i = 1; i <= nsamples; i++ )
     {    vec<double> x, y, tx;
          for ( int j = 0; j < lnbpfu[i].isize( ); j++ )
          {    bpfu[i].push( exp( lnbpfu[i][j].first ), lnbpfu[i][j].second );
               x.push_back( bpfu[i].back( ).first );
               y.push_back( bpfu[i].back( ).second );
               tx.push_back( stfu0[i][j].first );    }
          NaturalCubicSpline( x, y, Fbpfu[i] );
          NaturalCubicSpline( tx, y, Ftfu[i] );    }

     // Call peaks.  For each peak, find the valleys on both sides,
     // then estimate the area between them.  Valley finding is done badly.

     vec< vec<double> > pk(nsamples+1), pkat(nsamples+1), pkarea(nsamples+1), 
          tpk(nsamples+1), tpkat(nsamples+1);
     for ( int i = 0; i <= nsamples; i++ )
     {    CallPeaks( bpfu[i], pk[i], pkat[i], MIN_PEAK );
          CallPeaks( stfu0[i], tpk[i], tpkat[i], MIN_PEAK );
          ForceAssertEq( pkat[i].size( ), tpkat[i].size( ) );
          const piecewise_polynomial_function& F = Ftfu[i];
          for ( int j = 0; j < pk[i].isize( ); j++ )
          {    double t = tpkat[i][j];
               double v1 = t - 0.1, v2 = t + 0.1;
               double delta = 0.01;
               while( v1 >= stfu0[i].front( ).first - delta && F(v1-delta) < F(v1) )
               {    v1 -= delta;    }
               while( v2 <= stfu0[i].back( ).first - delta && F(v2+delta) < F(v2) 
                    && F(v2) >= 0 )
               {    v2 += delta;    }
               pkarea[i].push_back( F.Integrate( v1, v2 ) );    }    }

     // Report results for one lane.

     ForceAssertLe( LANE, nsamples );
     /*
     for ( int i = LANE; i <= LANE; i++ )
     {    for ( int j = 0; j < lnbpfu[i].isize( ); j++ )
          {    cout << bpfu[i][j].first << " " 
                    << lnbpfu[i][j].second << "\n";    }    }    
     for ( int i = LANE; i <= LANE; i++ )
     {    for ( int j = 0; j < tfu[i].isize( ); j++ )
          {    cout << tfu[i][j].first << " " 
                    << tfu[i][j].second << "\n";    }    }    
     */
     cout << "\nPeaks:\n";

     for ( int i = 0; i < pk[LANE].isize( ); i++ )
     {    double bp = pkat[LANE][i];
          cout << bp << " " << pkarea[LANE][i] / bp << "\n";    }    }

     /*
     vec<double> bp;
     bp.push_back(50);
     for ( int i = 0; i < pk[LANE].isize( ); i++ )
     {    cout << pkat[LANE][i] << " " 
               << 100.0 * pk[LANE][i] / pkat[LANE][i] << "\n";
          bp.push_back( pkat[LANE][i] );    }
     bp.push_back(10380);

     const vec< pair<double,double> >& xy = tfu[LANE];
     cout << "\n";
     int ipeak = 0;
     for ( int j = 0; j < xy.isize( ); j++ )
     {    if ( xy[j].second < MIN_PEAK ) continue;
          int k;
          for ( k = j + 1; k < xy.isize( ); k++ )
               if ( xy[k].second < xy[k-1].second ) break;
          if ( k - 3 >= 0 && k + 1 < xy.isize( ) && xy[k-3].second < xy[k-2].second 
               && xy[k-2].second < xy[k-1].second && xy[k].second < xy[k-1].second 
               && xy[k+1].second < xy[k].second )
          {    double a, b, c;
               vec<double> x, y;
               for ( int u = k-3; u <= k+1; u++ )
               {    x.push_back( xy[u].first ), y.push_back( xy[u].second );    }
               BestQuadratic( x, y, a, b, c );
               double t = -b / ( 2.0 * a );
               double p = a*t*t + b*t + c;
               cout << "peak at " << t << ", ";
               double area = pow( b*b - 4*a*c, 1.5 ) / (6.0*a*a);
               cout << "bp = " << bp[ipeak] << ", ";
               cout << "molarity = " << area / bp[ipeak] << "\n";
               ++ipeak;
               }
          j = k;    }

     }
     */
