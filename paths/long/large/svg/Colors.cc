///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ORIGIN OF THIS CODE
// modified superficially from stackoverflow.com
// http://stackoverflow.com/questions/1313/
// followup-finding-an-accurate-distance-between-colors
// and that was taken from
// http://svn.int64.org/viewvc/int64/colors/colors.js
// Copyright (c) 2011, Cory Nelson (phrosty@gmail.com).

#include <math.h>

#include "CoreTools.h"
#include "paths/long/large/svg/Colors.h"

#define REF_X 95.047; // Observer= 2°, Illuminant= D65
#define CV_PI   3.1415926535897932384626433832795

void bgr2xyz( const vec<int>& BGR, vec<double>& XYZ );
void xyz2lab( const vec<double>& XYZ, vec<double>& Lab );
void lab2lch( const vec<double>& Lab, vec<double>& LCH );
// double deltaE2000( const vec<int>& bgr1, const vec<int>& bgr2 );
double deltaE2000( const vec<double>& lch1, const vec<double>& lch2 );

void bgr2xyz( const vec<double>& BGR, vec<double>& XYZ )
{    XYZ.resize(3);
     double r = (double)BGR[2] / 255.0;
     double g = (double)BGR[1] / 255.0;
     double b = (double)BGR[0] / 255.0;
     if ( r > 0.04045 ) r = pow( ( r + 0.055 ) / 1.055, 2.4 );
     else r = r / 12.92;
     if ( g > 0.04045 ) g = pow( ( g + 0.055 ) / 1.055, 2.4 );
     else g = g / 12.92;
     if ( b > 0.04045 ) b = pow( ( b + 0.055 ) / 1.055, 2.4 );
     else b = b / 12.92;
     r *= 100.0;
     g *= 100.0;
     b *= 100.0;
     XYZ[0] = r * 0.4124 + g * 0.3576 + b * 0.1805;
     XYZ[1] = r * 0.2126 + g * 0.7152 + b * 0.0722;
     XYZ[2] = r * 0.0193 + g * 0.1192 + b * 0.9505;    }

void xyz2lab( const vec<double>& XYZ, vec<double>& Lab )
{    Lab.resize(3);
     double x = XYZ[0] / REF_X;
     double y = XYZ[1] / REF_X;
     double z = XYZ[2] / REF_X;
     if ( x > 0.008856 ) x = pow( x, .3333333333 );
     else x = ( 7.787 * x ) + ( 16.0 / 116.0 );
     if ( y > 0.008856 ) y = pow( y, .3333333333 );
     else y = ( 7.787 * y ) + ( 16.0 / 116.0 );
     if ( z > 0.008856 ) z = pow( z, .3333333333 );
     else z = ( 7.787 * z ) + ( 16.0 / 116.0 );
     Lab[0] = ( 116.0 * y ) - 16.0;
     Lab[1] = 500.0 * ( x - y );
     Lab[2] = 200.0 * ( y - z );    }

void lab2lch( const vec<double>& Lab, vec<double>& LCH )
{    LCH.resize(3);
     LCH[0] = Lab[0];
     LCH[1] = sqrt( ( Lab[1] * Lab[1] ) + ( Lab[2] * Lab[2] ) );
     LCH[2] = atan2( Lab[2], Lab[1] );    }

double deltaE2000_bgr( const vec<double>& bgr1, const vec<double>& bgr2 )
{    vec<double> xyz1, xyz2, lab1, lab2, lch1, lch2;
     bgr2xyz( bgr1, xyz1 );
     bgr2xyz( bgr2, xyz2 );
     xyz2lab( xyz1, lab1 );
     xyz2lab( xyz2, lab2 );
     lab2lch( lab1, lch1 );
     lab2lch( lab2, lch2 );
     return deltaE2000( lch1, lch2 );    }

double deltaE2000( const vec<double>& lch1, const vec<double>& lch2 )
{    double avg_L = ( lch1[0] + lch2[0] ) * 0.5;
     double delta_L = lch2[0] - lch1[0];
     double avg_C = ( lch1[1] + lch2[1] ) * 0.5;
     double delta_C = lch1[1] - lch2[1];
     double avg_H = ( lch1[2] + lch2[2] ) * 0.5;
     if ( fabs( lch1[2] - lch2[2] ) > CV_PI ) avg_H += CV_PI;
     double delta_H = lch2[2] - lch1[2];
     if ( fabs( delta_H ) > CV_PI )
     {    if ( lch2[2] <= lch1[2] ) delta_H += CV_PI * 2.0;
          else delta_H -= CV_PI * 2.0;    }
     delta_H = sqrt( lch1[1] * lch2[1] ) * sin( delta_H ) * 2.0;
     double T = 1.0 -
          0.17 * cos( avg_H - CV_PI / 6.0 ) +
          0.24 * cos( avg_H * 2.0 ) +
          0.32 * cos( avg_H * 3.0 + CV_PI / 30.0 ) -
          0.20 * cos( avg_H * 4.0 - CV_PI * 7.0 / 20.0 );
     double SL = avg_L - 50.0;
     SL *= SL;
     SL = SL * 0.015 / sqrt( SL + 20.0 ) + 1.0;
     double SC = avg_C * 0.045 + 1.0;
     double SH = avg_C * T * 0.015 + 1.0;
     double delta_Theta = avg_H / 25.0 - CV_PI * 11.0 / 180.0;
     delta_Theta = exp( delta_Theta * -delta_Theta ) * ( CV_PI / 6.0 );
     double RT = pow( avg_C, 7.0 );
     RT = sqrt( RT / ( RT + 6103515625.0 ) ) 
          * sin( delta_Theta ) * -2.0; // 6103515625 = 25^7
     delta_L /= SL;
     delta_C /= SC;
     delta_H /= SH;
     return sqrt( delta_L * delta_L + delta_C * delta_C 
          + delta_H * delta_H + RT * delta_C * delta_H );    }

double ColorDist( const vec<int>& rgbpc1, const vec<int>& rgbpc2 )
{    vec<double> bgr1(3), bgr2(3);
     for ( int j = 0; j < 3; j++ )
     {    bgr1[j] = 255.0 * rgbpc1[3-j-1] / 100.0;
          bgr2[j] = 255.0 * rgbpc2[3-j-1] / 100.0;    }
     double d = deltaE2000_bgr( bgr1, bgr2 );
     /*
     if ( d < 0 ) d = 0;
     if ( d > 100 ) d = 100;
     */
     return d;    }
