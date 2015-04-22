/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "dna/DNAHybridization.h"
#include "math/Functions.h"

float Tm(const String & s, float NaMol) {
  float w = 0, x = 0, y = 0, z = 0;
  for (int i=0; i != int(s.size()); ++i) {
    switch (toupper(s[i])) {
    case 'A': ++w; break;
    case 'C': ++z; break;
    case 'G': ++y; break;
    case 'T': ++x; break;
    }
  }
  const float N = (w+x+y+z);
  const float NaTerm = 16.6*log10(NaMol);
  if (N>70)
    return 81.5 + (41.0 * (y+z)/N) - (500.0/N) + NaTerm;
  if (N>20)
    return 81.5 + (41.0 * (y+z)/N) - (600.0/N) + NaTerm;
  return (w+x)*2 + (y+z)*4 - 16.6*log10(0.050) + NaTerm;
}

float Tm( const basevector & b, float NaMol) {
  return Tm(b.ToString(), NaMol);
}

// GetThermodynamicParametersDNA.
//
// Return nearest neighbor thermodynamic parameters, from table 1 in reference [1],
// Allawi and SantaLucia (1997).
//
// dH = enthalpy (kcal/mole) = dH^\circ
// dS = entropy = dS^\circ
// dG = Gibb's free energy at 37 degrees = dG^\circ_37
// dH, dS, dG are all {A,C,G,T} x {A,C,G,T} matrices.
// (d means delta)
//
// Although there are dG entries in the table, I've computed them here directly
// from dH and dS.  I checked a few entries and they agree with what's in the table.
// dG = dH - T dS (theoretical)
// dG = dH - (37 + 273.15) dS / 1000 (used).
// I'm not sure quite why the factor of 1000 is needed.

void GetThermodynamicParametersDNA( vec< vec<double> >& dH,
     vec< vec<double> >& dS, vec< vec<double> >& dG,
     double& dH_G_or_C_init, double& dH_A_or_T_init,
     double& dS_G_or_C_init, double& dS_A_or_T_init,
     double& dG_G_or_C_init, double& dG_A_or_T_init,
     double& dH_symmetry_correction, double& dS_symmetry_correction,
     double& dG_symmetry_correction )
{
     int A = 0, C = 1, G = 2, T = 3;
     dH.clear_and_resize(4);
     dS.clear_and_resize(4);
     dG.clear_and_resize(4);
     for ( int i = 0; i < 4; i++ )
     {    dH[i].resize(4), dS[i].resize(4), dG[i].resize(4);    }

     dH[A][A] = dH[T][T] = -7.9;
     dH[A][T] = -7.2;
     dH[T][A] = -7.2;
     dH[C][A] = dH[T][G] = -8.5;
     dH[G][T] = dH[A][C] = -8.4;
     dH[C][T] = dH[A][G] = -7.8;
     dH[G][A] = dH[T][C] = -8.2;
     dH[C][G] = -10.6;
     dH[G][C] = -9.8;
     dH[G][G] = dH[C][C] = -8.0;
     dS[A][A] = dS[T][T] = -22.2;
     dS[A][T] = -20.4;
     dS[T][A] = -21.3;
     dS[C][A] = dS[T][G] = -22.7;
     dS[G][T] = dS[A][C] = -22.4;
     dS[C][T] = dS[A][G] = -21.0;
     dS[G][A] = dS[T][C] = -22.2;
     dS[C][G] = -27.2;
     dS[G][C] = -24.4;
     dS[G][G] = dS[C][C] = -19.9;

     /*
     // Here are the parameters from SantaLucia and Hicks (2004).
     // They are almost the same.
     dH[A][A] = dH[T][T] = -7.6;
     dH[A][T] = -7.2;
     dH[T][A] = -7.2;
     dH[C][A] = dH[T][G] = -8.5;
     dH[G][T] = dH[A][C] = -8.4;
     dH[C][T] = dH[A][G] = -7.8;
     dH[G][A] = dH[T][C] = -8.2;
     dH[C][G] = -10.6;
     dH[G][C] = -9.8;
     dH[G][G] = dH[C][C] = -8.0;
     dS[A][A] = dS[T][T] = -21.3;
     dS[A][T] = -20.4;
     dS[T][A] = -21.3;
     dS[C][A] = dS[T][G] = -22.7;
     dS[G][T] = dS[A][C] = -22.4;
     dS[C][T] = dS[A][G] = -21.0;
     dS[G][A] = dS[T][C] = -22.2;
     dS[C][G] = -27.2;
     dS[G][C] = -24.4;
     dS[G][G] = dS[C][C] = -19.9;
     */

     double temp = 37.0 + 273.15;
     for ( int i = 0; i < 4; i++ )
     {    for ( int j = 0; j < 4; j++ )
          {    dG[i][j] = dH[i][j] - temp * dS[i][j] / 1000.0;    }    }
     dH_G_or_C_init = 0.1;
     dH_A_or_T_init = 2.3;
     dS_G_or_C_init = -2.8;
     dS_A_or_T_init = 4.1;
     dG_G_or_C_init = dH_G_or_C_init - temp * dS_G_or_C_init / 1000.0;
     dG_A_or_T_init = dH_A_or_T_init - temp * dS_A_or_T_init / 1000.0;
     dH_symmetry_correction = 0.0;
     dS_symmetry_correction = -1.4;
     dG_symmetry_correction
          = dH_symmetry_correction - temp * dS_symmetry_correction / 1000.0;    }

// GetLockedThermodynamicParametersDNA.  Return the 32 nearest-neighbor
// thermodynamic parameter set for LNA incorporation, from the left half of Table 4
// in reference [3].

void GetLockedThermodynamicParametersDNA( vec< vec<double> >& ddH_left,
     vec< vec<double> >& ddS_left, vec< vec<double> >& ddG_left,
     vec< vec<double> >& ddH_right, vec< vec<double> >& ddS_right,
     vec< vec<double> >& ddG_right )
{    int A = 0, C = 1, G = 2, T = 3;
     ddH_left.clear_and_resize(4);
     ddS_left.clear_and_resize(4);
     ddG_left.clear_and_resize(4);
     ddH_right.clear_and_resize(4);
     ddS_right.clear_and_resize(4);
     ddG_right.clear_and_resize(4);
     for ( int i = 0; i < 4; i++ )
     {    ddH_left[i].resize(4), ddS_left[i].resize(4), ddG_left[i].resize(4);
          ddH_right[i].resize(4), ddS_right[i].resize(4), ddG_right[i].resize(4);   }
     ddH_left[A][A] = 0.707;
     ddH_left[A][C] = 1.131;
     ddH_left[A][G] = 0.264;
     ddH_left[A][T] = 2.282;
     ddH_left[C][A] = 1.049;
     ddH_left[C][C] = 2.096;
     ddH_left[C][G] = 0.785;
     ddH_left[C][T] = 0.708;
     ddH_left[G][A] = 3.162;
     ddH_left[G][C] = -0.360;
     ddH_left[G][G] = -2.844;
     ddH_left[G][T] = -0.212;
     ddH_left[T][A] = -0.046;
     ddH_left[T][C] = 1.893;
     ddH_left[T][G] = -1.540;
     ddH_left[T][T] = 1.528;
     ddS_left[A][A] = 2.477;
     ddS_left[A][C] = 4.064;
     ddS_left[A][G] = 2.613;
     ddS_left[A][T] = 7.457;
     ddS_left[C][A] = 4.320;
     ddS_left[C][C] = 7.996;
     ddS_left[C][G] = 3.709;
     ddS_left[C][T] = 4.175;
     ddS_left[G][A] = 10.544;
     ddS_left[G][C] = -0.251;
     ddS_left[G][G] = -6.680;
     ddS_left[G][T] = 0.073;
     ddS_left[T][A] = 1.562;
     ddS_left[T][C] = 6.685;
     ddS_left[T][G] = -3.044;
     ddS_left[T][T] = 5.298;
     ddG_left[A][A] = -0.092;
     ddG_left[A][C] = -0.122;
     ddG_left[A][G] = -0.561;
     ddG_left[A][T] = -0.007;
     ddG_left[C][A] = -0.270;
     ddG_left[C][C] = -0.457;
     ddG_left[C][G] = -0.332;
     ddG_left[C][T] = -0.666;
     ddG_left[G][A] = -0.072;
     ddG_left[G][C] = -0.414;
     ddG_left[G][G] = -0.700;
     ddG_left[G][T] = -0.194;
     ddG_left[T][A] = -0.563;
     ddG_left[T][C] = -0.208;
     ddG_left[T][G] = -0.548;
     ddG_left[T][T] = -0.130;
     ddH_right[A][A] = 0.992;
     ddH_right[A][C] = 2.890;
     ddH_right[A][G] = -1.200;
     ddH_right[A][T] = 1.816;
     ddH_right[C][A] = 1.358;
     ddH_right[C][C] = 2.063;
     ddH_right[C][G] = -0.276;
     ddH_right[C][T] = -1.671;
     ddH_right[G][A] = 0.444;
     ddH_right[G][C] = -0.925;
     ddH_right[G][G] = -0.943;
     ddH_right[G][T] = -0.635;
     ddH_right[T][A] = 1.591;
     ddH_right[T][C] = 0.609;
     ddH_right[T][G] = 2.165;
     ddH_right[T][T] = 2.326;
     ddS_right[A][A] = 4.065;
     ddS_right[A][C] = 10.576;
     ddS_right[A][G] = -1.826;
     ddS_right[A][T] = 6.863;
     ddS_right[C][A] = 4.367;
     ddS_right[C][C] = 7.565;
     ddS_right[C][G] = -0.718;
     ddS_right[C][T] = -4.070;
     ddS_right[G][A] = 2.898;
     ddS_right[G][C] = -1.111;
     ddS_right[G][G] = -0.933;
     ddS_right[G][T] = -0.342;
     ddS_right[T][A] = 5.281;
     ddS_right[T][C] = 3.169;
     ddS_right[T][G] = 7.163;
     ddS_right[T][T] = 8.051;
     ddG_right[A][A] = -0.396;
     ddG_right[A][C] = -0.390;
     ddG_right[A][G] = -0.603;
     ddG_right[A][T] = -0.309;
     ddG_right[C][A] = 0.046;
     ddG_right[C][C] = -0.404;
     ddG_right[C][G] = -0.003;
     ddG_right[C][T] = -0.409;
     ddG_right[G][A] = -0.437;
     ddG_right[G][C] = -0.535;
     ddG_right[G][G] = -0.666;
     ddG_right[G][T] = -0.520;
     ddG_right[T][A] = 0.004;
     ddG_right[T][C] = -0.396;
     ddG_right[T][G] = -0.106;
     ddG_right[T][T] = -0.212;    }

// ThermodynamicSumsDNA.  Compute dH_sum, dS_sum, dG_sum.  We allow some nucleotides
// to be locked, as specified by the variable "locked".  The following assumptions
// are enforced:
// (a) no two locked bases in a row;
// (b) no locked base at the beginning or end, or adjacent to those positions.
// These are consistent with the assumptions in [3].

void ThermodynamicSumsDNA( const String& S, double& dH_sum, double& dS_sum,
     double& dG_sum, const Bool include_symmetry_correction = True,
     const Bool include_initiation_terms = True,
     const vec<Bool>& locked = vec<Bool>( ) )
{    vec< vec<double> > dH, dS, dG;
     double dH_G_or_C_init, dH_A_or_T_init;
     double dS_G_or_C_init, dS_A_or_T_init;
     double dG_G_or_C_init, dG_A_or_T_init;
     double dH_symmetry_correction, dS_symmetry_correction, dG_symmetry_correction;
     GetThermodynamicParametersDNA( dH, dS, dG, dH_G_or_C_init, dH_A_or_T_init,
          dS_G_or_C_init, dS_A_or_T_init, dG_G_or_C_init, dG_A_or_T_init,
          dH_symmetry_correction, dS_symmetry_correction, dG_symmetry_correction );
     dH_sum = dS_sum = dG_sum = 0.0;
     if (include_symmetry_correction)
     {    dH_sum += dH_symmetry_correction;
          dS_sum += dS_symmetry_correction;
          dG_sum += dG_symmetry_correction;    }
     if (include_initiation_terms)
     {    if ( S[0] == 'A' || S[0] == 'T' )
          {    dH_sum += dH_A_or_T_init;
               dS_sum += dS_A_or_T_init;
               dG_sum += dG_A_or_T_init;    }
          else
          {    dH_sum += dH_G_or_C_init;
               dS_sum += dS_G_or_C_init;
               dG_sum += dG_G_or_C_init;    }
          if ( S[ (int) S.size( ) - 1 ] == 'A' || S[ (int) S.size( ) - 1 ] == 'T' )
          {    dH_sum += dH_A_or_T_init;
               dS_sum += dS_A_or_T_init;
               dG_sum += dG_A_or_T_init;    }
          else
          {    dH_sum += dH_G_or_C_init;
               dS_sum += dS_G_or_C_init;
               dG_sum += dG_G_or_C_init;    }    }
     basevector b;
     b.SetFromString(S);
     for ( int i = 0; i < S.isize( ) - 1; i++ )
     {    dH_sum += dH[ b[i] ][ b[i+1] ];
          dS_sum += dS[ b[i] ][ b[i+1] ];
          dG_sum += dG[ b[i] ][ b[i+1] ];    }

     // Handle locked bases.

     if ( locked.empty( ) || Sum(locked) == 0 ) return;
     for ( int i = 1; i < locked.isize( ); i++ )
          if ( locked[i] ) ForceAssert( !locked[i-1] );
     ForceAssertGe( locked.isize( ), 5 );
     ForceAssert( !locked[0] && !locked[1] );
     ForceAssert( !locked[ locked.isize( ) - 1 ] );
     ForceAssert( !locked[ locked.isize( ) - 2 ] );
     vec< vec<double> > ddH_left, ddS_left, ddG_left;
     vec< vec<double> > ddH_right, ddS_right, ddG_right;
     GetLockedThermodynamicParametersDNA( ddH_left, ddS_left, ddG_left,
          ddH_right, ddS_right, ddG_right );
     for ( int i = 0; i < locked.isize( ); i++ )
     {    if ( !locked[i] ) continue;
          dH_sum += ddH_left[ b[i] ][ b[i+1] ] + ddH_right[ b[i-1] ][ b[i] ];
          dS_sum += ddS_left[ b[i] ][ b[i+1] ] + ddS_right[ b[i-1] ][ b[i] ];
          dG_sum += ddG_left[ b[i] ][ b[i+1] ]
               + ddG_right[ b[i-1] ][ b[i] ];    }    }

void VerifyDNA( const String& S )
{ String::const_iterator end = S.end();
  for (String::const_iterator itr = S.begin(); itr != end; ++itr )
          ForceAssert( Base::isCanonicalBase(*itr) ); }

double FreeEnergyDNA( const String& S )
{    VerifyDNA(S);
     double dH_sum, dS_sum, dG_sum;
     ThermodynamicSumsDNA( S, dH_sum, dS_sum, dG_sum, False, False );
     return dG_sum;    }

double FreeEnergyDNA_IDT_Heterodimer( const String& S )
{    VerifyDNA(S);
     double dG_sum = 0.0;
     basevector b;
     b.SetFromString(S);
     vec< vec<double> > dG;
     dG.clear_and_resize(4);
     for ( int i = 0; i < 4; i++ )
          dG[i].resize(4);
     int A = 0, C = 1, G = 2, T = 3;

     // I have absolutely no idea where the following parameters come from.
     // (I got them by reverse engineering.)

     dG[A][A] = dG[T][T] = -1.94;
     dG[A][T] = -1.47;
     dG[T][A] = -0.96;
     dG[C][A] = dG[T][G] = -1.95;
     dG[G][T] = dG[A][C] = -1.34;
     dG[C][T] = dG[A][G] = -1.60;
     dG[G][A] = dG[T][C] = -1.57;
     dG[C][G] = -3.61;
     dG[G][C] = -3.14;
     dG[G][G] = dG[C][C] = -3.07;

     for ( int i = 0; i < (int) S.size( ) - 1; i++ )
          dG_sum += dG[ b[i] ][ b[i+1] ];
     return dG_sum;    }

double FreeEnergyDNA_IDT_Heterodimer( const String& S, const String& T )
{    int nS = S.size( ), nT = T.size( );
     String Trc;
     StringReverseComplement( T, Trc );
     int max_match = 0;
     double min_dG = 0.0;
     for ( int pass = 1; pass <= 2; pass++ )
     {    for ( int offset = - nT + 2; offset <= nS - 2; offset++ )
          {    int spos_start, spos_stop, tpos_start, tpos_stop;
               if ( offset >= 0 )
               {    spos_start = offset;
                    spos_stop = Min( nS, offset + nT );
                    tpos_start = 0;
                    tpos_stop = Min( nT, nS - offset );    }
               else
               {    spos_start = 0;
                    spos_stop = Min( nS, offset + nT );
                    tpos_start = -offset;
                    tpos_stop = Min( nT, -offset + nS );    }
               for ( int js = spos_start; js < spos_stop; js++ )
               {    int jt = js + tpos_start - spos_start;
                    if ( S[js] != Trc[jt] ) continue;
                    int ks;
                    for ( ks = js + 1; ks < spos_stop; ks++ )
                    {    int kt = ks + tpos_start - spos_start;
                         if ( S[ks] != Trc[kt] ) break;    }
                    int match = ks - js;
                    if ( pass == 1 ) max_match = Max( max_match, match );
                    else if ( match >= 2 && max_match - match <= 2 )
                    {    String M;
                         for ( int u = js; u < ks; u++ )
                              M += S[u];
                         double dG = FreeEnergyDNA_IDT_Heterodimer(M);
                         min_dG = Min( min_dG, dG );    }
                    js = ks;    }    }    }
     return min_dG;    }

double FreeEnergyDNA_IDT_Heterodimer_SecondBest( const String& S, const String& T )
{    int nS = S.size( ), nT = T.size( );
     String Trc;
     StringReverseComplement( T, Trc );
     int max_match = 0;
     vec<double> dGs;
     for ( int pass = 1; pass <= 2; pass++ )
     {    for ( int offset = - nT + 2; offset <= nS - 2; offset++ )
          {    int spos_start, spos_stop, tpos_start, tpos_stop;
               if ( offset >= 0 )
               {    spos_start = offset;
                    spos_stop = Min( nS, offset + nT );
                    tpos_start = 0;
                    tpos_stop = Min( nT, nS - offset );    }
               else
               {    spos_start = 0;
                    spos_stop = Min( nS, offset + nT );
                    tpos_start = -offset;
                    tpos_stop = Min( nT, -offset + nS );    }
               for ( int js = spos_start; js < spos_stop; js++ )
               {    int jt = js + tpos_start - spos_start;
                    if ( S[js] != Trc[jt] ) continue;
                    int ks;
                    for ( ks = js + 1; ks < spos_stop; ks++ )
                    {    int kt = ks + tpos_start - spos_start;
                         if ( S[ks] != Trc[kt] ) break;    }
                    int match = ks - js;
                    if ( pass == 1 ) max_match = Max( max_match, match );
                    else if ( match >= 2 )
                    {    String M;
                         for ( int u = js; u < ks; u++ )
                              M += S[u];
                         double dG = FreeEnergyDNA_IDT_Heterodimer(M);
                         dGs.push_back(dG);    }
                    js = ks;    }    }    }
     Sort(dGs);
     if ( dGs.size( ) >= 2 ) return dGs[1];
     else return 0.0;    }

double Tm_NearestNeighbor( const String& S, const double SMol, const double NaMol,
     const vec<Bool>& locked )
{
     // Allow for + symbols.

     if ( S.Contains( "+" ) )
     {    ForceAssert( locked.empty( ) );
          ForceAssert( !S.Contains( "++" ) );
          ForceAssert( !S.Contains( "+", -1 ) );
          String Sx;
          vec<Bool> lockedx;
          for ( int i = 0; i < S.isize( ); i++ )
          {    if ( S[i] != '+' )
               {    Sx += S[i];
                    lockedx.push_back(False);    }
               else
               {    Sx += S[i+1];
                    lockedx.push_back(True);
                    ++i;    }    }
          return Tm_NearestNeighbor( Sx, SMol, NaMol, lockedx );    }

     // Validate sequence.

     VerifyDNA(S);

     // Compute thermodynamic sums.

     double dH_sum, dS_sum, dG_sum;
     ThermodynamicSumsDNA( S, dH_sum, dS_sum, dG_sum, True, True, locked );

     // Compute melting temperature based on nearest-neighbor model.

     const double ideal_gas_const = 1.987; // calories per Kelvin per mole
     const double Kelvin_to_Celsius = 273.15;
     double temp = 1000.0 *
          dH_sum / ( dS_sum + ideal_gas_const * log(SMol) ) - Kelvin_to_Celsius;

     // Correct for Na concentration, following [2].

     int GC = 0;
     for ( int i = 0; i < (int) S.size( ); i++ )
          if ( S[i] == 'G' || S[i] == 'C' ) ++GC;
     double GC_fract = double(GC) / double( S.size( ) );
     double lnNaMol = log(NaMol);
     temp = -Kelvin_to_Celsius + 1.0 /
          ( 1.0/(temp+Kelvin_to_Celsius)
               + ( 4.29 * GC_fract - 3.95 ) * 0.00001 * lnNaMol
               + 9.40 * 0.000001 * lnNaMol * lnNaMol );
     return temp;    }
