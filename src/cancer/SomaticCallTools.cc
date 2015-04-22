/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// require <cstddef> to use gmp in GCC-4.9
#include <cstddef>
#include <gmpxx.h>

#include "CoreTools.h"
#include "cancer/SomaticCallTools.h"
#include "random/Bernoulli.h"
#include "util/AlleleLikelihood.h"

// MakeDepend: library GMPXX

typedef mpf_class big_float;

Bool SomaticMutation( 

     // inputs:

     const int chr, const int pos, const vecbasevector& ref,
     const vec<unsigned char>* TQ, const vec<unsigned char>* NQ, 
     const vec<unsigned char>* TQ_rc, const vec<unsigned char>* NQ_rc, 
     const unsigned char refbase, 

     // heuristic constants:

     const int min_mutant_sum_pretest,
     const int min_mutant_sum,
     const double tumor_threshold,
     const double normal_threshold,

     // input and output (if != 4, assume altbase is given):

     unsigned char& altbase, 

     // outputs (apart from return value):

     String& info

     )
{    
     // Test sum at alternate base.

     vec<int> Tsums( 4, 0 );
     for ( int k = 0; k < 4; k++ )
     {    for ( int j = 0; j < TQ[k].isize( ); j++ )
               Tsums[k] += TQ[k][j];    }
     int Trefsum = Tsums[refbase], Taltsum = 0;
     if ( altbase == 4 )
     {    for ( int k = 0; k < 4; k++ )
          {    if ( k != refbase && Tsums[k] > Taltsum )
               {    Taltsum = Tsums[k];
                    altbase = k;    }    }    }
     else Taltsum = Tsums[altbase];
     if ( Taltsum < min_mutant_sum_pretest ) return False;

     // Compute asymmetry score.

     int fw = 0, rc = 0, fw_alt = 0, rc_alt = 0;
     for ( int pass = 1; pass <= 2; pass++ )
     {    const vec<unsigned char>*& Q_rc = ( pass == 1 ? TQ_rc : NQ_rc );
          for ( int j = 0; j < 4; j++ )
          {    for ( int k = 0; k < Q_rc[j].isize( ); k++ )
               {    if ( Q_rc[j][k] ) ++rc;
                    else ++fw;
                    if ( j == altbase ) 
                    {    if ( Q_rc[j][k] ) ++rc_alt;
                         else ++fw_alt;    }    }    }    }
     double asymmetry_score1 = -log10( BinomialSum( fw_alt + rc_alt, fw_alt,
          double(fw) / double(fw+rc) ) );
     double asymmetry_score2 = -log10( BinomialSum( fw_alt + rc_alt, rc_alt,
          double(rc) / double(fw+rc) ) );
     double asymmetry_score = Max( asymmetry_score1, asymmetry_score2 );

     // Compute genotype likelihoods.

     genotype_likelihoods NG;
     for ( int j = 0; j < 4; j++ )
     {    for ( int k = 0; k < NQ[j].isize( ); k++ )
          {    genotype_likelihoods NGk( as_base(refbase), as_base(j), NQ[j][k] );
               NG += NGk;    }    }
     genotype_likelihoods TG;
     for ( int j = 0; j < 4; j++ )
     {    for ( int k = 0; k < TQ[j].isize( ); k++ )
          {    genotype_likelihoods TGk( as_base(refbase), as_base(j), TQ[j][k] );
               TG += TGk;    }    }

     int Tcov = 0, Tsum = 0, Ncov = 0, Nsum = 0;
     for ( int j = 0; j < 4; j++ )
     {    for ( int k = 0; k < TQ[j].isize( ); k++ )
          {    Tcov++;
               Tsum += TQ[j][k];    }
          for ( int k = 0; k < NQ[j].isize( ); k++ )
          {    Ncov++;
               Nsum += NQ[j][k];    }    }

     vec<String> alleles;
     vec<double> likelihoods;
     NG.GetSorted( alleles, likelihoods );

     String refref(2);
     refref[0] = as_base(refbase), refref[1] = as_base(refbase);
     vec<unsigned char> all(2);
     all[0] = as_base(refbase), all[1] = as_base(altbase);
     Sort(all);
     String altref(2);
     altref[0] = all[0], altref[1] = all[1];
     String altalt(2);
     altalt[0] = as_base(altbase), altalt[1] = as_base(altbase);

     int AA = Position( alleles, refref );
     int AC = Position( alleles, altref );
     int CC = Position( alleles, altalt );
     double AC_CC = log10( pow( 10, likelihoods[AC] ) + pow( 10, likelihoods[CC] ) );
     double normal_score = likelihoods[AA] - AC_CC;

     if ( normal_score < normal_threshold ) return False;

     TG.GetSorted( alleles, likelihoods );
     AA = Position( alleles, refref );
     AC = Position( alleles, altref );
     CC = Position( alleles, altalt );
     AC_CC = log10( pow( 10, likelihoods[AC] ) + pow( 10, likelihoods[CC] ) );
     double tumor_score = AC_CC - likelihoods[AA];

     double Tq = double(Tsum)/double(Tcov); // mean observed quality in tumor
     double Tc = double(Taltsum)/Tq;
     int Tc0 = int(floor(Tc)), Tc1 = int(ceil(Tc));
     double Tqc;
     if ( Tc0 > 0 )
     {    double Tqc0 = 0, Tqc1 = 0;
          for ( int pass = 1; pass <= 2; pass++ )
          {    int Tc = ( pass == 1 ? Tc0 : Tc1 ); 
               big_float Tp = pow( 10, -(Tq/10) ) / 3;
               big_float Tx = 3 * ( 1 - BinomialSum( Tcov, Tc-1, Tp ) );
               ( pass == 1 ? Tqc0 : Tqc1 ) = -10 * log10( Tx.get_d( ) );    }
          Tqc = Tqc0 + (Tc - Tc0) * (Tqc1 - Tqc0);    }
     else Tqc = Taltsum;

     int Naltsum = 0;
     for ( int j = 0; j < NQ[altbase].isize( ); j++ )
          Naltsum += NQ[altbase][j];

     double Talt_score = Taltsum, Nalt_score = Naltsum;
     for ( int i = 0; i < 4; i++ )
     {    if ( i == altbase ) continue;
          for ( int j = 0; j < TQ[i].isize( ); j++ )
          {    double q = TQ[i][j];
               Talt_score -= q / ( 3.0 * pow( 10.0, q/10.0 ) );    }    
          for ( int j = 0; j < NQ[i].isize( ); j++ )
          {    double q = NQ[i][j];
               Nalt_score -= q / ( 3.0 * pow( 10.0, q/10.0 ) );    }    }
     Talt_score = Max( 0.0, Talt_score / 10.0 );
     Nalt_score = Max( 0.0, Nalt_score / 10.0 );

     big_float Talt_count = TQ[altbase].size( );
     int Nalt_count = NQ[altbase].size( );
     big_float Tcount 
          = TQ[0].size( ) + TQ[1].size( ) + TQ[2].size( ) + TQ[3].size( );
     int Ncount = NQ[0].size( ) + NQ[1].size( ) + NQ[2].size( ) + NQ[3].size( );

     big_float inv_ratio = Talt_count / Tcount;
     double invisible;
     if ( Talt_count < Tcount )
     {    big_float inv_sum = BinomialSum( Ncount, Nalt_count, inv_ratio );
          invisible = -log10( inv_sum.get_d( ) );    }
     else
     {    if ( Nalt_count < Ncount ) invisible = 1000.0; // infinite
          else invisible = 0.0;    }

     if ( Talt_score < tumor_threshold ) return False;
     if ( asymmetry_score + Nalt_score > normal_threshold ) return False;
     if ( invisible < 3.0 ) return False;

     double Nq = double(Nsum)/double(Ncov); // mean observed quality in tumor
     double Nc = double(Naltsum)/Nq;
     int Nc0 = int(floor(Nc)), Nc1 = int(ceil(Nc));
     double Nqc;
     if ( Nc0 > 0 )
     {    double Nqc0 = 0, Nqc1 = 0;
          for ( int pass = 1; pass <= 2; pass++ )
          {    int Nc = ( pass == 1 ? Nc0 : Nc1 ); 
               big_float Np = pow( 10, -(Nq/10) ) / 3;
               big_float Nx = 3 * ( 1 - BinomialSum( Ncov, Nc-1, Np ) );
               ( pass == 1 ? Nqc0 : Nqc1 ) = -10 * log10( Nx.get_d( ) );    }
          Nqc = Nqc0 + (Nc - Nc0) * (Nqc1 - Nqc0);    }
     else Nqc = Naltsum;

     // Get sequence context.

     String context; // note sloppy def: other bases may have changed
     if ( pos-10 >= 0 && pos+10 < ref[chr].isize( ) )
     {    basevector c;
          c.SetToSubOf( ref[chr], pos-10, 21 );
          c.Set( 10, altbase );
          context = c.ToString( );    }
     else context = "near_end";
    
     // Report success.

     int ocount = 3;
     ostringstream out;
     out << "reference_fasta_record_" << chr << ":" << pos;
     #define OPARAM( NAME, DEF )                                            \
          out << " " << ocount << "." << #NAME << " " << DEF; ocount += 2;
     OPARAM( refbase, as_base(refbase) )
     OPARAM( altbase, as_base(altbase) )
     OPARAM( tumor_context, context )
     OPARAM( T_A, TQ[0].size( ) )
     OPARAM( T_C, TQ[1].size( ) )
     OPARAM( T_G, TQ[2].size( ) )
     OPARAM( T_T, TQ[3].size( ) )
     OPARAM( N_A, NQ[0].size( ) )
     OPARAM( N_C, NQ[1].size( ) )
     OPARAM( N_G, NQ[2].size( ) )
     OPARAM( N_T, NQ[3].size( ) )
     OPARAM( Tscore, tumor_score )
     OPARAM( Tcov, Tcov )
     OPARAM( Tsum, Tsum )
     OPARAM( Tref_sum, Trefsum )
     OPARAM( Talt_sum, Taltsum )
     OPARAM( Talt_sum_adj, Tqc )
     OPARAM( Nscore, normal_score )
     OPARAM( Ncov, Ncov )
     OPARAM( Nalt_sum, Naltsum )
     OPARAM( Nalt_sum_adj, Nqc )
     OPARAM( Talt_score, Talt_score )
     OPARAM( Nalt_score, Nalt_score )
     OPARAM( asymmetry_score, asymmetry_score )
     info = out.str( );
     return True;    }
