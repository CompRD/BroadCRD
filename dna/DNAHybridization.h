/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Stuff relating to melting temperature of DNA, etc.

#ifndef DNA_HYBRIDIZATION_H
#define DNA_HYBRIDIZATION_H

#include "Basevector.h"
#include "CoreTools.h"

//  =============================================================
/// Calculate melting temperature of a DNA string.  Default salt
/// concentration is 0.05 M.
/// 
/// \function Tm
/// 
/// Formulas from
/// http://www.basic.northwestern.edu/biotools/oligocalc.html and from
/// Andi Gnirke and from M.L.M. Anderson, _Nucleic Acid Hybridization_.
///
/// For long molecules, > 70 bases, with GC content 30-75%:
/// Tm= 81.5 + (41 * (yG+zC)/(wA+xT+yG+zC)) - (500/(wA+xT+yG+zC)) + 16.6*log10([Na+])
///
/// For molecules 20 to 70 bases:
/// Tm= 81.5 + (41 * (yG+zC)/(wA+xT+yG+zC)) - (600/(wA+xT+yG+zC)) + 16.6*log10([Na+])
/// For oligos < 20 bases:
///
///  Tm= (wA+xT)*2 + (yG+zC)*4 - 16.6*log10(0.050) + 16.6*log10([Na+])

/// Unused alternative formula for oligos >14 < 50:
///
/// Tm= 100.5 + (41 * (yG+zC)/(wA+xT+yG+zC)) - (820/(wA+xT+yG+zC)) + 16.6*log10([Na+])
///
//  ==============================================================

float Tm( const String & s, float NaMol = 0.05 );
float Tm( const basevector & b, float NaMol = 0.05 );

// Tm_NearestNeighbor: compute melting temperature of a given DNA sequence S 
// at molarity SMol (default 0.25 microM), where Na+ concentration is NaMol
// (default 50 mM), using a nearest-neighbor model, as used by
//
// http://www.idtdna.com/ANALYZER/Applications/OligoAnalyzer
//
// and described in
//
// [1] Allawi, H. T. and SantaLucia J. Jr.  Thermodynamics and NMR of internal
// G.T mismatches in DNA.  Biochemistry 36 (1997), 10581-10594.
//
// [2] Owczarzy, et al.  Effects of sodium ions on DNA duplex oligomers: improved
// predictions of melting temperatures.  Biochemistry 43 (2004), 3537-3554.
//
// [3] Patricia M. McTigue, Raymond J. Peterson, Jason D. Kuhn.  Sequence-dependent
// thermodynamic parameters for locked nucleic acid (LNA)-DNA duplex formation.
// Biochemistry 43 (2004), 5388-5405.
// (See also www.celadonlabs.com/Software/ModChem/default.aspx.)
//
// This method should be quite good for short sequences.
//
// The results are close to (but not identical to) the results given by the IDT
// OligoAnalyzer.  I don't know what the discrepancy is due to.
//
// We allow some nucleotides to be locked, as specified by the vector "locked".
// The following assumptions are enforced: 
// (a) no two locked bases in a row;
// (b) no locked base at the beginning or end, or adjacent to those positions.
// These are consistent with the assumptions in [3].
//
// Alternatively, if S contains + symbols, they will be interpreted as instructions
// to "lock the following nucleotide", and the vector "locked" will be generated
// accordingly.

double Tm_NearestNeighbor( const String& S, 
     const double SMol = 0.00000025, const double NaMol = 0.05,
     const vec<Bool>& locked = vec<Bool>( ) );

inline double Tm_NearestNeighbor( const basevector& b, 
     const double SMol = 0.00000025, const double NaMol = 0.05,
     const vec<Bool>& locked = vec<Bool>( ) )
{    return Tm_NearestNeighbor( b.ToString( ), SMol, NaMol, locked );    }

// FreeEnergyDNA.  Compute the delta_G that you get if you ignore the symmetry and
// initiation terms.

double FreeEnergyDNA( const String& S );

// FreeEnergyDNA_IDT_Heterodimer(S).  Compute the delta_G that you get if you ignore 
// the symmetry and initiation terms, using the nearest neighbor coefficients that
// are implied by the idtdna website computation for heterodimers 
// (at least as of Sept. 6, 2006).  Found by reverse engineering the web site.

double FreeEnergyDNA_IDT_Heterodimer( const String& S );

// FreeEnergyDNA_IDT_Heterodimer( S, T ).  Find the longest perfect match M between
// S and T, and then find all perfect matches P whose length is >= |M| - 2, and
// compute FreeEnergyDNA_IDT_Heterodimer(P).  Return the max of all such values.
//
// This is very similar to what the idtdna website returns if you ask it to compute
// delta G for a "heterodimer".  Note that the site does not use the oligo 
// concentration or salt concentration.
//
// Note that this must be a very bad way to compute the actual maximum delta G,
// because it only looks at perfect matches.

double FreeEnergyDNA_IDT_Heterodimer( const String& S, const String& T );

// Same but return second best free energy:

double FreeEnergyDNA_IDT_Heterodimer_SecondBest( const String& S, const String& T );

#endif
