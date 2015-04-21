///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file CleanFasta.h
 * \author tsharpe
 * \date Feb 16, 2012
 *
 * \brief
 */
#ifndef CLEANFASTA_H_
#define CLEANFASTA_H_

#include "String.h"

/// Copies a FASTA file, removing all N's from the sequence.
void CleanFasta( String const& fastaIn, String const& fastaOut );

#endif /* CLEANFASTA_H_ */
