/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// GetCodingMask.  Build a vecbitvector corresponding to human build 18, with a 
// given bit turned on if it corresponds to a coding base.  This file has
// hard-coded paths!
//
// *********************************************************************************
//
// The source file Ensembl48.genes was generated as follows:
//
// Go to www.ensembl.org, choose 'Mine Ensembl with Biomart', Database 'Ensembl
// 47', 'Homo sapiens genes', 'Attributes', 'Structures', and then choose the
// fields you want -
// 
// You must checkmark the boxes in the same order as appears below, as it
// determines the format of the data.
// 
// Ensembl Gene ID
// chromosome
// biotype
// Ensembl Transcript ID
// transcript start
// transcript end
// strand
// external gene id
// description
// exon start
// exon end
// ensembl exon id
// coding start
// coding end
// 
// Hit the Results button then choose
// 
// compressed file, TSV
//
// *********************************************************************************

#include "Basevector.h"
#include "Bitvector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "TokenizeString.h"

int main( )
{
     // HARD-CODED PATHS!  WARNING!

     String output_file = "/wga/scr3/jaffe/human_coding.vecbitvector";
     String gene_list = "/seq/solexa/Ensembl48.genes";
     String genome_fasta = "/seq/references/Homo_sapiens_assembly18/v0/"
          "Homo_sapiens_assembly18.fasta";

     // Build the file.

     fast_ifstream in(gene_list);
     vecbasevector genome;
     vecString genome_names;
     FetchReads( genome, genome_names, genome_fasta );
     vecbitvector coding;
     Mimic( genome, coding );
     String line, chr;
     vec<String> fields;
     vec<char> seps;
     seps.push_back( '\t' );
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( !line.Contains( "protein_coding" ) ) continue;
          TokenizeStrictly( line, seps, fields );
          if ( fields[12].size( ) == 0 || fields[13].size( ) == 0 ) continue;
          int start = fields[12].Int( ) - 1, stop = fields[13].Int( );
          size_t t = ~0UL;
          for ( size_t i = 0; i < genome_names.size( ); i++ )
          {    if ( genome_names[i] == "chr" + fields[1] )
               {    t = i;
                    break;    }    }
          if ( t != ~0UL )
          {    for ( int j = start; j < stop; j++ )
                    coding[t].Set( j, True );    }    }
     cout << "coding fraction = " << Coverage(coding) << endl;
     coding.WriteAll(output_file);    }
