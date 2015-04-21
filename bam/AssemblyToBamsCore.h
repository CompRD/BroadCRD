///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BAM__ASSEMBLY_TO_BAMS_CORE__H
#define BAM__ASSEMBLY_TO_BAMS_CORE__H

#include "String.h"

/**
 * CheckPerlVersion
 *
 * Make sure perl v5.10 is in the path (a warning is sent to cout if
 * not).
 */
bool CheckPerlVersion( const String &tmp_dir );

/**
 * CheckJavaVersion
 *
 * Make sure java 1.6 is in the path (a warning is sent to cout if
 * not).
 */
bool CheckJavaVersion( const String &tmp_dir );

/**
 * PrepareReference
 *
 * Generate various dictionary and lookup files needed by different
 * aligners. Output is saved in the subdir <base_dir>/reference.
 *
 * ref_fasta: full path name of input fasta
 * base_dir: output will be saved in the subdir <base_dir>/reference
 * createSeqDict_jar: full path name to picard's CreateSequenceDictionary.jar
 * bwa_bin: full path name to bwa binary
 * FORCE: force-regenerate all output files
 */
void PrepareReference( const String &ref_fasta,
		       const String &createSeqDict_jar,
		       const String &bwa_bin,
		       const String &base_dir,
		       const bool FORCE = false );

/**
 * AlignReadsToReference
 *
 * Use bwa to align reads o reference, and generate .bam and .bam.bai
 * files.
 *
 * in_bam_file: input bam file
 * base_dir, head_name: define output dir (<base_dir>/<head_name>)
 * picard_dir: where picard tools are
 * bwa_bin: full path name to bwa binary
 * FORCE: force-regenerate all output files
 */
void AlignReadsToReference( const String &in_bam_file,
			    const String &base_dir,
			    const String &head_name,
			    const String &picard_dir,
			    const String &bwa_bin,
			    const bool FORCE = false );

/**
 * AlignContigsToReference
 *
 * Align a fastb object (contigs) onto a reference, and generate
 * corresponding bam file.
 *
 * in_fastb_file: full path name of input fastb
 * base_dir, head_name: define output dir (<base_dir>/<head_name>)
 * FORCE: force-regenerate all output files
 * FW_ONLY: only consider fw aligns
 */
void AlignContigsToReference( const String &in_fastb_file,
			      const String &base_dir,
			      const String &head_name,
			      const bool FORCE,
			      const bool FW_ONLY = false );

#endif
