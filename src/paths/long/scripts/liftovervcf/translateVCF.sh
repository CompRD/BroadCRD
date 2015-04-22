#!/usr/bin/tcsh
###############################################################################
##                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     ##
##       This software and its documentation are copyright (2013) by the     ##
##   Broad Institute.  All rights are reserved.  This software is supplied   ##
##   without any warranty or guaranteed support whatsoever. The Broad        ##
##   Institute is not responsible for its use, misuse, or functionality.     ##
###############################################################################

# Script to translate VCF files computed using hg18 to hg19 co-ordinates
#
# Usage: translateVCF.sh filename
#
# Where filename is the filename of the vcf file to convert, minus the .vcf extension.
#
# Requires:
#   GATK toolkit
#   GATK chain files, available from:
#     /humgen/gsa-hpprojects/GATK/data/Liftover_Chain_Files/
#   GATK liftOverVCF.pl script (modified - included here)
#   GATK sortByRef.pl (included here)
#   Hg18 fasta, dict and index files
#   Hg19 fasta, dict and index files
#
# The liftOverVCF.pl from GATK doesn't work with their latest toolkit. The version
# included here has been modified to work.
#
# GATK is rather fussy about the format of the VCF file:
#   The columns must be tab separated
#   The chromosomes must be in the same order as those in the reference
#   The chromosome names must match
# This script ensures this by first preparing the vcf file.
#
# The script requires that the GATK jar and perl scripts are in the same directory
# You will need to alter the paths to the reference and chain files below in the
# call to liftOverVCF.pl
# This script does not clean up after itself.

set vcf_in = $1

# Reorder the VCF by chromosome number - numeric not alpahabetic order
grep '^#' $vcf_in.vcf > ordered.vcf
grep -v -E '^X|^Y|^#' $vcf_in.vcf | sort -n -k1 -k2 >> ordered.vcf
grep -E '^X|^Y' $vcf_in.vcf | sort -k1,1d -k2,2n >> ordered.vcf

# Replace chromosome numbers with names in the form ChrN
cat ordered.vcf | sed 's/^[0-9]/chr&/' | sed 's/^X/chr&/' > renamed.vcf

# Replace spaces with tabs
cat renamed.vcf | tr --squeeze-repeat " " "\t" > fixed.vcf

# Run GATK conversion script
liftOverVCF.pl -vcf fixed.vcf -chain hg18ToHg19.broad.over.chain -out $vcf_in.hg19.vcf -gatk . -newRef hg19/Homo_sapiens_assembly19 -oldRef hg18/Homo_sapiens_assembly18 >& $vcf_in.log

