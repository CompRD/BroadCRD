===============================================
How to convert BAM files to ALLPATHS-ready data
===============================================


AllPaths requires two types of read libraries as input: a fragment
library and a jumping library.  For each library quality scores and
pairing information is also necessary.  Examples would be:

  <DATA_DIR>/frag_reads_orig.fastb
  <DATA_DIR>/frag_reads_orig.qualb
  <DATA_DIR>/frag_reads_orig.pairs

  <DATA_DIR>/jump_reads_orig.fastb
  <DATA_DIR>/jump_reads_orig.qualb
  <DATA_DIR>/jump_reads_orig.pairs


These files can be automatically generated from BAM files by running
the script 'BAMsToAllPathsInputs.pl'.  This script takes as input two
comma-separated-values (.csv) files

  ./in_groups.csv
  ./in_libs.csv

which describe the locations and library information of the various
BAM files available.


-------------
in_groups.csv
-------------

Each line in './in_groups.csv' provides, for each BAM file, the
following information:


  group_name:    a UNIQUE nickname for this specific BAM file.

  library_name:  the library to which the BAM file belongs.

  file_name:     the absolute path to the BAM file.



Example (NOTE: column order is not important):


                     file_name, library_name, group_name
/picard/Solexa-11541/302GJ.bam, Solexa-11541,      302GJ
/picard/Solexa-11542/303GJ.bam, Solexa-11542,      303GJ



-----------
in_libs.csv
-----------

Each line in './in_libs.csv' describes a specific library.  The
specific fields are:


  library_name:     the same as in 'in_groups.csv'.

  project_name:     a string naming the project.

  organism_name:    the organism.

  type:             'fragment', 'jumping', 'EcoP15', etc.

  paired:           '1': paired reads.
                    '0': UNpaired reads.

  frag_size:        average number of bases in the fragments.

  frag_stddev:      estimated standard deviation of the fragments
                    sizes.

  insert_size:      average number of bases in the inserts. 
                    NOTE: if this field is empty this is a FRAGMENT
                          library, otherwise it is a JUMPING library.

  insert_stddev:    estimated standard deviation of the inserts sizes.
                    NOTE: if this field is empty this is a FRAGMENT
                          library, otherwise it is a JUMPING library.

  read_orientation: 'inward' or 'outward'.
                    NOTE: 'outward' reads will be reverted.

  genomic_start:    index of the FIRST genomic base in the reads.
                    NOTE: if non-zero, reads will be trimmed.

  genomic_end:      index of the LAST genomic base in the reads.
                      NOTE: if non-zero, reads will be trimmed.



Example (NOTE: all the fields should be on a single line; that makes
the lines too long to show here, hence the '...'):


  library_name, project_name, organism_name,     type, paired, ...
  Solexa-11541,      Awesome,        E.coli, fragment,      1, ...
  Solexa-11542,      Awesome,        E.coli,  jumping,      1, ...


       ... frag_size, frag_stddev, insert_size, insert_stddev, ... 
       ...       180,          10,            ,              , ...
       ...       300,          20,        3000,           500, ...


                  ... read_orientation, genomic_start, genomic_end
                  ...           inward,             0,         100
                  ...          outward,             0,         100 


-----------------------
BAMsToAllPathsInputs.pl
-----------------------

Simplest example of 'BAMsToAllPathsInputs.pl':

  
  BAMsToAllPathsInputs.pl \
    DATA_DIR=./projects/EColi/data\
    PICARD_TOOLS_DIR=/picard/bin


where DATA_DIR is the ALLPATHS DATA directory (see the AllPaths manual
for details) where the input reads will be placed, and PICARD_TOOLS_DIR 
is the path to the Picard tools (http://picard.sourceforge.net/).
'BAMsToAllPathsInputs.pl' uses the Picard tools to process BAM files.

There are a few other options that can be specified:


  IN_GROUPS_CSV          use a file other than ./in_groups.csv.
 
  IN_LIBS_CSV            use a file other than ./in_libs.csv.

  GENOME_SIZE            estimated genome size for the purpose of
                         coverage estimation.

  FRAG_COVERAGE          fragment library desired coverage, e.g. 45

  JUMP_COVERAGE          jumping library desired coverage, e.g. 45

  INCLUDE_NON_PF_READS   '1':(default) include non-PF reads. 
                         '0': include only PF reads.

  HOSTS                  e.g. '2,3.host2,4.host3'
                         when converting several BAM files, fork:
                           2 processes on the localhost;
                           3 processes on host2;
                           4 processes on host3;
                         NOTE: forking to remote hosts requires
                         password-less ssh access (using ssh-agent /
                         ssh-add for example). 
















