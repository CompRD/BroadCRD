#!/bin/csh

#Download and assemble the jumping library.

  set current=`pwd`
  #/prodinfo/prodapps/arachneExtractor/ArachneExtractor D322 all palvarez
  mkdir -p qual
  mkdir -p fasta
  mkdir -p traceinfo
  mv *fasta.* fasta
  mv *qual.* qual
  mv *xml* traceinfo
  rm -f /wga/dev1/WGAdata/projects/jumping
  ln -s $current /wga/dev1/WGAdata/projects/jumping 
  cp ../D322b/genome.size .
  cp ../D322b/nhaplotypes .
  cp ../D322b/reads_config.xml .
  cp ../D322b/organism .
  cp ../D322b/vector.fasta .
  cp /wga/dev1/WGAdata/dtds/configuration.dtd traceinfo
  echo  "WIBR	1000	F	GGTTTAAACGAATTCGCCCCT	0" > insert.sites
  echo  "WIBR	1000	R	GGTTTAAACGAATTCGCCCCT	0" >> insert.sites
  pushd ~palvarez/projects/Arachne
  Assemblez DATA=projects/jumping RUN=run1 STOP=PartitionInput
  popd
