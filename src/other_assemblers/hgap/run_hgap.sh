#!/bin/bash -l

echo "Setting up environment"
use .smrtanalysis-2.0.1-centos
source $SEYMOUR_HOME/etc/setup.sh

picard_dir=/wga/scr4/picard/pacbio
working_dir=/wga/scr4/iainm/hgap

echo "Preparing list of pacbio files to assemble"
cat $picard_dir/*/input.fofn > $working_dir/pacbio_h5_filelist

echo "Running fofnToSmrtpipeInput.py"
fofnToSmrtpipeInput.py $working_dir/pacbio_h5_filelist > $working_dir/input.xml

cd $working_dir

echo "Running hgap"
smrtpipe.py --params=settings.xml xml:input.xml
