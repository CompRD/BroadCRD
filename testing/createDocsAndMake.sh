#!/bin/tcsh

# We are using a specific version which I compiled.
# This turns out to be easier to maintain than if IT does it.
set path=( /wga/scr2/Binaries/ohsumit/binutils-2.17.50/bin /wga/scr2/Binaries/ohsumit/m4-1.4/bin /wga/scr2/Binaries/ohsumit/flex-2.5.33/bin /wga/scr2/Binaries/ohsumit/bison-2.3/bin /wga/scr2/Binaries/ohsumit/doxygen-1.5.4/bin $path )

# Check out a clean copy of Arachne.
cd /wga/dev/Nightly/doxygen
rm -rf *
cvs -q co Arachne > /dev/null
echo Source checked out of repository.

# Regenerate the doxygen docs.
#  The Doxygen configuration file is located in the root directory of Arachne
cd Arachne
( doxygen > /wga/dev/Nightly/doxygen/doxygen.log ) >& /wga/dev/Nightly/doxygen/doxygen.errs
echo Documentation rebuilt.

# Check for compilation errors.
make -k OPTIM=none > /dev/null
echo Compilation tested.

