#!/bin/sh


RELEASE_DIR=Arachne

cvs -q checkout -r release-2_0_1 -P -d ArachneSrc Arachne >> /dev/null

# create the subdirs in Arachne
DIRS="Arachne_src/CVS Arachne_src/doc"
for DIR in $DIRS ; do
	echo "$DIR"
	mkdir -p "$DIR"
done

# copy the files that we're releasing
cd ArachneSrc
cp `make find_dependencies | grep -v xerces_include | grep -v '\.o'` ../Arachne_src

# there are some files that are not included in the above list, so we'll do
# them manually
FILES="Makefile Makefile_g++ MakeDepend.cc \
		ViewAlignments.cc libxerces-c1_5_1.so Time doc/Manual.html \
		doc/ReleaseNotes.html \
        LICENSE.txt LICENSE.xerces.txt CVS/Entries xerces_include"

for FILE in $FILES ; do
    cp -fR "$FILE" "../Arachne_src/$FILE"
done

cd ..
rm -rf ArachneSrc

# remove CVS directories from xerces_include and its subdirs
(cd Arachne_src ; rm -rf `find xerces_include -name CVS -print` ; cd ..)

# package up the source for distribution
(tar vcf Arachne_src.tar Arachne_src; gzip -vf Arachne_src.tar)

# compile Arachne to make sure everything kosher
(cd Arachne_src ; make DOWNLOAD_DIR=.. bin_download data_download )
