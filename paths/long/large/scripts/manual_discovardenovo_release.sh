#! /bin/bash

# Script to manually create a discovar release and post it to
# the FTP site under the 'experimental' subdirectory.
 
# INSTRUCTIONS
#
# To avoid potential environment issues, run this script as dexter with a clean environment
#
# First remove the directory /tmp/manual_release, if it exists.
#
# Then let R be the revision id that you want to push out.  Run this:
#
#      sudo -u dexter env -i paths/long/large/scripts/manual_discovardenovo_release.sh R
#
# To push out the latest version, omit R.
#
# There should be a tarball in
# /web/ftp/distribution/crd/nightly/discovardenovo/experimental
# Otherwise take a look at /tmp/manual_release/test_log.
#
# If it succeeds, then tell the user that they can download this revision from
# our experimental directory at:
# ftp://ftp.broadinstitute.org/pub/crd/nightly/discovardenovo/experimental/
#       discovardenovo-R.tar.gz

# Halt on any error
set -o errexit

if [[ $# -gt 2 ]]; then
  echo "Function: Creates a manual release of the discovar de novo package" 1>&2
  echo "Usage: `basename $0` [revision]" 1>&2
  exit 1
fi

username=$(whoami)
if [ "$username" != "dexter" ]; then
    echo -e "\nWARNING: Best run as dexter with a clean environment. E.g."
    echo -e "sudo -u dexter env -i $0 [revision] \n"
fi

tmp_dir="/tmp"
revision=""
if [ $# -gt 0 ]; then
    echo "Using SVN revision $1"
    revision=$1
else
    echo "Using latest SVN revision"
fi

base_dir=$tmp_dir/manual_release
package_dir=$base_dir/packaging
test_dir=$base_dir/testing

# setup working directory

if  [ "${base_dir:0:1}" != "/" ] ; then
    echo "You must supply an absolute path for the target directory"
    exit 1
fi

if [ -e $base_dir ]; then
  echo "$base_dir exists" 1>&2
  exit 1
fi

mkdir $base_dir
if [ $? -ne 0 ]; then
  echo "Error making directory $base_dir" 1>&2
  exit 1
fi

pushd $base_dir

echo -e "\n--------------------------------------------------"
echo -e "Retrieving packaging and test scripts"
echo -e "--------------------------------------------------\n"

# get the latest version of the general manual release code
svn export https://svn.broadinstitute.org/comprd/trunk/packard/manual_release.sh

# get the latest version of the dexter discovar package testing code
svn export https://svn.broadinstitute.org/comprd/trunk/dexter/scripts/discovardenovo_package_test.sh

echo -e "\n--------------------------------------------------"
echo -e "Running packaging script"
echo -e "--------------------------------------------------\n"

# package up code
./manual_release.sh $package_dir discovardenovo $revision &> package_log

# find the package name
package_name=$(find $package_dir -maxdepth 1 -name "discovardenovo-*.tar.gz" -printf '%P\n')

echo -e "\n--------------------------------------------------"
echo -e "Running test script"
echo -e "--------------------------------------------------\n"

# test the package
./discovardenovo_package_test.sh $test_dir $package_dir/$package_name &> test_log 

echo -e "\n--------------------------------------------------"
echo -e "Pushing to FTP site"
echo -e "--------------------------------------------------\n"

cp $package_dir/$package_name /web/ftp/pub/crd/nightly/discovardenovo/experimental

echo "$package_name released."

echo -e "\n--------------------------------------------------"
echo -e "Cleanup"
echo -e "--------------------------------------------------\n"

echo -e "\nRelease may be found at:\n"
echo -e "ftp://ftp.broadinstitute.org/pub/crd/nightly/discovardenovo/experimental"
echo -e "/discovardenovo-$revision.tar.gz\n"

rm -rf $base_dir
