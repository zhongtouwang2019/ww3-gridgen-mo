#!/bin/bash

set -eu

WRKDIR=$1
MYDIR=$PWD

cp smcProps.exe $WRKDIR

cd $WRKDIR
echo 'Launching smcProps.exe for propagation tests'
smcProps.exe >smcProps_output.txt

if [ ! -d smcProps ]; then
  mkdir smcProps
fi
mv Cn*.d smcProps
cd smcProps
echo 'Propagation test data file list' > ../cfiles.txt
ls Cn*.d >> ../cfiles.txt

cd $MYDIR

exit 0
