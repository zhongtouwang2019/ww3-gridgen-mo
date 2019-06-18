#!/bin/bash

set -eu

WRKDIR=$1
MYDIR=$PWD

cp genSides.exe $WRKDIR
cp countijsdnew $WRKDIR

cd $WRKDIR
echo 'Launching genSides.exe to generate face arrays'
#genSides.exe >>genSides_output.txt 2>&1
genSides.exe >genSides_output.txt

echo 'Sorting the face arrays'
countijsdnew ww3GISide.d ww3GJSide.d

#echo 'Tidying up'
#rm genSides.exe
#rm countijsdnew
cd $MYDIR

exit 0
