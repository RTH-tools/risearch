#!/bin/bash

WD=`pwd`

cd RIsearch1/src || exit 1
mkdir -p ../bin
make clean || exit 1
make || exit 1
cd $WD

cd RIsearch2
cd src || exit 1
mkdir -p ../bin
make clean || exit 1
make || exit 1
cd $WD

echo
echo '-------------------------------------------------'
echo 'You should now be able to run RIsearch1 and RIsearch2 by typing'
echo '-------------------------------------------------'

echo "$PWD"/RIsearch1/bin/RISearch
echo "$PWD"/RIsearch2/bin/risearch2.x


