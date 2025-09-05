#!/bin/bash

WD=`pwd`

aclocal
autoheader
autoconf
automake --add-missing --copy --force

cd RIsearch2
rm -rf libdivsufsort-2.0.1/build || exit 1
mkdir -p libdivsufsort-2.0.1/build || exit 1
cd libdivsufsort-2.0.1/build || exit 1
cmake -DCMAKE_BUILD_TYPE="Release" -DBUILD_DIVSUFSORT64:BOOL=ON -DUSE_OPENMP:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=OFF .. || exit 1
make || exit 1
cd - || exit 1
cd $WD

echo
echo '-------------------------------------------------'
echo 'You should now be able to run'
echo './configure && make && make install'
echo '-------------------------------------------------'

echo "$PWD"/RIsearch1/bin/RISearch
echo "$PWD"/RIsearch2/bin/risearch2.x



