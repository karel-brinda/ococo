#! /usr/bin/env bash

make clean
find . | grep CMakeCache.txt | xargs rm -f

cmake . -DCMAKE_BUILD_TYPE=Release -DOCOCO32=OFF
make -j 3 || make
mv ococo ococo16

make clean

cmake . -DCMAKE_BUILD_TYPE=Release -DOCOCO32=ON
make -j 3 || make
mv ococo ococo32

