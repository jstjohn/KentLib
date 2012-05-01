#!/bin/bash

KENTSRC=$HOME/kent/src

pushd inc
for F in `ls *.h`
do
    cp ${KENTSRC}/inc/${F} ${F}
done

popd


pushd lib
for F in `ls *.c`
do
    cp ${KENTSRC}/lib/${F} ${F}
done
popd

make clean
make CC=gcc44
