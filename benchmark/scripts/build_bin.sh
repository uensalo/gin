#!/bin/bash
cd ../bin/
rm -rf *
cmake -DCMAKE_BUILD_TYPE=Release ../..
make
cd ../scripts