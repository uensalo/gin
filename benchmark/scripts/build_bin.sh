#!/bin/bash
cd ../bin/
rm -rf *
cmake -DCMAKE_BUILD_TYPE=Release ../..
make
cd ../bin-sdsl/
cmake -DCMAKE_BUILD_TYPE=Release ../.. -DBUILD_SDSL=ON
make
cd ../scripts
cd ../bin-fast/
cmake -DCMAKE_BUILD_TYPE=Release ../.. -DBUILD_DNA_FMI=ON
make
cd ../scripts
