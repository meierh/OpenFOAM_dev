#!/bin/bash

cd src/structure
rm -rf lnInclude Make/linux64GccDPInt32Opt
wmake -j 12 &

cd ../meshRefiner
rm -rf lnInclude Make/linux64GccDPInt32Opt
wmake -j 3 &

cd ../interaction
rm -rf lnInclude Make/linux64GccDPInt32Opt
wmake -j 12 &

cd ../optimizer
rm -rf lnInclude Make/linux64GccDPInt32Opt
wmake -j 3 &

cd ../icoImmersedBoundary
rm -rf lnInclude Make/linux64GccDPInt32Opt
wmake -j 4 &

cd ../icoAdjointImmersedBoundary
rm -rf lnInclude Make/linux64GccDPInt32Opt
wmake -j 8 &

cd ../../applications/solvers/icoIB
rm -rf Make/linux64GccDPInt32Opt
wmake -j 3 &

cd ../optimIB
rm -rf Make/linux64GccDPInt32Opt
wmake -j 3 &

cd ../fdIB
rm -rf Make/linux64GccDPInt32Opt
wmake -j 3 &

wait
cd ../../../

cd src/structure
wmake &

cd ../meshRefiner
wmake &

cd ../interaction
wmake &

cd ../optimizer
wmake &

cd ../icoImmersedBoundary
wmake &

cd ../icoAdjointImmersedBoundary
wmake &

cd ../../applications/solvers/icoIB
wmake &

cd ../optimIB
wmake &

cd ../fdIB
wmake &

wait
echo "Completed compilation"
