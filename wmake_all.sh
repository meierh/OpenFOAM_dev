#!/bin/bash

cd src/structure
wmake -j 12 &

cd ../meshRefiner
wmake -j 3 &

cd ../interaction
wmake -j 12 &

cd ../optimizer
wmake -j 3 &

cd ../icoImmersedBoundary
wmake -j 4 &

cd ../icoAdjointImmersedBoundary
wmake -j 8 &

cd ../../applications/solvers/icoIB
wmake -j 3 &

cd ../optimIB
wmake -j 3 &

cd ../fdIB
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
