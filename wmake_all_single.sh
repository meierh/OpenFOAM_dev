#!/bin/bash

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
cd ../setupOpenFOAM
wmake &

cd ../../applications/solvers/icoIB
wmake &
cd ../optimIB
wmake &
cd ../fdIB
wmake &

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
cd ../setupOpenFOAM
wmake &

cd ../../applications/solvers/icoIB
wmake &
cd ../optimIB
wmake &
cd ../fdIB
wmake &

wait
echo "Completed compilation"
