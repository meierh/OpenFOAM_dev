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

cd ../../applications/solvers/icoIB
wmake &

cd ../optimIB
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

cd ../../applications/solvers/icoIB
wmake &

cd ../optimIB
wmake &

wait
echo "Completed compilation"
