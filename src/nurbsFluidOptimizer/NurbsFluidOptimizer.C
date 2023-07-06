#include "NurbsFluidOptimizer.H"

Foam::NurbsFluidOptimizer::NurbsFluidOptimizer
(
    Foam::argList args
)
{
    fileName runDirectory = args.rootPath();
    fileName caseName = args.caseName();
    NurbsReader Reader(runDirectory,caseName);
    nurbsGeom = Reader.getNurbsData();
}

std::shared_ptr<std::vector<NurbsData>> Foam::NurbsFluidOptimizer::getNurbsParameters
(
)
{
    return nurbsGeom;
}
