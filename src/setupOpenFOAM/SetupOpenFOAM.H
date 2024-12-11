#ifndef SETUPOPENFOAM_H
#define SETUPOPENFOAM_H
#include "argList.H"
#include "pimpleSingleRegionControl.H"
#include <memory>

namespace Foam
{

std::unique_ptr<Time> createTime(const argList& args);
std::unique_ptr<fvMesh> createMesh(Time& runTime);

}
#endif
