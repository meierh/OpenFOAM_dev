#ifndef functionObjects_FSIFluidBoundaryConditions_H
#define functionObjects_FSIFluidBoundaryConditions_H

#include "fvMeshFunctionObject.H"

namespace Foam
{
namespace functionObjects
{
class FSIFluidBoundaryConditions : public fvMeshFunctionObject
{
    const word& name,
    const Time& runTime,
    const dictionary& dict
};
}
}

#endif
