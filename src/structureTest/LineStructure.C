#include "LineStructure.H"

Foam::LineStructure::LineStructure
(
    Time& runTime,
    const dimensionedScalar& alpha,
    const volScalarField& T,
    const volScalarField& p,
    const volVectorField& U,
    dynamicRefineFvMesh& mesh,
    const dimensionedScalar nu
):
Structure(runTime,alpha,T,p,U,mesh,nu)
{
}
