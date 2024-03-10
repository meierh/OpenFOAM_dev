#ifndef LINESTRUCTURE_H
#define LINESTRUCTURE_H
#include "Structure.H"

namespace Foam
{
    
class LineStructure : public Structure
{
    public:
        explicit LineStructure
        (
            Time& runTime,
            const dimensionedScalar& alpha,
            const volScalarField& T,
            const volScalarField& p,
            const volVectorField& U,
            dynamicRefineFvMesh& mesh,
            const dimensionedScalar nu
        );
        LineStructure(const LineStructure&) = delete;
                
    private:
     
};
}
#endif
