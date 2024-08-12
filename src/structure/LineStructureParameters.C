#include "LineStructure.H"

Foam::LineStructureParameter::LineStructureParameter
():
dimension(-1)
{}

Foam::LineStructureParameter::LineStructureParameter
(
    NurbsCoeffReference coeff
):
dimension(-1)
{
    dimension = coeff.dimension;
    coeffs.push_back(coeffRef);
}


Foam::LineStructureParameter::LineStructureParameter
(
    const std::vector<NurbsCoeffReference>& coeffs
):
coeffs(coeffs),
dimension(-1)
{
    for(const NurbsCoeffReference& coeffRef : coeffs)
    {
        if(dimension==-1)
            dimension = coeffRef.dimension;
        else if(dimension!=coeffRef.dimension)
            FatalErrorInFunction<<"Coefficients of different dimensions can not be one parameter!"<<exit(FatalError);
    }
}

void Foam::LineStructureParameter::addCoeff
(
    NurbsCoeffReference coeffRef
)
{
    if(dimension==-1)
        dimension = coeffRef.dimension;
    else if(dimension!=coeffRef.dimension)
        FatalErrorInFunction<<"Coefficients of different dimensions can not be one parameter!"<<exit(FatalError);
    coeffs.push_back(coeffRef);
}

void Foam::LineStructureParameters::collectParameters
(
    const LineStructure* structure
)
{
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        label nbrCoeffs = structure->numberCoeffs(rodNumber);
        for(label coeffInd=0; coeffInd<nbrCoeffs; coeffInd++)
        {
            std::array<LineStructureParameter,3> thriDimPara;
            for(dim=0; dim<3; dim++)
            {
                LineStructureParameter para(NurbsCoeffReference{rodNumber,coeffInd,dim});
                thriDimPara[dim] = para;
            }
            threeDimParameters.push_back(thriDimPara);
        }
    }
}
