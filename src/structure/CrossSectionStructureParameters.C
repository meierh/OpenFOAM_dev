#include "CrossSectionStructure.H"

Foam::CrossSectionStructureParameter::CrossSectionStructureParameter
(
    CrossSectionCoeffReference coeff
)
{
    coeffs.push_back(coeffRef);
}


Foam::CrossSectionStructureParameter::CrossSectionStructureParameter
(
    const std::vector<CrossSectionCoeffReference>& coeffs
):
coeffs(coeffs)
{}

void Foam::CrossSectionStructureParameter::addCoeff
(
    CrossSectionCoeffReference coeffRef
)
{
    coeffs.push_back(coeffRef);
}

void Foam::CrossSectionStructureParameters::collectParameters
(
    const CrossSectionStructure* structure
)
{
    LineStructureParameters::collectParameters(structure);
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        const CrossSection& crossSec = structure->getRodCrossSections()[rodNumber];
        label nbrFourierCoeffs = crossSec.numberFourierCoeff();
        for(label fourCoeffInd=0; fourCoeffInd<nbrFourierCoeffs; fourCoeffInd++)
        {
            label nbrNurbsCoeffs = crossSec.numberNurbsCoeffs(fourCoeffInd);
            for(label coeffInd=0; coeffInd<nbrNurbsCoeffs; coeffInd++)
            {
                parameters.push_back(CrossSectionStructureParameter{rodNumber,fourCoeffInd,coeffInd});
            }
        }
    }
}
