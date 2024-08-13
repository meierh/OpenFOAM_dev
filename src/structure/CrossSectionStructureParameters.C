#include "CrossSectionStructure.H"

Foam::CrossSectionCoeffReference::CrossSectionCoeffReference
(
    label rodNumber,
    label fourierCoeffNumber,
    label coeffNumber
):
rodNumber(rodNumber),
fourierCoeffNumber(fourierCoeffNumber),
coeffNumber(coeffNumber)
{}

Foam::CrossSectionStructureParameter::CrossSectionStructureParameter
(
    CrossSectionCoeffReference coeffRef
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
    for(label rodNumber=0; rodNumber<structure->getNumberRods(); rodNumber++)
    {
        const CrossSection& crossSec = structure->getRodCrossSections()[rodNumber];
        label nbrFourierCoeffs = crossSec.numberFourierCoeff();
        for(label fourCoeffInd=0; fourCoeffInd<nbrFourierCoeffs; fourCoeffInd++)
        {
            label nbrNurbsCoeffs = crossSec.numberNurbsCoeffs(fourCoeffInd);
            for(label coeffInd=0; coeffInd<nbrNurbsCoeffs; coeffInd++)
            {
                CrossSectionStructureParameter para(CrossSectionCoeffReference(rodNumber,fourCoeffInd,coeffInd));
                parameters.push_back(para);
            }
        }
    }
}
