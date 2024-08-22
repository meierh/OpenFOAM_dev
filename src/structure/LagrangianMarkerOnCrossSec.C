#include "CrossSectionStructure.H"

Foam::LagrangianMarkerOnCrossSec::LagrangianMarkerOnCrossSec
(
    const Structure& structure,
    const fvMesh& mesh,
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* baseRod,
    const scalar markerParameter,
    const CrossSection* baseCrossSec,
    const scalar markerAngle,
    const scalar radiusFrac
):
LagrangianMarker(structure,mesh,rodNumber,baseRod),
markerAngle(markerAngle),
markerRadiusFrac(radiusFrac),
baseCrossSec(baseCrossSec)
{
    this->markerParameter = markerParameter;
    evaluateMarker();
}

void Foam::LagrangianMarkerOnCrossSec::evaluateMarker()
{
    markerPosition = CrossSectionStructure::evaluateRodCircumPos
    (
        baseRod,markerParameter,baseCrossSec,markerAngle,markerRadiusFrac
    );
    markerCell = mesh.findCell(markerPosition);
    computeSupport();
    Pair<vector> h = minMaxNeighbourWidth(directSupport);
    h_plus = h.first();
    h_minus = h.second();
    dilation = dilationFactors(h);
    checkDirectSupport();
    reduceSupport();
}
