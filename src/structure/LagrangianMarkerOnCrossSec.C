#include "CrossSectionStructure.H"

Foam::LagrangianMarkerOnCrossSec::LagrangianMarkerOnCrossSec
(
    const Structure& structure,
    const dynamicRefineFvMesh& mesh,
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* baseRod,
    const scalar markerParameter,
    const CrossSection* baseCrossSec,
    const scalar markerAngle,
    const scalar radiusFrac
):
structure(structure),
mesh(mesh),
rodNumber(rodNumber),
baseRod(baseRod),
markerParameter(markerParameter),
baseCrossSec(baseCrossSec),
markerAngle(markerAngle),
markerRadiusFrac(radiusFrac)
{
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
    minMaxSupportWidth();
    dilationFactors();
}
