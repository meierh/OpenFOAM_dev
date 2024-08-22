#include "MeshRefiner.H"

Foam::MeshRefiner::MeshRefiner
(
    volScalarField& doRefine,
    const LineStructure& structure
):
doRefine(doRefine),
structure(structure)
{
}

void MeshRefiner::setRefinementCells()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    for(label markerInd=0; markerInd<markers.size(); markerInd++)
    {
        const LagrangianMarker* oneMarker = markers[markerInd];       
        label markerCellInd = oneMarker->getMarkerCell();
        scalar charLen = oneMarker->getMarkerCharacLen();
        const cell& markerCell = fvMesh.cells()[markerCellInd];
        scalar markerCellLen = Structure::initialSpacingFromMesh(mesh,markerCellInd);
    }
}
