#include "MeshRefiner.H"

Foam::MeshRefiner::MeshRefiner
(
    volScalarField& doRefine,
    const LineStructure& structure,
    fvMesh& mesh,
    Time& runTime
):
doRefine(doRefine),
structure(structure),
mesh(mesh),
runTime(runTime),
markerRefinementInfo
(
    IOobject
    (
        "markerRefine",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
)
{
}

void Foam::MeshRefiner::adaptMesh()
{
    bool refined;
    do
    {
        setMarkerRefinementInfo();
        for(label cellInd=0; cellInd<markerRefinementInfo.size(); cellInd++)
        {
            doRefine[cellInd] = markerRefinementInfo[cellInd];
        }
        refined = mesh.update();
    }
    while(refined);
}

void Foam::MeshRefiner::setMarkerRefinementInfo()
{
    markerRefinementInfo = -1;
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    for(uint markerInd=0; markerInd<markers.size(); markerInd++)
    {
        const LagrangianMarker* oneMarker = markers[markerInd];       
        label markerCellInd = oneMarker->getMarkerCell();
        scalar charLen = oneMarker->getMarkerCharacLen();
        const cell& markerCell = mesh.cells()[markerCellInd];
        scalar markerCellLen = Structure::initialSpacingFromMesh(mesh,markerCellInd);
        if(charLen<markerCellLen)
            markerRefinementInfo[markerCellInd] = 1;
        else
            markerRefinementInfo[markerCellInd] = 0;
    }
}
