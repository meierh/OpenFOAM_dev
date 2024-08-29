#include "MeshRefiner.H"

Foam::MeshRefiner::MeshRefiner
(
    volScalarField& doRefine,
    const LineStructure& structure,
    fvMesh& mesh,
    Time& runTime
):
mesh(mesh),
runTime(runTime),
structure(structure),
doRefine(doRefine),
fieldRefineDemands("fieldDemands",doRefine),
markerRefineDemands("markerDemands",doRefine)
{}

void Foam::MeshRefiner::adaptMesh()
{
    fieldRefinement();
    markerRefinement(UNREFINE);

    for(label cellInd=0; cellInd<fieldRefineDemands.size(); cellInd++)
    {
        scalar fieldRefValue = fieldRefineDemands[cellInd];
        const scalar markerRefValue = markerRefineDemands[cellInd];
        scalar& refineValue = doRefine[cellInd];
        refineValue = refinementDemandMerge(fieldRefValue,markerRefValue);
    }

    bool refined = mesh.update();
    while(refined)
    {
        markerRefinement(MUSTKEEP);
        for(label cellInd=0; cellInd<markerRefineDemands.size(); cellInd++)
        {
            doRefine[cellInd] = markerRefineDemands[cellInd];
        }
        refined = mesh.update();
    }
}

void Foam::MeshRefiner::fieldRefinement()
{
    dimensionSet dimensions = fieldRefineDemands.dimensions();
    Foam::dimensioned<Foam::scalar> val("fieldRefinement",dimensions,-0.25);
    fieldRefineDemands = val;
}

void Foam::MeshRefiner::markerRefinement(scalar defaultValue)
{
    dimensionSet dimensions = fieldRefineDemands.dimensions();
    Foam::dimensioned<Foam::scalar> val("fieldRefinement",dimensions,defaultValue);
    markerRefineDemands = val;

    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    for(uint markerInd=0; markerInd<markers.size(); markerInd++)
    {
        const LagrangianMarker* oneMarker = markers[markerInd];
        label markerCellInd = oneMarker->getMarkerCell();
        scalar charLen = oneMarker->getMarkerCharacLen();
        scalar markerCellLen = Structure::initialSpacingFromMesh(mesh,markerCellInd);
        if(charLen<markerCellLen)
            markerRefineDemands[markerCellInd] = 1;
        else
            markerRefineDemands[markerCellInd] = 0;
    }
}

Foam::scalar Foam::MeshRefiner::refinementDemandMerge
(
    Foam::scalar fieldDemand,
    Foam::scalar markerDemand
)
{
    scalar refineValue;

    if(fieldDemand>REFINE_LIM) // refine
    {
        refineValue = fieldDemand;
    }
    else if(fieldDemand>MUSTKEEP_LIM) // mustkeep
    {
        if(markerDemand>REFINE_LIM)
            refineValue = markerDemand;
        else
            refineValue = fieldDemand;
    }
    else if(fieldDemand>DONTCARE_LIM) // dontcare
    {
        refineValue = markerDemand;
    }
    else // unrefine
    {
        if(markerDemand>MUSTKEEP_LIM)
            refineValue = markerDemand;
        else
            refineValue = fieldDemand;
    }

    return refineValue;
}

