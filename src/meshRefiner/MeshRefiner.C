#include "MeshRefiner.H"

Foam::MeshRefiner::MeshRefiner
(
    fvMesh& mesh,
    LineStructure& structure,
    volScalarField& doRefine,
    dictionary& dynamicMeshDict
):
mesh(mesh),
structure(structure),
doRefine(doRefine),
dynamicMeshDict(dynamicMeshDict),
fieldRefineDemands("fieldDemands",doRefine),
markerRefineDemands("markerDemands",doRefine)
{
    dictionary& topoChangerDict = dynamicMeshDict.subDict("topoChanger");
    if(topoChangerDict.found("markerCellFactor"))
    {
        //Set refine/stay/unrefine field values
        if(topoChangerDict.found("upperRefineLevel")) topoChangerDict.set("upperRefineLevel",1.5);
        else topoChangerDict.add("upperRefineLevel",1.5);
        
        if(topoChangerDict.found("lowerRefineLevel")) topoChangerDict.set("lowerRefineLevel",0.5);
        else topoChangerDict.add("lowerRefineLevel",0.5);
        
        if(topoChangerDict.found("unrefineLevel")) topoChangerDict.set("unrefineLevel",-0.5);
        else topoChangerDict.add("unrefineLevel",-0.5);
        
        
        ITstream topoChangerFactorStream = topoChangerDict.lookup("markerCellFactor");
        token topoChangerFactorToken;
        topoChangerFactorStream.read(topoChangerFactorToken);
        if(!topoChangerFactorToken.isScalar())
        {
            Info<<"topoChangerFactorToken:"<<topoChangerFactorToken<<Foam::endl;
            FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChanger/markerCellFactor -- must be scalar"<<exit(FatalError);
        }
        scalar topoChangerFactorScalar = topoChangerFactorToken.scalarToken();
        if(topoChangerFactorScalar<=0)
            FatalErrorInFunction<<"Invalid topoChanger/markerCellFactor is"<<topoChangerFactorScalar<<" -- valid {]0,inf[}"<<exit(FatalError);
        markerCharLengthToCellSizeFactor = topoChangerFactorScalar;
    }

}

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
    mesh.update();
    
    refineMeshAndMarkers();
}

void Foam::MeshRefiner::refineMeshOnStaticMarkers()
{
    bool refined = true;
    while(refined)
    {
        Info<<"------------Refine mesh iteration-------------"<<Foam::endl;
        markerRefinement(MUSTKEEP);
        for(label cellInd=0; cellInd<markerRefineDemands.size(); cellInd++)
        {
            doRefine[cellInd] = markerRefineDemands[cellInd];
        }
        refined = mesh.update();
    }
}

void Foam::MeshRefiner::refineMeshAndMarkers()
{
    bool refined = true;
    while(refined)
    {
        Info<<"------------Refine mesh and marker iteration-------------"<<Foam::endl;
        markerRefinement(MUSTKEEP);
        for(label cellInd=0; cellInd<markerRefineDemands.size(); cellInd++)
        {
            doRefine[cellInd] = markerRefineDemands[cellInd];
        }
        refined = mesh.update();
        structure.settleIntoRefinedMesh();
        structure.refineMarkersOnRefinedMesh();
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

