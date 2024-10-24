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
    
    if(topoChangerDict.found("field")) topoChangerDict.set("field",doRefine.name());
    else topoChangerDict.add("field",doRefine.name());
        
    //Set refine/stay/unrefine field values
    if(topoChangerDict.found("upperRefineLevel")) topoChangerDict.set("upperRefineLevel",1.5);
    else topoChangerDict.add("upperRefineLevel",1.5);
    
    if(topoChangerDict.found("lowerRefineLevel")) topoChangerDict.set("lowerRefineLevel",0.5);
    else topoChangerDict.add("lowerRefineLevel",0.5);
    
    if(topoChangerDict.found("unrefineLevel")) topoChangerDict.set("unrefineLevel",-0.5);
    else topoChangerDict.add("unrefineLevel",-0.5);
    
    if(topoChangerDict.found("markerCellFactor"))
    {
        ITstream topoChangerFactorStream = topoChangerDict.lookup("markerCellFactor");
        token topoChangerFactorToken;
        topoChangerFactorStream.read(topoChangerFactorToken);
        if(!topoChangerFactorToken.isScalar())
        {
            Info<<"topoChangerFactorToken:"<<topoChangerFactorToken<<Foam::nl;
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
    FatalErrorInFunction<<"Not in use anymore"<<exit(FatalError);
    /*
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
    */
}

bool Foam::MeshRefiner::refineMeshOnStaticMarkers()
{   
    bool meshWasRefined = false;
    bool refined = true;
    while(refined)
    {
        Info<<"------------refineMeshOnStaticMarkers-------------"<<Foam::nl;
        auto start = std::chrono::system_clock::now();
        markerRefinement(MUSTKEEP);
        for(label cellInd=0; cellInd<markerRefineDemands.size(); cellInd++)
        {
            doRefine[cellInd] = markerRefineDemands[cellInd];
        }
        refined = applyMeshAdaption();
        meshWasRefined |= refined;
        auto end = std::chrono::system_clock::now();
        Info<<"------------------------------------------------- took:"<<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<" milliseconds"<<Foam::nl;
    }    
    return meshWasRefined;
}

bool Foam::MeshRefiner::refineMeshAndMarkers
(
    bool preRefinedMesh
)
{    
    if(preRefinedMesh)
        structure.refineMarkersOnRefinedMesh();
            
    bool meshWasRefined = false;
    bool refined = true;
    while(refined)
    {
        Info<<"------------refineMeshAndMarkers-------------"<<Foam::nl;
        auto start = std::chrono::system_clock::now();
        markerRefinement(MUSTKEEP);
        for(label cellInd=0; cellInd<markerRefineDemands.size(); cellInd++)
        {
            doRefine[cellInd] = markerRefineDemands[cellInd];
        }
        refined = applyMeshAdaption();
        Info<<"Refined:"<<refined<<Foam::endl;
        meshWasRefined |= refined;
        if(refined)
            structure.refineMarkersOnRefinedMesh();
        auto end = std::chrono::system_clock::now();
        Info<<"------------------------------------------------- took:"<<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<" milliseconds"<<Foam::nl;
    }
        
    return meshWasRefined;
}

void Foam::MeshRefiner::fieldRefinement()
{
    dimensionSet dimensions = fieldRefineDemands.dimensions();
    Foam::dimensioned<Foam::scalar> val("fieldRefinement",dimensions,-0.25);
    fieldRefineDemands = val;
}

void Foam::MeshRefiner::markerRefinement
(
    scalar defaultValue
)
{    
    dimensionSet dimensions = fieldRefineDemands.dimensions();
    Foam::dimensioned<Foam::scalar> val("fieldRefinement",dimensions,defaultValue);
    markerRefineDemands = val;
       
    
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    for(uint markerInd=0; markerInd<markers.size(); markerInd++)
    {
        const LagrangianMarker* oneMarker = markers[markerInd];
        label markerCellInd = oneMarker->getMarkerCell();
        if(markerCellInd!=-1)
        {
            scalar charLen = oneMarker->getMarkerCharacLen();
            
            DynamicList<label> neighbours;
            Structure::neighbourCells(mesh,markerCellInd,neighbours);
            neighbours.append(markerCellInd);
            for(label cell : neighbours)
            {
                scalar cellLen = Structure::spacingFromMesh(mesh,cell);
                if(charLen*markerCharLengthToCellSizeFactor < cellLen)
                {
                    Info<<"Refine because:"<<charLen<<"*"<<markerCharLengthToCellSizeFactor<<" < "<<cellLen<<Foam::endl;
                    markerRefineDemands[cell] = std::max(REFINE,markerRefineDemands[cell]);
                }
                else
                    markerRefineDemands[cell] = std::max(MUSTKEEP,markerRefineDemands[cell]);
            }
        }
    }    
}

bool Foam::MeshRefiner::applyMeshAdaption()
{
    bool refined = mesh.update();
    if(refined)
    {
        structure.computeHaloData();
        structure.settleIntoRefinedMesh();
    }
    return refined;
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

