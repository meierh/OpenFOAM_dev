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
dynamicMeshDict(dynamicMeshDict),
fieldRefineDemands("fieldDemands",doRefine),
markerRefineDemands("markerDemands",doRefine),
doRefine(doRefine),
cellLevel("cellLevel",doRefine)
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
    
    if(topoChangerDict.found("fluidRefineOnIB"))
    {
        ITstream fluidRefineOnIBStream = topoChangerDict.lookup("fluidRefineOnIB");
        token fluidRefineOnIBToken;
        fluidRefineOnIBStream.read(fluidRefineOnIBToken);
        if(!fluidRefineOnIBToken.isWord())
        {
            Info<<"fluidRefineOnIBToken:"<<fluidRefineOnIBToken<<Foam::nl;
            FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChanger/fluidRefineOnIB -- must be word"<<exit(FatalError);
        }
        word fluidRefineOnIBWord = fluidRefineOnIBToken.wordToken();
        if(fluidRefineOnIBWord=="yes")
            fluidRefineOnIB = true;
        else if(fluidRefineOnIBWord=="no")
            fluidRefineOnIB = false;
        else
            FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChanger/fluidRefineOnIB -- must be {yes,no}"<<exit(FatalError);
    }
    
    const fvMeshTopoChanger& topoCh = mesh.topoChanger();
    const fvMeshTopoChanger* topoChPtr = &topoCh;
    refinerPtr = dynamic_cast<const fvMeshTopoChangers::refiner*>(topoChPtr);
}

bool Foam::MeshRefiner::refineMeshAndMarkers
(
    bool preRefinedMesh
)
{    
    if(preRefinedMesh)
        structure.refineMarkersOnRefinedMesh();
    
    dimensionSet dimensions = fieldRefineDemands.dimensions();
    Foam::dimensioned<Foam::scalar> val("fieldRefinement",dimensions,DONTCARE);
    doRefine = val;
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
        meshWasRefined |= refined;
        if(refined)
            structure.refineMarkersOnRefinedMesh();
        auto end = std::chrono::system_clock::now();
        Info<<"------------------------------------------------- took:"<<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<" milliseconds"<<Foam::nl;
    }
    mesh.write();
        
    return meshWasRefined;
}

bool Foam::MeshRefiner::refineMeshOnFluid()
{
    fieldRefinement();
    bool meshChange = false;
    List<bool> markerCell(fieldRefineDemands.size(),false);
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    for(const LagrangianMarker* oneMarker : markers)
    {
        label markerCellInd = oneMarker->getMarkerCell();
        markerCell[markerCellInd] = true;
    }
    
    const labelList& cellLevel = refinerPtr->meshCutter().cellLevel();
       
    std::unordered_map<label,label> levelToCount;
    for(label cLevel : cellLevel)
    {
        levelToCount[cLevel]++;
    }
    for(auto iter=levelToCount.begin(); iter!=levelToCount.end(); iter++)
        Info<<"MeshRefiner::level count:"<<iter->first<<" : "<<iter->second<<Foam::endl;
    
    
    std::unordered_map<label,label> levelToRefineCount;
    dimensionSet dimensions = fieldRefineDemands.dimensions();
    Foam::dimensioned<Foam::scalar> val("fieldRefinement",dimensions,DONTCARE);
    doRefine = val;
    for(label cellInd=0; cellInd<fieldRefineDemands.size(); cellInd++)
    {
        scalar refineDemand = fieldRefineDemands[cellInd];
        bool isMarkerCell = markerCell[cellInd];
        
        if(isMarkerCell)
        {
            if(refineDemand==REFINE)
            {
                if(!fluidRefineOnIB)
                    refineDemand = MUSTKEEP;
            }
            if(refineDemand==UNREFINE)
            {
                refineDemand = MUSTKEEP;
            }
        }
        
        if(refineDemand==REFINE || refineDemand==UNREFINE)
        {
            meshChange = true;
        }        
        doRefine[cellInd] = refineDemand;
        
        if(doRefine[cellInd]>REFINE_LIM)
        {
            levelToRefineCount[cellLevel[cellInd]]++;
        }
    }
    
    for(auto iter=levelToRefineCount.begin(); iter!=levelToRefineCount.end(); iter++)
        Info<<"MeshRefiner::level fluid refine count:"<<iter->first<<" : "<<iter->second<<Foam::endl;
    
    for(label cellInd=0; cellInd<doRefine.size(); cellInd++)
    {
        label level = cellLevel[cellInd];
        scalar refCommand = doRefine[cellInd];
        if(level==1)
        {
            if(refCommand>REFINE_LIM)
                FatalErrorInFunction<<"Not further refine"<<exit(FatalError);
        }
    }
    
    
    bool refined=false;
        
    if(meshChange)
    {
        refined = applyMeshAdaption();
        mesh.write();
    }
    
    return refined;
}

void Foam::MeshRefiner::fieldRefinement()
{
    dimensionSet dimensions = fieldRefineDemands.dimensions();
    Foam::dimensioned<Foam::scalar> val("fieldRefinement",dimensions,DONTCARE);
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
        structure.computeMeshSetup();
        structure.settleIntoRefinedMesh();
    }
    //const labelList& cellLevel = refinerPtr->meshCutter().cellLevel();
    //Info<<"cellLevel:"<<cellLevel<<Foam::nl;
    return refined;
}

Foam::scalar Foam::MeshRefiner::refinementDemandMerge
(
    Foam::scalar fieldDemand,
    Foam::scalar markerDemand
)
{
    FatalErrorInFunction<<"No longer in use"<<exit(FatalError);
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

