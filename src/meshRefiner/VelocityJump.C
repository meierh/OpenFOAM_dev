#include "VelocityJump.H"

#include <typeinfo>

Foam::VelocityJump::VelocityJump
(
    fvMesh& mesh,
    LineStructure& structure,
    volScalarField& doRefine,
    dictionary& dynamicMeshDict,
    volVectorField& velocity
):
VelocityRefiner(mesh,structure,doRefine,dynamicMeshDict,velocity),
maxCellJump("maxCellJump",doRefine)
{
    dictionary& topoChangerDict = dynamicMeshDict.subDict("topoChanger");
    dictionary& fluidCriterionDict = topoChangerDict.subDict("fluidCriterion");
    
    ITstream refineStream = fluidCriterionDict.lookup("refineJumpLimit");
    token refineToken;
    refineStream.read(refineToken);
    if(!refineToken.isScalar())
        FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChangerDict/fluidCriterion/refineJumpLimit -- must be scalar"<<exit(FatalError);
    refineJumpLimit = refineToken.scalarToken();
    
    ITstream unrefineStream = fluidCriterionDict.lookup("unrefineJumpLimit");
    token unrefineToken;
    unrefineStream.read(unrefineToken);
    if(!unrefineToken.isScalar())
        FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChangerDict/fluidCriterion/unrefineJumpLimit -- must be scalar"<<exit(FatalError);
    unrefineJumpLimit = unrefineToken.scalarToken();
    
    if(refineJumpLimit <= unrefineJumpLimit)
        FatalErrorInFunction<<"Error refineJumpLimit <= unrefineJumpLimit"<<exit(FatalError);
    
    if(fluidCriterionDict.found("timeOffset"))
    {
        ITstream timeOffsetStream = fluidCriterionDict.lookup("timeOffset");
        token timeOffsetToken;
        timeOffsetStream.read(timeOffsetToken);
        if(!timeOffsetToken.isScalar())
            FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChangerDict/fluidCriterion/timeOffset -- must be scalar"<<exit(FatalError);
        timeOffset = timeOffsetToken.scalarToken();
        timeOffsetSet = true;
    }
    else
        timeOffsetSet = false;
    
    if(fluidCriterionDict.found("maxRefinement"))
    {
        ITstream maxRefinementStream = fluidCriterionDict.lookup("maxRefinement");
        token maxRefinementToken;
        maxRefinementStream.read(maxRefinementToken);
        if(!maxRefinementToken.isLabel())
            FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChangerDict/fluidCriterion/maxRefinement -- must be scalar"<<exit(FatalError);
        maxRefinement = maxRefinementToken.labelToken();
        maxRefinementSet = true;
    }
    else
        maxRefinementSet = false;
}

void Foam::VelocityJump::fieldRefinement()
{
    dimensionSet dimensions = fieldRefineDemands.dimensions();
    Foam::dimensioned<Foam::scalar> val("fieldRefinement",dimensions,MUSTKEEP);
    fieldRefineDemands = val;
    
    maxCellJump = dimensioned<scalar>("maxCellJump",dimensions,-1);
    
    if(timeOffsetSet)
    {
        if(mesh.time().value() < timeOffset)
            return;
    }
    
    const labelList& cellLevel = refinerPtr->meshCutter().cellLevel();
    
    Info<<"VelocityJump::fieldRefinement: "<<"maxRefinementSet:"<<maxRefinementSet<<" maxRefinement:"<<maxRefinement<<Foam::nl;
    
    std::unordered_map<label,label> levelToCount;
    for(label cLevel : cellLevel)
    {
        levelToCount[cLevel]++;
    }
    
    for(auto iter=levelToCount.begin(); iter!=levelToCount.end(); iter++)
        Info<<"level count:"<<iter->first<<" : "<<iter->second<<Foam::endl;
    
    std::unordered_map<label,label> levelToRefineCount;   
    const cellList& cells = mesh.cells();
    if(cells.size()!=cellLevel.size())
        FatalErrorInFunction<<"Mismatch in cell number"<<exit(FatalError);
    const faceList& faces = mesh.faces();
    const labelList& owners = mesh.owner();
    const labelList& neighbours = mesh.neighbour();    
    for(label cellInd=0; cellInd<velocity.size(); cellInd++)
    {
        const cell& oneCell = cells[cellInd];
        scalar maxJump = 0;
        bool cellJumpSeen = false;
        for(label cellFaceInd=0; cellFaceInd<oneCell.size(); cellFaceInd++)
        {
            label faceInd = oneCell[cellFaceInd];
            if(faceInd>=neighbours.size())
                continue;
            
            label faceNeiCellInd = (owners[faceInd]==cellInd) ? neighbours[faceInd] : owners[faceInd];            
            if(faceNeiCellInd==-1)
                continue;
            
            label opposingFaceLabel = oneCell.opposingFaceLabel(faceInd,faces);
            if(opposingFaceLabel<0)
                continue;
            if(opposingFaceLabel>=neighbours.size())
                continue;
            
            label oppoFaceNeiCellInd = (owners[opposingFaceLabel]==cellInd) ? neighbours[opposingFaceLabel] : owners[opposingFaceLabel];
            
            if(oppoFaceNeiCellInd==-1)
                continue;

            vector velocityFaceSide =  velocity[faceNeiCellInd];
            vector velocityOppoFaceSide =  velocity[oppoFaceNeiCellInd];            
            vector diffVelocity = velocityFaceSide-velocityOppoFaceSide;
            

            for(label dim=0; dim<3; dim++)
                if(maxJump<std::abs(diffVelocity[dim]))
                    maxJump = std::abs(diffVelocity[dim]);
            cellJumpSeen = true;
        }
        
        maxCellJump[cellInd] = maxJump;
        
        fieldRefineDemands[cellInd] = DONTCARE;
        if(cellJumpSeen)
        {
            if(maxJump>refineJumpLimit)
            {
                if(maxRefinementSet && cellLevel[cellInd]<maxRefinement)
                {
                    levelToRefineCount[cellLevel[cellInd]]++;
                    fieldRefineDemands[cellInd] = REFINE;
                }
                else
                    fieldRefineDemands[cellInd] = MUSTKEEP;
            }
            else if(maxJump<unrefineJumpLimit)
                fieldRefineDemands[cellInd] = UNREFINE;
            else
                fieldRefineDemands[cellInd] = MUSTKEEP;
        }
    }
    
    for(auto iter=levelToRefineCount.begin(); iter!=levelToRefineCount.end(); iter++)
        Info<<"level fluid refine count:"<<iter->first<<" : "<<iter->second<<Foam::endl;
    
    label fluidRefineDemandCount=0;
    for(scalar refnDem : fieldRefineDemands)
        if(refnDem>REFINE_LIM)
            fluidRefineDemandCount++;
    Info<<"fluidRefineDemandCount:"<<fluidRefineDemandCount<<Foam::endl;    
    //maxCellJump.write();
}
