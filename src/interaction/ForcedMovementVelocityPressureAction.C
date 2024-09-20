#include "ForcedMovementVelocityPressureAction.H"
#include "codeStream.H"

Foam::ForcedMovementVelocityPressureAction::ForcedMovementVelocityPressureAction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volVectorField& input_U,
    volVectorField& output_Uf,
    const IOdictionary& structureDict,
    std::shared_ptr<MeshRefiner> refinement_
):
VelocityPressureForceInteraction(mesh,structure,input_U,output_Uf,refinement_),
structureDict(structureDict)
{
    Info<<"Created ForcedMovementVelocityPressureAction"<<Foam::endl;
}

void Foam::ForcedMovementVelocityPressureAction::preSolveMovement()
{
    moveRodsAndMarkers();
}

std::unique_ptr<Foam::List<Foam::List<Foam::vector>>> Foam::ForcedMovementVelocityPressureAction::getDeformation()
{   
    const dictionary& rodMovementFieldDict = structureDict.subDict("rodMovementField");
    List<keyType> rodMovementFieldKeys = rodMovementFieldDict.keys();
    
    auto movementListPtr = std::make_unique<List<List<vector>>>(rodMovementFieldKeys.size());
    List<List<vector>>& movementList = *movementListPtr;
    
    if(rodMovementFieldKeys.size()!=structure.getNumberRods())
    {
        Info<<"rodMovementFieldKeys:"<<rodMovementFieldKeys<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in movement field to rod number!"<<exit(FatalError);
    }
    
    label rodNumber = 0;
    for(keyType oneRodMoveFieldKey : rodMovementFieldKeys)
    {
        const dictionary& oneRodMovementDict = rodMovementFieldDict.subDict(oneRodMoveFieldKey);
        const entry* cdstr = oneRodMovementDict.lookupEntryPtr("move",false,false);
        ITstream stream = cdstr->stream();
        List<vector> oneRodMovementData(stream);        
        movementList[rodNumber] = oneRodMovementData;
    }
    return movementListPtr;
}

Foam::vector Foam::ForcedMovementVelocityPressureAction::getVelocity
(
    const LagrangianMarker* marker
)
{
    return marker->getMarkerVelocity();
}
