#include "ForcedMovementVelocityPressureAction.H"

Foam::ForcedMovementVelocityPressureAction::ForcedMovementVelocityPressureAction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volVectorField& input_U,
    volVectorField& output_Uf,
    const IOdictionary& structureDict
):
VelocityPressureForceInteraction(mesh,structure,input_U,output_Uf),
structureDict(structureDict)
{
    Info<<"Created ForcedMovementVelocityPressureAction"<<Foam::endl;
    readMovementFromDict();
}

void Foam::ForcedMovementVelocityPressureAction::preSolveMovement()
{
    std::unique_ptr<List<List<vector>>> deformationCoeffsPtr = readMovementFromDict();
    List<List<vector>> deformationCoeffs = *deformationCoeffsPtr;
    
    
}

Foam::vector Foam::ForcedMovementVelocityPressureAction::getVelocity
(
    const LagrangianMarker* marker
)
{
    return vector(0,0,0);
}

std::unique_ptr<Foam::List<Foam::List<Foam::vector>>> Foam::ForcedMovementVelocityPressureAction::readMovementFromDict()
{
    const dictionary& rodMovementFieldDict = structureDict.subDict("rodMovementField");
    List<keyType> rodMovementFieldKeys = rodMovementFieldDict.keys();
    
    auto movementListPtr = std::make_unique<List<List<vector>>>(rodMovementFieldKeys.size());
    List<List<vector>>& movementList = *movementListPtr;
    
    if(rodMovementFieldKeys.size()!=structure.getNumberRods())
        FatalErrorInFunction<<"Mismatch in movement field to rod number!"<<exit(FatalError);
    
    label rodNumber = 0;
    for(keyType oneRodMoveFieldKey : rodMovementFieldKeys)
    {
        const dictionary& oneRodMovementDict = rodMovementFieldDict.subDict(oneRodMoveFieldKey);
        List<vector> oneRodMovementData = oneRodMovementDict.lookup("move");
        
        if(oneRodMovementData.size()!=structure.numberCoeffs(rodNumber))
            FatalErrorInFunction<<"Wrong number of coefficient data at "<<oneRodMovementDict.name()<<Foam::endl
            <<"Expected "<<structure.numberCoeffs(rodNumber)<<" but got "<<oneRodMovementData.size()<<exit(FatalError);
        
        movementList[rodNumber] = oneRodMovementData;
    }
    return movementListPtr;
}
