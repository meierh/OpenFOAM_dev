#include "ForcedMovementVelocityPressureAction.H"
#include "codeStream.H"

Foam::ForcedMovementVelocityPressureAction::ForcedMovementVelocityPressureAction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volVectorField& input_U,
    volVectorField& output_Uf,
    const IOdictionary& structureDict,
    std::shared_ptr<MeshRefiner> refinement_,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
VelocityPressureForceInteraction(mesh,structure,input_U,output_Uf,structureDict,refinement_,modusFieldToMarker,modusMarkerToField)
{
    Info<<"Created ForcedMovementVelocityPressureAction"<<Foam::endl;
}

void Foam::ForcedMovementVelocityPressureAction::preSolveMovement()
{
    std::unique_ptr<List<List<vector>>> allRodsDeformation = readDeformationDict();
    structure.setDeformation(*allRodsDeformation);
    moveMarkers();
}

void Foam::ForcedMovementVelocityPressureAction::preSolveMarkerMeshAdaption()
{
    meshMarkerAdaptation();
};


std::unique_ptr<Foam::List<Foam::List<Foam::vector>>> Foam::ForcedMovementVelocityPressureAction::readDeformationDict()
{   
    Info<<"Foam::ForcedMovementVelocityPressureAction::readDeformationDict()"<<Foam::endl;
    const dictionary& rodMovementFieldDict = structureDict.subDict("rodMovementField");
    List<keyType> rodMovementFieldKeys = rodMovementFieldDict.keys();
    Info<<"rodMovementFieldKeys:"<<rodMovementFieldKeys<<Foam::endl;
    
    auto movementListPtr = std::make_unique<List<List<vector>>>(rodMovementFieldKeys.size());
    List<List<vector>>& movementList = *movementListPtr;
    
    if(rodMovementFieldKeys.size()!=structure.getNumberRods())
    {
        Info<<"rodMovementFieldKeys:"<<rodMovementFieldKeys<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in movement field to rod number!"<<exit(FatalError);
    }
    
    for(label rodNumber=0; rodNumber<structure.getNumberRods(); rodNumber++)
    {
        keyType oneRodMoveFieldKey = rodMovementFieldKeys[rodNumber];
        Info<<"rodNumber:"<<rodNumber<<" : "<<oneRodMoveFieldKey<<Foam::endl;
        const dictionary& oneRodMovementDict = rodMovementFieldDict.subDict(oneRodMoveFieldKey);
        const entry* cdstr = oneRodMovementDict.lookupEntryPtr("move",false,false);
        ITstream stream = cdstr->stream();
        List<vector> oneRodMovementData(stream);
        Info<<"oneRodMovementData:"<<oneRodMovementData<<Foam::endl;
        const gsNurbs<scalar>& deformation = structure.getDeformation(rodNumber);
        label nbrDefParameters = deformation.coefs().rows();
        if(nbrDefParameters!=oneRodMovementData.size())
            FatalErrorInFunction<<"Mismatch in given deformation parameter size! Is "<<oneRodMovementData.size()<<" but should be "<<nbrDefParameters<<exit(FatalError);
        if(deformation.coefs().cols()!=3)
            FatalErrorInFunction<<"Deformation has non vector type parameter"<<exit(FatalError);
        movementList[rodNumber].setSize(nbrDefParameters);
        for(label coeffI=0; coeffI<movementList[rodNumber].size(); coeffI++)
        {
            movementList[rodNumber][coeffI] = oneRodMovementData[coeffI];
        }
    }
    Info<<"movementList:"<<movementList<<Foam::endl;
    return movementListPtr;
}

Foam::vector Foam::ForcedMovementVelocityPressureAction::getVelocity
(
    const LagrangianMarker* marker
)
{
    return marker->getMarkerVelocity();
}
