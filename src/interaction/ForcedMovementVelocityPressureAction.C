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
VelocityPressureForceInteraction(mesh,structure,input_U,output_Uf,refinement_,modusFieldToMarker,modusMarkerToField),
structureDict(structureDict)
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
    const dictionary& rodMovementFieldDict = structureDict.subDict("rodMovementField");
    List<keyType> rodMovementFieldKeys = rodMovementFieldDict.keys();
    
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
        const dictionary& oneRodMovementDict = rodMovementFieldDict.subDict(oneRodMoveFieldKey);
        const entry* cdstr = oneRodMovementDict.lookupEntryPtr("move",false,false);
        ITstream stream = cdstr->stream();
        List<vector> oneRodMovementData(stream);
        const gsNurbs<scalar>& deformation = structure.getDeformation(rodNumber);
        gsMatrix<scalar> defCoeffs;
        Structure::fitNurbsCoeffsToPoints(oneRodMovementData,deformation,defCoeffs);
        movementList[rodNumber].setSize(defCoeffs.rows());
        for(label coeffI=0; coeffI<movementList[rodNumber].size(); coeffI++)
        {
            movementList[rodNumber][coeffI][0] = defCoeffs(coeffI,0);
            movementList[rodNumber][coeffI][1] = defCoeffs(coeffI,1);
            movementList[rodNumber][coeffI][2] = defCoeffs(coeffI,2);
        }
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
