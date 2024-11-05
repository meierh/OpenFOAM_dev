#include "VelocityPressureForceMarkerInteraction.H"


Foam::VelocityPressureForceInteraction::VelocityPressureForceInteraction
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
FieldMarkerStructureInteraction(mesh,structure,structureDict,modusFieldToMarker,modusMarkerToField),
input_U(input_U),
output_Uf(output_Uf),
refinement_(refinement_)
{
    if(structureDict.found("recordRodForce"))
    {
        ITstream recordRodForceStream = structureDict.lookup("recordRodForce");
        token recordRodForceToken;
        recordRodForceStream.read(recordRodForceToken);
        if(!recordRodForceToken.isString())
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/recordRodForce -- must be string"<<exit(FatalError);
        recordRodForceFileName = recordRodForceToken.stringToken();
        if(Pstream::master())
            recordRodForceFile = std::make_unique<std::ofstream>(recordRodForceFileName);
        printSummedRodForces = true;
    }
    if(structureDict.found("recordRodMoment"))
    {
        ITstream recordRodMomentStream = structureDict.lookup("recordRodMoment");
        token recordRodMomentToken;
        recordRodMomentStream.read(recordRodMomentToken);
        if(!recordRodMomentToken.isString())
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/recordRodMoment -- must be string"<<exit(FatalError);
        recordRodMomentFileName = recordRodMomentToken.stringToken();
        if(Pstream::master())
            recordRodMomentFile = std::make_unique<std::ofstream>(recordRodMomentFileName);
        printSummedRodMoments = true;
    }
}

void Foam::VelocityPressureForceInteraction::solve()
{
    interpolateFluidVelocityToMarkers();
    computeCouplingForceOnMarkers();
    computeRodForceMoment();
    interpolateFluidForceField();
    
    if(printSummedRodForces)
    {
        vector sumForcesVal = sumForces();
        if(Pstream::master())
        {
            (*recordRodForceFile)<<mesh.time().value()<<":  "<<sumForcesVal[0]<<" "<<sumForcesVal[1]<<" "<<sumForcesVal[2]<<std::endl;
        }
    }
    
    if(printSummedRodMoments)
    {
        vector sumMomentsVal = sumMoments();
        if(Pstream::master())
        {
            (*recordRodMomentFile)<<mesh.time().value()<<":  "<<sumMomentsVal[0]<<" "<<sumMomentsVal[1]<<" "<<sumMomentsVal[2]<<std::endl;
        }
    }
}

void Foam::VelocityPressureForceInteraction::store()
{
    std::tuple<DynamicList<vector>,DynamicList<vector>,DynamicList<vector>,DynamicList<vector>>& markerValues = storage[mesh.time().value()];
    
    std::get<0>(markerValues) = markerFluidVelocity;
    std::get<1>(markerValues) = markerCouplingForce;
    std::get<2>(markerValues) = rodForce;
    std::get<3>(markerValues) = rodMoment;
}

void Foam::VelocityPressureForceInteraction::setToTime(scalar time)
{
    if(storage.find(time)==storage.end())
        FatalErrorInFunction<<"No data exists for "<<time<<exit(FatalError);
    
    std::tuple<DynamicList<vector>,DynamicList<vector>,DynamicList<vector>,DynamicList<vector>>& markerValues = storage[time];
    
    markerFluidVelocity = std::get<0>(markerValues);
    markerCouplingForce = std::get<1>(markerValues);
    rodForce = std::get<2>(markerValues);
    rodMoment = std::get<3>(markerValues);
}

void Foam::VelocityPressureForceInteraction::interpolateFluidVelocityToMarkers()
{
    fieldToMarker<vector>(input_U,markerFluidVelocity);
}

void Foam::VelocityPressureForceInteraction::computeCouplingForceOnMarkers()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    
    if(markerFluidVelocity.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of markerFluidVelocity and markers"<<exit(FatalError);
    
    scalar deltaT = mesh.time().deltaTValue();
    Info<<"Foam::VelocityPressureForceInteraction::computeCouplingForceOnMarkers: "<<deltaT<<Foam::nl;
    
    markerCouplingForce.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarkerPtr = markers[markerInd];
        vector markerVelocity = getVelocity(oneMarkerPtr);
        vector velocity = markerFluidVelocity[markerInd];
        markerCouplingForce[markerInd] = (markerVelocity-velocity)/deltaT;
    }
}

void Foam::VelocityPressureForceInteraction::interpolateFluidForceField()
{
    markerToField<vector>(markerCouplingForce,output_Uf);
}

void Foam::VelocityPressureForceInteraction::moveMarkers()
{  
    List<bool> prevRodInMesh = structure.getRodInMesh();
    structure.moveMarkersOnRodMovement();
    const List<bool>& rodInMesh = structure.getRodInMesh();
    for(label rodNumber=0; rodNumber<prevRodInMesh.size(); rodNumber++)
    {
        if(!prevRodInMesh[rodNumber])
        {
            if(rodInMesh[rodNumber])
            {
                structure.createMarkersOnRod(rodNumber);
            }
        }
    }
}

Foam::vector Foam::VelocityPressureForceInteraction::sumForces
(
    std::function<bool(LagrangianMarker)> condition
)
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    
    if(rodForce.size()!=static_cast<label>(markers.size()))
    {
        Info<<"rodForce.size():"<<rodForce.size()<<Foam::endl;
        Info<<"markers.size():"<<markers.size()<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in size of rodForce and markers"<<exit(FatalError);
    }
    
    vector result = Foam::zero();
    for(std::size_t i=0; i<markers.size(); i++)
    {
        LagrangianMarker* oneMarker = markers[i];
        if(condition(*oneMarker))
            result += rodForce[i];
    }
    Pstream::gather(result,std::plus<vector>());
    Pstream::scatter(result);
    
    return result;
}

Foam::vector Foam::VelocityPressureForceInteraction::sumMoments
(
    std::function<bool(LagrangianMarker)> condition
)
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    
    if(rodMoment.size()!=static_cast<label>(markers.size()))
    {
        Info<<"rodMoment.size():"<<rodMoment.size()<<Foam::endl;
        Info<<"markers.size():"<<markers.size()<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in size of rodMoment and markers"<<exit(FatalError);
    }
    
    vector result = Foam::zero();
    for(std::size_t i=0; i<markers.size(); i++)
    {
        LagrangianMarker* oneMarker = markers[i];
        if(condition(*oneMarker))
            result += rodMoment[i];
    }
    Pstream::gather(result,std::plus<vector>());
    Pstream::scatter(result);
    
    return result;
}

void Foam::VelocityPressureForceInteraction::meshMarkerAdaptation()
{
    if(refinement_)
    {
        Info<<"|||||||||||||||||||||||||Do refinement|||||||||||||||||||||||||"<<Foam::endl;
        //bool meshRefined = refinement_->refineMeshOnStaticMarkers();
        //meshRefined = refinement_->refineMeshAndMarkers(meshRefined);
        refinement_->refineMeshAndMarkers();
        Info<<"||||||||||||||||||||||||Done refinement||||||||||||||||||||||||"<<Foam::endl;
    }
    structure.finalizeMarkers(true);
}

void Foam::VelocityPressureForceInteraction::computeRodForceMoment()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();

    if(markerCouplingForce.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of makerCouplingForce and markers"<<exit(FatalError);
    
    rodForce.resize(markers.size());
    rodMoment.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarker = markers[markerInd];
        scalar volume = oneMarker->getMarkerVolume();
        scalar rho = 1.225;
        
        rodForce[markerInd] = rho*volume*markerCouplingForce[markerInd];
        
        vector basePnt;
        Structure::rodEval(oneMarker->getBaseRod(),oneMarker->getMarkerParameter(),basePnt);
        vector vectorToMarker = oneMarker->getMarkerPosition()-basePnt;
        vector momentum = vectorToMarker^markerCouplingForce[markerInd];
        rodMoment[markerInd] = rho*momentum*volume;
    }
}
