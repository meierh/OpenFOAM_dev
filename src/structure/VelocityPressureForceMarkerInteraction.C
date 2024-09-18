#include "VelocityPressureForceMarkerInteraction.H"


Foam::VelocityPressureForceInteraction::VelocityPressureForceInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volVectorField& input_U,
    volVectorField& output_Uf,
    std::shared_ptr<MeshRefiner> refinement_
):
FieldMarkerStructureInteraction(mesh,structure),
input_U(input_U),
output_Uf(output_Uf),
refinement_(refinement_)
{
}

void Foam::VelocityPressureForceInteraction::solve()
{
    interpolateFluidVelocityToMarkers();
    computeCouplingForceOnMarkers();
    computeRodForceMoment();
    interpolateFluidForceField();
}

void Foam::VelocityPressureForceInteraction::store()
{
    std::tuple<DynamicList<vector>,DynamicList<vector>,DynamicList<vector>,DynamicList<vector>>& markerValues = storage[mesh.time().value()];
    
    std::get<0>(markerValues) = markerFluidVelocity;
    std::get<1>(markerValues) = makerCouplingForce;
    std::get<2>(markerValues) = rodForce;
    std::get<3>(markerValues) = rodMoment;
}

void Foam::VelocityPressureForceInteraction::setToTime(scalar time)
{
    if(storage.find(time)==storage.end())
        FatalErrorInFunction<<"No data exists for "<<time<<exit(FatalError);
    
    std::tuple<DynamicList<vector>,DynamicList<vector>,DynamicList<vector>,DynamicList<vector>>& markerValues = storage[time];
    
    markerFluidVelocity = std::get<0>(markerValues);
    markerFluidVelocity = std::get<1>(markerValues);
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
    
    makerCouplingForce.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarkerPtr = markers[markerInd];
        vector markerVelocity = getVelocity(oneMarkerPtr);
        vector velocity = markerFluidVelocity[markerInd];
        makerCouplingForce[markerInd] = (markerVelocity-velocity)/deltaT;
    }
}

void Foam::VelocityPressureForceInteraction::interpolateFluidForceField()
{
    markerToField<vector>(makerCouplingForce,output_Uf);
}

void Foam::VelocityPressureForceInteraction::moveRodsAndMarkers()
{
    List<bool> prevRodInMesh = structure.getRodInMesh();
    std::unique_ptr<List<List<vector>>> defPtr = getDeformation();
    if(defPtr)
    {
        structure.setDeformation(*defPtr);
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
        if(refinement_)
        {
            refinement_->refineMeshOnStaticMarkers();
            refinement_->refineMeshAndMarkers();
        }
        structure.finalizeMarkers();    
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

void Foam::VelocityPressureForceInteraction::computeRodForceMoment()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();

    if(makerCouplingForce.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of makerCouplingForce and markers"<<exit(FatalError);
    
    rodForce.resize(markers.size());
    rodMoment.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarker = markers[markerInd];
        scalar volume = oneMarker->getMarkerVolume();
        scalar rho = 1.225;
        
        rodForce[markerInd] = rho*volume*makerCouplingForce[markerInd];
        
        vector basePnt;
        Structure::rodEval(oneMarker->getBaseRod(),oneMarker->getMarkerParameter(),basePnt);
        vector vectorToMarker = oneMarker->getMarkerPosition()-basePnt;
        vector momentum = vectorToMarker^makerCouplingForce[markerInd];
        rodMoment[markerInd] = rho*momentum*volume;
    }
}

void Foam::VelocityPressureForceInteraction::assignForceOnRod()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    
    computeRodForceMoment();
    List<std::multimap<scalar,vector>> forces(rodMesh->m_Rods.size());
    List<std::multimap<scalar,vector>> moments(rodMesh->m_Rods.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarker = markers[markerInd];
        label rodNumber = oneMarker->getRodNumber();
        scalar markerParameter = oneMarker->getMarkerParameter();
        forces[rodNumber].insert({markerParameter,rodForce[markerInd]});
        moments[rodNumber].insert({markerParameter,rodMoment[markerInd]});        
    }
    List<std::map<scalar,vector>> forcesComb(rodMesh->m_Rods.size());
    List<std::map<scalar,vector>> momentsComb(rodMesh->m_Rods.size());
    for(label rodInd=0; rodInd<forcesComb.size(); rodInd++)
    {
        for(auto iterForce=forces[rodInd].begin(); iterForce!=forces[rodInd].end(); iterForce++)
        {
            forcesComb[rodInd][iterForce->first] += iterForce->second;
        }
        for(auto iterMom=moments[rodInd].begin(); iterMom!=moments[rodInd].end(); iterMom++)
        {
            momentsComb[rodInd][iterMom->first] += iterMom->second;
        }
    }
    
    List<List<scalar>> parametersList(forcesComb.size());
    List<List<vector>> forcesList(forcesComb.size());
    List<List<vector>> momentsList(forcesComb.size());
    for(label rodInd=0; rodInd<forcesComb.size(); rodInd++)
    {
        parametersList[rodInd].resize(forcesComb[rodInd].size());
        forcesList[rodInd].resize(forcesComb[rodInd].size());
        momentsList[rodInd].resize(forcesComb[rodInd].size());
        label index = 0;
        for(auto iterForce=forcesComb[rodInd].begin(); iterForce!=forcesComb[rodInd].end(); iterForce++)
        {
            parametersList[rodInd][index] = iterForce->first;
            forcesList[rodInd][index] = iterForce->second;
            index++;
        }
        index = 0;
        for(auto iterMom=momentsComb[rodInd].begin(); iterMom!=momentsComb[rodInd].end(); iterMom++)
        {
            momentsList[rodInd][index] = iterMom->second;
            index++;
        }
    }
}

Foam::vector Foam::VelocityPressureForceInteraction::getVelocity
(
    const LagrangianMarker* marker
)
{
    return vector(0,0,0);
}
