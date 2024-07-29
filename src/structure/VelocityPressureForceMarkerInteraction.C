#include "VelocityPressureForceMarkerInteraction.H"


Foam::VelocityPressureForceInteraction::VelocityPressureForceInteraction
(
    dynamicRefineFvMesh& mesh,
    LineStructure& structure,
    volVectorField& input_U,
    volVectorField& output_Uf
):
FieldMarkerStructureInteraction(mesh,structure),
input_U(input_U),
output_Uf(output_Uf)
{
}

void Foam::VelocityPressureForceInteraction::interpolateFluidVelocityToMarkers()
{
    fieldToMarker<vector>(input_U,markerFluidVelocity);
}

void Foam::VelocityPressureForceInteraction::computeCouplingForceOnMarkers()
{
    if(markerFluidVelocity.size()!=markers.size())
        FatalErrorInFunction<<"Mismatch in size of markerFluidVelocity and markers"<<exit(FatalError);
    
    scalar deltaT = mesh.time().deltaTValue();
    
    makerCouplingForce.resize(markers.size());
    for(label rodInd=0; rodInd<markers.size(); rodInd++)
    {
        for(label markerInd=0; markerInd<markers.size(); markerInd++)
        {
            LagrangianMarker* oneMarkerPtr = markers[rodInd][markerInd];
            vector markerVelocity = oneMarkerPtr->getMarkerVelocity();
            vector fluidVelocity = markerFluidVelocity[rodInd][markerInd];
            makerCouplingForce[markerInd] = (markerVelocity-fluidVelocity)/deltaT;
        }
    }
}

void Foam::VelocityPressureForceInteraction::interpolateFluidForceField()
{
    scalar deltaT = mesh.time().deltaTValue();
    markerToField<vector>(makerCouplingForce,output_Uf);
}

vector Foam::VelocityPressureForceInteraction::sumForces
(
    std::function<bool(LagrangianMarker)> condition
)
{
    if(rodForce.size()!=markers.size())
    {
        Info<<"rodForce.size():"<<rodForce.size()<<Foam::endl;
        Info<<"markers.size():"<<markers.size()<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in size of rodForce and markers"<<exit(FatalError);
    }
    
    vector result = Foam::zero();
    for(label i=0; i<markers.size(); i++)
    {
        LagrangianMarker* oneMarker = markers[i];
        if(condition(*oneMarker))
            result += rodForce[i];
    }
    return result;
}

void Foam::VelocityPressureForceInteraction::computeRodForceMoment()
{
    if(makerCouplingForce.size()!=markers.size())
        FatalErrorInFunction<<"Mismatch in size of makerCouplingForce and markers"<<exit(FatalError);
    
    rodForce.resize(markers.size());
    rodMoment.resize(markers.size());
    for(label markerInd=0; markerInd<markers.size(); markerInd++)
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
    computeRodForceMoment();
    List<std::multimap<scalar,vector>> forces(rodMesh->m_Rods.size());
    List<std::multimap<scalar,vector>> moments(rodMesh->m_Rods.size());
    for(label markerInd=0; markerInd<markers.size(); markerInd++)
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
