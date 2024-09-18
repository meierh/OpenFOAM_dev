#include "FSIAction.H"

Foam::FSIAction::FSIAction
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
}

void Foam::FSIAction::computeRodForceMoment()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();

    if(makerCouplingForce.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of makerCouplingForce and markers"<<exit(FatalError);
    
    List<List<vector>> globalRodForce(Pstream::nProcs());
    List<List<vector>> globalRodMoment(Pstream::nProcs());
    List<List<label>> globalRodNumber(Pstream::nProcs());
    List<List<scalar>> globalParameter(Pstream::nProcs());  

    List<label>& rodNumber = globalRodNumber[Pstream::myProcNo()];
    List<scalar>& parameter = globalParameter[Pstream::myProcNo()];
    
    rodForce.resize(markers.size());
    rodMoment.resize(markers.size());
    rodNumber.resize(markers.size());
    parameter.resize(markers.size());
    
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
        
        rodNumber[markerInd] = oneMarker->getRodNumber();
        parameter[markerInd] = oneMarker->getMarkerParameter();
    }
    globalRodForce[Pstream::myProcNo()] = rodForce;
    globalRodMoment[Pstream::myProcNo()] = rodMoment;
    
    Pstream::gatherList(globalRodForce);
    Pstream::gatherList(globalRodMoment);
    Pstream::gatherList(globalRodNumber);
    Pstream::gatherList(globalParameter);
    
    if(Pstream::master())
    {
        std::pair<gsNurbs<scalar>,gsNurbs<scalar>> forcesMoment = assignForceOnRod(globalRodForce,globalRodMoment,globalRodNumber,globalParameter);
    }
    
    
    //Gather force and moment
}

Foam::vector Foam::FSIAction::getVelocity
(
    const LagrangianMarker* marker
)
{
    return marker->getMarkerVelocity();
}

std::pair<gsNurbs<Foam::scalar>,gsNurbs<Foam::scalar>> Foam::FSIAction::assignForceOnRod
(
    const List<List<vector>>& globalRodForce,
    const List<List<vector>>& globalRodMoment,
    const List<List<label>>& globalRodNumber,
    const List<List<scalar>>& globalParameter
)
{
    List<std::multimap<scalar,vector>> forces(structure.getNumberRods());
    List<std::multimap<scalar,vector>> moments(structure.getNumberRods());
    for(label proc=0; proc<Pstream::nProcs(); proc++)
    {
        const List<vector>& rodForce = globalRodForce[proc];
        const List<vector>& rodMoment = globalRodMoment[proc];
        const List<label>& rodNumber = globalRodNumber[proc];
        const List<scalar>& parameter = globalParameter[proc];
        
        for(label markerInd=0; markerInd<rodNumber.size(); markerInd++)
        {
            forces[rodNumber[markerInd]].insert({parameter[markerInd],rodForce[markerInd]});
            moments[rodNumber[markerInd]].insert({parameter[markerInd],rodMoment[markerInd]});
        }
    }
    
    List<std::map<scalar,vector>> forcesComb(structure.getNumberRods());
    List<std::map<scalar,vector>> momentsComb(structure.getNumberRods());
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
    
    List<List<scalar>> parametersList(structure.getNumberRods());
    List<List<vector>> forcesList(structure.getNumberRods());
    List<List<vector>> momentsList(structure.getNumberRods());
    for(label rodInd=0; rodInd<structure.getNumberRods(); rodInd++)
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
    
    //Transfer to nurbs
    FatalErrorInFunction<<"Incomplete implementation"<<exit(FatalError);
    
    return {};
}
