#include "TemperatureInteraction.H"

Foam::TemperatureInteraction::TemperatureInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volScalarField& input_T,
    volScalarField& output_Tf,
    const IOdictionary& structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
FieldMarkerStructureInteraction(mesh,structure,modusFieldToMarker,modusMarkerToField),
structureDict(structureDict),
input_T(input_T),
output_Tf(output_Tf)
{
}

void Foam::TemperatureInteraction::solve()
{
    interpolateTemperatureToMarkers();
    computeCouplingHeatingOnMarkers();
    computeRodHeating();
    interpolateHeatingField();
}

void Foam::TemperatureInteraction::store()
{
    std::tuple<DynamicList<scalar>,DynamicList<scalar>,DynamicList<scalar>>& markerValues = storage[mesh.time().value()];
    
    std::get<0>(markerValues) = markerFluidTemperature;
    std::get<1>(markerValues) = makerCouplingHeating;
    std::get<2>(markerValues) = rodHeating;
}

void Foam::TemperatureInteraction::setToTime(scalar time)
{
    if(storage.find(time)==storage.end())
        FatalErrorInFunction<<"No data exists for "<<time<<exit(FatalError);
    
    std::tuple<DynamicList<scalar>,DynamicList<scalar>,DynamicList<scalar>>& markerValues = storage[time];
    
    markerFluidTemperature = std::get<0>(markerValues);
    makerCouplingHeating = std::get<1>(markerValues);
    rodHeating = std::get<2>(markerValues);
}

void Foam::TemperatureInteraction::interpolateTemperatureToMarkers()
{
    fieldToMarker<scalar>(input_T,markerFluidTemperature);
}

void Foam::TemperatureInteraction::computeCouplingHeatingOnMarkers()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(markerFluidTemperature.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of markerTemperature and markers"<<exit(FatalError);
    
    scalar deltaT = mesh.time().deltaTValue();
    
    makerCouplingHeating.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarkerPtr = markers[markerInd];
        scalar markerTemperature = getTemperature(oneMarkerPtr);
        scalar temperature = markerFluidTemperature[markerInd];
        makerCouplingHeating[markerInd] = (markerTemperature-temperature)/deltaT;
    }
}

void Foam::TemperatureInteraction::interpolateHeatingField()
{
    scalar deltaT = mesh.time().deltaTValue();
    markerToField<scalar>(makerCouplingHeating,output_Tf);
    for(scalar Tf : output_Tf)
        Tf = Tf*deltaT;
}

void Foam::TemperatureInteraction::computeRodHeating()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(makerCouplingHeating.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of makerCouplingHeating and markers"<<exit(FatalError);
    
    rodHeating.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarker = markers[markerInd];
        scalar volume = oneMarker->getMarkerVolume();
        scalar cp = 1.0035;
        scalar rho = 1.225;
        
        rodHeating[markerInd] = rho*cp*volume*makerCouplingHeating[markerInd];
    }
}

Foam::scalar Foam::TemperatureInteraction::sumHeating
(
    std::function<bool(LagrangianMarker)> condition
)
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(rodHeating.size()!=static_cast<label>(markers.size()))
    {
        Info<<"rodHeating.size():"<<rodHeating.size()<<Foam::endl;
        Info<<"markers.size():"<<markers.size()<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in size of rodHeating and markers"<<exit(FatalError);
    }
    
    scalar result = Foam::zero();
    for(std::size_t i=0; i<markers.size(); i++)
    {
        LagrangianMarker* oneMarker = markers[i];
        if(condition(*oneMarker))
            result += rodHeating[i];
    }
    return result;
}

Foam::scalar Foam::TemperatureInteraction::getTemperature
(
    const LagrangianMarker* marker
)
{
    return 0;
}
