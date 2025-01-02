#include "FixedTemperatureAction.H"

Foam::FixedTemperatureAction::FixedTemperatureAction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volScalarField& input_T,
    volScalarField& output_Tf,
    const IOdictionary& structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
TemperatureInteraction(mesh,structure,input_T,output_Tf,structureDict,modusFieldToMarker,modusMarkerToField),
fixedTemperature(readTemperatureFromDict())
{
}

Foam::scalar Foam::FixedTemperatureAction::getTemperature
(
    const LagrangianMarker* marker
)
{
    return fixedTemperature;
}

Foam::scalar Foam::FixedTemperatureAction::readTemperatureFromDict()
{   
    ITstream rodTemperatureStream = structureDict.lookup("rodTemperature");
    token rodTemperatureToken;
    rodTemperatureStream.read(rodTemperatureToken);
    if(!rodTemperatureToken.isScalar())
        FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodTemperature -- must be scalar"<<exit(FatalError);
    scalar temperature = rodTemperatureToken.scalarToken();
    return temperature;
}
