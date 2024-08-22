#include "FixedTemperatureAction.H"

Foam::FixedTemperatureAction::FixedTemperatureAction
(
    fvMesh& mesh,
    LineStructure& structure,
    volScalarField& input_T,
    volScalarField& output_Tf,
    scalar temperature
):
TemperatureInteraction(mesh,structure,input_T,output_Tf),
fixedTemperature(temperature)
{
}

Foam::scalar Foam::FixedTemperatureAction::getTemperature
(
    const LagrangianMarker* marker
)
{
    return fixedTemperature;
}
