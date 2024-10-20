#include "SensitivityTemperatureInteraction.H"

Foam::SensitivityTemperatureInteraction::SensitivityTemperatureInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    const TemperatureInteraction& primalInteraction,
    volScalarField& input_adj_T,
    volScalarField& output_adj_Tf,
    const IOdictionary& structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
SensitivityInteraction(mesh,structure,structureDict,modusFieldToMarker,modusMarkerToField),
primalInteraction(primalInteraction),
heatingDerivativeField(primalInteraction.getTemperatureField()),
input_adj_T(input_adj_T),
output_adj_Tf(output_adj_Tf)
{}

Foam::scalar Foam::SensitivityTemperatureInteraction::computeSensitivity
(
    const Parameter& para
)
{
    scalar sensitivity = integrateTemperatureForcingSensitivity(para)+integrateTemperatureSensitivity(para);
    Pstream::gather(sensitivity,std::plus<scalar>());
    return sensitivity;
}

Foam::scalar Foam::SensitivityTemperatureInteraction::integrateTemperatureForcingSensitivity
(
    const Parameter& par
)
{
    const DynamicList<scalar>& F_T = primalInteraction.getMarkerCouplingHeating();
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(F_T.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in marker size!"<<exit(FatalError);
    
    deriveParamMarkerToField<scalar>(F_T,heatingDerivativeField,par);
    if(heatingDerivativeField.size()!=output_adj_Tf.size())
        FatalErrorInFunction<<"Mismatch in field size!"<<exit(FatalError);    
    
    scalar sumTemperatureForcingSensitivity = 0;
    for(label cellInd=0; cellInd<heatingDerivativeField.size(); cellInd++)
    {
        sumTemperatureForcingSensitivity +=  output_adj_Tf[cellInd] * -1 * heatingDerivativeField[cellInd];
    }
    return sumTemperatureForcingSensitivity;
}

Foam::scalar Foam::SensitivityTemperatureInteraction::integrateTemperatureSensitivity
(
    const Parameter& par
)
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(makerCouplingAdjointHeating.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in marker size!"<<exit(FatalError);
    
    DynamicList<scalar> heatingDerivationMarkers(markers.size());
    deriveParamFieldToMarker<scalar>(primalInteraction.getTemperatureField(),heatingDerivationMarkers,par); 
    
    scalar sumTemperatureSensitivity = 0;
    for(std::size_t cellInd=0; cellInd<markers.size(); cellInd++)
    {
        sumTemperatureSensitivity +=  makerCouplingAdjointHeating[cellInd] * -1 * heatingDerivationMarkers[cellInd];
    }
    return sumTemperatureSensitivity;
}
