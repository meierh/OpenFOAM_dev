#include "SensitivityVelocityPressure.H"

Foam::SensitivityVelocityPressure::SensitivityVelocityPressure
(
    const fvMesh& mesh,
    LineStructure& structure,
    const VelocityPressureForceInteraction& primalInteraction,
    const IOdictionary& structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
SensitivityInteraction(mesh,structure,structureDict,modusFieldToMarker,modusMarkerToField),
primalInteraction(primalInteraction),
forcingDerivativeField(primalInteraction.getReferenceInOutField())
{}

Foam::scalar Foam::SensitivityVelocityPressure::computeSensitivity
(
    const Parameter& para
)
{
    Info<<"computeSensitivity:"<<para.to_string()<<Foam::nl;
    vector sensitivityVector = integrateVelocityForcingSensitivity(para)+integrateVelocitySensitivity(para);
    scalar sensitivity = sensitivityVector[0]+sensitivityVector[1]+sensitivityVector[2];
    Pstream::gather<scalar>(sensitivity,std::plus<scalar>());
    Pstream::scatter<scalar>(sensitivity);
    return sensitivity;
}

Foam::vector Foam::SensitivityVelocityPressure::integrateVelocityForcingSensitivity
(
    const Parameter& par
)
{
    const DynamicList<vector>& F_U = primalInteraction.getMarkerCouplingForce();
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(F_U.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in marker size!"<<exit(FatalError);
    
    deriveParamMarkerToField<vector>(F_U,forcingDerivativeField,par);
    if(forcingDerivativeField.size()!=output_adj_Uf.size())
        FatalErrorInFunction<<"Mismatch in field size!"<<exit(FatalError);    
    
    vector sumVelocityForcingSensitivity = vector(0,0,0);
    for(label cellInd=0; cellInd<forcingDerivativeField.size(); cellInd++)
    {
        for(label dim=0; dim<3; dim++)
            sumVelocityForcingSensitivity[dim] +=  output_adj_Uf[cellInd][dim] * -1 * forcingDerivativeField[cellInd][dim];
    }
    return sumVelocityForcingSensitivity;
}

Foam::vector Foam::SensitivityVelocityPressure::integrateVelocitySensitivity
(
    const Parameter& par
)
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(makerCouplingAdjointForce.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in marker size!"<<exit(FatalError);
    
    DynamicList<vector> velocityDerivationMarkers(markers.size());
    //deriveParamFieldToMarker<vector>(primalInteraction.getVelocityField(),velocityDerivationMarkers,par); 
    
    vector sumTemperatureSensitivity = vector(0,0,0);
    for(std::size_t cellInd=0; cellInd<markers.size(); cellInd++)
    {
        for(label dim=0; dim<3; dim++)
            sumTemperatureSensitivity[dim] +=  makerCouplingAdjointForce[cellInd][dim] * -1 * velocityDerivationMarkers[cellInd][dim];
    }
    return sumTemperatureSensitivity;
}
