#include "SensitivityTemperatureInteraction.H"

Foam::SensitivityTemperatureInteraction::SensitivityTemperatureInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    const TemperatureInteraction& primalInteraction,
    volScalarField& adj_T,
    volScalarField& adj_fT,
    const IOdictionary& structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
SensitivityInteraction(mesh,structure,structureDict,modusFieldToMarker,modusMarkerToField),
primalInteraction(primalInteraction),
heatingDerivativeField(primalInteraction.getTemperatureField()),
adj_T(adj_T),
adj_fT(adj_fT)
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

void Foam::SensitivityTemperatureInteraction::solve(scalar timeStep)
{
    interpolateAdjTemperatureToMarkers();
    computeAdjCouplingHeatingOnMarkers(timeStep);
    interpolateAdjHeatingField();
}

void Foam::SensitivityTemperatureInteraction::recomputeMarkerValues()
{
    interpolateAdjTemperatureToMarkers();
    
    scalar virtualAdjTemperatureTimestep=0.1;
    IOobject fvSolutionIO("fvSolution","system",mesh.time(),IOobject::MUST_READ,IOobject::NO_WRITE);
    if(!fvSolutionIO.filePath("",true).empty())
    {
        IOdictionary fvSolutionDict(fvSolutionIO);
        dictionary& adj_pimpleDict = fvSolutionDict.subDict("adj_PIMPLE");

        ITstream virtualTemperatureTimestepStream = adj_pimpleDict.lookup("virtualTemperatureTimestep");
        token virtualTemperatureTimestepToken;
        virtualTemperatureTimestepStream.read(virtualTemperatureTimestepToken);
        if(!virtualTemperatureTimestepToken.isScalar())
            FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/virtualMomentumTimestep -- must be scalar"<<exit(FatalError);
        virtualAdjTemperatureTimestep = virtualTemperatureTimestepToken.scalarToken();
    }
    else
        FatalErrorInFunction<<"Missing file in system/fvSolution"<<exit(FatalError);
    computeAdjCouplingHeatingOnMarkers(virtualAdjTemperatureTimestep);
    
    interpolateAdjHeatingField();
}

void Foam::SensitivityTemperatureInteraction::interpolateAdjTemperatureToMarkers()
{
    // lambda_F_Tfluid = int 1/rho lambda_T_fluid delta(x-X) dOmega
    fieldToMarker<scalar>(adj_T,markerFluidAdjointTemperature);
    //markerFluidAdjointTemperature /= rho;
}

void Foam::SensitivityTemperatureInteraction::computeAdjCouplingHeatingOnMarkers(scalar timeStep)
{
    // lambda_T_fluidm = - lambda_F_Tfluid / delta t
    
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    
    if(markerFluidAdjointTemperature.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of markerFluidAdjointTemperature and markers"<<exit(FatalError);
    
    Info<<"Foam::SensitivityTemperatureInteraction::computeAdjCouplingHeatingOnMarkers: "<<timeStep<<Foam::nl;
    
    makerCouplingAdjointHeating.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        //LagrangianMarker* oneMarkerPtr = markers[markerInd];
        scalar markerAdjTemperature = 0;
        scalar adj_temperature = markerFluidAdjointTemperature[markerInd];
        makerCouplingAdjointHeating[markerInd] = (markerAdjTemperature-adj_temperature)/timeStep;
    }
}

void Foam::SensitivityTemperatureInteraction::interpolateAdjHeatingField()
{
    markerToField<scalar>(makerCouplingAdjointHeating,adj_fT);
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
    if(heatingDerivativeField.size()!=adj_fT.size())
        FatalErrorInFunction<<"Mismatch in field size!"<<exit(FatalError);    
    
    Field<scalar> temperatureForcingSensitivity(mesh.cells().size(),Foam::zero());
    for(label cellInd=0; cellInd<heatingDerivativeField.size(); cellInd++)
    {
        temperatureForcingSensitivity[cellInd] =  adj_T[cellInd] * -1 * heatingDerivativeField[cellInd];
    }
        
    return integrateField(temperatureForcingSensitivity);
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
    
    DynamicList<scalar> temperatureSensitivity(markers.size(),Foam::zero());
    for(std::size_t cellInd=0; cellInd<markers.size(); cellInd++)
    {
        temperatureSensitivity[cellInd] =  makerCouplingAdjointHeating[cellInd] * -1 * heatingDerivationMarkers[cellInd];
    }
    
    return integrateMarkers(temperatureSensitivity);
}
