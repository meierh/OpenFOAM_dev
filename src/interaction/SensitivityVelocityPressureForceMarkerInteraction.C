#include "SensitivityVelocityPressureForceMarkerInteraction.H"

Foam::SensitivityVelocityPressureForceInteraction::SensitivityVelocityPressureForceInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    const VelocityPressureForceInteraction& primalInteraction,
    volVectorField& adj_U,
    volVectorField& adj_fU,
    const IOdictionary& structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
SensitivityInteraction(mesh,structure,structureDict,modusFieldToMarker,modusMarkerToField),
primalInteraction(primalInteraction),
forcingDerivativeField(primalInteraction.getReferenceInOutField()),
adj_U(adj_U),
adj_fU(adj_fU)
{}

Foam::scalar Foam::SensitivityVelocityPressureForceInteraction::computeSensitivity
(
    const Parameter& para
)
{
    vector velocityForcingSensitivity = integrateVelocityForcingSensitivity(para);
    Info<<"velocityForcingSensitivity:"<<velocityForcingSensitivity<<Foam::endl;
    vector velocitySensitivity = integrateVelocitySensitivity(para);
    Info<<"velocitySensitivity:"<<velocitySensitivity<<Foam::endl;
    vector sensitivityVector = velocityForcingSensitivity+velocitySensitivity;
    Info<<"sensitivityVector:"<<sensitivityVector<<Foam::endl;
    scalar sensitivity = sensitivityVector[0]+sensitivityVector[1]+sensitivityVector[2];
    Info<<"sensitivity:"<<sensitivity<<Foam::endl;
    Pstream::gather<scalar>(sensitivity,std::plus<scalar>());
    Pstream::scatter<scalar>(sensitivity);
    return sensitivity;
}

void Foam::SensitivityVelocityPressureForceInteraction::solve(scalar timeStep)
{
    interpolateAdjVelocityToMarkers();
    computeAdjCouplingForceOnMarkers(timeStep);
    interpolateAdjFluidForceField();
}

void Foam::SensitivityVelocityPressureForceInteraction::recomputeMarkerValues()
{
    interpolateAdjVelocityToMarkers();
    
    scalar virtualAdjMomentumTimestep=0.1;
    IOobject fvSolutionIO("fvSolution","system",mesh.time(),IOobject::MUST_READ,IOobject::NO_WRITE);
    if(!fvSolutionIO.filePath("",true).empty())
    {
        IOdictionary fvSolutionDict(fvSolutionIO);
        dictionary& adj_pimpleDict = fvSolutionDict.subDict("adj_PIMPLE");

        ITstream virtualMomentumTimestepStream = adj_pimpleDict.lookup("virtualMomentumTimestep");
        token virtualMomentumTimestepToken;
        virtualMomentumTimestepStream.read(virtualMomentumTimestepToken);
        if(!virtualMomentumTimestepToken.isScalar())
            FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/virtualMomentumTimestep -- must be scalar"<<exit(FatalError);
        virtualAdjMomentumTimestep = virtualMomentumTimestepToken.scalarToken();
    }
    else
        FatalErrorInFunction<<"Missing file in system/fvSolution"<<exit(FatalError);
    computeAdjCouplingForceOnMarkers(virtualAdjMomentumTimestep);
    
    interpolateAdjFluidForceField();
}

void Foam::SensitivityVelocityPressureForceInteraction::interpolateAdjVelocityToMarkers()
{
    // lambda_Fu = int 1/rho lambda_u delta(x-X) dOmega
    fieldToMarker<vector>(adj_U,markerFluidAdjointVelocity);
    //markerFluidAdjointVelocity /= rho;
}

void Foam::SensitivityVelocityPressureForceInteraction::computeAdjCouplingForceOnMarkers(scalar timeStep)
{
    // lambda_Um = - lambda_Fu / delta t
    
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    
    if(markerFluidAdjointVelocity.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of markerFluidAdjointVelocity and markers"<<exit(FatalError);
    
    Info<<"Foam::SensitivityVelocityPressureForceInteraction::computeCouplingForceOnMarkers: "<<timeStep<<Foam::nl;
    
    makerCouplingAdjointForce.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        //LagrangianMarker* oneMarkerPtr = markers[markerInd];
        vector markerAdjVelocity = vector(0,0,0);
        vector adj_velocity = markerFluidAdjointVelocity[markerInd];
        makerCouplingAdjointForce[markerInd] = (markerAdjVelocity-adj_velocity)/timeStep;
    }
    
    vector sum = vector(0,0,0);
    for(vector val : makerCouplingAdjointForce)
        sum += val;
    Info<<"Foam::SensitivityVelocityPressureForceInteraction::computeCouplingForceOnMarkers::makerCouplingAdjointForce = "<<sum<<Foam::nl;
}

void Foam::SensitivityVelocityPressureForceInteraction::interpolateAdjFluidForceField()
{
    adj_fU = Foam::zero();
    markerToField<vector>(makerCouplingAdjointForce,adj_fU);
}

Foam::vector Foam::SensitivityVelocityPressureForceInteraction::integrateVelocityForcingSensitivity
(
    const Parameter& par
)
{   
    const DynamicList<vector>& F_U = primalInteraction.getMarkerCouplingForce();
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(F_U.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in marker size!"<<exit(FatalError);
    
    deriveParamMarkerToField<vector>(F_U,forcingDerivativeField,par);
    if(forcingDerivativeField.size()!=adj_fU.size())
        FatalErrorInFunction<<"Mismatch in field size!"<<exit(FatalError);    
    
    Field<vector> velocityForcingSensitivity(mesh.cells().size(),Foam::zero());
    for(label cellInd=0; cellInd<forcingDerivativeField.size(); cellInd++)
    {
        for(label dim=0; dim<3; dim++)
            velocityForcingSensitivity[cellInd][dim] =  adj_U[cellInd][dim] * -1 * forcingDerivativeField[cellInd][dim];
    }

    return integrateField<vector>(velocityForcingSensitivity);
}

Foam::vector Foam::SensitivityVelocityPressureForceInteraction::integrateVelocitySensitivity
(
    const Parameter& par
)
{   
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(makerCouplingAdjointForce.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in marker size!"<<exit(FatalError);
    
    DynamicList<vector> velocityDerivationMarkers(markers.size());
    deriveParamFieldToMarker<vector>(primalInteraction.getVelocityField(),velocityDerivationMarkers,par); 
    
    DynamicList<vector> velocitySensitivity(markers.size(),Foam::zero());
    for(std::size_t cellInd=0; cellInd<markers.size(); cellInd++)
    {
        for(label dim=0; dim<3; dim++)
            velocitySensitivity[cellInd][dim] =  makerCouplingAdjointForce[cellInd][dim] * -1 * velocityDerivationMarkers[cellInd][dim];
    }
        
    return integrateMarkers<vector>(velocitySensitivity);
}
