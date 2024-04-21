#include "MarkerImplementation.H"

Foam::VelocityPressureForceInteraction::VelocityPressureForceInteraction
(
    dynamicRefineFvMesh& mesh,
    volVectorField& input_U,
    volScalarField& input_p,
    volVectorField& output_Uf
):
FieldMarkerStructureInteraction(mesh),
input_U(input_U),
input_p(input_p),
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
    for(label markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarkerPtr = markers[markerInd];
        vector markerVelocity = oneMarkerPtr->getMarkerVelocity();
        vector fluidVelocity = markerFluidVelocity[markerInd];
        makerCouplingForce[markerInd] = (markerVelocity-fluidVelocity)/deltaT;
        //Info<<markerInd<<" velocity:"<<fluidVelocity<<" couplingForce:"<<makerCouplingForce[markerInd]<<Foam::endl;
    }
}

void Foam::VelocityPressureForceInteraction::interpolateFluidForceField()
{
    scalar deltaT = mesh.time().deltaTValue();
    markerToField<vector>(makerCouplingForce,output_Uf);
    for(vector uf : output_Uf)
        uf = uf*deltaT;
}
