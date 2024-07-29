#include "TestInteraction.H"

Foam::TestInteraction::TestInteraction
(
    dynamicRefineFvMesh& mesh,
    volScalarField& testField
):
FieldMarkerStructureInteraction(mesh),
testField(testField)
{}

void Foam::TestInteraction::printSupportToField()
{
    testField = Foam::zero();
    for(label markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarkerPtr = markers[markerInd];
        for(auto cellIter=oneMarkerPtr->getSupportCells().begin();
            cellIter!=oneMarkerPtr->getSupportCells().end();
            cellIter++)
        {
            testField[std::get<2>(*cellIter)] += 1;
        }
    }
}
