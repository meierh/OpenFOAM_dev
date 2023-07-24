#include "TfsiTapeStorage.H"

Foam::TfsiTapeStorage::TfsiTapeStorage
(   
    cutCellFvMesh& mesh 
):
mesh(mesh)
{
}

void Foam::TfsiTapeStorage::registerField
(
    volScalarField& field
)
{
    scalarTape.push_back(Tape<volScalarField>(field));
}

void Foam::TfsiTapeStorage::registerField
(
    volVectorField& field
)
{
    vectorTape.push_back(Tape<volVectorField>(field));
}

void Foam::TfsiTapeStorage::record
(
    scalar timeValue
)
{
    for(Tape<volScalarField>& oneField : scalarTape)
    {
        oneField.push(timeValue);
    }
    for(Tape<volVectorField>& oneField : vectorTape)
    {
        oneField.push(timeValue);
    }
    
    verticeTape.push(mesh.points())
    times.push(time);
}

bool Foam::TfsiTapeStorage::rewind
(
    scalar timeValue
)
{
    if(times.size()>0)
    {
        for(Tape<volScalarField>& oneField : scalarTape)
        {
            oneField.pop(timeValue);
        }
        for(Tape<volVectorField>& oneField : vectorTape)
        {
            oneField.pop(timeValue);
        }
    
        verticeTape.push(mesh.points())
        times.push(time);
        
        return true;
    }
    else
    {
        return false;
    }
}
