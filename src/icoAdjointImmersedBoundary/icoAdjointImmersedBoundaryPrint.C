#include "icoAdjointImmersedBoundary.H"

void Foam::printMinMax(surfaceScalarField const& field)
{
    if(field.size()<1)
        FatalErrorInFunction<<"Error"<<exit(FatalError);
    scalar minVal = field[0];
    scalar maxVal = field[0];
    for(label i=0; i<field.size(); i++)
    {
        minVal = std::min(minVal,field[i]);
        maxVal = std::max(maxVal,field[i]);
    }
    Info<<"  Min: "<<minVal<<" Max: "<<maxVal<<Foam::nl;
}

void Foam::printMinMax(volScalarField const& field)
{
    if(field.size()<1)
        FatalErrorInFunction<<"Error"<<exit(FatalError);
    scalar minVal = field[0];
    scalar maxVal = field[0];
    for(label i=0; i<field.size(); i++)
    {
        minVal = std::min(minVal,field[i]);
        maxVal = std::max(maxVal,field[i]);
    }
    Info<<"  Min: "<<minVal<<" Max: "<<maxVal<<Foam::nl;
}

void Foam::printMinMax(volVectorField const& field)
{
    if(field.size()<1)
        FatalErrorInFunction<<"Error"<<exit(FatalError);
    vector minVal = field[0];
    vector maxVal = field[0];
    for(label i=0; i<field.size(); i++)
    {
        for(label d=0; d<3; d++)
        {
            minVal[d] = std::min(minVal[d],field[i][d]);
            maxVal[d] = std::max(maxVal[d],field[i][d]);
        }
    }
    Info<<"  Min: "<<minVal<<" Max: "<<maxVal<<Foam::nl;
}
