#include "LineStructure.H"

Foam::ParameterVariation* Foam::ParameterVariation::variator = nullptr;;

Foam::ParameterVariation::ParameterVariation
(
    LineStructure* structure,
    const Parameter& para
): 
structure(structure),
para(para),
values(structure->getParameterValue(para))
{}

Foam::ParameterVariation::~ParameterVariation()
{
    reset();
}

Foam::ParameterVariation& Foam::ParameterVariation::getParameterVariator()
{
    if(variator)
        return *variator;
    else
        FatalErrorInFunction<<"ParameterVariation not set"<<exit(FatalError);
    return *variator;
}

Foam::ParameterVariation& Foam::ParameterVariation::createParameterVariator
(
    LineStructure* structure,
    const Parameter& para
)
{
    if(variator!=nullptr)
        delete variator;
    variator = new ParameterVariation(structure,para);
    return *variator;
}

void Foam::ParameterVariation::reset()
{
    structure->setParameterValue(para,values);
}

void Foam::ParameterVariation::vary(scalar epsilon)
{
    List<scalar> new_values = values;
    for(scalar& val : new_values)
        val += epsilon;
    Info<<"Vary from "<<structure->getParameterValue(para)<<Foam::endl;
    Info<<para<<Foam::endl;
    
    structure->setParameterValue(para,new_values);
    
    Info<<" to "<<structure->getParameterValue(para)<<Foam::endl;
}
