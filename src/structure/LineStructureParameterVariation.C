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
{
    Info<<"Create Parameter Variation"<<Foam::endl;
    Info<<para;
    Info<<"Create Parameter Variation done"<<Foam::endl;
}

Foam::ParameterVariation::~ParameterVariation()
{
    Info<<"Clear Parameter Variation"<<Foam::endl;
    Info<<para;
    reset();
    Info<<"Clear Parameter Variation done"<<Foam::endl;
}

const Foam::Parameter& Foam::ParameterVariation::getParameter()
{
    return para;
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
    {
        Info<<"Remove ";
        Info<<variator->getParameter()<<Foam::endl;
        delete variator;
    }
    variator = new ParameterVariation(structure,para);
    return *variator;
}

void Foam::ParameterVariation::reset()
{
    Info<<"Reset parameter";
    Info<<para<<" from "<<structure->getParameterValue(para)<<" to "<<values<<Foam::endl;
    structure->setParameterValue(para,values);
}

void Foam::ParameterVariation::vary(scalar epsilon)
{
    List<scalar> new_values = values;
    for(scalar& val : new_values)
        val += epsilon;
    
    Info<<para<<Foam::endl;
    Info<<"Vary from "<<structure->getParameterValue(para)<<" to aim "<<new_values<<Foam::endl;
    structure->setParameterValue(para,new_values);
    Info<<" and is "<<structure->getParameterValue(para)<<Foam::endl;
}
