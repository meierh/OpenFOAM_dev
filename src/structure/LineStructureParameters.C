#include "LineStructure.H"

Foam::NurbsCoeffReference::NurbsCoeffReference
(
    label rodNumber,
    label coeffNumber,
    label dimension
):
rodNumber(rodNumber),
coeffNumber(coeffNumber),
dimension(dimension)
{}

Foam::Parameter::Parameter
():
parameterType(Type::None),
dimension(-1)
{}

Foam::Parameter::Parameter
(
    NurbsCoeffReference coeff
):
parameterType(Type::Rod),
dimension(-1),
valid(true)
{
    dimension = coeff.dimension;
    nurbsCoeffs.push_back(coeff);
}


Foam::Parameter::Parameter
(
    const std::vector<NurbsCoeffReference>& coeffs
):
nurbsCoeffs(coeffs),
dimension(-1),
valid(true)
{
    for(const NurbsCoeffReference& coeffRef : nurbsCoeffs)
    {
        if(dimension==-1)
            dimension = coeffRef.dimension;
        else if(dimension!=coeffRef.dimension)
            FatalErrorInFunction<<"Coefficients of different dimensions can not be one parameter!"<<exit(FatalError);
    }
}

void Foam::Parameter::addCoeff
(
    NurbsCoeffReference coeffRef
)
{
    valid=true;
    if(parameterType==Type::None)
        parameterType=Type::Rod;
    else if(parameterType==Type::CrossSection)
        FatalErrorInFunction<<"Can not add nurbs references to cross section parameter!"<<exit(FatalError);
    if(dimension==-1)
        dimension = coeffRef.dimension;
    else if(dimension!=coeffRef.dimension)
        FatalErrorInFunction<<"Coefficients of different dimensions can not be one parameter!"<<exit(FatalError);
    nurbsCoeffs.push_back(coeffRef);
}

Foam::Parameter::Parameter
(
    CrossSectionCoeffReference coeff
):
parameterType(Type::CrossSection),
dimension(-1),
valid(true)
{
    crossSecCoeffs.push_back(coeff);
}

Foam::Parameter::Parameter
(
    const std::vector<CrossSectionCoeffReference>& coeffs
):
parameterType(Type::CrossSection),
dimension(-1),
crossSecCoeffs(coeffs),
valid(true)
{}

void Foam::Parameter::addCoeff
(
    CrossSectionCoeffReference coeffRef
)
{
    valid = true;
    if(parameterType==Type::None)
        parameterType=Type::CrossSection;
    else if(parameterType==Type::Rod)
        FatalErrorInFunction<<"Can not add cross section references to nurbs parameter!"<<exit(FatalError);
    crossSecCoeffs.push_back(coeffRef);
}

std::string Foam::Parameter::to_string() const
{
    std::string result = "Parameter type:";
    switch (parameterType)
    {
        case Rod:
        {
            result+="Rod";
            result+=" / dim:"+std::to_string(dimension)+" / ";
            for(const NurbsCoeffReference& nurbsCoef : nurbsCoeffs)
                result+="("+std::to_string(nurbsCoef.rodNumber)+","+std::to_string(nurbsCoef.coeffNumber)+","+std::to_string(nurbsCoef.dimension)+") ";
            result+="/ valid:"+std::to_string(valid);
            break;
        }
        case CrossSection:
        {
            result+="CrossSection";
            result+=" / ";
            for(const CrossSectionCoeffReference& crossSecCoef : crossSecCoeffs)
                result+="("+std::to_string(crossSecCoef.rodNumber)+","+std::to_string(crossSecCoef.phase)+","+std::to_string(crossSecCoef.fourierCoeffNumber)+","+std::to_string(crossSecCoef.coeffNumber)+") ";
            result+="/ valid:"+std::to_string(valid);
            break;
        }
        case None:
        {
            result+="None";        
            break;
        }
        default:
            FatalErrorInFunction<<"Invalid"<<exit(FatalError);
    }
    return result;
}

void Foam::LineStructureParameters::collectParameters
(
    const LineStructure* structure
)
{
    for(label rodNumber=0; rodNumber<structure->getNumberRods(); rodNumber++)
    {
        label nbrCoeffs = structure->numberCurveCoeffs(rodNumber);
        for(label coeffInd=0; coeffInd<nbrCoeffs; coeffInd++)
        {
            std::array<Parameter,3> thriDimPara;
            for(label dim=0; dim<3; dim++)
            {
                Parameter para(NurbsCoeffReference(rodNumber,coeffInd,dim));
                thriDimPara[dim] = para;
            }
            threeDimParameters.push_back(thriDimPara);
        }
    }
}
