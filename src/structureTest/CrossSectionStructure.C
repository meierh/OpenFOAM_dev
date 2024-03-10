#include "CrossSectionStructure.H"

Foam::CrossSection::CrossSection
(
    std::vector<gsNurbs<double>> sinusCoeffs,
    std::vector<gsNurbs<double>> cosinusCoeffs
):
a_k(sinusCoeffs),
b_k(cosinusCoeffs),
numberCoeffs(sinusCoeffs.size())
{
    if(b_k.size()!=numberCoeffs)
        FatalErrorInFunction<<"Invalid coefficents of Cross Section"<<exit(FatalError);
    if(numberCoeffs<1)
        FatalErrorInFunction<<"Number of coefficients must be larger than 0"<<exit(FatalError);
    
    for(const gsNurbs<double>& oneCoeff : a_k)
        if(oneCoeff.domainDim()!=1)
            FatalErrorInFunction<<"Coeff Nurbs must have a Basis of dimension 1"<<exit(FatalError);
    for(const gsNurbs<double>& oneCoeff : b_k)
        if(oneCoeff.domainDim()!=1)
            FatalErrorInFunction<<"Coeff Nurbs must have a Basis of dimension 1"<<exit(FatalError);
    for(const gsNurbs<double>& oneCoeff : a_k)
        if(oneCoeff.size()!=1)
            FatalErrorInFunction<<"Coeff Nurbs must have a dimension of 1"<<exit(FatalError);
    for(const gsNurbs<double>& oneCoeff : b_k)
        if(oneCoeff.size()!=1)
            FatalErrorInFunction<<"Coeff Nurbs must have a dimension of 1"<<exit(FatalError);
}

std::function<double(double)> Foam::CrossSection::getEvalOnPoint(gsMatrix<double,1,1> parameter)
{
    std::shared_ptr<std::vector<double>> aCoeffs;
    for(label coeffI=0; coeffI<a_k.size(); coeffI++)
    {
        gsMatrix<double> aCoeffI = a_k[coeffI].eval(parameter);
        aCoeffs->push_back(aCoeffI.at(0));
    }
    std::shared_ptr<std::vector<double>> bCoeffs;
    for(label coeffI=0; coeffI<b_k.size(); coeffI++)
    {
        gsMatrix<double> bCoeffI = b_k[coeffI].eval(parameter);
        bCoeffs->push_back(bCoeffI.at(0));
    }
    return [num_Coeff=numberCoeffs,a_k=aCoeffs,b_k=bCoeffs](double rad)
    {
        double value = (*a_k)[0]/2;
        for(label coeffI=1; coeffI<num_Coeff; coeffI++)
        {
            value += (*a_k)[coeffI]*cos(coeffI*rad)+(*b_k)[coeffI]*sin(coeffI*rad);
        }
        return value;
    };
}

double Foam::CrossSection::operator()(gsMatrix<double,1,1> parameter,double rad)
{
    std::function<double(double)> evalOnPoint = getEvalOnPoint(parameter);
    return evalOnPoint(rad);
}

Foam::CrossSectionStructure::CrossSectionStructure
(
    Time& runTime,
    const dimensionedScalar& alpha,
    const volScalarField& T,
    const volScalarField& p,
    const volVectorField& U,
    dynamicRefineFvMesh& mesh,
    const dimensionedScalar nu,
    std::vector<CrossSection> rodCrossSection
):
LineStructure(runTime,alpha,T,p,U,mesh,nu),
rodCrossSection(rodCrossSection)
{
}
