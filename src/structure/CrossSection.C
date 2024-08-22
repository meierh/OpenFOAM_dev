#include "CrossSectionStructure.H"

Foam::CrossSection::CrossSection
(
    scalar radius
)
{
    std::tuple<gsNurbs<scalar>,std::vector<gsNurbs<scalar>>,std::vector<gsNurbs<scalar>>,gsNurbs<scalar>> nurbs;
    nurbs = constCrossSec(2*radius,{},{});
    init(std::get<0>(nurbs),std::get<1>(nurbs),std::get<2>(nurbs),std::get<3>(nurbs));
}

Foam::CrossSection::CrossSection
(
    scalar a_0,
    std::vector<scalar> a_k,
    std::vector<scalar> b_k,
    scalar phaseShift
)        
{
    std::tuple<gsNurbs<scalar>,std::vector<gsNurbs<scalar>>,std::vector<gsNurbs<scalar>>,gsNurbs<scalar>> nurbs;
    nurbs = constCrossSec(a_0,a_k,b_k,phaseShift);
    init(std::get<0>(nurbs),std::get<1>(nurbs),std::get<2>(nurbs),std::get<3>(nurbs));
}

Foam::CrossSection::CrossSection
(
    scalar a_0,
    std::vector<scalar> a_k,
    std::vector<scalar> b_k,
    gsNurbs<scalar> phaseShift
)        
{
    std::tuple<gsNurbs<scalar>,std::vector<gsNurbs<scalar>>,std::vector<gsNurbs<scalar>>,gsNurbs<scalar>> nurbs;
    nurbs = twistedConstCrossSec(a_0,a_k,b_k,phaseShift);
    init(std::get<0>(nurbs),std::get<1>(nurbs),std::get<2>(nurbs),std::get<3>(nurbs));
}

Foam::CrossSection::CrossSection
(
    gsNurbs<scalar> a_0,
    std::vector<gsNurbs<scalar>> a_k,
    std::vector<gsNurbs<scalar>> b_k,
    gsNurbs<scalar> phaseShift
)
{
    init(a_0,a_k,b_k,phaseShift);
}

std::tuple<gsNurbs<Foam::scalar>,std::vector<gsNurbs<Foam::scalar>>,std::vector<gsNurbs<Foam::scalar>>,gsNurbs<Foam::scalar>> Foam::CrossSection::constCrossSec
(
    scalar a_0,
    std::vector<scalar> a_k,
    std::vector<scalar> b_k,
    scalar phaseShift
)
{        
    gsNurbs<scalar> n_a0 = createConstNurbs(a_0);
    
    std::vector<gsNurbs<scalar>> n_ak;
    for(scalar a : a_k)
        n_ak.push_back(createConstNurbs(a));
    
    std::vector<gsNurbs<scalar>> n_bk;
    for(scalar b : b_k)
        n_bk.push_back(createConstNurbs(b));
    
    gsNurbs<scalar> pS = createConstNurbs(phaseShift);
    
    return std::make_tuple(n_a0,n_ak,n_bk,pS);
}

std::tuple<gsNurbs<Foam::scalar>,std::vector<gsNurbs<Foam::scalar>>,std::vector<gsNurbs<Foam::scalar>>,gsNurbs<Foam::scalar>> Foam::CrossSection::twistedConstCrossSec
(
    scalar a_0,
    std::vector<scalar> a_k,
    std::vector<scalar> b_k,
    gsNurbs<scalar> phaseShift
)
{       
    gsNurbs<scalar> n_a0 = createConstNurbs(a_0);
    
    std::vector<gsNurbs<scalar>> n_ak;
    for(scalar a : a_k)
        n_ak.push_back(createConstNurbs(a));
    
    std::vector<gsNurbs<scalar>> n_bk;
    for(scalar b : b_k)
        n_bk.push_back(createConstNurbs(b));

    return std::make_tuple(n_a0,n_ak,n_bk,phaseShift);
}

void Foam::CrossSection::init
(
    gsNurbs<scalar> a_0,
    std::vector<gsNurbs<scalar>> a_k,
    std::vector<gsNurbs<scalar>> b_k,
    gsNurbs<scalar> phaseShift
)
{
    this->a_0 = a_0;
    this->a_k = a_k;
    this->b_k = b_k;
    this->phaseShift = phaseShift;
    
    if(this->a_k.size() != this->b_k.size())
        FatalErrorInFunction<<"Coefficient number mismatch!"<<exit(FatalError);
    numberCoeffs = a_k.size();
    
    domStart = a_0.domainStart();
    domEnd = a_0.domainEnd();
    
    for(const gsNurbs<scalar>& oneCoeff : a_k)
        if(oneCoeff.domainDim()!=1)
            FatalErrorInFunction<<"Coeff Nurbs must have a Basis of dimension 1"<<exit(FatalError);
    for(const gsNurbs<scalar>& oneCoeff : b_k)
        if(oneCoeff.domainDim()!=1)
            FatalErrorInFunction<<"Coeff Nurbs must have a Basis of dimension 1"<<exit(FatalError);
    if(phaseShift.domainDim()!=1)
        FatalErrorInFunction<<"Coeff Nurbs must have a Basis of dimension 1"<<exit(FatalError);
    
    for(const gsNurbs<scalar>& oneCoeff : a_k)
        if(oneCoeff.size()!=1)
            FatalErrorInFunction<<"Coeff Nurbs must have a dimension of 1"<<exit(FatalError);
    for(const gsNurbs<scalar>& oneCoeff : b_k)
        if(oneCoeff.size()!=1)
            FatalErrorInFunction<<"Coeff Nurbs must have a dimension of 1"<<exit(FatalError);
    if(phaseShift.size()!=1)
        FatalErrorInFunction<<"Coeff Nurbs must have a Basis of dimension 1"<<exit(FatalError);
    
    for(const gsNurbs<scalar>& oneCoeff : a_k)
    {
        if(oneCoeff.domainStart()!=domStart)
            FatalErrorInFunction<<"Mismatch in domain range"<<exit(FatalError);
        if(oneCoeff.domainEnd()!=domEnd)
            FatalErrorInFunction<<"Mismatch in domain range"<<exit(FatalError);        
    }
    for(const gsNurbs<scalar>& oneCoeff : b_k)
    {
        if(oneCoeff.domainStart()!=domStart)
            FatalErrorInFunction<<"Mismatch in domain range"<<exit(FatalError);
        if(oneCoeff.domainEnd()!=domEnd)
            FatalErrorInFunction<<"Mismatch in domain range"<<exit(FatalError);
    }
    if(phaseShift.domainStart()!=domStart)
        FatalErrorInFunction<<"Mismatch in domain range"<<exit(FatalError);
    if(phaseShift.domainEnd()!=domEnd)
        FatalErrorInFunction<<"Mismatch in domain range"<<exit(FatalError);
}

Foam::scalar Foam::CrossSection::operator()(scalar parameter,scalar angle) const
{
    //Info<<"par:"<<parameter<<" angle:"<<angle<<Foam::endl;
    std::function<scalar(scalar)> evalOnPoint = getEvalOnPoint(parameter);
    return evalOnPoint(angle);
}

std::function<Foam::scalar(Foam::scalar)> Foam::CrossSection::getEvalOnPoint(scalar parameter) const
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    //Info<<"parameter:"<<parameter<<Foam::endl;

    scalar a0Coeff;
    gsMatrix<scalar> a0CoeffI = a_0.eval(parMat);
    a0Coeff = a0CoeffI.at(0);
    
    //Info<<"a0Coeff:"<<a0Coeff<<Foam::endl;
    
    auto aCoeffs = std::make_shared<std::vector<scalar>>();
    for(uint coeffI=0; coeffI<a_k.size(); coeffI++)
    {
        gsMatrix<scalar> aCoeffI = a_k[coeffI].eval(parMat);
        aCoeffs->push_back(aCoeffI.at(0));
    }
    
    auto bCoeffs = std::make_shared<std::vector<scalar>>();
    for(uint coeffI=0; coeffI<b_k.size(); coeffI++)
    {
        gsMatrix<scalar> bCoeffI = b_k[coeffI].eval(parMat);
        bCoeffs->push_back(bCoeffI.at(0));
    }
    
    scalar phShift;
    gsMatrix<scalar> phaseShiftM = phaseShift.eval(parMat);
    phShift = phaseShiftM.at(0);
    
    //Info<<"parameter:"<<parameter<<Foam::endl;
    return [num_Coeff=numberCoeffs,a_0=a0Coeff,a_k=aCoeffs,b_k=bCoeffs,pS=phShift](scalar rad)
    {
        scalar value = a_0/2;
        for(uint coeffI=0; coeffI<num_Coeff; coeffI++)
        {
            label k=coeffI+1;
            value += (*a_k)[coeffI]*std::cos(k*rad+pS);
            value += (*b_k)[coeffI]*std::sin(k*rad+pS);
        }
        return value;
    };
}

Foam::scalar Foam::CrossSection::deriv_angle(scalar parameter,scalar angle) const
{
    //Info<<"par:"<<parameter<<" angle:"<<angle<<Foam::endl;
    std::function<scalar(scalar)> derivOnPoint = getDerivAngleOnPoint(parameter);
    return derivOnPoint(angle);
}

std::function<Foam::scalar(Foam::scalar)> Foam::CrossSection::getDerivAngleOnPoint(scalar parameter) const
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    //Info<<"parameter:"<<parameter<<Foam::endl;

    scalar a0Coeff;
    gsMatrix<scalar> a0CoeffI = a_0.eval(parMat);
    a0Coeff = a0CoeffI.at(0);
    
    //Info<<"a0Coeff:"<<a0Coeff<<Foam::endl;
    
    auto aCoeffs = std::make_shared<std::vector<scalar>>();
    for(uint coeffI=0; coeffI<a_k.size(); coeffI++)
    {
        gsMatrix<scalar> aCoeffI = a_k[coeffI].eval(parMat);
        aCoeffs->push_back(aCoeffI.at(0));
    }
    
    auto bCoeffs = std::make_shared<std::vector<scalar>>();
    for(uint coeffI=0; coeffI<b_k.size(); coeffI++)
    {
        gsMatrix<scalar> bCoeffI = b_k[coeffI].eval(parMat);
        bCoeffs->push_back(bCoeffI.at(0));
    }
    
    scalar phShift;
    gsMatrix<scalar> phaseShiftM = phaseShift.eval(parMat);
    phShift = phaseShiftM.at(0);
    
    //Info<<"parameter:"<<parameter<<Foam::endl;
    return [num_Coeff=numberCoeffs,a_0=a0Coeff,a_k=aCoeffs,b_k=bCoeffs,pS=phShift](scalar angle)
    {
        scalar value = 0;
        for(uint coeffI=0; coeffI<num_Coeff; coeffI++)
        {
            label k=coeffI+1;
            value += (*a_k)[coeffI]*std::sin(k*angle+pS)*(-k);
            value += (*b_k)[coeffI]*std::cos(k*angle+pS)*k;
        }
        return value;
    };
}

void Foam::CrossSection::delete_coeffDerivedCurves()
{
    for(std::vector<gsNurbs<scalar>*>* coeffDerivedCurve : coeffDerivedCurves)
    {
        if(coeffDerivedCurve!=nullptr)
        {
            for(gsNurbs<scalar>* coeffOneCoeffDerivedCurve : *coeffDerivedCurve)
            {
                if(coeffOneCoeffDerivedCurve!=nullptr)
                {
                    delete coeffOneCoeffDerivedCurve;
                }
            }
            delete coeffDerivedCurve;
        }
    }
};


const gsNurbs<Foam::scalar>* Foam::CrossSection::getCurvePtr
(
    label fourierCoeffNumber
) const
{
    const gsNurbs<scalar>* curve;
    if(fourierCoeffNumber==0)
    {
        curve = &a_0;
    }
    else
    {
        fourierCoeffNumber -= 1;
        if(fourierCoeffNumber%2==0)
        {
            fourierCoeffNumber /= 2;
            curve = &(a_k[fourierCoeffNumber]);
        }
        else
        {
            fourierCoeffNumber -= 1;
            fourierCoeffNumber /= 2;
            curve = &(b_k[fourierCoeffNumber]);
        }
    }
    return curve;
}

gsNurbs<Foam::scalar>* Foam::CrossSection::getCurvePtr
(
    label fourierCoeffNumber
)
{
    gsNurbs<scalar>* curve;
    if(fourierCoeffNumber==0)
    {
        curve = &a_0;
    }
    else
    {
        fourierCoeffNumber -= 1;
        if(fourierCoeffNumber%2==0)
        {
            fourierCoeffNumber /= 2;
            curve = &(a_k[fourierCoeffNumber]);
        }
        else
        {
            fourierCoeffNumber -= 1;
            fourierCoeffNumber /= 2;
            curve = &(b_k[fourierCoeffNumber]);
        }
    }
    return curve;
}

Foam::label Foam::CrossSection::numberFourierCoeff() const
{
    return numberCoeffs+1;
}

Foam::label Foam::CrossSection::numberNurbsCoeffs
(
    label fourierCoeffNumber
) const
{
    const gsNurbs<scalar>* curve = getCurvePtr(fourierCoeffNumber);
    return curve->coefs().cols();
}

Foam::scalar Foam::CrossSection::evalRadiusDerivCoeff
(
    label fourierCoeffNumber,
    label derivCoeffNumber,
    scalar parameter
)
{
    std::vector<gsNurbs<scalar>*>* coeffDerivedCurve = coeffDerivedCurves[fourierCoeffNumber];
    if(coeffDerivedCurve==nullptr)
    {
        coeffDerivedCurve = new std::vector<gsNurbs<scalar>*>();
        coeffDerivedCurve->resize(numberNurbsCoeffs(fourierCoeffNumber));
    }
    gsNurbs<scalar>* coeffOneCoeffDerivedCurve = (*coeffDerivedCurve)[derivCoeffNumber];
    if(coeffOneCoeffDerivedCurve==nullptr)
    {
        coeffOneCoeffDerivedCurve = new gsNurbs<scalar>();
        *coeffOneCoeffDerivedCurve = *(getCurvePtr(fourierCoeffNumber));
        gsMatrix<scalar>& coeffs = coeffOneCoeffDerivedCurve->coefs();
        for(label col=0; col<coeffs.cols(); col++)
        {
            scalar replVal = col==derivCoeffNumber?1:0;
            for(label row=0; row<coeffs.rows(); row++)
            {
                coeffs(col,row) = replVal;
            }
        }
    }
    
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    const gsNurbs<scalar>& coeffDerivCurve = *coeffOneCoeffDerivedCurve;
    gsMatrix<scalar> coeffDerivEval;
    coeffDerivCurve.eval_into(parMat,coeffDerivEval);    
    return coeffDerivEval(0,0);
}

void Foam::CrossSection::setNurbsCoeff
(
    label fourierCoeffNumber,
    label derivCoeffNumber,
    scalar value
)
{
    if(fourierCoeffNumber<0 || fourierCoeffNumber>=numberFourierCoeff())
        FatalErrorInFunction<<"Invalid fourierCoeffNumber"<<exit(FatalError);
    gsNurbs<scalar>* curve = getCurvePtr(fourierCoeffNumber);
    gsMatrix<scalar>& coeffs = curve->coefs();
    if(derivCoeffNumber<0 || derivCoeffNumber>=coeffs.cols())
        FatalErrorInFunction<<"Invalid derivCoeffNumber"<<exit(FatalError);
    coeffs(derivCoeffNumber,0) = value;
    
    std::vector<gsNurbs<scalar>*>* coeffDerivedCurve = coeffDerivedCurves[fourierCoeffNumber];
    if(coeffDerivedCurve!=nullptr)
    {
        for(gsNurbs<scalar>* coeffOneCoeffDerivedCurve : *coeffDerivedCurve)
        {
            if(coeffOneCoeffDerivedCurve!=nullptr)
            {
                delete coeffOneCoeffDerivedCurve;
            }
        }
        delete coeffDerivedCurve;
    }
}

Foam::scalar Foam::CrossSection::lowerLimitRadius(scalar parameter) const
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    scalar a0Coeff;
    gsMatrix<scalar> a0CoeffI = a_0.eval(parMat);
    a0Coeff = a0CoeffI.at(0);
        
    auto aCoeffs = std::make_shared<std::vector<scalar>>();
    for(uint coeffI=0; coeffI<a_k.size(); coeffI++)
    {
        gsMatrix<scalar> aCoeffI = a_k[coeffI].eval(parMat);
        aCoeffs->push_back(aCoeffI.at(0));
    }
    
    auto bCoeffs = std::make_shared<std::vector<scalar>>();
    for(uint coeffI=0; coeffI<b_k.size(); coeffI++)
    {
        gsMatrix<scalar> bCoeffI = b_k[coeffI].eval(parMat);
        bCoeffs->push_back(bCoeffI.at(0));
    }
    
    scalar minValue = a0Coeff/2;
    for(scalar ak : *aCoeffs)
        minValue -= std::abs(ak);
    for(scalar bk : *bCoeffs)
        minValue -= std::abs(bk);
    return minValue;
}

Foam::scalar Foam::CrossSection::upperLimitRadius(scalar parameter) const
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    scalar a0Coeff;
    gsMatrix<scalar> a0CoeffI = a_0.eval(parMat);
    a0Coeff = a0CoeffI.at(0);
        
    auto aCoeffs = std::make_shared<std::vector<scalar>>();
    for(uint coeffI=0; coeffI<a_k.size(); coeffI++)
    {
        gsMatrix<scalar> aCoeffI = a_k[coeffI].eval(parMat);
        aCoeffs->push_back(aCoeffI.at(0));
    }
    
    auto bCoeffs = std::make_shared<std::vector<scalar>>();
    for(uint coeffI=0; coeffI<b_k.size(); coeffI++)
    {
        gsMatrix<scalar> bCoeffI = b_k[coeffI].eval(parMat);
        bCoeffs->push_back(bCoeffI.at(0));
    }
    
    scalar maxValue = a0Coeff/2;
    for(scalar ak : *aCoeffs)
        maxValue += std::abs(ak);
    for(scalar bk : *bCoeffs)
        maxValue += std::abs(bk);
    return maxValue;
}

std::pair<Foam::scalar,Foam::scalar> Foam::CrossSection::nurbsBounds
(
    gsNurbs<scalar> curve
) const
{
    if(curve.coefs().cols()<1)
        FatalErrorInFunction<<"Size mismatch"<<exit(FatalError);
    
    scalar min = curve.coefs()(0,0);
    scalar max = curve.coefs()(0,0);
    for(label i=0; i<curve.coefs().cols(); i++)
    {
        min = std::min(min,curve.coefs()(i,0));
        max = std::max(max,curve.coefs()(i,0));
    }
    return {min,max};
}

std::pair<Foam::scalar,Foam::scalar> Foam::CrossSection::nurbsBounds
(
    gsNurbs<scalar> curve,
    scalar start,
    scalar end
) const
{
    std::unordered_set<label> knotSet;
    for(scalar knot : curve.knots())
        knotSet.insert(knot);

    if(knotSet.find(start)==knotSet.end())
        curve.insertKnot(start);
    if(knotSet.find(end)==knotSet.end())
        curve.insertKnot(end);

    label degree = curve.knots().degree();
    label knot_start = -1;
    label knot_end = -1;
    for(std::size_t knotI=0; knotI<curve.knots().size()-1; knotI++)
    {
        scalar knot_i0 = curve.knots()[knotI];
        scalar knot_i1 = curve.knots()[knotI+1];
        
        if(knot_i0==start && knot_i0!=knot_i1)
        {
            if(knot_start!=-1)
                FatalErrorInFunction<<"Duplicate assigned"<<exit(FatalError);
            knot_start=static_cast<label>(knotI);
        }
        if(knot_i1==end && knot_i0!=knot_i1)
        {
            if(knot_end!=-1)
                FatalErrorInFunction<<"Duplicate assigned"<<exit(FatalError);
            knot_end=static_cast<label>(knotI+1);
        }
    }
    if(knot_start==-1 || knot_end==-1)
        FatalErrorInFunction<<"Not assigned"<<exit(FatalError);

    label coeff_start = knot_start-degree;
    label coeff_end = knot_end-1;

    if(curve.coefs().rows()!=1)
        FatalErrorInFunction<<"Rows number out of range"<<exit(FatalError);
    if(coeff_start<0 || coeff_start>=curve.coefs().cols())
        FatalErrorInFunction<<"Coeff start out of range"<<exit(FatalError);
    if(coeff_end<0 || coeff_end>=curve.coefs().cols())
        FatalErrorInFunction<<"Coeff start out of range"<<exit(FatalError);

    gsMatrix<scalar> coeffs(coeff_end-coeff_start+1,1);
    label ind=0;
    for(label c_s=coeff_start; c_s<coeff_end+1; c_s++,ind++)
    {
        coeffs(ind,0) = curve.coefs()(c_s,0);
    }
    if(ind!=coeffs.cols())
        FatalErrorInFunction<<"Size mismatch"<<exit(FatalError);
    
    scalar min = coeffs(0,0);
    scalar max = coeffs(0,0);
    for(label i=0; i<coeffs.cols(); i++)
    {
        min = std::min(min,coeffs(i,0));
        max = std::max(max,coeffs(i,0));
    }
    return {min,max};
}

std::pair<Foam::scalar,Foam::scalar> Foam::CrossSection::radiusBounds() const
{
    auto minMax = nurbsBounds(a_0);
    scalar min = minMax.first/2;
    scalar max = minMax.second/2;
    
    for(const gsNurbs<scalar>& ak : a_k)
    {
        auto minMax = nurbsBounds(ak);
        min += minMax.first;
        max += minMax.second;
    }
    for(const gsNurbs<scalar>& bk : b_k)
    {
        auto minMax = nurbsBounds(bk);
        min += minMax.first;
        max += minMax.second;
    }
    return {min,max};
}

std::pair<Foam::scalar,Foam::scalar> Foam::CrossSection::radiusBounds
(
    scalar start,
    scalar end
) const
{
    auto minMax = nurbsBounds(a_0,start,end);
    scalar min = minMax.first/2;
    scalar max = minMax.second/2;
    
    for(const gsNurbs<scalar>& ak : a_k)
    {
        auto minMax = nurbsBounds(ak,start,end);
        min += minMax.first;
        max += minMax.second;
    }
    for(const gsNurbs<scalar>& bk : b_k)
    {
        auto minMax = nurbsBounds(bk,start,end);
        min += minMax.first;
        max += minMax.second;
    }
    
    return {min,max};
}

gsNurbs<Foam::scalar> Foam::CrossSection::createConstNurbs(scalar coeff) const
{
    std::vector<scalar> knotContainer = {0,0,1,1};
    gsKnotVector<scalar> cKnots(knotContainer,1);
    gsMatrix<scalar> cWeight(2,1); cWeight.at(0) = 1; cWeight.at(1) = 1;
    gsMatrix<scalar> cCoeff(2,1); cCoeff.at(0) = coeff; cCoeff.at(1) = coeff;
    return gsNurbs<scalar>(cKnots,cWeight,cCoeff);
};
