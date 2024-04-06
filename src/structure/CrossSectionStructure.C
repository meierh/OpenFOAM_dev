#include "CrossSectionStructure.H"


Foam::CrossSection::CrossSection
(
    scalar radius
)
{
    std::tuple<gsNurbs<scalar>,std::vector<gsNurbs<scalar>>,std::vector<gsNurbs<scalar>>> nurbs;
    nurbs = constCrossSec(2*radius,{},{});
    init(std::get<0>(nurbs),std::get<1>(nurbs),std::get<2>(nurbs));
}

Foam::CrossSection::CrossSection
(
    scalar a_0,
    std::vector<scalar> a_k,
    std::vector<scalar> b_k
)        
{
    std::tuple<gsNurbs<scalar>,std::vector<gsNurbs<scalar>>,std::vector<gsNurbs<scalar>>> nurbs;
    nurbs = constCrossSec(a_0,a_k,b_k);
    init(std::get<0>(nurbs),std::get<1>(nurbs),std::get<2>(nurbs));
}

Foam::CrossSection::CrossSection
(
    gsNurbs<scalar> a_0,
    std::vector<gsNurbs<scalar>> a_k,
    std::vector<gsNurbs<scalar>> b_k
)
{
    init(a_0,a_k,b_k);
}

std::tuple<gsNurbs<scalar>,std::vector<gsNurbs<scalar>>,std::vector<gsNurbs<scalar>>> Foam::CrossSection::constCrossSec
(
    scalar a_0,
    std::vector<scalar> a_k,
    std::vector<scalar> b_k
)
{
    auto createConstNurbs = [](scalar coeff)
    {
        std::vector<scalar> knotContainer = {0,0,1,1};
        gsKnotVector<scalar> cKnots(knotContainer,1);
        gsMatrix<scalar> cWeight(2,1); cWeight.at(0) = 1; cWeight.at(1) = 1;
        gsMatrix<scalar> cCoeff(2,1); cCoeff.at(0) = coeff; cCoeff.at(1) = coeff;
        return gsNurbs<scalar>(cKnots,cWeight,cCoeff);
    };
        
    gsNurbs<scalar> n_a0 = createConstNurbs(a_0);
    
    std::vector<gsNurbs<scalar>> n_ak;
    for(scalar a : a_k)
        n_ak.push_back(createConstNurbs(a));
    
    std::vector<gsNurbs<scalar>> n_bk;
    for(scalar b : b_k)
        n_bk.push_back(createConstNurbs(b));

    return std::make_tuple(n_a0,n_ak,n_bk);
}

void Foam::CrossSection::init
(
    gsNurbs<scalar> a_0,
    std::vector<gsNurbs<scalar>> a_k,
    std::vector<gsNurbs<scalar>> b_k
)
{
    this->a_0 = a_0;
    this->a_k = a_k;
    this->b_k = b_k;
    
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
    for(const gsNurbs<scalar>& oneCoeff : a_k)
        if(oneCoeff.size()!=1)
            FatalErrorInFunction<<"Coeff Nurbs must have a dimension of 1"<<exit(FatalError);
    for(const gsNurbs<scalar>& oneCoeff : b_k)
        if(oneCoeff.size()!=1)
            FatalErrorInFunction<<"Coeff Nurbs must have a dimension of 1"<<exit(FatalError);
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
}

scalar Foam::CrossSection::operator()(scalar parameter,scalar rad) const
{
    //Info<<"par:"<<parameter<<" rad:"<<rad<<Foam::endl;
    std::function<scalar(scalar)> evalOnPoint = getEvalOnPoint(parameter);
    return evalOnPoint(rad);
}

std::function<scalar(scalar)> Foam::CrossSection::getEvalOnPoint(scalar parameter) const
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    //Info<<"parameter:"<<parameter<<Foam::endl;

    scalar a0Coeff;
    gsMatrix<scalar> a0CoeffI = a_0.eval(parMat);
    a0Coeff = a0CoeffI.at(0);
    
    //Info<<"a0Coeff:"<<a0Coeff<<Foam::endl;
    
    auto aCoeffs = std::make_shared<std::vector<scalar>>();
    for(label coeffI=0; coeffI<a_k.size(); coeffI++)
    {
        gsMatrix<scalar> aCoeffI = a_k[coeffI].eval(parMat);
        aCoeffs->push_back(aCoeffI.at(0));
    }
    
    auto bCoeffs = std::make_shared<std::vector<scalar>>();
    for(label coeffI=0; coeffI<b_k.size(); coeffI++)
    {
        gsMatrix<scalar> bCoeffI = b_k[coeffI].eval(parMat);
        bCoeffs->push_back(bCoeffI.at(0));
    }
    //Info<<"parameter:"<<parameter<<Foam::endl;
    return [num_Coeff=numberCoeffs,a_0=a0Coeff,a_k=aCoeffs,b_k=bCoeffs](scalar rad)
    {
        scalar value = a_0/2;
        for(label coeffI=0; coeffI<num_Coeff; coeffI++)
        {
            label k=coeffI+1;
            value += (*a_k)[coeffI]*std::cos(k*rad)+(*b_k)[coeffI]*std::sin(k*rad);
        }
        return value;
    };
}

scalar Foam::CrossSectionStructure::distance
(
    const LagrangianMarkerOnCrossSec& A,
    const LagrangianMarkerOnCrossSec& B
)
{
    if(A.baseRod!=B.baseRod)
        FatalErrorInFunction<<"Distance can not be computed between points on different rods"<<exit(FatalError);
    
    using node = std::tuple<scalar,scalar,vector>;
    
    auto restrictAngle = [](scalar angle)
    {
        while(angle<0)
        {
            angle += 2*Foam::constant::mathematical::pi;
        }
        while(angle>2*Foam::constant::mathematical::pi)
        {
            angle -= 2*Foam::constant::mathematical::pi;
        }
        return angle;
    };
    
    node start = {A.markerParameter,A.markerAngle,A.markerPosition};
    std::get<1>(start) = restrictAngle(std::get<1>(start));
    node end = {B.markerParameter,B.markerAngle,B.markerPosition};
    std::get<1>(end) = restrictAngle(std::get<1>(end));
    
    std::list<node> nodes;
    nodes.push_back(start);
    nodes.push_back(end);
    
    scalar distance = std::numeric_limits<scalar>::max();
    bool errorSufficient = false;
    while(!errorSufficient)
    {
        auto node0 = nodes.begin();
        auto node1 = ++nodes.begin();
        for( ; node1!=nodes.end() ; )
        {
            node insertNode;
            std::get<0>(insertNode) = (std::get<0>(*node0)+std::get<0>(*node1))/2;
            
            scalar anglNode0 = std::get<1>(*node0);
            scalar anglNode1 = std::get<1>(*node1);
            scalar anglNodeInter;
            scalar distInner = std::abs(anglNode0-anglNode1);
            if(anglNode0<anglNode1)
            {
                scalar distToUp = 2*Foam::constant::mathematical::pi-anglNode1;
                scalar outerDist = distToUp+anglNode0;
                if(outerDist>distInner)
                    anglNodeInter = (anglNode0+anglNode1)/2;
                else
                {
                    scalar halfOuterDist = outerDist/2;
                    if(distToUp<anglNode0)
                        anglNodeInter = anglNode0-halfOuterDist;
                    else
                        anglNodeInter = anglNode1+halfOuterDist;
                }
            }
            else
            {
                scalar distToUp = 2*Foam::constant::mathematical::pi-anglNode0;
                scalar outerDist = distToUp+anglNode1;
                if(outerDist>distInner)
                    anglNodeInter = (anglNode0+anglNode1)/2;
                else
                {
                    scalar halfOuterDist = outerDist/2;
                    if(distToUp<anglNode1)
                        anglNodeInter = anglNode1-halfOuterDist;
                    else
                        anglNodeInter = anglNode0+halfOuterDist;
                }
            }
            std::get<1>(insertNode) = anglNodeInter;
            std::get<2>(insertNode) = evaluateRodCircumPos(A.baseRod,std::get<0>(insertNode),
                                                           A.baseCrossSec,std::get<1>(insertNode));
            auto inserted = nodes.insert(node1,insertNode);
            node0 = node1;
            node1++;
        }
        
        scalar newDist = 0;
        node0 = nodes.begin();
        node1 = ++nodes.begin();
        for( ; node1!=nodes.end() ; )
        {
            vector connec = std::get<2>(*node1) - std::get<2>(*node0);
            newDist += std::sqrt(connec&connec);
            node0++;
            node1++;
        }
        
        scalar change = std::abs(newDist-distance);
        scalar changePerc = change/newDist;
        distance = newDist;
        if(nodes.size()>10 && changePerc<1e-3)
            errorSufficient=true;
    }
    return distance;
}

Foam::CrossSectionStructure::CrossSectionStructure
(
    dynamicRefineFvMesh& mesh,
    const dimensionedScalar& alpha,
    const volScalarField& T,
    const volScalarField& p,
    const volVectorField& U,
    const dimensionedScalar nu,
    std::vector<CrossSection> rodCrossSection
):
LineStructure(mesh,alpha,T,p,U,nu,{}),
rodCrossSection(rodCrossSection)
{
}

void Foam::CrossSectionStructure::transferMarkers(FieldMarkerStructureInteraction& connector)
{
    Info<<"Transfer markers"<<Foam::endl;
    if(!myMesh)
        FatalErrorInFunction<<"Rod Mesh not set!"<<exit(FatalError);
    if(myMesh->m_nR!=rodMarkers.size())
    {
        rodMarkers.resize(myMesh->m_nR);
    }
    connector.markers.resize(0);
    
    if(rodCrossSection.size()!=rodMarkers.size())
        FatalErrorInFunction<<"Mismatch in size of rodCrossSection and rodMarkers"<<exit(FatalError);
    
    for(uint rodIndex=0; rodIndex<rodMarkers.size(); rodIndex++)
    {
        Info<<"  Transfer markers from "<<rodIndex<<Foam::endl;
        std::unique_ptr<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>& oneRodMarkers = rodMarkers[rodIndex];
        if(!oneRodMarkers)
        {
            myMesh->m_Rods[rodIndex]->m_Rot.setCoefs(myMesh->m_Rods[rodIndex]->m_Rot_coefs.transpose());
            oneRodMarkers = constructMarkerSet(myMesh->m_Rods[rodIndex],rodCrossSection[rodIndex]);
        }
        
        for(std::vector<LagrangianMarkerOnCrossSec>& circMarker : *oneRodMarkers)
        {
            for(LagrangianMarkerOnCrossSec& marker : circMarker)
            {
                connector.markers.push_back(&marker);
            }
        }
    }
    Info<<"Transfer markers done"<<Foam::endl;
}

Foam::LagrangianMarkerOnCrossSec Foam::CrossSectionStructure::createLagrangianMarkerOnCrossSec
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    const CrossSection& oneCrossSec,
    scalar angle
)
{
    LagrangianMarkerOnCrossSec marker
    (
        parameter,
        evaluateRodCircumPos(oneRod,parameter,oneCrossSec,angle),
        oneRod,
        angle,
        oneCrossSec
    );
    label cellOfMarker = mesh.findCell(marker.markerPosition);
    if(cellOfMarker!=-1)
        getSupportDomain(cellOfMarker,marker.supportCells);
    marker.markerCell = -1;
    for(auto iter=marker.supportCells.begin(); iter!=marker.supportCells.end(); iter++)
        if(mesh.pointInCell(marker.markerPosition,*iter))
        {
            if(marker.markerCell!=-1)
                FatalErrorInFunction<<"Double assignment marker cell"<<exit(FatalError);
            marker.markerCell = *iter;
        }       
    return marker;
}

Foam::scalar Foam::CrossSectionStructure::evaluateCircumArcLen
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameterA,
    scalar parameterB,
    const CrossSection& oneCrossSec,
    scalar angleA,
    scalar angleB
)
{
    vector parAVec = evaluateRodCircumPos(oneRod,parameterA,oneCrossSec,angleA);
    vector parBVec = evaluateRodCircumPos(oneRod,parameterB,oneCrossSec,angleB);
    vector connec = parAVec-parBVec;
    return Foam::mag(connec);
}

Foam::vector Foam::CrossSectionStructure::evaluateRodCircumPos
(
    const ActiveRodMesh::rodCosserat* oneRod,
    Foam::scalar parameter,
    const CrossSection& oneCrossSec,
    Foam::scalar angle
)
{
    vector d1,d2,d3,r;
    rodEval(oneRod,parameter,d1,d2,d3,r);    
    scalar radius = oneCrossSec(parameter,angle);
    vector coordXDir = std::cos(angle)*radius*d2;
    vector coordYDir = std::sin(angle)*radius*d3; 
    return r+coordXDir+coordYDir;
}

std::unique_ptr<std::vector<std::vector<LagrangianMarkerOnCrossSec>>> Foam::CrossSectionStructure::constructMarkerSet
(
    const ActiveRodMesh::rodCosserat* oneRod,
    const CrossSection& oneCrossSec
)
{
    LagrangianMarkerOnCrossSec marker1 = createLagrangianMarkerOnCrossSec(oneRod,0,oneCrossSec,0.985247);
    LagrangianMarkerOnCrossSec marker2 = createLagrangianMarkerOnCrossSec(oneRod,0,oneCrossSec,0.985295);
    
    scalar distOnNurbs = distance(marker1,marker2);
    scalar distDirect = Foam::mag(marker1.markerPosition-marker2.markerPosition);
    
    Info<<distOnNurbs<<"|"<<distDirect<<"/ "<<marker1.markerPosition<<"->"<<marker2.markerPosition<<Foam::endl;
    
    //FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);    

    
    
    Info<<"constructMarkerSet"<<Foam::endl;
    std::list<std::list<LagrangianMarkerOnCrossSec>> markers;
    
    // Insert start marker
    scalar startPar = oneRod->m_Curve.domainStart();
    markers.push_back(createInitialCircumMarkers(oneRod,startPar,oneCrossSec));
    constructMarkerSetCircumferential(oneRod,oneCrossSec,markers.back());

    // Insert end marker
    scalar endPar = oneRod->m_Curve.domainEnd();
    markers.push_back(createInitialCircumMarkers(oneRod,endPar,oneCrossSec));
    constructMarkerSetCircumferential(oneRod,oneCrossSec,markers.back());
        
    for(auto mark : markers)
        Info<<mark.front().markerParameter<<Foam::endl;
    
    bool refined=true;
    label refinementCount = 0;
    while(refined)
    {
        refined = false;
        bool cond = true;
        auto markersIter0 = markers.begin();
        auto markersIter1 = ++markers.begin();
        for( ; markersIter1!=markers.end() ; )
        {
            bool subdivide = doSubdivision(*markersIter0, *markersIter1);
            //Info<<"subdivide:"<<subdivide<<Foam::endl;
            if(refinementCount<minRefinement)
                subdivide = true;
            if(subdivide)
            {
                scalar middlePar = 0.5*(markersIter0->front().markerParameter+markersIter1->front().markerParameter);
                //Info<<markersIter0->front().markerParameter<<"->"<<markersIter1->front().markerParameter<<" | middlePar:"<<middlePar<<Foam::endl;
                auto inserted = markers.insert(markersIter1,createInitialCircumMarkers(oneRod,middlePar,oneCrossSec));
                constructMarkerSetCircumferential(oneRod,oneCrossSec,*inserted);
                refined=true;
            }
            markersIter0 = markersIter1;
            markersIter1++;
        }
        refinementCount++;
        Info<<"refinementCount:"<<refinementCount<<Foam::endl;
        if(refinementCount==20)
            FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);    
    }
    
    
    std::unique_ptr<std::vector<std::vector<LagrangianMarkerOnCrossSec>>> markersPtr(new std::vector<std::vector<LagrangianMarkerOnCrossSec>>());
    for(auto iter=markers.begin(); iter!=markers.end(); iter++)
    {
        markersPtr->push_back(std::vector<LagrangianMarkerOnCrossSec>());
        for(auto iterSub=iter->begin(); iterSub!=iter->end(); iterSub++)
            markersPtr->back().push_back(*iterSub);
    }
    return markersPtr;
}

void Foam::CrossSectionStructure::constructMarkerSetCircumferential
(
    const ActiveRodMesh::rodCosserat* oneRod,
    const CrossSection& oneCrossSec,
    std::list<LagrangianMarkerOnCrossSec>& circleMarkers
)
{    
    //Info<<"constructMarkerSetCircumferential"<<Foam::endl;
    
    if(circleMarkers.size()!=2)
        FatalErrorInFunction<<"Must start with 2 markers"<< exit(FatalError);
    scalar parameter = circleMarkers.front().markerParameter;
    
    /*
    vector last;
    scalar pi = constant::mathematical::pi;
    for(scalar rad=-1*pi; rad<3*pi; rad+=0.1)
    {
        vector pos = evaluateRodCircumPos(oneRod,parameter,oneCrossSec,rad);
        Info<<"para:"<<parameter<<"  rad:"<<rad<<"  "<<pos<<"  distToLast:"<<((last-pos)&(last-pos))<<Foam::endl;
        last=pos;
    }
    
    for(auto mark : circleMarkers)
        Info<<mark.markerParameter<<"-"<<mark.markerAngle<<Foam::endl;
    */
    
    bool refined=true;
    label refinementCount = 0;
    while(refined)
    {
        refined = false;
        bool cond = true;
        auto markersIter0 = circleMarkers.begin();
        auto markersIter1 = ++circleMarkers.begin();
        for( ; markersIter1!=circleMarkers.end() ; )
        {
            bool subdivide = doSubdivisionAngle(*markersIter0, *markersIter1);
            if(refinementCount<minRefinement)
                subdivide = true;
            if(subdivide)
            {
                scalar middleAngle = 0.5*(markersIter0->markerAngle+markersIter1->markerAngle);
                //Info<<parameter<<"\t"<<markersIter0->markerAngle<<markersIter0->markerPosition<<"->"<<markersIter1->markerAngle<<markersIter1->markerPosition<<" | middleAngle:"<<middleAngle<<"  ||||||||"<<markersIter0->supportCells.size()<<"|"<<markersIter1->supportCells.size()<<Foam::endl;
                auto inserted = circleMarkers.insert
                    (
                        markersIter1,
                        createLagrangianMarkerOnCrossSec(oneRod,parameter,oneCrossSec,middleAngle)
                    );
                refined=true;
            }
            markersIter0 = markersIter1;
            markersIter1++;
        }
        bool subdivide = doSubdivisionAngle(circleMarkers.back(),circleMarkers.front());
        if(refinementCount<minRefinement)
            subdivide = true;
        if(subdivide)
        {
            scalar middleAngle = 0.5*(circleMarkers.back().markerAngle+2*Foam::constant::mathematical::pi);
            //Info<<parameter<<"\t"<<circleMarkers.back().markerAngle<<circleMarkers.back().markerPosition<<"->"<<2*Foam::constant::mathematical::pi<<circleMarkers.front().markerPosition<<" | middleAngle:"<<middleAngle<<"  ||||||||"<<circleMarkers.back().supportCells.size()<<"|"<<circleMarkers.front().supportCells.size()<<Foam::endl;
            auto inserted = circleMarkers.insert
                (
                    markersIter1,
                    createLagrangianMarkerOnCrossSec(oneRod,parameter,oneCrossSec,middleAngle)
                );
            refined=true;
        }
        //Info<<" angleRefinementCount:"<<refinementCount<<"  refined:"<<refined<<Foam::endl;
        refinementCount++;
        /*
        if(refinementCount==6)
        {
            LagrangianMarkerOnCrossSec marker1 = createLagrangianMarkerOnCrossSec(oneRod,0,oneCrossSec,6.18501);
            LagrangianMarkerOnCrossSec marker2 = createLagrangianMarkerOnCrossSec(oneRod,0,oneCrossSec,6.28319);
    
            scalar distOnNurbs = distance(marker1,marker2);
            scalar distDirect = Foam::mag(marker1.markerPosition-marker2.markerPosition);
    
            Info<<"XXXxxxXXX "<<distOnNurbs<<"|"<<distDirect<<"/ "<<marker1.markerPosition<<"->"<<marker2.markerPosition<<Foam::endl;
            
            
            FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);
        }
        */
    }
}

bool Foam::CrossSectionStructure::doSubdivision
(
    const std::list<LagrangianMarkerOnCrossSec>& smallerSide,
    const std::list<LagrangianMarkerOnCrossSec>& largerSide
)
{
    std::unordered_set<label> smallSideCompleteSupport;
    for(auto smIter=smallerSide.begin(); smIter!=smallerSide.end(); smIter++)
        smallSideCompleteSupport.insert(smIter->supportCells.begin(),smIter->supportCells.end());

    //Info<<"     suppSizesmallerSide:"<<smallSideCompleteSupport.size()<<Foam::endl;    
    
    std::unordered_set<label> largeSideCompleteSupport;
    for(auto laIter=largerSide.begin(); laIter!=largerSide.end(); laIter++)
        largeSideCompleteSupport.insert(laIter->supportCells.begin(),laIter->supportCells.end());
    
    //Info<<"     suppSizelargerSide:"<<largeSideCompleteSupport.size()<<Foam::endl;    
    
    bool allSmallOverlapLarge = true;
    for(auto smIter=smallerSide.begin(); smIter!=smallerSide.end(); smIter++)
    {
        for(auto smSuppCellIter=smIter->supportCells.begin();
            smSuppCellIter!=smIter->supportCells.end();
            smSuppCellIter++)
        {
            if(largeSideCompleteSupport.find(*smSuppCellIter)==largeSideCompleteSupport.end())
                allSmallOverlapLarge = false;
        }
    }
    //Info<<"     allSmallOverlapLarge:"<<allSmallOverlapLarge<<Foam::endl;
    
    bool allLargeOverlapSmall = true;
    for(auto laIter=largerSide.begin(); laIter!=largerSide.end(); laIter++)
    {
        for(auto laSuppCellIter=laIter->supportCells.begin();
            laSuppCellIter!=laIter->supportCells.end();
            laSuppCellIter++)
        {
            if(smallSideCompleteSupport.find(*laSuppCellIter)==smallSideCompleteSupport.end())
                allLargeOverlapSmall = false;
        }
    }
    //Info<<"     allLargeOverlapSmall:"<<allLargeOverlapSmall<<Foam::endl;
        
    if(allSmallOverlapLarge && allLargeOverlapSmall)
    {
        return false;
    }
    return true;
}

bool Foam::CrossSectionStructure::doSubdivisionAngle
(
    const LagrangianMarkerOnCrossSec& smallerSide,
    const LagrangianMarkerOnCrossSec& largerSide
)
{
    bool supportDomainOverlap = false;
    for(auto iterSm=smallerSide.supportCells.begin();
        iterSm!=smallerSide.supportCells.end(); iterSm++)
    {
        for(auto iterLa=largerSide.supportCells.begin();
            iterLa!=largerSide.supportCells.end(); iterLa++)
        {
            
            if(*iterSm == *iterLa)
                supportDomainOverlap = true;
        }
    }
    bool smallerSideHasSupp = smallerSide.supportCells.size()!=0;
    bool largerSideHasSupp = largerSide.supportCells.size()!=0;
    bool bothSidesHaveSupp = largerSideHasSupp&&smallerSideHasSupp;
    
    scalar distOnNurbs = distance(smallerSide,largerSide);
    scalar distDirect = Foam::mag(smallerSide.markerPosition-largerSide.markerPosition);
    
    scalar smallerSideDomainMinSize = std::numeric_limits<scalar>::max();
    if(smallerSideHasSupp)
        smallerSideDomainMinSize = supportDomainMinSize(smallerSide.supportCells);
    
    scalar largerSideDomainMinSize = std::numeric_limits<scalar>::max();
    if(largerSideHasSupp)
        largerSideDomainMinSize = supportDomainMinSize(largerSide.supportCells);
    
    scalar overallDomainMinSide = std::min(smallerSideDomainMinSize,largerSideDomainMinSize);
    
    //Info<<" doSubdivision:"<<distOnNurbs<<"|"<<distDirect<<"/"<<overallDomainMinSide<<smallerSide.markerPosition<<"->"<<largerSide.markerPosition<<" -- "<<distance(smallerSide,largerSide)<<Foam::endl;
    
    if(bothSidesHaveSupp)
    {
        if(distOnNurbs<overallDomainMinSide)
        {
            if(!supportDomainOverlap)
            {
                Info<<Foam::endl;
                std::cout<<"smallerSide.markerParameter:"<<smallerSide.markerParameter<<std::endl;
                std::cout<<"largerSide.markerParameter:"<<largerSide.markerParameter<<std::endl;
                Info<<"distOnNurbs:"<<distOnNurbs<<Foam::endl;
                Info<<"distDirect:"<<distDirect<<Foam::endl;
                Info<<"smallerSideDomainMinSize:"<<smallerSideDomainMinSize<<Foam::endl;
                Info<<"largerSideDomainMinSize:"<<largerSideDomainMinSize<<Foam::endl;
                Info<<"overallDomainMinSide:"<<overallDomainMinSide<<Foam::endl;
                Info<<"smallerSideHasSupp:"<<smallerSideHasSupp<<Foam::endl;
                Info<<"largerSideHasSupp:"<<largerSideHasSupp<<Foam::endl;
                FatalErrorInFunction<<"Dist matches but no overlap"<<exit(FatalError);
            }
            return false;
        }
        else
            return true;
    }
    else if(smallerSideHasSupp || largerSideHasSupp)
    {
        if(distOnNurbs<overallDomainMinSide)
            return false;
        else
            return true;
    }
    else
    {
        if(distOnNurbs<overallDomainMinSide)
            return false;
        else
            return true;        
    }
}

std::list<Foam::LagrangianMarkerOnCrossSec> Foam::CrossSectionStructure::createInitialCircumMarkers
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    const CrossSection& oneCrossSec
)
{
    std::list<LagrangianMarkerOnCrossSec> initialCircle;
    initialCircle.push_back(createLagrangianMarkerOnCrossSec(oneRod,parameter,oneCrossSec,0));
    initialCircle.push_back(createLagrangianMarkerOnCrossSec(oneRod,parameter,oneCrossSec,Foam::constant::mathematical::pi));
    return initialCircle;
}
