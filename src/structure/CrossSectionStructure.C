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

scalar Foam::CrossSectionStructure::restrictAngle
(
    scalar angle
)
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
}

scalar Foam::CrossSectionStructure::distance
(
    const LagrangianMarkerOnCrossSec& A,
    const LagrangianMarkerOnCrossSec& B
)
{
    if(A.baseRod!=B.baseRod)
        FatalErrorInFunction<<"Distance can not be computed between points on different rods"<<exit(FatalError);
    
    //markerParameter,markerAngle,radiusFrac,markerPosition
    using node = std::tuple<scalar,scalar,scalar,vector>;

    node start = {A.markerParameter,A.markerAngle,A.markerRadiusFrac,A.markerPosition};
    std::get<1>(start) = restrictAngle(std::get<1>(start));
    node end = {B.markerParameter,B.markerAngle,B.markerRadiusFrac,B.markerPosition};
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
            std::get<2>(insertNode) = (std::get<2>(*node0)+std::get<2>(*node1))/2;
            
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
            std::get<3>(insertNode) = CrossSectionStructure::evaluateRodCircumPos
                                    (A.baseRod,std::get<0>(insertNode),
                                     A.baseCrossSec,std::get<1>(insertNode),
                                     std::get<2>(insertNode));
            auto inserted = nodes.insert(node1,insertNode);
            node0 = node1;
            node1++;
        }
        
        scalar newDist = 0;
        node0 = nodes.begin();
        node1 = ++nodes.begin();
        for( ; node1!=nodes.end() ; )
        {
            vector connec = std::get<3>(*node1) - std::get<3>(*node0);
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
        std::unique_ptr<std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>>& oneRodMarkers = rodMarkers[rodIndex];
        if(!oneRodMarkers)
        {
            myMesh->m_Rods[rodIndex]->m_Rot.setCoefs(myMesh->m_Rods[rodIndex]->m_Rot_coefs.transpose());
            oneRodMarkers = constructMarkerSet(myMesh->m_Rods[rodIndex],rodCrossSection[rodIndex]);
        }
        
        for(std::vector<std::vector<LagrangianMarkerOnCrossSec>>& radMarker : *oneRodMarkers)
        {
            for(std::vector<LagrangianMarkerOnCrossSec>& circMarker : radMarker)
            {
                for(LagrangianMarkerOnCrossSec& marker : circMarker)
                {
                    connector.markers.push_back(&marker);
                }
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
    scalar angle,
    scalar radiusFrac
)
{
    LagrangianMarkerOnCrossSec marker
    (
        parameter,
        CrossSectionStructure::evaluateRodCircumPos(oneRod,parameter,oneCrossSec,angle,radiusFrac),
        oneRod,
        angle,
        oneCrossSec,
        radiusFrac
    );
    label cellOfMarker = mesh.findCell(marker.markerPosition);
    if(cellOfMarker!=-1)
        getSupportDomain(cellOfMarker,marker.supportCells);
    marker.markerCell = mesh.findCell(marker.markerPosition);
    /*
    for(auto iter=marker.supportCells.begin(); iter!=marker.supportCells.end(); iter++)
        if(mesh.pointInCell(marker.markerPosition,*iter))
        {
            Info<<*iter<<":"<<marker.markerPosition<<Foam::endl;
            if(marker.markerCell!=-1)
                FatalErrorInFunction<<"Double assignment marker cell"<<exit(FatalError);
            marker.markerCell = *iter;
        }
    */
    return marker;
}

Foam::scalar Foam::CrossSectionStructure::evaluateCircumArcLen
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameterA,
    scalar parameterB,
    const CrossSection& oneCrossSec,
    scalar angleA,
    scalar angleB,
    scalar radiusFracA,
    scalar radiusFracB
)
{
    vector parAVec = CrossSectionStructure::evaluateRodCircumPos(oneRod,parameterA,oneCrossSec,angleA,radiusFracA);
    vector parBVec = CrossSectionStructure::evaluateRodCircumPos(oneRod,parameterB,oneCrossSec,angleB,radiusFracB);
    vector connec = parAVec-parBVec;
    return Foam::mag(connec);
}

Foam::vector Foam::CrossSectionStructure::evaluateRodCircumPos
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    const CrossSection& oneCrossSec,
    scalar angle,
    scalar radiusFrac
)
{
    vector d1,d2,d3,r;
    rodEval(oneRod,parameter,d1,d2,d3,r);
    //Info<<Foam::endl<<"parameter:"<<parameter<<Foam::endl;
    scalar radius = oneCrossSec(parameter,angle)*radiusFrac;
    //Info<<"radius:"<<radius<<Foam::endl;
    vector coordXDir = std::cos(angle)*radius*d1;
    //Info<<"coordXDir:"<<coordXDir<<"  angle:"<<angle<<"  "<<d2<<Foam::endl;
    vector coordYDir = std::sin(angle)*radius*d2;
    //Info<<"coordYDir:"<<coordYDir<<"  angle:"<<angle<<"  "<<d3<<Foam::endl;
    return r+coordXDir+coordYDir;
}

std::unique_ptr<std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>> Foam::CrossSectionStructure::constructMarkerSet
(
    const ActiveRodMesh::rodCosserat* oneRod,
    const CrossSection& oneCrossSec
)
{
    /*
    LagrangianMarkerOnCrossSec marker1 = createLagrangianMarkerOnCrossSec(oneRod,0,oneCrossSec,0.985247,1);
    LagrangianMarkerOnCrossSec marker2 = createLagrangianMarkerOnCrossSec(oneRod,0,oneCrossSec,0.985295,1);
    
    scalar distOnNurbs = distance(marker1,marker2);
    scalar distDirect = Foam::mag(marker1.markerPosition-marker2.markerPosition);
    
    Info<<distOnNurbs<<"|"<<distDirect<<"/ "<<marker1.markerPosition<<"->"<<marker2.markerPosition<<Foam::endl;
    */
    
    
    //Info<<"constructMarkerSet"<<Foam::endl;
    std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>> markers;
    markers.resize(2);
    
    // Insert start marker
    scalar startPar = oneRod->m_Curve.domainStart();
    constructMarkerSetRadial(oneRod,startPar,oneCrossSec,markers.front());

    // Insert end marker
    scalar endPar = oneRod->m_Curve.domainEnd();
    constructMarkerSetRadial(oneRod,endPar,oneCrossSec,markers.back());
    
    //Info<<"Start and end done"<<Foam::endl;
    
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
            bool subdivide = doSubdivision(markersIter0->front(), markersIter1->front());
            //Info<<"subdivide:"<<subdivide<<Foam::endl;
            if(refinementCount<minRefinement)
                subdivide = true;
            if(subdivide)
            {
                scalar middlePar = 0.5*(markersIter0->front().front().markerParameter+markersIter1->front().front().markerParameter);
                //Info<<markersIter0->front().markerParameter<<"->"<<markersIter1->front().markerParameter<<" | middlePar:"<<middlePar<<Foam::endl;
                auto inserted = markers.insert
                (
                    markersIter1,
                    {createInitialCircumMarkers(oneRod,middlePar,oneCrossSec,1)}
                );
                constructMarkerSetCircumferential(oneRod,oneCrossSec,inserted->front());
                refined=true;
            }
            markersIter0 = markersIter1;
            markersIter1++;
        }
        refinementCount++;
        //Info<<"refinementCount:"<<refinementCount<<Foam::endl;
        if(refinementCount==20)
            FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);    
    }
    
    std::unique_ptr<std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>> markersPtr(new std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>());
    for(auto iter=markers.begin(); iter!=markers.end(); iter++)
    {
        markersPtr->push_back(std::vector<std::vector<LagrangianMarkerOnCrossSec>>());
        for(auto iterSub=iter->begin(); iterSub!=iter->end(); iterSub++)
        {
            markersPtr->back().push_back(std::vector<LagrangianMarkerOnCrossSec>());
            for(auto iterSubSub=iterSub->begin(); iterSubSub!=iterSub->end(); iterSubSub++)
            {
                markersPtr->back().back().push_back(*iterSubSub);
            }
        }
    }
    
    std::function<vector(scalar,scalar,scalar)> getPosition = 
    [rod=oneRod,crossSec=oneCrossSec](scalar para, scalar radFrac, scalar angle)
    {
        return CrossSectionStructure::evaluateRodCircumPos(rod,para,crossSec,angle,radFrac);
    };
    
    scalar sumVolume = 0;
    std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>& markersRef = *markersPtr;
    for(label rodWiseIndex=0; rodWiseIndex<markers.size(); rodWiseIndex++)
    {
        //Info<<"rodWiseIndex:"<<rodWiseIndex<<Foam::endl;
        std::vector<std::vector<LagrangianMarkerOnCrossSec>>& radialMarkers = markersRef[rodWiseIndex];
        scalar rodWiseIndexParameter = radialMarkers[0][0].markerParameter;
        
        scalar prevParameter;
        if(rodWiseIndex==0)
            prevParameter = rodWiseIndexParameter;
        else
            prevParameter = (markersRef[rodWiseIndex-1][0][0].markerParameter+rodWiseIndexParameter)/2;
        
        scalar subseqParameter;
        if(rodWiseIndex==markers.size()-1)
            subseqParameter = rodWiseIndexParameter;
        else
            subseqParameter = (markersRef[rodWiseIndex+1][0][0].markerParameter+rodWiseIndexParameter)/2;
        
        for(label radialIndex=0; radialIndex<radialMarkers.size(); radialIndex++)
        {
            //Info<<"\tradialIndex:"<<radialIndex<<Foam::endl;
            std::vector<LagrangianMarkerOnCrossSec>& circumMarkers = radialMarkers[radialIndex];
            scalar radialIndexFrac = circumMarkers[0].markerRadiusFrac;
            
            scalar outerRadiusFrac, innerRadiusFrac;
            if(radialMarkers.size()==1)
            {
                outerRadiusFrac = 1;
                innerRadiusFrac = 0;
            }
            else 
            {
                if(radialIndex==0)
                    outerRadiusFrac = radialIndexFrac;
                else
                    outerRadiusFrac = (radialMarkers[radialIndex-1][0].markerRadiusFrac+radialIndexFrac)/2;
                
                if(radialIndex==radialMarkers.size()-1)
                    innerRadiusFrac = radialIndexFrac;
                else
                    innerRadiusFrac = (radialMarkers[radialIndex+1][0].markerRadiusFrac+radialIndexFrac)/2;
            }
            
            for(label circIndex=0; circIndex<circumMarkers.size(); circIndex++)
            {
                //Info<<"\t\tcircIndex:"<<circIndex<<"/"<<circumMarkers.size()<<Foam::endl;

                LagrangianMarkerOnCrossSec& singleMarker = circumMarkers[circIndex];
                scalar circIndexAngle = radialMarkers[radialIndex][0].markerRadiusFrac;

                scalar lowerAngle, upperAngle;
                if(circIndex==0)
                    lowerAngle = circumMarkers.back().markerAngle-2*Foam::constant::mathematical::pi;
                else
                    lowerAngle = circumMarkers[circIndex-1].markerAngle;
                lowerAngle += circIndexAngle;
                lowerAngle /= 2;
                lowerAngle = CrossSectionStructure::restrictAngle(lowerAngle);
                
                if(circIndex==circumMarkers.size()-1)
                    upperAngle = circumMarkers[0].markerAngle+2*Foam::constant::mathematical::pi;
                else
                    upperAngle = circumMarkers[circIndex+1].markerAngle;
                upperAngle += circIndexAngle;
                upperAngle /= 2;
                upperAngle = CrossSectionStructure::restrictAngle(upperAngle);
                
                /*
                Info<<"prevParameter:"<<prevParameter<<"  subseqParameter:"<<subseqParameter<<Foam::endl;
                Info<<"innerRadiusFrac:"<<innerRadiusFrac<<"  outerRadiusFrac:"<<outerRadiusFrac<<Foam::endl;
                Info<<"lowerAngle:"<<lowerAngle<<"  upperAngle:"<<upperAngle<<Foam::endl;
                */
                
                std::array<std::array<std::array<label,2>,2>,2> paraRadAnglePnts;
                pointField points(8);
                points[0] = getPosition(prevParameter,innerRadiusFrac,lowerAngle);
                paraRadAnglePnts[0][0][0] = 0;
                points[1] = getPosition(prevParameter,innerRadiusFrac,upperAngle);
                paraRadAnglePnts[0][0][1] = 1;
                points[2] = getPosition(prevParameter,outerRadiusFrac,lowerAngle);
                paraRadAnglePnts[0][1][0] = 2;
                points[3] = getPosition(prevParameter,outerRadiusFrac,upperAngle);
                paraRadAnglePnts[0][1][1] = 3;
                points[4] = getPosition(subseqParameter,innerRadiusFrac,lowerAngle);
                paraRadAnglePnts[1][0][0] = 4;
                points[5] = getPosition(subseqParameter,innerRadiusFrac,upperAngle);
                paraRadAnglePnts[1][0][1] = 5;
                points[6] = getPosition(subseqParameter,outerRadiusFrac,lowerAngle);
                paraRadAnglePnts[1][1][0] = 6;
                points[7] = getPosition(subseqParameter,outerRadiusFrac,upperAngle);
                paraRadAnglePnts[1][1][1] = 7;
        
                face prevFace({0,1,3,2});
                face subseqFace({4,5,7,6});
                face innerFace({0,1,5,4});
                face outerFace({2,3,7,6});
                face lowerFace({0,2,6,4});
                face upperFace({1,3,7,5});
                
                List<face> faces(12);
                faces[0] = {prevFace[0],prevFace[1],prevFace[3]};
                faces[1] = {prevFace[1],prevFace[2],prevFace[3]};
                
                faces[2] = {subseqFace[0],subseqFace[1],subseqFace[3]};
                faces[3] = {subseqFace[1],subseqFace[2],subseqFace[3]};
                
                faces[4] = {innerFace[0],innerFace[1],innerFace[3]};
                faces[5] = {innerFace[1],innerFace[2],innerFace[3]};
                
                faces[6] = {outerFace[0],outerFace[1],outerFace[3]};
                faces[7] = {outerFace[1],outerFace[2],outerFace[3]};
                
                faces[8] = {lowerFace[0],lowerFace[1],lowerFace[3]};
                faces[9] = {lowerFace[1],lowerFace[2],lowerFace[3]};
                
                faces[10] = {upperFace[0],upperFace[1],upperFace[3]};
                faces[11] = {upperFace[1],upperFace[2],upperFace[3]};
                
                List<label> cellList(12);
                for(label ind=0; ind<cellList.size(); ind++)
                    cellList[ind] = ind;
                cell thisCell(cellList);
                
                singleMarker.markerVolume = thisCell.mag(points,faces);
                sumVolume+=singleMarker.markerVolume;
            }
        }
    }
    return markersPtr;
}

void Foam::CrossSectionStructure::constructMarkerSetRadial
(
    const ActiveRodMesh::rodCosserat* oneRod,
    const scalar parameter,
    const CrossSection& oneCrossSec,
    std::list<std::list<LagrangianMarkerOnCrossSec>>& radialMarkers
)
{
    if(!radialMarkers.empty())
        FatalErrorInFunction<<"List must be empty"<< exit(FatalError);
    
    //Info<<"constructMarkerSetRadial:"<<parameter<<Foam::endl;
    
    // Insert outer circum markers
    radialMarkers.push_back(createInitialCircumMarkers(oneRod,parameter,oneCrossSec,1));
    constructMarkerSetCircumferential(oneRod,oneCrossSec,radialMarkers.back());
    
    //Info<<"Done outer cirum markers"<<Foam::endl;

    // Insert inner circum marker
    LagrangianMarkerOnCrossSec marker
    (
        parameter,
        evaluateRodPos(oneRod,parameter),
        oneRod,
        0,
        oneCrossSec,
        0
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
    radialMarkers.push_back({marker});
    
    //Info<<"Done inner cirum markers"<<Foam::endl;

    
    bool refined=true;
    label refinementCount = 0;
    while(refined)
    {
        refined = false;
        bool cond = true;
        auto markersIter0 = radialMarkers.begin();
        auto markersIter1 = ++radialMarkers.begin();
        for( ; markersIter1!=radialMarkers.end() ; )
        {
            bool subdivide = doSubdivision(*markersIter0, *markersIter1);
            //Info<<"\t\tsubdivide:"<<subdivide<<Foam::endl;
            if(refinementCount<minRefinement)
                subdivide = true;
            if(subdivide)
            {
                scalar middleRadiusFrac = 0.5*(markersIter0->front().markerRadiusFrac+markersIter1->front().markerRadiusFrac);
                //Info<<"\t\t"<<markersIter0->front().markerRadiusFrac<<"->"<<markersIter1->front().markerRadiusFrac<<" | middleRadiusFrac:"<<middleRadiusFrac<<Foam::endl;
                auto inserted = radialMarkers.insert
                (
                    markersIter1,
                    createInitialCircumMarkers(oneRod,parameter,oneCrossSec,middleRadiusFrac)
                );
                constructMarkerSetCircumferential(oneRod,oneCrossSec,*inserted,middleRadiusFrac);
                refined=true;
            }
            markersIter0 = markersIter1;
            markersIter1++;
        }
        refinementCount++;
        //Info<<"\tconstructMarkerSetRadial refinementCount:"<<refinementCount<<Foam::endl;
        if(refinementCount==7)
            FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);    
    }
}

void Foam::CrossSectionStructure::constructMarkerSetCircumferential
(
    const ActiveRodMesh::rodCosserat* oneRod,
    const CrossSection& oneCrossSec,
    std::list<LagrangianMarkerOnCrossSec>& circleMarkers,
    scalar radiusFrac
)
{    
    //Info<<"constructMarkerSetCircumferential:"<<radiusFrac<<Foam::endl;
    
    if(circleMarkers.size()!=2)
        FatalErrorInFunction<<"Must start with 2 markers"<< exit(FatalError);
    scalar parameter = circleMarkers.front().markerParameter;
    
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
            bool subdivide = doSubdivisionCircumferential(*markersIter0, *markersIter1);
            if(refinementCount<minRefinement)
                subdivide = true;
            if(subdivide)
            {
                scalar middleAngle = 0.5*(markersIter0->markerAngle+markersIter1->markerAngle);
                //Info<<parameter<<"\t"<<markersIter0->markerAngle<<markersIter0->markerPosition<<"->"<<markersIter1->markerAngle<<markersIter1->markerPosition<<" | middleAngle:"<<middleAngle<<"  ||||||||"<<markersIter0->supportCells.size()<<"|"<<markersIter1->supportCells.size()<<Foam::endl;
                auto inserted = circleMarkers.insert
                    (
                        markersIter1,
                        createLagrangianMarkerOnCrossSec(oneRod,parameter,oneCrossSec,middleAngle,radiusFrac)
                    );
                refined=true;
            }
            markersIter0 = markersIter1;
            markersIter1++;
        }
        bool subdivide = doSubdivisionCircumferential(circleMarkers.back(),circleMarkers.front());
        if(refinementCount<minRefinement)
            subdivide = true;
        if(subdivide)
        {
            scalar middleAngle = 0.5*(circleMarkers.back().markerAngle+2*Foam::constant::mathematical::pi);
            //Info<<parameter<<"\t"<<circleMarkers.back().markerAngle<<circleMarkers.back().markerPosition<<"->"<<2*Foam::constant::mathematical::pi<<circleMarkers.front().markerPosition<<" | middleAngle:"<<middleAngle<<"  ||||||||"<<circleMarkers.back().supportCells.size()<<"|"<<circleMarkers.front().supportCells.size()<<Foam::endl;
            auto inserted = circleMarkers.insert
                (
                    markersIter1,
                    createLagrangianMarkerOnCrossSec(oneRod,parameter,oneCrossSec,middleAngle,radiusFrac)
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
    std::unordered_set<label> smallSideCells;
    for(auto smIter=smallerSide.begin(); smIter!=smallerSide.end(); smIter++)
        smallSideCells.insert(smIter->markerCell);
    smallSideCells.erase(-1);
    
    std::unordered_set<label> largeSideCompleteSupport;
    for(auto laIter=largerSide.begin(); laIter!=largerSide.end(); laIter++)
        largeSideCompleteSupport.insert(laIter->supportCells.begin(),laIter->supportCells.end());
    std::unordered_set<label> largeSideCells;
    for(auto laIter=largerSide.begin(); laIter!=largerSide.end(); laIter++)
        largeSideCells.insert(laIter->markerCell);
    largeSideCells.erase(-1);

    bool allSmallOverlapLarge = true;
    for(auto laIter=largeSideCells.begin(); laIter!=largeSideCells.end(); laIter++)
    {
        label laCell = *laIter;
        if(smallSideCompleteSupport.find(laCell)==smallSideCompleteSupport.end())
            allSmallOverlapLarge = false;
    }
    //Info<<"\t\t\t"<<smallerSide.size()<<"     allSmallOverlapLarge:"<<allSmallOverlapLarge<<"  smallSideCompleteSupportSize:"<<smallSideCompleteSupport.size()<<"  smallSideCellsSize:"<<smallSideCells.size()<<"  (p:"<<smallerSide.front().markerParameter<<",r:"<<smallerSide.front().markerRadiusFrac<<")"<<Foam::endl;
    
    bool allLargeOverlapSmall = true;
    for(auto smIter=smallSideCells.begin(); smIter!=smallSideCells.end(); smIter++)
    {
        label smCell = *smIter;
        if(largeSideCompleteSupport.find(smCell)==largeSideCompleteSupport.end())
            allLargeOverlapSmall = false;
    }
    //Info<<"\t\t\t"<<largerSide.size()<<"     allLargeOverlapSmall:"<<allLargeOverlapSmall<<"  largeSideCompleteSupportSize:"<<largeSideCompleteSupport.size()<<"  largeSideCellsSize:"<<largeSideCells.size()<<"  (p:"<<largerSide.front().markerParameter<<",r:"<<largerSide.front().markerRadiusFrac<<")"<<Foam::endl;
        
    if(allSmallOverlapLarge && allLargeOverlapSmall)
        return false;
    if(allSmallOverlapLarge && largeSideCompleteSupport.size()==0)
        return false;
    if(smallSideCompleteSupport.size()==0 && allLargeOverlapSmall)
        return false;
    if(smallSideCompleteSupport.size()==0 && largeSideCompleteSupport.size()==0)
        return false;
    
    return true;
}

bool Foam::CrossSectionStructure::doSubdivisionCircumferential
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
    const CrossSection& oneCrossSec,
    scalar radiusFrac
)
{
    //Info<<" createInitialCircumMarkers:"<<parameter<<" - "<<radiusFrac<<Foam::endl;
    
    std::list<LagrangianMarkerOnCrossSec> initialCircle;
    initialCircle.push_back(createLagrangianMarkerOnCrossSec(oneRod,parameter,oneCrossSec,0,radiusFrac));
    initialCircle.push_back(createLagrangianMarkerOnCrossSec(oneRod,parameter,oneCrossSec,Foam::constant::mathematical::pi,radiusFrac));
    return initialCircle;
}
