#include "LineStructure.H"

std::string Foam::LagrangianMarker::to_string()
{
    return std::to_string(markerParameter);
}

vector Foam::LagrangianMarker::getMarkerVelocity()
{
    return vector(0,0,0);
}

scalar Foam::LagrangianMarker::getMarkerTemperature()
{
    return 1;
}

scalar Foam::LineStructure::distance
(
    const LagrangianMarker& A,
    const LagrangianMarker& B
)
{
    if(A.baseRod!=B.baseRod)
        FatalErrorInFunction<<"Distance can not be computed between points on different rods"<<exit(FatalError);
        
    using node = std::pair<scalar,vector>;
    
    node start = {A.markerParameter,A.markerPosition};
    node end = {B.markerParameter,B.markerPosition};
    
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
            insertNode.first = (node0->first+node1->first)/2;
            insertNode.second = evaluateRodPos(A.baseRod,insertNode.first);
            auto inserted = nodes.insert(node1,insertNode);
            node0 = node1;
            node1++;
        }
        
        scalar newDist = 0;
        node0 = nodes.begin();
        node1 = ++nodes.begin();
        for( ; node1!=nodes.end() ; )
        {
            vector connec = node1->second - node0->second;
            newDist += std::sqrt(connec&connec);
            node0++;
            node1++;
        }
        
        scalar change = std::abs(newDist-distance);
        scalar changePerc = change/newDist;
        distance = newDist;
        if(nodes.size()>10 && changePerc<1e-3)
            errorSufficient=true;
        
        //Info<<" distance:"<<distance<<Foam::endl;
    }
    //Info<<"     distance("<<A.markerPosition<<","<<B.markerPosition<<"):"<<distance<<Foam::endl;
    
    return distance;
}

Foam::FieldMarkerStructureInteraction::FieldMarkerStructureInteraction
(
    dynamicRefineFvMesh& mesh,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
mesh(mesh),
h(std::cbrt(mesh.cells()[0].mag(mesh.points(),mesh.faces()))),
modusFieldToMarker(modusFieldToMarker),
modusMarkerToField(modusMarkerToField)
{
}

Foam::scalar Foam::FieldMarkerStructureInteraction::phiFunction(Foam::scalar r)
{
    //Info<<"phiFunction -- r:"<<r<<Foam::endl;
    
    scalar abs_r = std::abs(r);
    if(abs_r < 0.5)
    {
        scalar result = (1.0/3.0)*(1+std::sqrt(-3*(abs_r*abs_r)+1));
        //Info<<"abs_r:"<<abs_r<<" -> "<<result<<Foam::endl;
        return result;
    }
    else if(abs_r <= 1.5)
    {
        scalar result = (1.0/6.0)*(5.0-3.0*abs_r-std::sqrt(-3.0*((1-abs_r)*(1-abs_r))+1));
        //Info<<"abs_r:"<<abs_r<<" -> "<<result<<Foam::endl;
        return result;        
    }
    else
    {
        //Info<<"abs_r:"<<abs_r<<" -> "<<0<<Foam::endl;
        return 0;
    }
}

Foam::scalar Foam::FieldMarkerStructureInteraction::deltaDirac
(
    Foam::vector X,
    Foam::vector x,
    Foam::scalar h
)
{   
    vector sigma_d;
    for(label dim=0; dim<3; dim++)
    {
        scalar X_i_x_i = X[dim]-x[dim];
        scalar r = X_i_x_i / h;
        sigma_d[dim] = phiFunction(r);
    }
    scalar deltaDir = sigma_d[0]*sigma_d[1]*sigma_d[2];
    deltaDir /= (h*h*h);
    //Info<<h<<"  X:"<<X<<"  x:"<<x<<" -> "<<deltaDir<<Foam::endl;
    return deltaDir;
}

void Foam::FieldMarkerStructureInteraction::scatterNurbs
(
    std::pair<gsNurbs<scalar>,label> in,
    std::pair<gsNurbs<scalar>,label>& out
)
{
    List<List<scalar>> nurbsData(7);
    // degree
    // gsKnotVector
    // weights
    // cols // nbrControlPoints
    // rows // dimensionControlPoints
    // coefficients
    // number
    if(Pstream::master())
    {
        gismo::gsKnotVector<scalar> knots = in.first.knots();
        label degree = knots.degree();
        nurbsData[0] = List<scalar>(1);
        nurbsData[0][0] = degree;
        nurbsData[1] = List<scalar>(knots.size());
        for(label i=0;i<knots.size();i++)
            nurbsData[1][i] = knots[i];
        
        gismo::gsMatrix<scalar> weights = in.first.weights();
        if(weights.rows()!=1)
            FatalErrorInFunction<<"Wrong weight dimensions"<< exit(FatalError);
        nurbsData[2] = List<scalar>(weights.cols());
        for(label i=0;i<weights.cols();i++)
            nurbsData[2][i] = weights(0,i);

        gismo::gsMatrix<scalar> coefs = in.first.coefs();
        nurbsData[3] = List<scalar>(1);
        nurbsData[3][0] = coefs.cols();
        nurbsData[4] = List<scalar>(1);
        nurbsData[4][0] = coefs.rows();
        nurbsData[5] = List<scalar>(coefs.cols()*coefs.rows());
        label index = 0;
        for(int n=0;n<coefs.cols();n++)
        {
            for(int d=0;d<coefs.rows();d++)
            {
                nurbsData[5][index] = coefs(d,n);
                index++;
            }
        }
        
        nurbsData[6] = List<scalar>(1);
        nurbsData[6][0] = in.second;
    }
    Pstream::scatter(nurbsData);
    if(Pstream::master())
    {
        out = in;
    }
    else
    {
        std::vector<double> knotContainer(nurbsData[1].size());
        for(label i=0;i<nurbsData[1].size();i++)
            knotContainer[i] = nurbsData[1][i];
        gismo::gsKnotVector<scalar> knots(knotContainer,nurbsData[0][0]);

        gismo::gsMatrix<scalar> weights(nurbsData[2].size(),1);
        for(label i=0;i<nurbsData[2].size();i++)
            weights(i,0) = nurbsData[2][i];

        gismo::gsMatrix<scalar> coefs(nurbsData[4][0],nurbsData[3][0]);
        nurbsData[5] = List<scalar>(coefs.cols()*coefs.rows());
        label index = 0;
        for(int n=0;n<coefs.cols();n++)
        {
            for(int d=0;d<coefs.rows();d++)
            {
                coefs(d,n) = nurbsData[5][index];;
                index++;
            }
        }
       
        out.first = gismo::gsNurbs<double>(knots,weights,coefs);
        
        out.second = nurbsData[6][0];
    }
}

Foam::LineStructure::LineStructure
(
    dynamicRefineFvMesh& mesh,
    /*
    const dimensionedScalar& alpha,
    const volScalarField& T,
    const volScalarField& p,
    const volVectorField& U,
    const dimensionedScalar nu,
    */
    const List<scalar> crossSecArea
):
Structure(mesh,mesh.time()/*,alpha,T,p,U,nu*/),
crossSecArea(crossSecArea)
{
    
}

void Foam::LineStructure::transferMarkers(FieldMarkerStructureInteraction& connector)
{
    if(!myMesh)
        FatalErrorInFunction<<"Rod Mesh not set!"<<exit(FatalError);
    if(myMesh->m_nR!=rodMarkers.size())
    {
        rodMarkers.resize(myMesh->m_nR);
    }
    connector.markers.resize(0);
    
    if(crossSecArea.size()!=rodMarkers.size())
        FatalErrorInFunction<<"Mismatch in size of crossSecArea and rodMarkers"<<exit(FatalError);
    
    for(uint rodIndex=0; rodIndex<rodMarkers.size(); rodIndex++)
    {
        std::unique_ptr<std::vector<LagrangianMarker>>& oneRodMarkers = rodMarkers[rodIndex];
        if(!oneRodMarkers)
        {
            myMesh->m_Rods[rodIndex]->m_Rot.setCoefs(myMesh->m_Rods[rodIndex]->m_Rot_coefs.transpose());
            oneRodMarkers = constructMarkerSet(rodIndex,myMesh->m_Rods[rodIndex],crossSecArea[rodIndex]);
        }
        for(LagrangianMarker& marker : *oneRodMarkers)
        {
            connector.markers.push_back(&marker);
        }
    }
}

void Foam::LineStructure::getSupportDomain
(
    label cellInd,
    std::unordered_set<label>& neighbourhoodCells,
    label iterations
)
{
    const cellList& cellList = mesh.cells();
    const faceList& facesList = mesh.faces();
    if(cellInd<0 || cellInd>=cellList.size())
        FatalErrorInFunction<<"Invalid cell index"<< exit(FatalError);
    neighbourhoodCells.insert(cellInd);
    
    for(label iter=0; iter<iterations; iter++)
    {
        DynamicList<label> newCells;
        for(auto cellIter=neighbourhoodCells.begin(); cellIter!=neighbourhoodCells.end(); cellIter++)
        {
            const cell& thisCell = cellList[*cellIter];
            labelList cellVertices = thisCell.labels(facesList);
            for(label vertice : cellVertices)
            {
                labelList verticeCells = mesh.pointCells(vertice);
                newCells.append(verticeCells);
            }
        }
        neighbourhoodCells.insert(newCells.begin(),newCells.end());
    }    
}

std::pair<double,double> Foam::LineStructure::minMaxSpan
(
    const cell& thisCell
)
{
    std::vector<double> spans;
    cellDistances(thisCell,spans);
    std::sort(spans.begin(),spans.end());
    if(spans.size()<24)
        FatalErrorInFunction<<"Cell can not have less than 24 spans"<< exit(FatalError);
    return {spans.front(),spans.back()};
}

void Foam::LineStructure::cellDistances
(
    const cell& thisCell,
    std::vector<double>& spans
)
{
    const pointField& cellPoints = thisCell.points(mesh.faces(),mesh.points());
    for(label pntI=0; pntI<cellPoints.size(); pntI++)
    {
        for(label pntJ=0; pntJ<cellPoints.size(); pntJ++)
        {
            if(pntI!=pntJ)
            {
                vector distanceVec = cellPoints[pntI]-cellPoints[pntJ];
                scalar distanceVal = Foam::mag(distanceVec);
                spans.push_back(distanceVal);
            }
        }
    }    
}

scalar Foam::LineStructure::supportDomainMinSize
(
    label cellInd
)
{
    //Info<<"     supportDomainMinSize:"<<cellInd<<Foam::endl;
    const cellList& cells = mesh.cells();
    
    double minSuppSize = std::numeric_limits<double>::max(); 
    std::unordered_set<label> neighbourhoodCells;
    getSupportDomain(cellInd,neighbourhoodCells);
    //Info<<"     neighbourhoodCells.size():"<<neighbourhoodCells.size()<<Foam::endl;
    for(auto cellIter=neighbourhoodCells.begin(); cellIter!=neighbourhoodCells.end(); cellIter++)
    {
        //Info<<"     cell:"<<*cellIter<<Foam::endl;
        const cell& suppCell = cells[*cellIter];
        double minSpan = minMaxSpan(suppCell).first;
        minSuppSize = std::min<double>(minSuppSize,minSpan);
    }
    //Info<<"     supportDomainMinSize done"<<Foam::endl;
    return minSuppSize;
}

scalar Foam::LineStructure::supportDomainMinSize
(
    const std::unordered_set<label>& supportDomainCells
)
{    
    const cellList& cells = mesh.cells();    
    double minSuppSize = std::numeric_limits<double>::max(); 
    for(auto cellIter=supportDomainCells.begin(); cellIter!=supportDomainCells.end(); cellIter++)
    {
        const cell& suppCell = cells[*cellIter];
        double minSpan = minMaxSpan(suppCell).first;
        minSuppSize = std::min<double>(minSuppSize,minSpan);
    }
    return minSuppSize;
}

std::unique_ptr<std::vector<LagrangianMarker>> Foam::LineStructure::constructMarkerSet
(
    label rodNumber,
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar crossSecArea
)
{
    Info<<"constructMarkerSet"<<Foam::endl;
    
    std::list<LagrangianMarker> markers;
    
    // Insert start marker
    scalar startPar = oneRod->m_Curve.domainStart();
    markers.push_back(createLagrangianMarker(startPar,rodNumber,oneRod));
    Info<<"Start marker created: "<<startPar<<Foam::endl;
    Info<<"markers.back():"<<markers.back().markerParameter<<" - "<<markers.back().markerPosition<<Foam::endl;
    
    // Insert end marker
    scalar endPar = oneRod->m_Curve.domainEnd();
    markers.push_back(createLagrangianMarker(endPar,rodNumber,oneRod));
    Info<<"End marker created: "<<endPar<<Foam::endl;
    Info<<"markers.back():"<<markers.back().markerParameter<<" - "<<markers.back().markerPosition<<Foam::endl;
    
    bool refined=true;
    label refinementCount = 0;
    while(refined)
    {
        //Info<<"Refine "<<refinementCount<<" size:"<<markers.size()<<Foam::endl;
        refined = false;
        bool cond = true;
        auto markersIter0 = markers.begin();
        auto markersIter1 = ++(markers.begin());
        for( ; markersIter1!=markers.end() ; )
        {
            bool subdivide = doSubdivision(*markersIter0, *markersIter1);
            if(refinementCount<minRefinement)
                subdivide = true;
            if(subdivide)
            {
                scalar middlePar = 0.5*(markersIter0->markerParameter+markersIter1->markerParameter);
                auto inserted = markers.insert(markersIter1,createLagrangianMarker(middlePar,rodNumber,oneRod));
                refined=true;
                //Info<<"   Subdivision:"<<markersIter0->markerParameter<<"--"<<markersIter1->markerParameter<<" / "<<middlePar<<Foam::endl;
            }
            else
            {
                //Info<<"   Non Subdivision:"<<markersIter0->markerParameter<<"--"<<markersIter1->markerParameter<<Foam::endl;
            }
            markersIter0 = markersIter1;
            markersIter1++;
        }
        refinementCount++;
    }
        
    std::list<LagrangianMarker>::iterator iterPrev = markers.end();
    std::list<LagrangianMarker>::iterator iterNext;
    for(auto iter=markers.begin(); iter!=markers.end(); iter++)
    {
        iterNext = iter;
        iterNext++;
        
        scalar spanStart;
        if(iterPrev!=markers.end())
        {
            spanStart = iterPrev->markerParameter;
            spanStart = spanStart + iter->markerParameter;
            spanStart /= 2;
        }
        else
            spanStart = iter->markerParameter;

        scalar spanEnd;
        if(iterNext!=markers.end())
        {            
            spanEnd = iterNext->markerParameter;
            spanEnd = spanEnd + iter->markerParameter;
            spanEnd /= 2;
        }
        else
            spanEnd = iter->markerParameter;
        
        if(iterPrev==markers.end() && iterNext==markers.end())
            FatalErrorInFunction<<"Marker with no predecessor and no succesor"<<exit(FatalError);
        
        scalar span = spanEnd-spanStart;
        iter->markerVolume = span*crossSecArea;
        
        iterPrev = iter;
    }

    std::unique_ptr<std::vector<LagrangianMarker>> markersPtr(new std::vector<LagrangianMarker>());
    for(auto iter=markers.begin(); iter!=markers.end(); iter++)
    {
        Info<<iter->markerParameter<<" -- "<<iter->markerPosition<<" "<<iter->markerCell<<Foam::endl;
        if(iter->markerCell!=-1)
            markersPtr->push_back(*iter);
    }
    
    Info<<"Created "<<markersPtr->size()<<" markers!"<<Foam::endl;
    return markersPtr;
}

Foam::vector Foam::LineStructure::evaluateRodPos
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter
)
{
    vector r;
    rodEval(oneRod,parameter,r);
    return r;
}

/*
Foam::vector Foam::LineStructure::evaluateRodDeriv
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter
)
{
    gsMatrix<double,1,1> parMat(1,1);
    parMat.at(0) = parameter;
    
    gsMatrix<double> basePntVec;
    gsNurbs<double> base = oneRod->m_Curve;
    base.deriv_into(parMat,basePntVec);

    gsMatrix<double> defPntVec;
    gsNurbs<double> def = oneRod->m_Def;    
    def.deriv_into(parMat,defPntVec);
    
    gsMatrix<double> resultMat = basePntVec+defPntVec;
    return vector(resultMat(0,0),resultMat(0,1),resultMat(0,2));    
}
*/

Foam::LagrangianMarker Foam::LineStructure::createLagrangianMarker
(
    scalar markerParameter,
    label rodNumber,
    const ActiveRodMesh::rodCosserat* oneRod
)
{
    //Info<<" createLagrangianMarker"<<Foam::endl;
    LagrangianMarker marker(markerParameter,evaluateRodPos(oneRod,markerParameter),rodNumber,oneRod);
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

Foam::scalar Foam::LineStructure::evaluateRodArcLen
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parA,
    scalar parB
)
{
    vector parAVec = evaluateRodPos(oneRod,parA);
    vector parBVec = evaluateRodPos(oneRod,parB);
    vector connec = parAVec-parBVec;
    return Foam::mag(connec);
}

bool Foam::LineStructure::doSubdivision
(
    const LagrangianMarker& smallerSide,
    const LagrangianMarker& largerSide
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
