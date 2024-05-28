#include "LineStructure.H"

Foam::LineStructure::LineStructure
(
    dynamicRefineFvMesh& mesh,
    const List<scalar> crossSecArea,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
Structure(mesh,mesh.time()),
modusFieldToMarker(modusFieldToMarker),
modusMarkerToField(modusMarkerToField),
crossSecArea(crossSecArea),
initialSpacing(initialSpacingFromMesh())
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
        //Info<<"rodIndex:"<<rodIndex<<" initialSpacing:"<<initialSpacing<<Foam::endl;
        std::unique_ptr<std::vector<LagrangianMarker>>& oneRodMarkers = rodMarkers[rodIndex];
        if(!oneRodMarkers)
        {
            myMesh->m_Rods[rodIndex]->m_Rot.setCoefs(myMesh->m_Rods[rodIndex]->m_Rot_coefs.transpose());
            oneRodMarkers = constructMarkerSet(rodIndex,myMesh->m_Rods[rodIndex],crossSecArea[rodIndex],initialSpacing);
        }
        for(LagrangianMarker& marker : *oneRodMarkers)
        {
            connector.markers.push_back(&marker);
        }
    }
}

std::unique_ptr<std::vector<scalar>> Foam::LineStructure::createSpacedPointsOnRod
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar spacing
)
{
    auto pointsPtr = std::unique_ptr<std::vector<scalar>>(new std::vector<scalar>());
    std::vector<scalar>& pointsVec = *pointsPtr;
    
    std::list<scalar> points;
    points.push_back(oneRod->m_Curve.domainStart());
    points.push_back(oneRod->m_Curve.domainEnd());
    
    bool refined=true;
    while(refined)
    {
        refined = false;
        bool cond = true;
        scalar summedDist = 0;
        auto pntsIter0 = points.begin();
        auto pntsIter1 = ++(points.begin());
        for( ; pntsIter1!=points.end() ; )
        {
            scalar dist = distance(oneRod,*pntsIter0,*pntsIter1);
            summedDist+=dist;
            if(dist>spacing)
            {
                scalar middlePar = 0.5*(*pntsIter0 + *pntsIter1);
                auto inserted = points.insert(pntsIter1,middlePar);
                refined=true;
            }
            pntsIter0 = pntsIter1;
            pntsIter1++;
        }
    }
    
    for(scalar para : points)
    {
        pointsVec.push_back(para);
    }
    
    return pointsPtr;
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

/*
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
*/

scalar Foam::LineStructure::supportDomainMinSize
(
    const DynamicList<label>& supportDomainCells
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

std::unique_ptr<std::vector<LagrangianMarker>> Foam::LineStructure::constructMarkerSet
(
    label rodNumber,
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar crossSecArea,
    scalar initialSpacing,
    std::pair<bool,scalar> refineSpacing,
    bool reInitialize
)
{
    //Create initial spacing
    std::unique_ptr<std::vector<scalar>>& spacedPoints = initialRodPoints[rodNumber];
    if(reInitialize || !spacedPoints)
    {
        spacedPoints = createSpacedPointsOnRod(oneRod,initialSpacing);
        if(spacedPoints->size()<2)
            FatalErrorInFunction<<"Initial spaced points must be at least two"<< exit(FatalError);
    }

    //Create initial markers
    std::list<LagrangianMarker> markers;
    for(scalar point : *spacedPoints)
        markers.push_back(LagrangianMarker(mesh,rodNumber,oneRod,point));

    //Marker refinement
    bool refined=true;
    while(refined)
    {
        refined = false;
        bool cond = true;
        auto markersIter0 = markers.begin();
        auto markersIter1 = ++(markers.begin());
        for( ; markersIter1!=markers.end() ; )
        {
            scalar markers0Para = markersIter0->getMarkerParameter();
            scalar markers0CellSpacing = markersIter0->getMarkerCellMinSpacing();
            bool markers0InCell = (markersIter0->getMarkerCell()!=-1);
            
            scalar markers1Para = markersIter0->getMarkerParameter();
            scalar markers1CellSpacing = markersIter0->getMarkerCellMinSpacing();
            bool markers1InCell = (markersIter1->getMarkerCell()!=-1);
            
            scalar dist = distance(oneRod,markers0Para,markers1Para);
            
            bool subdivide = false;
            if(refineSpacing.first)
            {
                if(dist>refineSpacing.second)
                    subdivide=true;
            }
            if(markers0InCell || markers1InCell)
            {
                scalar minSpacing = std::min(markers0CellSpacing,markers1CellSpacing);
                if(dist>minSpacing)
                    subdivide=true;
            }
        
            if(subdivide)
            {
                scalar middlePar = 0.5*(markers0Para+markers1Para);
                LagrangianMarker middleMarker(mesh,rodNumber,oneRod,middlePar);
                auto inserted = markers.insert(markersIter1,middleMarker);
                refined=true;
            }
            markersIter0 = markersIter1;
            markersIter1++;
        }
    }

    // Compute marker volume
    std::list<LagrangianMarker>::iterator iterPrev = markers.end();
    std::list<LagrangianMarker>::iterator iterNext;
    for(auto iter=markers.begin(); iter!=markers.end(); iter++)
    {
        iterNext = iter;
        iterNext++;
        
        scalar spanStart;
        if(iterPrev!=markers.end())
        {
            spanStart = iterPrev->getMarkerParameter();
            spanStart = spanStart + iter->getMarkerParameter();
            spanStart /= 2;
        }
        else
            spanStart = iter->getMarkerParameter();

        scalar spanEnd;
        if(iterNext!=markers.end())
        {            
            spanEnd = iterNext->getMarkerParameter();
            spanEnd = spanEnd + iter->getMarkerParameter();
            spanEnd /= 2;
        }
        else
            spanEnd = iter->getMarkerParameter();
        
        if(iterPrev==markers.end() && iterNext==markers.end())
            FatalErrorInFunction<<"Marker with no predecessor and no succesor"<<exit(FatalError);
        
        scalar span = spanEnd-spanStart;
        if(modusMarkerToField==markerMeshType::NonUniform)
            iter->setMarkerVolume(span);
        else
            iter->setMarkerVolume(span*crossSecArea);
        
        iterPrev = iter;
    }

    //Transfer to std::vector
    std::unique_ptr<std::vector<LagrangianMarker>> markersPtr(new std::vector<LagrangianMarker>());
    for(auto iter=markers.begin(); iter!=markers.end(); iter++)
    {
        Info<<iter->getMarkerParameter()<<" -- "<<iter->getMarkerPosition()<<" "<<iter->getMarkerCell()<<Foam::endl;
        if(iter->getMarkerCell()!=-1)
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

scalar Foam::LineStructure::distance
(
    const ActiveRodMesh::rodCosserat* oneRod,   
    scalar parA,
    scalar parB
)
{
    std::function<scalar(scalar)> curveLen = 
    [rodPtr = oneRod](scalar parameter)
    {
        vector derivRod = rodDerivEval(rodPtr,parameter);
        return std::sqrt(derivRod&derivRod);
    };    
    return integrateRodwise<scalar>(oneRod,parA,parB,curveLen);
}

scalar Foam::LineStructure::distance
(
    const LagrangianMarker& A,
    const LagrangianMarker& B
)
{
    if(A.getBaseRod()!=B.getBaseRod())
        FatalErrorInFunction<<"Distance can not be computed between points on different rods"<<exit(FatalError);
    
    return distance(A.getBaseRod(),A.getMarkerParameter(),B.getMarkerParameter());
}

void Foam::LineStructure::extractStructureHaloMarkers
(
    const std::vector<LagrangianMarker*> markers 
)
{
    const DynamicList<label>& haloCells = getSelfHaloCells();
    haloMarkers.resize(haloCells.size());
    const std::unordered_map<label,DynamicList<label>>& haloCellsMap = getHaloCellMap();
    const std::unordered_map<label,label> selfHaloListIndexMap = getSelfHaloListIndexMap();
    for(const LagrangianMarker* markerPtr : markers)
    {
        label cell = markerPtr->getMarkerCell();
        if(haloCellsMap.find(cell)!=haloCellsMap.end())
        {
            auto iter = selfHaloListIndexMap.find(cell);
            if(iter==selfHaloListIndexMap.end())
                FatalErrorInFunction<<"Can not happen"<<exit(FatalError);
            label index = iter->second;
            haloMarkers[index].push_back(markerPtr);
        }
    }
}


/*
bool Foam::LineStructure::doSubdivision
(
    const LagrangianMarker& smallerSide,
    const LagrangianMarker& largerSide
)
{
    FatalErrorInFunction<<"Not implemented"<<exit(FatalError);
    
    bool supportDomainOverlap = false;
    for(auto iterSm=smallerSide.getSupportCells().begin();
        iterSm!=smallerSide.getSupportCells().end(); iterSm++)
    {
        for(auto iterLa=largerSide.getSupportCells().begin();
            iterLa!=largerSide.getSupportCells().end(); iterLa++)
        {
            
            if(*iterSm == *iterLa)
                supportDomainOverlap = true;
        }
    }
    bool smallerSideHasSupp = smallerSide.getSupportCells().size()!=0;
    bool largerSideHasSupp = largerSide.getSupportCells().size()!=0;
    bool bothSidesHaveSupp = largerSideHasSupp&&smallerSideHasSupp;
    
    scalar distOnNurbs = distance(smallerSide,largerSide);
    scalar distDirect = Foam::mag(smallerSide.markerPosition-largerSide.markerPosition);
    
    scalar smallerSideDomainMinSize = std::numeric_limits<scalar>::max();
    if(smallerSideHasSupp)
        smallerSideDomainMinSize = supportDomainMinSize(smallerSide.getSupportCells());
    
    scalar largerSideDomainMinSize = std::numeric_limits<scalar>::max();
    if(largerSideHasSupp)
        largerSideDomainMinSize = supportDomainMinSize(largerSide.getSupportCells());
    
    scalar overallDomainMinSide = std::min(smallerSideDomainMinSize,largerSideDomainMinSize);
    
    if(bothSidesHaveSupp)
    {
        if(distOnNurbs<overallDomainMinSide)
        {
            if(!supportDomainOverlap)
            {
                Info<<Foam::endl;
                std::cout<<"smallerSide.markerParameter:"<<smallerSide.getMarkerParameter()<<std::endl;
                std::cout<<"largerSide.markerParameter:"<<largerSide.getMarkerParameter()<<std::endl;
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
*/

/*
std::unique_ptr<std::pair<
    std::unordered_multimap<label,std::unordered_multimap<label,std::vector<vector>>>,
    std::unordered_multimap<label,std::unordered_multimap<label,std::vector<scalar>>>
>> Foam::LineStructure::shareMarkerData
(
    const std::vector<LagrangianMarker*>& markers,
    const std::vector<std::vector<vector>>& markersVecValues,
    const std::vector<std::vector<scalar>>& markersScalValues    
)
{
    const polyBoundaryMesh& boundaries = mesh.boundaryMesh();
    
    List<DynamicList<label>> procMarkersCell(Pstream::nProcs());
    List<DynamicList<vector>> procMarkersPosition(Pstream::nProcs());
    List<DynamicList<scalar>> procMarkersVolume(Pstream::nProcs());
    List<DynamicList<DynamicList<vector>>> procMarkersVecValues(Pstream::nProcs());
    List<DynamicList<DynamicList<scalar>>> procMarkersScalValues(Pstream::nProcs());

    for(label markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* marker = markers[markerInd];
        label markerCell = marker->getMarkerCell();
        if(selfHaloSet.find(markerCell)!=selfHaloSet.end())
        {
            procMarkersCell[Pstream::myProcNo()].append(markerCell);
            procMarkersPosition[Pstream::myProcNo()].append(marker->getMarkerPosition());
            procMarkersVolume[Pstream::myProcNo()].append(marker->getMarkerVolume());
            if(markersVecValues.size()>0)
            {
                procMarkersVecValues[Pstream::myProcNo()].append(DynamicList<vector>());
                for(vector vecVal : markersVecValues[markerInd])
                    procMarkersVecValues[Pstream::myProcNo()].back().append(vecVal);
            }
            if(markersScalValues.size()>0)
            {
                procMarkersScalValues[Pstream::myProcNo()].append(DynamicList<scalar>());
                for(scalar scalVal : markersScalValues[markerInd])
                    procMarkersScalValues[Pstream::myProcNo()].back().append(scalVal);
            }
        }
    }
    
    
    Pstream::gather(procMarkersCell);
    Pstream::gather(procMarkersPosition);
    Pstream::gather(procMarkersVolume);
    Pstream::gather(procMarkersVecValues);
    Pstream::gather(procMarkersScalValues);
    
    Pstream::scatter(procMarkersCell);
    Pstream::scatter(procMarkersPosition);
    Pstream::scatter(procMarkersVolume);
    Pstream::scatter(procMarkersVecValues);
    Pstream::scatter(procMarkersScalValues);
    
    globalMarkerData.clear();
    
    auto result = std::unique_ptr<std::pair<
        std::unordered_multimap<label,std::unordered_multimap<label,std::vector<vector>>>,
        std::unordered_multimap<label,std::unordered_multimap<label,std::vector<scalar>>> 
    >>
    (
        new std::pair<
            std::unordered_multimap<label,std::unordered_multimap<label,std::vector<vector>>>,
            std::unordered_multimap<label,std::unordered_multimap<label,std::vector<scalar>>> 
        >();
    );
    
    for(label patchIndex=0; patchIndex<boundaries.size(); patchIndex++)
    {
        const polyPatch& patch = boundaries[patchIndex];
        if(isA<processorPolyPatch>(patch))
        {
            const processorPolyPatch* pPP = dynamic_cast<const processorPolyPatch*>(&patch);
            label neighborProcess = pPP->neighbProcNo();
            
            for(label markerListInd=0; markerListInd<procMarkersCell[neighborProcess].size(); markerListInd++)
            {
                label cellInd = procMarkersCell[neighborProcess][markerListInd];
                vector pos = procMarkersPosition[neighborProcess][markerListInd];
                scalar vol = procMarkersVolume[neighborProcess][markerListInd];
                
                globalMarkerData[neighborProcess][cellInd].append(pos,vol);
                
                DynamicList<vector>& oneMarkerVecValues = procMarkersVecValues[neighborProcess][markerListInd];
                result->first[neighborProcess][cellInd].append(std::vector<vector>());
                for(vector& vec : oneMarkerVecValues)
                    result->first[neighborProcess][cellInd].back().push_back(vec);
                
                DynamicList<scalar>& oneMarkerScalValues = procMarkersScalValues[neighborProcess][markerListInd];
                result->second[neighborProcess][cellInd].append(std::vector<scalar>());
                for(scalar& scal : oneMarkerScalValues)
                    result->first[neighborProcess][cellInd].back().push_back(scal);                
            }
        }
    }
    return result;    
}
*/

scalar Foam::LineStructure::initialSpacingFromMesh()
{
    const cell& oneCell = mesh.cells()[0];
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    const pointField& oneCellPoints = oneCell.points(faces,points);
    vector oneCellCentre = oneCell.centre(points,faces);
    scalar minHalfDist = std::numeric_limits<scalar>::max();
    for(vector pnt : oneCellPoints)
    {
        vector conn = pnt-oneCellCentre;
        for(label dim=0; dim<3; dim++)
        {
            scalar connD = std::abs(conn[dim]);
            if(connD<minHalfDist)
                minHalfDist=connD;
        }
    }
    return 2*minHalfDist;
}
