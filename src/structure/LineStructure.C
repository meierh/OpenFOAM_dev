#include "LineStructure.H"
#include <Structure.H>

Foam::LineStructure::LineStructure
(
    const fvMesh& mesh,
    const List<scalar> crossSecArea,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
Structure(mesh,mesh.time()),
modusFieldToMarker(modusFieldToMarker),
modusMarkerToField(modusMarkerToField),
crossSecArea(crossSecArea)
{
    initialize();
}

Foam::LineStructure::LineStructure
(
    const fvMesh& mesh,
    const scalar crossSecArea,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
Structure(mesh,mesh.time()),
modusFieldToMarker(modusFieldToMarker),
modusMarkerToField(modusMarkerToField),
crossSecArea(List<scalar>(myMesh->m_nR,crossSecArea))
{
    initialize();
}

Foam::LineStructure::LineStructure
(
    const fvMesh& mesh,
    const std::shared_ptr<IOdictionary> structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
Structure(mesh,structureDict,mesh.time()),
modusFieldToMarker(modusFieldToMarker),
modusMarkerToField(modusMarkerToField)
{
    readRodPntsToMeshSpacingDict(*structureDict);
    initialize();
}

Foam::LineStructure::LineStructure
(
    const fvMesh& mesh,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
Structure(mesh,mesh.time()),
modusFieldToMarker(modusFieldToMarker),
modusMarkerToField(modusMarkerToField)
{}

Foam::LineStructure::LineStructure
(
    const fvMesh& mesh,
    const std::shared_ptr<IOdictionary> structureDict,
    bool parentConstructor,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
Structure(mesh,structureDict,mesh.time()),
modusFieldToMarker(modusFieldToMarker),
modusMarkerToField(modusMarkerToField)
{}

void Foam::LineStructure::finalizeMarkers
(
    bool doRefine
)
{
    refineEvaluateReduceCollect(doRefine);
    markerWeighting();
}

void Foam::LineStructure::refineEvaluateReduceCollect
(
    bool doRefine
)
{
    auto start = std::chrono::system_clock::now();
    
    auto t1 = std::chrono::system_clock::now();
    if(doRefine)
        refineMarkers();
    auto t2 = std::chrono::system_clock::now();
    if(doRefine)
        Info<<"refineEvaluateReduceCollect / refineMarkers took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
    else
        Info<<"refineEvaluateReduceCollect / -- took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
    
    if(doRefine)
        setMarkerVolume();
    t1 = std::chrono::system_clock::now();
    if(doRefine)
        Info<<"refineEvaluateReduceCollect / setMarkerVolume took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count()<<" milliseconds"<<Foam::nl;    
    else
        Info<<"refineEvaluateReduceCollect / -- took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count()<<" milliseconds"<<Foam::nl;
    
    evaluateMarkerMeshRelation();
    t2 = std::chrono::system_clock::now();
    Info<<"refineEvaluateReduceCollect / evaluateMarkerMeshRelation took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
    
    reduceMarkers();
    t1 = std::chrono::system_clock::now();
    Info<<"refineEvaluateReduceCollect / reduceMarkers took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count()<<" milliseconds"<<Foam::nl;
    
    collectMarkers();
    t2 = std::chrono::system_clock::now();
    Info<<"refineEvaluateReduceCollect / collectMarkers: "<<collectedMarkers.size()<<" took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
    
    auto end = std::chrono::system_clock::now();
    Info<<"refineEvaluateReduceCollect took "<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()<<" milliseconds"<<Foam::nl;
}

void Foam::LineStructure::markerWeighting()
{
    auto start = std::chrono::system_clock::now();
    
    auto t2 = std::chrono::system_clock::now();
    computeMarkerCellWeights();
    auto t1 = std::chrono::system_clock::now();
    Info<<"markerWeighting / computeMarkerCellWeights took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count()<<" milliseconds"<<Foam::nl;
    
    collectHaloMarkers();
    t2 = std::chrono::system_clock::now();
    Info<<"markerWeighting / collectHaloMarkers took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
    
    exchangeHaloMarkersData();
    t1 = std::chrono::system_clock::now();
    Info<<"markerWeighting / exchangeHaloMarkersData took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count()<<" milliseconds"<<Foam::nl;
    
    computeMarkerWeights();
    t2 = std::chrono::system_clock::now();
    Info<<"markerWeighting / computeMarkerWeights took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
    
    exchangeHaloMarkersWeight();
    t1 = std::chrono::system_clock::now();
    Info<<"exchangeHaloMarkersWeight / computeMarkerWeights took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count()<<" milliseconds"<<Foam::nl;
    
    auto end = std::chrono::system_clock::now();
    Info<<"markerWeighting took "<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()<<" milliseconds"<<Foam::nl;
}

void Foam::LineStructure::moveMarkersOnRodMovement()
{
    auto t1 = std::chrono::system_clock::now();
    evaluateMarkerMeshRelation();
    auto t2 = std::chrono::system_clock::now();
    Info<<"evaluateMarkerMeshRelation took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
}

Foam::vector Foam::LineStructure::evaluateRodVelocity
(
    label rodNumber,
    scalar parameter,
    scalar angle,
    scalar radiusFrac
)
{
    scalar currentTime = mesh.time().value();
    vector currentPosition = evaluateRodPos(Rods[rodNumber],parameter);
    
    const std::pair<gsNurbs<scalar>,scalar>* prevDef = readPrevRodDeformation(rodNumber);
    const std::pair<gsNurbs<scalar>,scalar>* prevRot = readPrevRodRotation(rodNumber);
    
    if(prevDef==nullptr || prevRot==nullptr)
        return vector(0,0,0);
    
    scalar prevTime = prevDef->second;
    if(prevTime!=prevRot->second)
        FatalErrorInFunction<<"Mismatch in prev time stamps"<<exit(FatalError);
    scalar deltaT = currentTime-prevTime;
    
    vector previousPosition;
    rodEval(Rods[rodNumber]->m_Curve,prevDef->first,parameter,previousPosition);
    
    return (currentPosition-previousPosition)/deltaT;
}

void Foam::LineStructure::settleIntoRefinedMesh()
{
    auto start = std::chrono::system_clock::now();
    
    auto t1 = std::chrono::system_clock::now();
    evaluateMarkerMeshRelation();
    auto t2 = std::chrono::system_clock::now();
    Info<<"settleIntoRefinedMesh / evaluateMarkerMeshRelation took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
    
    auto end = std::chrono::system_clock::now();
    Info<<"settleIntoRefinedMesh took "<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()<<" milliseconds"<<Foam::nl;
}

void Foam::LineStructure::refineMarkersOnRefinedMesh()
{
    auto start = std::chrono::system_clock::now();
    
    auto t1 = std::chrono::system_clock::now();
    refineMarkers();
    auto t2 = std::chrono::system_clock::now();
    Info<<"refineMarkersOnRefinedMesh / refineMarkers took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
    
    setMarkerVolume();
    t1 = std::chrono::system_clock::now();
    Info<<"refineMarkersOnRefinedMesh / setMarkerVolume took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count()<<" milliseconds"<<Foam::nl;    
    
    evaluateMarkerMeshRelation();
    t2 = std::chrono::system_clock::now();
    Info<<"refineMarkersOnRefinedMesh / evaluateMarkerMeshRelation took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
    
    collectMarkers();
    t1 = std::chrono::system_clock::now();
    Info<<"refineMarkersOnRefinedMesh / collectMarkers: "<<collectedMarkers.size()<<" took "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count()<<" milliseconds"<<Foam::nl;
    
    auto end = std::chrono::system_clock::now();
    Info<<"refineMarkersOnRefinedMesh took "<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()<<" milliseconds"<<Foam::nl;
}

void Foam::LineStructure::to_string()
{
    Pout<<"-------LineStructure::Markers-------"<<Foam::nl;
    int count=0;
    for(std::unique_ptr<std::list<LagrangianMarker>>& oneRodMarkers : rodMarkersList)
    {
        Pout<<"Rod "<<count<<Foam::nl;
        if(oneRodMarkers)
        {
            for(const LagrangianMarker& marker : *oneRodMarkers)
            {
                Pout<<"\t"<<marker.to_string()<<Foam::nl;
            }
        }
    }
    Pout<<"------------------------------------"<<Foam::nl;
}

void Foam::LineStructure::setNurbsParameters
(
    label rodNumber,
    label derivCoeffNumber,
    label dimension,
    scalar value
)
{
    FatalErrorInFunction<<"Not implemented"<<exit(FatalError);
    //setNurbsCoeff(rodNumber,derivCoeffNumber,dimension,value);
}

Foam::vector Foam::LineStructure::dXdParam
(
    const LagrangianMarker* marker,
    const Parameter& par
)
{
    return dXdParam(marker->getRodNumber(),marker->getMarkerParameter(),0,0,par);
}

Foam::vector Foam::LineStructure::dXdParam
(
    label rodNumber,
    scalar rodParameter,
    scalar angle,
    scalar radiusFrac,
    const Parameter& par
)
{
    if(!(par.isValid()))
        FatalErrorInFunction<<"Invalid parameter here!"<<exit(FatalError);
    if(par.getType()!=Parameter::Type::Rod)
        FatalErrorInFunction<<"Invalid type of parameter here!"<<exit(FatalError);
    vector rodDerive(0,0,0);
    label dimension = par.getDimension();
    const std::vector<NurbsCoeffReference>& nurbsCoeffs = par.getNurbsCoeffs();
    for(const NurbsCoeffReference& ref : nurbsCoeffs)
    {
        if(rodNumber==ref.rodNumber)
        {
            if(dimension!=ref.dimension)
                FatalErrorInFunction<<"Dimension mismatch!"<<exit(FatalError);
            
            vector d1dC,d2dC,d3dC,rdC;
            rodEvalDerivCoeff(rodNumber,ref.coeffNumber,dimension,rodParameter,d1dC,d2dC,d3dC,rdC);
            rodDerive += rdC;
        }
    }
    return rodDerive;
}

void Foam::LineStructure::initialize()
{
    check();
    auto t1 = std::chrono::system_clock::now();
    
    createSpacingPoints();
    auto t2 = std::chrono::system_clock::now();
    Info<<"createSpacingPoints took:"<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
    
    createMarkersFromSpacedPoints();
    t1 = std::chrono::system_clock::now();
    Info<<"createMarkersFromSpacedPoints took:"<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count()<<" milliseconds"<<Foam::nl;
    
    finalizeMarkers();
}

void Foam::LineStructure::check()
{
    if(!myMesh)
        FatalErrorInFunction<<"Rod Mesh not set!"<<exit(FatalError);
    if(crossSecArea.size()!=myMesh->m_nR)
        FatalErrorInFunction<<"Mismatch in size of crossSecArea and rodMarkersList"<<exit(FatalError);
    if(static_cast<int>(myMesh->m_Rods.size())!=myMesh->m_nR)
        FatalErrorInFunction<<"Mismatch in size of m_Rods and m_nR"<<exit(FatalError);
}

void Foam::LineStructure::store()
{
    Structure::store();
    
    std::vector<LagrangianMarker>& timeMarkers = storage[mesh.time().value()];
    for(const LagrangianMarker* marker : collectedMarkers)
    {
        timeMarkers.push_back(*marker);
    }
}

void Foam::LineStructure::setToTime(scalar time)
{
    Structure::setToTime(time);
    
    if(storage.find(time)==storage.end())
        FatalErrorInFunction<<"No data exists for "<<time<<exit(FatalError);
    
    std::vector<LagrangianMarker>& timeMarkers = storage[time];
    collectedMarkers.clear();
    for(LagrangianMarker& marker : timeMarkers)
    {
        collectedMarkers.push_back(&marker);
    }
}

void Foam::LineStructure::createSpacingPoints()
{
    status.execValid(status.initialPoints);
    setupRodBoundingBoxTree();
    LineStructure::initialRodPoints.resize(myMesh->m_nR);
    for(int rodIndex=0; rodIndex<myMesh->m_nR; rodIndex++)
    {
        createSpacedPointsOnRod(rodIndex,initialMeshSpacing);
    }
    status.executed(status.initialPoints);
}

void Foam::LineStructure::createSpacedPointsOnRod
(
    label rodNumber,
    scalar spacing
)
{   
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    auto pointsPtr = std::unique_ptr<std::vector<scalar>>(new std::vector<scalar>());
    std::vector<scalar>& pointsVec = *pointsPtr;
    
    std::list<scalar> points;
    points.push_back(oneRod->m_Curve.domainStart());
    points.push_back(oneRod->m_Curve.domainEnd());
    
    bool refined=true;
    while(refined)
    {
        refined = false;
        //bool cond = true;
        scalar summedDist = 0;
        auto pntsIter0 = points.begin();
        auto pntsIter1 = ++(points.begin());
        for( ; pntsIter1!=points.end() ; )
        {
            scalar dist = distance(oneRod,*pntsIter0,*pntsIter1);
            summedDist+=dist;
            if(dist > spacing*iniRodPntsDistToMeshSpacing)
            {
                scalar middlePar = 0.5*(*pntsIter0 + *pntsIter1);
                points.insert(pntsIter1,middlePar);
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
    
    initialRodPoints[rodNumber] = std::move(pointsPtr);
}

void Foam::LineStructure::createMarkersFromSpacedPoints()
{
    status.execValid(status.markers);
    rodMarkersList.resize(myMesh->m_nR);
    for(uint rodIndex=0; rodIndex<rodMarkersList.size(); rodIndex++)
    {
        createMarkersFromSpacedPointsOnRod(rodIndex);
    }
    status.executed(status.markers);
}

void Foam::LineStructure::createMarkersFromSpacedPointsOnRod
(
    label rodNumber
)
{
    if(rodInMesh[rodNumber])
    {
        const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
        auto markersPtr = std::unique_ptr<std::list<LagrangianMarker>>(new std::list<LagrangianMarker>());
        std::list<LagrangianMarker>& markers = *markersPtr;
        markers.clear();
        for(scalar point : *(initialRodPoints[rodNumber]))
            markers.push_back(LagrangianMarker(*this,mesh,rodNumber,oneRod,point));
        rodMarkersList[rodNumber] = std::move(markersPtr);
    }
}

void Foam::LineStructure::refineMarkers
(
    bool useMarkerCharLenSpacing,
    std::pair<bool,scalar> forcedSpacing
)
{
    status.execValid(status.markersRefined);
    for(uint rodIndex=0; rodIndex<rodMarkersList.size(); rodIndex++)
    {
        refineMarkersOnRod(rodIndex,useMarkerCharLenSpacing,forcedSpacing);
    }
    status.executed(status.markersRefined);
}

void Foam::LineStructure::refineMarkersOnRod
(
    label rodNumber,
    bool useMarkerCharLenSpacing,
    std::pair<bool,scalar> forcedSpacing
)
{
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    if(rodMarkersList[rodNumber])
    {
        std::list<LagrangianMarker>& markers = *(rodMarkersList[rodNumber]);
        bool refined=true;
        while(refined)
        {
            refined = false;
            //bool cond = true;
            auto markersIter0 = markers.begin();
            auto markersIter1 = ++(markers.begin());
            for( ; markersIter1!=markers.end() ; )
            {
                scalar markers0Para = markersIter0->getMarkerParameter();
                scalar markers0CellSpacing = markersIter0->getMarkerCellMinSpacing();
                scalar markers0CharacSpacing = markersIter0->getMarkerCharacLen();
                bool markers0InCell = (markersIter0->getMarkerCell()!=-1);
                
                scalar markers1Para = markersIter1->getMarkerParameter();
                scalar markers1CellSpacing = markersIter1->getMarkerCellMinSpacing();
                scalar markers1CharacSpacing = markersIter1->getMarkerCharacLen();
                bool markers1InCell = (markersIter1->getMarkerCell()!=-1);
                
                scalar dist = distance(oneRod,markers0Para,markers1Para);
                
                bool subdivide = false;
                if(forcedSpacing.first)
                {
                    if(dist>forcedSpacing.second)
                        subdivide=true;
                }
                if(markers0InCell || markers1InCell)
                {
                    scalar minSpacing = std::min(markers0CellSpacing,markers1CellSpacing);
                    if(useMarkerCharLenSpacing)
                    {
                        if(markers0InCell)
                            minSpacing = std::min(minSpacing,markers0CharacSpacing*rodPntDistToMarkerCharLen);
                        if(markers1InCell)
                            minSpacing = std::min(minSpacing,markers1CharacSpacing*rodPntDistToMarkerCharLen);
                    }
                    if(dist > minSpacing*refnRodMarkersDistToMeshSpacing)
                        subdivide=true;
                }
            
                if(subdivide)
                {
                    scalar middlePar = 0.5*(markers0Para+markers1Para);
                    LagrangianMarker middleMarker(*this,mesh,rodNumber,oneRod,middlePar);
                    markers.insert(markersIter1,middleMarker);
                    refined=true;
                }
                markersIter0 = markersIter1;
                markersIter1++;
            }
        }
    }
}

void Foam::LineStructure::setMarkerVolume()
{
    status.execValid(status.markersVolume);
    for(uint rodIndex=0; rodIndex<rodMarkersList.size(); rodIndex++)
    {
        setMarkerVolumeOnRod(rodIndex);
    }
    status.executed(status.markersVolume);
}

void Foam::LineStructure::setMarkerVolumeOnRod
(
    label rodNumber
)
{
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    if(rodMarkersList[rodNumber])
    {
        std::list<LagrangianMarker>& markers = *(rodMarkersList[rodNumber]);
        
        std::function<bool(scalar)> InMesh = 
        [&,rod=oneRod](scalar parameter)
        {
            vector position = evaluateRodPos(rod,parameter);
            label posCell = mesh.findCell(position);
            if(posCell==-1)
                return false;
            else
                return true;
        };

        std::list<LagrangianMarker>::iterator iterPrev = markers.end();
        std::list<LagrangianMarker>::iterator iterNext;
        for(auto iter=markers.begin(); iter!=markers.end(); iter++)
        {
            scalar span = 0;
            if(iter->getMarkerCell()!=-1)
            {
                iterNext = iter;
                iterNext++;
                scalar spanStart;
                if(iterPrev!=markers.end())
                {
                    if(iterPrev->getMarkerCell()==-1)
                    {
                        scalar prevPara = iterPrev->getMarkerParameter();
                        scalar thisPara = iter->getMarkerParameter();
                        if(!InMesh(thisPara))
                            FatalErrorInFunction<<"Both sides out of mesh"<<exit(FatalError);
                        scalar threshold = spacingFromMesh(mesh,iter->getMarkerCell())/100;
                        spanStart = bisectionBinary(prevPara,thisPara,InMesh,threshold);
                    }
                    else
                    {
                        spanStart = iterPrev->getMarkerParameter();
                        spanStart = spanStart + iter->getMarkerParameter();
                        spanStart /= 2;
                    }
                }
                else
                    spanStart = iter->getMarkerParameter();

                scalar spanEnd;
                if(iterNext!=markers.end())
                {
                    if(iterNext->getMarkerCell()==-1)
                    {
                        scalar nextPara = iterNext->getMarkerParameter();
                        scalar thisPara = iter->getMarkerParameter();
                        if(!InMesh(thisPara))
                            FatalErrorInFunction<<"Both sides out of mesh"<<exit(FatalError);
                        scalar threshold = spacingFromMesh(mesh,iter->getMarkerCell())/100;
                        spanEnd = bisectionBinary(thisPara,nextPara,InMesh,threshold);
                    }
                    else
                    {
                        spanEnd = iterNext->getMarkerParameter();
                        spanEnd = spanEnd + iter->getMarkerParameter();
                        spanEnd /= 2;
                    }
                }
                else
                    spanEnd = iter->getMarkerParameter();
                
                if(iterPrev==markers.end() && iterNext==markers.end())
                    FatalErrorInFunction<<"Marker with no predecessor and no succesor"<<exit(FatalError);
                
                span = LineStructure::distance(oneRod,spanStart,spanEnd);
                //Pout<<"Cell:"<<iter->getMarkerCell()<<":"<<iter->getMarkerParameter()<<"  start:"<<spanStart<<"  end:"<<spanEnd<<"  span:"<<span<<Foam::nl;
            }
            if(modusMarkerToField==markerMeshType::NonUniform)
            {

                iter->setMarkerVolume(span);
            }
            else
                iter->setMarkerVolume(span*crossSecArea[rodNumber]);
            
            iterPrev = iter;
        }
    }
}

void Foam::LineStructure::evaluateMarkerMeshRelation()
{
    status.execValid(status.markerMesh);
    for(std::unique_ptr<std::list<LagrangianMarker>>& singleRodMarkers : rodMarkersList)
        if(singleRodMarkers)
            evaluateMarkerMeshRelation(*singleRodMarkers);
    status.executed(status.markerMesh);
}

void Foam::LineStructure::evaluateMarkerMeshRelation
(
    std::list<LagrangianMarker>& markerList
)
{
    for(LagrangianMarker& marker : markerList)
        marker.evaluateMarker();
}

void Foam::LineStructure::reduceMarkers()
{
    std::cout<<"LineStructure::reduceMarkers()"<<std::endl;

    status.execValid(status.markersReduction);
    std::vector<MarkerReference<LagrangianMarker>> allMarkers;
    for(std::unique_ptr<std::list<LagrangianMarker>>& singleRodMarkersPtr : rodMarkersList)
    {
        if(singleRodMarkersPtr)
        {
            std::list<LagrangianMarker>* singleRodMarkers = &(*singleRodMarkersPtr);
            for(auto iter=singleRodMarkers->begin(); iter!=singleRodMarkers->end(); iter++)
            {
                allMarkers.push_back(MarkerReference<LagrangianMarker>(iter,singleRodMarkers));
            }
        }
    }
    reduceMarkers(allMarkers);
    removeOverlapMarkers();    
    status.executed(status.markersReduction);
}

void Foam::LineStructure::removeOverlapMarkers()
{}

void Foam::LineStructure::collectMarkers()
{
    status.execValid(status.markersCollected);    
    collectedMarkers.resize(0);
    for(std::size_t rodNumber=0; rodNumber<rodMarkersList.size(); rodNumber++)
    {
        std::unique_ptr<std::list<LagrangianMarker>>& oneRodMarkers = rodMarkersList[rodNumber];
        if(rodInMesh[rodNumber])
        {
            if(oneRodMarkers)
            {
                for(LagrangianMarker& oneMarker : *oneRodMarkers)
                {
                    collectedMarkers.push_back(&oneMarker);
                }
            }
            else
            {
                meshBoundingBox.print();
                Pout<<"rodInMesh["<<rodNumber<<"]:"<<rodInMesh[rodNumber]<<Foam::nl;
                Pout<<"rodTrees["<<rodNumber<<"]:"; rodTrees[rodNumber].printRoot();
                FatalErrorInFunction<<"Rod with no markers given but in mesh"<<exit(FatalError);
            }
        }
        else
        {
            if(oneRodMarkers)
                FatalErrorInFunction<<"Rod out of mesh but markers given"<<exit(FatalError);
        }
    }
    status.executed(status.markersCollected);
}

void Foam::LineStructure::collectHaloMarkers()
{
    status.execValid(status.markersHaloCollect);
    haloCellsRodMarkersList.clear();
    haloCellsRodMarkersList.resize(getGlobalHaloCellList_Sorted(Pstream::myProcNo()).size());
    const std::unordered_map<label,label>& selfHaloCellToIndex = getHaloCellToIndexMap(Pstream::myProcNo());
    for(uint localMarkerInd=0; localMarkerInd<collectedMarkers.size(); localMarkerInd++)
    {
        LagrangianMarker* marker = collectedMarkers[localMarkerInd];
        label cellOfMarker = marker->getMarkerCell();
        auto iter = selfHaloCellToIndex.find(cellOfMarker);
        if(iter!=selfHaloCellToIndex.end())
        {
            label index = iter->second;
            std::pair<LagrangianMarker*,label> markerData(marker,localMarkerInd);
            haloCellsRodMarkersList[index].push_back(markerData);
        }
    }
    status.executed(status.markersHaloCollect);
}

void Foam::LineStructure::exchangeHaloMarkersData()
{
    status.execValid(status.markersHaloExchange);

    globHaloMarkers = GlobalHaloMarkers(this,haloCellsRodMarkersList);
    for(uint haloCellInd=0; haloCellInd<haloCellsRodMarkersList.size(); haloCellInd++)
    {
        std::vector<std::pair<LagrangianMarker*,label>>& localHaloCellRodMarkers = haloCellsRodMarkersList[haloCellInd];
        for(uint haloMarkerInd=0; haloMarkerInd<haloCellsRodMarkersList[haloCellInd].size(); haloMarkerInd++)
        {
            std::pair<LagrangianMarker*,label> markerData = localHaloCellRodMarkers[haloMarkerInd];
            LagrangianMarker* marker = markerData.first;
            //Pout<<"Append:"<<marker->to_string()<<Foam::nl;
            label markerInd = markerData.second;
            
            vector position = marker->getMarkerPosition();
            scalar volume = marker->getMarkerVolume();
            label index = markerInd;
            vector dilation = marker->getDilation();
            DynamicList<Pair<label>> supportCellsIndices;
            DynamicList<vector> supportCellsCentre;
            DynamicList<scalar> supportCellsVolume;
            label nurbsInd = marker->getRodNumber();
            scalar para = marker->getMarkerParameter();
            scalar angle = marker->getMarkerAngle();
            scalar radiusFrac = marker->getMarkerRadiusFrac();
            
            const DynamicList<Pair<label>>& supportCells = marker->getSupportCells();
            for(Pair<label> oneSupport : supportCells)
            {
                vector oneSuppCellCentre;
                scalar oneSuppCellVolume;
                marker->getCellData(oneSupport,oneSuppCellCentre,oneSuppCellVolume);
                supportCellsIndices.append({oneSupport.first(),oneSupport.second()});
                supportCellsCentre.append(oneSuppCellCentre);
                supportCellsVolume.append(oneSuppCellVolume);
            }
            globHaloMarkers.appendMarkerData
            (
                haloCellInd,
                {position,volume,index,dilation,supportCellsIndices,supportCellsCentre,supportCellsVolume,marker->getCorrParaB(),nurbsInd,para,angle,radiusFrac}
            );
        }
    }
    globHaloMarkers.communicate();

    status.executed(status.markersHaloExchange);
}

void Foam::LineStructure::computeMarkerCellWeights()
{
    status.execValid(status.markersCellWeight);
    for(LagrangianMarker* marker : collectedMarkers)
    {
        marker->compNonUniformCorrWeights();
    }
    status.executed(status.markersCellWeight);
}

std::unique_ptr<Foam::LineStructure::LinearSystem> Foam::LineStructure::computeMarkerEpsilonMatrix()
{   
    label locProcMarkerNbr = collectedMarkers.size();

    List<label> globalMarkerNumber(Pstream::nProcs());
    globalMarkerNumber[Pstream::myProcNo()] = locProcMarkerNbr;
    Pstream::gatherList(globalMarkerNumber);
    Pstream::scatterList(globalMarkerNumber);
        
    label globNumOfMarkers = 0;
    for(label proc=0; proc<globalMarkerNumber.size(); proc++)
        globNumOfMarkers += globalMarkerNumber[proc];
    label smProcsNumOfMarkers = 0;
    for(label proc=0; proc<Pstream::myProcNo(); proc++)
        smProcsNumOfMarkers += globalMarkerNumber[proc];
    List<label> processGlobalMarkerOffset(Pstream::nProcs());
    processGlobalMarkerOffset[0] = 0;
    for(label proc=1; proc<processGlobalMarkerOffset.size(); proc++)
        processGlobalMarkerOffset[proc] = processGlobalMarkerOffset[proc-1]+globalMarkerNumber[proc-1];
            
    auto result = std::unique_ptr<LinearSystem>(new LinearSystem());
    std::get<0>(*result) = CSR_Matrix_par(locProcMarkerNbr,smProcsNumOfMarkers,globNumOfMarkers,globNumOfMarkers);
    CSR_Matrix_par& A = std::get<0>(*result);    
    std::get<1>(*result) = Vector_par(locProcMarkerNbr,smProcsNumOfMarkers,globNumOfMarkers);
    Vector_par& b = std::get<1>(*result);
    
    std::unordered_map<Pair<label>,DynamicList<label>,foamPairHash<label>> cellToMarkerInfluence;
    for(uint I=0; I<collectedMarkers.size(); I++)
    {
        const LagrangianMarker& markerI = *(collectedMarkers[I]);
        //cellToMarkerI.insert({markerI.getMarkerCell(),I});
        const DynamicList<Pair<label>>& supportI = markerI.getSupportCells();
        for(const Pair<label>& node : supportI)
        {
            cellToMarkerInfluence[node].append(I);
        }
    }
        
    struct VectorHash
    {
        std::size_t operator()(const vector& vec) const noexcept
        {
            std::size_t vec0 = std::hash<scalar>{}(vec[0]);
            std::size_t vec1 = std::hash<scalar>{}(vec[1]);
            std::size_t vec2 = std::hash<scalar>{}(vec[2]);
            return vec0 ^ vec1 ^ vec2;
        }
    };
    using CellSet = std::unordered_map<vector,scalar,VectorHash>;
    for(uint I=0; I<collectedMarkers.size(); I++)
    {
        const LagrangianMarker& markerI = *(collectedMarkers[I]);
        vector XI = markerI.getMarkerPosition();
        vector dilationI = markerI.getDilation();
        scalar dilationIMax = std::max<scalar>(dilationI[0],dilationI[1]);
        dilationIMax = std::max<scalar>(dilationIMax,dilationI[2]);
        //scalar markerIVol = markerI.getMarkerVolume();
        const DynamicList<Pair<label>>& supportI = markerI.getSupportCells();
        CellSet markerISupportMap;
        for(const Pair<label>& suppCell : supportI)
        {
            vector centre;
            scalar vol;
            markerI.getCellData(suppCell,centre,vol);
            markerISupportMap[centre] = vol;
        }

        std::unordered_set<label> localNeighbourMarkerK;
        for(const Pair<label>& suppCell : supportI)
        {
            auto iter = cellToMarkerInfluence.find(suppCell);
            if(iter!=cellToMarkerInfluence.end())
            {
                for(scalar markerK : iter->second)
                    localNeighbourMarkerK.insert(markerK);
            }
        }

        std::map<label,scalar> rowEntries;
        // Compute local matrix entries
        for(auto iterK=localNeighbourMarkerK.begin(); iterK!=localNeighbourMarkerK.end(); iterK++)
        {
            uint K = *iterK;
        //for(uint K=0; K<collectedMarkers.size(); K++)
        //{
            /*
            if(Pstream::myProcNo()==0)
                Pout<<"I:"<<I<<" K:"<<K<<Foam::nl;
            */

            const LagrangianMarker& markerK = *(collectedMarkers[K]);
            vector XK = markerK.getMarkerPosition();
            vector dilationK = markerK.getDilation();
            scalar dilationKMax = std::max<scalar>(dilationK[0],dilationK[1]);
            dilationKMax = std::max<scalar>(dilationKMax,dilationK[2]);
            scalar markerKVol = markerK.getMarkerVolume();
            const DynamicList<Pair<label>>& supportK = markerK.getSupportCells();
            
            vector distVec = XI-XK;
            scalar distMag = std::sqrt(distVec&distVec);
            
            scalar matrixEntry = 0;
            if(distMag < 10*dilationIMax || distMag < 10*dilationKMax)
            {
                DynamicList<std::pair<vector,scalar>> combinedIKSupport;
                for(const Pair<label>& suppCell : supportK)
                {
                    vector centre;
                    scalar vol;
                    markerK.getCellData(suppCell,centre,vol);
                    if(markerISupportMap.find(centre)!=markerISupportMap.end())
                    {
                        combinedIKSupport.append({centre,vol});
                    }
                }
                for(std::pair<vector,scalar> suppCell : combinedIKSupport)
                {
                    scalar weightI = markerI.correctedDeltaDirac(XI,suppCell.first);
                    scalar weightK = markerK.correctedDeltaDirac(XK,suppCell.first);                       
                    matrixEntry += weightK*weightI*suppCell.second;
                    /*
                    if(Pstream::myProcNo()==0)
                        Pout<<"    ("<<weightI<<","<<weightK<<") "<<"XI:"<<XI<<" XK:"<<XK<<" suppC:"<<suppCell.first<<" suppV:"<<suppCell.second<<Foam::nl;
                    */

                }
                matrixEntry *= markerKVol;
            }
            
            auto iter = rowEntries.find(smProcsNumOfMarkers+K);
            if(iter!=rowEntries.end())
                FatalErrorInFunction<<"Multiple writes in matrix"<<exit(FatalError);
            rowEntries[smProcsNumOfMarkers+K] = matrixEntry;
            /*
            if(Pstream::myProcNo()==0)
                Pout<<"Inner ("<<matrixEntry<<","<<(smProcsNumOfMarkers+K)<<")   markerKVol:"<<markerKVol<<Foam::nl;
            */
        }
                
        // Compute foreign local matrix entries
        std::unordered_set<label> neighborProcesses;
        for(const Pair<label>& tupl : supportI)
        {
            if(!(tupl.first()==Pstream::myProcNo()))
            {
                label process = tupl.first();
                label cellInd = tupl.second();
                neighborProcesses.insert(process);
                const std::unordered_map<label,label>& processHaloCellToIndexMap = getHaloCellToIndexMap(process);
                auto cellIndexIter = processHaloCellToIndexMap.find(cellInd);
                if(cellIndexIter==processHaloCellToIndexMap.end())
                    FatalErrorInFunction<<"Cell in process halo not found"<<exit(FatalError);
                label index = cellIndexIter->second;
                const DynamicList<CellDescription>& processHaloCellList = getGlobalHaloCellList_Sorted(process);
                const CellDescription& oneCell = processHaloCellList[index];
                if(oneCell.index!=cellInd)
                    FatalErrorInFunction<<"Cell index data mismatch"<<exit(FatalError);
            }
        }
        using LM=LagrangianMarker;
        for(auto iter=neighborProcesses.begin(); iter!=neighborProcesses.end(); iter++)
        {
            label process = *iter;
            label offset = processGlobalMarkerOffset[process];
            /*
            if(Pstream::myProcNo()==0)
                Pout<<"neighbour:"<<process<<Foam::nl;
            */
            label neighProcHaloCellNum = globHaloMarkers.size_haloCells(process);
            for(label haloCellInd=0; haloCellInd<neighProcHaloCellNum; haloCellInd++)
            {
                //Pout<<"neighbour:"<<process<<"  haloCellInd:"<<haloCellInd<<"/"<<neighProcHaloCellNum<<Foam::nl;
                label haloCellMarkerNum = globHaloMarkers.size_cellMarkers(process,haloCellInd);
                for(label haloCellMarkerInd=0; haloCellMarkerInd<haloCellMarkerNum; haloCellMarkerInd++)
                {
                    //Pout<<"neighbour:"<<process<<"  haloCellInd:"<<haloCellInd<<"  haloCellMarkerInd:"<<haloCellMarkerInd<<"/"<<haloCellMarkerNum<<Foam::nl;
                    
                    std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,FixedList<scalar,10>,label,scalar,scalar,scalar> markerData;
                    markerData = globHaloMarkers.getMarkerData(process,haloCellInd,haloCellMarkerInd);
                    
                    vector markerKposition = std::get<0>(markerData);
                    scalar markerKvolume = std::get<1>(markerData);
                    label markerKindex = std::get<2>(markerData);
                    vector markerKdilation = std::get<3>(markerData);
                    scalar dilationKMax = std::max<scalar>(markerKdilation[0],markerKdilation[1]);
                    dilationKMax = std::max<scalar>(dilationKMax,markerKdilation[2]);
                    
                    vector XK = markerKposition;
                    vector distVec = XI-XK;
                    scalar distMag = std::sqrt(distVec&distVec);
                    scalar matrixEntry = 0;
                    
                    
                    if(distMag < 10*dilationIMax || distMag < 10*dilationKMax)
                    {
                        DynamicList<vector>& markerKSuppCellCentres = std::get<5>(markerData);
                        DynamicList<scalar>& markerKSuppCellVol = std::get<6>(markerData);
                        FixedList<scalar,10>& markerKb = std::get<7>(markerData);
                        
                        DynamicList<std::pair<vector,scalar>> combinedIKSupport;
                        for(label suppKInd=0; suppKInd<markerKSuppCellCentres.size(); suppKInd++)
                        {
                            vector centre = markerKSuppCellCentres[suppKInd];
                            scalar vol = markerKSuppCellVol[suppKInd];
                            if(markerISupportMap.find(centre)!=markerISupportMap.end())
                            {
                                combinedIKSupport.append({centre,vol});
                            }
                        }
                        
                        for(std::pair<vector,scalar> suppCell : combinedIKSupport)
                        {
                            //Pout<<"suppCell:"<<suppCell.first<<"/"<<suppCell.second<<Foam::nl;
                            scalar weightI = markerI.correctedDeltaDirac(XI,suppCell.first);
                            scalar weightK = LM::correctedDeltaDirac
                            (
                                XK,suppCell.first,markerKdilation,markerKb
                            );
                            matrixEntry += weightK*weightI*suppCell.second;
                        }
                        matrixEntry *= markerKvolume;
                    }
                    
                    auto iter = rowEntries.find(offset+markerKindex);
                    if(iter!=rowEntries.end())
                        FatalErrorInFunction<<"Multiple writes in matrix"<<exit(FatalError);
                    rowEntries[offset+markerKindex] = matrixEntry;
                    /*
                    if(Pstream::myProcNo()==0)
                        Pout<<"Outer ("<<matrixEntry<<","<<(offset+markerKindex)<<") "<<Foam::nl;
                    */
                }
            }
        }
        
        DynamicList<std::pair<scalar,label>> matrixRow;
        for(auto iter=rowEntries.begin(); iter!=rowEntries.end(); iter++)
        {
            uint K = iter->first;
            scalar matrixEntry = iter->second;
            matrixRow.append({matrixEntry,K});
            //Check for small diagonal entries
            if(K==I+smProcsNumOfMarkers && std::abs(matrixEntry)<1e-5)
            {
                Pout<<"I:"<<I<<" K:"<<K<<" matrixEntry:"<<matrixEntry<<Foam::nl;
                Pout<<"markerI:"<<markerI.to_string()<<Foam::nl;
                FatalErrorInFunction<<"Near zero diagonal element"<<exit(FatalError);
            }
        }
        //Check for non existing diagonal entries
        if(rowEntries.find(I+smProcsNumOfMarkers)==rowEntries.end())
        {
            Pout<<"matrixRow:";
            for(auto entr : matrixRow)
                Pout<<"("<<entr.first<<","<<entr.second<<") ";
            Pout<<Foam::nl;
            FatalErrorInFunction<<"Missing therfore zero diagonal element"<<exit(FatalError);
        }
        A.addRow(matrixRow);
        b[I] = 1;
        /*
        if(true)
        {
            for(auto pair : matrixRow)
                Pout<<"("<<pair.first<<","<<pair.second<<") ";
            Pout<<Foam::nl;
        }
        */
    }

    return result;
}

void Foam::LineStructure::computeMarkerWeights()
{
    status.execValid(status.markersWeight);
    
    std::unique_ptr<Foam::LineStructure::LinearSystem> system = computeMarkerEpsilonMatrix();
    CSR_Matrix_par& A = std::get<0>(*system);
    Vector_par& ones = std::get<1>(*system);
    //Pout<<"Matrix assembly took :"<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<Foam::nl;
        
    switch (solutionStrategy)
    {
        case SystemSolve::Raw:
        {
            break;
        }
        case SystemSolve::RowAequilibration:
        {
            Vector_par rowSum = A*ones;
            for(label localRow=0; localRow<rowSum.getLocalSize().second; localRow++)
            {
                if(rowSum[localRow]==0)
                    FatalErrorInFunction<<"Matrix has a zero in the diagonal"<<exit(FatalError);
                rowSum[localRow] = 1.0/rowSum[localRow];
            }
            CSR_DiagMatrix_par P(rowSum);
            A = P*A;
            ones = P*ones;
            break;
        }
        case SystemSolve::ColAequilibration:
        {
            FatalErrorInFunction<<"Not implemented!"<<exit(FatalError);
            break;
        }
        case SystemSolve::GeoRescaled:
        {
            FatalErrorInFunction<<"Not matching!"<<exit(FatalError);
            break;
        }
        case SystemSolve::Jacobi:
        {
            CSR_DiagMatrix_par diagA = A.diagonalMatrix();
            //Pout<<diagA.to_string()<<Foam::nl;
            Vector_par diagAVec = diagA*ones;
            //Pout<<diagAVec.to_string()<<Foam::nl;
            for(label localRow=0; localRow<diagAVec.getLocalSize().second; localRow++)
            {
                if(diagAVec[localRow]==0)
                    FatalErrorInFunction<<"Matrix has a zero in the diagonal"<<exit(FatalError);
                diagAVec[localRow] = 1.0/diagAVec[localRow];
            }
            CSR_DiagMatrix_par P(diagAVec);
            A = P*A;
            ones = P*ones;        
            break;
        }
        default:
            FatalErrorInFunction<<"Invalid option"<<exit(FatalError);
    }
    BiCGSTAB solver(A);
    Vector_par eps = solver.solve(ones);

    //Pout<<"Linear system solution took :"<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count()<<" milliseconds"<<Foam::nl;
        
    for(uint I=0; I<collectedMarkers.size(); I++)
    {
        collectedMarkers[I]->setMarkerWeight(eps[I]);
    }
    status.executed(status.markersWeight);
}

void Foam::LineStructure::exchangeHaloMarkersWeight()
{    
    status.execValid(status.markersWeightExchange);
    
    for(uint haloCellInd=0; haloCellInd<haloCellsRodMarkersList.size(); haloCellInd++)
    {
        std::vector<std::pair<LagrangianMarker*,label>>& localHaloCellRodMarkers = haloCellsRodMarkersList[haloCellInd];
        for(uint haloMarkerInd=0; haloMarkerInd<haloCellsRodMarkersList[haloCellInd].size(); haloMarkerInd++)
        {
            std::pair<LagrangianMarker*,label> markerData = localHaloCellRodMarkers[haloMarkerInd];
            LagrangianMarker* marker = markerData.first;
            scalar weight = marker->getMarkerWeight();
            globHaloMarkers.insertMarkerWeight(haloCellInd,haloMarkerInd,weight);
        }
    }
    globHaloMarkers.communicateWeight();
    
    status.executed(status.markersWeightExchange);
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

Foam::vector Foam::LineStructure::derivateRodPos
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter
)
{
    vector dd1dp,dd2dp,dd3dp,drdp;
    rodDerivEval(oneRod,parameter,dd1dp,dd2dp,dd3dp,drdp);
    return drdp;
}

Foam::vector Foam::LineStructure::derivate2RodPos
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter
)
{
    vector d2d1dp,d2d2dp,d2d3dp,d2rdp;
    rodDeriv2Eval(oneRod,parameter,d2d1dp,d2d2dp,d2d3dp,d2rdp);
    return d2rdp;
}

Foam::scalar Foam::LineStructure::distance
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

Foam::scalar Foam::LineStructure::distance
(
    const LagrangianMarker& A,
    const LagrangianMarker& B
)
{
    if(A.getBaseRod()!=B.getBaseRod())
        FatalErrorInFunction<<"Distance can not be computed between points on different rods"<<exit(FatalError);
    
    return distance(A.getBaseRod(),A.getMarkerParameter(),B.getMarkerParameter());
}

Foam::scalar Foam::LineStructure::bisectionBinary
(
    scalar startValue,
    scalar endValue,
    std::function<bool(scalar)> criterion,
    scalar threshold
)
{
    bool startCriterion = criterion(startValue);
    bool endCriterion = criterion(endValue);
    if(startCriterion==endCriterion)
        FatalErrorInFunction<<"Invalid bisection start values"<<exit(FatalError);
    
    scalar dist;
    do
    {
        scalar center = 0.5*(startValue+endValue);
        bool centerCriterion = criterion(center);
        if(centerCriterion==startCriterion)
        {
            startValue = center;
        }
        else if(centerCriterion==endCriterion)
        {
            endValue = center;
        }
        dist = std::abs(endValue-startValue);
    }
    while(dist>=threshold);
    
    return 0.5*(startValue+endValue);
}

void Foam::LineStructure::readRodPntsToMeshSpacingDict
(
    const IOdictionary& structureDict
)
{
    ITstream iniSpacingFactorStream = structureDict.lookup("iniPntDistToCellSpacing");
    token iniSpacingFactorToken;
    iniSpacingFactorStream.read(iniSpacingFactorToken);
    if(!iniSpacingFactorToken.isScalar())
    {
        Pout<<"iniSpacingFactorToken:"<<iniSpacingFactorToken<<Foam::endl;
        Pout<<"iniSpacingFactorToken:"<<iniSpacingFactorToken.typeName()<<Foam::endl;
        FatalErrorInFunction<<"Invalid entry in constant/structureDict/iniPntDistToCellSpacing -- must be scalar"<<exit(FatalError);
    }
    iniRodPntsDistToMeshSpacing = iniSpacingFactorToken.scalarToken();
    
    ITstream refnSpacingFactorStream = structureDict.lookup("refnPntDistToCellSpacing");
    token refnSpacingFactorToken;
    refnSpacingFactorStream.read(refnSpacingFactorToken);
    if(!refnSpacingFactorToken.isScalar())
    {
        Pout<<"refnSpacingFactorToken:"<<refnSpacingFactorToken<<Foam::endl;
        Pout<<"refnSpacingFactorToken:"<<refnSpacingFactorToken.typeName()<<Foam::endl;
        FatalErrorInFunction<<"Invalid entry in constant/structureDict/refnPntDistToCellSpacing -- must be scalar"<<exit(FatalError);
    }
    refnRodMarkersDistToMeshSpacing = refnSpacingFactorToken.scalarToken();
    
    ITstream pntDistToMarkerCharLenStream = structureDict.lookup("pntDistToMarkerCharLen");
    token pntDistToMarkerCharLenToken;
    pntDistToMarkerCharLenStream.read(pntDistToMarkerCharLenToken);
    if(!pntDistToMarkerCharLenToken.isScalar())
    {
        Pout<<"pntDistToMarkerCharLenToken:"<<pntDistToMarkerCharLenToken<<Foam::endl;
        Pout<<"pntDistToMarkerCharLenToken:"<<pntDistToMarkerCharLenToken.typeName()<<Foam::endl;
        FatalErrorInFunction<<"Invalid entry in constant                                                            /structureDict/pntDistToMarkerCharLen -- must be scalar"<<exit(FatalError);
    }
    rodPntDistToMarkerCharLen = pntDistToMarkerCharLenToken.scalarToken();
    
    Info<<"iniRodPntsDistToMeshSpacing:"<<iniRodPntsDistToMeshSpacing<<Foam::nl;
    Info<<"refnRodMarkersDistToMeshSpacing:"<<refnRodMarkersDistToMeshSpacing<<Foam::nl;
    Info<<"rodPntDistToMarkerCharLen:"<<rodPntDistToMarkerCharLen<<Foam::nl;
}

Foam::BoundingBox Foam::LineStructure::computeBox
(
    label rodNumber
) const
{
    const ActiveRodMesh::rodCosserat* rod = Rods[rodNumber];
    const gsNurbs<scalar>& curve = rod->m_Curve;
    BoundingBox curve_box = BoundingBox::boundsOfNurbs(curve);
    const gsNurbs<scalar>& deformation = rod->m_Def;
    BoundingBox def_box = BoundingBox::boundsOfNurbs(deformation);
    return curve_box+def_box;
}

Foam::BoundingBox Foam::LineStructure::computeBox
(
    label rodNumber,
    scalar parStart,
    scalar parEnd
) const
{
    const ActiveRodMesh::rodCosserat* rod = Rods[rodNumber];

    if(parStart<rod->m_Curve.domainStart())
        FatalErrorInFunction<<"Out of range"<<exit(FatalError);
    if(parEnd<rod->m_Curve.domainEnd())
        FatalErrorInFunction<<"Out of range"<<exit(FatalError);
    if(parEnd<parStart)
        FatalErrorInFunction<<"End smaller than start"<<exit(FatalError);
    
    const gsNurbs<scalar>& curve = rod->m_Curve;
    const gsNurbs<scalar>& def = rod->m_Def;

    return BoundingBox::boundsOfNurbs(curve,parStart,parEnd)+BoundingBox::boundsOfNurbs(def,parStart,parEnd);
}

Foam::scalar Foam::LineStructure::characteristicSize
(
    label rodNumber,
    scalar par
) const
{
    FatalErrorInFunction<<"Not yet implemented!"<<exit(FatalError);
    return 0;
}

void Foam::LineStructure::buildTrees()
{
    rodTrees.resize(nR);
    for(label rodI=0; rodI<nR; rodI++)
        buildTreeOnRod(rodI);
}

void Foam::LineStructure::buildTreeOnRod
(
    label rodNumber
)
{
    const ActiveRodMesh::rodCosserat* rod = Rods[rodNumber];
    
    BoundingBoxTree& rodTree = rodTrees[rodNumber];
    std::unique_ptr<BoundingBoxTree::Node>& root = rodTree.getRoot();
    root = std::make_unique<BoundingBoxTree::Node>();
    
    root->key = 0.5*(rod->m_Curve.domainStart()+rod->m_Curve.domainEnd());
    root->value = computeBox(rodNumber);
    root->leftChild = std::make_unique<BoundingBoxTree::Node>();
    root->rightChild = std::make_unique<BoundingBoxTree::Node>();

    std::function<void(std::unique_ptr<BoundingBoxTree::Node>&,std::pair<scalar,scalar>)> recursiveBuildTree =
    [&](std::unique_ptr<BoundingBoxTree::Node>& node, std::pair<scalar,scalar> bound)
    {
        scalar center = 0.5*(bound.first+bound.second);
        node->key = center;
        node->value = computeBox(rodNumber,bound.first,bound.second);
        scalar boxSize = node->value.innerSize();
        scalar charSize = characteristicSize(rodNumber,center);
        if(charSize<boxSize)
        {
            node->leftChild = std::make_unique<BoundingBoxTree::Node>();
            recursiveBuildTree(node->leftChild,{bound.first,center});
            node->rightChild = std::make_unique<BoundingBoxTree::Node>();
            recursiveBuildTree(node->rightChild,{center,bound.second});
        }        
    };
}

void Foam::LineStructure::setupRodBoundingBoxTree()
{
    buildTrees();
    rodInMesh.resize(nR);
    for(label rodI=0; rodI<nR; rodI++)
    {
        if(meshBoundingBox.boundingBoxOverlap(rodTrees[rodI].rootBox()))
            rodInMesh[rodI] = true;
        else
            rodInMesh[rodI] = false;
    }
}

void Foam::LineStructure::setDeformation
(
    const List<List<vector>>& deformationCoeffs
)
{
    Structure::setDeformation(deformationCoeffs);
    setupRodBoundingBoxTree();
}

Foam::vector Foam::LineStructure::evaluateRodCircumPos
(
    label rodNumber,
    scalar parameter,
    scalar angle,
    scalar radiusFrac,
    scalar var_para,
    scalar var_radius
)
{
    FatalErrorInFunction<<"Invalid call"<<exit(FatalError);
    return vector(0,0,0);
}
        
Foam::Pair<Foam::vector> Foam::LineStructure::derivateRodCircumPos
(
    label rodNumber,
    scalar parameter,
    scalar angle,
    scalar radiusFrac
)
{
    FatalErrorInFunction<<"Invalid call"<<exit(FatalError);
    return Pair<vector>(vector(0,0,0),vector(0,0,0));
}
    
Foam::Pair<Foam::vector> Foam::LineStructure::derivate2RodCircumPos
(
    label rodNumber,
    scalar parameter,
    scalar angle,
    scalar radiusFrac
)
{
    FatalErrorInFunction<<"Invalid call"<<exit(FatalError);
    return Pair<vector>(vector(0,0,0),vector(0,0,0));
}

Foam::LineStructure::GlobalHaloMarkers::GlobalHaloMarkers
(
    const LineStructure* structure,
    // [haloCells] -> [marker of cell] -> ptr,index
    std::vector<std::vector<std::pair<LagrangianMarker*,label>>>& selfhaloCellsRodMarkersList
):
structure(structure),
selfhaloCellsRodMarkersList(&selfhaloCellsRodMarkersList),
nProcs(Pstream::nProcs()),
broadcasted(false),
broadcastedWeights(false)
{
    label selfHaloCellSize = this->selfhaloCellsRodMarkersList->size();
    
    globalHaloCellsMarkerPos.clear();
    globalHaloCellsMarkerPos.setSize(Pstream::nProcs());
    globalHaloCellsMarkerPos[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerPos[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsMarkerVolume.clear();
    globalHaloCellsMarkerVolume.setSize(Pstream::nProcs());
    globalHaloCellsMarkerVolume[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerVolume[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsLocalIndex.clear();
    globalHaloCellsLocalIndex.setSize(Pstream::nProcs());
    globalHaloCellsLocalIndex[Pstream::myProcNo()].clear();
    globalHaloCellsLocalIndex[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsMarkerDilation.clear();
    globalHaloCellsMarkerDilation.setSize(Pstream::nProcs());
    globalHaloCellsMarkerDilation[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerDilation[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsMarkerSupportCellIndices.clear();
    globalHaloCellsMarkerSupportCellIndices.setSize(Pstream::nProcs());
    globalHaloCellsMarkerSupportCellIndices[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerSupportCellIndices[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsMarkerSupportCellCentres.clear();
    globalHaloCellsMarkerSupportCellCentres.setSize(Pstream::nProcs());
    globalHaloCellsMarkerSupportCellCentres[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerSupportCellCentres[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsMarkerSupportCellVolume.clear();
    globalHaloCellsMarkerSupportCellVolume.setSize(Pstream::nProcs());
    globalHaloCellsMarkerSupportCellVolume[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerSupportCellVolume[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsMarkerb.clear();
    globalHaloCellsMarkerb.setSize(Pstream::nProcs());
    globalHaloCellsMarkerb[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerb[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsMarkerWeight.clear();
    globalHaloCellsMarkerWeight.setSize(Pstream::nProcs());
    globalHaloCellsMarkerWeight[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerWeight[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsMarkerNurbsInd.clear();
    globalHaloCellsMarkerNurbsInd.setSize(Pstream::nProcs());
    globalHaloCellsMarkerNurbsInd[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerNurbsInd[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsMarkerParameter.clear();
    globalHaloCellsMarkerParameter.setSize(Pstream::nProcs());
    globalHaloCellsMarkerParameter[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerParameter[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsMarkerAngle.clear();
    globalHaloCellsMarkerAngle.setSize(Pstream::nProcs());
    globalHaloCellsMarkerAngle[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerAngle[Pstream::myProcNo()].setSize(selfHaloCellSize);
    
    globalHaloCellsMarkerRadiusFrac.clear();
    globalHaloCellsMarkerRadiusFrac.setSize(Pstream::nProcs());
    globalHaloCellsMarkerRadiusFrac[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerRadiusFrac[Pstream::myProcNo()].setSize(selfHaloCellSize);

    procHaloCellsSize.setSize(nProcs);
    procHaloCellsSize[Pstream::myProcNo()] = selfHaloCellSize;
    structure->exchangeBetweenAll(procHaloCellsSize);

    procHaloCellMarkerSize.setSize(nProcs);
    const std::unordered_set<label>& neighbours = structure->getNeighbourProcesses();
    for(label proc=0; proc<nProcs; proc++)
    {
        if(neighbours.find(proc)!=neighbours.end() || proc==Pstream::myProcNo())
            procHaloCellMarkerSize[proc].setSize(procHaloCellsSize[proc],0);
    }
}

void Foam::LineStructure::GlobalHaloMarkers::appendMarkerData
(
    label haloCellInd,
    std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,FixedList<scalar,10>,label,scalar,scalar,scalar> data
)
{
    if(structure==nullptr)
        FatalErrorInFunction<<"Structure nullptr"<<exit(FatalError);
    if(broadcasted)
        FatalErrorInFunction<<"Already broadcasted"<<exit(FatalError);
    
    label proc = Pstream::myProcNo();
    if(haloCellInd>=procHaloCellsSize[proc])
        FatalErrorInFunction<<"Out of range value for halocCellInd"<<exit(FatalError);
    label hCI = haloCellInd;
    
    globalHaloCellsMarkerPos[proc][hCI].append(std::get<0>(data));
    globalHaloCellsMarkerVolume[proc][hCI].append(std::get<1>(data));
    globalHaloCellsLocalIndex[proc][hCI].append(std::get<2>(data));
    globalHaloCellsMarkerDilation[proc][hCI].append(std::get<3>(data));
    globalHaloCellsMarkerSupportCellIndices[proc][hCI].append(std::get<4>(data));
    globalHaloCellsMarkerSupportCellCentres[proc][hCI].append(std::get<5>(data));
    globalHaloCellsMarkerSupportCellVolume[proc][hCI].append(std::get<6>(data));
    globalHaloCellsMarkerb[proc][hCI].append(std::get<7>(data));
    globalHaloCellsMarkerWeight[proc][hCI].append(-1);
    globalHaloCellsMarkerNurbsInd[proc][hCI].append(std::get<8>(data));
    globalHaloCellsMarkerParameter[proc][hCI].append(std::get<9>(data));
    globalHaloCellsMarkerAngle[proc][hCI].append(std::get<10>(data));
    globalHaloCellsMarkerRadiusFrac[proc][hCI].append(std::get<11>(data));
    
    procHaloCellMarkerSize[proc][hCI]++;
}

void Foam::LineStructure::GlobalHaloMarkers::insertMarkerWeight
(
    label haloCellInd,
    label cellMarkerInd,
    scalar weight
)
{
    if(structure==nullptr)
        FatalErrorInFunction<<"Structure nullptr"<<exit(FatalError);
    if(broadcastedWeights)
        FatalErrorInFunction<<"Already broadcasted weights"<<exit(FatalError);
    
    label proc = Pstream::myProcNo();
    if(haloCellInd>=procHaloCellsSize[proc])
        FatalErrorInFunction<<"Out of range value for halocCellInd"<<exit(FatalError);
    label hCI = haloCellInd;
    if(cellMarkerInd>=procHaloCellMarkerSize[proc][hCI])
        FatalErrorInFunction<<"Out of range value for marker ind"<<exit(FatalError);
    label cMI = cellMarkerInd;
    
    globalHaloCellsMarkerWeight[proc][hCI][cMI] = weight;
}

std::tuple<Foam::vector,Foam::scalar,Foam::label,Foam::vector,Foam::DynamicList<Foam::Pair<Foam::label>>,Foam::DynamicList<Foam::vector>,Foam::DynamicList<Foam::scalar>,Foam::FixedList<Foam::scalar,10>,Foam::label,Foam::scalar,Foam::scalar,Foam::scalar> Foam::LineStructure::GlobalHaloMarkers::getMarkerData
(
    label process,
    label haloCellInd,
    label cellMarkerInd
) const
{
    if(structure==nullptr)
        FatalErrorInFunction<<"Structure nullptr"<<exit(FatalError);
    if(!broadcasted)
        FatalErrorInFunction<<"Not broadcasted yet"<<Foam::nl;
    
    if(!(structure->neighbourDataExists(process)))
        FatalErrorInFunction<<"Tried to access non neighbour"<<exit(FatalError);

    label proc = process;
    label hCI = haloCellInd;
    label cMI = cellMarkerInd;

    check(proc,hCI,cMI);
    
    vector position = globalHaloCellsMarkerPos[proc][hCI][cMI];
    scalar volume = globalHaloCellsMarkerVolume[proc][hCI][cMI];
    label index = globalHaloCellsLocalIndex[proc][hCI][cMI];
    vector dilation = globalHaloCellsMarkerDilation[proc][hCI][cMI];
    DynamicList<Pair<label>> suppCellsIndices = globalHaloCellsMarkerSupportCellIndices[proc][hCI][cMI];
    DynamicList<vector> suppCellsCentre = globalHaloCellsMarkerSupportCellCentres[proc][hCI][cMI];
    DynamicList<scalar> suppCellsVolume = globalHaloCellsMarkerSupportCellVolume[proc][hCI][cMI];
    FixedList<scalar,10> b = globalHaloCellsMarkerb[proc][hCI][cMI];
    label nurbsInd = globalHaloCellsMarkerNurbsInd[proc][hCI][cMI];
    scalar parameter = globalHaloCellsMarkerParameter[proc][hCI][cMI];
    scalar angle = globalHaloCellsMarkerAngle[proc][hCI][cMI];
    scalar radiusFrac = globalHaloCellsMarkerRadiusFrac[proc][hCI][cMI];
    
    return {position,volume,index,dilation,suppCellsIndices,suppCellsCentre,suppCellsVolume,b,nurbsInd,parameter,angle,radiusFrac};
}

Foam::scalar Foam::LineStructure::GlobalHaloMarkers::getMarkerWeight
(
    label process,
    label haloCellInd,
    label cellMarkerInd
) const
{
    if(structure==nullptr)
        FatalErrorInFunction<<"Structure nullptr"<<exit(FatalError);
    if(!broadcastedWeights)
        FatalErrorInFunction<<"Not broadcasted yet"<<exit(FatalError);
    
    if(!(structure->neighbourDataExists(process)))
        FatalErrorInFunction<<"Tried to access non neighbour"<<exit(FatalError);
    
    if(haloCellInd>=procHaloCellsSize[process])
        FatalErrorInFunction<<"Out of range value for halocCellInd"<<exit(FatalError);
    label hCI = haloCellInd;
    if(cellMarkerInd>=procHaloCellMarkerSize[process][hCI])
    {
        Pout<<"process:"<<process<<Foam::nl;
        Pout<<"haloCellInd:"<<haloCellInd<<Foam::nl;
        Pout<<"cellMarkerInd:"<<cellMarkerInd<<Foam::nl;
        Pout<<"procHaloCellsSize["<<process<<"]:"<<procHaloCellsSize[process]<<Foam::nl;
        Pout<<"procHaloCellMarkerSize["<<process<<"]["<<hCI<<"]:"<<procHaloCellMarkerSize[process][hCI]<<Foam::nl;
        FatalErrorInFunction<<"Out of range value for marker ind"<<exit(FatalError);
    }
    label cMI = cellMarkerInd;
    
    return globalHaloCellsMarkerWeight[process][hCI][cMI];
}

Foam::label Foam::LineStructure::GlobalHaloMarkers::size_processes() const
{
    if(structure==nullptr)
        FatalErrorInFunction<<"Structure nullptr"<<exit(FatalError);
    return Pstream::nProcs();
}

Foam::label Foam::LineStructure::GlobalHaloMarkers::size_haloCells
(
    label process
) const
{
    if(structure==nullptr)
        FatalErrorInFunction<<"Structure nullptr"<<exit(FatalError);
    return procHaloCellsSize[process];
}

Foam::label Foam::LineStructure::GlobalHaloMarkers::size_cellMarkers
(
    label process,
    label haloCellInd
) const
{
    if(structure==nullptr)
        FatalErrorInFunction<<"Structure nullptr"<<exit(FatalError);
    if(!(structure->neighbourDataExists(process)))
        FatalErrorInFunction<<"Tried to access non neighbour"<<exit(FatalError);
    return procHaloCellMarkerSize[process][haloCellInd];
}

void Foam::LineStructure::GlobalHaloMarkers::communicate()
{
    if(structure==nullptr)
        FatalErrorInFunction<<"Structure nullptr"<<exit(FatalError);
    
    broadcasted = true;
    
    structure->exchangeBetweenAll(globalHaloCellsMarkerPos);
    structure->exchangeBetweenAll(globalHaloCellsMarkerVolume);
    structure->exchangeBetweenAll(globalHaloCellsLocalIndex);
    structure->exchangeBetweenAll(globalHaloCellsMarkerDilation);
    structure->exchangeBetweenAll(globalHaloCellsMarkerSupportCellIndices);
    structure->exchangeBetweenAll(globalHaloCellsMarkerSupportCellCentres);
    structure->exchangeBetweenAll(globalHaloCellsMarkerSupportCellVolume);
    structure->exchangeBetweenAll(globalHaloCellsMarkerb);
    structure->exchangeBetweenAll(globalHaloCellsMarkerWeight);
    structure->exchangeBetweenAll(globalHaloCellsMarkerNurbsInd);
    structure->exchangeBetweenAll(globalHaloCellsMarkerParameter);
    structure->exchangeBetweenAll(globalHaloCellsMarkerAngle);
    structure->exchangeBetweenAll(globalHaloCellsMarkerRadiusFrac);
       
    structure->exchangeBetweenAll(procHaloCellMarkerSize);

    //check();
}

void Foam::LineStructure::GlobalHaloMarkers::communicateWeight()
{
    if(structure==nullptr)
        FatalErrorInFunction<<"Structure nullptr"<<exit(FatalError);
        
    broadcastedWeights = true;
    
    structure->exchangeBetweenAll(globalHaloCellsMarkerWeight);
}

void Foam::LineStructure::printMarkerStructure()
{
    for(std::size_t rodNumber=0; rodNumber<rodMarkersList.size(); rodNumber++)
    {
        Pout<<"--------------------rodNumber:"<<rodNumber<<"----------------------"<<Foam::endl;
        std::list<LagrangianMarker>& oneRod = *(rodMarkersList[rodNumber]);
        
        for(auto iterTang = oneRod.begin(); iterTang!=oneRod.end(); iterTang++)
        {
            Pout<<"tang:"<<iterTang->getMarkerPosition()<<Foam::endl;
        }
    }
}

void Foam::LineStructure::writeCellMarkerCountField
(
    volScalarField& field
) const
{
    field = Foam::zero();
    for(const LagrangianMarker* markerPtr : collectedMarkers)
    {
        label markerCell = markerPtr->getMarkerCell();
        if(markerCell>=0 && markerCell<mesh.cells().size())
            field[markerCell] += 1;
    }
}

void Foam::LineStructure::writeCellMarkerCharacSizeField
(
    volScalarField& field
) const
{
    field = Foam::zero();
    std::unordered_map<label,DynamicList<scalar>> cellToCharacLen;
    for(const LagrangianMarker* markerPtr : collectedMarkers)
    {
        label markerCell = markerPtr->getMarkerCell();
        scalar markerCharacLen = markerPtr->getMarkerCharacLen();
        if(markerCell>=0 && markerCell<mesh.cells().size())
        {
            cellToCharacLen[markerCell].append(markerCharacLen);
        }
    }
    for(auto iter=cellToCharacLen.begin(); iter!=cellToCharacLen.end(); iter++)
    {
        const DynamicList<scalar>& characLens = iter->second;
        label markerCell = iter->first;
        scalar avgLens = 0;
        for(scalar len : characLens)
            avgLens += len;
        avgLens /= characLens.size();
        field[markerCell] = avgLens;
    }
}

void Foam::LineStructure::GlobalHaloMarkers::check()
{
    label numProcs = size_processes();
    if(numProcs!=globalHaloCellsMarkerPos.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerPos size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsMarkerVolume.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerVolume size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsLocalIndex.size())
        FatalErrorInFunction<<"globalHaloCellsLocalIndex size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsMarkerDilation.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerDilation size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsMarkerSupportCellIndices.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellIndices size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsMarkerSupportCellCentres.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellCentres size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsMarkerSupportCellVolume.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellVolume size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsMarkerb.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerb size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsMarkerWeight.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerWeight size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsMarkerNurbsInd.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerNurbsInd size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsMarkerParameter.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerParameter size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsMarkerAngle.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerAngle size error"<<exit(FatalError);
    if(numProcs!=globalHaloCellsMarkerRadiusFrac.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerRadiusFrac size error"<<exit(FatalError);
    
    for(label proc=0; proc<numProcs; proc++)
    {
        label numHaloCells = size_haloCells(proc);
        if(numHaloCells!=globalHaloCellsMarkerPos[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerPos["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsMarkerVolume[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerVolume["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsLocalIndex[proc].size())
            FatalErrorInFunction<<"globalHaloCellsLocalIndex["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsMarkerDilation[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerDilation["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsMarkerSupportCellIndices[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellIndices["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsMarkerSupportCellCentres[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellCentres["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsMarkerSupportCellVolume[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellVolume["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsMarkerb[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerb["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsMarkerWeight[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerWeight["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsMarkerNurbsInd[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerNurbsInd["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsMarkerParameter[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerParameter["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsMarkerAngle[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerAngle["<<proc<<"] size error"<<exit(FatalError);
        if(numHaloCells!=globalHaloCellsMarkerRadiusFrac[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerRadiusFrac["<<proc<<"] size error"<<exit(FatalError);
        
        for(label haloCellInd=0; haloCellInd<numHaloCells; haloCellInd++)
        {
            label numHaloCellMarkers = size_cellMarkers(proc,haloCellInd);
            if(numHaloCellMarkers!=globalHaloCellsMarkerPos[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerPos["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsMarkerPos[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsMarkerVolume[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerVolume["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsMarkerVolume[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsLocalIndex[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsLocalIndex["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsLocalIndex[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsMarkerDilation[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerDilation["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsMarkerDilation[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsMarkerSupportCellIndices[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellIndices["<<proc<<" size "<<globalHaloCellsMarkerSupportCellIndices[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsMarkerSupportCellCentres[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellCentres["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsMarkerSupportCellCentres[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsMarkerSupportCellVolume[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellVolume["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsMarkerSupportCellVolume[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsMarkerb[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerb["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsMarkerb[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsMarkerWeight[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerWeight["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsMarkerWeight[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsMarkerNurbsInd[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerNurbsInd["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsMarkerNurbsInd[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsMarkerParameter[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerParameter["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsMarkerParameter[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsMarkerAngle[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerAngle["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsMarkerAngle[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
            if(numHaloCellMarkers!=globalHaloCellsMarkerRadiusFrac[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerRadiusFrac["<<proc<<"]["<<haloCellInd<<"] size "<<globalHaloCellsMarkerRadiusFrac[proc][haloCellInd].size()<<" error:"<<numHaloCellMarkers<<exit(FatalError);
        }
    }    
}

void Foam::LineStructure::GlobalHaloMarkers::check
(
    label process,
    label haloCellInd,
    label cellMarkerInd
) const
{
    if(process<0)
        FatalErrorInFunction<<"Invalid process index"<<Foam::nl;
    if(process>=globalHaloCellsMarkerPos.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerPos size error"<<exit(FatalError);
    if(process>=globalHaloCellsMarkerVolume.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerVolume size error"<<exit(FatalError);
    if(process>=globalHaloCellsMarkerDilation.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerDilation size error"<<exit(FatalError);
    if(process>=globalHaloCellsMarkerSupportCellIndices.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellIndices size error"<<exit(FatalError);
    if(process>=globalHaloCellsMarkerSupportCellCentres.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellCentres size error"<<exit(FatalError);
    if(process>=globalHaloCellsMarkerSupportCellVolume.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellVolume size error"<<exit(FatalError);
    if(process>=globalHaloCellsMarkerb.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerb size error"<<exit(FatalError);
    if(process>=globalHaloCellsMarkerWeight.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerWeight size error"<<exit(FatalError);
    if(process>=globalHaloCellsMarkerNurbsInd.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerNurbsInd size error"<<exit(FatalError);
    if(process>=globalHaloCellsMarkerParameter.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerParameter size error"<<exit(FatalError);
    if(process>=globalHaloCellsMarkerAngle.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerAngle size error"<<exit(FatalError);
    if(process>=globalHaloCellsMarkerRadiusFrac.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerRadiusFrac size error"<<exit(FatalError);
    
    if(haloCellInd<0)
        FatalErrorInFunction<<"Invalid haloCell index"<<Foam::nl;
    if(haloCellInd>=globalHaloCellsMarkerPos[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerPos[process] size error"<<exit(FatalError);
    if(haloCellInd>=globalHaloCellsMarkerVolume[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerVolume[process] size error"<<exit(FatalError);
    if(haloCellInd>=globalHaloCellsMarkerDilation[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerDilation[process] size error"<<exit(FatalError);
    if(haloCellInd>=globalHaloCellsMarkerSupportCellIndices[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellIndices[process] size error"<<exit(FatalError);
    if(haloCellInd>=globalHaloCellsMarkerSupportCellCentres[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellCentres[process] size error"<<exit(FatalError);
    if(haloCellInd>=globalHaloCellsMarkerSupportCellVolume[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellVolume[process] size error"<<exit(FatalError);
    if(haloCellInd>=globalHaloCellsMarkerb[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerb[process] size error"<<exit(FatalError);
    if(haloCellInd>=globalHaloCellsMarkerWeight[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerWeight[process] size error"<<exit(FatalError);
    if(haloCellInd>=globalHaloCellsMarkerNurbsInd[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerNurbsInd[process] size error"<<exit(FatalError);
    if(haloCellInd>=globalHaloCellsMarkerParameter[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerParameter[process] size error"<<exit(FatalError);
    if(haloCellInd>=globalHaloCellsMarkerAngle[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerAngle[process] size error"<<exit(FatalError);
    if(haloCellInd>=globalHaloCellsMarkerRadiusFrac[process].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerRadiusFrac[process] size error"<<exit(FatalError);
    
    if(cellMarkerInd<0)
        FatalErrorInFunction<<"Invalid cellMarker index"<<Foam::nl;
    if(cellMarkerInd>=globalHaloCellsMarkerPos[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerPos[process][haloCellInd] size error"<<exit(FatalError);
    if(cellMarkerInd>=globalHaloCellsMarkerVolume[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerVolume[process][haloCellInd] size error"<<exit(FatalError);
    if(cellMarkerInd>=globalHaloCellsMarkerDilation[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerDilation[process][haloCellInd] size error"<<exit(FatalError);
    if(cellMarkerInd>=globalHaloCellsMarkerSupportCellIndices[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellIndices[process][haloCellInd] size error"<<exit(FatalError);
    if(cellMarkerInd>=globalHaloCellsMarkerSupportCellCentres[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellCentres[process][haloCellInd] size error"<<exit(FatalError);
    if(cellMarkerInd>=globalHaloCellsMarkerSupportCellVolume[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellVolume[process][haloCellInd] size error"<<exit(FatalError);
    if(cellMarkerInd>=globalHaloCellsMarkerb[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerb[process][haloCellInd] size error"<<exit(FatalError);
    if(cellMarkerInd>=globalHaloCellsMarkerWeight[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerWeight[process][haloCellInd] size error"<<exit(FatalError);
    if(cellMarkerInd>=globalHaloCellsMarkerNurbsInd[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerNurbsInd[process][haloCellInd] size error"<<exit(FatalError);
    if(cellMarkerInd>=globalHaloCellsMarkerParameter[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerParameter[process][haloCellInd] size error"<<exit(FatalError);
    if(cellMarkerInd>=globalHaloCellsMarkerAngle[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerAngle[process][haloCellInd] size error"<<exit(FatalError);
    if(cellMarkerInd>=globalHaloCellsMarkerRadiusFrac[process][haloCellInd].size())
        FatalErrorInFunction<<"globalHaloCellsMarkerRadiusFrac[process][haloCellInd] size error"<<exit(FatalError);
}

void Foam::LineStructure::parameterGradientCheck()
{
    Structure::parameterGradientCheck();
    Pout<<"LineStructure::parameterGradientCheck"<<Foam::nl;
    
    scalar nbrSteps = 20;
    scalar epsilon = 1e-4;
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        scalar domainStart = this->domainStart(rodNumber)+5*epsilon;
        scalar domainEnd = this->domainEnd(rodNumber)-5*epsilon;
        scalar delta = domainEnd-domainStart;
        scalar stepsize = delta/nbrSteps;
        for(scalar parameter=domainStart; parameter<=domainEnd; parameter+=stepsize)
        {
            //Pout<<"parameter:"<<parameter<<Foam::nl;
            
            // Compute gradients
            vector drdp = derivateRodPos(Rods[rodNumber],parameter);
            vector d2rdp = derivate2RodPos(Rods[rodNumber],parameter);
            
            // Lower value
            scalar lower_parameter = parameter-epsilon;
            vector lower_r = evaluateRodPos(Rods[rodNumber],lower_parameter);
            vector lower_drdp = derivateRodPos(Rods[rodNumber],lower_parameter);
            
            // Upper value
            scalar upper_parameter = parameter+epsilon;
            vector upper_r = evaluateRodPos(Rods[rodNumber],upper_parameter);
            vector upper_drdp = derivateRodPos(Rods[rodNumber],upper_parameter);
            
            vector fd_drdp = (upper_r-lower_r)/(2*epsilon);
            vector fd_d2rdp = (upper_drdp-lower_drdp)/(2*epsilon);
                        
            if(vectorDistance(fd_drdp,drdp)>epsilon)
            {
                Pout<<"parameter:"<<parameter<<Foam::nl;
                Pout<<"drdp:"<<drdp<<Foam::nl;
                Pout<<"lower_r:"<<lower_r<<Foam::nl;
                Pout<<"upper_r:"<<upper_r<<Foam::nl;
                Pout<<"fd_drdp:"<<fd_drdp<<Foam::nl;
                FatalErrorInFunction<<"Error"<<exit(FatalError);
            }
            if(vectorDistance(fd_d2rdp,d2rdp)>epsilon)
            {
                Pout<<"parameter:"<<parameter;//<<Foam::nl;
                Pout<<" d2rdp:"<<d2rdp;//<<Foam::nl;
                //Pout<<"lower_drdp:"<<lower_drdp<<Foam::nl;
                //Pout<<"upper_drdp:"<<upper_drdp<<Foam::nl;
                Pout<<" fd_d2rdp:"<<fd_d2rdp<<Foam::nl;
                //FatalErrorInFunction<<"Error"<<exit(FatalError);
            }      
        }
    }
}

void Foam::LineStructure::selfCheck()
{
    Structure::selfCheck();
}
