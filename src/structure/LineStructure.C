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
crossSecArea(crossSecArea)
{
    Info<<"modusFieldToMarker:"<<modusFieldToMarker<<Foam::endl;
    Info<<"modusMarkerToField:"<<modusMarkerToField<<Foam::endl;
    
    label nbrOfRods = myMesh->m_nR;
    rodMarkersList.resize(nbrOfRods);
    check();
    initialRodPoints.resize(nbrOfRods);
    createSpacingPoints();
    
    for(const auto& oneRodPoints : initialRodPoints)
        for(scalar pnt : *oneRodPoints)
            Info<<pnt<<Foam::endl;
    
    /*
    for(const auto& oneRodPoints : rodMarkersList)
        for(auto pnt : *oneRodPoints)
            Info<<pnt.to_string()<<Foam::endl;
    */
    //rodMarkersList.resize(nbrOfRods);
    
    const DynamicList<CellDescription>& selfHaloCellList = getHaloCellList(Pstream::myProcNo());
    
    //haloCellsRodMarkersList.resize(selfHaloCellList.size());
    
    /*
    globalHaloCellsMarkerPos.resize(Pstream::nProcs());
    globalHaloCellsMarkerPos[Pstream::myProcNo()].resize(selfHaloCellList.size());
    globalHaloCellsMarkerVolume.resize(Pstream::nProcs());
    globalHaloCellsMarkerVolume[Pstream::myProcNo()].resize(selfHaloCellList.size());
    globalHaloCellsMarkerDilation.resize(Pstream::nProcs());
    globalHaloCellsMarkerDilation[Pstream::myProcNo()].resize(selfHaloCellList.size());
    globalHaloCellsMarkerSupportCellCentres.resize(Pstream::nProcs());
    globalHaloCellsMarkerSupportCellCentres[Pstream::myProcNo()].resize(selfHaloCellList.size());
    globalHaloCellsMarkerSupportCellVolume.resize(Pstream::nProcs());
    globalHaloCellsMarkerSupportCellVolume[Pstream::myProcNo()].resize(selfHaloCellList.size());
    */
    
    createMarkersFromSpacedPoints();
    refineMarkers();
    setMarkerVolume();
    for(LagrangianMarker& marker : *(rodMarkersList[0]))
        Info<<marker.to_string()<<Foam::endl;
    
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[0];
    Info<<"dist:0.375->0.40625:"<<LineStructure::distance(oneRod,0.375,0.40625)<<Foam::endl;
    evaluateMarkerMeshRelation();
    reduceMarkers();
    collectMarkers();
    collectHaloMarkers();
    exchangeHaloMarkersData();

    computeMarkerCellWeights();    
    for(LagrangianMarker& marker : *(rodMarkersList[0]))
        Info<<marker.to_string()<<Foam::endl;
    
    //std::unique_ptr<LinearSystem> system = computeMarkerEpsilonMatrix();
    computeMarkerWeights();
    
    /*
    if(Pstream::myProcNo()==1)
    {
        //std::cout<<std::get<1>(*system)<<std::endl;
        //std::cout<<"neighborEps:"<<std::get<2>(*system).size()<<std::endl;
        
    }
    */
}

void Foam::LineStructure::check()
{
    if(!myMesh)
        FatalErrorInFunction<<"Rod Mesh not set!"<<exit(FatalError);
    if(myMesh->m_nR!=static_cast<int>(rodMarkersList.size()))
        FatalErrorInFunction<<"Mismatch in size of m_Rods and rodMarkersList"<<exit(FatalError);
    if(crossSecArea.size()!=static_cast<int>(rodMarkersList.size()))
        FatalErrorInFunction<<"Mismatch in size of crossSecArea and rodMarkersList"<<exit(FatalError);
    if(static_cast<int>(myMesh->m_Rods.size())!=myMesh->m_nR)
        FatalErrorInFunction<<"Mismatch in size of m_Rods and m_nR"<<exit(FatalError);
}

void Foam::LineStructure::connect
(
    FieldMarkerStructureInteraction& connector
)
{
    connector.markers.resize(rodMarkersList.size());
    connectedInteractions.push_back(&connector);
}

/*
void Foam::LineStructure::initializeMarkers()
{
    reInitializeMarkers(false,false);
}
*/

/*
void Foam::LineStructure::reInitializeMarkers
(
    bool keepMarkers,
    bool keepSeedPoints
)
{
    for(uint rodIndex=0; rodIndex<rodMarkers.size(); rodIndex++)
    {
        std::unique_ptr<std::vector<LagrangianMarker>>& oneRodMarkers = rodMarkers[rodIndex];
        if(!oneRodMarkers || !keepMarkers)
        {
            myMesh->m_Rods[rodIndex]->m_Rot.setCoefs(myMesh->m_Rods[rodIndex]->m_Rot_coefs.transpose());
            oneRodMarkers = constructMarkerSet
            (
                rodIndex,myMesh->m_Rods[rodIndex],crossSecArea[rodIndex],initialSpacing,keepSeedPoints
            );
        }
        onnector.markers[rodIndex].resize(0);
        for(LagrangianMarker* oneMarkerPtr : *oneRodMarkers)
        {
            onnector.markers[rodIndex].push_back(oneMarkerPtr);
        }
    }
}
*/

void Foam::LineStructure::createSpacingPoints()
{
    for(uint rodIndex=0; rodIndex<rodMarkersList.size(); rodIndex++)
    {
        createSpacedPointsOnRod(rodIndex,initialMeshSpacing);
    }
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
            if(dist>spacing)
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
    for(uint rodIndex=0; rodIndex<rodMarkersList.size(); rodIndex++)
    {
        createMarkersFromSpacedPointsOnRod(rodIndex);
    }
}

void Foam::LineStructure::createMarkersFromSpacedPointsOnRod
(
    label rodNumber
)
{
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    auto markersPtr = std::unique_ptr<std::list<LagrangianMarker>>(new std::list<LagrangianMarker>());
    std::list<LagrangianMarker>& markers = *markersPtr;
    markers.clear();
    for(scalar point : *(initialRodPoints[rodNumber]))
        markers.push_back(LagrangianMarker(*this,mesh,rodNumber,oneRod,point));
    rodMarkersList[rodNumber] = std::move(markersPtr);
}

void Foam::LineStructure::refineMarkers
(
    std::pair<bool,scalar> forcedSpacing
)
{
    for(uint rodIndex=0; rodIndex<rodMarkersList.size(); rodIndex++)
    {
        refineMarkersOnRod(rodIndex,forcedSpacing);
    }
}

void Foam::LineStructure::refineMarkersOnRod
(
    label rodNumber,
    std::pair<bool,scalar> forcedSpacing
)
{
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
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
            bool markers0InCell = (markersIter0->getMarkerCell()!=-1);
            
            scalar markers1Para = markersIter0->getMarkerParameter();
            scalar markers1CellSpacing = markersIter0->getMarkerCellMinSpacing();
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
                if(dist>minSpacing)
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

void Foam::LineStructure::setMarkerVolume()
{
   for(uint rodIndex=0; rodIndex<rodMarkersList.size(); rodIndex++)
   {
       setMarkerVolumeOnRod(rodIndex);
   }
}

void Foam::LineStructure::setMarkerVolumeOnRod
(
    label rodNumber
)
{
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    std::list<LagrangianMarker>& markers = *(rodMarkersList[rodNumber]);

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
        
        scalar span = LineStructure::distance(oneRod,spanStart,spanEnd);
        Info<<"dist:"<<spanStart<<"->"<<spanEnd<<":"<<LineStructure::distance(oneRod,spanStart,spanEnd)<<Foam::endl;
        Info<<"span:"<<span<<Foam::endl;
        if(modusMarkerToField==markerMeshType::NonUniform)
        {
            Info<<"Set NonUniform"<<Foam::endl;
            iter->setMarkerVolume(span);
        }
        else
            iter->setMarkerVolume(span*crossSecArea[rodNumber]);
        Info<<"getMarkerVolume:"<<iter->getMarkerVolume()<<Foam::endl;
        
        iterPrev = iter;
    }
}

void Foam::LineStructure::evaluateMarkerMeshRelation()
{
    for(std::unique_ptr<std::list<LagrangianMarker>>& singleRodMarkers : rodMarkersList)
        evaluateMarkerMeshRelation(*singleRodMarkers);
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
    std::vector<MarkerReference<LagrangianMarker>> allMarkers;
    for(std::unique_ptr<std::list<LagrangianMarker>>& singleRodMarkersPtr : rodMarkersList)
    {
        std::list<LagrangianMarker>* singleRodMarkers = &(*singleRodMarkersPtr);
        for(auto iter=singleRodMarkers->begin(); iter!=singleRodMarkers->end(); iter++)
        {
            allMarkers.push_back(MarkerReference<LagrangianMarker>(iter,singleRodMarkers));
        }
    }
    reduceMarkers(allMarkers);
}

template<class MARKER>
void Foam::LineStructure::maxOneMarkerPerOctant
(
    std::vector<MarkerReference<MARKER>>& allMarkers,
    std::unordered_map<label,std::vector<MarkerReference<MARKER>*>>& cellToMarker
)
{
    using SingleMarker = MarkerReference<MARKER>;
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    for(auto iter=cellToMarker.begin(); iter!=cellToMarker.end(); iter++)
    {
        label cellInd = iter->first;
        if(cellInd>=0 && cellInd<mesh.cells().size())
        {
            const cell& oneCell = cells[cellInd];
            vector cellCentre = oneCell.centre(points,faces);
            std::pair<
                std::pair<
                    std::pair<
                        DynamicList<SingleMarker*>,
                        DynamicList<SingleMarker*>
                    >,
                    std::pair<
                        DynamicList<SingleMarker*>,
                        DynamicList<SingleMarker*>
                    >
                >,
                std::pair<
                    std::pair<
                        DynamicList<SingleMarker*>,
                        DynamicList<SingleMarker*>
                    >,
                    std::pair<
                        DynamicList<SingleMarker*>,
                        DynamicList<SingleMarker*>
                    >
                >
            > xyz_split;
            for(SingleMarker* markerPtr : iter->second)
            {
                vector markerPos = markerPtr->getMarker().getMarkerPosition();
                std::pair<
                    std::pair<
                        DynamicList<SingleMarker*>,
                        DynamicList<SingleMarker*>
                    >,
                    std::pair<
                        DynamicList<SingleMarker*>,
                        DynamicList<SingleMarker*>
                    >
                >* x_side;
                if(markerPos[0]<cellCentre[0])
                    x_side = &(xyz_split.first);
                else
                    x_side = &(xyz_split.second);
                
                std::pair<
                    DynamicList<SingleMarker*>,
                    DynamicList<SingleMarker*>
                >* xy_side;
                if(markerPos[1]<cellCentre[1])
                    xy_side = &(x_side->first);
                else
                    xy_side = &(x_side->second);
                
                DynamicList<SingleMarker*>* xyz_side;
                if(markerPos[2]<cellCentre[2])
                    xyz_side = &(xy_side->first);
                else
                    xyz_side = &(xy_side->second);
                
                xyz_side->append(markerPtr);
            }

            std::vector<DynamicList<SingleMarker*>*> subCells = 
            {
                &(xyz_split.first.first.first),
                &(xyz_split.first.first.second),
                &(xyz_split.first.second.first),
                &(xyz_split.first.second.second),
                &(xyz_split.second.first.first),
                &(xyz_split.second.first.second),
                &(xyz_split.second.second.first),
                &(xyz_split.second.second.first)
            };
            
            for(DynamicList<SingleMarker*>* subCell : subCells)
            {
                if(subCell->size()>0)
                {
                    vector averagePosition = vector(0,0,0);
                    for(SingleMarker* marker : *subCell)
                        averagePosition += marker->getMarker().getMarkerPosition();
                    averagePosition /= subCell->size();
                    
                    scalar summedVolume = 0;
                    for(SingleMarker* marker : *subCell)
                        summedVolume += marker.getMarkerVolume();
                    
                    SingleMarker* optMarker;
                    scalar optMarkerDistToAvgPos = std::numeric_limits<scalar>::max();
                    for(SingleMarker* marker : *subCell)
                    {
                        vector distVec = marker->getMarker().getMarkerPosition() - averagePosition;
                        scalar dist = std::sqrt(distVec&distVec);
                        if(dist<optMarkerDistToAvgPos)
                        {
                            optMarker = marker;
                            optMarkerDistToAvgPos = dist;
                        }
                    }
                    optMarker->setMarkerVolume(summedVolume);
                    for(SingleMarker* marker : *subCell)
                        if(marker!=optMarker)
                            marker->deleteMarker();
                }
            }
        }
        else
        {
            for(SingleMarker* markerPtr : iter->second)
            {
                markerPtr->deleteMarker();
            }
        }
    }
}

template<class MARKER>
void Foam::LineStructure::maxOneMarkerPerCell
(
    std::vector<MarkerReference<MARKER>>& allMarkers,
    std::unordered_map<label,std::vector<MarkerReference<MARKER>*>>& cellToMarker
)
{
    using SingleMarker = MarkerReference<MARKER>;
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    for(auto iter=cellToMarker.begin(); iter!=cellToMarker.end(); iter++)
    {
        label cellInd = iter->first;
        if(cellInd>=0 && cellInd<mesh.cells().size())
        {
            const cell& oneCell = cells[cellInd];
            vector cellCentre = oneCell.centre(points,faces);
            std::vector<SingleMarker*>& oneCellMarkers = iter->second;
            
            if(oneCellMarkers.size()>0)
            {
                vector averagePosition = vector(0,0,0);
                for(SingleMarker* marker : oneCellMarkers)
                    averagePosition += marker->getMarker().getMarkerPosition();
                averagePosition /= oneCellMarkers.size();
                
                scalar summedVolume = 0;
                for(SingleMarker* marker : oneCellMarkers)
                    summedVolume += marker->getMarker().getMarkerVolume();
                
                SingleMarker* optMarker;
                scalar optMarkerDistToAvgPos = std::numeric_limits<scalar>::max();
                for(SingleMarker* marker : oneCellMarkers)
                {
                    vector distVec = marker->getMarker().getMarkerPosition() - averagePosition;
                    scalar dist = std::sqrt(distVec&distVec);
                    if(dist<optMarkerDistToAvgPos)
                    {
                        optMarker = marker;
                        optMarkerDistToAvgPos = dist;
                    }
                }
                optMarker->getMarker().setMarkerVolume(summedVolume);
                for(SingleMarker* marker : oneCellMarkers)
                    if(marker!=optMarker)
                        marker->deleteMarker();
            }
        }
        else
        {
            for(SingleMarker* markerPtr : iter->second)
            {
                markerPtr->deleteMarker();
            }
        }
    }
}

template<class MARKER>
void Foam::LineStructure::oneCenteredMarker
(
    std::vector<MarkerReference<MARKER>>& allMarkers,
    std::unordered_map<label,std::vector<MarkerReference<MARKER>*>>& cellToMarker
)
{
    using SingleMarker = MarkerReference<MARKER>;
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    for(auto iter=cellToMarker.begin(); iter!=cellToMarker.end(); iter++)
    {
        label cellInd = iter->first;
        if(cellInd>=0 && cellInd<mesh.cells().size())
        {
            const cell& oneCell = cells[cellInd];
            vector cellCentre = oneCell.centre(points,faces);
            std::vector<SingleMarker*>& oneCellMarkers = iter->second;
            
            scalar summedVolume = 0;
            for(SingleMarker* marker : oneCellMarkers)
                summedVolume += marker->getMarker().getMarkerVolume();
            
            if(oneCellMarkers.size()>0)
            {                        
                SingleMarker* optMarker = nullptr;
                scalar optMarkerDistToAvgPos = std::numeric_limits<scalar>::max();
                for(SingleMarker* marker : oneCellMarkers)
                {
                    vector distVec = marker->getMarker().getMarkerPosition() - cellCentre;
                    scalar dist = std::sqrt(distVec&distVec);
                    if(dist<optMarkerDistToAvgPos)
                    {
                        optMarker = marker;
                        optMarkerDistToAvgPos = dist;
                    }
                }
                optMarker->getMarker().setMarkerVolume(summedVolume);
                for(SingleMarker* marker : oneCellMarkers)
                    if(marker!=optMarker)
                        marker->deleteMarker();
            }
        }
        else
        {
            for(SingleMarker* markerPtr : iter->second)
            {
                markerPtr->deleteMarker();
            }
        }
    }
}

void Foam::LineStructure::collectMarkers()
{
    collectedMarkers.resize(0);
    for(std::unique_ptr<std::list<LagrangianMarker>>& oneRodMarkers : rodMarkersList)
    {
        for(LagrangianMarker& oneMarker : *oneRodMarkers)
        {
            collectedMarkers.push_back(&oneMarker);
        }
    }
}

void Foam::LineStructure::collectHaloMarkers()
{
    haloCellsRodMarkersList.clear();
    haloCellsRodMarkersList.resize(getHaloCellList(Pstream::myProcNo()).size());
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
}

void Foam::LineStructure::exchangeHaloMarkersData()
{
    globHaloMarkers = GlobalHaloMarkers(haloCellsRodMarkersList);
    //const cellList& cells = mesh.cells();
    //const faceList& faces = mesh.faces();
    //const pointField& points = mesh.points();
        
    for(uint haloCellInd=0; haloCellInd<haloCellsRodMarkersList.size(); haloCellInd++)
    {
        std::vector<std::pair<LagrangianMarker*,label>>& localHaloCellRodMarkers = haloCellsRodMarkersList[haloCellInd];
        for(uint haloMarkerInd=0; haloMarkerInd<haloCellsRodMarkersList[haloCellInd].size(); haloMarkerInd++)
        {
            std::pair<LagrangianMarker*,label> markerData = localHaloCellRodMarkers[haloMarkerInd];
            LagrangianMarker* marker = markerData.first;
            label markerInd = markerData.second;
            
            vector position = marker->getMarkerPosition();
            scalar volume = marker->getMarkerVolume();
            label index = markerInd;
            vector dilation = marker->getDilation();
            DynamicList<Pair<label>> supportCellsIndices;
            DynamicList<vector> supportCellsCentre;
            DynamicList<scalar> supportCellsVolume;
            
            const DynamicList<std::tuple<bool,label,label>>& supportCells = marker->getSupportCells();
            for(std::tuple<bool,label,label> oneSupport : supportCells)
            {
                vector oneSuppCellCentre;
                scalar oneSuppCellVolume;
                marker->getCellData(oneSupport,oneSuppCellCentre,oneSuppCellVolume);
                supportCellsIndices.append({std::get<1>(oneSupport),std::get<2>(oneSupport)});
                supportCellsCentre.append(oneSuppCellCentre);
                supportCellsVolume.append(oneSuppCellVolume);
            }
            List<scalar> b(10);
            for(label i=0; i<10; i++)
                b[i] = marker->getCorrParaB()[i];
            globHaloMarkers.appendMarkerData
            (
                haloCellInd,
                {position,volume,index,dilation,supportCellsIndices,supportCellsCentre,supportCellsVolume,b}
            );
        }
    }
    globHaloMarkers.communicate();
}

void Foam::LineStructure::computeMarkerCellWeights()
{
    for(LagrangianMarker* marker : collectedMarkers)
        marker->compNonUniformCorrWeights();
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
        const DynamicList<std::tuple<bool,label,label>>& supportI = markerI.getSupportCells();
        CellSet markerISupportMap;
        for(const std::tuple<bool,label,label>& suppCell : supportI)
        {
            vector centre;
            scalar vol;
            markerI.getCellData(suppCell,centre,vol);
            markerISupportMap[centre] = vol;
        }
        
        std::map<label,scalar> rowEntries;
        
        // Compute local matrix entries
        for(uint K=0; K<collectedMarkers.size(); K++)
        {
            const LagrangianMarker& markerK = *(collectedMarkers[K]);
            vector XK = markerK.getMarkerPosition();
            vector dilationK = markerK.getDilation();
            scalar dilationKMax = std::max<scalar>(dilationK[0],dilationK[1]);
            dilationKMax = std::max<scalar>(dilationKMax,dilationK[2]);
            scalar markerKVol = markerK.getMarkerVolume();
            const DynamicList<std::tuple<bool,label,label>>& supportK = markerK.getSupportCells();
            
            vector distVec = XI-XK;
            scalar distMag = std::sqrt(distVec&distVec);
            
            scalar matrixEntry = 0;
            if(distMag < 10*dilationIMax || distMag < 10*dilationKMax)
            {
                DynamicList<std::pair<vector,scalar>> combinedIKSupport;
                for(const std::tuple<bool,label,label>& suppCell : supportK)
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
                }
                matrixEntry *= markerKVol;
            }
            
            auto iter = rowEntries.find(smProcsNumOfMarkers+K);
            if(iter!=rowEntries.end())
                FatalErrorInFunction<<"Multiple writes in matrix"<<exit(FatalError);
            rowEntries[smProcsNumOfMarkers+K] = matrixEntry;            
        }
        
        // Compute foreign local matrix entries
        std::unordered_set<label> neighborProcesses;
        for(const std::tuple<bool,label,label>& tupl : supportI)
        {
            if(!std::get<0>(tupl))
            {
                label process = std::get<1>(tupl);
                label cellInd = std::get<2>(tupl);
                neighborProcesses.insert(process);
                const std::unordered_map<label,label>& processHaloCellToIndexMap = getHaloCellToIndexMap(process);
                auto cellIndexIter = processHaloCellToIndexMap.find(cellInd);
                if(cellIndexIter==processHaloCellToIndexMap.end())
                    FatalErrorInFunction<<"Cell in process halo not found"<<exit(FatalError);
                label index = cellIndexIter->second;
                const DynamicList<CellDescription>& processHaloCellList = getHaloCellList(process);
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
            
            label neighProcHaloCellNum = globHaloMarkers.size_haloCells(process);
            for(label haloCellInd=0; haloCellInd<neighProcHaloCellNum; haloCellInd++)
            {
                label haloCellMarkerNum = globHaloMarkers.size_cellMarkers(process,haloCellInd);               
                for(label haloCellMarkerInd=0; haloCellMarkerInd<haloCellMarkerNum; haloCellMarkerInd++)
                {
                    std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,List<scalar>> markerData;
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
                        List<scalar>& markerKb = std::get<7>(markerData);
                        std::array<scalar,10> b;
                        for(label i=0; i<markerKb.size(); i++)
                            b[i] = markerKb[i];
                        
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
                            scalar weightI = markerI.correctedDeltaDirac(XI,suppCell.first);
                            scalar weightK = LM::correctedDeltaDirac
                            (
                                XK,suppCell.first,markerKdilation,b
                            );
                            matrixEntry += weightK*weightI*suppCell.second;
                        }
                        matrixEntry *= markerKvolume;
                    }
                    
                    auto iter = rowEntries.find(offset+markerKindex);
                    if(iter!=rowEntries.end())
                        FatalErrorInFunction<<"Multiple writes in matrix"<<exit(FatalError);
                    rowEntries[offset+markerKindex] = matrixEntry;
                }
            }
        }
        
        DynamicList<std::pair<scalar,label>> matrixRow;
        for(auto iter=rowEntries.begin(); iter!=rowEntries.end(); iter++)
        {
            label K = iter->first;
            scalar matrixEntry = iter->second;
            matrixRow.append({matrixEntry,K});
        }
        A.addRow(matrixRow);
        b[I] = 1;
    }
    return result;
}

void Foam::LineStructure::computeMarkerWeights()
{
    std::unique_ptr<Foam::LineStructure::LinearSystem> system = computeMarkerEpsilonMatrix();
    CSR_Matrix_par& A = std::get<0>(*system);
    Vector_par& ones = std::get<1>(*system);
    Info<<"A:"<<A.to_string()<<Foam::endl;
    Info<<"ones:"<<ones.to_string()<<Foam::endl;
    
    BiCGSTAB solver(A);
    Vector_par eps = solver.solve(ones);
    Info<<"eps:"<<eps.to_string()<<Foam::endl;
    
    FatalErrorInFunction<<"Temp stop"<<exit(FatalError);
    for(uint I=0; I<collectedMarkers.size(); I++)
    {
        collectedMarkers[I]->setMarkerWeight(eps[I]);
    }
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

Foam::LineStructure::GlobalHaloMarkers::GlobalHaloMarkers
(
    // [haloCells] -> [marker of cell] -> ptr,index
    std::vector<std::vector<std::pair<LagrangianMarker*,label>>>& selfhaloCellsRodMarkersList
):
selfhaloCellsRodMarkersList(&selfhaloCellsRodMarkersList),
broadcasted(false)
{
    globalHaloCellsMarkerPos.clear();
    globalHaloCellsMarkerPos.setSize(Pstream::nProcs());
    globalHaloCellsMarkerPos[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerPos[Pstream::myProcNo()].setSize(this->selfhaloCellsRodMarkersList->size());
    
    globalHaloCellsMarkerVolume.clear();
    globalHaloCellsMarkerVolume.setSize(Pstream::nProcs());
    globalHaloCellsMarkerVolume[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerVolume[Pstream::myProcNo()].setSize(this->selfhaloCellsRodMarkersList->size());
    
    globalHaloCellsLocalIndex.clear();
    globalHaloCellsLocalIndex.setSize(Pstream::nProcs());
    globalHaloCellsLocalIndex[Pstream::myProcNo()].clear();
    globalHaloCellsLocalIndex[Pstream::myProcNo()].setSize(this->selfhaloCellsRodMarkersList->size());
    
    globalHaloCellsMarkerDilation.clear();
    globalHaloCellsMarkerDilation.setSize(Pstream::nProcs());
    globalHaloCellsMarkerDilation[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerDilation[Pstream::myProcNo()].setSize(this->selfhaloCellsRodMarkersList->size());
    
    globalHaloCellsMarkerSupportCellIndices.clear();
    globalHaloCellsMarkerSupportCellIndices.setSize(Pstream::nProcs());
    globalHaloCellsMarkerSupportCellIndices[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerSupportCellIndices[Pstream::myProcNo()].setSize(this->selfhaloCellsRodMarkersList->size());
    
    globalHaloCellsMarkerSupportCellCentres.clear();
    globalHaloCellsMarkerSupportCellCentres.setSize(Pstream::nProcs());
    globalHaloCellsMarkerSupportCellCentres[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerSupportCellCentres[Pstream::myProcNo()].setSize(this->selfhaloCellsRodMarkersList->size());
    
    globalHaloCellsMarkerSupportCellVolume.clear();
    globalHaloCellsMarkerSupportCellVolume.setSize(Pstream::nProcs());
    globalHaloCellsMarkerSupportCellVolume[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerSupportCellVolume[Pstream::myProcNo()].setSize(this->selfhaloCellsRodMarkersList->size());
    
    globalHaloCellsMarkerb.clear();
    globalHaloCellsMarkerb.setSize(Pstream::nProcs());
    globalHaloCellsMarkerb[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerb[Pstream::myProcNo()].setSize(this->selfhaloCellsRodMarkersList->size());
}

void Foam::LineStructure::GlobalHaloMarkers::appendMarkerData
(
    label haloCellInd,
    std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,List<scalar>> data
)
{
    if(broadcasted)
        FatalErrorInFunction<<"Already broadcasted"<<Foam::endl;
    
    label proc = Pstream::myProcNo();
    label hCI = haloCellInd;
    globalHaloCellsMarkerPos[proc][hCI].append(std::get<0>(data));
    globalHaloCellsMarkerVolume[proc][hCI].append(std::get<1>(data));
    globalHaloCellsLocalIndex[proc][hCI].append(std::get<2>(data));
    globalHaloCellsMarkerDilation[proc][hCI].append(std::get<3>(data));
    globalHaloCellsMarkerSupportCellIndices[proc][hCI].append(std::get<4>(data));
    globalHaloCellsMarkerSupportCellCentres[proc][hCI].append(std::get<5>(data));
    globalHaloCellsMarkerSupportCellVolume[proc][hCI].append(std::get<6>(data));
    globalHaloCellsMarkerb[proc][hCI].append(std::get<7>(data));
}

std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,List<scalar>> Foam::LineStructure::GlobalHaloMarkers::getMarkerData
(
    label process,
    label haloCellInd,
    label cellMarkerInd
) const
{
    if(!broadcasted)
        FatalErrorInFunction<<"Not broadcasted yet"<<Foam::endl;
        
    label proc = process;
    label hCI = haloCellInd;
    label cMI = cellMarkerInd;
    vector position = globalHaloCellsMarkerPos[proc][hCI][cMI];
    scalar volume = globalHaloCellsMarkerVolume[proc][hCI][cMI];
    label index = globalHaloCellsLocalIndex[proc][hCI][cMI];
    vector dilation = globalHaloCellsMarkerDilation[proc][hCI][cMI];
    DynamicList<Pair<label>> suppCellsIndices = globalHaloCellsMarkerSupportCellIndices[proc][hCI][cMI];
    DynamicList<vector> suppCellsCentre = globalHaloCellsMarkerSupportCellCentres[proc][hCI][cMI];
    DynamicList<scalar> suppCellsVolume = globalHaloCellsMarkerSupportCellVolume[proc][hCI][cMI];
    List<scalar> b = globalHaloCellsMarkerb[proc][hCI][cMI];
    return {position,volume,index,dilation,suppCellsIndices,suppCellsCentre,suppCellsVolume,b};
}

label Foam::LineStructure::GlobalHaloMarkers::size_processes() const
{
    return Pstream::nProcs();
}

label Foam::LineStructure::GlobalHaloMarkers::size_haloCells
(
    label process
) const
{
    return globalHaloCellsMarkerPos[process].size();
}

label Foam::LineStructure::GlobalHaloMarkers::size_cellMarkers
(
    label process,
    label haloCellInd
) const
{
    return globalHaloCellsMarkerPos[process][haloCellInd].size();
}

void Foam::LineStructure::GlobalHaloMarkers::communicate()
{
    broadcasted = true;
    
    Pstream::gatherList(globalHaloCellsMarkerPos);
    Pstream::gatherList(globalHaloCellsMarkerVolume);
    Pstream::gatherList(globalHaloCellsMarkerDilation);
    Pstream::gatherList(globalHaloCellsMarkerSupportCellIndices);
    Pstream::gatherList(globalHaloCellsMarkerSupportCellCentres);
    Pstream::gatherList(globalHaloCellsMarkerSupportCellVolume);
    Pstream::gatherList(globalHaloCellsMarkerb);
    
    Pstream::scatterList(globalHaloCellsMarkerPos);
    Pstream::scatterList(globalHaloCellsMarkerVolume);
    Pstream::scatterList(globalHaloCellsMarkerDilation);
    Pstream::scatterList(globalHaloCellsMarkerSupportCellIndices);
    Pstream::scatterList(globalHaloCellsMarkerSupportCellCentres);
    Pstream::scatterList(globalHaloCellsMarkerSupportCellVolume);
    Pstream::scatterList(globalHaloCellsMarkerb);
}

template<typename T> std::unique_ptr<List<List<DynamicList<T>>>> Foam::LineStructure::GlobalHaloMarkers::broadcasHaloMarkerField
(
    const List<T>& collectedMarkerField
)
{
    if(selfhaloCellsRodMarkersList==nullptr)
        FatalErrorInFunction<<"Invalid GlobalHaloMarkers object"<<Foam::endl;
    
    auto resultPtr = std::make_unique<List<List<DynamicList<T>>>>(Pstream::nProcs());
    List<List<DynamicList<T>>>& result = *resultPtr;
    List<DynamicList<T>>& myProcResult = result[Pstream::myProcNo()];
    myProcResult.setSize(selfhaloCellsRodMarkersList->size());
    for(uint haloCellInd=0; haloCellInd<selfhaloCellsRodMarkersList->size(); haloCellInd++)
    {
        std::vector<std::pair<LagrangianMarker*,label>>& oneHaloCellMarkerList = (*selfhaloCellsRodMarkersList)[haloCellInd];
        DynamicList<T> oneHaloCellMarkerField = myProcResult[haloCellInd];
        for(uint markerInd=0; markerInd<oneHaloCellMarkerList.size(); markerInd++)
        {
            std::pair<LagrangianMarker*,label>& oneMarker = oneHaloCellMarkerList[markerInd];
            label index = oneMarker.second;
            oneHaloCellMarkerField.append(collectedMarkerField[index]);
        }
    }
    Pstream::gatherList(result);
    Pstream::scatterList(result);
    return resultPtr;
}

template<typename T> std::unique_ptr<List<List<DynamicList<scalar>>>> Foam::LineStructure::broadcastHaloMarkerFields
(
    const List<T>& collectedMarkerField
)
{
    if(collectedMarkerField.size()!=collectedMarkers.size())
        FatalErrorInFunction<<"Marker Data size mismatch"<< exit(FatalError);
    return globHaloMarkers.broadcasHaloMarkerField<T>(collectedMarkerField);
}

/*
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
*/

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
