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
    
    haloCellsRodMarkersList.resize(selfHaloCellList.size());
    
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
    
    createMarkersFromSpacedPoints();
    refineMarkers();
    setMarkerVolume();
    for(LagrangianMarker& marker : *(rodMarkersList[0]))
        Info<<marker.to_string()<<Foam::endl;
    
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[0];
    Info<<"dist:0.375->0.40625:"<<LineStructure::distance(oneRod,0.375,0.40625)<<Foam::endl;
    evaluateMarkerMeshRelation();
    reduceMarkers();
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
    if(myMesh->m_nR!=rodMarkersList.size())
        FatalErrorInFunction<<"Mismatch in size of m_Rods and rodMarkersList"<<exit(FatalError);
    if(crossSecArea.size()!=rodMarkersList.size())
        FatalErrorInFunction<<"Mismatch in size of crossSecArea and rodMarkersList"<<exit(FatalError);
    if(myMesh->m_Rods.size()!=myMesh->m_nR)
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
                auto inserted = markers.insert(markersIter1,middleMarker);
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
    haloCellsRodMarkersList.resize(selfHaloCellList.size());
    for(label localMarkerInd=0; localMarkerInd<collectedMarkers.size(); localMarkerInd++)
    {
        LagrangianMarker* marker = collectedMarkers[localMarkerInd];
        label cellOfMarker = marker->getMarkerCell();
        auto iter = selfHaloCellToIndex.find(cellOfMarker);
        if(iter!=selfHaloCellToIndex.end())
        {
            label index = iter->second;
            haloCellsRodMarkersList[index].push_back({&marker,localMarkerInd});
        }
    }
}

void Foam::LineStructure::exchangeHaloMarkersData()
{
    globHaloMarkers = GlobalHaloMarkers(haloCellsRodMarkersList);
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
        
    for(label haloCellInd=0; haloCellInd<haloCellsRodMarkersList.size(); haloCellInd++)
    {
        std::vector<std::pair<LagrangianMarker*,label>>& localHaloCellRodMarkers = haloCellsRodMarkersList[haloCellInd];
        for(label haloMarkerInd=0; haloMarkerInd<haloCellsRodMarkersList[haloCellInd].size(); haloMarkerInd++)
        {
            std::pair<LagrangianMarker*,label> markerData = localHaloCellRodMarkers[haloMarkerInd];
            LagrangianMarker* marker = markerData.first;
            label markerInd = markerData.second;
            
            vector position = marker->getMarkerPosition();
            scalar volume = marker->getMarkerVolume();
            label index = markerInd;
            vector dilation = marker->getDilation();
            DynamicList<vector> supportCellsCentre;
            DynamicList<scalar> supportCellsVolume;
            
            const DynamicList<std::tuple<bool,label,label>>& supportCells = marker->getSupportCells();
            for(std::tuple<bool,label,label> oneSupport : supportCells)
            {
                vector oneSuppCellCentre;
                scalar oneSuppCellVolume;
                marker->getCellData(oneSupport,oneSuppCellCentre,oneSuppCellVolume);
                supportCellsCentre.append(oneSuppCellCentre);
                supportCellsVolume.append(oneSuppCellVolume);
            }
            globHaloMarkers.appendMarkerData
            (
                haloCellInd,
                {position,volume,index,vector,supportCellsCentre,supportCellsVolume}
            );
        }
    }
    GlobalHaloMarkers.communicate();
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
        processGlobalMarkerOffset[proc] = processGlobalMarkerOffset[proc-1]+globNumOfMarkers[proc-1];
    
    auto result = std::unique_ptr<LinearSystem>(new LinearSystem());
    std::get<0>(*result) = CSR_Matrix_par(locProcMarkerNbr,smProcsNumOfMarkers,globNumOfMarkers,globNumOfMarkers);
    CSR_Matrix_par& A = std::get<0>(*result);
    std::get<1>(*result) = Vector_par(locProcMarkerNbr,smProcsNumOfMarkers,globNumOfMarkers);
    Vector_par& b = std::get<1>(*result);
    
    for(label I=0; I<markers.size(); I++)
    {
        const LagrangianMarker& markerI = *(markers[I]);
        vector XI = markerI.getMarkerPosition();
        vector dilationI = markerI.getDilation();
        scalar dilationIMax = std::max<scalar>(dilationI[0],dilationI[1]);
        dilationIMax = std::max<scalar>(dilationIMax,dilationI[2]);
        scalar markerIVol = markerI.getMarkerVolume();
        const DynamicList<std::tuple<bool,label,label>>& supportI = markerI.getSupportCells();
        DynamicList<vector> markerISuppCellCentres;
        DynamicList<scalar> markerISuppCellVol;
        for(std::tuple<bool,label,label>& suppCell : supportI)
        {
            vector centre;
            scalar vol;
            markerI.getCellData(suppCell,centre,vol);
            markerISuppCellCentres.append(centre);
            markerISuppCellVol.append(vol);
        }
        
        std::map<label,scalar> rowEntries;
        
        // Compute in matrix entries
        for(label K=0; K<markers.size(); K++)
        {
            const LagrangianMarker& markerK = *(markers[K]);
            vector XK = markerK.getMarkerPosition();
            vector dilationK = markerK.getDilation();
            scalar dilationKMax = std::max<scalar>(dilationK[0],dilationK[1]);
            dilationKMax = std::max<scalar>(dilationKMax,dilationK[2]);
            scalar markerKVol = markerK.getMarkerVolume();
            const DynamicList<std::tuple<bool,label,label>>& supportK = markerK.getSupportCells();
            DynamicList<vector> markerKSuppCellCentres;
            DynamicList<scalar> markerKSuppCellVol;
            for(std::tuple<bool,label,label>& suppCell : supportK)
            {
                vector centre;
                scalar vol;
                markerK.getCellData(suppCell,centre,vol);
                markerKSuppCellCentres.append(centre);
                markerKSuppCellVol.append(vol);
            }

            vector distVec = XI-XK;
            scalar distMag = std::sqrt(distVec&distVec);
            
            scalar matrixEntry = 0;
            if(distMag < 10*dilationIMax || distMag < 10*dilationKMax)
            {
                for(label suppIInd=0; suppIInd<supportI.size(); suppIInd++)
                {
                    scalar weightI = markerI.correctedDeltaDirac(XI,markerISuppCellCentres[suppIInd));
                    scalar weightK = markerK.correctedDeltaDirac(XK,markerISuppCellCentres[suppIInd));                       
                    matrixEntry += weightK*weightI*markerISuppCellVol[suppIInd];
                }
                for(label suppKInd=0; suppKInd<supportK.size(); suppKInd++)
                {
                    scalar weightI = markerI.correctedDeltaDirac(XI,markerKSuppCellCentres[suppKInd));
                    scalar weightK = markerK.correctedDeltaDirac(XK,markerKSuppCellCentres[suppKInd));                       
                    matrixEntry += weightK*weightI*markerKSuppCellVol[suppKInd];
                }
                matrixEntry *= markerKVol;
            }
            
            auto iter = rowEntries.find(smProcsNumOfMarkers+K);
            if(iter!=rowEntries.end())
                FatalErrorInFunction<<"Multiple writes in matrix"<<exit(FatalError);
            rowEntries[smProcsNumOfMarkers+K] = matrixEntry;            
        }
        
        GlobalHaloMarkers& globHaloMarkers = globHaloMarkers;
        
        //Compute additional right hand side entries
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
                    std::tuple<vector,scalar,label,vector,DynamicList<vector>,DynamicList<scalar>> markerData;
                    markerData = globHaloMarkers.getMarkerData(process,haloCellInd,haloCellMarkerInd);
                    
                    vector markerKposition = std::get<0>(markerData);
                    scalar markerKvolume = std::get<1>(markerData);
                    label markerKindex = std::get<2>(markerData);
                    vector markerKdilation = std::get<3>(markerData);
                    DynamicList<vector>& markerKSuppCellCentres = std::get<4>(markerData);
                    DynamicList<scalar>& markerKSuppCellVol = std::get<5>(markerData);
                    List<scalar>& markerKb = std::get<6>(markerData);
                    std::array<scalar,10> b;
                    for(label i=0; i<markerKb.size(); i++)
                        b[i] = markerKb[i];
                    
                    vector XK = markerKposition;
                    scalar markerKVol = markerKvolume;
                    
                    vector distVec = XI-XK;
                    scalar distMag = std::sqrt(distVec&distVec);
                    
                    scalar matrixEntry = 0;
                    if(distMag < 10*dilationIMax || distMag < 10*dilationKMax)
                    {        
                        for(label suppIInd=0; suppIInd<supportI.size(); suppIInd++)
                        {
                            scalar weightI = markerI.correctedDeltaDirac(XI,markerISuppCellCentres[suppIInd));
                            scalar weightK = LM::correctedDeltaDirac
                            (
                                XK,markerISuppCellCentres[suppIInd),markerKdilation,b
                            );                       
                            matrixEntry += weightK*weightI*markerISuppCellVol[suppIInd];
                        }
                        for(label suppKInd=0; suppKInd<supportK.size(); suppKInd++)
                        {
                            scalar weightI = markerI.correctedDeltaDirac(XI,markerKSuppCellCentres[suppKInd));
                            scalar weightK = LM::correctedDeltaDirac
                            (
                                XK,markerKSuppCellCentres[suppKInd),markerKdilation,b
                            );     
                            matrixEntry += weightK*weightI*markerKSuppCellVol[suppKInd];
                        }
                        matrixEntry *= markerKVol;
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
    std::unique_ptr<LinearSystem> system = computeMarkerEpsilonMatrix();
    
    Pout<<"Computed system matrix"<<Foam::endl;
    std::cout<<std::get<1>(*system)<<std::endl;
    
    std::vector<LagrangianMarker*>& markers = std::get<0>(*system);
    const gismo::gsMatrix<scalar>& A = std::get<1>(*system);
    // <markerI,DynList<tuple<procK,haloCellIndK,markerIndK,aIK>>>
    const std::unordered_map<label,DynamicList<std::tuple<label,label,label,scalar>>>& bEntry = std::get<2>(*system);
    Pout<<"Decomposed system matrix"<<Foam::endl;
    
    const DynamicList<CellDescription>& selfHaloCellList = getHaloCellList(Pstream::myProcNo());
    List<List<DynamicList<scalar>>> globalHaloMarkerEpsValues(Pstream::nProcs());
    globalHaloMarkerEpsValues[Pstream::myProcNo()].resize(selfHaloCellList.size());
    
    Pout<<"Created eps vector"<<Foam::endl;

    // markerPtr -> {cellInd, list index}
    std::unordered_map<const LagrangianMarker*,std::pair<label,label>> markerPtrHaloCellIndex;
    List<List<DynamicList<scalar>>> globalHaloCellsMarkerEpsilon(Pstream::nProcs());
    globalHaloCellsMarkerEpsilon[Pstream::myProcNo()].resize(haloCellsRodMarkersList.size());
    List<DynamicList<scalar>>& localHaloCellsMarkerEpsilon = globalHaloCellsMarkerEpsilon[Pstream::myProcNo()];
    
    Pout<<"Loop"<<Foam::endl;
    
    Pout<<"globalHaloCellsMarkerEpsilon.size():"<<globalHaloCellsMarkerEpsilon.size()<<Foam::endl;
    Pout<<"localHaloCellsMarkerEpsilon.size():"<<localHaloCellsMarkerEpsilon.size()<<Foam::endl;

    Pout<<"haloCellsRodMarkersList.size():"<<haloCellsRodMarkersList.size()<<Foam::endl;
        
    for(label haloCellInd=0; haloCellInd<haloCellsRodMarkersList.size(); haloCellInd++)
    {
        Pout<<"haloCellInd:"<<haloCellInd<<Foam::endl;
        DynamicList<scalar>& cellMarkersEpsilon = localHaloCellsMarkerEpsilon[haloCellInd];
        std::vector<const LagrangianMarker*>& cellRodMarkerList = haloCellsRodMarkersList[haloCellInd];
        
        cellMarkersEpsilon.resize(cellRodMarkerList.size());
        Pout<<"--haloCellInd:"<<haloCellInd<<Foam::endl;
        
        for(label haloCellMarkerInd=0; haloCellMarkerInd<cellMarkersEpsilon.size(); haloCellMarkerInd++)
        {
            cellMarkersEpsilon[haloCellMarkerInd] = 0;
            const LagrangianMarker* oneMarker = cellRodMarkerList[haloCellMarkerInd];
            markerPtrHaloCellIndex[oneMarker] = {haloCellInd,haloCellMarkerInd};
        }
    }
    
    Pout<<"Halo marker to map:"<<Pstream::myProcNo()<<"  "<<markerPtrHaloCellIndex.size()<<Foam::endl;
        
    // markerI -> pair<haloCellIndK,markerIndK>>
    std::unordered_map<label,std::pair<label,label>> markerIIndHaloCellIndex;
    for(label I=0; I<markers.size(); I++)
    {
        const LagrangianMarker* markerPtr = markers[I];
        auto iter = markerPtrHaloCellIndex.find(markerPtr);
        if(iter!=markerPtrHaloCellIndex.end())
        {
            markerIIndHaloCellIndex[I] = iter->second;
        }
    }
    
    Pout<<"Halo marker I to map:"<<Pstream::myProcNo()<<"  "<<markerIIndHaloCellIndex.size()<<Foam::endl;
    
    
    //[i] -> {markerI,{haloCellIndK,markerIndK}}
    List<std::pair<label,std::pair<label,label>>> I_haloMarkerCellIndex(markerIIndHaloCellIndex.size());
    label i=0;
    for(auto iter=markerIIndHaloCellIndex.begin(); iter!=markerIIndHaloCellIndex.end(); iter++,i++)
    {
        I_haloMarkerCellIndex[i] = {iter->first,iter->second};
    }
    
    /*
    Pout<<"I:"<<0<<" "<<I_haloMarkerCellIndex[0].first<<"|"<<I_haloMarkerCellIndex[0].second.first<<"|"<<I_haloMarkerCellIndex[0].second.second<<Foam::endl;
    Pout<<"I:"<<1<<" "<<I_haloMarkerCellIndex[1].first<<"|"<<I_haloMarkerCellIndex[1].second.first<<"|"<<I_haloMarkerCellIndex[1].second.second<<Foam::endl;
    Pout<<"I:"<<2<<" "<<I_haloMarkerCellIndex[2].first<<"|"<<I_haloMarkerCellIndex[2].second.first<<"|"<<I_haloMarkerCellIndex[2].second.second<<Foam::endl;
    Pout<<"Halo marker I list:"<<Pstream::myProcNo()<<"  "<<I_haloMarkerCellIndex.size()<<Foam::endl;
    */
        
    gismo::gsMatrix<scalar> ones(A.cols(),1),eps(A.cols(),1);
    Info<<"Start solving equential out"<<Foam::endl;
    if(Pstream::parRun())
    {
        for(label col=0; col<A.cols(); col++)
            eps(col,0) = 0;
        
        for(scalar tolerance = 1e-1; tolerance>=1e-10; tolerance/=10)
        {
            Pout<<"tolerance:"<<tolerance<<Foam::endl;
            bool shortSolution = false;
            while(!shortSolution)
            {
                Pout<<"  Start subiteration"<<Foam::endl;
                gismo::gsGMRes WM(A,tolerance);
                
                // gather and scatter eps marker values
                for(const std::pair<label,std::pair<label,label>>& markerInd : I_haloMarkerCellIndex)
                {
                    label I = markerInd.first;
                    label haloCellInd = markerInd.second.first;
                    label haloCellMarkerInd = markerInd.second.second;
                    scalar markerEps = eps(I,0);
                    globalHaloCellsMarkerEpsilon[Pstream::myProcNo()][haloCellInd][haloCellMarkerInd] = markerEps;
                }
                Pstream::gatherList(globalHaloCellsMarkerEpsilon);
                Pstream::scatterList(globalHaloCellsMarkerEpsilon);
                Pout<<"\t Exchanged"<<Foam::endl;

                for(label col=0; col<A.cols(); col++)
                {
                    ones(col,0) = 1;
                    auto iter = bEntry.find(col);
                    if(iter!=bEntry.end())
                    {
                        for(const std::tuple<label,label,label,scalar>& oneHaloMarkerK : iter->second)
                        {
                            label procK = std::get<0>(oneHaloMarkerK);
                            label haloCellIndK = std::get<1>(oneHaloMarkerK);
                            label markerIndK = std::get<2>(oneHaloMarkerK);
                            if( procK>=globalHaloCellsMarkerEpsilon.size() ||
                                haloCellIndK>=globalHaloCellsMarkerEpsilon[procK].size() ||
                                markerIndK>=globalHaloCellsMarkerEpsilon[procK][haloCellIndK].size()
                              )
                            {
                                Pout<<"procK:"<<procK<<"/"<<globalHaloCellsMarkerEpsilon.size()<<Foam::endl;
                                Pout<<"haloCellIndK:"<<haloCellIndK<<"/"<<globalHaloCellsMarkerEpsilon[procK].size()<<Foam::endl;
                                Pout<<"markerIndK:"<<markerIndK<<"/"<<globalHaloCellsMarkerEpsilon[procK][haloCellIndK].size()<<Foam::endl;
                                FatalErrorInFunction<<"Invalid indices of foreign marker"<<exit(FatalError);
                            }
                            scalar aIK = std::get<3>(oneHaloMarkerK);
                            scalar eps = globalHaloCellsMarkerEpsilon[procK][haloCellIndK][markerIndK];
                            ones(col,0) -= aIK*eps;
                        }
                    }
                }
                Pout<<"\t Added right side"<<Foam::endl;
                // Communicate and modify ones
                std::cout<<"A:"<<A<<std::endl;
                std::cout<<"ones:"<<ones<<std::endl;
                std::cout<<"eps"<<eps<<std::endl;
                WM.solve(ones,eps);
                if(WM.iterations()<=2)
                    shortSolution=true;
                Pout<<"\t Solved"<<Foam::endl;
            }
        }
    }
    else
    {
        Info<<"Sequential out"<<Foam::endl;
        for(label col=0; col<A.cols(); col++)
        {
            ones(col,0) = 1;
        }
        linearSolve(A,eps,ones,1e-60);        
        FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
    }
    
    for(const cell& oneCell : mesh.cells())
        Info<<"oneCell:"<<oneCell.centre(mesh.points(),mesh.faces())<<" -> "<<oneCell.mag(mesh.points(),mesh.faces())<<Foam::endl;
    Info<<Foam::endl;
    Info<<"para / cell /   volume   /      position       /       dilation"<<Foam::endl;
    for(auto markers : (*rodMarkersList[0]))
        Info<<markers.getMarkerParameter()<<"  /  "<<markers.getMarkerCell()<<"  /  "<<markers.getMarkerVolume()<<"  /  "<<markers.getMarkerPosition()<<"  /  "<<markers.getDilation()<<Foam::endl;
    
    
    for(label I=0; I<markers.size(); I++)
    {
        
    }
}

void Foam::LineStructure::computeMarkerCellWeights()
{
    std::vector<LagrangianMarker*> allMarkers;
    for(std::unique_ptr<std::list<LagrangianMarker>>& singleRodMarkersPtr : rodMarkersList)
    {
        for(LagrangianMarker& marker : *singleRodMarkersPtr)
        {
            allMarkers.push_back(&marker);
        }
    }
    computeMarkerCellWeights(allMarkers);
}

void Foam::LineStructure::computeMarkerCellWeights
(
    const std::vector<LagrangianMarker*>& markers
)
{
    Info<<"computeMarkerCellWeights"<<Foam::endl;
    for(LagrangianMarker* markerPtr : markers)
    {
        markerPtr->compNonUniformCorrWeights();
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

/*
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
*/

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

Foam::GlobalHaloMarkers::GlobalHaloMarkers
(
    // [haloCells] -> [marker of cell] -> ptr,index
    std::vector<std::vector<std::pair<LagrangianMarker*,label>>>& selfhaloCellsRodMarkersList
):
selfhaloCellsRodMarkersList(selfhaloCellsRodMarkersList),
broadcasted(false)
{
    globalHaloCellsMarkerPos.clear();
    globalHaloCellsMarkerPos.setSize(Pstream::nProcs());
    globalHaloCellsMarkerPos[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerPos[Pstream::myProcNo()].setSize(selfhaloCellsRodMarkersList);
    
    globalHaloCellsMarkerVolume.clear();
    globalHaloCellsMarkerVolume.setSize(Pstream::nProcs());
    globalHaloCellsMarkerVolume[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerVolume[Pstream::myProcNo()].setSize(selfhaloCellsRodMarkersList);
    
    globalHaloCellsLocalIndex.clear();
    globalHaloCellsLocalIndex.setSize(Pstream::nProcs());
    globalHaloCellsLocalIndex[Pstream::myProcNo()].clear();
    globalHaloCellsLocalIndex[Pstream::myProcNo()].setSize(selfhaloCellsRodMarkersList);
    
    globalHaloCellsMarkerDilation.clear();
    globalHaloCellsMarkerDilation.setSize(Pstream::nProcs());
    globalHaloCellsMarkerDilation[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerDilation[Pstream::myProcNo()].setSize(selfhaloCellsRodMarkersList);
    
    globalHaloCellsMarkerSupportCellCentres.clear();
    globalHaloCellsMarkerSupportCellCentres.setSize(Pstream::nProcs());
    globalHaloCellsMarkerSupportCellCentres[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerSupportCellCentres[Pstream::myProcNo()].setSize(selfhaloCellsRodMarkersList);
    
    globalHaloCellsMarkerSupportCellVolume.clear();
    globalHaloCellsMarkerSupportCellVolume.setSize(Pstream::nProcs());
    globalHaloCellsMarkerSupportCellVolume[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerSupportCellVolume[Pstream::myProcNo()].setSize(selfhaloCellsRodMarkersList);
    
    globalHaloCellsMarkerb.clear();
    globalHaloCellsMarkerb.setSize(Pstream::nProcs());
    globalHaloCellsMarkerb[Pstream::myProcNo()].clear();
    globalHaloCellsMarkerb[Pstream::myProcNo()].setSize(selfhaloCellsRodMarkersList);
}

void appendMarkerData
(
    label haloCellInd,
    std::tuple<vector,scalar,label,vector,DynamicList<vector>,DynamicList<scalar>,List<scalar>> data
)
{
    label proc = Pstream::myProcNo();
    label hCI = haloCellInd;
    globalHaloCellsMarkerPos[proc][hCI].append(std::get<0>(data));
    globalHaloCellsMarkerVolume[proc][hCI].append(std::get<1>(data));
    globalHaloCellsLocalIndex[proc][hCI].append(std::get<2>(data));
    globalHaloCellsMarkerDilation[proc][hCI].append(std::get<3>(data));
    globalHaloCellsMarkerSupportCellCentres[proc][hCI].append(std::get<4>(data));
    globalHaloCellsMarkerSupportCellVolume[proc][hCI].append(std::get<5>(data));
    globalHaloCellsMarkerb[proc][hCI].append(std::get<6>(data));
}

std::tuple<vector,scalar,label,vector,DynamicList<vector>,DynamicList<scalar>> getMarkerData
(
    label process,
    label haloCellInd,
    label cellMarkerInd
)
{
    label proc = process;
    label hCI = haloCellInd;
    label cMI = cellMarkerInd;
    vector position = globalHaloCellsMarkerPos[proc][hCI][cMI];
    scalar volume = globalHaloCellsMarkerVolume[proc][hCI][cMI];
    label index = globalHaloCellsLocalIndex[proc][hCI][cMI];
    vector dilation = globalHaloCellsMarkerDilation[proc][hCI][cMI];
    DynamicList<vector> suppCellsCentre = globalHaloCellsMarkerSupportCellCentres[proc][hCI][cMI];
    DynamicList<scalar> suppCellsVolume = globalHaloCellsMarkerSupportCellVolume[proc][hCI][cMI];
    List<scalar> b = globalHaloCellsMarkerb[proc][hCI][cMI];
    return {position,volume,index,dilation,suppCellsCentre,suppCellsVolume,b};
}

label size_processes()
{
    return Pstream::nProcs();
}

label size_haloCells(label process)
{
    return globalHaloCellsMarkerPos[process].size();
}

label size_cellMarkers(label process, label haloCellInd)
{
    return globalHaloCellsMarkerPos[process][haloCellInd].size();
}

void Foam::GlobalHaloMarkers::communicate()
{
    broadcasted = true;
    
    Pstream::gatherList(globalHaloCellsMarkerPos);
    Pstream::gatherList(globalHaloCellsMarkerVolume);
    Pstream::gatherList(globalHaloCellsMarkerDilation);
    Pstream::gatherList(globalHaloCellsMarkerSupportCellCentres);
    Pstream::gatherList(globalHaloCellsMarkerSupportCellVolume);
    Pstream::gatherList(globalHaloCellsMarkerb);
    
    Pstream::scatterList(globalHaloCellsMarkerPos);
    Pstream::scatterList(globalHaloCellsMarkerVolume);
    Pstream::scatterList(globalHaloCellsMarkerDilation);
    Pstream::scatterList(globalHaloCellsMarkerSupportCellCentres);
    Pstream::scatterList(globalHaloCellsMarkerSupportCellVolume);
    Pstream::scatterList(globalHaloCellsMarkerb);
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
