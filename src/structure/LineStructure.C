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
    /*
    Info<<"modusFieldToMarker:"<<modusFieldToMarker<<Foam::endl;
    Info<<"modusMarkerToField:"<<modusMarkerToField<<Foam::endl;
    */
    
    label nbrOfRods = myMesh->m_nR;
    rodMarkersList.resize(nbrOfRods);
    check();
    initialRodPoints.resize(nbrOfRods);
    createSpacingPoints();
    
    /*
    for(const auto& oneRodPoints : initialRodPoints)
        for(scalar pnt : *oneRodPoints)
            Info<<pnt<<Foam::endl;
    */
    
    const DynamicList<CellDescription>& selfHaloCellList = getHaloCellList(Pstream::myProcNo());
    
    createMarkersFromSpacedPoints();
    refineMarkers();
    setMarkerVolume();
    
    /*
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[0];
    Info<<"dist:0.375->0.40625:"<<LineStructure::distance(oneRod,0.375,0.40625)<<Foam::endl;
    */
    evaluateMarkerMeshRelation();
    reduceMarkers();
    collectMarkers();
    collectHaloMarkers();
    exchangeHaloMarkersData();

    computeMarkerCellWeights();    
    for(LagrangianMarker& marker : *(rodMarkersList[0]))
        Pout<<marker.to_string()<<Foam::endl;
    
    computeMarkerWeights();
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
        if(modusMarkerToField==markerMeshType::NonUniform)
        {
            iter->setMarkerVolume(span);
        }
        else
            iter->setMarkerVolume(span*crossSecArea[rodNumber]);
        
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
    
    Info<<"A:"<<A.to_string()<<Foam::endl;
    
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
    Pout<<"collectedMarkers.size():"<<collectedMarkers.size()<<Foam::endl;
    for(uint I=0; I<collectedMarkers.size(); I++)
    {
        Pout<<"I:"<<I<<Foam::endl;
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
            Pout<<"I:"<<I<<" K:"<<K<<Foam::endl;

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
            Pout<<"neighbour:"<<process<<Foam::endl;
            label neighProcHaloCellNum = globHaloMarkers.size_haloCells(process);
            for(label haloCellInd=0; haloCellInd<neighProcHaloCellNum; haloCellInd++)
            {
                Pout<<"neighbour:"<<process<<"  haloCellInd:"<<haloCellInd<<"/"<<neighProcHaloCellNum<<Foam::endl;
                label haloCellMarkerNum = globHaloMarkers.size_cellMarkers(process,haloCellInd);               
                for(label haloCellMarkerInd=0; haloCellMarkerInd<haloCellMarkerNum; haloCellMarkerInd++)
                {
                    Pout<<"neighbour:"<<process<<"  haloCellInd:"<<haloCellInd<<"  haloCellMarkerInd:"<<haloCellMarkerInd<<"/"<<haloCellMarkerNum<<Foam::endl;
                    
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
    Pout<<"A:"<<A.to_string()<<Foam::endl;
    
    return result;
}

void Foam::LineStructure::computeMarkerWeights()
{
    Info<<"Compute marker weights"<<Foam::endl;
    std::unique_ptr<Foam::LineStructure::LinearSystem> system = computeMarkerEpsilonMatrix();
    Info<<"Computed marker weights"<<Foam::endl;
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
    
    nProcs = Pstream::nProcs();
    
    procHaloCellsSize[Pstream::myProcNo()] = selfHaloCellSize;
    Pstream::gatherList(procHaloCellsSize);
    Pstream::scatterList(procHaloCellsSize);
    
    procHaloCellMarkerSize.setSize(nProcs);
    for(label proc=0; proc<nProcs; procs++)
        procHaloCellMarkerSize[procs].setSize(procHaloCellsSize[proc]);
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
    if(procHaloCellsSize[proc]>=haloCellInd)
        FatalErrorInFunction<<"Out of range value for halocCellInd"<<Foam::endl;
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
    Pout<<"-----------------------getMarkerData--------------------------"<<Foam::endl;
    if(!broadcasted)
        FatalErrorInFunction<<"Not broadcasted yet"<<Foam::endl;

    label proc = process;
    label hCI = haloCellInd;
    label cMI = cellMarkerInd;

    Pout<<"("<<proc<<","<<hCI<<","<<cMI<<")"<<Foam::endl;
    
    vector position = globalHaloCellsMarkerPos[proc][hCI][cMI];
    Pout<<"---"<<Foam::endl;
    scalar volume = globalHaloCellsMarkerVolume[proc][hCI][cMI];
    Pout<<"---|---"<<Foam::endl;
    Pout<<"globalHaloCellsLocalIndex.size():"<<globalHaloCellsLocalIndex.size()<<Foam::endl;
    Pout<<"globalHaloCellsLocalIndex["<<proc<<"].size():"<<globalHaloCellsLocalIndex[proc].size()<<Foam::endl;
    Pout<<"globalHaloCellsLocalIndex["<<proc<<"]["<<hCI<<"].size():"<<globalHaloCellsLocalIndex[proc][hCI].size()<<Foam::endl;
    label index = globalHaloCellsLocalIndex[proc][hCI][cMI];
    Pout<<"---|---|---"<<Foam::endl;
    vector dilation = globalHaloCellsMarkerDilation[proc][hCI][cMI];
    Pout<<"---|---|---|---"<<Foam::endl;
    DynamicList<Pair<label>> suppCellsIndices = globalHaloCellsMarkerSupportCellIndices[proc][hCI][cMI];
    Pout<<"---|---|---|---|---"<<Foam::endl;
    DynamicList<vector> suppCellsCentre = globalHaloCellsMarkerSupportCellCentres[proc][hCI][cMI];
    Pout<<"---|---|---|---|---|---"<<Foam::endl;
    DynamicList<scalar> suppCellsVolume = globalHaloCellsMarkerSupportCellVolume[proc][hCI][cMI];
    Pout<<"---|---|---|---|---|---|---"<<Foam::endl;
    List<scalar> b = globalHaloCellsMarkerb[proc][hCI][cMI];
    Pout<<"---|---|---|---|---|---|---|---"<<Foam::endl;
    
    Pout<<"||||||||||||||||||getMarkerData|||||||||||||||||||||"<<Foam::endl;
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
    return procHaloCellsSize[process];
}

label Foam::LineStructure::GlobalHaloMarkers::size_cellMarkers
(
    label process,
    label haloCellInd
) const
{
    return procHaloCellMarkerSize[process][haloCellInd];
}

void Foam::LineStructure::GlobalHaloMarkers::communicate()
{
    check();
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
    check();
}

void Foam::LineStructure::GlobalHaloMarkers::check()
{
    label numProcs = size_processes();
    if(numProcs!=globalHaloCellsMarkerPos.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerPos size error"<<Foam::endl;
    if(numProcs!=globalHaloCellsMarkerVolume.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerVolume size error"<<Foam::endl;
    if(numProcs!=globalHaloCellsMarkerDilation.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerDilation size error"<<Foam::endl;
    if(numProcs!=globalHaloCellsMarkerSupportCellIndices.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellIndices size error"<<Foam::endl;
    if(numProcs!=globalHaloCellsMarkerSupportCellCentres.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellCentres size error"<<Foam::endl;
    if(numProcs!=globalHaloCellsMarkerSupportCellVolume.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellVolume size error"<<Foam::endl;
    if(numProcs!=globalHaloCellsMarkerb.size())
        FatalErrorInFunction<<"globalHaloCellsMarkerb size error"<<Foam::endl;
    
    for(label proc=0; proc<numProcs; proc++)
    {
        label numHaloCells = size_haloCells(proc);
        if(numHaloCells!=globalHaloCellsMarkerPos[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerPos["<<proc<<"] size error"<<Foam::endl;
        if(numHaloCells!=globalHaloCellsMarkerVolume[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerVolume["<<proc<<"] size error"<<Foam::endl;
        if(numHaloCells!=globalHaloCellsMarkerDilation[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerDilation["<<proc<<"] size error"<<Foam::endl;
        if(numHaloCells!=globalHaloCellsMarkerSupportCellIndices[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellIndices["<<proc<<"] size error"<<Foam::endl;
        if(numHaloCells!=globalHaloCellsMarkerSupportCellCentres[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellCentres["<<proc<<"] size error"<<Foam::endl;
        if(numHaloCells!=globalHaloCellsMarkerSupportCellVolume[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellVolume["<<proc<<"] size error"<<Foam::endl;
        if(numHaloCells!=globalHaloCellsMarkerb[proc].size())
            FatalErrorInFunction<<"globalHaloCellsMarkerb["<<proc<<"] size error"<<Foam::endl;

        for(label haloCellInd=0; haloCellInd<numHaloCells; haloCellInd++)
        {
            label numHaloCellMarkers = size_cellMarkers(proc,haloCellInd);
            if(numHaloCellMarkers!=globalHaloCellsMarkerPos[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerPos["<<proc<<"]["<<haloCellInd<<"] size error"<<Foam::endl;
            if(numHaloCellMarkers!=globalHaloCellsMarkerVolume[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerVolume["<<proc<<"]["<<haloCellInd<<"] size error"<<Foam::endl;
            if(numHaloCellMarkers!=globalHaloCellsMarkerDilation[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerDilation["<<proc<<"]["<<haloCellInd<<"] size error"<<Foam::endl;
            if(numHaloCellMarkers!=globalHaloCellsMarkerSupportCellIndices[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellIndices["<<proc<<"]["<<haloCellInd<<"] size error"<<Foam::endl;
            if(numHaloCellMarkers!=globalHaloCellsMarkerSupportCellCentres[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellCentres["<<proc<<"]["<<haloCellInd<<"] size error"<<Foam::endl;
            if(numHaloCellMarkers!=globalHaloCellsMarkerSupportCellVolume[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerSupportCellVolume["<<proc<<"]["<<haloCellInd<<"] size error"<<Foam::endl;
            if(numHaloCellMarkers!=globalHaloCellsMarkerb[proc][haloCellInd].size())
                FatalErrorInFunction<<"globalHaloCellsMarkerb["<<proc<<"]["<<haloCellInd<<"] size error"<<Foam::endl;
        }
    }
}
