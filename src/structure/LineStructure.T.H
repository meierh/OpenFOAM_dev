template<class MARKER>
void Foam::LineStructure::reduceMarkers
(
    std::vector<MarkerReference<MARKER>>& allMarkers
)
{
    using SingleMarker = MarkerReference<MARKER>;
    
    std::unordered_map<label,std::vector<SingleMarker*>> cellToMarker;
    for(SingleMarker& markerData : allMarkers)
    {
        const MARKER& marker = markerData.getMarker();
        label cellInd = marker.getMarkerCell();
        cellToMarker[cellInd].push_back(&markerData);
    }
    oneCenteredMarker(allMarkers,cellToMarker);
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

template<typename T> std::unique_ptr<List<List<DynamicList<T>>>> Foam::LineStructure::broadcastHaloMarkerField
(
    const List<T>& collectedMarkerField
)
{
    if(collectedMarkerField.size()!=static_cast<label>(collectedMarkers.size()))
        FatalErrorInFunction<<"Marker Data size mismatch"<< exit(FatalError);
    return globHaloMarkers.broadcasHaloMarkerField<T>(collectedMarkerField);
}

template<typename T>
T Foam::LineStructure::integrateRodwise
(
    const ActiveRodMesh::rodCosserat* oneRod,
    std::function<T(scalar)> function
)
{
    return integrateRodwise<T>(oneRod,oneRod->m_Curve.knots().first(),oneRod->m_Curve.knots().last(),function);
}

template<typename T>
T Foam::LineStructure::integrateRodwise
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parA,
    scalar parB,
    std::function<T(scalar)> function
)
{
    if(parA>parB)
        FatalErrorInFunction<<"Integration from higher to lower value"<<exit(FatalError);
    
    std::set<scalar> knotSet;
    for(scalar knot : oneRod->m_Curve.knots())
        knotSet.insert(knot);
    
    label degree = oneRod->m_Curve.knots().degree();
    auto start = knotSet.insert(parA);
    auto end = knotSet.insert(parB);
    
    label numberOfAbcissa = (degree+1)*(degree+1);
    
    T totalValue = Foam::zero();
    
    std::set<scalar>::iterator currNode = start.first;
    std::set<scalar>::iterator nextNode = ++(start.first);
    for(; currNode!=knotSet.end() && nextNode!=knotSet.end() && currNode!=end.first ; )
    {
        scalar startPar = *currNode;
        scalar endPar = *nextNode;
        scalar dist = endPar-startPar;
        scalar step = dist/numberOfAbcissa;
        scalar initialStep = startPar+step/2;
        T summedValue = Foam::zero();
        for(label i=0;i<numberOfAbcissa;i++)
        {
            summedValue += function(initialStep)*step;
            initialStep+=step;
        }
        totalValue+=summedValue;
        
        currNode=nextNode;
        nextNode++;
    }
    return totalValue;
}
