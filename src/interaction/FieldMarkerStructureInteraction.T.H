template<typename T>
void Foam::FieldMarkerStructureInteraction::fieldToMarker
(
    const Field<T>& fieldData,
    DynamicList<T>& markerData
)
{
    if(fieldData.size()!=mesh.cells().size())
    {
        Info<<"fieldData.size():"<<fieldData.size()<<Foam::nl;
        Info<<"mesh.size():"<<mesh.cells().size()<<Foam::nl;
        FatalErrorInFunction<<"Mismatch in cell number"<<exit(FatalError);
    }
    
    std::function<T(const LagrangianMarker*,std::function<T(Pair<label>)> valueFunction)> convolution;
    convolution = [this](const LagrangianMarker* marker, std::function<T(Pair<label>)> valueFunction)
    {
        std::function<scalar(vector,vector)> weightFunction;
        if(modusFieldToMarker==markerMeshType::Uniform)
        {
            weightFunction = [&marker] (vector X, vector x)
            {
                return marker->deltaDirac(X,x);
            };
        }
        else
        {
            weightFunction = [&marker] (vector X, vector x)
            {
                return marker->correctedDeltaDirac(X,x);
            };
        }
        return marker->convolute<T>(weightFunction,valueFunction);
    };
    
    integrateFieldOverMarker(fieldData,markerData,convolution);
    
    /*
    LineStructure& structure = this->structure;
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    
    std::unique_ptr<List<List<T>>> haloFieldsPtr = structure.broadcastHaloFields(fieldData);
    const List<List<T>>& haloFields = *haloFieldsPtr;

    markerData.resize(markers.size());
    std::function<T(Pair<label>)> valueFunction = 
    [&fieldData,&haloFields,&structure](Pair<label> cell)
    {
        const Pair<label>& suppCellData = cell;
        label proc = suppCellData.first();
        label cellInd = suppCellData.second();
        if(proc==Pstream::myProcNo())
        {
            return fieldData[cellInd];
        }
        else
        {
            const std::unordered_map<label,label>& procCellToInd = structure.getHaloCellToIndexMap(proc);
            auto iterIndex = procCellToInd.find(cellInd);
            if(iterIndex!=procCellToInd.end())
                FatalErrorInFunction<<"Halo cell does not exist"<<exit(FatalError);
            label haloIndex = iterIndex->second;
            return haloFields[proc][haloIndex];                    
        }
    };
    for(label index=0; index<markerData.size(); index++)
    {
        const LagrangianMarker& oneMarker = *(markers[index]);
        std::function<scalar(vector,vector)> weightFunction;
        if(modusFieldToMarker==markerMeshType::Uniform)
        {
            weightFunction = [&oneMarker] (vector X, vector x)
            {
                return oneMarker.deltaDirac(X,x);
            };
        }
        else
        {
            weightFunction = [&oneMarker] (vector X, vector x)
            {
                return oneMarker.correctedDeltaDirac(X,x);
            };
        }
        markerData[index] = oneMarker.convolute<T>(weightFunction,valueFunction);
    }
    */
}

template<typename T>
void Foam::FieldMarkerStructureInteraction::integrateFieldOverMarker
(
    const Field<T>& fieldData,
    DynamicList<T>& markerData,
    std::function<T(const LagrangianMarker*, std::function<T(Pair<label>)>)> integration
)
{
    if(fieldData.size()!=mesh.cells().size())
        FatalErrorInFunction<<"Mismatch in cell number"<<exit(FatalError);
    
    LineStructure& structure = this->structure;
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    
    std::unique_ptr<List<List<T>>> haloFieldsPtr = structure.broadcastHaloFields(fieldData);
    const List<List<T>>& haloFields = *haloFieldsPtr;

    markerData.resize(markers.size());
    std::function<T(Pair<label>)> valueFunction = 
    [&fieldData,&haloFields,&structure](Pair<label> cell)
    {
        const Pair<label>& suppCellData = cell;
        label proc = suppCellData.first();
        label cellInd = suppCellData.second();
        if(proc==Pstream::myProcNo())
        {
            return fieldData[cellInd];
        }
        else
        {
            const std::unordered_map<label,label>& procCellToInd = structure.getHaloCellToIndexMap(proc);
            auto iterIndex = procCellToInd.find(cellInd);
            if(iterIndex==procCellToInd.end())
            {
                Pout<<"proc:"<<proc<<Foam::nl;
                Pout<<"cellInd:"<<cellInd<<Foam::nl;
                //Pout<<"iterIndex->second:"<<iterIndex->second<<Foam::nl;
                Pout<<"haloFields.size():"<<haloFields.size()<<Foam::nl;
                Pout<<"haloFields["<<proc<<"].size():"<<haloFields[proc].size()<<Foam::nl;
                FatalErrorInFunction<<"Halo cell does not exist"<<exit(FatalError);
            }
            label haloIndex = iterIndex->second;
            if(haloIndex<0 || haloIndex>=haloFields[proc].size())
            {
                Pout<<"haloIndex:"<<haloIndex<<Foam::nl;
                Pout<<"haloFields.size():"<<haloFields.size()<<Foam::nl;
                Pout<<"haloFields["<<proc<<"].size():"<<haloFields[proc].size()<<Foam::nl;
                FatalErrorInFunction<<"Halo cell out of range"<<exit(FatalError);
            }
            return haloFields[proc][haloIndex];                    
        }
    };
    for(label index=0; index<markerData.size(); index++)
    {
        const LagrangianMarker* oneMarker = markers[index];
        markerData[index] = integration(oneMarker,valueFunction);

        /*
        const LagrangianMarker& oneMarker = *(markers[index]);
        std::function<scalar(vector,vector)> weightFunction;
        if(modusFieldToMarker==markerMeshType::Uniform)
        {
            weightFunction = [&oneMarker] (vector X, vector x)
            {
                return oneMarker.deltaDirac(X,x);
            };
        }
        else
        {
            weightFunction = [&oneMarker] (vector X, vector x)
            {
                return oneMarker.correctedDeltaDirac(X,x);
            };
        }
        markerData[index] = oneMarker.convolute<T>(weightFunction,valueFunction);
        */
    }
}

template<typename T>
void Foam::FieldMarkerStructureInteraction::markerToField
(
    const DynamicList<T>& markerData,
    Field<T>& fieldData
)
{
    if(fieldData.size()!=mesh.cells().size())
        FatalErrorInFunction<<"Mismatch in cell number"<<exit(FatalError);
    
    fieldData = Foam::zero();
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    const std::vector<LagrangianMarker*>& markers = this->structure.getCollectedMarkers();
    
    // Spread marker field of own markers 
    for(label index=0; index<markerData.size(); index++)
    {
        const LagrangianMarker& oneMarker = *(markers[index]);
        const vector& oneMarkerPos = oneMarker.getMarkerPosition();
        const DynamicList<Pair<label>>& supportCells = oneMarker.getSupportCells();

        for(label suppInd=0; suppInd<supportCells.size(); suppInd++)
        {
            const Pair<label>& oneSuppCell = supportCells[suppInd];
            label proc = oneSuppCell.first();
            label cellInd = oneSuppCell.second();
            if(proc==Pstream::myProcNo())
            {
                vector cellCentre = cells[cellInd].centre(points,faces);
                scalar factor;
                if(modusFieldToMarker==markerMeshType::Uniform)
                    factor = oneMarker.deltaDirac(oneMarkerPos,cellCentre);
                else
                    factor = oneMarker.correctedDeltaDirac(oneMarkerPos,cellCentre);                
                factor *= oneMarker.getMarkerVolume();
                factor *= oneMarker.getMarkerWeight();
                fieldData[cellInd] += markerData[index]*factor;
            }
        }
    }
        
    // Read broadcasted marker weights
    /*
    std::unique_ptr<List<List<DynamicList<scalar>>>>& haloMarkerWeightsPtr = structure.get_haloMarkerWeights();
    if(!haloMarkerWeightsPtr)
        FatalErrorInFunction<<"Halo Marker Values not broadcasted before markerToField"<<exit(FatalError);
    const List<List<DynamicList<scalar>>>& haloMarkerWeights = *haloMarkerWeightsPtr;
    */
        
    // Broadcast marker field
    if(markerData.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Marker Data size mismatch"<< exit(FatalError);
    std::unique_ptr<List<List<DynamicList<T>>>> broadcastedHaloMarkerFieldPtr; 
    broadcastedHaloMarkerFieldPtr = structure.broadcastHaloMarkerField<T>(markerData);
    List<List<DynamicList<T>>> broadcastedHaloMarkerField = *broadcastedHaloMarkerFieldPtr;
        
    const LineStructure::GlobalHaloMarkers& gHM = structure.get_globalHaloMarkers();
    using LM=LagrangianMarker;
    const std::unordered_set<label>& neighbourProcs = structure.getNeighbourProcesses();
        
    for(label proc : neighbourProcs)
    {
        //const List<DynamicList<scalar>>& thisProcHaloMarkerWeights = haloMarkerWeights[proc];
        const List<DynamicList<T>>& thisProcHaloMarkeField = broadcastedHaloMarkerField[proc];
        label numHaloCells = gHM.size_haloCells(proc);
        if(numHaloCells!=thisProcHaloMarkeField.size())
            FatalErrorInFunction<<"Halo cell number mismatch"<<exit(FatalError);        
        for(label haloCellInd=0; haloCellInd<numHaloCells; haloCellInd++)
        {
            label numHaloCellMarkers = gHM.size_cellMarkers(proc,haloCellInd);
            //const DynamicList<scalar>& thisHaloCellMarkerWeights = thisProcHaloMarkerWeights[haloCellInd];
            const DynamicList<T>& thisHaloCellMarkeField = thisProcHaloMarkeField[haloCellInd];
            if(numHaloCellMarkers!=thisHaloCellMarkeField.size())
            {
                Pout<<"proc:"<<proc<<Foam::nl;
                Pout<<"haloCellInd:"<<haloCellInd<<Foam::nl;
                Pout<<"numHaloCells:"<<numHaloCells<<Foam::nl;
                Pout<<"thisProcHaloMarkeField.size():"<<thisProcHaloMarkeField.size()<<Foam::nl;
                Pout<<"numHaloCellMarkers:"<<numHaloCellMarkers<<Foam::nl;
                Pout<<"thisHaloCellMarkeField.size():"<<thisHaloCellMarkeField.size()<<Foam::nl;
                FatalErrorInFunction<<"Marker number mismatch"<<exit(FatalError);
            }
            for(label markerInd=0; markerInd<numHaloCellMarkers; markerInd++)
            {
                std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,FixedList<scalar,10>,label,scalar,scalar,scalar> haloMarkerData;
                haloMarkerData = gHM.getMarkerData(proc,haloCellInd,markerInd);

                vector position = std::get<0>(haloMarkerData);
                scalar volume = std::get<1>(haloMarkerData);
                //label index = std::get<2>(haloMarkerData);
                vector dilation = std::get<3>(haloMarkerData);
                FixedList<scalar,10>& markerKb = std::get<7>(haloMarkerData);
                scalar weight = gHM.getMarkerWeight(proc,haloCellInd,markerInd);
                T fieldValue = thisHaloCellMarkeField[markerInd];

                DynamicList<Pair<label>>& suppCellIndices =  std::get<4>(haloMarkerData);
                for(label suppInd=0; suppInd<suppCellIndices.size(); suppInd++)
                {
                    label proc = suppCellIndices[suppInd].first();
                    if(proc==Pstream::myProcNo())
                    {
                        label cellInd = suppCellIndices[suppInd].second();
                        vector cellCentre = cells[cellInd].centre(points,faces);
                     
                        scalar factor;
                        if(modusFieldToMarker==markerMeshType::Uniform)
                            factor = LM::deltaDirac
                            (
                                position,cellCentre,dilation
                            );
                        else
                            factor = LM::correctedDeltaDirac
                            (
                                position,cellCentre,dilation,markerKb
                            );
                        
                        factor *= volume;
                        factor *= weight;
                        fieldData[cellInd] += fieldValue*factor;
                    }
                }
            }
        }
    }
}

template<typename T>
T Foam::FieldMarkerStructureInteraction::integrateField
(
    const Field<T>& field
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();

    if(field.size()!=cells.size())
        FatalErrorInFunction<<"Mismatch in cell number"<<exit(FatalError);

    T sum = Foam::zero();
    for(label cellInd=0; cellInd<cells.size(); cellInd++)
    {
        sum += field[cellInd]*cells[cellInd].mag(points,faces);
    }
    
    return sum;
}

template<typename T> 
T Foam::FieldMarkerStructureInteraction::integrateMarkers
(
    const DynamicList<T>& markerVals
)
{
    const cellList& cells = mesh.cells();
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    
    if(markerVals.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in marker number"<<exit(FatalError);
    
    Field<T> meshField(cells.size());
    markerToField<T>(markerVals,meshField);
    return integrateField<T>(meshField);
}
