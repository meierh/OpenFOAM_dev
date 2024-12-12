template<typename T> void Foam::SensitivityInteraction::deriveParamFieldToMarker
(
    const GeometricField<T,fvPatchField,volMesh>& fieldData,
    DynamicList<T>& markerData,
    const Parameter& parameter
)
{
    std::function<T(const LagrangianMarker*,std::function<T(Pair<label>)>)> convolution;
    convolution = [this,&parameter](const LagrangianMarker* marker, std::function<T(Pair<label>)> valueFunction)
    {
        vector dX_dParam = structure.dXdParam(marker,parameter);
        
        std::function<scalar(vector,vector)> weightFunction;
        if(modusFieldToMarker==markerMeshType::Uniform)
        {
            weightFunction = [&marker,&dX_dParam] (vector X, vector x)
            {
                vector ddeltaDirac_dX = marker->ddeltaDirac_dX(X,x);
                return ddeltaDirac_dX & dX_dParam;
            };
        }
        else
        {
            weightFunction = [&marker,&dX_dParam] (vector X, vector x)
            {
                vector dcorrectedDeltaDirac_dX = marker->dcorrectedDeltaDirac_dX(X,x);
                return dcorrectedDeltaDirac_dX & dX_dParam;
            };
        }
        return marker->convolute<T>(weightFunction,valueFunction);
    };
    
    integrateFieldOverMarker(fieldData,markerData,convolution);
}

template<typename T> void Foam::SensitivityInteraction::deriveParamMarkerToField
(
    const DynamicList<T>& markerData,
    GeometricField<T,fvPatchField,volMesh>& fieldData,
    const Parameter& parameter
)
{
    fieldData = Foam::zero();
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    const std::vector<LagrangianMarker*>& markers = this->structure.getCollectedMarkers();
    
    if(markers.size()!=markerData.size())
        FatalErrorInFunction<<"Mismatch in markers.size():"<<markers.size()<<" and markerData.size():"<<markerData.size()<<exit(FatalError);
    
    // Spread marker field of own markers 
    for(label index=0; index<markerData.size(); index++)
    {
        const LagrangianMarker* oneMarkerPtr = markers[index];
        const LagrangianMarker& oneMarker = *(oneMarkerPtr);
        const vector& oneMarkerPos = oneMarker.getMarkerPosition();
        const DynamicList<Pair<label>>& supportCells = oneMarker.getSupportCells();
        
        vector dX_dParam = structure.dXdParam(oneMarkerPtr,parameter);

        for(label suppInd=0; suppInd<supportCells.size(); suppInd++)
        {
            const Pair<label>& oneSuppCell = supportCells[suppInd];
            label proc = oneSuppCell.first();
            label cellInd = oneSuppCell.second();
            if(proc==Pstream::myProcNo())
            {
                vector cellCentre = cells[cellInd].centre(points,faces);
                vector X = oneMarkerPos;
                vector x = cellCentre;
                scalar factor;
                if(modusFieldToMarker==markerMeshType::Uniform)
                {
                    vector ddeltaDirac_dX = oneMarker.ddeltaDirac_dX(X,x);
                    factor = ddeltaDirac_dX & dX_dParam;
                }
                else
                {
                    vector dcorrectedDeltaDirac_dX = oneMarker.dcorrectedDeltaDirac_dX(X,x);
                    factor = dcorrectedDeltaDirac_dX & dX_dParam;   
                }
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
        if(/*numHaloCells!=thisProcHaloMarkerWeights.size() ||*/ numHaloCells!=thisProcHaloMarkeField.size())
            FatalErrorInFunction<<"Halo cell number mismatch"<<exit(FatalError);        
        for(label haloCellInd=0; haloCellInd<numHaloCells; haloCellInd++)
        {
            label numHaloCellMarkers = gHM.size_cellMarkers(proc,haloCellInd);
            //const DynamicList<scalar>& thisHaloCellMarkerWeights = thisProcHaloMarkerWeights[haloCellInd];
            const DynamicList<T>& thisHaloCellMarkeField = thisProcHaloMarkeField[haloCellInd];
            if(/*numHaloCellMarkers!=thisHaloCellMarkerWeights.size() ||*/ numHaloCellMarkers!=thisHaloCellMarkeField.size())
                FatalErrorInFunction<<"Marker number mismatch"<<exit(FatalError);
            for(label markerInd=0; markerInd<numHaloCellMarkers; markerInd++)
            {
                std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,FixedList<scalar,10>,label,scalar,scalar,scalar> haloMarkerData;
                haloMarkerData = gHM.getMarkerData(proc,haloCellInd,markerInd);

                vector position = std::get<0>(haloMarkerData);
                scalar volume = std::get<1>(haloMarkerData);
                //label index = std::get<2>(haloMarkerData);
                vector dilation = std::get<3>(haloMarkerData);
                FixedList<scalar,10>& markerKb = std::get<7>(haloMarkerData);
                label rodNumber = std::get<8>(haloMarkerData);
                scalar param = std::get<9>(haloMarkerData);
                scalar angle = std::get<10>(haloMarkerData);
                scalar radiusFrac = std::get<11>(haloMarkerData);
                
                vector dX_dParam = structure.dXdParam(rodNumber,param,angle,radiusFrac,parameter);
                
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
                        {
                            vector ddeltaDirac_dX = LM::ddeltaDirac_dX(position,cellCentre,dilation);
                            factor = ddeltaDirac_dX & dX_dParam;
                        }
                        else
                        {
                            vector dcorrectedDeltaDirac_dX = LM::dcorrectedDeltaDirac_dX(position,cellCentre,dilation,markerKb);
                            factor = dcorrectedDeltaDirac_dX & dX_dParam;
                        }
                        
                        factor *= volume;
                        factor *= weight;
                        fieldData[cellInd] += fieldValue*factor;
                    }
                }
            }
        }
    }
}
