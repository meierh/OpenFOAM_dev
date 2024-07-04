#include "LineStructure.H"

Foam::FieldMarkerStructureInteraction::FieldMarkerStructureInteraction
(
    dynamicRefineFvMesh& mesh,
    LineStructure& structure,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
mesh(mesh),
structure(structure),
h(std::cbrt(mesh.cells()[0].mag(mesh.points(),mesh.faces()))),
modusFieldToMarker(modusFieldToMarker),
modusMarkerToField(modusMarkerToField)
{
    Info<<"Completed FieldMarkerStructureInteraction setup"<<Foam::endl;
}

template<typename T>
void Foam::FieldMarkerStructureInteraction::fieldToMarker
(
    const GeometricField<T,fvPatchField,volMesh>& fieldData,
    DynamicList<T>& markerData
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    Structure& structure = this->structure;
    
    std::unique_ptr<List<List<T>>> haloFieldsPtr = structure.broadcastHaloFields(fieldData);
    const List<List<T>>& haloFields = *haloFieldsPtr;

    markerData.resize(markers.size());
    std::function<T(std::tuple<bool,label,label>)> valueFunction = 
    [&fieldData,&haloFields,&structure](std::tuple<bool,label,label> cell)
    {
        const std::tuple<bool,label,label>& suppCellData = cell;
        if(std::get<0>(suppCellData))
        {
            label cellInd = std::get<2>(suppCellData);
            return fieldData[cellInd];
        }
        else
        {
            label proc = std::get<1>(suppCellData);
            label cellInd = std::get<2>(suppCellData);
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
}

template<typename T>
void Foam::FieldMarkerStructureInteraction::markerToField
(
    const DynamicList<T>& markerData,
    GeometricField<T,fvPatchField,volMesh>& fieldData
)
{
    fieldData = Foam::zero();
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    // Spread marker field of own markers 
    for(label index=0; index<markerData.size(); index++)
    {
        T value = Foam::zero();
        const LagrangianMarker& oneMarker = *(markers[index]);
        const vector& oneMarkerPos = oneMarker.getMarkerPosition();
        const DynamicList<std::tuple<bool,label,label>>& supportCells = oneMarker.getSupportCells();

        for(label suppInd=0; suppInd<supportCells.size(); suppInd++)
        {
            const std::tuple<bool,label,label>& oneSuppCell = supportCells[suppInd];
            label cellInd = std::get<2>(oneSuppCell);
            if(std::get<0>(oneSuppCell))
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
    std::unique_ptr<List<List<DynamicList<scalar>>>>& haloMarkerWeightsPtr = structure.get_haloMarkerWeights();
    if(!haloMarkerWeightsPtr)
        FatalErrorInFunction<<"Halo Marker Values not broadcasted before markerToField"<<exit(FatalError);
    const List<List<DynamicList<scalar>>>& haloMarkerWeights = *haloMarkerWeightsPtr;
    
    // Broadcast marker field
    if(markerData.size()!=markers.size())
        FatalErrorInFunction<<"Marker Data size mismatch"<< exit(FatalError);
    std::unique_ptr<List<List<DynamicList<T>>>> broadcastedHaloMarkerFieldPtr = structure.broadcastHaloMarkerFields(markerData);
    List<List<DynamicList<T>>> broadcastedHaloMarkerField = *broadcastedHaloMarkerFieldPtr;
    
    const LineStructure::GlobalHaloMarkers& gHM = structure.get_globalHaloMarkers();
    
    using LM=LagrangianMarker;
    const std::unordered_set<label>& neighbourProcs = structure.getNeighbourProcesses();
    for(label proc : neighbourProcs)
    {
        const List<DynamicList<scalar>>& thisProcHaloMarkerWeights = haloMarkerWeights[proc];
        const List<DynamicList<T>>& thisProcHaloMarkeField = broadcastedHaloMarkerField[proc];
        label numHaloCells = gHM.size_haloCells(proc);
        if(numHaloCells!=thisProcHaloMarkerWeights.size() || numHaloCells!=thisProcHaloMarkeField.size())
            FatalErrorInFunction<<"Halo cell number mismatch"<<exit(FatalError);        
        for(label haloCellInd=0; haloCellInd<numHaloCells; haloCellInd++)
        {
            label numHaloCellMarkers = gHM.size_cellMarkers(proc,haloCellInd);
            const DynamicList<scalar>& thisHaloCellMarkerWeights = thisProcHaloMarkerWeights[haloCellInd];
            const DynamicList<T>& thisHaloCellMarkeField = thisProcHaloMarkeField[haloCellInd];
            if(numHaloCellMarkers!=thisHaloCellMarkerWeights.size() || numHaloCellMarkers!=thisHaloCellMarkeField.size())
                FatalErrorInFunction<<"Marker number mismatch"<<exit(FatalError);
            for(label markerInd=0; markerInd<numHaloCellMarkers; markerInd++)
            {
                std::tuple<vector,scalar,label,vector,DynamicList<Pair<label>>,DynamicList<vector>,DynamicList<scalar>,List<scalar>> haloMarkerData;
                haloMarkerData = gHM.getMarkerData(proc,haloCellInd,markerInd);

                vector position = std::get<0>(haloMarkerData);
                scalar volume = std::get<1>(haloMarkerData);
                label index = std::get<2>(haloMarkerData);
                vector dilation = std::get<3>(haloMarkerData);
                FixedList<scalar,10>& markerKb = std::get<7>(haloMarkerData);
                scalar weight = thisHaloCellMarkerWeights[markerInd];
                T fieldValue = thisHaloCellMarkeField[markerInd];

                DynamicList<Pair<label>>& suppCellIndices =  std::get<4>(haloMarkerData);
                for(label suppInd=0; suppInd<suppCellIndices.size(); suppInd++)
                {
                    label proc = suppCellIndices[suppInd].first;
                    if(proc==Pstream::myProcNo())
                    {
                        label cellInd = suppCellIndices[suppInd].second;
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
        for(uint i=0;i<knots.size();i++)
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
