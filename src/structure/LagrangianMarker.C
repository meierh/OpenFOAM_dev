#include "LineStructure.H"

Foam::LagrangianMarker::LagrangianMarker
(
    const Structure& structure,
    const dynamicRefineFvMesh& mesh,
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* baseRod,
    const scalar markerParameter
):
structure(structure),
mesh(mesh),
rodNumber(rodNumber),
baseRod(baseRod),
markerParameter(markerParameter)
{
    evaluateMarker();
}

Foam::LagrangianMarker::LagrangianMarker
(    
    const Structure& structure,
    const dynamicRefineFvMesh& mesh,
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* baseRod
):
structure(structure),
mesh(mesh),
rodNumber(rodNumber),
baseRod(baseRod),
markerParameter(0)
{}

vector Foam::LagrangianMarker::getMarkerVelocity()
{
    return vector(0,0,0);
}

scalar Foam::LagrangianMarker::getMarkerTemperature()
{
    return 1;
}

scalar Foam::LagrangianMarker::getMarkerCellVolume()
{
    if(markerCell<0)
        return std::numeric_limits<scalar>::max();
    else
        return mesh.cells()[markerCell].mag(mesh.points(),mesh.faces());
}

std::pair<scalar,scalar> Foam::LagrangianMarker::getMarkerCellSpacing()
{
    if(markerCell<0)
    {
        //scalar maxLen = std::numeric_limits<scalar>::max();
        return {std::numeric_limits<scalar>::max(),std::numeric_limits<scalar>::min()};
    }
    else
    {
        const edgeList cellEdges = mesh.cells()[markerCell].edges(mesh.faces());
        scalar minLen = std::numeric_limits<scalar>::max();
        scalar maxLen = std::numeric_limits<scalar>::min();
        for(edge oneEdge : cellEdges)
        {
            minLen = std::min(oneEdge.mag(mesh.points()),minLen);
            maxLen = std::min(oneEdge.mag(mesh.points()),maxLen);
        }
        return {minLen,maxLen};
    }
}

scalar Foam::LagrangianMarker::getMarkerCellMinSpacing()
{
    return getMarkerCellSpacing().first;
}

std::string Foam::LagrangianMarker::to_string() const
{
    return std::to_string(markerParameter)+" cell:"+std::to_string(markerCell)+" pos:("+
           std::to_string(markerPosition[0])+","+std::to_string(markerPosition[1])+","+
           std::to_string(markerPosition[2])+")"+" vol:"+std::to_string(markerVolume);
}

void Foam::LagrangianMarker::evaluateMarker()
{
    markerPosition = LineStructure::evaluateRodPos(baseRod,markerParameter);
    markerCell = mesh.findCell(markerPosition);
    computeSupport();
    minMaxSupportWidth();
    dilationFactors();
}

void Foam::LagrangianMarker::computeSupport
(
    label iterations
)
{
    //Info<<"Compute support"<<Foam::endl;
    const cellList& cellList = mesh.cells();
    const faceList& facesList = mesh.faces();
    const labelList& owners = mesh.owner();
    const labelList& neighbours = mesh.neighbour();
    const pointField& points = mesh.points();
    //faceInd -> [{proc,cell}]
    const std::unordered_map<label,DynamicList<std::pair<label,label>>>& patchFaceToCell = structure.getPatchFaceToCellMap();

    struct FirstHash
    {
        label operator()(const std::pair<label, label> &p) const
        {
            return std::hash<label>{}(p.first);
        }
    };
    std::unordered_set<label> treatedCell;
    //pair{iteration,cellInd}
    std::unordered_set<std::pair<label,label>,FirstHash> ownSupportCells;
    //process -> {cellInd}
    std::unordered_map<label,std::unordered_set<label>> foreignHaloCells;
    std::unordered_set<label> directNeighborFaces;
    std::unordered_set<label> directNonNeighborFaces;
    //Info<<"-----------------------markerCell----------------:"<<markerCell<<Foam::endl;
    if(markerCell!=-1)
    {
        if(markerCell<0 || markerCell>=cellList.size())
            FatalErrorInFunction<<"Invalid cell index"<< exit(FatalError);
        ownSupportCells.insert({-1,markerCell});
        treatedCell.insert(markerCell);
        
        for(label iter=0; iter<iterations; iter++)
        {
            std::vector<std::pair<label,label>> addedOwnSupport;
            for(auto cellIter=ownSupportCells.begin(); cellIter!=ownSupportCells.end(); cellIter++)
            {
                if(cellIter->first==iter-1)
                {
                    label cellInd = cellIter->second;
                    //Info<<"iteration:"<<iter<<"  base:"<<cellInd<<Foam::endl;
                    const cell& thisCell = cellList[cellInd];
                    for(label faceInd : thisCell)
                    {
                        bool interProcessBound = false;
                        label neighborCell=-1;
                        if(owners[faceInd]==cellInd)
                        {
                            if(faceInd<neighbours.size())
                            {
                                neighborCell = neighbours[faceInd];
                            }
                            else
                            {
                                auto iter = patchFaceToCell.find(faceInd);
                                if(iter!=patchFaceToCell.end())
                                    interProcessBound = true;
                            }
                        }
                        else
                            neighborCell = owners[faceInd];
                        
                        //Info<<" faceInd:"<<faceInd<<"  |"<<neighborCell<<"/"<<interProcessBound<<Foam::endl;
                        
                        if(interProcessBound)
                        {
                            auto iter = patchFaceToCell.find(faceInd);
                            if(iter!=patchFaceToCell.end())
                            {
                                const DynamicList<std::pair<label,label>>& haloCells = iter->second;
                                for(const std::pair<label,label>& haloCell : haloCells)
                                {
                                    label process = haloCell.first;
                                    label cellInd = haloCell.second;
                                    foreignHaloCells[process].insert(cellInd);
                                }
                            }
                            directNeighborFaces.insert(faceInd);
                        }
                        else
                        {
                            if(neighborCell!=-1)
                            {
                                if(treatedCell.find(neighborCell)==treatedCell.end())
                                {
                                    addedOwnSupport.push_back({iter,neighborCell});
                                    treatedCell.insert(neighborCell);
                                }
                                directNeighborFaces.insert(faceInd);
                            }
                            else
                                directNonNeighborFaces.insert(faceInd);
                        }
                    }
                }
            }
            ownSupportCells.insert(addedOwnSupport.begin(),addedOwnSupport.end());
        }
    }
    //Info<<"ownSupportCells.size():"<<ownSupportCells.size()<<"  iterations:"<<iterations<<Foam::endl;
    this->supportCells.resize(0);
    for(auto iterCells=ownSupportCells.begin(); iterCells!=ownSupportCells.end(); iterCells++)
    {
        std::tuple<bool,label,label> cellData = {true,Pstream::myProcNo(),iterCells->second};
        this->supportCells.append(cellData);
    }
    for(auto iterProc=foreignHaloCells.begin(); iterProc!=foreignHaloCells.end(); iterProc++)
    {
        label processNo = iterProc->first;
        for(auto iterCell=iterProc->second.begin(); iterCell!=iterProc->second.end(); iterCell++)
        {
            label cellInd = *iterCell;
            std::tuple<bool,label,label> cellData = {false,processNo,cellInd};
            this->supportCells.append(cellData);
        }
    }
    
    //Info<<"markerParameter:"<<markerParameter<<"  markerPosition:"<<markerPosition<<"  "<<Foam::endl;
    
    List<List<bool>> directionsExist(3,List<bool>(2,false));
    for(label faceInd : directNeighborFaces)
    {
        const face& thisFace = facesList[faceInd];
        label owner = owners[faceInd];
        vector faceNormal = thisFace.normal(points);
        //Info<<"faceNormal:"<<faceNormal<<Foam::endl;
        faceNormal /= std::sqrt(faceNormal&faceNormal);
        if(owner!=markerCell)
        {
            if(markerCell==-1)
                FatalErrorInFunction<<"Invalid cell can not have neighbours"<<Foam::endl;
            faceNormal *= -1;
        }
        //Info<<"faceNormal:"<<faceNormal<<Foam::endl;
        label dir = -1;
        label posNeg = 1;
        if(std::abs(faceNormal[0]) >= 0.9)
        {
            dir = 0;
            if(faceNormal[dir]<0)
                posNeg = 0;
        }
        else if(std::abs(faceNormal[1]) >= 0.9)
        {
            dir = 1;
            if(faceNormal[dir]<0)
                posNeg = 0;
        }
        else if(std::abs(faceNormal[2]) >= 0.9)
        {
            dir = 2;
            if(faceNormal[dir]<0)
                posNeg = 0;
        }
        else
            FatalErrorInFunction<<"Faces have to be cartesian"<<exit(FatalError);

        directionsExist[dir][posNeg] = true;
    }
    //Info<<"directionsExist:"<<directionsExist<<Foam::endl;
    for(label dim=0; dim<3; dim++)
    {
        existingDims[dim] = directionsExist[dim][0] && directionsExist[dim][1];
    }
    
    /*
    for(auto neighTupl : supportCells)
        Info<<"     "<<std::get<2>(neighTupl)<<Foam::endl;
    */
    
}

void Foam::LagrangianMarker::computeSupportDimensions()
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    DynamicList<vector> centeredSupport;
    
    for(label i=0; i<supportCells.size(); i++)
    {
        const std::tuple<bool,label,label>& suppCellData = supportCells[i];
        vector cellCentre;
        if(std::get<0>(suppCellData))
            cellCentre = cells[std::get<2>(suppCellData)].centre(points,faces);
        else
        {
            label neighProcess = std::get<1>(suppCellData);
            const std::unordered_map<label,label>& neighborHaloCellToIndexMap = structure.getHaloCellToIndexMap(neighProcess);
            auto iter = neighborHaloCellToIndexMap.find(std::get<2>(suppCellData));
            if(iter==neighborHaloCellToIndexMap.end())
                FatalErrorInFunction<<"Halo cell does not exist"<<exit(FatalError);
            label index = iter->second;
            cellCentre = structure.getHaloCellList(neighProcess)[index].centre;
        }
        centeredSupport.append(cellCentre-markerPosition);
    }
    
    
}

void Foam::LagrangianMarker::minMaxSupportWidth()
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    vector minSpan;
    vector maxSpan;
    //Info<<"supportCells.size():"<<supportCells.size()<<Foam::endl;
    
    if(supportCells.size()>1)
    {
        scalar min = std::numeric_limits<scalar>::min();
        scalar max = std::numeric_limits<scalar>::max();
        minSpan = vector(max,max,max);
        maxSpan = vector(min,min,min);
        
        for(label i=0; i<supportCells.size(); i++)
        {
            const std::tuple<bool,label,label>& suppCellDataA = supportCells[i];
            vector cellCentreA;
            if(std::get<0>(suppCellDataA))
                cellCentreA = cells[std::get<2>(suppCellDataA)].centre(points,faces);
            else
            {
                label neighProcess = std::get<1>(suppCellDataA);
                const std::unordered_map<label,label>& neighborHaloCellToIndexMap = structure.getHaloCellToIndexMap(neighProcess);
                auto iter = neighborHaloCellToIndexMap.find(std::get<2>(suppCellDataA));
                if(iter==neighborHaloCellToIndexMap.end())
                    FatalErrorInFunction<<"Halo cell does not exist"<<exit(FatalError);
                label index = iter->second;
                cellCentreA = structure.getHaloCellList(neighProcess)[index].centre;
            }
            
            for(label j=0; j<supportCells.size(); j++)
            {
                if(i!=j)
                {
                    const std::tuple<bool,label,label>& suppCellDataB = supportCells[j];
                    vector cellCentreB;
                    if(std::get<0>(suppCellDataB))
                        cellCentreB = cells[std::get<2>(suppCellDataB)].centre(points,faces);
                    else
                    {
                        label neighProcess = std::get<1>(suppCellDataB);
                        const std::unordered_map<label,label>& neighborHaloCellToIndexMap = structure.getHaloCellToIndexMap(neighProcess);
                        auto iter = neighborHaloCellToIndexMap.find(std::get<2>(suppCellDataB));
                        if(iter==neighborHaloCellToIndexMap.end())
                            FatalErrorInFunction<<"Halo cell does not exist"<<exit(FatalError);
                        label index = iter->second;
                        cellCentreB = structure.getHaloCellList(neighProcess)[index].centre;
                    }
                                    
                    vector conn = cellCentreA-cellCentreB;
                    for(label d=0; d<3; d++)
                    {
                        conn[d] = std::abs(conn[d]);
                        minSpan[d] = std::min(minSpan[d],conn[d]);
                        maxSpan[d] = std::max(maxSpan[d],conn[d]);
                    }
                }
            }
        }
    }
    else
    {
        scalar min = structure.initialMeshSpacing;
        scalar max = structure.initialMeshSpacing;
        minSpan = vector(max,max,max);
        maxSpan = vector(min,min,min);
    }
    h_plus = maxSpan;
    h_minus = minSpan;
    //Info<<"maxSpan:"<<maxSpan<<Foam::endl;
    //Info<<"minSpan:"<<minSpan<<Foam::endl;
}

void Foam::LagrangianMarker::dilationFactors()
{
    vector minSpan = h_minus;
    vector maxSpan = h_plus;
    //Info<<"     h_minus:"<<h_minus<<Foam::endl;
    //Info<<"     h_plus:"<<h_plus<<Foam::endl;
    
    //const cellList& cells = mesh.cells();
    //const faceList& faces = mesh.faces();
    //const pointField& points = mesh.points();
    
    scalar eps = 0.1*std::sqrt(minSpan&minSpan);
    dilation = 5.0/6.0 * maxSpan + 1.0/6.0 * minSpan + vector(eps,eps,eps);
    //Info<<" dilation:"<<dilation<<Foam::endl;
}

void Foam::LagrangianMarker::getCellData
(
    const std::tuple<bool,label,label>& cell,
    vector& cellCentre,
    scalar& cellVolume
) const
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    const std::tuple<bool,label,label>& suppCellData = cell;
    if(std::get<0>(suppCellData))
    {
        label cellInd = std::get<2>(suppCellData);
        cellCentre = cells[cellInd].centre(points,faces);
        cellVolume = cells[cellInd].mag(points,faces);
    }
    else
    {
        label neighProcess = std::get<1>(suppCellData);
        const std::unordered_map<label,label>& neighborHaloCellToIndexMap = structure.getHaloCellToIndexMap(neighProcess);
        auto iter = neighborHaloCellToIndexMap.find(std::get<2>(suppCellData));
        if(iter==neighborHaloCellToIndexMap.end())
            FatalErrorInFunction<<"Halo cell does not exist"<<exit(FatalError);
        label index = iter->second;
        cellCentre = structure.getHaloCellList(neighProcess)[index].centre;
        cellVolume = structure.getHaloCellList(neighProcess)[index].volume;
    }
}

scalar Foam::LagrangianMarker::computeMoment
(
    vector indices,
    std::function<scalar(vector,vector)> deltaFunction
) const
{
    const LagrangianMarker& marker = *this;
    std::function<scalar(std::tuple<bool,label,label>)> valueFunction = 
    [marker,indices](std::tuple<bool,label,label> cell)
    {
        vector X = marker.getMarkerPosition();
        vector cellCentre;
        scalar cellVolume;
        marker.getCellData(cell,cellCentre,cellVolume);
        vector x = cellCentre;
        vector conn = x-X;
        vector coeff(1,1,1);
        for(label dim=0; dim<3; dim++)
        {
            for(label index=0; index<indices[dim]; index++)
            {
                coeff[dim]*=conn[dim];
            }
        }
        return coeff[0]*coeff[1]*coeff[2];
    };
    return convolute<scalar>(deltaFunction,valueFunction);
}

scalar Foam::LagrangianMarker::computeCorrectedMoment
(
    vector indices
) const
{
    const LagrangianMarker& marker = *this;
    std::function<scalar(vector,vector)> deltaFunction = 
    [marker](vector X, vector x)
    {
        return marker.correctedDeltaDirac(X,x);
    };
    return computeMoment(indices,deltaFunction);
}

scalar Foam::LagrangianMarker::computeMoment
(
    vector indices
) const
{
    const LagrangianMarker& marker = *this;
    std::function<scalar(vector,vector)> deltaFunction = 
    [marker](vector X, vector x)
    {
        return marker.deltaDirac(X,x);
    };
    return computeMoment(indices,deltaFunction);
}

/*
scalar Foam::LagrangianMarker::computeMoment
(
    vector indices,
    std::function<scalar(vector,vector)> deltaFunction
) const
{
    vector X = getMarkerPosition();
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    const DynamicList<std::tuple<bool,label,label>>& supportCells = getSupportCells();
    
    scalar moment = 0;
    for(label i=0; i<supportCells.size(); i++)
    {
        const std::tuple<bool,label,label>& suppCellData = supportCells[i];
        vector cellCentre;
        scalar cellVolume;
        if(std::get<0>(suppCellData))
        {
            label cellInd = std::get<2>(suppCellData);
            cellCentre = cells[cellInd].centre(points,faces);
            cellVolume = cells[cellInd].mag(points,faces);
        }
        else
        {
            label neighProcess = std::get<1>(suppCellData);
            const std::unordered_map<label,label>& neighborHaloCellToIndexMap = structure.getHaloCellToIndexMap(neighProcess);
            auto iter = neighborHaloCellToIndexMap.find(std::get<2>(suppCellData));
            if(iter==neighborHaloCellToIndexMap.end())
                FatalErrorInFunction<<"Halo cell does not exist"<<exit(FatalError);
            label index = iter->second;
            cellCentre = structure.getHaloCellList(neighProcess)[index].centre;
            cellVolume = structure.getHaloCellList(neighProcess)[index].volume;
        }
        
        vector x = cellCentre;
        vector conn = x-X;
        vector coeff(1,1,1);
        for(label dim=0; dim<3; dim++)
        {
            for(label index=0; index<indices[dim]; index++)
            {
                coeff[dim]*=conn[dim];
            }
        }
        scalar dirac = deltaFunction(X,x);
        //Info<<conn<<"|"<<dirac<<Foam::endl;
        moment += coeff[0]*coeff[1]*coeff[2]*dirac*cellVolume;
    }
    return moment;
};
*/


/*
scalar Foam::LagrangianMarker::computeCorrectedMoment
(
    vector indices
) const
{
    vector X = getMarkerPosition();
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    const DynamicList<std::tuple<bool,label,label>>& supportCells = getSupportCells();
    
    scalar moment = 0;
    for(label i=0; i<supportCells.size(); i++)
    {
        const std::tuple<bool,label,label>& suppCellData = supportCells[i];
        vector cellCentre;
        scalar cellVolume;
        if(std::get<0>(suppCellData))
        {
            label cellInd = std::get<2>(suppCellData);
            cellCentre = cells[cellInd].centre(points,faces);
            cellVolume = cells[cellInd].mag(points,faces);
        }
        else
        {
            label neighProcess = std::get<1>(suppCellData);
            const std::unordered_map<label,label>& neighborHaloCellToIndexMap = structure.getHaloCellToIndexMap(neighProcess);
            auto iter = neighborHaloCellToIndexMap.find(std::get<2>(suppCellData));
            if(iter==neighborHaloCellToIndexMap.end())
                FatalErrorInFunction<<"Halo cell does not exist"<<exit(FatalError);
            label index = iter->second;
            cellCentre = structure.getHaloCellList(neighProcess)[index].centre;
            cellVolume = structure.getHaloCellList(neighProcess)[index].volume;
        }
        
        vector x = cellCentre;
        vector conn = x-X;
        vector coeff(1,1,1);
        for(label dim=0; dim<3; dim++)
        {
            for(label index=0; index<indices[dim]; index++)
            {
                coeff[dim]*=coeff[dim];
            }
        }
        scalar dirac = correctedDeltaDirac(X,x,getDilation(),b);
        moment += coeff[0]*coeff[1]*coeff[2]*dirac*cellVolume;
    }
    return moment;
};
*/

std::unique_ptr<gismo::gsMatrix<scalar>> Foam::LagrangianMarker::computeMomentMatrix3D() const
{
    std::array<vector,10> gen3DIndices = 
    {
        vector(0,0,0),
        vector(1,0,0),vector(0,1,0),vector(0,0,1),
        vector(1,1,0),vector(0,1,1),vector(1,0,1),vector(2,0,0),vector(0,2,0),vector(0,0,2)        
    };
    
    List<List<vector>> indices(10,List<vector>(10));
    for(uint i=0; i<gen3DIndices.size(); i++)
    {
        vector indi = gen3DIndices[i];
        for(uint j=0; j<gen3DIndices.size(); j++)
        {
            vector indj = gen3DIndices[j];
            vector index = indi+indj;
            indices[i][j] = index;
        }
    }

    auto moments3DPtr = std::unique_ptr<gismo::gsMatrix<scalar>>(new gismo::gsMatrix<scalar>(10,10));
    gismo::gsMatrix<scalar>& moments3D = *moments3DPtr;
    for(uint i=0; i<gen3DIndices.size(); i++)
    {
        for(uint j=0; j<gen3DIndices.size(); j++)
        {
            moments3D(i,j) = computeMoment(indices[i][j]);
        }
    }

    return moments3DPtr;
}

std::unique_ptr<gismo::gsMatrix<scalar>> Foam::LagrangianMarker::computeMomentMatrix2D
(
    std::array<label,2> dim
) const
{
    std::array<std::array<label,2>,6> gen2DIndices;
    gen2DIndices[0] = {0,0}; gen2DIndices[1] = {1,0}; gen2DIndices[2] = {0,1};
    gen2DIndices[3] = {1,1}; gen2DIndices[4] = {2,0}; gen2DIndices[5] = {0,2};
    
    if(dim[0]!=0 && dim[0]!=1 && dim[0]!=2) 
        FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    if(dim[1]!=0 && dim[1]!=1 && dim[1]!=2)
        FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    if(dim[0]==dim[1] || dim[0]>dim[1])
        FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);

    std::array<vector,6> gen3DIndices;
    for(uint i=0; i<gen3DIndices.size(); i++)
    {
        std::array<label,2> singleEntry = gen2DIndices[i];
        vector indices3D(0,0,0);
        indices3D[dim[0]] = singleEntry[0];
        indices3D[dim[1]] = singleEntry[1];
        gen3DIndices[i] = indices3D;
    }

    List<List<vector>> indices(6,List<vector>(6));
    for(uint i=0; i<gen3DIndices.size(); i++)
    {
        vector indi = gen3DIndices[i];
        for(uint j=0; j<gen3DIndices.size(); j++)
        {
            vector indj = gen3DIndices[j];
            vector index = indi+indj;
            indices[i][j] = index;
        }
    }
    
    auto moments2DPtr = std::unique_ptr<gismo::gsMatrix<scalar>>(new gismo::gsMatrix<scalar>(6,6));
    gismo::gsMatrix<scalar>& moments2D = *moments2DPtr;
    for(uint i=0; i<gen3DIndices.size(); i++)
    {
        for(uint j=0; j<gen3DIndices.size(); j++)
        {
            moments2D(i,j) = computeMoment(indices[i][j]);
        }
    }

    return moments2DPtr;
}

std::unique_ptr<gismo::gsMatrix<scalar>> Foam::LagrangianMarker::computeMomentMatrix1D
(
    label dim
) const
{
    if(dim!=0 && dim!=1 && dim!=2)
        FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    
    auto moments1DPtr = std::unique_ptr<gismo::gsMatrix<scalar>>(new gismo::gsMatrix<scalar>(3,3));
    gismo::gsMatrix<scalar>& moments1D = *moments1DPtr;
    
    //Initial
    moments1D(0,0) = computeMoment(vector(0,0,0));
    
    vector indices1(0,0,0);
    indices1[dim] = 1;    
    vector indices2(0,0,0);
    indices2[dim] = 2;    
    vector indices3(0,0,0);
    indices3[dim] = 3;    
    vector indices4(0,0,0);
    indices4[dim] = 4;
           
    //Diagonal
    moments1D(1,1) = computeMoment(indices2);
    moments1D(2,2) = computeMoment(indices4);
    
    //Lower / Upper triangle
    moments1D(0,1) = moments1D(1,0) = computeMoment(indices1);
    moments1D(0,2) = moments1D(2,0) = computeMoment(indices2);
    moments1D(1,2) = moments1D(2,1) = computeMoment(indices3);
    
    return moments1DPtr;
}

std::unique_ptr<std::array<scalar,10>> Foam::LagrangianMarker::rescalingDiagonal3D() const
{
    auto diagPtr = std::unique_ptr<std::array<scalar,10>>(new std::array<scalar,10>());
    std::array<scalar,10>& diag = *diagPtr;
    vector dilation = getDilation();
    
    diag[0] = 1;
    diag[1] = 1/dilation[0];
    diag[2] = 1/dilation[1];
    diag[3] = 1/dilation[2];
    diag[4] = 1/(dilation[0]*dilation[1]);
    diag[5] = 1/(dilation[1]*dilation[2]);
    diag[6] = 1/(dilation[0]*dilation[2]);
    diag[7] = 1/(dilation[0]*dilation[0]);
    diag[8] = 1/(dilation[1]*dilation[1]);
    diag[9] = 1/(dilation[2]*dilation[2]);
    return diagPtr;
}

std::unique_ptr<std::array<scalar,6>> Foam::LagrangianMarker::rescalingDiagonal2D
(
    label normalDim
) const
{
    if(normalDim!=0 && normalDim!=1 && normalDim!=2)
        FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    
    auto diagPtr = std::unique_ptr<std::array<scalar,6>>(new std::array<scalar,6>());
    std::array<scalar,6>& diag = *diagPtr;
    vector dilation = getDilation();
    
    diag[0] = 1;
    
    if(normalDim==0)
    {
        diag[1] = 1/dilation[1];
        diag[2] = 1/dilation[2];
        diag[3] = 1/(dilation[1]*dilation[2]);
        diag[4] = 1/(dilation[1]*dilation[1]);
        diag[5] = 1/(dilation[2]*dilation[2]);
    }
    else if(normalDim==1)
    {
        diag[1] = 1/dilation[0];
        diag[2] = 1/dilation[2];
        diag[3] = 1/(dilation[0]*dilation[2]);
        diag[4] = 1/(dilation[0]*dilation[0]);
        diag[5] = 1/(dilation[2]*dilation[2]);
    }
    else
    {
        diag[1] = 1/dilation[0];
        diag[2] = 1/dilation[1];
        diag[3] = 1/(dilation[0]*dilation[1]);
        diag[4] = 1/(dilation[0]*dilation[0]);
        diag[5] = 1/(dilation[1]*dilation[1]);
    }
    return diagPtr;
}

std::unique_ptr<std::array<scalar,3>> Foam::LagrangianMarker::rescalingDiagonal1D
(
    label dim
) const
{
    if(dim!=0 && dim!=1 && dim!=2)
        FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    
    auto diagPtr = std::unique_ptr<std::array<scalar,3>>(new std::array<scalar,3>());
    std::array<scalar,3>& diag = *diagPtr;
    vector dilation = getDilation();
    
    diag[0] = 1;
    diag[1] = 1/dilation[dim];
    diag[2] = 1/(dilation[dim]*dilation[dim]);
    return diagPtr;
}

void Foam::LagrangianMarker::computeCorrectionWeights()
{
    label numberDims = 0;
    DynamicList<label> listDims;
    for(label dim=0; dim<existingDims.size(); dim++)
    {
        if(existingDims[dim])
        {
            numberDims++;
            listDims.append(dim);
        }
    }
    WeightSystemSolve solStrategy = WeightSystemSolve::Raw;
    
    for(label i=0; i<10; i++)
        b[i] = 0;
        
    if(numberDims==0)
    /*
        b0 +
        0    + 0      + 0 +
        0    + 0      + 0 +
        0    + 0      + 0
    */
    {
        scalar m000 = computeMoment(vector(0,0,0));
        b[0] = 1.0/m000;
    }
    else if(numberDims==1)
    {        
        if(listDims[0]==0)
        /*
            b0 +
            (x-s)b1    + 0      + 0 +
            0          + 0      + 0 +
            (x-s)²b7   + 0      + 0
        */
        {
            gismo::gsMatrix<scalar> x;
            solveWeights<3>(rescalingDiagonal1D(0),computeMomentMatrix1D(0),x,solStrategy);
            b[0] = x(0,0);
            b[1] = x(1,0);
            b[7] = x(2,0);            
        }
        else if(listDims[0]==1)
        /*
            b0 +
            0  + (y-t)b2      + 0 +
            0  + 0            + 0 +
            0  + (y-t)²b8     + 0
        */
        {
            gismo::gsMatrix<scalar> x;
            solveWeights<3>(rescalingDiagonal1D(1),computeMomentMatrix1D(1),x,solStrategy);
            b[0] = x(0,0);
            b[2] = x(1,0);
            b[8] = x(2,0);
        }
        else if(listDims[0]==2)
        /*
            b0 +
            0  + 0  + (z-v)b3 +
            0  + 0  +      0  +
            0  + 0  + (z-v)²b9
        */
        {
            gismo::gsMatrix<scalar> x;
            solveWeights<3>(rescalingDiagonal1D(2),computeMomentMatrix1D(2),x,solStrategy);
            b[0] = x(0,0);
            b[3] = x(1,0);
            b[9] = x(2,0);
        }
        else
            FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    }
    else if(numberDims==2)
    {        
        if(listDims[0]==0 && listDims[1]==1)
        /*
            b0 +
            (x-s)b1      + (y-t)b2      + 0 +
            (x-s)(y-t)b4 + 0            + 0 +
            (x-s)²b7     + (y-t)²b8     + 0
        */
        {
            gismo::gsMatrix<scalar> x;
            solveWeights<6>(rescalingDiagonal2D(2),computeMomentMatrix2D({0,1}),x,solStrategy);
            b[0] = x(0,0);
            b[1] = x(1,0);
            b[2] = x(2,0);
            b[4] = x(3,0);
            b[7] = x(4,0);
            b[8] = x(5,0);
        }
        else if(listDims[0]==1 && listDims[1]==2)
        /*
            b0 +
            0       + (y-t)b2      + (z-v)b3 +
            0       + (y-t)(z-v)b5 + 0 +
            0       + (y-t)²b8     + (z-v)²b9
        */
        {
            gismo::gsMatrix<scalar> x;
            solveWeights<6>(rescalingDiagonal2D(0),computeMomentMatrix2D({1,2}),x,solStrategy);           
            b[0] = x(0,0);
            b[2] = x(1,0);
            b[3] = x(2,0);
            b[5] = x(3,0);
            b[8] = x(4,0);
            b[9] = x(5,0);
        }
        else if(listDims[0]==0 && listDims[1]==2)        
        /*
            b0 +
            (x-s)b1      + 0    + (z-v)b3 +
            0            + 0    + (z-v)(x-s)b6 +
            (x-s)²b7     + 0    + (z-v)²b9
        */
        {
            gismo::gsMatrix<scalar> x;
            solveWeights<6>(rescalingDiagonal2D(1),computeMomentMatrix2D({0,2}),x,solStrategy);
            b[0] = x(0,0);
            b[1] = x(1,0);
            b[3] = x(2,0);
            b[6] = x(3,0);
            b[7] = x(4,0);
            b[9] = x(5,0);
        }
        else
            FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    }
    else if(numberDims==3)
    /*
        b0 +
        (x-s)b1      + (y-t)b2      + (z-v)b3 +
        (x-s)(y-t)b4 + (y-t)(z-v)b5 + (z-v)(x-s)b6 +
        (x-s)²b7     + (y-t)²b8     + (z-v)²b9
    */
    {
        gismo::gsMatrix<scalar> x;
        solveWeights<10>(rescalingDiagonal3D(),computeMomentMatrix3D(),x,solStrategy);
        for(label i=0; i<10; i++)
            b[i] = x(i,0);
    }
}

scalar Foam::LagrangianMarker::deltaDirac
(
    vector X,
    vector x,
    scalar h
)
{
    return deltaDirac(X,x,vector(h,h,h));
}

scalar Foam::LagrangianMarker::deltaDirac
(
    vector X,
    vector x,
    vector h
)
{
    vector sigma_d;
    for(label dim=0; dim<3; dim++)
    {
        scalar X_i_x_i = X[dim]-x[dim];
        scalar r = X_i_x_i / h[dim];
        sigma_d[dim] = phiFunction(r);
    }
    scalar deltaDir = sigma_d[0]*sigma_d[1]*sigma_d[2];
    deltaDir /= (h[0]*h[1]*h[2]);
    //Info<<h<<"  X:"<<X<<"  x:"<<x<<" -> "<<deltaDir<<Foam::endl;
    return deltaDir;
}

scalar Foam::LagrangianMarker::correctedDeltaDirac
(
    vector X,
    vector x
) const
{
    return correctedDeltaDirac(X,x,dilation,b);
}

scalar Foam::LagrangianMarker::deltaDirac
(
    vector X,
    vector x
) const
{
    return deltaDirac(X,x,dilation);
}

scalar Foam::LagrangianMarker::correctedDeltaDirac
(
    vector X,
    vector x,
    scalar h,
    const std::array<scalar,10>& b
)
{
    return correctedDeltaDirac(X,x,vector(h,h,h),b);
}

scalar Foam::LagrangianMarker::correctedDeltaDirac
(
    vector X,
    vector x,
    vector h,
    const std::array<scalar,10>& b
)
{
    vector conn = x-X;
    scalar correctionFactor =   b[0] + 
                                conn[0]*b[1] + conn[1]*b[2] + conn[2]*b[3] +
                                conn[0]*conn[1]*b[4] + conn[1]*conn[2]*b[5] + conn[2]*conn[0]*b[6] +
                                conn[0]*conn[0]*b[7] + conn[1]*conn[1]*b[8] + conn[2]*conn[2]*b[9];
    return correctionFactor*deltaDirac(X,x,h);
}

scalar Foam::LagrangianMarker::phiFunction
(
    scalar r
)
{
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

std::unique_ptr<gismo::gsMatrix<scalar>> Foam::LagrangianMarker::computeCorrectedMomentMatrix() const
{
    auto moments3DPtr = std::unique_ptr<gismo::gsMatrix<scalar>>(new gismo::gsMatrix<scalar>(10,10));
    gismo::gsMatrix<scalar>& moments3D = *moments3DPtr;
    
    //Diagonal
    moments3D(0,0) = computeCorrectedMoment(vector(0,0,0));
    moments3D(1,1) = computeCorrectedMoment(vector(2,0,0));
    moments3D(2,2) = computeCorrectedMoment(vector(0,2,0));
    moments3D(3,3) = computeCorrectedMoment(vector(0,0,2));
    moments3D(4,4) = computeCorrectedMoment(vector(2,2,0));
    moments3D(5,5) = computeCorrectedMoment(vector(0,2,2));
    moments3D(6,6) = computeCorrectedMoment(vector(2,0,2));
    moments3D(7,7) = computeCorrectedMoment(vector(4,0,0));
    moments3D(8,8) = computeCorrectedMoment(vector(0,4,0));
    moments3D(9,9) = computeCorrectedMoment(vector(0,0,4));
    
    //Row / Column 0
    moments3D(0,1) = moments3D(1,0) = computeCorrectedMoment(vector(1,0,0));
    moments3D(0,2) = moments3D(2,0) = computeCorrectedMoment(vector(0,1,0));
    moments3D(0,3) = moments3D(3,0) = computeCorrectedMoment(vector(0,0,1));
    moments3D(0,4) = moments3D(4,0) = computeCorrectedMoment(vector(1,1,0));
    moments3D(0,5) = moments3D(5,0) = computeCorrectedMoment(vector(0,1,1));
    moments3D(0,6) = moments3D(6,0) = computeCorrectedMoment(vector(1,0,1));
    moments3D(0,7) = moments3D(7,0) = computeCorrectedMoment(vector(2,0,0));
    moments3D(0,8) = moments3D(8,0) = computeCorrectedMoment(vector(0,2,0));
    moments3D(0,9) = moments3D(9,0) = computeCorrectedMoment(vector(0,0,2));

    //Row / Column 1
    moments3D(1,2) = moments3D(2,1) = computeCorrectedMoment(vector(1,1,0));
    moments3D(1,3) = moments3D(3,1) = computeCorrectedMoment(vector(1,0,1));
    moments3D(1,4) = moments3D(4,1) = computeCorrectedMoment(vector(2,1,0));
    moments3D(1,5) = moments3D(5,1) = computeCorrectedMoment(vector(1,1,1));
    moments3D(1,6) = moments3D(6,1) = computeCorrectedMoment(vector(2,0,1));
    moments3D(1,7) = moments3D(7,1) = computeCorrectedMoment(vector(3,0,0));
    moments3D(1,8) = moments3D(8,1) = computeCorrectedMoment(vector(1,2,0));
    moments3D(1,9) = moments3D(9,1) = computeCorrectedMoment(vector(1,0,2));
    
    //Row / Column 2
    moments3D(2,3) = moments3D(3,2) = computeCorrectedMoment(vector(0,1,1));
    moments3D(2,4) = moments3D(4,2) = computeCorrectedMoment(vector(1,2,0));
    moments3D(2,5) = moments3D(5,2) = computeCorrectedMoment(vector(0,2,1));
    moments3D(2,6) = moments3D(6,2) = computeCorrectedMoment(vector(1,1,1));
    moments3D(2,7) = moments3D(7,2) = computeCorrectedMoment(vector(2,1,0));
    moments3D(2,8) = moments3D(8,2) = computeCorrectedMoment(vector(0,3,0));
    moments3D(2,9) = moments3D(9,2) = computeCorrectedMoment(vector(0,1,2));
    
    //Row / Column 3
    moments3D(3,4) = moments3D(4,3) = computeCorrectedMoment(vector(1,1,1));
    moments3D(3,5) = moments3D(5,3) = computeCorrectedMoment(vector(0,1,2));
    moments3D(3,6) = moments3D(6,3) = computeCorrectedMoment(vector(1,0,2));
    moments3D(3,7) = moments3D(7,3) = computeCorrectedMoment(vector(2,0,1));
    moments3D(3,8) = moments3D(8,3) = computeCorrectedMoment(vector(0,2,1));
    moments3D(3,9) = moments3D(9,3) = computeCorrectedMoment(vector(0,0,3));

    //Row / Column 4
    moments3D(3,4) = moments3D(4,3) = computeCorrectedMoment(vector(1,1,1));
    moments3D(3,5) = moments3D(5,3) = computeCorrectedMoment(vector(0,1,2));
    moments3D(3,6) = moments3D(6,3) = computeCorrectedMoment(vector(1,0,2));
    moments3D(3,7) = moments3D(7,3) = computeCorrectedMoment(vector(2,0,1));
    moments3D(3,8) = moments3D(8,3) = computeCorrectedMoment(vector(0,2,1));
    moments3D(3,9) = moments3D(9,3) = computeCorrectedMoment(vector(0,0,3));

    //Row / Column 5
    moments3D(4,5) = moments3D(5,4) = computeCorrectedMoment(vector(1,2,1));
    moments3D(4,6) = moments3D(6,4) = computeCorrectedMoment(vector(2,1,1));
    moments3D(4,7) = moments3D(7,4) = computeCorrectedMoment(vector(3,1,0));
    moments3D(4,8) = moments3D(8,4) = computeCorrectedMoment(vector(1,3,0));
    moments3D(4,9) = moments3D(9,4) = computeCorrectedMoment(vector(1,1,2));
    
    //Row / Column 6
    moments3D(5,6) = moments3D(6,5) = computeCorrectedMoment(vector(1,1,2));
    moments3D(5,7) = moments3D(7,5) = computeCorrectedMoment(vector(2,1,1));
    moments3D(5,8) = moments3D(8,5) = computeCorrectedMoment(vector(0,3,1));
    moments3D(5,9) = moments3D(9,5) = computeCorrectedMoment(vector(0,1,3));
    
    //Row / Column 7
    moments3D(6,7) = moments3D(7,6) = computeCorrectedMoment(vector(3,0,1));
    moments3D(6,8) = moments3D(8,6) = computeCorrectedMoment(vector(1,2,1));
    moments3D(6,9) = moments3D(9,6) = computeCorrectedMoment(vector(1,0,3));
    
    //Row / Column 8
    moments3D(7,8) = moments3D(8,7) = computeCorrectedMoment(vector(2,2,0));
    moments3D(7,9) = moments3D(9,7) = computeCorrectedMoment(vector(2,0,2));
    
    //Row / Column 9
    moments3D(8,9) = moments3D(9,8) = computeCorrectedMoment(vector(0,2,2));
    
    return moments3DPtr;
}
