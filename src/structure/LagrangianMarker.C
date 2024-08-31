#include "LineStructure.H"

Foam::LagrangianMarker::LagrangianMarker
(
    const Structure& structure,
    const fvMesh& mesh,
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
    const fvMesh& mesh,
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* baseRod
):
structure(structure),
mesh(mesh),
rodNumber(rodNumber),
baseRod(baseRod),
markerParameter(0)
{}

Foam::vector Foam::LagrangianMarker::getMarkerVelocity()
{
    return vector(0,0,0);
}

Foam::scalar Foam::LagrangianMarker::getMarkerTemperature()
{
    return 1;
}

Foam::scalar Foam::LagrangianMarker::getMarkerCellVolume()
{
    if(markerCell<0)
        return std::numeric_limits<scalar>::max();
    else
        return mesh.cells()[markerCell].mag(mesh.points(),mesh.faces());
}

std::pair<Foam::scalar,Foam::scalar> Foam::LagrangianMarker::getMarkerCellSpacing()
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

Foam::scalar Foam::LagrangianMarker::getMarkerCellMinSpacing()
{
    return getMarkerCellSpacing().first;
}

std::string Foam::LagrangianMarker::to_string() const
{
    return std::to_string(markerParameter)+" cell:"+std::to_string(markerCell)+" pos:("+
           std::to_string(markerPosition[0])+","+std::to_string(markerPosition[1])+","+
           std::to_string(markerPosition[2])+")"+" vol:"+std::to_string(markerVolume)+" dilation:("+
           std::to_string(dilation[0])+","+std::to_string(dilation[1])+","+
           std::to_string(dilation[2])+")"+ " existingDims:("+std::to_string(existingDims[0])+","+std::to_string(existingDims[1])+","+std::to_string(existingDims[2])+")";
}

void Foam::LagrangianMarker::total_print() const
{
    const cellList& cellList = mesh.cells();
    const faceList& facesList = mesh.faces();
    const pointField& points = mesh.points();
    
    Info<<"h_plus:"<<h_plus<<Foam::endl;
    Info<<"h_minus:"<<h_minus<<Foam::endl;
    Info<<"dilation:"<<dilation<<Foam::endl;
    
    Pair<vector> h = minMaxNeighbourWidth(directSupport);
    Info<<"h:"<<h<<Foam::endl;
    
    auto printCell = [&](Pair<label> cell)
    {
        Info<<"cell:"<<cell.second();
        if(cell.first()!=Pstream::myProcNo())
            Info<<Foam::endl;
        else
        {
            vector centre = cellList[cell.second()].centre(points,facesList);
            vector dist = centre - markerPosition;
            Info<<" : "<<"dist:"<<dist<<" -- "<<deltaDirac(markerPosition,centre)<<Foam::endl;
        }
    };
    
    Pout<<"\t directSupport:"<<Foam::endl;
    for(auto cell : directSupport)
    {
        printCell(cell);
    }
    Pout<<"\t fullSupport:"<<Foam::endl;
    for(auto cell : fullSupport)
    {
        printCell(cell);
    }
}

void Foam::LagrangianMarker::evaluateMarker()
{
    markerPosition = LineStructure::evaluateRodPos(baseRod,markerParameter);
    markerCell = mesh.findCell(markerPosition);
    computeSupport();
    Pair<vector> h = minMaxNeighbourWidth(directSupport);
    h_plus = h.first();
    h_minus = h.second();
    dilation = dilationFactors(h);
    checkDirectSupport();
    reduceSupport();
}

void Foam::LagrangianMarker::computeSupport
(
    label iterations
)
{
    const cellList& cellList = mesh.cells();
    const faceList& facesList = mesh.faces();
    const labelList& owners = mesh.owner();
    const pointField& points = mesh.points();
    
    const List<List<Pair<label>>>& localMeshGraph = structure.getMeshGraph();
    std::unordered_set<Pair<label>,foamPairHash<label>> direct;
    std::unordered_set<Pair<label>,foamPairHash<label>> full;
    if(markerCell!=-1)
    {
        if(markerCell<0 || markerCell>=cellList.size())
            FatalErrorInFunction<<"Invalid cell index"<< exit(FatalError);
        
        direct.insert({Pstream::myProcNo(),markerCell});
        full.insert({Pstream::myProcNo(),markerCell});
        DynamicList<Pair<label>> frontNodes;
        for(const Pair<label>& edge : localMeshGraph[markerCell])
        {
            if(edge.first()!=-1 && edge.second()!=-1)
            {
                direct.insert({edge.first(),edge.second()});
                full.insert({edge.first(),edge.second()});
                frontNodes.append({edge.first(),edge.second()});
            }
        }
        
        for(label iteration=1; iteration<iterations; iteration++)
        {
            DynamicList<Pair<label>> newFront;
            for(const Pair<label>& node : frontNodes)
            {
                label proc = node.first();
                label cellInd = node.second();
                
                if(proc!=-1 && cellInd!=-1)
                {
                    const List<Pair<label>>* nodeNeighbours;
                    if(proc==Pstream::myProcNo())
                    {
                        nodeNeighbours = &(localMeshGraph[cellInd]);
                    }
                    else
                    {
                        const List<List<Pair<label>>>& procHaloMeshGraph = structure.getHaloMeshGraph(proc);
                        const std::unordered_map<label,label>& procHaloCellToIndex = structure.getHaloCellToIndexMap(proc);
                        auto iter = procHaloCellToIndex.find(cellInd);
                        if(iter==procHaloCellToIndex.end())
                            FatalErrorInFunction<<"Support iteration depth mismatch!"<<exit(FatalError);
                        label procHaloIndex = iter->second;
                        const List<Pair<label>>& procHaloCellGraph = procHaloMeshGraph[procHaloIndex];
                        nodeNeighbours = &procHaloCellGraph;
                    }
                                        
                    for(const Pair<label>& node : *nodeNeighbours)
                    {
                        if(node.first()!=-1 && node.second()!=-1)
                        {
                            if(full.find(node)==full.end())
                            {
                                newFront.append(node);
                                full.insert(node);
                            }
                        }
                    }
                }
            }
            frontNodes = newFront;
            newFront.clear();
        }
    }
    
    directSupport.resize(0);
    for(auto iterCells=direct.begin(); iterCells!=direct.end(); iterCells++)
    {
        directSupport.append(*iterCells);
    }
    
    fullSupport.resize(0);
    for(auto iterCells=full.begin(); iterCells!=full.end(); iterCells++)
    {
        fullSupport.append(*iterCells);
    }

    List<List<bool>> directionsExist(3,List<bool>(2,false));
    if(markerCell!=-1)
    {
        const cell& mCell = cellList[markerCell];
        const List<Pair<label>>& mCellGraph = localMeshGraph[markerCell];
        if(mCell.size()!=mCellGraph.size())
            FatalErrorInFunction<<"Incompatible sizes"<<exit(FatalError);
        //Info<<"to_string:"<<to_string()<<Foam::endl;
        //Info<<"mCell:"<<mCell<<Foam::endl;
        for(label cellFaceInd=0; cellFaceInd<mCell.size(); cellFaceInd++)
        {
            label faceInd = mCell[cellFaceInd];
            const face& thisFace = facesList[faceInd];
            label owner = owners[faceInd];
            Pair<label> edge = mCellGraph[cellFaceInd];
            if(edge.first()!=-1 && edge.second()!=-1)
            {
                vector faceNormal = thisFace.normal(points);
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
                //Info<<faceInd<<" \t owner:"<<owner<<" nei:"<<edge<<"  normal:"<<faceNormal<<Foam::endl;
                
            }
        }
    }
    //Info<<"directionsExist:"<<directionsExist<<Foam::endl;
    for(label dim=0; dim<3; dim++)
    {
        existingDims[dim] = directionsExist[dim][0] && directionsExist[dim][1];
    }
}

Foam::Pair<Foam::vector> Foam::LagrangianMarker::minMaxNeighbourWidth
(
    const List<Pair<label>>& support
) const
{   
    vector minSpan;
    vector maxSpan;
        
    if(support.size()>1 && markerCell!=-1)
    {
        scalar min = std::numeric_limits<scalar>::min();
        scalar max = std::numeric_limits<scalar>::max();
        scalar cellSpacing = structure.initialSpacingFromMesh(mesh,markerCell);
        minSpan = vector(max,max,max);
        maxSpan = vector(min,min,min);
        Vector<bool> minSet(false,false,false);
        Vector<bool> maxSet(false,false,false);
        for(label i=0; i<support.size(); i++)
        {
            const Pair<label>& neiCellDataA = support[i];           
            if(neiCellDataA.first()!=Pstream::myProcNo() || neiCellDataA.second()!=markerCell)
            {
                vector cellCentreN;
                scalar cellVolN;
                getCellData(neiCellDataA,cellCentreN,cellVolN);
                vector centreToNeighbour = cellCentreN-markerPosition;
                        
                for(label d=0; d<3; d++)
                {
                    centreToNeighbour[d] = std::abs(centreToNeighbour[d]);
                    if(centreToNeighbour[d]<=minSpan[d])
                    {
                        minSpan[d] = centreToNeighbour[d];
                        minSet[d] = true;
                    }
                    if(centreToNeighbour[d]>=maxSpan[d])
                    {
                        maxSpan[d] = centreToNeighbour[d];
                        maxSet[d] = true;
                    }
                }
            }
        }
        for(label d=0; d<3; d++)
        {
            if(!minSet[d])
                minSpan[d] = cellSpacing;
            if(!maxSet[d])
                maxSpan[d] = cellSpacing;
        }
    }
    else
    {
        scalar min = 2*structure.initialMeshSpacing;
        scalar max = structure.initialMeshSpacing;
        minSpan = vector(max,max,max);
        maxSpan = vector(min,min,min);
    }
    return Pair<vector>(maxSpan,minSpan);
}

Foam::vector Foam::LagrangianMarker::dilationFactors
(
    Pair<vector> h
) const
{
    vector maxSpan = h.first();
    vector minSpan = h.second();
    if(markerCell!=-1)
    {
        
        scalar maxLen = std::sqrt(maxSpan&maxSpan);
        scalar minLen = std::sqrt(minSpan&minSpan);
        scalar lenFrac = minLen/maxLen;
        vector dilation;
        if(lenFrac<0.5)
            dilation = 0.8*maxSpan;
        else
            dilation = maxSpan;
        scalar eps = 0.05*minLen;
        dilation = dilation + vector(eps,eps,eps);
        if(dilation[0]==0 || dilation[1]==0 || dilation[2]==0)
            FatalErrorInFunction<<"Invalid dilation"<<exit(FatalError);
        return dilation;
    }
    else
    {
        scalar max = std::numeric_limits<scalar>::max()/1e100;
        return vector(max,max,max);
    }
}

void Foam::LagrangianMarker::checkDirectSupport() const
{
    for(const Pair<label>& cell : directSupport)
    {
        vector x;
        scalar vol;
        getCellData(cell,x,vol);
        scalar value = deltaDirac(markerPosition,x);
        if(value<=0)
            FatalErrorInFunction<<"Direct support out of function range!"<<exit(FatalError);
    }
}

void Foam::LagrangianMarker::reduceSupport()
{
    DynamicList<Pair<label>> fullSupport;
    for(const Pair<label>& cell : this->fullSupport)
    {
        vector x;
        scalar vol;
        getCellData(cell,x,vol);
        scalar value = deltaDirac(markerPosition,x);
        if(value>0)
            fullSupport.append(cell);
    }
    this->fullSupport = fullSupport;
}

void Foam::LagrangianMarker::getCellData
(
    const Pair<label>& cell,
    vector& cellCentre,
    scalar& cellVolume
) const
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    const Pair<label>& suppCellData = cell;
    if(suppCellData.first()==Pstream::myProcNo())
    {
        label cellInd = suppCellData.second();
        if(cellInd<0 || cellInd>=cells.size())
            FatalErrorInFunction<<"Out of range cell ind"<<exit(FatalError);
        cellCentre = cells[cellInd].centre(points,faces);
        cellVolume = cells[cellInd].mag(points,faces);
    }
    else
    {
        label neighProcess = suppCellData.first();
        const DynamicList<Structure::CellDescription>& neighHaloCells = structure.getHaloCellList(neighProcess);
        const std::unordered_map<label,label>& neighborHaloCellToIndexMap = structure.getHaloCellToIndexMap(neighProcess);
        auto iter = neighborHaloCellToIndexMap.find(suppCellData.second());
        if(iter==neighborHaloCellToIndexMap.end())
            FatalErrorInFunction<<"Halo cell does not exist"<<exit(FatalError);
        label index = iter->second;
        if(index<0 || index>=neighHaloCells.size())
            FatalErrorInFunction<<"Out of bounds cell number!"<<exit(FatalError);
        cellCentre = neighHaloCells[index].centre;
        cellVolume = neighHaloCells[index].volume;
    }
}

Foam::scalar Foam::LagrangianMarker::computeMoment
(
    vector indices,
    std::function<scalar(vector,vector)> deltaFunction
) const
{
    const LagrangianMarker& marker = *this;
    std::function<scalar(Pair<label>)> valueFunction = 
    [marker,indices](Pair<label> cell)
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
    
    const List<Pair<label>>* momentsCells;
    if(momentsSupportType==SupportType::Direct)
        momentsCells = &directSupport;
    else
        momentsCells = &fullSupport;
    
    return convolute<scalar>(deltaFunction,valueFunction,*momentsCells);
}

Foam::scalar Foam::LagrangianMarker::computeCorrectedMoment
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

Foam::scalar Foam::LagrangianMarker::computeMoment
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

std::unique_ptr<gismo::gsMatrix<Foam::scalar>> Foam::LagrangianMarker::computeMomentMatrix3D() const
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

std::unique_ptr<gismo::gsMatrix<Foam::scalar>> Foam::LagrangianMarker::computeMomentMatrix2D
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

std::unique_ptr<gismo::gsMatrix<Foam::scalar>> Foam::LagrangianMarker::computeMomentMatrix1D
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

std::unique_ptr<gismo::gsMatrix<Foam::scalar>> Foam::LagrangianMarker::rescalingDiagonal3D() const
{
    auto diagPtr = std::make_unique<gismo::gsMatrix<scalar>>(10,10);
    gismo::gsMatrix<scalar>& diag = *diagPtr;
    for(label row=0; row<diag.rows(); row++)
        for(label col=0; col<diag.cols(); col++)
            diag(row,col) = 0;
    
    vector dilation = getDilation();
    
    diag(0,0) = 1;
    diag(1,1) = 1/dilation[0];
    diag(2,2) = 1/dilation[1];
    diag(3,3) = 1/dilation[2];
    diag(4,4) = 1/(dilation[0]*dilation[1]);
    diag(5,5) = 1/(dilation[1]*dilation[2]);
    diag(6,6) = 1/(dilation[0]*dilation[2]);
    diag(7,7) = 1/(dilation[0]*dilation[0]);
    diag(8,8) = 1/(dilation[1]*dilation[1]);
    diag(9,9) = 1/(dilation[2]*dilation[2]);
    return diagPtr;
}

std::unique_ptr<gismo::gsMatrix<Foam::scalar>> Foam::LagrangianMarker::rescalingDiagonal2D
(
    label normalDim
) const
{
    if(normalDim!=0 && normalDim!=1 && normalDim!=2)
        FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    
    auto diagPtr = std::make_unique<gismo::gsMatrix<scalar>>(6,6);
    gismo::gsMatrix<scalar>& diag = *diagPtr;
    for(label row=0; row<diag.rows(); row++)
        for(label col=0; col<diag.cols(); col++)
            diag(row,col) = 0;
    
    vector dilation = getDilation();
    
    diag(0,0) = 1;
    
    if(normalDim==0)
    {
        diag(1,1) = 1/dilation[1];
        diag(2,2) = 1/dilation[2];
        diag(3,3) = 1/(dilation[1]*dilation[2]);
        diag(4,4) = 1/(dilation[1]*dilation[1]);
        diag(5,5) = 1/(dilation[2]*dilation[2]);
    }
    else if(normalDim==1)
    {
        diag(1,1) = 1/dilation[0];
        diag(2,2) = 1/dilation[2];
        diag(3,3) = 1/(dilation[0]*dilation[2]);
        diag(4,4) = 1/(dilation[0]*dilation[0]);
        diag(5,5) = 1/(dilation[2]*dilation[2]);
    }
    else
    {
        diag(1,1) = 1/dilation[0];
        diag(2,2) = 1/dilation[1];
        diag(3,3) = 1/(dilation[0]*dilation[1]);
        diag(4,4) = 1/(dilation[0]*dilation[0]);
        diag(5,5) = 1/(dilation[1]*dilation[1]);
    }
    return diagPtr;
}

std::unique_ptr<gismo::gsMatrix<Foam::scalar>> Foam::LagrangianMarker::rescalingDiagonal1D
(
    label dim
) const
{
    if(dim!=0 && dim!=1 && dim!=2)
        FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    
    auto diagPtr = std::make_unique<gismo::gsMatrix<scalar>>(3,3);
    gismo::gsMatrix<scalar>& diag = *diagPtr;
    for(label row=0; row<diag.rows(); row++)
        for(label col=0; col<diag.cols(); col++)
            diag(row,col) = 0;
    
    vector dilation = getDilation();
    
    diag(0,0) = 1;
    diag(1,1) = 1/dilation[dim];
    diag(2,2) = 1/(dilation[dim]*dilation[dim]);
    return diagPtr;
}

std::unique_ptr<Foam::Pair<gismo::gsMatrix<Foam::scalar>>> Foam::LagrangianMarker::computeMomentsMatrix
(
    Vector<bool> dims
) const
{
    return computeMomentsMatrix(dims,solutionStrategy);
}

std::unique_ptr<Foam::Pair<gismo::gsMatrix<Foam::scalar>>> Foam::LagrangianMarker::computeMomentsMatrix
(
    Vector<bool> dims,
    SystemSolve solutionStrategy
) const
{
    label numberDims = 0;
    DynamicList<label> listDims;
    for(label dim=0; dim<dims.size(); dim++)
    {
        if(dims[dim])
        {
            numberDims++;
            listDims.append(dim);
        }
    }

    auto system = std::make_unique<Pair<gismo::gsMatrix<scalar>>>();
    
    std::unique_ptr<gismo::gsMatrix<scalar>> matrixPtr;
    gismo::gsMatrix<scalar>& rhs = system->second();
    if(numberDims==0)
    {
        matrixPtr = std::make_unique<gismo::gsMatrix<scalar>>(1,1);
        gismo::gsMatrix<scalar>& matrix = *matrixPtr;
        matrix(0,0) = computeMoment(vector(0,0,0));
        rhs = gismo::gsMatrix<scalar>(1,1);
    }
    else if(numberDims==1)
    {
        rhs = gismo::gsMatrix<scalar>(3,1);
        if(listDims[0]==0)
        {
            matrixPtr = computeMomentMatrix1D(0);          
        }
        else if(listDims[0]==1)
        {
            matrixPtr = computeMomentMatrix1D(1);
        }
        else if(listDims[0]==2)
        {
            matrixPtr = computeMomentMatrix1D(2);
        }
        else
            FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    }
    else if(numberDims==2)
    {
        rhs = gismo::gsMatrix<scalar>(6,1);
        if(listDims[0]==0 && listDims[1]==1)
        {
            matrixPtr = computeMomentMatrix2D({0,1});
        }
        else if(listDims[0]==1 && listDims[1]==2)
        {
            matrixPtr = computeMomentMatrix2D({1,2});
        }
        else if(listDims[0]==0 && listDims[1]==2)
        {
            matrixPtr = computeMomentMatrix2D({0,2});
        }
        else
            FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    }
    else if(numberDims==3)
    {
        rhs = gismo::gsMatrix<scalar>(10,1);
        matrixPtr = computeMomentMatrix3D();
    }
    system->first() = *matrixPtr;
    for(label row=0; row<rhs.rows(); row++)
        rhs(row,0) = 0;
    rhs(0,0) = 1;
    
    gismo::gsMatrix<scalar> P(matrixPtr->rows(),matrixPtr->cols());
    for(label row=0; row<P.rows(); row++)
        for(label col=0; col<P.cols(); col++)
            P(row,col) = 0;
    
    switch (solutionStrategy)
    {
        case SystemSolve::Raw:
            for(label row=0; row<system->first().rows(); row++)
                P(row,row) = 1;
            break;
        case SystemSolve::RowAequilibration:
            for(label row=0; row<matrixPtr->rows(); row++)
            {
                scalar rowSum = 0;
                for(label col=0; col<system->first().cols(); col++)
                {
                    rowSum += system->first()(row,col);
                }
                P(row,row) = 1/rowSum;
            }
            break;
        case SystemSolve::ColAequilibration:
            for(label col=0; col<matrixPtr->cols(); col++)
            {
                scalar colSum = 0;
                for(label row=0; row<system->first().rows(); row++)
                {
                    colSum += system->first()(row,col);
                }
                P(col,col) = 1/colSum;
            }
            break;
        case SystemSolve::GeoRescaled:
            if(numberDims==0)
                P(0,0) = 1;
            else if(numberDims==1)
            {
                if(listDims[0]==0)
                    P = *rescalingDiagonal1D(0);
                else if(listDims[0]==1)
                    P = *rescalingDiagonal1D(1);
                else if(listDims[0]==2)
                    P = *rescalingDiagonal1D(2);
                else
                    FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
            }
            else if(numberDims==2)
            {
                if(listDims[0]==0 && listDims[1]==1)
                    P = *rescalingDiagonal2D(2);
                else if(listDims[0]==1 && listDims[1]==2)
                    P = *rescalingDiagonal2D(0);
                else if(listDims[0]==0 && listDims[1]==2)
                    P = *rescalingDiagonal2D(1);
                else
                    FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
            }
            else if(numberDims==3)
                P = *rescalingDiagonal3D();
            break;
        case SystemSolve::Jacobi:
            for(label col=0; col<matrixPtr->cols(); col++)
                P(col,col) = 1.0/system->first()(col,col);
            break;
        default:
            FatalErrorInFunction<<"Invalid option"<<exit(FatalError);
    }

    //std::cout<<"P:"<<std::endl<<P<<std::endl;
    
    system->first() = P*system->first();
    system->second() = P*system->second();

    return system;
}

bool Foam::LagrangianMarker::checkSolvability
(
    std::unique_ptr<Pair<gismo::gsMatrix<scalar>>>& system,
    Vector<bool> dimensions
)
{
    system = computeMomentsMatrix(dimensions);
    scalar condBaseMatrix = condition(system->first());
    if(condBaseMatrix<conditionThreshold)
    {
        //Info<<dimensions<<" :"<<condBaseMatrix<<" < "<<conditionThreshold<<Foam::endl;
        return true;
    }
    else
    {
        //Info<<dimensions<<" :"<<condBaseMatrix<<" >= "<<conditionThreshold<<Foam::endl;
        return false;
    }
}

void Foam::LagrangianMarker::searchValidConvolutionSetup
(
    std::unique_ptr<Pair<gismo::gsMatrix<scalar>>>& system
)
{   
    // Try to expand the support
    Pair<vector> newh = minMaxNeighbourWidth(fullSupport);
    h_plus = newh.first();
    h_minus = newh.second();
    dilation = dilationFactors(newh);
    computeSupport(supportWidth+2);
    if(checkSolvability(system,existingDims))
        return;
    else
        evaluateMarker();
    
    // Remove dimensions
    Vector<bool> dimensions;
    auto bitAnd = [](Vector<bool> pattern, Vector<bool> directions)
    {
        for(label dim=0; dim<directions.size(); dim++)
        {
            directions[dim] = pattern[dim] && directions[dim];
        }
        return directions;
    };
    
    // Remove a single dimension
    dimensions = bitAnd(existingDims,Vector<bool>(0,1,1));
    if(checkSolvability(system,dimensions))
    {
        existingDims = dimensions;
        return;
    }
    dimensions = bitAnd(existingDims,Vector<bool>(1,0,1));
    if(checkSolvability(system,dimensions))
    {
        existingDims = dimensions;
        return;
    }
    dimensions = bitAnd(existingDims,Vector<bool>(1,1,0));
    if(checkSolvability(system,dimensions))
    {
        existingDims = dimensions;
        return;
    }
    
    // Remove two dimensions
    dimensions = bitAnd(existingDims,Vector<bool>(1,0,0));
    if(checkSolvability(system,dimensions))
    {
        existingDims = bitAnd(existingDims,dimensions);
        return;
    }
    dimensions = bitAnd(existingDims,Vector<bool>(0,1,0));
    if(checkSolvability(system,dimensions))
    {
        existingDims = bitAnd(existingDims,dimensions);
        return;
    }
    dimensions = bitAnd(existingDims,Vector<bool>(0,0,1));
    if(checkSolvability(system,dimensions))
    {
        existingDims = bitAnd(existingDims,dimensions);
        return;
    }
    
    // Remove three dimensions
    dimensions = bitAnd(existingDims,Vector<bool>(0,0,0));
    if(checkSolvability(system,dimensions))
    {
        existingDims = bitAnd(existingDims,dimensions);
        return;
    }
    
    FatalErrorInFunction<<"Can not happen"<<exit(FatalError);
}

Foam::Vector<bool> Foam::LagrangianMarker::analyseMomentsMatrix
(
    Vector<bool> dimensions,
    const gismo::gsMatrix<scalar>& momentsMatrix
) const
{
    gismo::gsMatrix<scalar> evals,evecs;
    scalar determinant,condition;
    invertableInfo(momentsMatrix,evals,evecs,determinant,condition);
    
    std::cout<<"eval:"<<evals<<std::endl;
    std::cout<<"evecs:"<<std::endl<<evecs<<std::endl;
    std::cout<<"det:"<<determinant<<std::endl;
    std::cout<<"cond:"<<condition<<std::endl;
    
    return {false,false,false};
}

void Foam::LagrangianMarker::computeCorrectionWeights()
{
    //Pout<<Foam::endl<<"----------- computeCorrectionWeights ---------------"<<Foam::endl;
    //Info<<"existingDims:"<<existingDims<<Foam::endl;
    
    std::unique_ptr<Pair<gismo::gsMatrix<scalar>>> system;
    if(checkSolvability(system,existingDims))
    {
    }
    else
    {
        searchValidConvolutionSetup(system);
    }
    
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

    label expectedMatrixDim;
    if(numberDims==0)
        expectedMatrixDim = 1;
    else if(numberDims==1)
        expectedMatrixDim = 3;
    else if(numberDims==2)
        expectedMatrixDim = 6;
    else if(numberDims==3)
        expectedMatrixDim = 10;

    if(system->first().cols()!=expectedMatrixDim || system->first().rows()!=expectedMatrixDim)
    {
        Info<<"exp: "<<expectedMatrixDim<<Foam::endl;
        Info<<"shape: "<<system->first().cols()<<"/"<<system->first().rows()<<Foam::endl;
        FatalErrorInFunction<<"Matrix dimension wrong!"<<exit(FatalError);
    }
    if(system->second().cols()!=1 || system->second().rows()!=expectedMatrixDim)
        FatalErrorInFunction<<"RHS dimension wrong!"<<exit(FatalError);
    
    /*
    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"A:"<<std::endl<<system->first()<<std::endl;
    */
    gismo::gsMatrix<scalar> x;
    try
    {
        linearSolve(system->first(),x,system->second());
    }
    catch(...)
    {
        Info<<"condition:"<<condition(system->first())<<Foam::endl;
        Info<<"determinant:"<<determinant(system->first())<<Foam::endl;
        gismo::gsMatrix<scalar> eval,evec;
        eig(system->first(),eval,evec);
        std::cout<<eval<<std::endl;
        FatalErrorInFunction<<"Linear solver failed"<<exit(FatalError);
    }
    
    /*
    std::cout<<"x:"<<std::endl<<x<<std::endl;
    std::cout<<"b:"<<std::endl<<system->second()<<std::endl;
    std::cout<<"-------------------"<<std::endl;
    */
    
    if(x.cols()!=1 || x.rows()!=expectedMatrixDim)
        FatalErrorInFunction<<"RHS dimension wrong!"<<exit(FatalError);
    
    for(label i=0; i<10; i++)
        b[i] = 0;

    /*
    std::cout<<"-------------------"<<std::endl;
    std::unique_ptr<Pair<gismo::gsMatrix<scalar>>> checkSystem = computeMomentsMatrix(existingDims,SystemSolve::Raw);
    std::cout<<"A:"<<std::endl<<checkSystem->first()<<std::endl;
    gismo::gsMatrix<scalar> x2;
    linearSolve(checkSystem->first(),x2,checkSystem->second());
    std::cout<<"x:"<<std::endl<<x2<<std::endl;
    std::cout<<"b:"<<std::endl<<checkSystem->second()<<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;
    */
    
    /*
    std::unique_ptr<Pair<gismo::gsMatrix<scalar>>> checkSystem = computeMomentsMatrix(existingDims,WeightSystemSolve::Raw);
    gismo::gsMatrix<scalar> rhs_sol = checkSystem->first()*x;
    std::cout<<rhs_sol<<std::endl;
    */
    
    if(existingDims == Vector<bool>(0,0,0))
    /*
        b0 +
        0    + 0      + 0 +
        0    + 0      + 0 +
        0    + 0      + 0
    */
    {
        b[0] = x(0,0);
    }
    else if(existingDims == Vector<bool>(1,0,0))
    /*
        b0 +
        (x-s)b1    + 0      + 0 +
        0          + 0      + 0 +
        (x-s)²b7   + 0      + 0
    */
    {
        b[0] = x(0,0);
        b[1] = x(1,0);
        b[7] = x(2,0);            
    }
    else if(existingDims == Vector<bool>(0,1,0))
    /*
        b0 +
        0  + (y-t)b2      + 0 +
        0  + 0            + 0 +
        0  + (y-t)²b8     + 0
    */
    {
        b[0] = x(0,0);
        b[2] = x(1,0);
        b[8] = x(2,0);
    }
    else if(existingDims == Vector<bool>(0,0,1))
    /*
        b0 +
        0  + 0  + (z-v)b3 +
        0  + 0  +      0  +
        0  + 0  + (z-v)²b9
    */
    {
        b[0] = x(0,0);
        b[3] = x(1,0);
        b[9] = x(2,0);
    }
    else if(existingDims == Vector<bool>(0,1,1))
    /*
        b0 +
        0       + (y-t)b2      + (z-v)b3 +
        0       + (y-t)(z-v)b5 + 0 +
        0       + (y-t)²b8     + (z-v)²b9
    */
    {         
        b[0] = x(0,0);
        b[2] = x(1,0);
        b[3] = x(2,0);
        b[5] = x(3,0);
        b[8] = x(4,0);
        b[9] = x(5,0);
    }
    else if(existingDims == Vector<bool>(1,0,1))
    /*
        b0 +
        (x-s)b1      + 0    + (z-v)b3 +
        0            + 0    + (z-v)(x-s)b6 +
        (x-s)²b7     + 0    + (z-v)²b9
    */
    {
        b[0] = x(0,0);
        b[1] = x(1,0);
        b[3] = x(2,0);
        b[6] = x(3,0);
        b[7] = x(4,0);
        b[9] = x(5,0);
    }
    else if(existingDims == Vector<bool>(1,1,0))
    /*
        b0 +
        (x-s)b1      + (y-t)b2      + 0 +
        (x-s)(y-t)b4 + 0            + 0 +
        (x-s)²b7     + (y-t)²b8     + 0
    */
    {
        b[0] = x(0,0);
        b[1] = x(1,0);
        b[2] = x(2,0);
        b[4] = x(3,0);
        b[7] = x(4,0);
        b[8] = x(5,0);
    }
    else if(existingDims == Vector<bool>(1,1,1))
    /*
        b0 +
        (x-s)b1      + (y-t)b2      + (z-v)b3 +
        (x-s)(y-t)b4 + (y-t)(z-v)b5 + (z-v)(x-s)b6 +
        (x-s)²b7     + (y-t)²b8     + (z-v)²b9
    */
    {
        for(label i=0; i<10; i++)
            b[i] = x(i,0);
    }
}

Foam::scalar Foam::LagrangianMarker::deltaDirac
(
    vector X,
    vector x,
    scalar h
)
{
    return deltaDirac(X,x,vector(h,h,h));
}

Foam::scalar Foam::LagrangianMarker::deltaDirac
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
    return deltaDir;
}

Foam::vector Foam::LagrangianMarker::ddeltaDirac_dX
(
    vector X,
    vector x,
    vector h
)
{
    vector ddeltaDir_dX;
    for(label dim=0; dim<3; dim++)
    {
        scalar X_i_x_i = X[dim]-x[dim];
        scalar r = X_i_x_i / h[dim];
        ddeltaDir_dX[dim] = dphiFunction_dr(r);
        ddeltaDir_dX[dim] /= h[dim];
    }
    ddeltaDir_dX /= (h[0]*h[1]*h[2]);
    return ddeltaDir_dX;
}

Foam::scalar Foam::LagrangianMarker::correctedDeltaDirac
(
    vector X,
    vector x,
    scalar h,
    const FixedList<scalar,10>& b
)
{
    return correctedDeltaDirac(X,x,vector(h,h,h),b);
}

Foam::scalar Foam::LagrangianMarker::correctedDeltaDirac
(
    vector X,
    vector x,
    vector h,
    const FixedList<scalar,10>& b
)
{
    vector conn = x-X;
    scalar correctionFactor =   b[0] + 
                                conn[0]*b[1] + conn[1]*b[2] + conn[2]*b[3] +
                                conn[0]*conn[1]*b[4] + conn[1]*conn[2]*b[5] + conn[2]*conn[0]*b[6] +
                                conn[0]*conn[0]*b[7] + conn[1]*conn[1]*b[8] + conn[2]*conn[2]*b[9];
    return correctionFactor*deltaDirac(X,x,h);
}

Foam::vector Foam::LagrangianMarker::dcorrectedDeltaDirac_dX
(
    vector X,
    vector x,
    vector h,
    const FixedList<scalar,10>& b
)
{
    vector conn = x-X;
    scalar dcorrectionFactor_dX0 = -1*b[1]+ -1*conn[1]*b[4]+ conn[2]*-1*b[6]+ -2*conn[0]*b[7];
    scalar dcorrectionFactor_dX1 = -1*b[2]+ -1*conn[0]*b[4]+ -1*conn[2]*b[5]+ -2*conn[1]*b[8];
    scalar dcorrectionFactor_dX2 = -1*b[3]+ conn[1]*-1*b[5]+ -1*conn[0]*b[6]+ -2*conn[2]*b[9];
    vector ddeltaDirac_dX = LagrangianMarker::ddeltaDirac_dX(X,x,h);
    ddeltaDirac_dX[0] *= dcorrectionFactor_dX0;
    ddeltaDirac_dX[1] *= dcorrectionFactor_dX1;
    ddeltaDirac_dX[2] *= dcorrectionFactor_dX2;
    return ddeltaDirac_dX;
}

Foam::scalar Foam::LagrangianMarker::phiFunction
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

Foam::scalar Foam::LagrangianMarker::dphiFunction_dr
(
    scalar r
)
{
    scalar abs_r = std::abs(r);
    if(abs_r < 0.5)
    {
        scalar result = -1*r / std::sqrt(-3*r*r+1);
        return result;
    }
    else if(abs_r <= 1.5)
    {
        scalar result = -0.5*(r/abs_r) - 0.5*((1-abs_r)*r/abs_r) / std::sqrt(-3.0*((1-abs_r)*(1-abs_r))+1);
        return result;
    }
    else
    {
        return 0;
    }
}

Foam::scalar Foam::LagrangianMarker::correctedDeltaDirac
(
    vector X,
    vector x
) const
{
    return correctedDeltaDirac(X,x,dilation,b);
}

Foam::scalar Foam::LagrangianMarker::deltaDirac
(
    vector X,
    vector x
) const
{
    return deltaDirac(X,x,dilation);
}

Foam::vector Foam::LagrangianMarker::dcorrectedDeltaDirac_dX
(
    vector X,
    vector x
) const
{
    return dcorrectedDeltaDirac_dX(X,x,dilation,b);
}

Foam::vector Foam::LagrangianMarker::ddeltaDirac_dX
(
    vector X,
    vector x
) const
{
    return ddeltaDirac_dX(X,x,dilation);
}


std::unique_ptr<gismo::gsMatrix<Foam::scalar>> Foam::LagrangianMarker::computeCorrectedMomentMatrix() const
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
