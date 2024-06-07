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
markerParameter(markerParameter),
markerPosition(LineStructure::evaluateRodPos(baseRod,markerParameter)),
markerCell(mesh.findCell(markerPosition))
{
    computeSupport();
    minMaxSupportWidth();
    dilationFactors();
}

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
        scalar maxLen = std::numeric_limits<scalar>::max();
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

std::string Foam::LagrangianMarker::to_string()
{
    return std::to_string(markerParameter);
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
    const cellList& cellList = mesh.cells();
    const faceList& facesList = mesh.faces();
    const labelList& owners = mesh.owner();
    const labelList& neighbours = mesh.neighbour();
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
    std::unordered_set<std::pair<label,label>,FirstHash> supportCells;
    //process -> {cellInd}
    std::unordered_map<label,std::unordered_set<label>> foreignHaloCells;
    this->supportCells.resize(0);
    if(markerCell!=-1)
    {
        if(markerCell<0 || markerCell>=cellList.size())
            FatalErrorInFunction<<"Invalid cell index"<< exit(FatalError);
        supportCells.insert({-1,markerCell});
        treatedCell.insert(markerCell);
        
        for(label iter=0; iter<iterations; iter++)
        {
            for(auto cellIter=supportCells.begin(); cellIter!=supportCells.end(); cellIter++)
            {
                if(cellIter->first==iter-1)
                {
                    label cellInd = cellIter->second;
                    const cell& thisCell = cellList[cellInd];
                    for(label faceInd : thisCell)
                    {
                        label owner = owners[faceInd];
                        if(treatedCell.find(owner)==treatedCell.end())
                        {
                            supportCells.insert({iter,owner});
                            treatedCell.insert(owner);
                            if(faceInd<neighbours.size())
                            {
                                auto iter = patchFaceToCell.find(faceInd);
                                if(iter!=patchFaceToCell.end())
                                {
                                    const DynamicList<std::pair<label,label>>& haloCells = iter->second;
                                    for(const std::pair<label,label> haloCell : haloCells)
                                    {
                                        label process = haloCell.first;
                                        label cellInd = haloCell.second;
                                        foreignHaloCells[process].insert(cellInd);
                                    }
                                }
                            }
                        }         
                    }
                }
            }
        }
    }
    for(auto iterCells=supportCells.begin(); iterCells!=supportCells.end(); iterCells++)
    {
        this->supportCells.append({true,Pstream::myProcNo(),iterCells->second});
    }
    for(auto iterProc=foreignHaloCells.begin(); iterProc!=foreignHaloCells.end(); iterProc++)
    {
        label processNo = iterProc->first;
        for(auto iterCell=iterProc->second.begin(); iterCell!=iterProc->second.end(); iterCell++)
        {
            label cellInd = *iterCell;
            this->supportCells.append({false,processNo,cellInd});
        }
    }
}

void Foam::LagrangianMarker::minMaxSupportWidth()
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    scalar min = std::numeric_limits<scalar>::max();
    scalar max = std::numeric_limits<scalar>::min();
    vector minSpan(max,max,max);
    vector maxSpan(min,min,min);

    /*
    for(label i=0; i<supportCells.size(); i++)
    {
        vector cellCentreA = cells[supportCells[i]].centre(points,faces);
        for(label j=0; j<supportCells.size(); j++)
        {
            if(i!=j)
            {
                vector cellCentreB = cells[supportCells[j]].centre(points,faces);
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
    */
    h_plus = maxSpan;
    h_minus = minSpan;
}

void Foam::LagrangianMarker::dilationFactors()
{
    vector minSpan = h_minus;
    vector maxSpan = h_plus;
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    scalar eps = 0.1*std::sqrt(minSpan&minSpan);
    dilation = 5.0/6.0 * maxSpan + 1.0/6.0 * minSpan + vector(eps,eps,eps);
}

scalar Foam::LagrangianMarker::computeMoment
(
    vector indices
)
{
    vector X = getMarkerPosition();
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    const DynamicList<std::tuple<bool,label,label>>& supportCells = getSupportCells();
    
    scalar moment = 0;
    for(label i=0; i<supportCells.size(); i++)
    {
        label cellInd = std::get<2>(supportCells[i]);
        vector x = cells[cellInd].centre(points,faces);
        vector conn = x-X;
        vector coeff(1,1,1);
        for(label dim=0; dim<3; dim++)
        {
            for(label index=0; index<indices[dim]; index++)
            {
                coeff[dim]*=coeff[dim];
            }
        }
        scalar volume = cells[cellInd].mag(points,faces);
        scalar dirac = deltaDirac(X,x,getDilation());
        moment += coeff[0]*coeff[1]*coeff[2]*dirac*volume;
    }
    return moment;
};

std::unique_ptr<gismo::gsMatrix<scalar>> Foam::LagrangianMarker::computeMomentMatrix()
{
    auto moments3DPtr = std::unique_ptr<gismo::gsMatrix<scalar>>(new gismo::gsMatrix<scalar>(10,10));
    gismo::gsMatrix<scalar>& moments3D = *moments3DPtr;
    
    //Diagonal
    moments3D(0,0) = computeMoment(vector(0,0,0));
    moments3D(1,1) = computeMoment(vector(2,0,0));
    moments3D(2,2) = computeMoment(vector(0,2,0));
    moments3D(3,3) = computeMoment(vector(0,0,2));
    moments3D(4,4) = computeMoment(vector(2,2,0));
    moments3D(5,5) = computeMoment(vector(0,2,2));
    moments3D(6,6) = computeMoment(vector(2,0,2));
    moments3D(7,7) = computeMoment(vector(4,0,0));
    moments3D(8,8) = computeMoment(vector(0,4,0));
    moments3D(9,9) = computeMoment(vector(0,0,4));
    
    //Row / Column 0
    moments3D(0,1) = moments3D(1,0) = computeMoment(vector(1,0,0));
    moments3D(0,2) = moments3D(2,0) = computeMoment(vector(0,1,0));
    moments3D(0,3) = moments3D(3,0) = computeMoment(vector(0,0,1));
    moments3D(0,4) = moments3D(4,0) = computeMoment(vector(1,1,0));
    moments3D(0,5) = moments3D(5,0) = computeMoment(vector(0,1,1));
    moments3D(0,6) = moments3D(6,0) = computeMoment(vector(1,0,1));
    moments3D(0,7) = moments3D(7,0) = computeMoment(vector(2,0,0));
    moments3D(0,8) = moments3D(8,0) = computeMoment(vector(0,2,0));
    moments3D(0,9) = moments3D(9,0) = computeMoment(vector(0,0,2));

    //Row / Column 1
    moments3D(1,2) = moments3D(2,1) = computeMoment(vector(1,1,0));
    moments3D(1,3) = moments3D(3,1) = computeMoment(vector(1,0,1));
    moments3D(1,4) = moments3D(4,1) = computeMoment(vector(2,1,0));
    moments3D(1,5) = moments3D(5,1) = computeMoment(vector(1,1,1));
    moments3D(1,6) = moments3D(6,1) = computeMoment(vector(2,0,1));
    moments3D(1,7) = moments3D(7,1) = computeMoment(vector(3,0,0));
    moments3D(1,8) = moments3D(8,1) = computeMoment(vector(1,2,0));
    moments3D(1,9) = moments3D(9,1) = computeMoment(vector(1,0,2));
    
    //Row / Column 2
    moments3D(2,3) = moments3D(3,2) = computeMoment(vector(0,1,1));
    moments3D(2,4) = moments3D(4,2) = computeMoment(vector(1,2,0));
    moments3D(2,5) = moments3D(5,2) = computeMoment(vector(0,2,1));
    moments3D(2,6) = moments3D(6,2) = computeMoment(vector(1,1,1));
    moments3D(2,7) = moments3D(7,2) = computeMoment(vector(2,1,0));
    moments3D(2,8) = moments3D(8,2) = computeMoment(vector(0,3,0));
    moments3D(2,9) = moments3D(9,2) = computeMoment(vector(0,1,2));
    
    //Row / Column 3
    moments3D(3,4) = moments3D(4,3) = computeMoment(vector(1,1,1));
    moments3D(3,5) = moments3D(5,3) = computeMoment(vector(0,1,2));
    moments3D(3,6) = moments3D(6,3) = computeMoment(vector(1,0,2));
    moments3D(3,7) = moments3D(7,3) = computeMoment(vector(2,0,1));
    moments3D(3,8) = moments3D(8,3) = computeMoment(vector(0,2,1));
    moments3D(3,9) = moments3D(9,3) = computeMoment(vector(0,0,3));

    //Row / Column 4
    moments3D(3,4) = moments3D(4,3) = computeMoment(vector(1,1,1));
    moments3D(3,5) = moments3D(5,3) = computeMoment(vector(0,1,2));
    moments3D(3,6) = moments3D(6,3) = computeMoment(vector(1,0,2));
    moments3D(3,7) = moments3D(7,3) = computeMoment(vector(2,0,1));
    moments3D(3,8) = moments3D(8,3) = computeMoment(vector(0,2,1));
    moments3D(3,9) = moments3D(9,3) = computeMoment(vector(0,0,3));

    //Row / Column 5
    moments3D(4,5) = moments3D(5,4) = computeMoment(vector(1,2,1));
    moments3D(4,6) = moments3D(6,4) = computeMoment(vector(2,1,1));
    moments3D(4,7) = moments3D(7,4) = computeMoment(vector(3,1,0));
    moments3D(4,8) = moments3D(8,4) = computeMoment(vector(1,3,0));
    moments3D(4,9) = moments3D(9,4) = computeMoment(vector(1,1,2));
    
    //Row / Column 6
    moments3D(5,6) = moments3D(6,5) = computeMoment(vector(1,1,2));
    moments3D(5,7) = moments3D(7,5) = computeMoment(vector(2,1,1));
    moments3D(5,8) = moments3D(8,5) = computeMoment(vector(0,3,1));
    moments3D(5,9) = moments3D(9,5) = computeMoment(vector(0,1,3));
    
    //Row / Column 7
    moments3D(6,7) = moments3D(7,6) = computeMoment(vector(3,0,1));
    moments3D(6,8) = moments3D(8,6) = computeMoment(vector(1,2,1));
    moments3D(6,9) = moments3D(9,6) = computeMoment(vector(1,0,3));
    
    //Row / Column 8
    moments3D(7,8) = moments3D(8,7) = computeMoment(vector(2,2,0));
    moments3D(7,9) = moments3D(9,7) = computeMoment(vector(2,0,2));
    
    //Row / Column 9
    moments3D(8,9) = moments3D(9,8) = computeMoment(vector(0,2,2));
    
    return moments3DPtr;
}

std::unique_ptr<std::array<scalar,10>> Foam::LagrangianMarker::rescalingDiagonal()
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

void Foam::LagrangianMarker::computeCorrectionWeights()
{
    std::unique_ptr<gismo::gsMatrix<scalar>> M = computeMomentMatrix();
    std::unique_ptr<std::array<scalar,10>> H = rescalingDiagonal();
    gismo::gsMatrix<scalar> e(10,1),c(10,1);
    for(label i=0; i<10; i++)
    {
        (*M)(i,i) = (*M)(i,i) * (*H)[i];
        e(i,0) = 0;
        c(i,0) = 0;
    }
    e(0,0) = 1;

    gismo::gsGMRes HM(*M);
    HM.solve(e,c);
    
    for(label i=0; i<10; i++)
        b[i] = c(i,0) * (*H)[i];
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

scalar Foam::LagrangianMarker::deltaDirac
(
    vector X,
    vector x
)
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
