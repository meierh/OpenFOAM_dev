#include "LineStructure.H"

Foam::LagrangianMarker::LagrangianMarker
(
    const dynamicRefineFvMesh& mesh,
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* baseRod,
    const scalar markerParameter
):
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
    // std::pair<level,cell>
    
    struct FirstHash
    {
        label operator()(const std::pair<label, label> &p) const
        {
            return std::hash<label>{}(p.first);
        }
    };
    std::unordered_set<std::pair<label,label>,FirstHash> supportCells;
    std::unordered_set<std::pair<label,label>,FirstHash> supportFaces;
    this->supportCells.resize(0);
    if(markerCell!=-1)
    {
        if(markerCell<0 || markerCell>=cellList.size())
            FatalErrorInFunction<<"Invalid cell index"<< exit(FatalError);
        supportCells.insert({-1,markerCell});
        
        for(label iter=0; iter<iterations; iter++)
        {
            //DynamicList<label> newCells;
            for(auto cellIter=supportCells.begin(); cellIter!=supportCells.end(); cellIter++)
            {
                label cellInd = cellIter->first;
                const cell& thisCell = cellList[cellInd];
                for(label faceInd : thisCell)
                {
                    
                    supportCells.insert({owners[faceInd],iter});
                    if(faceInd<neighbours.size())
                        supportCells.insert({neighbours[faceInd],iter});                    
                }
            }
            //supportCells.insert(newCells.begin(),newCells.end());
        }
    }
    for(auto iterCells=supportCells.begin(); iterCells!=supportCells.end(); iterCells++)
    {
        //this->supportCells.append(*iterCells);
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
