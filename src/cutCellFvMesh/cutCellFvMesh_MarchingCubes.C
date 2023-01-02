#include "cutCellFvMesh.H"

void Foam::cutCellFvMesh::MarchingCubes::setBit(std::uint8_t& bitField, std::uint8_t bitIndex)
{
    if(bitIndex>7)
        FatalErrorInFunction<<"Bit Index must be inside [0,7]"<< exit(FatalError);

    std::uint8_t mask = 1;
    mask = mask << bitIndex;
    
    bitField = bitField | mask;
}

std::uint8_t Foam::cutCellFvMesh::MarchingCubes::getBit(std::uint8_t bitField, std::uint8_t bitIndex)
{
    if(bitIndex>7)
        FatalErrorInFunction<<"Bit Index must be inside [0,7]"<< exit(FatalError);

    std::uint8_t mask = 1;
    mask = mask << bitIndex;
    
    return bitField & mask;
}

Foam::cutCellFvMesh::MarchingCubes::MarchingCubes
(
    cutCellFvMesh& mesh
):
mesh(mesh)
{
    /*
    classifyVertex();
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const cellShapeList& meshShapes = mesh.cellShapes();
    cubeConfiguration.resize(cells.size());
    for(int i=0;i<cells.size();i++)
    {
        cell thisCell = cells[i];
        cellShape thisCellShape = meshShapes[i];
        cellModel thisCellModel = thisCellShape.model();
    }
    for(int i=0;i<cells.size();i++)
    {
        cell thisCell = cells[i];
        labelList labels = thisCell.labels(faces);
        cubeConfiguration[i] = 0;
        for(int j=0;j<labels.size();j++)
        {
            if(posVertice[i]==1)
                setBit(cubeConfiguration[i],j);
        }
    }
    */
    Info<<"Created cubeConfig"<<Foam::endl;
}

void Foam::cutCellFvMesh::MarchingCubes::classifyVertex()
{    
    scalarList& pointDist = mesh.pointDist;
    posVertice.resize(pointDist.size());
    for(int i=0;i<pointDist.size();i++)
    {
        posVertice[i] = (pointDist[i]>0)?1:0;
    }
}
