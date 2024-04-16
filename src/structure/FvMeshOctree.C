#include "FvMeshOctree.H"

Foam::FvMeshOctree::FvMeshOctree
(
    const fvMesh& mesh
):
mesh(mesh)
{
}

void Foam::FvMeshOctree::computingBoundingBox()
{
    const pointField& meshPoints = mesh.points();
    std::pair<vector,vector> boundingBox = topBottomFromPoints(meshPoints);
    boundingBoxTop = boundingBox.first;
    boundingBoxBottom = boundingBox.second;
}

void Foam::FvMeshOctree::constructOctree
(
    const fvMesh& mesh
)
{
    root = std::unique_ptr<OctreeNode>(new OctreeNode);
    for(label cellInd=0; cellInd<mesh.cells().size(); cellInd++)
        root->cells.append(cellInd);
}

void Foam::FvMeshOctree::createNode
(
    const fvMesh& mesh,
    std::unique_ptr<OctreeNode>& node,
    uint maxNumberOfCellsInNode
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    node->seperationCoord = vector(-1,-1,-1);
    if(node->cells.size()>maxNumberOfCellsInNode)
    {
    
        vector top;
        vector bottom;
        DynamicList<vector> pointList;
        for(label cellInd : node->cells)
        {
            const cell& oneCell = cells[cellInd];
            const pointField& cellPoints = oneCell.points(faces,points);
            for(const point& onePoint : cellPoints)
            {
                pointList.append(onePoint);
            }
        }
        
        std::pair<vector,vector> boundingBox = topBottomFromPoints(pointList);
        vector boundingBoxCenter = (boundingBox.first + boundingBox.second)/2;
        node->seperationCoord = boundingBoxCenter;
        node->hasChilds = true;
        for(label childX=0; childX<2; childX++)
        {
            for(label childY=0; childY<2; childY++)
            {
                for(label childZ=0; childZ<2; childZ++)
                {
                    node->childs[childX][childY][childZ] = std::unique_ptr<OctreeNode>();
                }
            }
        }
        
        DynamicList<label> childNodeCells[2][2][2];
        for(label cellInd : node->cells)
        {
            DynamicList<std::tuple<label,label,label>> bins;
            assignCellToSide(cellInd,boundingBoxCenter,bins);
            for(std::tuple<label,label,label>& bin : bins)
            {
                label x = std::get<0>(bin);
                label y = std::get<0>(bin);
                label z = std::get<0>(bin);
                childNodeCells[x][y][z].append(cellInd);
            }
        }
        
        for(label childX=0; childX<2; childX++)
        {
            for(label childY=0; childY<2; childY++)
            {
                for(label childZ=0; childZ<2; childZ++)
                {
                    node->childs[childX][childY][childZ]->cells.append(childNodeCells[childX][childY][childZ]);
                }
            }
        }
        
        for(label childX=0; childX<2; childX++)
        {
            for(label childY=0; childY<2; childY++)
            {
                for(label childZ=0; childZ<2; childZ++)
                {
                    createNode(mesh,node->childs[childX][childY][childZ],maxNumberOfCellsInNode);
                }
            }
        }        
    }
}

std::pair<vector,vector> Foam::FvMeshOctree::topBottomFromPoints
(
    const pointField& points
) const
{
    DynamicList<vector> pointList;
    for(vector onePoint : points)
        pointList.append(onePoint);
    return topBottomFromPoints(pointList);
}

std::pair<vector,vector> Foam::FvMeshOctree::topBottomFromPoints
(
    const DynamicList<vector>& points
) const
{
    std::pair<vector,vector> boundingBox;
    vector& top = boundingBox.first;
    vector& bottom = boundingBox.second;
    for(const point& onePoint : points)
    {
        for(label dim=0; dim<3; dim++)
        {
            if(onePoint[dim]>top[dim])
                top[dim] = onePoint[dim];
            if(onePoint[dim]<bottom[dim])
                bottom[dim] = onePoint[dim];
        }
    }
}

void Foam::FvMeshOctree::assignCellToSide
(
    label cellInd,
    vector splitVec,
    DynamicList<std::tuple<label,label,label>>& bins
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    const cell& oneCell = cells[cellInd];
    const pointField& cellPoints = oneCell.points(faces,points);
    
    std::tuple<DynamicList<label>,DynamicList<label>,DynamicList<label>> sides = assignPointsToSide(cellPoints,splitVec);
    
    for(label sideX : std::get<0>(sides))
    {
        for(label sideY : std::get<1>(sides))
        {
            for(label sideZ : std::get<2>(sides))
            {
                bins.append({sideX,sideY,sideZ});
            }
        }
    }
}

std::tuple<DynamicList<label>,DynamicList<label>,DynamicList<label>> Foam::FvMeshOctree::assignPointsToSide
(
    const pointField& points,
    vector splitVec
) const
{
    std::tuple<DynamicList<label>,DynamicList<label>,DynamicList<label>> cellOnSide;
    DynamicList<label>& cellsInViaX = std::get<0>(cellOnSide);
    DynamicList<label>& cellsInViaY = std::get<1>(cellOnSide);
    DynamicList<label>& cellsInViaZ = std::get<2>(cellOnSide);
    
    std::array<DynamicList<label>*,3> cellsVia = {  &cellsInViaX,
                                                    &cellsInViaY,
                                                    &cellsInViaZ };
    for(const vector& onePoint : points)
    {
        for(label dim=0; dim<3; dim++)
        {
            if(onePoint[dim]<splitVec[dim])
            {
                cellsVia[dim]->append(0);
            }
            else
            {
                cellsVia[dim]->append(1);
            }
        }
    }
    return cellOnSide;
}

void Foam::FvMeshOctree::doOctreeInference(vector point, DynamicList<label>& cells)
{
    OctreeNode* node = root.get();
    while(node->hasChilds)
    {
        std::array<label,3> binCoord;
        for(label dim=0; dim<3; dim++)
        {
            if(point[dim]<node->seperationCoord[dim])
            {
                binCoord[dim]=0;
            }
            else
            {
                binCoord[dim]=1;
            }
        }
        node = node->childs[binCoord[0]][binCoord[1]][binCoord[2]].get();
    }
    for(label cellInd : node->cells)
    {
        cells.append(cellInd);
    }
}

