#include "cutCellPolyMesh.H"

Foam::cutCellPolyMesh::cutCellPolyMesh
(
    const IOobject& io,
    std::function<scalar(const vector)> levelSet,
    cutStatus state
)  
:   polyMesh(io)
{
    //
    {
        const polyBoundaryMesh& boundMesh = this->boundaryMesh();
        patchStarts = labelList(boundMesh.size());
        patchSizes = labelList(boundMesh.size());
        for(int i=0;i<boundMesh.size();i++)
        {
            patchStarts[i] = boundMesh[i].start();
            patchSizes[i] = boundMesh[i].faceCentres().size();
        }
        if(patchStarts.last() != this->nFaces() || patchSizes.last() != 0)
        {
                FatalErrorInFunction
                << " The cutCellPolyMesh must be defined with an empty boundary patch at the end."
                << " Starting at: "<<this->nFaces()<<" (curr:"<<patchStarts.last()<<") and with the size: 0 ("<<patchSizes.last()<<")"
                << abort(FatalError);
        }
    }
    
    this->levelSet = levelSet;
    
    projectLevelSet();
    
    newMeshPoints();
    //printAddedPoints();
    newMeshEdges();
    //edgesToSide();
    newMeshFaces();
    cutOldFaces();
    
    faceList faces(0);
    labelList owner(0);
    labelList neighbour(0);
        
    if(state == internalCut)
    {
        createNewMeshData();
    
        faces.append(addedCutFaces);
        faces.append(splitAndUnsplitFacesInterior);
        faces.append(splitAndUnsplitFacesBoundary);
    
        owner.append(addedCutFaceOwner);
        owner.append(splitAndUnsplitFaceInteriorOwner);
        owner.append(splitAndUnsplitFaceBoundaryOwner);
    
        neighbour.append(addedCutFaceNeighbor);
        neighbour.append(splitAndUnsplitFaceInteriorNeighbor);
        neighbour.append(splitAndUnsplitFaceBoundaryNeighbor);
    }
    else if(state == delNegMesh)
    {
        createNewMeshData_cutNeg();

        faces.append(splitAndUnsplitFacesInterior);
        faces.append(splitAndUnsplitFacesBoundary);
        faces.append(addedCutFaces);
    
        owner.append(splitAndUnsplitFaceInteriorOwner);
        owner.append(splitAndUnsplitFaceBoundaryOwner);
        owner.append(addedCutFaceOwner);
        
        neighbour.append(splitAndUnsplitFaceInteriorNeighbor);
        neighbour.append(splitAndUnsplitFaceBoundaryNeighbor);
        neighbour.append(addedCutFaceNeighbor);
        
        patchStarts[patchStarts.size()-1] = (patchStarts.last()+patchSizes.last());
        patchSizes[patchSizes.size()-1] = (addedCutFaces.size());
        
        Info<<"--"<<endl;
        for(int i=0;i<patchStarts.size();i++)
        {
            Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
        }

    }
    printMesh();
    
    const pointField& oldPoints = this->points();
    const faceList& oldFaceList = this->faces();
    const cellList& oldCells = this->cells();
    
    oldCellVolume = scalarList(oldCells.size());
    for(int i=0;i<oldCells.size();i++)
    {
        oldCellVolume[i] = oldCells[i].mag(oldPoints,oldFaceList);
    }
    
    resetPrimitives(Foam::clone(newMeshPoints_),
                    Foam::clone(faces),
                    Foam::clone(owner),
                    Foam::clone(neighbour),
                    patchSizes,
                    patchStarts,
                    true);    
    
    this->write();
    printMesh();
}

void Foam::cutCellPolyMesh::pointsToSide
(
)
{
    labelList pointsToSide(newMeshPoints_.size());
    scalar lvlSet;
    for(int i=0;i<newMeshPoints_.size();i++)
    {
        lvlSet = pointDist[i];
        if(lvlSet > 0)
            pointsToSide[i] = 1;
        else if(lvlSet < 0)
            pointsToSide[i] = -1;
        else
            pointsToSide[i] = 0;
        
        //Info<<points[i]<<"\t"<<lvlSet<<"\t"<<pointsToSide[i]<<endl;
    }
    this->pointsToSide_ = pointsToSide;
}

void Foam::cutCellPolyMesh::pointsToSide
(
    const pointField& points
)
{
    labelList pointsToSide(points.size());
    scalar lvlSet;
    for(int i=0;i<points.size();i++)
    {
        lvlSet = pointDist[i];
        if(lvlSet > 0)
            pointsToSide[i] = 1;
        else if(lvlSet < 0)
            pointsToSide[i] = -1;
        else
            pointsToSide[i] = 0;
        
        //Info<<points[i]<<"\t"<<lvlSet<<"\t"<<pointsToSide[i]<<endl;
    }
    this->pointsToSide_ = pointsToSide;
}

void Foam::cutCellPolyMesh::edgesToSide
(
)
{
    const edgeList& edges = this->edges();
    labelList edgesToSide(edges.size());
    
    for(int i=0;i<edges.size();i++)
    {
        scalar lvlSetStart = pointDist[edges[i].start()];
        scalar lvlSetEnd = pointDist[edges[i].end()];

        if(lvlSetStart != 0 && lvlSetEnd != 0)
        {
            if(lvlSetStart > 0 && lvlSetEnd > 0)
                edgesToSide[i] = +1;
            else if(lvlSetStart < 0 && lvlSetEnd < 0)
                edgesToSide[i] = -1;
            else
                edgesToSide[i] = 0;
        }
        else if(lvlSetStart+lvlSetEnd == 0)
        {
            edgesToSide[i] = 0;
        }
        else
        {
            if(lvlSetStart+lvlSetEnd > 0)
                edgesToSide[i] = +1;
            else
                edgesToSide[i] = -1;
        }
    }
    this->edgesToSide_ = edgesToSide;
}

void Foam::cutCellPolyMesh::edgesToSide
(
    const edgeList& edges
)
{
    labelList edgesToSide(edges.size());
    
    for(int i=0;i<edges.size();i++)
    {
        scalar lvlSetStart = pointDist[edges[i].start()];
        scalar lvlSetEnd = pointDist[edges[i].end()];

        if(lvlSetStart != 0 && lvlSetEnd != 0)
        {
            if(lvlSetStart > 0 && lvlSetEnd > 0)
                edgesToSide[i] = +1;
            else if(lvlSetStart < 0 && lvlSetEnd < 0)
                edgesToSide[i] = -1;
            else
                edgesToSide[i] = 0;
        }
        else if(lvlSetStart+lvlSetEnd == 0)
        {
            edgesToSide[i] = 0;
        }
        else
        {
            if(lvlSetStart+lvlSetEnd > 0)
                edgesToSide[i] = +1;
            else
                edgesToSide[i] = -1;
        }
    }
    this->edgesToSide_ = edgesToSide;
}

void Foam::cutCellPolyMesh::facesToSide
(
)
{
    const faceList& faces = this->faces();
    labelList facesToSide(faces.size());
    
    for(int i=0;i<faces.size();i++)
    {
        bool nullExist = false;
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<faces[i].size();k++)
        {
            if(pointDist[faces[i][k]] == 0)
                nullExist = true;
            else if(pointDist[faces[i][k]] > 0)
                posExist = true;
            else
                negExist = true;
        }
        
        if(nullExist)
        {
            if(posExist && negExist)
                facesToSide[i] = 0;
            else if(posExist)
                facesToSide[i] = +1;
            else if(negExist)
                facesToSide[i] = -1;
            else
                facesToSide[i] = 0;
        }
        else
        {
            if(posExist && negExist)
                facesToSide[i] = 0;
            else if(posExist)
                facesToSide[i] = +1;
            else if(negExist)
                facesToSide[i] = -1;
            else
            {
                FatalErrorInFunction
                << "A face cannot have neither postive,"
                << " negative nor null points"
                << abort(FatalError);
            }
        }
    }
    this->facesToSide_ = facesToSide;
}

void Foam::cutCellPolyMesh::facesToSide
(
    const faceList& faces
)
{
    labelList facesToSide(faces.size());
    
    for(int i=0;i<faces.size();i++)
    {
        bool nullExist = false;
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<faces[i].size();k++)
        {
            if(pointDist[faces[i][k]] == 0)
                nullExist = true;
            else if(pointDist[faces[i][k]] > 0)
                posExist = true;
            else
                negExist = true;
        }
        
        if(nullExist)
        {
            if(posExist && negExist)
                facesToSide[i] = 0;
            else if(posExist)
                facesToSide[i] = +1;
            else if(negExist)
                facesToSide[i] = -1;
            else
                facesToSide[i] = 0;
        }
        else
        {
            if(posExist && negExist)
                facesToSide[i] = 0;
            else if(posExist)
                facesToSide[i] = +1;
            else if(negExist)
                facesToSide[i] = -1;
            else
            {
                FatalErrorInFunction
                << "A face cannot have neither postive,"
                << " negative nor null points"
                << abort(FatalError);
            }
        }
    }
    this->facesToSide_ = facesToSide;
}

void Foam::cutCellPolyMesh::cellsToSide
(
)
{
    const pointField& points = this->points();
    const faceList& faces = this->faces();
    const cellList& cells = this->cells();
    labelList cellsToSide(cells.size());
    
    for(int i=0;i<cells.size();i++)
    {
        pointField cellPoints = cells[i].points(faces,points);
        labelList cellLabels = cells[i].labels(faces);

        
        bool nullExist = false;
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<cellLabels.size();k++)
        {
            if(pointDist[cellLabels[k]] == 0)
                nullExist = true;
            else if(pointDist[cellLabels[k]] > 0)
                posExist = true;
            else
                negExist = true;
        }
        
        if(nullExist)
        {
            if(posExist && negExist)
                cellsToSide[i] = 0;
            else if(posExist)
                cellsToSide[i] = +1;
            else if(negExist)
                cellsToSide[i] = -1;
            else
                cellsToSide[i] = 0;
        }
        else
        {
            if(posExist && negExist)
                cellsToSide[i] = 0;
            else if(posExist)
                cellsToSide[i] = +1;
            else if(negExist)
                cellsToSide[i] = -1;
            else
            {
                FatalErrorInFunction
                << "A face cannot have neither postive,"
                << " negative nor null points"
                << abort(FatalError);
            }
        }
    }
    this->cellsToSide_ = cellsToSide;
}

void Foam::cutCellPolyMesh::cellsToSide
(
    const cellList& cells
)
{
    const pointField& points = this->points();
    const faceList& faces = this->faces();
    labelList cellsToSide(cells.size());
    
    for(int i=0;i<cells.size();i++)
    {
        pointField cellPoints = cells[i].points(faces,points);
        labelList cellLabels = cells[i].labels(faces);
        
        bool nullExist = false;
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<cellLabels.size();k++)
        {
            if(pointDist[cellLabels[k]] == 0)
                nullExist = true;
            else if(pointDist[cellLabels[k]] > 0)
                posExist = true;
            else
                negExist = true;
        }
        
        if(nullExist)
        {
            if(posExist && negExist)
                cellsToSide[i] = 0;
            else if(posExist)
                cellsToSide[i] = +1;
            else if(negExist)
                cellsToSide[i] = -1;
            else
                cellsToSide[i] = 0;
        }
        else
        {
            if(posExist && negExist)
                cellsToSide[i] = 0;
            else if(posExist)
                cellsToSide[i] = +1;
            else if(negExist)
                cellsToSide[i] = -1;
            else
            {
                FatalErrorInFunction
                << "A face cannot have neither postive,"
                << " negative nor null points"
                << abort(FatalError);
            }
        }
    }
    this->cellsToSide_ = cellsToSide;
}

void Foam::cutCellPolyMesh::projectLevelSet()
{
    const pointField& points = this->points();
    pointDist = scalarList(points.size());
    
    for(int i=0;i<points.size();i++)
    {
        pointDist[i] = levelSet(points[i]);
    }
}

void Foam::cutCellPolyMesh::newMeshPoints
(
)
{
    //Info<<"Starting adding Points"<<endl;
    const cellList& meshCells = this->cells();
    const pointField& basisPoints = this->points();
    const faceList& basisFaces = this->faces();
    const edgeList& basisEdges = this->edges();

    nbrOfPrevPoints = basisPoints.size();
    
    newMeshPoints_ = pointField(0);
    newMeshPoints_.append(basisPoints);
    
    pointsToSide(newMeshPoints_);
    
    pointToEgde_ = labelList(basisPoints.size());
    for(int i=0;i<basisPoints.size();i++)
    {
        pointToEgde_[i] = -1;
    }
    
    edgeToPoint_ = labelList(basisEdges.size());
    
    pointToFaces_ = labelListList(basisPoints.size());
    labelListList pointFaces = this->pointFaces();
    labelListList edgeFaces = this->edgeFaces();
    for(int i=0;i<basisPoints.size();i++)
    {
        if(pointsToSide_[i] == 0)
        {
            pointToFaces_[i] = pointFaces[i];
        }
    }
        
    faceToPoints_ = labelListList(basisFaces.size());
    for(int i=0;i<basisPoints.size();i++)
    {
        if(pointToFaces_[i].size() != 0)
        {
            for(int k=0;k<pointToFaces_[i].size();k++)
            {
                if(faceToPoints_[pointToFaces_[i][k]].size() == 0)
                {
                    faceToPoints_[pointToFaces_[i][k]] = labelList(0);
                    faceToPoints_[pointToFaces_[i][k]].append(i);
                }
                else
                {
                    faceToPoints_[pointToFaces_[i][k]].append(i);
                }
            }
        }
    }
    
    pointToCells_ = labelListList(basisPoints.size());
    labelListList pointCells = this->pointCells();
    labelListList edgeCells = this->edgeCells();
    for(int i=0;i<basisPoints.size();i++)
    {
        if(pointsToSide_[i] == 0)
        {
            pointToCells_[i] = pointCells[i];
        }
    }
    
    cellToPoints_ = labelListList(meshCells.size());
    for(int i=0;i<basisPoints.size();i++)
    {
        if(pointToCells_[i].size() != 0)
        {
            for(int k=0;k<pointToCells_[i].size();k++)
            {
                if(cellToPoints_[pointToCells_[i][k]].size() == 0)
                {
                    cellToPoints_[pointToCells_[i][k]] = labelList(0);
                    cellToPoints_[pointToCells_[i][k]].append(i);
                }
                else
                {
                    cellToPoints_[pointToCells_[i][k]].append(i);
                }
            }
        }
    }
    
    
    
    label pos,neg;
    
    for(int i=0;i<basisEdges.size();i++)
    {
        label startLabel = basisEdges[i].start();
        label endLabel = basisEdges[i].end();        
        pos = 0;
        neg = 0;
        scalar phiStart = pointDist[startLabel];
        scalar phiEnd = pointDist[endLabel];
        
        //Info<<phiStart<<" -> "<<phiEnd<<endl;
        
        if(phiStart>0 || phiEnd>0)
            pos = +1;
        if(phiStart<0 || phiEnd<0)
            neg = -1;
        
        if(pos == +1 && neg == -1)
        {
            Info<<basisPoints[startLabel]<<" -> "<<basisPoints[endLabel]<<endl;
            Info<<phiStart<<" -> "<<phiEnd<<endl;
            vector startToEnd = basisEdges[i].vec(basisPoints);
            Info<<"startToEnd: "<<basisEdges[i].vec(basisPoints)<<endl;;
            scalar norm_startToEnd = basisEdges[i].mag(basisPoints);
            Info<<"norm_startToEnd: "<<norm_startToEnd<<endl;;
            scalar distPhi = std::abs(phiEnd-phiStart);
            Info<<"distPhi: "<<distPhi<<endl;;
            scalar norm_phiStart = std::abs(phiStart);
            Info<<"norm_phiStart: "<<norm_phiStart<<endl;;
            scalar scalePoint = phiStart / (phiStart - phiEnd);
            Info<<"scalePoint: "<<scalePoint<<endl;;
            vector newPoint = basisPoints[startLabel] + scalePoint * startToEnd;
            
            Info<<"Added: "<<newPoint<<endl<<endl;

            newMeshPoints_.append(newPoint);
            
            pointsToSide_.append(0);
            
            pointToEgde_.append(i);
            
            edgeToPoint_[i] = newMeshPoints_.size()-1;
            
            pointToFaces_.append(edgeFaces[i]);
            
            for(int k=0;k<edgeFaces[i].size();k++)
            {
                label faceLabel = edgeFaces[i][k];
                if(faceToPoints_[faceLabel].size() == 0)
                {
                    faceToPoints_[faceLabel] = labelList(0);
                    faceToPoints_[faceLabel].append(newMeshPoints_.size()-1);
                }
                else
                {
                    faceToPoints_[faceLabel].append(newMeshPoints_.size()-1);
                }
            }
            
            pointToCells_.append(edgeCells[i]);
            
            for(int k=0;k<edgeCells[i].size();k++)
            {
                label cellLabel = edgeCells[i][k];
                if(cellToPoints_[cellLabel].size() == 0)
                {
                    cellToPoints_[cellLabel] = labelList(0);
                    cellToPoints_[cellLabel].append(newMeshPoints_.size()-1);
                }
                else
                {
                    cellToPoints_[cellLabel].append(newMeshPoints_.size()-1);
                }
            }
        }
        else
        {
            edgeToPoint_[i] = -1;
        }
    }
}

void Foam::cutCellPolyMesh::printAddedPoints
(
)
{
    Info<<"---------------------------------------AddedPoints-----------------------------------"<<endl;
    for(int k=nbrOfPrevPoints;k<newMeshPoints_.size();k++)
    {
        Info<<"Point added: "<<k<<"-"<<newMeshPoints_[k]<<"\t"<<"at cells: ";
        for(int l=0;l<pointToCells_[k].size();l++)
            Info<<pointToCells_[k][l]<<" ";
        Info<<"and at faces: ";
        for(int l=0;l<pointToFaces_[k].size();l++)
            Info<<pointToFaces_[k][l]<<" ";
        Info<<"and at edge: ";
        Info<<pointToEgde_[k]<<endl;
    }
    Info<<"--------------------------------------CelltoPoints-----------------------------------"<<endl;
    for(int k=0;k<cellToPoints_.size();k++)
    {
        Info<<"Cell "<<k<<" has added Points: ";
        for(int j=0;j<cellToPoints_[k].size();j++)
        {
            int index = cellToPoints_[k][j];
            Info<<index<<"-"<<newMeshPoints_[index]<<" ";
        }
        Info<<endl;
    }
    Info<<"--------------------------------------FacetoPoints-----------------------------------"<<endl;
    for(int k=0;k<faceToPoints_.size();k++)
    {
        if(faceToPoints_[k].size() == 0)
            continue;
        Info<<"Face "<<k<<" has added Points: ";
        for(int j=0;j<faceToPoints_[k].size();j++)
        {
            int index = faceToPoints_[k][j];
            Info<<index<<"-"<<newMeshPoints_[index]<<" ";
        }
        Info<<endl;
    }
    Info<<"--------------------------------------EdgetoPoints-----------------------------------"<<endl;
    for(int k=0;k<edgeToPoint_.size();k++)
    {
        int index;
        if((index = edgeToPoint_[k])==-1)
            continue;
        Info<<"Edge "<<k<<" has added Point: ";
        Info<<index<<"-"<<newMeshPoints_[index]<<" ";
        Info<<endl;
    }
    
    Info<<"--------------------------------------PointsToSide-----------------------------------"<<endl;
    for(int k=0;k<newMeshPoints_.size();k++)
    {
        Info<<"Point "<<k<<" "<<newMeshPoints_[k]<<" Side: "<<pointsToSide_[k];
        Info<<endl;
    }
    
    /*
    Info<<"---------------------------------------AddedEdges------------------------------------"<<endl;
    for(int k=0;k<addedMeshItems.addedEdges.size();k++)
    {
        Info<<"Edge added: -"<<k<<"-"<< newMeshPoints_[addedMeshItems.addedEdges[k].start()]
        <<"->"<<
        newMeshPoints_[addedMeshItems.addedEdges[k].end()]<<"\t"<<"at cells: ";
        for(int l=0;l<addedMeshItems.addedEdgeToOldCells[k].size();l++)
            Info<<addedMeshItems.addedEdgeToOldCells[k][l]<<" ";
        Info<<"and at faces: ";
        Info<<addedMeshItems.addedEdgeToOldFace[k]<<endl;
    }        
    Info<<"---------------------------------------CelltoEdge------------------------------------"<<endl;
    for(int k=0;k<addedMeshItems.oldCellsToAddedEdges.size();k++)
    {
        Info<<"Cell "<<k<<" has added Edges: ";
        for(int j=0;j<addedMeshItems.oldCellsToAddedEdges[k].size();j++)
        {
            int index = addedMeshItems.oldCellsToAddedEdges[k][j];
            Info<<" -"<<index<<"-" <<newMeshPoints_[addedMeshItems.addedEdges[index].start()]<<"->"<<
            newMeshPoints_[addedMeshItems.addedEdges[index].end()];
        }
        Info<<endl;
    }
    Info<<"---------------------------------------FacetoEdge------------------------------------"<<endl;
    for(int k=0;k<addedMeshItems.oldFacesToAddedEdge.size();k++)
    {
        if(addedMeshItems.oldFacesToAddedEdge[k] == -1)
            continue;
        Info<<"Face "<<k<<" has added Edge: ";
        int index = addedMeshItems.oldFacesToAddedEdge[k];
        Info<<"-"<<index<<"-"<<
        newMeshPoints_[addedMeshItems.addedEdges[index].start()]<<"->"<<
        newMeshPoints_[addedMeshItems.addedEdges[index].end()];
        Info<<endl;
    }
    Info<<"---------------------------------------AddedFaces------------------------------------"<<endl;
    for(int k=0;k<addedMeshItems.addedFaces.size();k++)
    {
        Info<<"Face added: -"<<k<<"-";
        pointField facePoints = addedMeshItems.addedFaces[k].points(newMeshPoints_);
        for(int j=0;j<facePoints.size();j++)
            Info<<facePoints[j]<<"->";
        Info<<"\t"<<"with face centre at: "<<addedMeshItems.addedFaces[k].centre(newMeshPoints_);
        Info<<"at cell: "<<addedMeshItems.addedFaceToOldCells[k]<<endl;
    }   
    Info<<"---------------------------------------CelltoFace------------------------------------"<<endl;
    for(int k=0;k<addedMeshItems.oldCellsToAddedFace.size();k++)
    {
        Info<<"Cell "<<k<<" has added Face: -"<<addedMeshItems.oldCellsToAddedFace[k]<<"-";
        Info<<endl;
    }
    Info<<"---------------------------------------CelltoSide------------------------------------"<<endl;
    for(int i=0;i<addedMeshItems.oldPointsAtCutCellsToSide.size();i++)
    {
        Info<<"Point "<<i<<" is on side: "<<addedMeshItems.oldPointsAtCutCellsToSide[i]<<endl;
    }
    */    
}

void Foam::cutCellPolyMesh::newMeshEdges
(
)
{
    Info<<"Starting adding Edges"<<endl;
    const cellList& meshCells = this->cells();
    //const pointField& basisPoints = this->points();
    const faceList& basisFaces = this->faces();
    const edgeList& basisEdges = this->edges();
    
    nbrOfPrevEdges = basisEdges.size();
    
    newMeshEdges_ = edgeList(0);
    newMeshEdges_.append(basisEdges);

    edgesToSide(newMeshEdges_);
    
    /*
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        Info<<"Edge:"<<i<<" Side:"<<edgesToSide_[i]<<endl;
    }
    */
    
    Info<<"Put edges to side"<<endl;
    
    edgeToFaces_ = labelListList(basisEdges.size());
    labelListList edgeFaces = this->edgeFaces();
    for(int i=0;i<basisEdges.size();i++)
    {
        label startPoint = basisEdges[i].start();
        label endPoint = basisEdges[i].end();
        
        if(pointsToSide_[startPoint] == 0 && pointsToSide_[endPoint] == 0)
        {
            edgeToFaces_[i] = edgeFaces[i];
        }
    }
    
    /*
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        Info<<"Edge:"<<i;
        Info<<" from "<<newMeshEdges_[i].start()<<newMeshPoints_[newMeshEdges_[i].start()]<<"->";
        Info<<newMeshEdges_[i].end()<<newMeshPoints_[newMeshEdges_[i].end()];
        for(int k=0;k<edgeToFaces_[i].size();k++)
        {
            Info<<" face:"<<edgeToFaces_[i][k];
        }
        Info<<" Side:"<<edgesToSide_[i]<<endl;
    }
    */
    
    labelList thisFacePoints;
    for(int i=0;i<basisFaces.size();i++)
    {
        thisFacePoints = faceToPoints_[i];
        if(thisFacePoints.size() > 2)
        {
            bool allPointsOld = true;
            for(int k=0;k<thisFacePoints.size();k++)
            {
                if(thisFacePoints[k] >= nbrOfPrevPoints)
                    allPointsOld = false;
            }
            if(!allPointsOld)
            {
                FatalErrorInFunction
                << "A face cannot have "<< thisFacePoints.size()
                << " cut points while one or more "
                << "cut points are not old points! "
                << abort(FatalError);
            }
            bool allPointsBelongToOldEdges = true;
            for(int k=0;k<thisFacePoints.size();k++)
            {
                label localIndx = basisFaces[i].which(thisFacePoints[k]);
                label nextGlobalIndx = basisFaces[i].nextLabel(localIndx);
                label prevGlobalIndx = basisFaces[i].prevLabel(localIndx);
                if(pointsToSide_[nextGlobalIndx] != 0 && pointsToSide_[prevGlobalIndx] != 0)
                {
                    allPointsBelongToOldEdges = false;
                }
            }
            if(!allPointsBelongToOldEdges)
            {
                FatalErrorInFunction
                << "Face has "<< thisFacePoints.size()
                << " cut points while one or more "
                << "cut points are connected to the other ones! "
                << abort(FatalError);
            }            
        }
        else if(thisFacePoints.size() == 2)
        {
            label pt0 = thisFacePoints[0];
            label pt1 = thisFacePoints[1];
            // two old points
            if(pt0 < nbrOfPrevPoints && pt1 < nbrOfPrevPoints)
            {
                //Two not connected points
                label localIndx = basisFaces[i].which(pt0);
                label nextGlobalIndx = basisFaces[i].nextLabel(localIndx);
                label prevGlobalIndx = basisFaces[i].prevLabel(localIndx);
                //that dont lie on one edge
                if(nextGlobalIndx != pt1 && prevGlobalIndx != pt1)
                {
                    newMeshEdges_.append(edge(pt0,pt1));
                    edgesToSide_.append(0);
                    edgeToFaces_.append(labelList(0));
                    edgeToFaces_[edgeToFaces_.size()-1].append(i);

                }
            }
            //one or more new point
            else
            {
                newMeshEdges_.append(edge(pt0,pt1));
                edgesToSide_.append(0);
                edgeToFaces_.append(labelList(0));
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
            }
        }
    }
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        Info<<"Edge:"<<i;
        Info<<" from "<<newMeshEdges_[i].start()<<newMeshPoints_[newMeshEdges_[i].start()]<<"->";
        Info<<newMeshEdges_[i].end()<<newMeshPoints_[newMeshEdges_[i].end()];
        for(int k=0;k<edgeToFaces_[i].size();k++)
        {
            Info<<" face:"<<edgeToFaces_[i][k];
        }
        Info<<" Side:"<<edgesToSide_[i]<<endl;
    }
    Info<<"edge to face done"<<endl;
    
    faceToEdges_ = labelListList(basisFaces.size());
    Info<<faceToEdges_.size()<<endl;
    Info<<newMeshEdges_.size()<<endl;
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        //Info<<"Edge "<<i<<endl;
        if(edgeToFaces_[i].size() != 0)
        {
            //Info<<"\t"<<"edgeToFaces.size():"<<edgeToFaces_[i].size()<<endl;
            for(int k=0;k<edgeToFaces_[i].size();k++)
            {
                /*
                if(faceToEdges_[edgeToFaces_[i][k]].size() == 0)
                {
                    faceToEdges_[edgeToFaces_[i][k]] = labelList(0);
                    faceToEdges_[edgeToFaces_[i][k]].append(i);
                }
                else
                {
                    faceToEdges_[edgeToFaces_[i][k]].append(i);
                }
                */
                //Info<<"\t\t"<<"add edge "<<i<<" to face: "<<edgeToFaces_[i][k]<<endl;
                faceToEdges_[edgeToFaces_[i][k]].append(i);
            }
        }
    }
    Info<<"face to edge done"<<endl;
    
    Info<<faceToEdges_.size()<<endl;
    
    /*
    for(int i=0;i<faceToEdges_.size();i++)
    {
        Info<<"Face:"<<i<<" Edges:";
        for(int k=0;k<faceToEdges_[i].size();k++)
        {
            Info<<faceToEdges_[i][k]<<"-";   
        }
        Info<<endl;
    }
    */

    
    edgeToCells_ = labelListList(newMeshEdges_.size());
    labelList owner = this->faceOwner();
    labelList neighbour = this->faceNeighbour();
    labelListList edgeCells = this->edgeCells();
    
    Info<<"newMeshEdges_: "<<newMeshEdges_.size()<<endl;
    Info<<"edgeToCells_: "<<edgeToCells_.size()<<endl;
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        Info<<"Edge: "<<i<<endl;
        if(edgeToFaces_[i].size() != 0)
        {
            Info<<"edgeToFaces_["<<i<<"]: "<<edgeToFaces_[i].size()<<endl;
            if(i<nbrOfPrevEdges)
            {
                edgeToCells_[i] = Foam::clone(edgeCells[i]);
            }
            else
            {
                if(edgeToFaces_[i].size() != 1)
                {
                    FatalErrorInFunction
                    << "Added Edge has  "<< edgeToFaces_[i].size()
                    << " neighboring faces instead of one neighboring face! "
                    << abort(FatalError);
                }
                label thisFace = edgeToFaces_[i][0];
                Info<<"thisFace: "<<thisFace<<endl;
                Info<<"edgeToCells_[i].size(): "<<edgeToCells_[i].size()<<endl;
                Info<<"Owner of thisFace: "<<owner[thisFace]<<endl;
                edgeToCells_[i].append(owner[thisFace]);
                Info<<"Peter Pan"<<endl;
                if(thisFace < neighbour.size())
                    edgeToCells_[i].append(neighbour[thisFace]);
            }
        }
    }
    
    /*
    for(int i=0;i<edgeToCells_.size();i++)
    {
        Info<<"Edge "<<i; 
        for(int k=0;k<edgeToCells_[i].size();k++)
        {
            Info<<" cell:"<<edgeToCells_[i][k];
        }
        Info<<endl;
    }
    */
    
    Info<<"edge to cell done"<<endl;

    
    cellToEdges_ = labelListList(meshCells.size());
    Info<<"nCells: "<<meshCells.size()<<endl;
    Info<<"newMeshEdges_: "<<newMeshEdges_.size()<<endl;
    Info<<"edgeToCells_: "<<edgeToCells_.size()<<endl;
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        if(edgeToCells_[i].size() != 0)
        {
            Info<<"edgeToCells_["<<i<<"]: "<<edgeToCells_[i].size()<<endl;
            for(int k=0;k<edgeToCells_[i].size();k++)
            {
                Info<<"edgeToCells_["<<i<<"]["<<k<<"]:"<<edgeToCells_[i][k]<<endl;
                /*
                if(cellToEdges_[edgeToCells_[i][k]].size() == 0)
                {
                    cellToEdges_[edgeToCells_[i][k]] = labelList(0);
                    cellToEdges_[edgeToCells_[i][k]].append(i);
                }
                else
                {
                    cellToEdges_[edgeToCells_[i][k]].append(i);
                }
                */
                cellToEdges_[edgeToCells_[i][k]].append(i);
            }
        }
    }
    Info<<"cell to edge done"<<endl;
}

void Foam::cutCellPolyMesh::printAddedEdges
(
)
{
    Info<<"---------------------------------------AddedEdges------------------------------------"<<endl;
    for(int k=nbrOfPrevEdges;k<newMeshEdges_.size();k++)
    {
        Info<<"Edge added: -"<<k<<"-"<< newMeshPoints_[newMeshEdges_[k].start()]
        <<"->"<<
        newMeshPoints_[newMeshEdges_[k].end()]<<"\t"<<"at cells: ";
        for(int l=0;l<edgeToCells_[k].size();l++)
            Info<<edgeToCells_[k][l]<<" ";
        Info<<"and at faces: ";
        Info<<edgeToFaces_[k]<<endl;
    }        
    Info<<"---------------------------------------CelltoEdge------------------------------------"<<endl;
    for(int k=0;k<cellToEdges_.size();k++)
    {
        Info<<"Cell "<<k<<" has added Edges: ";
        for(int j=0;j<cellToEdges_[k].size();j++)
        {
            int index = cellToEdges_[k][j];
            Info<<" -"<<index<<"-" <<newMeshPoints_[newMeshEdges_[index].start()]<<"->"<<
            newMeshPoints_[newMeshEdges_[index].end()];
        }
        Info<<endl;
    }
    Info<<"---------------------------------------FacetoEdge------------------------------------"<<endl;
    for(int k=0;k<faceToEdges_.size();k++)
    {
        if(faceToEdges_[k].size() == 0)
            continue;
        Info<<"Face "<<k;
        labelList indexes = faceToEdges_[k];
        if(indexes.size() == 1)
        {
            Info<<" has added cut edge:" <<
            newMeshPoints_[newMeshEdges_[indexes[0]].start()]<<"->"<<
            newMeshPoints_[newMeshEdges_[indexes[0]].end()];
        }
        else
        {
            
            Info<<" has old cut edges:";
            for(int j=0;j<indexes.size();j++)
            {
                Info<<newMeshPoints_[newMeshEdges_[indexes[j]].start()]
                <<"->"<<newMeshPoints_[newMeshEdges_[indexes[j]].end()]<<endl;
            }
        }
        Info<<endl;
    }
    Info<<"--------------------------------------EdgesToSide-----------------------------------"<<endl;
    for(int k=0;k<newMeshEdges_.size();k++)
    {
        Info<<"Edge "<<k<<" "<<newMeshPoints_[newMeshEdges_[k].start()]<<"->"<<
        newMeshPoints_[newMeshEdges_[k].end()]<<" Side: "<<edgesToSide_[k];
        Info<<endl;
    }    
}

void Foam::cutCellPolyMesh::newMeshFaces
(
)
{
    Info<<"Starting adding Faces"<<endl;
    const cellList& meshCells = this->cells();
    //const pointField& basisPoints = this->points();
    const faceList& basisFaces = this->faces();
    //const edgeList& basisEdges = this->edges();
    
    nbrOfPrevFaces = basisFaces.size();
    
    newMeshFaces_ = faceList(0);
    newMeshFaces_.append(basisFaces);

    facesToSide(newMeshFaces_);
    
    /*
     * A existing face is a cut face if all its points are cut points.
     * Only if thats true the respective face is connected to its connected
     * cells via the faceToCells_ list. A non empty faceToCells_ for a face 
     * is the sign that the face is a cut face
     */
    faceToCells_ = labelListList(basisFaces.size());
    labelList owner = this->faceOwner();
    labelList neighbour = this->faceNeighbour();
    for(int i=0;i<basisFaces.size();i++)
    {
        labelList facePoints = basisFaces[i];
        
        Info<<"  Face:"<<i<<" Owner:"<<owner[i]<<" ";
        if(i < neighbour.size())
            Info<<" Neighbor:"<<neighbour[i]<<" ";
        for(int k=0;k<basisFaces[i].size();k++)
        {
            Info<<newMeshPoints_[basisFaces[i][k]]<<"->";
        }
        Info<<" with centre:"<<basisFaces[i].centre(newMeshPoints_);
        Info<<" and normal vector:"<<basisFaces[i].normal(newMeshPoints_);
        Info<<" and area:"<<basisFaces[i].mag(newMeshPoints_);

        /*
        if(facePoints.size() != 4)
        {
            FatalErrorInFunction
            << "A face cannot have "<< facePoints.size()
            << " cut points while one or more "
            << "cut points are not old points! "
            << abort(FatalError);
        }
        */
        bool isCutFace = true;
        for(int k=0;k<facePoints.size();k++)
        {
            if(pointsToSide_[facePoints[k]] != 0)
                isCutFace = false;
        }
        Info<<"\t"<<"is CutFace:"<<isCutFace<<endl;
        if(isCutFace)
        {
            faceToCells_[i].append(owner[i]);
            //Info<<"I want a neighbour"<<endl;
            if(i < neighbour.size())
                faceToCells_[i].append(neighbour[i]);
                //Info<<"I have a neighbour"<<endl;
        }        
    }

    /*
     * The following creates the added faces as well as the cellToFaces_ list 
     * and finishes the faceToCells_ list of the previous part with the new faces.
     * 
     */
    cellToFaces_ = labelListList(meshCells.size());
    for(int i=0;i<meshCells.size();i++)
    {    
        /* 1)
         * If a cell has no cut edges it is not cut by a face
         */
        if(cellToEdges_[i].size() == 0)
        {
            continue;
        }
        
        labelList facePoints;
        labelList cellCutEdgeList = cellToEdges_[i];
        
        /* 2)
         * If a cell has one cut edge it is not cut. It should be impossible that a
         * cell has two cut edges because the minimum cut face thinkable is a three edge
         * face. Because of that a failure abort is called if these states appear 
         */
        if(cellCutEdgeList.size() <= 2)
        {
            FatalErrorInFunction
            << "A cell cannot be cut by "<< cellCutEdgeList.size()
            << " edges! "
            << abort(FatalError);
        }
        
        labelListList cellCutEdgeFacesList = labelListList(cellCutEdgeList.size());
        for(int k=0;k<cellCutEdgeList.size();k++)
        {
            cellCutEdgeFacesList[k] = edgeToFaces_[cellCutEdgeList[k]];
        }       
        
        bool allCutEdgesOld = true;
        //label oldCutEdgeFace = ;
        for(int k=0;k<cellCutEdgeList.size();k++)
        {
            if(cellCutEdgeList[k] >= nbrOfPrevEdges)
            {
                allCutEdgesOld = false;
            }
        }
        
        if(allCutEdgesOld)
        {
            labelList equalFace(0);
            for(int a=0;a<cellCutEdgeFacesList[0].size();a++)
            {
                for(int b=0;b<cellCutEdgeFacesList[1].size();b++)
                {
                    if(cellCutEdgeFacesList[0][a] == cellCutEdgeFacesList[1][b])
                        equalFace.append(cellCutEdgeFacesList[0][a]);
                }
            }
            if(equalFace.size() >= 2)
            {
                FatalErrorInFunction
                << "Two cut edges can not have "<< equalFace.size()
                << " equal faces! "
                << abort(FatalError);
            }
            if(equalFace.size() == 0)
            {
                FatalErrorInFunction
                << "Cell is cut by two old edges that do not"
                << " share a face! This must not happen! "
                << abort(FatalError);
            }
            bool allCutEdgesSameFace = true;
            for(int k=2; k<cellCutEdgeFacesList.size();k++)
            {
                bool sameCutedgeFace = false;
                for(int a=0;a<cellCutEdgeFacesList[k].size();a++)
                {
                    if(cellCutEdgeFacesList[k][a] == equalFace[0])
                        sameCutedgeFace = true;
                }
                if(!sameCutedgeFace)
                    allCutEdgesSameFace = false;
            }
            if(!allCutEdgesSameFace)
            {
                FatalErrorInFunction
                << "Cell is cut by "<<cellCutEdgeFacesList.size()
                << " old edges that do not"
                << " share a face! This must not happen! "
                << abort(FatalError);
            }
            
            cellToFaces_[i] = equalFace;
        }
        
        /*
        for(int k=0;k<cellCutEdgeList.size();k++)
        {
            Info<<"  Edge:"<<cellCutEdgeList[k]<<endl;
            Info<<"    from Point "<<"-"<<addedEdges[cellCutEdgeList[k]].start()<<"-";
            Info<<addedPoints[addedEdges[cellCutEdgeList[k]].start()]<<"->";
            Info<<"-"<<addedEdges[cellCutEdgeList[k]].end()<<"-";
            Info<<addedPoints[addedEdges[cellCutEdgeList[k]].end()]<<endl;
        }
        
        Info<<addedEdges[cellCutEdgeList[0]].start()<<endl;
        Info<<addedEdges[cellCutEdgeList[0]].end()<<endl;
        */
//////////////////////////////////////////////////////////////////////////////////////////////
        facePoints.append(newMeshEdges_[cellCutEdgeList[0]].start());
        facePoints.append(newMeshEdges_[cellCutEdgeList[0]].end());
        label frontPoint = facePoints[facePoints.size()-1];
        label endPoint = facePoints[0];
        //Info<<i<<endl;        
        
        label currentEdge = cellCutEdgeList[0];
        bool closedFace = false;
        //Info<<"endPoint:"<<endPoint<<endl;
        //Info<<"Face: "<<facePoints[0]<<" "<<facePoints[1];
        //Info<<endl;
        int count = 0;
        while(!closedFace)
        {
            count++;
            for(int k=0;k<cellCutEdgeList.size();k++)
            {
                if(cellCutEdgeList[k] == currentEdge)
                    continue;
                /*
                Info<<"\tEdge:"<<cellCutEdgeList[k];
                Info<<" from Point "<<"-"<<addedEdges[cellCutEdgeList[k]].start()<<"-";
                Info<<addedPoints[addedEdges[cellCutEdgeList[k]].start()]<<"->";
                Info<<"-"<<addedEdges[cellCutEdgeList[k]].end()<<"-";
                Info<<addedPoints[addedEdges[cellCutEdgeList[k]].end()]<<endl;
                Info<<"frontPoint:"<<frontPoint<<endl;
                Info<<"currentEdge:"<<currentEdge<<endl;
                */
                
                if(newMeshEdges_[cellCutEdgeList[k]].start()==frontPoint)
                {
                    //Info<<"\tAdd edge:"<<cellCutEdgeList[k];
                    frontPoint = newMeshEdges_[cellCutEdgeList[k]].end();
                    if(frontPoint == endPoint)
                    {
                        closedFace = true;
                        /*
                        Info<<endl<<"----------------------Closed----------------------------------"<<endl;
                        for(int t=0;t<facePoints.size();t++)
                            Info<<facePoints[t]<<" ";
                        Info<<endl;
                        */
                        break;
                    }
                    facePoints.append(frontPoint);
                    currentEdge = cellCutEdgeList[k];
                    /*
                    Info<<" and Point:"<<facePoints[facePoints.size()-1]<<endl;
                    Info<<"--------------------------------------------------------------"<<endl;
                    */
                    continue;
                }
                if(newMeshEdges_[cellCutEdgeList[k]].end()  ==frontPoint)  
                {
                    //Info<<"\tAdd edge:"<<cellCutEdgeList[k];
                    frontPoint = newMeshEdges_[cellCutEdgeList[k]].start();
                    if(frontPoint == endPoint)
                    {
                        closedFace = true;
                        /*
                        Info<<endl<<"----------------------Closed----------------------------------"<<endl;
                        for(int t=0;t<facePoints.size();t++)
                            Info<<facePoints[t]<<" ";
                        Info<<endl;
                        */
                        break;
                    }
                    facePoints.append(frontPoint);
                    currentEdge = cellCutEdgeList[k];
                    /*
                    Info<<" and Point:"<<facePoints[facePoints.size()-1]<<endl;
                    Info<<"--------------------------------------------------------------"<<endl;
                    */
                    continue;
                }
            }
            if(count>10)
            {
                Info<<"---------------------------Unrealistic face with more than 10 edges";
                Info<<"---------------------------"<<endl;
                break;
            }                
        }
        //add to addedFaces
        face newFace(facePoints);
        //newFace = newFace.reverseFace();
        
        newMeshFaces_.append(newFace);
        facesToSide_.append(0);
        
        labelList newFaceCell(0);
        newFaceCell.append(i);
        faceToCells_.append(newFaceCell);
        
        //input to oldCellsToAddedFace
        cellToFaces_[i].append(newMeshFaces_.size()-1);
    }
    Info<<"End adding faces"<<endl;
}

void Foam::cutCellPolyMesh::printAddedFaces
(
)
{
    Info<<"---------------------------------------AddedFaces------------------------------------"<<endl; 
    for(int k=0;k<newMeshFaces_.size();k++)
    {
        if(faceToCells_[k].size() != 0)
        {
            Info<<"Face: -"<<k<<"-";
            for(int j=0;j<newMeshFaces_[k].size();j++)
            {
                Info<<newMeshPoints_[newMeshFaces_[k][j]]<<"-";
            }
            if(faceToCells_[k].size() == 1)
                Info<<" cuts cell"<<faceToCells_[k][0];
            else if(faceToCells_[k].size() > 1)
            {
                Info<<" is cut face and neighbors cell:";
                for(int l=0;l<faceToCells_[k].size();l++)
                {
                    Info<<" "<<faceToCells_[k][l];
                }
            }
            Info<<endl;
        }
    } 
    
    Info<<"---------------------------------------CelltoFace------------------------------------"<<endl;
    for(int k=0;k<cellToFaces_.size();k++)
    {
        Info<<"Cell "<<k<<" has added faces: ";
        for(int j=0;j<cellToFaces_[k].size();j++)
        {
            int index = cellToFaces_[k][j];
            Info<<"-"<<index;
        }
        Info<<endl;
    }
    Info<<"---------------------------------------FacetoEdge------------------------------------"<<endl;
    for(int k=0;k<faceToCells_.size();k++)
    {
        if(faceToCells_[k].size() == 0)
            continue;
        Info<<"Face "<<k;
        labelList indexes = faceToCells_[k];
        if(indexes.size() == 1)
        {
            Info<<" cuts cell: " <<indexes[0];
        }
        else
        {
            
            Info<<" is cut face and neighbors cells: ";
            for(int j=0;j<indexes.size();j++)
            {
                Info<<indexes[j]<<"-";
            }
        }
        Info<<endl;
    }
    
    Info<<"--------------------------------------FacesToSide-----------------------------------"<<endl;
    for(int k=0;k<newMeshFaces_.size();k++)
    {
        Info<<"Face "<<k<<" "<<" Side: "<<facesToSide_[k];
        Info<<endl;
    } 
}

void Foam::cutCellPolyMesh::cutOldFaces
(
)
{
    //const cellList& meshCells = this->cells();
    const faceList& meshFaces = this->faces();
    //const edgeList& meshEdges = this->edges();
    //const pointField& meshPoints = this->points();
    
    cellsToSide();  
    
/*    
    edgeList        cutEdges;
    labelListList   oldEdgesToCutEdges(meshEdges.size());
    labelList       cutEdgesToSide;
    
    for(int i=0;i<meshEdges.size();i++)
    {
        if(addedStruc.oldEdgesToAddedPoint[i] != -1)
        {
            scalar startPointInd = meshEdges[i].start();
            scalar startPointSide = addedStruc.oldPointsAtCutCellsToSide[startPointInd];
            scalar endPointInd = meshEdges[i].end();
            scalar endPointSide = addedStruc.oldPointsAtCutCellsToSide[endPointInd];
            scalar cutPointInd = addedStruc.oldEdgesToAddedPoint[i] + meshPoints.size();
            
            edge newEdgeS_M(startPointInd,cutPointInd);
            cutEdges.append(newEdgeS_M);
            oldEdgesToCutEdges[i].append(cutEdges.size()-1);
            cutEdgesToSide.append(startPointSide);
            
            edge newEdgeM_E(cutPointInd,endPointInd);
            cutEdges.append(newEdgeM_E);
            oldEdgesToCutEdges[i].append(cutEdges.size()-1);
            cutEdgesToSide.append(endPointSide);            
        }
    }
    
    cutCellsItems.cutEdges = cutEdges;
    cutCellsItems.oldEdgesToCutEdges = oldEdgesToCutEdges;
    cutCellsItems.cutEdgesToSide = cutEdgesToSide;
*/
    
    cutFaces_ = faceList(0);
    oldFacesToCutFaces_ = labelListList(meshFaces.size());
    cutFacesToSide_ = labelList(0);
    
    for(int i=0;i<meshFaces.size();i++)
    {
        if(faceToEdges_[i].size() == 1 && faceToEdges_[i][0] >= nbrOfPrevEdges)
        {
            edge        addedEdge = newMeshEdges_[faceToEdges_[i][0]];
            face        currFace = meshFaces[i];
            labelList   newFace1(0);
            scalar      newFace1Sign;
            labelList   newFace2(0);
            scalar      newFace2Sign;
            label       faceInd = 0;
            label       currPointIndex = 0;
            label       relPointIndex = 0; // zero at the first appended point
            label       firstFacePoint = currFace[faceInd];
            label       currPoint = firstFacePoint;
            label       firstCutPoint;
            label       secondCutPoint;
            
            Info<<endl<<"OrigFace "<<i<<" size:"<<currFace.size()<<"| ";
            for(int l=0;l<currFace.size();l++)
            {                
                label point = currFace[l];
                Info<<point<<newMeshPoints_[point]<<" ";
            }
            Info<<"Added Edge: "<<newMeshPoints_[addedEdge.start()]<<"->"<<newMeshPoints_[addedEdge.end()];
            Info<<endl;
            
            // Cycle to first non cut point
            while(pointsToSide_[currPoint]==0)
            {
                currPoint = currFace.nextLabel(currPointIndex);
                currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                relPointIndex++;
            }
            
            // Append fist point to first face
            newFace1.append(currPoint);
            newFace1Sign = pointsToSide_[currPoint];
            
            
                Info<<"1 Face 1 size:"<<newFace1.size()<<"| ";
                for(int l=0;l<newFace1.size();l++)
                {                
                    label point = newFace1[l];
                    Info<<point<<newMeshPoints_[point]<<" ";
                }
                Info<<relPointIndex<<endl;
            
            
            //Info<<"First part of face 1"<<endl;
            while(pointsToSide_[currPoint] == pointsToSide_[currFace.nextLabel(currPointIndex)])
            {
                currPoint = currFace.nextLabel(currPointIndex);
                currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                relPointIndex++;
                newFace1.append(currPoint);
    
                    Info<<"2 Face 1 size:"<<newFace1.size()<<"| ";
                    for(int l=0;l<newFace1.size();l++)
                    {                
                        label point = newFace1[l];
                        Info<<point<<newMeshPoints_[point]<<" ";
                    }
                    Info<<relPointIndex<<endl;
            }
//////////////////////////////////////////////////////////////////////////////////
//      Fallunterscheidung jump point old or addedCutFaceNeighbor
//////////////////////////////////////////////////////////////////////////////////
            label nextPoint = -1;
            
            // Cut point is old point
            if(pointsToSide_[currFace.nextLabel(currPointIndex)] == 0)
            {
                //Info<<newMeshPoints_[currFace[faceInd+1]];
                firstCutPoint = currFace.nextLabel(currPointIndex);
                if(relPointIndex+2 >= currFace.size())
                {
                    FatalErrorInFunction
                    << "First cut point is the last point of face. "
                    << "THis can not happen"
                    << abort(FatalError);
                }
                
                nextPoint = currFace.nextLabel((currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1);
                
                currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                relPointIndex++;
                
                secondCutPoint = -1;
                if(addedEdge.start() == firstCutPoint)
                    secondCutPoint = addedEdge.end();
                else if(addedEdge.end() == firstCutPoint)
                    secondCutPoint = addedEdge.start();
                else
                {
                    FatalErrorInFunction
                    << "Error: No Point matches in CutFace"
                    << abort(FatalError);
                }                
            }
            // Cut point is added point
            else if(pointsToSide_[currPoint] * pointsToSide_[currFace.nextLabel(currPointIndex)] < 0)
            {
                nextPoint = currFace.nextLabel(currPointIndex);
                labelList currPointEdges = this->pointEdges()[currPoint];
                labelList nextPointEdges = this->pointEdges()[nextPoint];
                /*
// Zusatzinformation Anfang
                for(int l=0;l<4;l++)
                {                
                    label point = currFace[l];
                    Info<<"Point "<<point<<"-"<<meshPoints[point];
                }
                Info<<endl;
            
                Info<<"currPoint "<<currPoint<<"-"<<meshPoints[currPoint]<<endl;
                Info<<"nextPoint "<<nextPoint<<"-"<<meshPoints[nextPoint]<<endl;
                
                Info<<"currPointEdges:";
                for(int k=0;k<currPointEdges.size();k++)
                {
                    Info<<currPointEdges[k]<<" ";
                }
                Info<<endl;
                Info<<"nextPointEdges:";
                for(int j=0;j<nextPointEdges.size();j++)
                {
                    Info<<nextPointEdges[j]<<" ";
                }
                Info<<endl;
//Zusatzinformation Ende
                */
                
                label sharedEdge = -1;
                for(int k=0;k<currPointEdges.size();k++)
                {
                    for(int j=0;j<nextPointEdges.size();j++)
                    {
                        if(currPointEdges[k] == nextPointEdges[j])
                            sharedEdge = currPointEdges[k];
                    }
                }
                //Info<<"SharedEdge: "<<sharedEdge<<endl;
            
                if(sharedEdge == -1)
                {
                    FatalErrorInFunction
                    << "Error: No shared Edge found in Cut Cell Method"
                    << abort(FatalError);
                }
                firstCutPoint = edgeToPoint_[sharedEdge];
                //Info<<firstCutPoint<<endl;
                
                secondCutPoint = -1;
                if(addedEdge.start() == firstCutPoint)
                    secondCutPoint = addedEdge.end();
                else if(addedEdge.end() == firstCutPoint)
                    secondCutPoint = addedEdge.start();
                else
                {
                    FatalErrorInFunction
                    << "Error: No Point matches in CutFace"
                    << abort(FatalError);
                }
                //Info<<firstCutPoint<<endl;
            }
            else
            {
                FatalErrorInFunction
                << "Next point is neither a cut point nor a point at the " 
                << "other side of the cut."
                << abort(FatalError);
            }
            
            //Info<<firstCutPoint<<endl;

            //Info<<"Add cut points"<<endl;
            newFace1.append(firstCutPoint);
            newFace2.append(firstCutPoint);
            
            
            Info<<"3 Face 1 size:"<<newFace1.size()<<"| ";
            for(int l=0;l<newFace1.size();l++)
            {                
                label point = newFace1[l];
                Info<<point<<newMeshPoints_[point]<<" ";
            }
            Info<<relPointIndex<<endl;
            
            
            currPoint = nextPoint;
            currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
            relPointIndex++;
            
            newFace2.append(currPoint);
            newFace2Sign = pointsToSide_[currPoint];
            
            
            Info<<"1 Face 2 size:"<<newFace2.size()<<"| ";
            for(int l=0;l<newFace2.size();l++)
            {                
                label point = newFace2[l];
                Info<<point<<newMeshPoints_[point]<<" ";
            }
            Info<<relPointIndex<<endl;
            
            //Info<<"Create face 2"<<endl;
            while(pointsToSide_[currPoint] == pointsToSide_[currFace.nextLabel(currPointIndex)])
            {
                Info<<newMeshPoints_[currPoint]<<pointsToSide_[currPoint]<<endl;
                Info<<newMeshPoints_[currFace.nextLabel(currPointIndex)]<<pointsToSide_[currFace.nextLabel(currPointIndex)]<<endl;
                //Info<<"Index:"<<faceInd<<endl;
                currPoint = currFace.nextLabel(currPointIndex);
                currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                relPointIndex++;
                newFace2.append(currPoint);
                
                
                    Info<<"2 Face 2 size:"<<newFace2.size()<<"| ";
                    for(int l=0;l<newFace2.size();l++)
                    {                
                        label point = newFace2[l];
                        Info<<point<<newMeshPoints_[point]<<" ";
                    }
                    Info<<relPointIndex<<endl;
                
                /*
                Info<<endl;
                Info<<currPoint<<":";
                Info<<"currPointSign:"<<addedStruc.oldPointsAtCutCellsToSide[currPoint]<<endl;
                Info<<faceInd+1<<"-"<<currFace[faceInd+1]<<":";
                Info<<"nextPointSign:"<<addedStruc.oldPointsAtCutCellsToSide[currFace[faceInd+1]]<<endl;
                */
            }
            //Info<<"Finished loop"<<endl;
            newFace2.append(secondCutPoint);
            newFace1.append(secondCutPoint);
            
            if(secondCutPoint < nbrOfPrevPoints)
            {
                currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                relPointIndex++;
            }
            
                Info<<"3 Face 2 size:"<<newFace2.size()<<"| ";
                for(int l=0;l<newFace2.size();l++)
                {                
                    label point = newFace2[l];
                    Info<<point<<newMeshPoints_[point]<<" ";
                }
                Info<<relPointIndex<<endl;
                Info<<"4 Face 1 size:"<<newFace1.size()<<"| ";
                for(int l=0;l<newFace1.size();l++)
                {                
                    label point = newFace1[l];
                    Info<<point<<newMeshPoints_[point]<<" ";
                }
                Info<<relPointIndex<<endl;

            
            Info<<"Finish face 1"<<endl;
            while(relPointIndex+1 < meshFaces[i].size())
            {
                currPoint = currFace.nextLabel(currPointIndex);
                currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                relPointIndex++;
                newFace1.append(currPoint);
                
                if(newFace1[0] == currPoint)
                {
                    FatalErrorInFunction
                    << "The starting point is added again. Something is wrong here."
                    << abort(FatalError);
                }                    
            }
            
            
            Info<<"5 Face 1 size:"<<newFace1.size()<<"| ";
            for(int l=0;l<newFace1.size();l++)
            {                
                label point = newFace1[l];
                Info<<point<<newMeshPoints_[point]<<" ";
            }
            Info<<relPointIndex<<endl;
            
            
            cutFaces_.append(face(newFace1));
            cutFacesToSide_.append(newFace1Sign);
            oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
            cutFaces_.append(face(newFace2));
            cutFacesToSide_.append(newFace2Sign);
            oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
            
            Info<<"6 Face 1 size:"<<newFace1.size()<<"| ";
            for(int l=0;l<newFace1.size();l++)
            {                
                label point = newFace1[l];
                Info<<point<<newMeshPoints_[point]<<" ";
            }
            Info<<relPointIndex<<endl;
            
            Info<<"4 Face 2 size:"<<newFace2.size()<<"| ";
            for(int l=0;l<newFace2.size();l++)
            {                
                label point = newFace2[l];
                Info<<point<<newMeshPoints_[point]<<" ";
            }
            Info<<relPointIndex<<endl;
        }
    }

/*
    pointField combinedPoints;
    combinedPoints.append(meshPoints);
    combinedPoints.append(addedStruc.addedPoints);

    Info<<"---------------------------------------CutEdges-----------------------------------"<<endl;
    for(int k=0;k<oldEdgesToCutEdges.size();k++)
    {
        Info<<"Edge "<<k<<" from "<<combinedPoints[meshEdges[k].start()]<<" to "<<
            combinedPoints[meshEdges[k].end()];
        if(oldEdgesToCutEdges[k].size()==0)
            Info<<" not cut"<<endl;
        else
        {
            Info<<" cut ";
            label firstind = oldEdgesToCutEdges[k][0];
            label secind = oldEdgesToCutEdges[k][1];
            edge edge1 = cutEdges[firstind];
            label sign1 = cutEdgesToSide[firstind];
            edge edge2 = cutEdges[secind];
            label sign2 = cutEdgesToSide[secind];
            Info<<sign1<<"-"<<combinedPoints[edge1.start()]<<"->"<<combinedPoints[edge1.end()]<<" ";
            Info<<sign2<<"-"<<combinedPoints[edge2.start()]<<"->"<<combinedPoints[edge2.end()]<<endl;
            
        }
    }
    
    Info<<"---------------------------------------CutFaces------------------------------------"<<endl;
    for(int k=0;k<oldFacesToCutFaces.size();k++)
    {
        Info<<"Face "<<k<<" from ";        
        pointField facePoints = meshFaces[k].points(combinedPoints);
        for(int j=0;j<facePoints.size();j++)
            Info<<facePoints[j]<<"->";
        
        if(oldFacesToCutFaces[k].size()==0)
            Info<<"\t"<<"not cut"<<endl;
        else
        {            
            Info<<"\t"<<"cut"<<endl;
            label firstind = oldFacesToCutFaces[k][0];
            label secind = oldFacesToCutFaces[k][1];
            label sign1 = cutFacesToSide[firstind];
            label sign2 = cutFacesToSide[secind];
            pointField facePoints1 = cutFaces[firstind].points(combinedPoints);
            pointField facePoints2 = cutFaces[secind].points(combinedPoints);
            
            Info<<"\t"<<sign1<<" ";
            for(int j=0;j<facePoints1.size();j++)
                Info<<facePoints1[j]<<"->";
            Info<<endl;
            
            Info<<"\t"<<sign2<<" ";
            for(int j=0;j<facePoints2.size();j++)
                Info<<facePoints2[j]<<"->";
            Info<<endl;            
        }
    }
    */
    
}

void Foam::cutCellPolyMesh::printCutFaces
(
)
{
    Info<<"------------------------------------OldFacetoCutFace---------------------------------"<<endl;
    for(int k=0;k<oldFacesToCutFaces_.size();k++)
    {
        Info<<"Face "<<k;
        labelList cutFaces = oldFacesToCutFaces_[k];
        if(cutFaces.size() == 2)
        {
            Info<<" is cut into"<<endl;
            Info<<cutFacesToSide_[cutFaces[0]]<<"\t";
            for(int j=0;j<cutFaces_[cutFaces[0]].size();j++)
            {
                Info<<newMeshPoints_[cutFaces_[cutFaces[0]][j]]<<"-";
            }
            Info<<endl;
            Info<<cutFacesToSide_[cutFaces[1]]<<"\t";
            for(int j=0;j<cutFaces_[cutFaces[1]].size();j++)
            {
                Info<<newMeshPoints_[cutFaces_[cutFaces[1]][j]]<<"-";
            }
            Info<<endl;
        }
        else
        {
            Info<<" is not splitted"<<endl;
        }
    }
    
    Info<<"-------------------------------------CellsToSide-----------------------------------"<<endl;
    for(int k=0;k<cellsToSide_.size();k++)
    {
        Info<<"Cell "<<k<<" "<<" Side: "<<cellsToSide_[k];
        Info<<endl;
    } 
}

void Foam::cutCellPolyMesh::createNewMeshData
(
)
{
    const cellList& meshCells = this->cells();
    const faceList& meshFaces = this->faces();
    //const edgeList& meshEdges = this->edges();
    const pointField& meshPoints = this->points();
    const labelList owner   = this->faceOwner();
    const labelList neighbour = this->faceNeighbour();    
    const polyBoundaryMesh& boundMesh = this->boundaryMesh();
    
    // Store old boundary patches
    patchStarts = labelList(boundMesh.size());
    patchSizes = labelList(boundMesh.size());
    for(int i=0;i<boundMesh.size();i++)
    {
        patchStarts[i] = boundMesh[i].start();
        patchSizes[i] = boundMesh[i].faceCentres().size();
    }
    
    Info<<"------------------------------------OldFaces------------------------------------"<<endl;
    for(int k=0;k<meshFaces.size();k++)
    {
        Info<<"Face: -"<<k<<"-";
        pointField facePoints = meshFaces[k].points(meshPoints);
        for(int j=0;j<facePoints.size();j++)
            Info<<facePoints[j]<<"->";
        Info<<"owner:"<<owner[k];
        //Info<<" neighbour:"<<neighbour[k];
        Info<<endl;
    }
    
    for(int i=0;i<boundMesh.size();i++)
    {        
        Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
    }
    Info<<"------------------------------------EndOldFaces---------------------------------"<<endl;
    
    oldSplittedCellToNewPlusCell = labelList(meshCells.size());
    oldSplittedCellToNewMinusCell = labelList(meshCells.size());
    for(int i=0;i<meshCells.size();i++)
    {
        oldSplittedCellToNewPlusCell[i] = -1;
        oldSplittedCellToNewMinusCell[i] = -1;
    }
    
    // Compute new cellIndexes for added cells
    labelList oldCellsToAddedMinusSideCellIndex(meshCells.size());
    label addedCellIndex = 0;
    for(int i=0;i<meshCells.size();i++)
    {
        if(cellToFaces_[i].size() == 1 && cellToFaces_[i][0] >= nbrOfPrevFaces)
        {
            oldCellsToAddedMinusSideCellIndex[i] = addedCellIndex+oldCellsToAddedMinusSideCellIndex.size();
            oldSplittedCellToNewPlusCell[i] = i;
            oldSplittedCellToNewMinusCell[i]= oldCellsToAddedMinusSideCellIndex[i];
            Info<<i<<"->"<<oldCellsToAddedMinusSideCellIndex[i]<<endl;
            addedCellIndex++;
        }
    }
    
/*
// Problem here
    pointField facePoints = addedStruc.addedFaces[addedStruc.oldCellsToAddedFace[0]].points(combinedPoints);
    for(int j=0;j<facePoints.size();j++)
        Info<<facePoints[j]<<"->";
    Info<<endl;
    
    facePoints = addedStruc.addedFaces[0].points(combinedPoints);
    for(int j=0;j<facePoints.size();j++)
        Info<<facePoints[j]<<"->";
    Info<<endl;
*/
    
    Info<<"Insert Split cell faces"<<endl;
    // Compute List of new faces splitting old cells
    label addedCutFacesNbr = 0;
    addedCutFaces = faceList(0);
    addedCutFaceNeighbor = labelList(0);
    addedCutFaceOwner = labelList(0);
    for(int i=0;i<cellToFaces_.size();i++)
    {
        if(cellToFaces_[i].size() == 1 && cellToFaces_[i][0] >= nbrOfPrevFaces)
        {
            face addedFace = newMeshFaces_[cellToFaces_[i][0]];
            
            labelList thisCellPointLabels = meshCells[i].labels(meshFaces);
            cell thisCell = meshCells[i];
            vector thisNormal = addedFace.normal(newMeshPoints_);
            Info<<"This Normal: "<<thisNormal<<endl;
            point thisCentre = addedFace.centre(newMeshPoints_);
            Info<<"This Centre: "<<thisCentre<<endl;
            
            label testInd = -1;
            for(int i=0;i<thisCellPointLabels.size();i++)
            {
                if(pointsToSide_[thisCellPointLabels[i]] == -1)
                {
                    testInd = thisCellPointLabels[i];
                    break;
                }
            }
            Info<<"test Point:"<<newMeshPoints_[testInd]<<endl;
            vector centreToPointInd = newMeshPoints_[testInd] - thisCentre;
            //centreToPointInd -= thisCentre;
            Info<<"centreToPointInd: "<<centreToPointInd<<endl;
            scalar dir = centreToPointInd && thisNormal;
            Info<<"dir: "<<dir<<endl;
            if(dir < 0)
                addedFace = addedFace.reverseFace();
            
            Info<<centreToPointInd<<endl;
            
            addedCutFaces.append(addedFace);
            addedCutFaceNeighbor.append(oldCellsToAddedMinusSideCellIndex[i]);
            addedCutFaceOwner.append(i);

            /*
            Info<<"+New: ";
            pointField facePoints = addedFace.points(combinedPoints);
            for(int j=0;j<facePoints.size();j++)
                Info<<facePoints[j]<<"->";
            Info<<"owner:"<<i;
            Info<<" neighbour:"<<oldCellsToAddedMinusSideCellIndex[i]<<endl;
            */
            
            addedCutFacesNbr++; 
        }
    }

    
    Info<<"Insert split faces interior"<<endl;
    // Compute the List of new faces resulting from the splitting of old faces
    label addedSplitCellsInteriorNbr = 0;
    splitAndUnsplitFacesInterior = faceList(0);
    splitAndUnsplitFaceInteriorNeighbor = labelList(0);
    splitAndUnsplitFaceInteriorOwner = labelList(0);    
    for(int i=0;i<neighbour.size();i++)
    {
        /*
        Info<<"-Old: ";
        pointField facePoints = meshFaces[i].points(combinedPoints);
        for(int j=0;j<facePoints.size();j++)
            Info<<facePoints[j]<<"->";
        Info<<"owner:"<<owner[i];
        Info<<" neighbour:"<<neighbour[i]<<endl;
        */
        
        if(faceToEdges_[i].size() == 1 && faceToEdges_[i][0] >= nbrOfPrevEdges)
        {
            if(oldFacesToCutFaces_[i].size() != 2)
            {
                FatalErrorInFunction
                << " Splitted interior cell is cut into"<<oldFacesToCutFaces_[i].size()
                << " faces instead of the expected 2."
                << abort(FatalError);
            }
            face face1      = cutFaces_[oldFacesToCutFaces_[i][0]];
            label signFace1 = cutFacesToSide_[oldFacesToCutFaces_[i][0]];
            face face2      = cutFaces_[oldFacesToCutFaces_[i][1]];
            //label signFace2 = cutFacesToSide_[oldFacesToCutFaces_[i][1]];
            
            face sameCellFace;
            face addedCellFace;
            if(signFace1 > 0)
            {
                sameCellFace = face1;
                addedCellFace = face2;
            }
            else
            {
                sameCellFace = face2;
                addedCellFace = face1;
            }
            splitAndUnsplitFacesInterior.append(sameCellFace);
            splitAndUnsplitFaceInteriorNeighbor.append(neighbour[i]);
            splitAndUnsplitFaceInteriorOwner.append(owner[i]);

            /*
            Info<<"+New: ";
            pointField facePoints = sameCellFace.points(combinedPoints);
            for(int j=0;j<facePoints.size();j++)
                Info<<facePoints[j]<<"->";
            Info<<"owner:"<<owner[i];
            Info<<" neighbour:"<<neighbour[i]<<endl;
            */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            splitAndUnsplitFacesInterior.append(addedCellFace);
            splitAndUnsplitFaceInteriorNeighbor.append(oldCellsToAddedMinusSideCellIndex[neighbour[i]]);
            splitAndUnsplitFaceInteriorOwner.append(oldCellsToAddedMinusSideCellIndex[owner[i]]);

            /*
            Info<<"+New: ";
            facePoints = addedCellFace.points(combinedPoints);
            for(int j=0;j<facePoints.size();j++)
                Info<<facePoints[j]<<"->";
            Info<<"owner:"<<oldCellsToAddedMinusSideCellIndex[owner[i]];
            Info<<" neighbour:"<<oldCellsToAddedMinusSideCellIndex[neighbour[i]]<<endl;
            */
            
            // -1 face + 2 faces -> +1 face
            addedSplitCellsInteriorNbr++;
        }
        else
        {
            splitAndUnsplitFacesInterior.append(meshFaces[i]);
            if(facesToSide_[i] == 1)
            {
                splitAndUnsplitFaceInteriorNeighbor.append(neighbour[i]);
                splitAndUnsplitFaceInteriorOwner.append(owner[i]);
            }
            else if(facesToSide_[i] == -1)
            {                
                if(cellsToSide_[neighbour[i]] == 0)
                    splitAndUnsplitFaceInteriorNeighbor.append(oldCellsToAddedMinusSideCellIndex[neighbour[i]]);
                else
                    splitAndUnsplitFaceInteriorNeighbor.append(neighbour[i]);
                
                if(cellsToSide_[owner[i]] == 0)
                    splitAndUnsplitFaceInteriorOwner.append(oldCellsToAddedMinusSideCellIndex[owner[i]]);
                else
                    splitAndUnsplitFaceInteriorOwner.append(owner[i]);                
            }
            else
            {
                FatalErrorInFunction
                << "A face with the side: "<<facesToSide_[i]<<" was not treated."
                << " This must not happen."
                << abort(FatalError);
            }
            
            /*
            // Modify !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Info<<"+New: ";
            pointField facePoints = meshFaces[i].points(combinedPoints);
            for(int j=0;j<facePoints.size();j++)
                Info<<facePoints[j]<<"->";
            Info<<"owner:"<<owner[i];
            Info<<" neighbour:"<<neighbour[i]<<endl;
            */
        }
        
        if(splitAndUnsplitFaceInteriorOwner[splitAndUnsplitFaceInteriorOwner.size()-1] == -1)
        {
            FatalErrorInFunction
            << " Owner of face must not be -1 as happend in face "<<i
            << abort(FatalError);
        }
    }

    
    for(int i=0;i<boundMesh.size();i++)
    {
        patchStarts[i] += addedSplitCellsInteriorNbr+addedCutFacesNbr;
    }
    
    
    Info<<"Insert split faces boundary"<<endl;
    label currBoundaryPatch = 0;
    label countOldBoundaryFaces = 0;
    label countNewBoundaryFaces = 0;
    splitAndUnsplitFacesBoundary = faceList(0);
    splitAndUnsplitFaceBoundaryNeighbor = labelList(0);
    splitAndUnsplitFaceBoundaryOwner = labelList(0);
    for(int i=neighbour.size();i<meshFaces.size();i++)
    {
        /*
        Info<<"-Old: ";
        pointField facePoints = meshFaces[i].points(combinedPoints);
        for(int j=0;j<facePoints.size();j++)
            Info<<facePoints[j]<<"->";
        Info<<"owner:"<<owner[i];
        Info<<" neighbour:"<<-1<<endl;
        */
        
        if(faceToEdges_[i].size() == 1 && faceToEdges_[i][0] >= nbrOfPrevEdges)
        {
            face face1      = cutFaces_[oldFacesToCutFaces_[i][0]];
            label signFace1 = cutFacesToSide_[oldFacesToCutFaces_[i][0]];
            face face2      = cutFaces_[oldFacesToCutFaces_[i][1]];
            //label signFace2 = cutFacesToSide_[oldFacesToCutFaces_[i][1]];
            
            face sameCellFace;
            face addedCellFace;
            if(signFace1 > 0)
            {
                sameCellFace = face1;
                addedCellFace = face2;
            }
            else
            {
                sameCellFace = face2;
                addedCellFace = face1;
            }
            splitAndUnsplitFacesBoundary.append(sameCellFace);
            splitAndUnsplitFaceBoundaryNeighbor.append(-1);
            splitAndUnsplitFaceBoundaryOwner.append(owner[i]);

            /*
            Info<<"+New: ";
            pointField facePoints = sameCellFace.points(combinedPoints);
            for(int j=0;j<facePoints.size();j++)
                Info<<facePoints[j]<<"->";
            Info<<"owner:"<<owner[i];
            Info<<" neighbour:"<<-1<<endl;
            */
            
            splitAndUnsplitFacesBoundary.append(addedCellFace);
            splitAndUnsplitFaceBoundaryNeighbor.append(-1);
            splitAndUnsplitFaceBoundaryOwner.append(oldCellsToAddedMinusSideCellIndex[owner[i]]);
            
            /*
            Info<<"+New: ";
            facePoints = addedCellFace.points(combinedPoints);
            for(int j=0;j<facePoints.size();j++)
                Info<<facePoints[j]<<"->";
            Info<<"owner:"<<oldCellsToAddedMinusSideCellIndex[owner[i]];
            Info<<" neighbour:"<<-1<<endl;
            */
            
            countNewBoundaryFaces += 2;
        }
        else
        {
            splitAndUnsplitFacesBoundary.append(meshFaces[i]);
            
            if(facesToSide_[i] == 1)
            {
                splitAndUnsplitFaceBoundaryNeighbor.append(-1);
                splitAndUnsplitFaceBoundaryOwner.append(owner[i]);
            }
            else if(facesToSide_[i] == -1)
            {
                splitAndUnsplitFaceBoundaryNeighbor.append(-1);
                
                if(cellsToSide_[owner[i]] == 0)
                    splitAndUnsplitFaceBoundaryOwner.append(oldCellsToAddedMinusSideCellIndex[owner[i]]);
                else
                    splitAndUnsplitFaceBoundaryOwner.append(owner[i]);                
            }
            else
            {
                FatalErrorInFunction
                << "A face with the side: "<<facesToSide_[i]<<" was not treated."
                << " This must not happen."
                << abort(FatalError);
            }
            
            /*
            // Modify !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
            Info<<"+New: ";
            pointField facePoints = meshFaces[i].points(combinedPoints);
            for(int j=0;j<facePoints.size();j++)
                Info<<facePoints[j]<<"->";
            Info<<"owner:"<<owner[i];
            Info<<" neighbour:"<<-1<<endl;
            */
            
            countNewBoundaryFaces += 1;
        }
        countOldBoundaryFaces++;
        
        if(patchSizes[currBoundaryPatch] <= countOldBoundaryFaces)
        {
            patchSizes[currBoundaryPatch] = countNewBoundaryFaces;
            currBoundaryPatch++;
            countNewBoundaryFaces = 0;
            countOldBoundaryFaces = 0;            
        }        
    }
    for(int i=1;i<patchStarts.size();i++)
    {
        Info<<i<<":"<<patchStarts[i-1]<<"+"<<patchSizes[i-1]<<"="<<patchStarts[i-1] + patchSizes[i-1]<<endl;
        patchStarts[i] = patchStarts[i-1] + patchSizes[i-1]; 
    }
    

    /*
    Info<<"---------------------------------------Faces------------------------------------"<<endl;
    Info<<"NbrFaces:"<<combinedFaces.size()<<endl;
    Info<<"NbrOwner:"<<combinedOwner.size()<<endl;
    Info<<"NbrNeighbor:"<<combinedNeighbor.size()<<endl;
    for(int k=0;k<combinedFaces.size();k++)
    {
        Info<<"Face: -"<<k<<"-";
        pointField facePoints = combinedFaces[k].points(combinedPoints);
        for(int j=0;j<facePoints.size();j++)
            Info<<facePoints[j]<<"->";
        Info<<"owner:"<<combinedOwner[k];
        Info<<" neighbour:"<<combinedNeighbor[k];
        Info<<endl;
    }
    */
    
    for(int i=0;i<boundMesh.size();i++)
    {
        Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
    }
}

void Foam::cutCellPolyMesh::createNewMeshData_cutNeg
(
)
{
    const cellList& meshCells = this->cells();
    const faceList& meshFaces = this->faces();
    //const edgeList& meshEdges = this->edges();
    const pointField& meshPoints = this->points();
    const labelList owner   = this->faceOwner();
    const labelList neighbour = this->faceNeighbour();    
    const polyBoundaryMesh& boundMesh = this->boundaryMesh();
    
    // Store old boundary patches
    patchStarts = labelList(boundMesh.size());
    patchSizes = labelList(boundMesh.size());
    for(int i=0;i<boundMesh.size();i++)
    {
        patchStarts[i] = boundMesh[i].start();
        patchSizes[i] = boundMesh[i].faceCentres().size();
    }
    
    Info<<"------------------------------------OldFaces------------------------------------"<<endl;
    for(int k=0;k<meshFaces.size();k++)
    {
        Info<<"Face: -"<<k<<"-";
        pointField facePoints = meshFaces[k].points(meshPoints);
        for(int j=0;j<facePoints.size();j++)
            Info<<facePoints[j]<<"->";
        Info<<"owner:"<<owner[k];
        //Info<<" neighbour:"<<neighbour[k];
        Info<<endl;
    }
    
    for(int i=0;i<boundMesh.size();i++)
    {        
        Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
    }
    Info<<"------------------------------------EndOldFaces---------------------------------"<<endl;
    
    oldSplittedCellToNewPlusCell = labelList(meshCells.size());
    oldSplittedCellToNewMinusCell = labelList(meshCells.size());
    for(int i=0;i<meshCells.size();i++)
    {
        oldSplittedCellToNewPlusCell[i] = -1;
        oldSplittedCellToNewMinusCell[i] = -1;
    }
    
    // Compute new cellIndexes for added cells
    labelList oldCellsToAddedMinusSideCellIndex(meshCells.size());
    label addedCellIndex = 0;
    for(int i=0;i<meshCells.size();i++)
    {
        if(cellToFaces_[i].size() == 1 && cellToFaces_[i][0] >= nbrOfPrevFaces)
        {
            oldCellsToAddedMinusSideCellIndex[i] = addedCellIndex+oldCellsToAddedMinusSideCellIndex.size();
            oldSplittedCellToNewPlusCell[i] = i;
            oldSplittedCellToNewMinusCell[i]= oldCellsToAddedMinusSideCellIndex[i];
            Info<<i<<"->"<<oldCellsToAddedMinusSideCellIndex[i]<<endl;
            addedCellIndex++;
        }
    }
    
/*
// Problem here
    pointField facePoints = addedStruc.addedFaces[addedStruc.oldCellsToAddedFace[0]].points(combinedPoints);
    for(int j=0;j<facePoints.size();j++)
        Info<<facePoints[j]<<"->";
    Info<<endl;
    
    facePoints = addedStruc.addedFaces[0].points(combinedPoints);
    for(int j=0;j<facePoints.size();j++)
        Info<<facePoints[j]<<"->";
    Info<<endl;
*/
    
    Info<<"Insert Split cell faces"<<endl;
    // Compute List of new faces splitting old cells
    //label addedCutFacesNbr = 0;
    addedCutFaces = faceList(0);
    addedCutFaceNeighbor = labelList(0);
    addedCutFaceOwner = labelList(0);
    for(int i=0;i<cellToFaces_.size();i++)
    {
        if(cellToFaces_[i].size() == 1 && cellToFaces_[i][0] >= nbrOfPrevFaces)
        {
            face addedFace = newMeshFaces_[cellToFaces_[i][0]];
            
            labelList thisCellPointLabels = meshCells[i].labels(meshFaces);
            cell thisCell = meshCells[i];
            vector thisNormal = addedFace.normal(newMeshPoints_);
            Info<<"This Normal: "<<thisNormal<<endl;
            point thisCentre = addedFace.centre(newMeshPoints_);
            Info<<"This Centre: "<<thisCentre<<endl;
            
            label testInd = -1;
            for(int i=0;i<thisCellPointLabels.size();i++)
            {
                if(pointsToSide_[thisCellPointLabels[i]] == -1)
                {
                    testInd = thisCellPointLabels[i];
                    break;
                }
            }
            Info<<"test Point:"<<newMeshPoints_[testInd]<<endl;
            vector centreToPointInd = newMeshPoints_[testInd] - thisCentre;
            //centreToPointInd -= thisCentre;
            Info<<"centreToPointInd: "<<centreToPointInd<<endl;
            scalar dir = centreToPointInd && thisNormal;
            Info<<"dir: "<<dir<<endl;
            if(dir < 0)
                addedFace = addedFace.reverseFace();
            
            Info<<centreToPointInd<<endl;
            
            addedCutFaces.append(addedFace);
            addedCutFaceNeighbor.append(-1);
            addedCutFaceOwner.append(i);

            /*
            Info<<"+New: ";
            pointField facePoints = addedFace.points(combinedPoints);
            for(int j=0;j<facePoints.size();j++)
                Info<<facePoints[j]<<"->";
            Info<<"owner:"<<i;
            Info<<" neighbour:"<<oldCellsToAddedMinusSideCellIndex[i]<<endl;
            */
            
        }
    }

    
    Info<<"Insert split faces interior"<<endl;
    // Compute the List of new faces resulting from the splitting of old faces
    label addedSplitCellsInteriorNbr = 0;
    splitAndUnsplitFacesInterior = faceList(0);
    splitAndUnsplitFaceInteriorNeighbor = labelList(0);
    splitAndUnsplitFaceInteriorOwner = labelList(0);    
    for(int i=0;i<neighbour.size();i++)
    {
        /*
        Info<<"-Old: ";
        pointField facePoints = meshFaces[i].points(combinedPoints);
        for(int j=0;j<facePoints.size();j++)
            Info<<facePoints[j]<<"->";
        Info<<"owner:"<<owner[i];
        Info<<" neighbour:"<<neighbour[i]<<endl;
        */
        
        if(faceToEdges_[i].size() == 1 && faceToEdges_[i][0] >= nbrOfPrevEdges)
        {
            if(oldFacesToCutFaces_[i].size() != 2)
            {
                FatalErrorInFunction
                << " Splitted interior cell is cut into"<<oldFacesToCutFaces_[i].size()
                << " faces instead of the expected 2."
                << abort(FatalError);
            }
            face face1      = cutFaces_[oldFacesToCutFaces_[i][0]];
            label signFace1 = cutFacesToSide_[oldFacesToCutFaces_[i][0]];
            face face2      = cutFaces_[oldFacesToCutFaces_[i][1]];
            //label signFace2 = cutFacesToSide_[oldFacesToCutFaces_[i][1]];
            
            face sameCellFace;
            face addedCellFace;
            if(signFace1 > 0)
            {
                sameCellFace = face1;
                addedCellFace = face2;
            }
            else
            {
                sameCellFace = face2;
                addedCellFace = face1;
            }
            splitAndUnsplitFacesInterior.append(sameCellFace);
            splitAndUnsplitFaceInteriorNeighbor.append(neighbour[i]);
            splitAndUnsplitFaceInteriorOwner.append(owner[i]);
            
            addedSplitCellsInteriorNbr++;
        }
        else
        {
            if(facesToSide_[i] == 1)
            {
                splitAndUnsplitFacesInterior.append(meshFaces[i]);
                splitAndUnsplitFaceInteriorNeighbor.append(neighbour[i]);
                splitAndUnsplitFaceInteriorOwner.append(owner[i]);
                
                addedSplitCellsInteriorNbr++;
            }
            else if(facesToSide_[i] != -1)
            {
                FatalErrorInFunction
                << "A face with the side: "<<facesToSide_[i]<<" was not treated."
                << " This must not happen."
                << abort(FatalError);
            }
            
        }
        
        if(splitAndUnsplitFaceInteriorOwner[splitAndUnsplitFaceInteriorOwner.size()-1] == -1)
        {
            FatalErrorInFunction
            << " Owner of face must not be -1 as happend in face "<<i
            << abort(FatalError);
        }
    }

    for(int i=0;i<boundMesh.size();i++)
    {
        patchStarts[i] += addedSplitCellsInteriorNbr-neighbour.size();
    }
    
    Info<<"Insert split faces boundary"<<endl;
    label currBoundaryPatch = 0;
    label countOldBoundaryFaces = 0;
    label countNewBoundaryFaces = 0;
    splitAndUnsplitFacesBoundary = faceList(0);
    splitAndUnsplitFaceBoundaryNeighbor = labelList(0);
    splitAndUnsplitFaceBoundaryOwner = labelList(0);
    for(int i=neighbour.size();i<meshFaces.size();i++)
    {
        Info<<"Boundary face "<<i;
        if(faceToEdges_[i].size() == 1 && faceToEdges_[i][0] >= nbrOfPrevEdges)
        {
            face face1      = cutFaces_[oldFacesToCutFaces_[i][0]];
            label signFace1 = cutFacesToSide_[oldFacesToCutFaces_[i][0]];
            face face2      = cutFaces_[oldFacesToCutFaces_[i][1]];
            //label signFace2 = cutFacesToSide_[oldFacesToCutFaces_[i][1]];
            
            face sameCellFace;
            face addedCellFace;
            if(signFace1 > 0)
            {
                sameCellFace = face1;
                addedCellFace = face2;
            }
            else
            {
                sameCellFace = face2;
                addedCellFace = face1;
            }
            splitAndUnsplitFacesBoundary.append(sameCellFace);
            splitAndUnsplitFaceBoundaryNeighbor.append(-1);
            splitAndUnsplitFaceBoundaryOwner.append(owner[i]);

            /*
            Info<<"+New: ";
            pointField facePoints = sameCellFace.points(combinedPoints);
            for(int j=0;j<facePoints.size();j++)
                Info<<facePoints[j]<<"->";
            Info<<"owner:"<<owner[i];
            Info<<" neighbour:"<<-1<<endl;
            */
            
            /*
            Info<<"+New: ";
            facePoints = addedCellFace.points(combinedPoints);
            for(int j=0;j<facePoints.size();j++)
                Info<<facePoints[j]<<"->";
            Info<<"owner:"<<oldCellsToAddedMinusSideCellIndex[owner[i]];
            Info<<" neighbour:"<<-1<<endl;
            */
            
            countNewBoundaryFaces++;
        }
        else
        {            
            if(facesToSide_[i] == 1)
            {
                countNewBoundaryFaces++;
                splitAndUnsplitFacesBoundary.append(meshFaces[i]);
                splitAndUnsplitFaceBoundaryNeighbor.append(-1);
                splitAndUnsplitFaceBoundaryOwner.append(owner[i]);
            }
            else if(facesToSide_[i] == 0)
            {
                splitAndUnsplitFacesBoundary.append(meshFaces[i]);
                splitAndUnsplitFaceBoundaryNeighbor.append(-1);
                splitAndUnsplitFaceBoundaryOwner.append(owner[i]);
            }
            else if(facesToSide_[i] != -1)
            {
                FatalErrorInFunction
                << "A face with the side: "<<facesToSide_[i]<<" was not treated."
                << " This must not happen."
                << abort(FatalError);
            }
            
            /*
            // Modify !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
            Info<<"+New: ";
            pointField facePoints = meshFaces[i].points(combinedPoints);
            for(int j=0;j<facePoints.size();j++)
                Info<<facePoints[j]<<"->";
            Info<<"owner:"<<owner[i];
            Info<<" neighbour:"<<-1<<endl;
            */
            
        }
        countOldBoundaryFaces++;
        
        if(patchSizes[currBoundaryPatch] <= countOldBoundaryFaces)
        {
            patchSizes[currBoundaryPatch] = countNewBoundaryFaces;
            currBoundaryPatch++;
            countNewBoundaryFaces = 0;
            countOldBoundaryFaces = 0;            
        }        
    }
    for(int i=1;i<patchStarts.size();i++)
    {
        Info<<i<<":"<<patchStarts[i-1]<<"+"<<patchSizes[i-1]<<"="<<patchStarts[i-1] + patchSizes[i-1]<<endl;
        patchStarts[i] = patchStarts[i-1] + patchSizes[i-1]; 
    }
    

    /*
    Info<<"---------------------------------------Faces------------------------------------"<<endl;
    Info<<"NbrFaces:"<<combinedFaces.size()<<endl;
    Info<<"NbrOwner:"<<combinedOwner.size()<<endl;
    Info<<"NbrNeighbor:"<<combinedNeighbor.size()<<endl;
    for(int k=0;k<combinedFaces.size();k++)
    {
        Info<<"Face: -"<<k<<"-";
        pointField facePoints = combinedFaces[k].points(combinedPoints);
        for(int j=0;j<facePoints.size();j++)
            Info<<facePoints[j]<<"->";
        Info<<"owner:"<<combinedOwner[k];
        Info<<" neighbour:"<<combinedNeighbor[k];
        Info<<endl;
    }
    */
    
    /*
    word boundName = "Cut Bound";
    labelList patches(addedCutFaces.size());
    label nbrFacesBeforeCutFaces =  splitAndUnsplitFacesInterior.size() 
                                    + splitAndUnsplitFacesBoundary.size();
    for(int i=0;i<patches.size();i++)
    {
        patches[i] = nbrFacesBeforeCutFaces + i;
    }
    boundMesh.setGroup(boundName,patches);
    */
    
    for(int i=0;i<boundMesh.size();i++)
    {
        Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
    }
    
    /*
    labelList posPoints(0);
    for(int i=0;i<pointsToSide_.size();i++)
    {
        if(pointsToSide_[i] >= 0)
            posPoints.append(
    }
    */
}

void Foam::cutCellPolyMesh::printNewMeshData
(
)
{
    const polyBoundaryMesh& boundMesh = this->boundaryMesh();
    const wordList& types = boundMesh.types();
    const wordList& phystypes = boundMesh.physicalTypes();
    label nextBound = patchStarts[0];
    label countFaces = 0;
    label boundaryIndex = -1;
    Info<<"---------------------------------------CutFaces------------------------------------"<<endl;
    for(int i=0;i<addedCutFaces.size();i++)
    {
        if(countFaces == nextBound)
        {
            boundaryIndex++;
            Info<<"Boundarypatch "<<boundaryIndex<<" patchStart "<<patchStarts[boundaryIndex]<<" patchSize "<<patchSizes[boundaryIndex]<<"  type:"<<types[boundaryIndex]<<"  phystype:"<<phystypes[boundaryIndex]<<endl;
            if(boundaryIndex != patchStarts.size()-1)
                nextBound = patchStarts[boundaryIndex+1];
        }
        Info<<"  Face:"<<i<<" Owner:"<<addedCutFaceOwner[i]<<" Neighbor:"<<addedCutFaceNeighbor[i]<<" ";
        for(int k=0;k<addedCutFaces[i].size();k++)
        {
            Info<<newMeshPoints_[addedCutFaces[i][k]]<<"->";
        }
        Info<<endl;
        countFaces++;
    }
    
    Info<<"-----------------------------splitAndUnsplitFacesInterior--------------------------"<<endl;
    for(int i=0;i<splitAndUnsplitFacesInterior.size();i++)
    {
        if(countFaces == nextBound)
        {
            boundaryIndex++;
            Info<<"Boundarypatch "<<boundaryIndex<<" patchStart "<<patchStarts[boundaryIndex]<<" patchSize "<<patchSizes[boundaryIndex]<<"  type:"<<types[boundaryIndex]<<"  phystype:"<<phystypes[boundaryIndex]<<endl;
            if(boundaryIndex != patchStarts.size()-1)
                nextBound = patchStarts[boundaryIndex+1];
        }
        Info<<"  Face:"<<i<<" Owner:"<<splitAndUnsplitFaceInteriorOwner[i]<<" Neighbor:"<<splitAndUnsplitFaceInteriorNeighbor[i]<<" ";
        for(int k=0;k<splitAndUnsplitFacesInterior[i].size();k++)
        {
            Info<<newMeshPoints_[splitAndUnsplitFacesInterior[i][k]]<<"->";
        }
        Info<<endl;
        countFaces++;
    }
    
    Info<<"-----------------------------splitAndUnsplitFacesBoundary--------------------------"<<endl;
    for(int i=0;i<splitAndUnsplitFacesBoundary.size();i++)
    {
        if(countFaces == nextBound)
        {
            boundaryIndex++;
            Info<<"Boundarypatch "<<boundaryIndex<<" patchStart "<<patchStarts[boundaryIndex]<<" patchSize "<<patchSizes[boundaryIndex]<<"  type:"<<types[boundaryIndex]<<"  phystype:"<<phystypes[boundaryIndex]<<endl;
            if(boundaryIndex != patchStarts.size()-1)
                nextBound = patchStarts[boundaryIndex+1];
        }
        Info<<"  Face:"<<i<<" Owner:"<<splitAndUnsplitFaceBoundaryOwner[i]<<" Neighbor:"<<splitAndUnsplitFaceBoundaryNeighbor[i]<<" ";
        for(int k=0;k<splitAndUnsplitFacesBoundary[i].size();k++)
        {
            Info<<newMeshPoints_[splitAndUnsplitFacesBoundary[i][k]]<<"->";
        }        
        Info<<endl;
        countFaces++;
    }  
    
}

void Foam::cutCellPolyMesh::printMesh
(
)
{
    Info<<"------------------------------------printMesh---------------------------------"<<endl;
    const cellList& meshCells = this->cells();
    const faceList& meshFaces = this->faces();
    const edgeList& meshEdges = this->edges();
    const pointField& meshPoints = this->points();
    const labelList owner   = this->faceOwner();
    const labelList neighbour = this->faceNeighbour();
    const polyBoundaryMesh& boundMesh = this->boundaryMesh();
    
    for(int i=0;i<meshCells.size();i++)
    {
        Info<<"Cell:"<<i<<" with "<<meshCells[i].nFaces()<<" |";
        for(int k=0;k<meshCells[i].size();k++)
        {
            Info<<meshCells[i][k]<<"->";
        }
        Info<<" with centre:"<<meshCells[i].centre(meshPoints,meshFaces);
        Info<<" and volume:"<<meshCells[i].mag(meshPoints,meshFaces)<<endl;
    }

    for(int i=0;i<meshFaces.size();i++)
    {
        Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
        if(i < neighbour.size())
            Info<<" Neighbor:"<<neighbour[i]<<" ";
        for(int k=0;k<meshFaces[i].size();k++)
        {
            Info<<meshPoints[meshFaces[i][k]]<<"->";
        }
        Info<<" with centre:"<<meshFaces[i].centre(meshPoints);
        Info<<" and normal vector:"<<meshFaces[i].normal(meshPoints);
        Info<<" and area:"<<meshFaces[i].mag(meshPoints)<<endl;
    }
    
    for(int i=0;i<meshEdges.size();i++)
    {
        Info<<"Edge:"<<i<< " Start:"<<meshPoints[meshEdges[i].start()]<<"-> End:"<<meshPoints[meshEdges[i].end()]<<endl;
    }
    for(int i=0;i<boundMesh.size();i++)
    {
        Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
    }
}

vector crossProd(const vector& v1, const vector& v2)
{
    return vector(  v1[1]*v2[2]-v1[2]*v2[1],
                    v1[2]*v2[0]-v1[0]*v2[2],
                    v1[0]*v2[1]-v1[1]*v2[0]
                 );
}


void Foam::cutCellPolyMesh::selfTestMesh()
{  
    Info<<"START MESH SELF TEST";
    const cellList& meshCells = this->cells();
    const faceList& meshFaces = this->faces();
    const edgeList& meshEdges = this->edges();
    const pointField& meshPoints = this->points();
    const labelList owner   = this->faceOwner();
    const labelList neighbour = this->faceNeighbour();    
    const polyBoundaryMesh& boundMesh = this->boundaryMesh();
    
    //Test faces
    for(int i=0;i<meshFaces.size();i++)
    {
        //Info<<"Test face "<<i<<endl;
        point centreFace = meshFaces[i].centre(meshPoints);
        vector normalFace = meshFaces[i].normal(meshPoints);
        scalar area = meshFaces[i].mag(meshPoints);
        
        //Test for face shape
        point curr,prev,next;
        vector edge1,edge2;
        List<vector> crossProds(meshFaces[i].size()+1);
        for(int k=0;k<meshFaces[i].size();k++)
        {
            prev = meshPoints[meshFaces[i].prevLabel(k)];
            curr = meshPoints[meshFaces[i][k]];
            next = meshPoints[meshFaces[i].nextLabel(k)];
            
            edge1 = curr-prev;
            edge2 = next-curr;
            crossProds[k] = crossProd(edge1,edge2);
        }
        crossProds[crossProds.size()-1] = crossProds[0];
        scalar res;
        for(int k=0;k<meshFaces[i].size();k++)
        {
            res = crossProds[k] && crossProds[k+1];
            if(res<0)
            {
                Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
                if(i < neighbour.size())
                    Info<<" Neighbor:"<<neighbour[i]<<" ";
                for(int k=0;k<meshFaces[i].size();k++)
                {
                    Info<<meshPoints[meshFaces[i][k]]<<"->";
                }
                Info<<" with centre:"<<centreFace;
                Info<<" and normal vector:"<<normalFace;
                Info<<" and area:"<<area<<endl;
                
                FatalErrorInFunction
                << "Face must not have a concave shape!"
                << abort(FatalError);
            }
        }
        
        //Test for centre point internal
        vector toCentre1,toCentre2;
        crossProds = List<vector>(meshFaces[i].size()+1);
        //Info<<"---------------------------------------------"<<endl;
        for(int k=0;k<meshFaces[i].size();k++)
        {
            curr = meshPoints[meshFaces[i][k]];
            next = meshPoints[meshFaces[i].nextLabel(k)];
            
            edge1 = next-curr;
            toCentre1 = centreFace-curr;
            toCentre2 = centreFace-next;
            
            //Info<<"Edge: "<<edge1<<endl;
            //Info<<"toCentre 1: "<<toCentre1<<endl;
            //Info<<"toCentre 2: "<<toCentre2<<endl;
            
            crossProds[k] = 0.5*(crossProd(edge1,toCentre1) + crossProd(edge1,toCentre2));
            //Info<<"crossProds: "<<crossProds[k]<<endl;

        }
        crossProds[crossProds.size()-1] = crossProds[0];
        for(int k=0;k<meshFaces[i].size();k++)
        {
            res = crossProds[k] && crossProds[k+1];
            if(res<0)
            {
                Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
                if(i < neighbour.size())
                    Info<<" Neighbor:"<<neighbour[i]<<" ";
                for(int k=0;k<meshFaces[i].size();k++)
                {
                    Info<<meshPoints[meshFaces[i][k]]<<"->";
                }
                Info<<" with centre:"<<centreFace;
                Info<<" and normal vector:"<<normalFace;
                Info<<" and area:"<<area<<endl;
                
                Info<<"CrossProds k: "<<crossProds[k]<<endl;
                Info<<"CrossProds k+1: "<<crossProds[k+1]<<endl;
                Info<<"res: "<<res<<endl;
                
                
                FatalErrorInFunction
                << "Face  has a centre thats not strictly inside!"
                << abort(FatalError);
            }
        }
        
        //Test if face has double points
        //collapse not allowed in this area
        for(int j = 0;j<meshFaces[i].size();j++)
        {
            for(int k=0;k<meshFaces[i].size();k++)
            {
                if(k==j)
                    continue;
                else if(meshFaces[i][j] == meshFaces[i][k])
                {
                    Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
                    if(i < neighbour.size())
                        Info<<" Neighbor:"<<neighbour[i]<<" ";
                    for(int s=0;s<meshFaces[i].size();s++)
                    {
                        Info<<meshPoints[meshFaces[i][s]]<<"->";
                    }
                    Info<<" with centre:"<<centreFace;
                    Info<<" and normal vector:"<<normalFace;
                    Info<<" and area:"<<area<<endl;
            
                    FatalErrorInFunction
                    << "Face had a double point!"
                    << abort(FatalError);
                }
            }
        }        
        
        //Test for face with negative area
        if(area<0)
        {
            Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
            if(i < neighbour.size())
                Info<<" Neighbor:"<<neighbour[i]<<" ";
            for(int k=0;k<meshFaces[i].size();k++)
            {
                Info<<meshPoints[meshFaces[i][k]]<<"->";
            }
            Info<<" with centre:"<<centreFace;
            Info<<" and normal vector:"<<normalFace;
            Info<<" and area:"<<area<<endl;
            
            FatalErrorInFunction
            << "Face has negative area!"
            << abort(FatalError); 
        }
        
        //Test if the specified owner of each face is the actual owner
        label ownerCell = owner[i];
        bool isOwnerCell = false;
        for(int k=0;k<meshCells[ownerCell].size();k++)
        {
            if(meshCells[ownerCell][k] == i)
            {
                isOwnerCell = true;
            }
        }
        if(!isOwnerCell)
        {
            Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
            if(i < neighbour.size())
                Info<<" Neighbor:"<<neighbour[i]<<" ";
            for(int k=0;k<meshFaces[i].size();k++)
            {
                Info<<meshPoints[meshFaces[i][k]]<<"->";
            }
            Info<<" with centre:"<<centreFace;
            Info<<" and normal vector:"<<normalFace;
            Info<<" and area:"<<area<<endl;
            
            FatalErrorInFunction
            << "Is listed as owned by "<<ownerCell<<" but this cell does not have this face!"
            << abort(FatalError); 
        }
        
        //Test if the specified neighbour of each face is the actual neighbour
        if(i<neighbour.size())
        {
            label neighbourCell = neighbour[i];
            bool isNeighbourCell = false;
            for(int k=0;k<meshCells[neighbourCell].size();k++)
            {
                if(meshCells[neighbourCell][k] == i)
                {
                    isNeighbourCell = true;
                }
            }
            if(!isNeighbourCell)
            {
                Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
                if(i < neighbour.size())
                    Info<<" Neighbor:"<<neighbour[i]<<" ";
                for(int k=0;k<meshFaces[i].size();k++)
                {
                    Info<<meshPoints[meshFaces[i][k]]<<"->";
                }
                Info<<" with centre:"<<centreFace;
                Info<<" and normal vector:"<<normalFace;
                Info<<" and area:"<<area<<endl;
            
                FatalErrorInFunction
                << "Is listed as neighbouring "<<neighbourCell<<" but this cell does not have this face!"
                << abort(FatalError); 
            }
        }
        
        
        vector faceCentreToOwnerCentre = meshCells[ownerCell].centre(meshPoints,meshFaces)-centreFace;
        if((faceCentreToOwnerCentre && normalFace)>=0)
        {
            Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
            if(i < neighbour.size())
                Info<<"Neighbor:"<<neighbour[i]<<" ";
            for(int k=0;k<meshFaces[i].size();k++)
            {
                Info<<meshPoints[meshFaces[i][k]]<<"->";
            }
            Info<<" with centre:"<<centreFace;
            Info<<" and normal vector:"<<normalFace;
            Info<<" and area:"<<area<<endl;
            
            FatalErrorInFunction
            <<"Normal vector is "<<normalFace<<" while faceCentreToOwnerCentre is "<<faceCentreToOwnerCentre<<"!"
            <<" They must have a opposite direction"
            << abort(FatalError); 
        }
        
        if(i<neighbour.size())
        {
            label neighbourCell = neighbour[i];
            vector centreToNeighbourCentre = meshCells[neighbourCell].centre(meshPoints,meshFaces)-centreFace;
            if((centreToNeighbourCentre && normalFace)<=0)
            {
                Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
                if(i < neighbour.size())
                    Info<<"Neighbor:"<<neighbour[i]<<" ";
                for(int k=0;k<meshFaces[i].size();k++)
                {
                    Info<<meshPoints[meshFaces[i][k]]<<"->";
                }
                Info<<" with centre:"<<centreFace;
                Info<<" and normal vector:"<<normalFace;
                Info<<" and area:"<<area<<endl;
            
                FatalErrorInFunction
                <<"Normal vector is "<<normalFace<<" while faceCentreToNeighbourCentre is "<<centreToNeighbourCentre<<"!"
                <<" They must have the same direction"
                << abort(FatalError);
            }
        }      
    }
    
    //Test cells
    //Test if cell centre is inside cell
    for(int i=0;i<meshCells.size();i++)
    {
        const point cellCentre = meshCells[i].centre(meshPoints, meshFaces);
        scalar mag = meshCells[i].mag(meshPoints, meshFaces);
        
        //Test for correct volume
        if(mag < 0)
        {
            Info<<"Cell:"<<i<<" with "<<meshCells[i].nFaces()<<" |";
            for(int k=0;k<meshCells[i].size();k++)
            {
                Info<<meshCells[i][k]<<"->";
            }
            Info<<" with centre:"<<cellCentre;
            Info<<" and volume:"<<mag<<endl;
            
            FatalErrorInFunction
            << "Cell cannot have Volume smaller than zero! "
            << abort(FatalError);
        }
        if(mag == 0 && cellsToSide_[i] != -1)
        {
            Info<<"Cell:"<<i<<" with "<<meshCells[i].nFaces()<<"faces |";
            for(int k=0;k<meshCells[i].size();k++)
            {
                Info<<meshCells[i][k]<<"->";
            }
            Info<<" with centre:"<<cellCentre;
            Info<<" and volume:"<<mag<<endl;
            
            FatalErrorInFunction
            << "Cell cannot have Volume equal zero while being on side:"<<cellsToSide_[i]
            << abort(FatalError);
        }
        
        
        //Test if centre is really inside cell
        for(int a=0;a<meshCells[i].size();a++)
        {
            label oneFaceInd = meshCells[i][a];
            vector thisFaceNormal = meshFaces[oneFaceInd].normal(meshPoints);
            if(owner[oneFaceInd] != i)
                thisFaceNormal = -1*thisFaceNormal;
            vector thisFaceCentre = meshFaces[oneFaceInd].centre(meshPoints);
            if(((thisFaceCentre-cellCentre) && thisFaceNormal) <= 0)
            {            
                Info<<"Cell:"<<i<<" with "<<meshCells[i].nFaces()<<"faces |";
                for(int k=0;k<meshCells[i].size();k++)
                {
                    Info<<meshCells[i][k]<<"->";
                }
                Info<<" with centre:"<<cellCentre;
                Info<<" and volume:"<<mag<<endl;
                
                FatalErrorInFunction
                << "Cell Face "<<a<<" has a normal "<<thisFaceNormal<<" but cellCentreToFaceCentre is "<<thisFaceCentre-cellCentre
                << abort(FatalError);                
            }
        }
    }
    Info<<": MESH IS CORRECT"<<endl;

}

