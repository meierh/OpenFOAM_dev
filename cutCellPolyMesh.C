#include "cutCellPolyMesh.H"

Foam::cutCellPolyMesh::cutCellPolyMesh
(
    const IOobject& io,
    std::function<scalar(const vector)> levelSet
)  
:   polyMesh(io)
{
    this->levelSet = levelSet;
    newMeshPoints();
    this->levelSet = levelSet;
    pointsToSide();
}

void Foam::cutCellPolyMesh::pointsToSide()
{
    const pointField& points = this->points();
    labelList pointsToSide(points.size());
    scalar lvlSet;
    for(int i=0;i<points.size();i++)
    {
        lvlSet = levelSet(points[i]);
        if(lvlSet > 0)
            pointsToSide[i] = 1;
        else if(lvlSet < 0)
            pointsToSide[i] = -1;
        else
            pointsToSide[i] = 0;
        
        Info<<points[i]<<"\t"<<lvlSet<<"\t"<<pointsToSide[i]<<endl;
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
        lvlSet = levelSet(points[i]);
        if(lvlSet > 0)
            pointsToSide[i] = 1;
        else if(lvlSet < 0)
            pointsToSide[i] = -1;
        else
            pointsToSide[i] = 0;
        
        Info<<points[i]<<"\t"<<lvlSet<<"\t"<<pointsToSide[i]<<endl;
    }
    this->pointsToSide_ = pointsToSide;
}

void Foam::cutCellPolyMesh::edgesToSide
(
)
{
    const pointField& points = this->points();
    const edgeList& edges = this->edges();
    labelList edgesToSide(edges.size());
    
    for(int i=0;i<edges.size();i++)
    {
        point start = points[edges[i].start()];
        point end = points[edges[i].end()];
        scalar lvlSetStart = levelSet(start);
        scalar lvlSetEnd = levelSet(end);

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
    const pointField& points = this->points();
    labelList edgesToSide(edges.size());
    
    for(int i=0;i<edges.size();i++)
    {
        point start = points[edges[i].start()];
        point end = points[edges[i].end()];
        scalar lvlSetStart = levelSet(start);
        scalar lvlSetEnd = levelSet(end);

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

void Foam::cutCellPolyMesh::facesToSide()
{
    const pointField& points = this->points();
    const faceList& faces = this->faces();
    labelList facesToSide(faces.size());
    
    for(int i=0;i<faces.size();i++)
    {
        bool nullExist = false;
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<faces[i].size();k++)
        {
            point thisPoint = points[faces[i][k]];
            if(levelSet(thisPoint) == 0)
                nullExist = true;
            else if(levelSet(thisPoint) > 0)
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
    const pointField& points = this->points();
    labelList facesToSide(faces.size());
    
    for(int i=0;i<faces.size();i++)
    {
        bool nullExist = false;
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<faces[i].size();k++)
        {
            point thisPoint = points[faces[i][k]];
            if(levelSet(thisPoint) == 0)
                nullExist = true;
            else if(levelSet(thisPoint) > 0)
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
        
        bool nullExist = false;
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<cellPoints.size();k++)
        {
            point thisPoint = cellPoints[k];
            if(levelSet(thisPoint) == 0)
                nullExist = true;
            else if(levelSet(thisPoint) > 0)
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
        
        bool nullExist = false;
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<cellPoints.size();k++)
        {
            point thisPoint = cellPoints[k];
            if(levelSet(thisPoint) == 0)
                nullExist = true;
            else if(levelSet(thisPoint) > 0)
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

void Foam::cutCellPolyMesh::newMeshPoints
(
)
{
    Info<<"Starting adding Points"<<endl;
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
        scalar phiStart = levelSet(basisPoints[startLabel]);
        scalar phiEnd = levelSet(basisPoints[endLabel]);
        
        if(phiStart>0 || phiEnd>0)
            pos = +1;
        if(phiStart<0 || phiEnd<0)
            neg = -1;
        
        if(pos == +1 && neg == -1)
        {          
            vector endToStart = basisEdges[i].vec(basisPoints);
            scalar norm_endToStart = basisEdges[i].mag(basisPoints);
            scalar distPhi = std::abs(phiEnd-phiStart);
            scalar norm_phiEnd = std::abs(phiEnd);
            scalar scalePoint = norm_phiEnd * distPhi / norm_endToStart; 
            vector newPoint = basisPoints[endLabel] + scalePoint * endToStart;

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
    
    printAddedElements();
}

void Foam::cutCellPolyMesh::printAddedElements
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
    
    faceToEdges_ = labelListList(basisFaces.size());
    for(int i=0;i<basisEdges.size();i++)
    {
        if(edgeToFaces_[i].size() != 0)
        {
            for(int k=0;k<edgeToFaces_[i].size();k++)
            {
                if(faceToEdges_[edgeToFaces_[i][k]].size() == 0)
                {
                    faceToEdges_[edgeToFaces_[i][k]] = labelList(0);
                    faceToEdges_[edgeToFaces_[i][k]].append(i);
                }
                else
                {
                    faceToEdges_[edgeToFaces_[i][k]].append(i);
                }
            }
        }
    }
    
    edgeToCells_ = labelListList(newMeshEdges_.size());
    labelList owner = this->faceOwner();
    labelList neighbour = this->faceNeighbour();
    labelListList edgeCells = this->edgeCells();
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        if(edgeToFaces_[i].size() != 0)
        {
            if(i<nbrOfPrevEdges)
            {
                edgeToCells_[i] = edgeCells[i];
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
                edgeToCells_[i].append(owner[thisFace]);
                edgeToCells_[i].append(neighbour[thisFace]);
            }
        }
    }
    
    cellToEdges_ = labelListList(meshCells.size());
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        if(edgeToCells_[i].size() != 0)
        {
            for(int k=0;k<edgeToCells_[i].size();k++)
            {
                if(cellToEdges_[edgeToCells_[i][k]].size() == 0)
                {
                    cellToEdges_[edgeToCells_[i][k]] = labelList(0);
                    cellToEdges_[edgeToCells_[i][k]].append(i);
                }
                else
                {
                    cellToEdges_[edgeToCells_[i][k]].append(i);
                }
            }
        }
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
    
    faceToCells_ = labelListList(basisFaces.size());
    labelList owner = this->faceOwner();
    labelList neighbour = this->faceNeighbour();
    for(int i=0;i<basisFaces.size();i++)
    {
        labelList facePoints = basisFaces[i];
        if(facePoints.size() != 4)
        {
            FatalErrorInFunction
            << "A face cannot have "<< facePoints.size()
            << " cut points while one or more "
            << "cut points are not old points! "
            << abort(FatalError);
        }
        bool isCutFace = true;
        for(int k=0;k<facePoints.size();k++)
        {
            if(pointsToSide_[facePoints[k]] != 0)
                isCutFace = false;
        }
        if(isCutFace)
        {
            faceToCells_[i].append(owner[i]);
            faceToCells_[i].append(neighbour[i]);
        }        
    }

    cellToFace_ = labelListList(meshCells.size());
    for(int i=0;i<meshCells.size();i++)
    {     
        if(cellToEdges_[i].size() == 0)
        {
            continue;
        }
        
        labelList facePoints;
        labelList cellCutEdgeList = cellToEdges_[i];
        
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
            
            cellToFace_[i] = equalFace[0];
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
        newMeshFaces_.append(newFace);
        facesToSide_.append(0);
        
        labelList newFaceCell(0);
        newFaceCell.append(i);
        faceToCells_.append(newFaceCell);
        
        //input to oldCellsToAddedFace
        cellToFace_[i] = newMeshFaces_.size()-1;
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
        if(faceToEdges_[i].size() == 1)
        {
            edge        addedEdge = newMeshEdges_[faceToEdges_[i][0]];
            face currFace = meshFaces[i];
            labelList   newFace1(0);
            scalar      newFace1Sign;
            labelList   newFace2(0);
            scalar      newFace2Sign;
            label       faceInd = 0;
            label       firstFacePoint = currFace[faceInd];
            label       currPoint = firstFacePoint;
            
            
            Info<<endl<<"OrigFace  size:"<<currFace.size()<<"| ";
            for(int l=0;l<currFace.size();l++)
            {                
                label point = currFace[l];
                Info<<point<<" "; //<<meshPoints[point];
            }
            Info<<endl;
            
            
            newFace1.append(firstFacePoint);
            newFace1Sign = pointsToSide_[firstFacePoint];
            
            /*
            Info<<"1 Face 1 size:"<<newFace1.size()<<"| ";
            for(int l=0;l<newFace1.size();l++)
            {                
                label point = newFace1[l];
                Info<<point<<" "; //<<meshPoints[point];
            }
            Info<<endl;
            */
            
            //Info<<"First part of face 1"<<endl;
            while(pointsToSide_[currPoint] ==
                  pointsToSide_[currFace[faceInd+1]])
            {
                currPoint = currFace[faceInd+1];
                newFace1.append(currPoint);
                faceInd++;
                
                /*
                Info<<"2 Face 1 size:"<<newFace1.size()<<"| ";
                for(int l=0;l<newFace1.size();l++)
                {                
                    label point = newFace1[l];
                    Info<<point<<" "; //<<meshPoints[point];
                }
                Info<<endl;
                */
            }
            label nextPoint = currFace[faceInd+1];
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
            label firstCutPoint = edgeToPoint_[sharedEdge];
            label secondCutPoint = -1;
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

            //Info<<"Add cut points"<<endl;
            newFace1.append(firstCutPoint);
            newFace2.append(firstCutPoint);
            
            /*
            Info<<"3 Face 1 size:"<<newFace1.size()<<"| ";
            for(int l=0;l<newFace1.size();l++)
            {                
                label point = newFace1[l];
                Info<<point<<" "; //<<meshPoints[point];
            }
            Info<<endl;
            */
            
            currPoint = nextPoint;
            faceInd++;
            
            newFace2.append(currPoint);
            newFace2Sign = pointsToSide_[currPoint];
            
            /*
            Info<<"1 Face 2 size:"<<newFace2.size()<<"| ";
            for(int l=0;l<newFace2.size();l++)
            {                
                label point = newFace2[l];
                Info<<point<<" "; //<<meshPoints[point];
            }
            Info<<endl;
            */
            
            //Info<<"Create face 2"<<endl;
            while(pointsToSide_[currPoint] ==
                  pointsToSide_[currFace[faceInd+1]])
            {
                //Info<<"Index:"<<faceInd<<endl;
                currPoint = meshFaces[i][faceInd+1];
                newFace2.append(currPoint);
                faceInd++;
                
                /*
                Info<<"2 Face 2 size:"<<newFace2.size()<<"| ";
                for(int l=0;l<newFace2.size();l++)
                {                
                    label point = newFace2[l];
                    Info<<point<<" "; //<<meshPoints[point];
                }
                */
                
                if(faceInd+1>=currFace.size())
                    break;
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
            
            /*
            Info<<"3 Face 2 size:"<<newFace2.size()<<"| ";
            for(int l=0;l<newFace2.size();l++)
            {                
                label point = newFace2[l];
                Info<<point<<" "; //<<meshPoints[point];
            }
            Info<<endl;
            */
        
            /*
            Info<<"4 Face 1 size:"<<newFace1.size()<<"| ";
            for(int l=0;l<newFace1.size();l++)
            {                
                label point = newFace1[l];
                Info<<point<<" "; //<<meshPoints[point];
            }
            Info<<endl;
            */
            
            //Info<<"Finish face 1"<<endl;
            while(++faceInd < meshFaces[i].size())
            {
                currPoint = meshFaces[i][faceInd];
                //Info<<currPoint<<"!="<<firstFacePoint<<endl;
                newFace1.append(currPoint);
            }
            
            /*
            Info<<"5 Face 1 size:"<<newFace1.size()<<"| ";
            for(int l=0;l<newFace1.size();l++)
            {                
                label point = newFace1[l];
                Info<<point<<" "; //<<meshPoints[point];
            }
            Info<<endl;
            */
            
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
                Info<<point<<" "; //<<meshPoints[point];
            }
            Info<<endl;
            
            Info<<"4 Face 2 size:"<<newFace2.size()<<"| ";
            for(int l=0;l<newFace2.size();l++)
            {                
                label point = newFace2[l];
                Info<<point<<" "; //<<meshPoints[point];
            }
            Info<<endl;
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
        if(cellToFace_[i].size() == 1 && cellToFace_[i][0] >= nbrOfPrevFaces)
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
    for(int i=0;i<cellToFace_.size();i++)
    {
        if(cellToFace_[i].size() == 1 && cellToFace_[i][0] >= nbrOfPrevFaces)
        {
            face addedFace = newMeshFaces_[cellToFace_[i][0]];
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
            if(facesToSide_[i] == 1)
            {
                splitAndUnsplitFacesInterior.append(meshFaces[i]);
                splitAndUnsplitFaceInteriorNeighbor.append(neighbour[i]);
                splitAndUnsplitFaceInteriorOwner.append(owner[i]);
            }
            else
            {
                splitAndUnsplitFacesInterior.append(meshFaces[i]);
                
                if(cellsToSide_[neighbour[i]] == 0)
                    splitAndUnsplitFaceInteriorNeighbor.append(oldCellsToAddedMinusSideCellIndex[neighbour[i]]);
                else
                    splitAndUnsplitFaceInteriorNeighbor.append(neighbour[i]);
                
                if(cellsToSide_[owner[i]] == 0)
                    splitAndUnsplitFaceInteriorOwner.append(oldCellsToAddedMinusSideCellIndex[owner[i]]);
                else
                    splitAndUnsplitFaceInteriorOwner.append(owner[i]);                
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
            if(cellsToSide_[i] == 1)
            {
                splitAndUnsplitFacesBoundary.append(meshFaces[i]);
                splitAndUnsplitFaceBoundaryNeighbor.append(-1);
                splitAndUnsplitFaceBoundaryOwner.append(owner[i]);
            }
            else
            {
                splitAndUnsplitFacesBoundary.append(meshFaces[i]);
                splitAndUnsplitFaceBoundaryNeighbor.append(-1);
                
                            
                if(cellsToSide_[owner[i]] == 0)
                    splitAndUnsplitFaceBoundaryOwner.append(oldCellsToAddedMinusSideCellIndex[owner[i]]);
                else
                    splitAndUnsplitFaceBoundaryOwner.append(owner[i]);                
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
