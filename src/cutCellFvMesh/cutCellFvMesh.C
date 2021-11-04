#include "cutCellFvMesh.H"

Foam::cutCellFvMesh::cutCellFvMesh
(
    const IOobject& io,
    std::function<scalar(const vector)> levelSet,
    cutStatus state
)  
:   fvMesh(io)
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
                << " The cutCellFvMesh must be defined with an empty boundary patch at the end."
                << " Starting at: "<<this->nFaces()<<" (curr:"<<patchStarts.last()<<") and with the size: 0 ("<<patchSizes.last()<<")"
                << exit(FatalError);
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

void Foam::cutCellFvMesh::pointsToSide
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

void Foam::cutCellFvMesh::pointsToSide
(
    const pointField& points
)
{
    pointsToSide_.setSize(points.size());
    scalar lvlSet;
    for(int i=0;i<points.size();i++)
    {
        lvlSet = pointDist[i];
        pointsToSide_[i] = (lvlSet > 0) * 1 + (lvlSet < 0) * -1 + 0;
        
        /*
        if(lvlSet > 0)
            pointsToSide_[i] = 1;
        else if(lvlSet < 0)
            pointsToSide_[i] = -1;
        else
            pointsToSide_[i] = 0;
        */
        
        //Info<<points[i]<<"\t"<<lvlSet<<"\t"<<pointsToSide[i]<<endl;
    }
}

void Foam::cutCellFvMesh::edgesToSide
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

void Foam::cutCellFvMesh::edgesToSide
(
    const edgeList& edges
)
{
    edgesToSide_.setSize(edges.size());
    //labelList edgesToSide(edges.size());
    
    for(int i=0;i<edges.size();i++)
    {
        scalar lvlSetStart = pointDist[edges[i].start()];
        scalar lvlSetEnd = pointDist[edges[i].end()];

        if(lvlSetStart != 0 && lvlSetEnd != 0)
        {
            if(lvlSetStart > 0 && lvlSetEnd > 0)
                edgesToSide_[i] = +1;
            else if(lvlSetStart < 0 && lvlSetEnd < 0)
                edgesToSide_[i] = -1;
            else
                edgesToSide_[i] = 0;
        }
        else if(lvlSetStart+lvlSetEnd == 0)
        {
            edgesToSide_[i] = 0;
        }
        else
        {
            if(lvlSetStart+lvlSetEnd > 0)
                edgesToSide_[i] = +1;
            else
                edgesToSide_[i] = -1;
        }
    }
    //this->edgesToSide_ = edgesToSide;
}

void Foam::cutCellFvMesh::facesToSide
(
)
{
    const faceList& faces = this->faces();
    facesToSide_.setCapacity(faces.size());
    
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
                facesToSide_[i] = 0;
            else if(posExist)
                facesToSide_[i] = +1;
            else if(negExist)
                facesToSide_[i] = -1;
            else
                facesToSide_[i] = 0;
        }
        else
        {
            if(posExist && negExist)
                facesToSide_[i] = 0;
            else if(posExist)
                facesToSide_[i] = +1;
            else if(negExist)
                facesToSide_[i] = -1;
            else
            {
                FatalErrorInFunction
                << "A face cannot have neither postive,"
                << " negative nor null points"
                << exit(FatalError);
            }
        }
    }
}

void Foam::cutCellFvMesh::facesToSide
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
                << exit(FatalError);
            }
        }
    }
    this->facesToSide_ = facesToSide;
}

void Foam::cutCellFvMesh::cellsToSide
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
                << exit(FatalError);
            }
        }
    }
    this->cellsToSide_ = cellsToSide;
}

void Foam::cutCellFvMesh::cellsToSide
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
                << exit(FatalError);
            }
        }
    }
    this->cellsToSide_ = cellsToSide;
}

void Foam::cutCellFvMesh::projectLevelSet()
{
    const pointField& points = this->points();
    pointDist = scalarList(points.size());
    
    for(int i=0;i<points.size();i++)
    {
        pointDist[i] = levelSet(points[i]);
    }
}

void Foam::cutCellFvMesh::newMeshPoints
(
)
{
    Info<<"-----------------------------------------------"<<endl;
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span;
    
    Info<<"Data 1 preparation ";
    t1 = std::chrono::high_resolution_clock::now();
    
    //Info<<"Starting adding Points"<<endl;
    const cellList& meshCells = this->cells();
    const pointField& basisPoints = this->points();
    const faceList& basisFaces = this->faces();
    const edgeList& basisEdges = this->edges();

    nbrOfPrevPoints = basisPoints.size();
    
    newMeshPointsInFunc.setCapacity(basisPoints.size()*2);
    //newMeshPoints_.setSize(basisPoints.size()*2);
    newMeshPointsInFunc.append(basisPoints);
    
    pointsToSide(basisPoints);
    //Info<<"Point to Side done"<<endl;
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Data 2 preparation ";
    t1 = std::chrono::high_resolution_clock::now();
    
    pointToEgde_.setSize(basisPoints.size());
    for(int i=0;i<basisPoints.size();i++)
    {
        pointToEgde_[i] = -1;
    }
    pointToEgde_.setCapacity(basisPoints.size()*2);
    
    edgeToPoint_.setSize(basisEdges.size());
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Data 3 preparation ";
    t1 = std::chrono::high_resolution_clock::now();
    
    pointToFaces_.setCapacity(basisPoints.size()*2);
    pointToFaces_.setSize(basisPoints.size());
    const labelListList& pointFaces = this->pointFaces();
    const labelListList& edgeFaces = this->edgeFaces();
    for(int i=0;i<basisPoints.size();i++)
    {
        if(pointsToSide_[i] == 0)
        {
            pointToFaces_[i] = pointFaces[i];
        }
    }
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t" << time_span.count() << " seconds."<<endl;
    //Info<<"Point to faces"<<endl;

    Info<<"Data 4 preparation ";
    t1 = std::chrono::high_resolution_clock::now();
    
    faceToPoints_.setSize(basisFaces.size());
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
    //Info<<"Point to faces done"<<endl;
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Data 5 preparation ";
    t1 = std::chrono::high_resolution_clock::now();
    
    pointToCells_.setCapacity(basisPoints.size()*2);
    pointToCells_.setSize(basisPoints.size());
    //Info<<"1"<<endl;
    const labelListList& pointCells = this->pointCells();
    //Info<<"2"<<endl;
    const labelListList& edgeCells = this->edgeCells();
    //Info<<"3"<<endl;
    for(int i=0;i<basisPoints.size();i++)
    {
        //Info<<pointsToSide_[i]<<endl;
        if(pointsToSide_[i] == 0)
        {
            pointToCells_[i] = pointCells[i];
        }
    }
    //Info<<"Point to cells "<<basisPoints.size()<<" done"<<endl;
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Data 6 preparation ";
    t1 = std::chrono::high_resolution_clock::now();
    
    cellToPoints_.setSize(meshCells.size());
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
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Point computing ";
    t1 = std::chrono::high_resolution_clock::now();
    
    label pos,neg;
    
    for(int i=0;i<basisEdges.size();i++)
    {
        label startLabel = basisEdges[i].start();
        label endLabel = basisEdges[i].end();        
        pos = 0;
        neg = 0;
        scalar phiStart = pointDist[startLabel];
        scalar phiEnd = pointDist[endLabel];
        
        if(phiStart>0 || phiEnd>0)
            pos = +1;
        if(phiStart<0 || phiEnd<0)
            neg = -1;
        
        if(pos == +1 && neg == -1)
        {
            //Info<<basisPoints[startLabel]<<" -> "<<basisPoints[endLabel]<<endl;
            Info<<phiStart<<" -> "<<phiEnd<<endl;
            vector startToEnd = basisEdges[i].vec(basisPoints);
            //Info<<"startToEnd: "<<basisEdges[i].vec(basisPoints)<<endl;;
            //scalar norm_startToEnd = basisEdges[i].mag(basisPoints);
            //Info<<"norm_startToEnd: "<<norm_startToEnd<<endl;;
            //scalar distPhi = std::abs(phiEnd-phiStart);
            //Info<<"distPhi: "<<distPhi<<endl;;
            //scalar norm_phiStart = std::abs(phiStart);
            //Info<<"norm_phiStart: "<<norm_phiStart<<endl;;
            scalar scalePoint = phiStart / (phiStart - phiEnd);
            //Info<<"scalePoint: "<<scalePoint<<endl;;
            vector newPoint = basisPoints[startLabel] + scalePoint * startToEnd;
            Info<<"--------------------------------------"<<endl;
            Info<<"dist: "<<distToNursOfEdge(basisPoints[startLabel],basisPoints[endLabel],11)<<endl;
            bool found;
            Info<<"Start: "<<basisPoints[startLabel]<<distToNurbs(basisPoints[startLabel],found)<<endl;
            Info<<"End: "<<basisPoints[endLabel]<<distToNurbs(basisPoints[endLabel],found)<<endl;
            label count=0;
            scalar distOfNew = distToNurbs(newPoint,found);
            if(!found)
                FatalErrorInFunction<<"Not found!"<< exit(FatalError);
            Info<<"startPoint: "<<basisPoints[startLabel]<<endl;
            Info<<"startToEnd:"<<startToEnd<<endl;
            Info<<"phiN: "<<distOfNew<<endl;
            Info<<"scalePointN: "<<scalePoint<<endl;
            Info<<"newPoint: "<<newPoint<<endl;
            Info<<"distOfNew: "<<distOfNew<<endl;
            if(std::abs(distOfNew)>10e-10)
            {
                scalar scalePointN;
                scalar scalePointO = scalePoint;
                scalar phiN;
                scalar phiO = distOfNew;
                Info<<"phiO:"<<phiO<<endl;
                Info<<"phiStart:"<<phiStart<<endl;
                if(std::signbit(phiO)!=std::signbit(phiStart))
                {
                    scalar m=(phiO-phiStart)/(scalePointO-0);
                    scalar b=phiO-m*scalePointO;
                    scalePointN = -b/m;
                }
                else if(std::signbit(phiO)!=std::signbit(phiEnd))
                {
                    scalar m=(phiO-phiEnd)/(scalePointO-1);
                    scalar b=phiO-m*scalePointO;
                    scalePointN = -b/m;
                }
                else
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);

                newPoint = basisPoints[startLabel] + scalePointN * startToEnd;
                phiN = distToNurbs(newPoint,found);
                Info<<"phiN: "<<phiN<<endl;
                Info<<"scalePointN: "<<scalePointN<<endl;
                Info<<"newPoint: "<<newPoint<<endl;
                if(!found)
                    FatalErrorInFunction<<"Not found!"<< exit(FatalError);
                
                while(std::abs(phiN)>10e-10)
                {
                    if(count>=90)
                    {
                        Info<<"Over 90"<<endl;
                        Info<<"phiN: "<<phiN<<endl;
                        Info<<"scalePointN: "<<scalePointN<<endl;
                        Info<<"newPoint: "<<newPoint<<endl;
                    }
                    if(count>100)
                    {
                        FatalErrorInFunction<<"Can not find a zero point!"<< exit(FatalError);
                    }
                    Info<<"scalePointO:"<<scalePointO<<endl;
                    Info<<"scalePointN:"<<scalePointN<<endl;
                    Info<<"phiO:"<<phiO<<endl;
                    Info<<"phiN:"<<phiN<<endl;
                    if(phiO!=phiN || scalePointO!=scalePointN)
                    {
                        scalar m=(phiO-phiN)/(scalePointO-scalePointN);
                        Info<<"m:"<<m<<endl;
                        scalar b=phiO-m*scalePointO;
                        Info<<"b:"<<b<<endl;
                        scalar x = -b/m;
                        Info<<"x:"<<x<<endl;
                        phiO = phiN;
                        scalePointO = scalePointN;
                        phiO = phiN;
                        scalePointN = x;
                    }
                    else
                    {
                        if(std::signbit(phiO)!=std::signbit(phiStart))
                        {
                            scalar m=(phiO-phiStart)/(scalePointO-0);
                            scalar b=phiO-m*scalePointO;
                            scalePointN = -b/m;
                        }
                        else if(std::signbit(phiO)!=std::signbit(phiEnd))
                        {
                            scalar m=(phiO-phiEnd)/(scalePointO-1);
                            scalar b=phiO-m*scalePointO;
                            scalePointN = -b/m;
                        }
                        else
                            FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                    }
                    newPoint = basisPoints[startLabel] + scalePointN * startToEnd;
                    phiN = distToNurbs(newPoint,found);
                    if(!found)
                        FatalErrorInFunction<<"Not found!"<< exit(FatalError);
                    count++;
                    Info<<"--dist:"<<phiN<<endl;
                }
                //FatalErrorInFunction<<"Temp Stop!"<< exit(FatalError);
            }
            Info<<"End"<<endl;
            
            if(startLabel==14315 || endLabel==14315)
            {
                Info<<"startLabel:"<<startLabel<<endl;
                Info<<"endLabel:"<<endLabel<<endl;
                Info<<phiStart<<" -> "<<phiEnd<<endl;
                Info<<"Added: "<<newPoint<<endl<<endl;
            }

            newMeshPointsInFunc.append(newPoint);
            
            pointsToSide_.append(0);
            
            pointToEgde_.append(i);
            
            edgeToPoint_[i] = newMeshPointsInFunc.size()-1;
                        
            DynamicList<label> insertedgeFaces;
            insertedgeFaces.setSize(edgeFaces[i].size());
            for(int j=0;j<edgeFaces[i].size();j++)
                insertedgeFaces[j] = edgeFaces[i][j];
            pointToFaces_.append(insertedgeFaces);
            
            for(int k=0;k<edgeFaces[i].size();k++)
            {
                label faceLabel = edgeFaces[i][k];
                if(faceToPoints_[faceLabel].size() == 0)
                {
                    faceToPoints_[faceLabel] = labelList(0);
                    faceToPoints_[faceLabel].append(newMeshPointsInFunc.size()-1);
                }
                else
                {
                    faceToPoints_[faceLabel].append(newMeshPointsInFunc.size()-1);
                }
            }
            
            DynamicList<label> insertedgeCells;        
            insertedgeCells.setSize(edgeCells[i].size());
            for(int j=0;j<edgeCells[i].size();j++)
                insertedgeCells[j] = edgeCells[i][j];
            pointToCells_.append(insertedgeCells);
            
            for(int k=0;k<edgeCells[i].size();k++)
            {
                label cellLabel = edgeCells[i][k];
                if(cellToPoints_[cellLabel].size() == 0)
                {
                    cellToPoints_[cellLabel] = labelList(0);
                    cellToPoints_[cellLabel].append(newMeshPointsInFunc.size()-1);
                }
                else
                {
                    cellToPoints_[cellLabel].append(newMeshPointsInFunc.size()-1);
                }
            }
        }
        else
        {
            edgeToPoint_[i] = -1;
        }
    }
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"End information ";
    t1 = std::chrono::high_resolution_clock::now();
    
    newMeshPoints_ = pointField(newMeshPointsInFunc.size());
    for(int i=0;i<newMeshPoints_.size();i++)
        newMeshPoints_[i] = newMeshPointsInFunc[i];
    newMeshPointsInFunc.setCapacity(0);
    
    pointToEgde_.setCapacity(pointToEgde_.size());
    pointToFaces_.setCapacity(pointToFaces_.size());
    pointToCells_.setCapacity(pointToFaces_.size());
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t\t" << time_span.count() << " seconds."<<endl;
    Info<<"-----------------------------------------------";

}

void Foam::cutCellFvMesh::printAddedPoints
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

scalar Foam::cutCellFvMesh::distToNurbs
(
    point pnt,
    bool& foundFlag
)
{
    scalar dist;
    std::unique_ptr<labelList> firstOrderNearNurbs = MainTree->nearNurbsCurves(pnt);
    if(firstOrderNearNurbs->size()==0)
    {
        Info<<"No firstOderNurbs found"<<endl;
        foundFlag = false;
        return -1;
    }
    DynamicList<scalar> distToNurbsSurface;
    DynamicList<scalar> paraToNurbsSurface;
    bool allOutSideNurbsBox = true;
    for(int k=0;k<firstOrderNearNurbs->size();k++)
    {
        label thisNurbs = (*firstOrderNearNurbs)[k];
        scalar thisNodePara = NurbsTrees[thisNurbs]->closestParaOnNurbsToPoint(pnt);
        //Info<<"\tIndex of nurbs:"<<thisNurbs<<" with para: "<<thisNodePara<<endl;
        if(thisNodePara < this->Curves[thisNurbs]->min_U())
        {
            dist = std::numeric_limits<scalar>::max();
            continue;
        }
        allOutSideNurbsBox = false;
        paraToNurbsSurface.append(thisNodePara);
        distToNurbsSurface.append(this->Curves[thisNurbs]->distanceToNurbsSurface(thisNodePara,pnt));
    }
    if(allOutSideNurbsBox)
    {
        foundFlag = false;
        return -1;
    }
    scalar minDistToNurbsSurface = std::numeric_limits<scalar>::max();
    for(int k=0;k<distToNurbsSurface.size();k++)
    {
        if(distToNurbsSurface[k] < minDistToNurbsSurface)
            minDistToNurbsSurface = distToNurbsSurface[k];
    }
    foundFlag = true;
    dist = minDistToNurbsSurface;
    return dist;
}

scalar Foam::cutCellFvMesh::nearestNurbsIndexPara
(
    point pnt,
    bool& foundFlag,
    label& nurbsInd,
    scalar& nurbsPara
)
{
    scalar dist;
    std::unique_ptr<labelList> firstOrderNearNurbs = MainTree->nearNurbsCurves(pnt);
    if(firstOrderNearNurbs->size()==0)
    {
        Info<<"No firstOderNurbs found"<<endl;
        foundFlag = false;
        return -1;
    }
    DynamicList<label> indNurbs;
    DynamicList<scalar> distToNurbsSurface;
    DynamicList<scalar> paraToNurbsSurface;
    bool allOutSideNurbsBox = true;
    for(int k=0;k<firstOrderNearNurbs->size();k++)
    {
        label thisNurbs = (*firstOrderNearNurbs)[k];
        scalar thisNodePara = NurbsTrees[thisNurbs]->closestParaOnNurbsToPoint(pnt);
        //Info<<"\tIndex of nurbs:"<<thisNurbs<<" with para: "<<thisNodePara<<endl;
        if(thisNodePara < this->Curves[thisNurbs]->min_U())
        {
            dist = std::numeric_limits<scalar>::max();
            continue;
        }
        allOutSideNurbsBox = false;
        indNurbs.append(thisNurbs);
        paraToNurbsSurface.append(thisNodePara);
        distToNurbsSurface.append(this->Curves[thisNurbs]->distanceToNurbsSurface(thisNodePara,pnt));
    }
    if(allOutSideNurbsBox)
    {
        foundFlag = false;
        return -1;
    }
    scalar minDistToNurbsSurface = std::numeric_limits<scalar>::max();
    label minDistNurbsInd = 0;
    scalar minDistPara = 0;
    Info<<"distToNurbsSurface:"<<distToNurbsSurface<<endl;
    Info<<"indNurbs:"<<indNurbs<<endl;
    Info<<"paraToNurbsSurface:"<<paraToNurbsSurface<<endl;
    for(int k=0;k<distToNurbsSurface.size();k++)
    {
        if(distToNurbsSurface[k] < minDistToNurbsSurface)
        {
            minDistToNurbsSurface = distToNurbsSurface[k];
            minDistNurbsInd = indNurbs[k];
            minDistPara = paraToNurbsSurface[k];
        }
    }
    foundFlag = true;
    nurbsInd = minDistNurbsInd;
    nurbsPara = minDistPara+minDistNurbsInd;
    return minDistToNurbsSurface;
}

List<scalar> Foam::cutCellFvMesh::distToNursOfEdge
(
    point startPoint,
    point endPoint,
    label nbrOfPoints
)
{
    List<scalar> distOfEdge(nbrOfPoints);
    vector connec = endPoint - startPoint;
    Info<<endl<<"connec:"<<connec<<endl;
    scalar stepSize = 1.0/static_cast<scalar>(nbrOfPoints-1);
    Info<<"stepSize:"<<stepSize<<endl;
    for(int i=0;i<nbrOfPoints;i++)
    {
        vector pnt = startPoint + connec*stepSize*i;
        Info<<"i:"<<i<<pnt<<endl;
        bool found;
        distOfEdge[i] = distToNurbs(pnt,found);
        if(!found)
            FatalErrorInFunction<<"Not found vector to Nurbs. Can not happen."<< exit(FatalError);
    }
    return distOfEdge;
}

vector Foam::cutCellFvMesh::vectorToNurbs
(
    point pnt,
    bool& foundFlag
)
{
    scalar dist;
    std::unique_ptr<labelList> firstOrderNearNurbs = MainTree->nearNurbsCurves(pnt);
    if(firstOrderNearNurbs->size()==0)
    {
        Info<<"No firstOderNurbs found"<<endl;
        foundFlag = false;
        return vector();
    }
    DynamicList<scalar> distToNurbsSurface;
    DynamicList<scalar> paraToNurbsSurface;
    DynamicList<label> nurbsInd;
    bool allOutSideNurbsBox = true;
    for(int k=0;k<firstOrderNearNurbs->size();k++)
    {
        label thisNurbs = (*firstOrderNearNurbs)[k];
        scalar thisNodePara = NurbsTrees[thisNurbs]->closestParaOnNurbsToPoint(pnt);
        //Info<<"\tIndex of nurbs:"<<thisNurbs<<" with para: "<<thisNodePara<<endl;
        if(thisNodePara < this->Curves[thisNurbs]->min_U())
        {
            dist = std::numeric_limits<scalar>::max();
            continue;
        }
        allOutSideNurbsBox = false;
        paraToNurbsSurface.append(thisNodePara);
        distToNurbsSurface.append(this->Curves[thisNurbs]->distanceToNurbsSurface(thisNodePara,pnt));
        nurbsInd.append(thisNurbs);
    }
    if(allOutSideNurbsBox)
    {
        foundFlag = false;
        return vector();
    }
    label minNurbs = -1;
    scalar minParaToNurbsSurface;
    scalar minDistToNurbsSurface = std::numeric_limits<scalar>::max();
    for(int k=0;k<distToNurbsSurface.size();k++)
    {
        if(distToNurbsSurface[k] < minDistToNurbsSurface)
        {
            minDistToNurbsSurface = distToNurbsSurface[k];
            minParaToNurbsSurface = paraToNurbsSurface[k];
            minNurbs = nurbsInd[k];
        }
    }
    foundFlag = true;
    dist = minDistToNurbsSurface;
    scalar para = minParaToNurbsSurface;
    vector nurbsPoint = this->Curves[minNurbs]->Curve_Derivative(0,para);
    vector pointToNurbsVector = nurbsPoint-pnt;
    return pointToNurbsVector;
}

List<vector> Foam::cutCellFvMesh::vectorsToNurbsOfEdge
(
    point startPoint,
    point endPoint,
    label nbrOfVectors
)
{
    List<vector> vectorsToNurbs(nbrOfVectors);
    vector connec = endPoint - startPoint;
    Info<<endl<<"connec:"<<connec<<endl;
    scalar stepSize = 1.0/static_cast<scalar>(nbrOfVectors-1);
    Info<<"stepSize:"<<stepSize<<endl;
    for(int i=0;i<nbrOfVectors;i++)
    {
        vector pnt = startPoint + connec*stepSize*i;
        Info<<"i:"<<i<<pnt<<endl;
        bool found;
        vectorsToNurbs[i] = vectorToNurbs(pnt,found);
        if(!found)
            FatalErrorInFunction<<"Not found vector to Nurbs. Can not happen."<< exit(FatalError);
    }
    return vectorsToNurbs;
}

label Foam::cutCellFvMesh::sideToNurbs(point pnt,bool& foundFlag)
{
    scalar dist = distToNurbs(pnt,foundFlag);
    
    label side;
    if(dist > 0)
        side = 1;
    else if(dist < 0)
        side = -1;
    else
        side = 0;
    return side;
}

scalar norm2(vector pnt)
{
    return std::sqrt(pnt.x()*pnt.x()+pnt.y()*pnt.y()+pnt.z()*pnt.z());
}

bool facesShareEdge
(
    const face& faceA,
    const face& faceB
)
{
    std::unordered_set<label> faceA_map;
    std::unordered_set<label> faceB_map;
    for(int i=0;i<faceA.size();i++)
        faceA_map.insert(faceA[i]);
    for(int i=0;i<faceB.size();i++)
        faceB_map.insert(faceB[i]);
    
    DynamicList<label> edge;
    for(const label& a :faceA_map)
        if(faceB_map.count(a)!=0)
            edge.append(a);
    
    if(edge.size()==0)
        FatalErrorInFunction<<"Can not happen. Faces must be connected!"<< exit(FatalError);
    if(edge.size()>=3)
    {
        Info<<"faceA:"<<faceA<<endl;
        Info<<"faceB:"<<faceB<<endl;
        FatalErrorInFunction<<"Can not happen. Faces are identical!"<< exit(FatalError);
    }
    if(edge.size()==1)
        return false;
    if(edge.size()==2)
    {
        List<label> localIndA(2);
        List<label> localIndB(2);
        for(int i=0;i<edge.size();i++)
        {
            localIndA[i] = faceA.which(edge[i]);
            localIndB[i] = faceB.which(edge[i]);
        }
        if(((localIndA[0]+1)%faceA.size()==localIndA[1] || (localIndA[1]+1)%faceA.size()==localIndA[0]) &&
            ((localIndB[0]+1)%faceB.size()==localIndB[1] || (localIndB[1]+1)%faceB.size()==localIndB[0]))
            return true;
        else
        {
            Info<<"faceA:"<<faceA<<endl;
            Info<<"faceB:"<<faceB<<endl;
            FatalErrorInFunction<<"Can not happen. Faces are identical!"<< exit(FatalError);
        }
    }
    return false;
}

void computeClosedFaceFront
(
    label centerPointInd,
    List<DynamicList<face>>& facesInCellsIn,
    List<DynamicList<std::unordered_set<label>>>& facesMapInCellsIn,
    DynamicList<DynamicList<face>>& closedFaceFrontOut,
    DynamicList<DynamicList<std::unordered_set<label>>>& closedFaceFrontMapOut
)
{
    if(facesInCellsIn.size()!=facesMapInCellsIn.size())
        FatalErrorInFunction<<"Unequal range of parameters!"<< exit(FatalError);
    
    //Remove faces not touching the centerPointInd
    List<DynamicList<face>> facesOfCells(facesInCellsIn.size());
    List<DynamicList<std::unordered_set<label>>> facesOfCellsMap(facesInCellsIn.size());
    for(int i=0;i<facesInCellsIn.size();i++)
    {
        for(int j=0;j<facesInCellsIn[i].size();j++)
        {
            if(facesMapInCellsIn[i][j].count(centerPointInd)!=0)
            {
                facesOfCells[i].append(facesInCellsIn[i][j]);
                facesOfCellsMap[i].append(facesMapInCellsIn[i][j]);
            }
        }
    }
    
    
    DynamicList<DynamicList<std::pair<label,label>>> uniqueCycles;
    DynamicList<std::unordered_map<label,label>> uniqueCyclesMap;    
    //Combine facesInCells to faceFronts
    enum fcState {OPEN,CLOSED,FAIL};    
    for(int i=0;i<facesOfCells.size();i++)
    {
        for(int j=0;j<facesOfCells[i].size();j++)
        {
            std::pair<label,label> startFace(i,j);
            DynamicList<std::unordered_set<label>> usedCells;
            usedCells.setSize(1);
            usedCells[0].insert(i);
            DynamicList<DynamicList<std::pair<label,label>>> faceCycle;
            faceCycle.setSize(1);
            faceCycle[0].append(startFace);
            DynamicList<fcState> state;
            state.setSize(1);
            state[0]=OPEN;
            bool allDone = false;
            
            while(!allDone)
            {
                label len = faceCycle.size();
                for(int k=0;k<len;k++)
                //Iterate over all cycles
                {
                    if(state[k] != OPEN)
                        continue;
                    
                    std::pair<label,label> currFace = faceCycle[k].last();
                    
                    // Test if cycle is closed
                    label equalPnts = 0;
                    for(int n=0;n<facesOfCells[startFace.first][startFace.second].size();n++)
                    {
                        if(facesOfCellsMap[currFace.first][currFace.second].count(facesOfCells[startFace.first][startFace.second][n])!=0)
                            equalPnts++;
                    }
                    if(faceCycle[k].size()>2 && equalPnts==2)
                    {
                        state[k] = CLOSED;
                        continue;
                    }
                    else if((faceCycle[k].size()==1 && equalPnts!=4) || (faceCycle[k].size()==2 && equalPnts!=2))
                        FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                    
                    DynamicList<std::pair<label,label>> addedFaces;
                    for(int l=0;l<facesOfCells.size();l++)
                    {
                        if(usedCells[k].count(l)!=0)
                            continue;
                        
                        for(int m=0;m<facesOfCells[l].size();m++)
                        {
                            label equalPnts = 0;
                            for(int n=0;n<facesOfCells[l][m].size();n++)
                            {
                                if(facesOfCellsMap[currFace.first][currFace.second].count(facesOfCells[l][m][n])!=0)
                                    equalPnts++;
                            }
                            if(equalPnts==2)
                            {
                                addedFaces.append(std::pair<label,label>(l,m));
                            }
                            else if(equalPnts!=1)
                                FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                        }
                    }
                    if(addedFaces.size()==0)
                    {
                        state[k] = FAIL;
                    }
                    else if(addedFaces.size()==1)
                    {
                        faceCycle[k].append(addedFaces[0]);
                    }
                    else
                    {
                        DynamicList<std::pair<label,label>> cpCycle = faceCycle[k];
                        fcState cpState = state[k];
                        faceCycle[k].append(addedFaces[0]);
                        for(int o=1;o<addedFaces.size();o++)
                        {
                            faceCycle.append(cpCycle);
                            faceCycle.last().append(addedFaces[o]);
                            state.append(cpState);
                        }
                    }
                }
                bool noOpen=true;
                for(int k=0;k<state.size();k++)
                {
                    if(state[k]==OPEN)
                        noOpen=false;
                }
            }
            for(int k=0;k<faceCycle.size();k++)
            //iterate across all found closed face cycles
            {
                if(state[k]==FAIL)
                    continue;
                if(state[k]==OPEN)
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);

                DynamicList<std::pair<label,label>> oneCylce = faceCycle[k];
                bool noMatch = true;
                for(int l=0;l<uniqueCyclesMap.size();l++)
                //test one Face cycle against one map
                {
                    bool allFaceMatch = true;
                    std::unordered_map<label,label> oneCycleMap = uniqueCyclesMap[l];
                    for(int m=0;m<oneCylce.size();m++)
                    //test one face against one cycleMap    
                    {
                        std::pair<label,label> oneFace = oneCylce[m];
                        auto item = oneCycleMap.find(oneFace.first);
                        if(!((item != oneCycleMap.end()) && (item->second==oneFace.second)))
                        //face does not exist in map
                        {                              
                            allFaceMatch = false;
                        }
                    }
                    if(allFaceMatch)
                        noMatch = false;
                }
                if(noMatch)
                {
                    uniqueCycles.append(faceCycle[k]);
                    std::unordered_map<label,label> cycleMap;
                    for(int l=0;l<faceCycle[k].size();l++)
                    {
                        cycleMap.insert(faceCycle[k][l]);
                    }
                    uniqueCyclesMap.append(cycleMap);
                }
            }
        }
    }
    
    if(uniqueCycles.size()==0)
        FatalErrorInFunction<<"There must be at least one face front!"<< exit(FatalError);
    for(int i=0;i<uniqueCycles.size();i++)
    {
        if(uniqueCycles[i].size()<3 || uniqueCycles[i].size()>8)
            FatalErrorInFunction<<"There must be between 3 and 8 faces in a front!"<< exit(FatalError);
    }
    
    for(int i=0;i<uniqueCycles.size();i++)
    {
        closedFaceFrontOut.append(DynamicList<face>());
        closedFaceFrontMapOut.append(DynamicList<std::unordered_set<label>>());
        for(int j=0;j<uniqueCycles[i].size();j++)
        {
            std::pair<label,label> face = uniqueCycles[i][j];
            closedFaceFrontOut.last().append(facesOfCells[face.first][face.second]);
            closedFaceFrontMapOut.last().append(facesOfCellsMap[face.first][face.second]);
        }
    }
    
    /*
    Info<<"centerPointInd:"<<centerPointInd<<endl;
    Info<<"facesInCellsIn:"<<facesInCellsIn<<endl;
        
    labelList offset(facesInCellsIn.size(),0);
    label counter = 0;
    for(int i=0;i<offset.size();i++)
    {
        offset[i] = counter;
        counter += facesInCellsIn[i].size();
    }    
    DynamicList<DynamicList<std::pair<label,label>>> faceFrontIndList;
    DynamicList<std::unordered_set<label>> faceFrontMap;
    
    Info<<"In Function----------------------------------"<<endl;
    Info<<"facesInCellsIn"<<facesInCellsIn<<endl;
    Info<<"centerPointInd"<<centerPointInd<<endl;
    
    for(int i=0;i<facesInCellsIn.size();i++)
    {
        for(int j=0;j<facesInCellsIn[i].size();j++)
        {
            if(facesMapInCellsIn[i][j].count(centerPointInd)==0)
                continue;
            
            std::pair<label,label> startFace(i,j);
            Info<<"startFace:("<<startFace.first<<"|"<<startFace.second<<")"<<endl;
            DynamicList<DynamicList<std::pair<label,label>>> faceConnectionsOfInitalFace;

            DynamicList<std::pair<label,label>> neighborFaces;
            for(int k=0;k<facesInCellsIn.size();k++)
            {
                if(k==i)
                    continue;
                for(int l=0;l<facesInCellsIn[k].size();l++)
                {
                    //Info<<"k:"<<k<<"  l:"<<l<<"   startFace.first:"<<startFace.first<<"   startFace.second:"<<startFace.second<<endl;
                    if(facesShareEdge(facesInCellsIn[k][l],facesInCellsIn[startFace.first][startFace.second]))
                    {
                        neighborFaces.append(std::pair<label,label>(k,l));
                    }
                }
            }
            
            Info<<"cell:"<<i<<" face:"<<j<<" <->";
            for(int l=0;l<neighborFaces.size();l++)
            {
                Info<<"  cell:"<<neighborFaces[l].first<<" face:"<<neighborFaces[l].second;
            }
            Info<<endl;
            
            DynamicList<label> neighborCells;   //list of all neighboring cells of this face
            std::unordered_map<label,label> neighborCellsMap; // map cell to neighborCells index
            for(int k=0;k<neighborFaces.size();k++)
            {
                if(neighborCellsMap.count(neighborFaces[k].first)==0)
                {
                    neighborCellsMap.insert(std::pair<label,label>(neighborFaces[k].first,neighborCells.size()));
                    neighborCells.append(neighborFaces[k].first);
                }
            }
            Info<<"neighborCells:"<<neighborCells<<endl;
            Info<<"neighborCellsMap:";
            for(auto& key:neighborCellsMap)
            {
                Info<<" ("<<key.first<<"|"<<key.second<<")";
            }
            Info<<endl;
            DynamicList<DynamicList<std::pair<label,label>>> cellwiseNeighborFaces;
            cellwiseNeighborFaces.setSize(neighborCells.size());
            for(int k=0;k<neighborFaces.size();k++)
            {
                cellwiseNeighborFaces[neighborCellsMap[neighborFaces[k].first]].append(neighborFaces[k]);
            }
            Info<<"cellwiseNeighborFaces:"<<endl;
            for(int k=0;k<cellwiseNeighborFaces.size();k++)
            {
                Info<<"\t\tk:"<<k<<" ";
                for(int l=0;l<cellwiseNeighborFaces[k].size();l++)
                {
                    Info<<"("<<cellwiseNeighborFaces[k][l].first<<"|"<<cellwiseNeighborFaces[k][l].second<<")";
                }
                Info<<endl;
            }
            if(facesInCellsIn.size()==1)
            {
                DynamicList<std::pair<label,label>> temp({startFace});                                        
                faceConnectionsOfInitalFace.append(temp);
            }
            else if(facesInCellsIn.size()==4 || facesInCellsIn.size()==8)
            {
                for(int k=0;k<cellwiseNeighborFaces.size();k++)
                {
                    for(int l=0;l<cellwiseNeighborFaces[k].size();l++)
                    {
                        std::pair<label,label> faceA = cellwiseNeighborFaces[k][l];
                        Info<<"faceA: ("<<faceA.first<<"|"<<faceA.second<<")"<<endl;
                        if(facesMapInCellsIn[faceA.first][faceA.second].count(centerPointInd)==0)
                            continue;
                        
                        for(int m=0;m<cellwiseNeighborFaces.size();m++)
                        {
                            if(k==m)
                                //Avoid building pairs of faces in same cell
                                continue;
                            
                            for(int n=0;n<cellwiseNeighborFaces[m].size();n++)
                            {
                                std::pair<label,label> faceB = cellwiseNeighborFaces[m][n];
                                Info<<"faceB: ("<<faceB.first<<"|"<<faceB.second<<")"<<endl;
                                if(facesMapInCellsIn[faceB.first][faceB.second].count(centerPointInd)==0)
                                    continue;
                                
                                for(int o=0;o<facesInCellsIn.size();o++)
                                {
                                    if(o==faceA.first || o==faceB.first || o==startFace.first)
                                        //Avoid searching closing face at cell of faceA, faceB or startFace
                                        continue;
                                    
                                    for(int p=0;p<facesInCellsIn[o].size();p++)
                                    {
                                        std::pair<label,label> posFourthFace(o,p);
                                        if(facesShareEdge(facesInCellsIn[posFourthFace.first][posFourthFace.second],facesInCellsIn[faceA.first][faceA.second]) &&
                                        facesShareEdge(facesInCellsIn[posFourthFace.first][posFourthFace.second],facesInCellsIn[faceB.first][faceB.second]))
                                        {
                                            DynamicList<std::pair<label,label>> temp({startFace,faceA,faceB,posFourthFace});                                        
                                            faceConnectionsOfInitalFace.append(temp);
                                            Info<<"Added:";
                                            for(int x=0;x<temp.size();x++)
                                            {
                                                Info<<" ("<<temp[x].first<<"|"<<temp[x].second<<")"<<endl;
                                            }
                                            Info<<endl;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(facesInCellsIn.size()==4)
                {
                    for(int k=0;k<cellwiseNeighborFaces.size();k++)
                    {
                        for(int l=0;l<cellwiseNeighborFaces[k].size();l++)
                        {
                            std::pair<label,label> faceA = cellwiseNeighborFaces[k][l];
                            DynamicList<std::pair<label,label>> temp({startFace,faceA});                                        
                            faceConnectionsOfInitalFace.append(temp);
                        }
                    }
                }
            }
            else
                FatalErrorInFunction<<"Can not happen. Cell count must be 1, 4 or 8!"<< exit(FatalError);
            for(int k=0;k<faceConnectionsOfInitalFace.size();k++)
            {
                bool faceConnectionIsOld = false;
                for(int l=0;l<faceFrontMap.size();l++)
                {
                    bool faceFrontEquals = true;
                    for(int m=0;m<faceConnectionsOfInitalFace[k].size();m++)
                    {
                        std::pair<label,label> face = faceConnectionsOfInitalFace[k][m];
                        label faceNbr = offset[face.first]+face.second;
                        if(faceFrontMap[l].count(faceNbr)==0)
                            faceFrontEquals = false;
                    }
                    if(faceFrontMap[l].size()!=faceConnectionsOfInitalFace[k].size())
                        faceFrontEquals = false;
                    if(faceFrontEquals)
                        faceConnectionIsOld = true;
                }
                if(!faceConnectionIsOld)
                {
                    faceFrontIndList.append(faceConnectionsOfInitalFace[k]);
                    faceFrontMap.append(std::unordered_set<label>());
                    for(int m=0;m<faceConnectionsOfInitalFace[k].size();m++)
                    {
                        std::pair<label,label> face = faceConnectionsOfInitalFace[k][m];
                        label faceNbr = offset[face.first]+face.second;
                        faceFrontMap.last().insert(faceNbr);
                    }
                }
            }
        }
    }
    Info<<"faceFrontIndList.size:"<<faceFrontIndList.size()<<endl;
    for(int i=0;i<faceFrontIndList.size();i++)
    {
        closedFaceFrontOut.append(DynamicList<face>());
        closedFaceFrontMapOut.append(DynamicList<std::unordered_set<label>>());
        for(int j=0;j<faceFrontIndList[i].size();j++)
        {
            std::pair<label,label> face = faceFrontIndList[i][j];
            closedFaceFrontOut.last().append(facesInCellsIn[face.first][face.second]);
            closedFaceFrontMapOut.last().append(facesMapInCellsIn[face.first][face.second]);
        }
    }
    */
}

List<bool> pointInFaceFront
(
    DynamicList<DynamicList<face>>& closedFaceFrontIn,
    DynamicList<DynamicList<std::unordered_set<label>>>& closedFaceFrontMapIn,
    label pointInd
)
{
    if(closedFaceFrontIn.size()!=closedFaceFrontMapIn.size())
        FatalErrorInFunction<<"Unequal range of parameters!"<< exit(FatalError);

    List<bool> testPointInside(closedFaceFrontIn.size());
    for(int i=0;i<closedFaceFrontMapIn.size();i++)
    {
        bool pointInside = false;
        for(int j=0;j<closedFaceFrontMapIn[i].size();j++)
        {
            if(closedFaceFrontMapIn[i][j].count(pointInd)!=0)
                pointInside = true;
        }
        testPointInside[i] = pointInside;
    }    
    return testPointInside;
}

scalar computeFaceFrontAngle
(
    DynamicList<DynamicList<face>>& oneZeroPointClosedFaces,
    DynamicList<DynamicList<std::unordered_set<label>>>& oneZeroPointClosedFaceMap,
    edge innerEdge,
    label centralPoint,
    pointField& points
)
{
    label outerPoint = innerEdge.otherVertex(centralPoint);
    DynamicList<DynamicList<face>> edgeMatchingZeroPointClosedFaces;
    if(oneZeroPointClosedFaceMap.size()!=2 || oneZeroPointClosedFaces.size()!=2)
        FatalErrorInFunction<<"Four edge face with a point that has other than two face fronts! Can not happen! May be wrong limitation but must be at least equal or more than two!"<< exit(FatalError);
    for(int i=0;i<oneZeroPointClosedFaceMap.size();i++)
    {
        bool faceFrontInside = true;
        label centralPointInsideCnt = 0;
        label outerPointInsideCnt = 0;
        for(int j=0;j<oneZeroPointClosedFaceMap[i].size();j++)
        {
            if(oneZeroPointClosedFaceMap[i][j].count(centralPoint)==1)
                centralPointInsideCnt++;
            if(oneZeroPointClosedFaceMap[i][j].count(outerPoint)==1)
                outerPointInsideCnt++;
        }
        if(centralPointInsideCnt!=4)
            FatalErrorInFunction<<"Face front must have central Point in all four faces!"<< exit(FatalError);
        if(outerPointInsideCnt==2)
            edgeMatchingZeroPointClosedFaces.append(oneZeroPointClosedFaces[i]);
        else if(outerPointInsideCnt==0)
        {}
        else
            FatalErrorInFunction<<"Outer Point must be either twice or not in face front!"<< exit(FatalError);
    }
    if(edgeMatchingZeroPointClosedFaces.size()!=1)
        FatalErrorInFunction<<"There must exist exactly one face front. May be wrong limitation but must be at least equal or more than two!"<< exit(FatalError);
    DynamicList<face> twoOuterFaces;
    for(int i=0;i<edgeMatchingZeroPointClosedFaces[0].size();i++)
    {
        if(edgeMatchingZeroPointClosedFaces[0][i].which(outerPoint)==-1)
            twoOuterFaces.append(edgeMatchingZeroPointClosedFaces[0][i]);
    }
    if(twoOuterFaces.size()!=2)
        FatalErrorInFunction<<"There must exist two faces that do no contain the outer Point!"<< exit(FatalError);
    
    DynamicList<label> sameOuterPoint;
    for(int i=0;i<twoOuterFaces[0].size();i++)
    {
        for(int j=0;j<twoOuterFaces[1].size();i++)
        {
            if((twoOuterFaces[0][i]==twoOuterFaces[1][j]) && twoOuterFaces[0][i]!=centralPoint)
                sameOuterPoint.append(twoOuterFaces[0][i]);
        }
    }
    if(sameOuterPoint.size()!=1)
        FatalErrorInFunction<<"Other than one Outer Point! Can not happen!"<< exit(FatalError);
    
    label connectVertextInd = centralPoint;
    if(connectVertextInd==-1)
        FatalErrorInFunction<<"No connect point! Can not happen!"<< exit(FatalError);
    label otherVertexInnerInd = outerPoint;
    label otherVertexOuterInd = sameOuterPoint[0];
    
    point connectVertex = points[connectVertextInd];
    point otherVertexInner = points[otherVertexInnerInd];
    point otherVertexOuter = points[otherVertexOuterInd];
    
    vector innerVector = connectVertex-otherVertexInner;
    vector outerVector = otherVertexOuter-connectVertex;
    
    scalar normInnerVector = norm2(outerVector);
    scalar normOuterVector = norm2(outerVector);
    scalar scalProd = innerVector & outerVector;
    
    scalar pseudoAngle = scalProd/(normInnerVector*normOuterVector);
    if(pseudoAngle>1 || pseudoAngle<-1)
        FatalErrorInFunction<<"Pseudo angle must be inside [-1,1]! Can not happen!"<< exit(FatalError);
    return pseudoAngle;
}

void Foam::cutCellFvMesh::newMeshEdges
(
)
{
    //Info<<"Starting adding Edges"<<endl;
    const cellList& meshCells = this->cells();
    const pointField& basisPoints = this->points();
    const faceList& basisFaces = this->faces();
    const edgeList& basisEdges = this->edges();
    const labelList& cellOwner = this->owner();
    const labelList& cellNeighbor = this->neighbour();
    const labelListList& faceToEdge = this->faceEdges();
    const labelListList& cellToEdge = this->cellEdges();
    const labelListList& pointToEdge = this->pointEdges();
    
    nbrOfPrevEdges = basisEdges.size();
    
    newMeshEdges_.setCapacity(basisEdges.size()*2);
    newMeshEdges_.append(basisEdges);

    edgesToSide(newMeshEdges_);
    
    edgesToSide_.setCapacity(basisEdges.size()*2);
    
    /*
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        Info<<"Edge:"<<i<<" Side:"<<edgesToSide_[i]<<endl;
    }
    */
    
    //Info<<"Put edges to side"<<endl;
    List<bool> problematicFace(basisFaces.size(),false);
    List<label> problematicFacePoints(basisFaces.size(),-1);
    List<label> problematicFaceNewPoints(basisFaces.size(),-1);
    
    edgeToFaces_.setSize(basisEdges.size());
    const labelListList& edgeFaces = this->edgeFaces();
    for(int i=0;i<basisEdges.size();i++)
    {
        label startPoint = basisEdges[i].start();
        label endPoint = basisEdges[i].end();
        
        if(pointsToSide_[startPoint] == 0 && pointsToSide_[endPoint] == 0)
        {
            edgeToFaces_[i] = edgeFaces[i];
            for(int j=0;j<edgeFaces[i].size();j++)
            {
                problematicFace[edgeFaces[i][j]] = true;
            }
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
        if(thisFacePoints.size() > 4 || basisFaces[i].size() > 4)
        {
            FatalErrorInFunction
            << "Five vertice face not implemented!"
            << exit(FatalError);
        }
        else if(thisFacePoints.size() == 4)
        {
            bool allPointsOld = true;
            label newPointsNum = 0;
            for(int k=0;k<thisFacePoints.size();k++)
            {
                if(thisFacePoints[k] >= nbrOfPrevPoints)
                {
                    allPointsOld = false;
                    newPointsNum++;
                }
            }
            if(allPointsOld)
            {
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
                    << exit(FatalError);
                }    
            }
            else if(newPointsNum == 1)
            {
                FatalErrorInFunction
                << "Face with four cut point but three of them are old ones and one is new! This can not happen "<<exit(FatalError);
            }
            else if(newPointsNum == 2)
            {
                FatalErrorInFunction
                << "Face with four cut point but two of them are old ones and two are new! This can not happen "<<exit(FatalError);
            }
            else if(newPointsNum == 3)
            {
                FatalErrorInFunction
                << "Face with four cut point but one is an old ones and three are new! This can not happen "<<exit(FatalError);
            }
            else if(newPointsNum == 4)
            {
                labelList edgIndx(4);
                edgIndx[0] = pointToEgde_[thisFacePoints[0]];
                edgIndx[1] = pointToEgde_[thisFacePoints[1]];
                edgIndx[2] = pointToEgde_[thisFacePoints[2]];
                edgIndx[3] = pointToEgde_[thisFacePoints[3]];
                
                edgeList edgItem(4);
                edgItem[0] = newMeshEdges_[edgIndx[0]];
                edgItem[1] = newMeshEdges_[edgIndx[1]];
                edgItem[2] = newMeshEdges_[edgIndx[2]];
                edgItem[3] = newMeshEdges_[edgIndx[3]];
                
                labelList facePoints(4,-1);
                
                labelList edgeOrder(4,-1);
                edgeOrder[0] = 0;
                for(int a_1=1;a_1<edgeOrder.size();a_1++)
                {
                    edge pivEdg = edgItem[edgeOrder[a_1-1]];
                    for(int a_2=0;a_2<edgItem.size();a_2++)
                    {
                        if(edgeOrder[0]==a_2||edgeOrder[1]==a_2||edgeOrder[2]==a_2||edgeOrder[3]==a_2)
                            continue;
                        label connectPt = pivEdg.commonVertex(edgItem[a_2]);
                        if(connectPt == -1)
                            continue;
                        facePoints[a_1-1] = connectPt;
                        edgeOrder[a_1] = a_2;
                        break;
                    }
                    if(facePoints[a_1-1] == -1)
                    {
                        Info<<endl;
                        Info<<"basisFaces["<<i<<"]: "<<basisFaces[i]<<endl;
                        Info<<"nbrOfPrevPoints:"<<nbrOfPrevPoints<<endl;
                        Info<<"nbrOfPrevEdges:"<<nbrOfPrevEdges<<endl;
                        Info<<"thisFacePoints:"<<thisFacePoints<<endl;
                        Info<<"edgIndx:"<<edgIndx<<endl;
                        Info<<"pivEdg:"<<pivEdg<<endl;
                        Info<<"edgItem:"<<edgItem<<endl;
                        Info<<"a_1:"<<a_1<<endl;
                        Info<<"facePoints: "<<facePoints<<endl;
                        FatalErrorInFunction<<"No connection point found!"<<exit(FatalError);
                    }
                }
                if((facePoints[edgeOrder.size()-1]=edgItem[edgeOrder[0]].commonVertex(edgItem[edgeOrder[3]]))==-1)
                    FatalErrorInFunction<<"No connection point found!"<<exit(FatalError);
                
                point faceCentre = basisFaces[i].centre(basisPoints);
                bool foundFlag;
                label centrePtToSide = sideToNurbs(faceCentre,foundFlag);
                if(!foundFlag)
                    FatalErrorInFunction<<"No Nurbs near point!"<<exit(FatalError);
                
                DynamicList<label> facePlusPoints;
                DynamicList<label> faceMinusPoints;
                for(int a=0;a<facePoints.size();a++)
                {
                    if(pointsToSide_[facePoints[a]]>0)
                        facePlusPoints.append(a);
                    else if(pointsToSide_[facePoints[a]]<0)
                        faceMinusPoints.append(a);
                    else
                    {
                        Info<<endl;
                        Info<<nbrOfPrevPoints<<endl;
                        Info<<thisFacePoints<<endl;
                        Info<<facePoints<<endl;
                        FatalErrorInFunction<<"Zero point can not be here!"<<exit(FatalError);
                    }
                }
                
                if(facePlusPoints.size()!=2||faceMinusPoints.size()!=2)
                    FatalErrorInFunction<<"Wrong number of plus and minus face points!"<<exit(FatalError);
                if((facePlusPoints[0]+1)%4==facePlusPoints[1] || (facePlusPoints[1]+1)%4==facePlusPoints[0])
                    FatalErrorInFunction<<"Plus point assignment must be wrong!"<<exit(FatalError);
                if((faceMinusPoints[0]+1)%4==faceMinusPoints[1] || (faceMinusPoints[1]+1)%4==faceMinusPoints[0])
                    FatalErrorInFunction<<"Minus point assignment must be wrong!"<<exit(FatalError);
                
                labelList pt(4);
                pt[0] = thisFacePoints[edgeOrder[0]];
                pt[1] = thisFacePoints[edgeOrder[1]];
                pt[2] = thisFacePoints[edgeOrder[2]];
                pt[3] = thisFacePoints[edgeOrder[3]];
                
                /*
                pointField pnt(4);
                pnt[0] = newMeshPoints_[pt[0]];
                pnt[1] = newMeshPoints_[pt[1]];
                pnt[2] = newMeshPoints_[pt[2]];
                pnt[3] = newMeshPoints_[pt[3]];

                pointField pnt_(4);
                pnt_[0] = 0.5*(pnt[0]+pnt[1]);
                pnt_[1] = 0.5*(pnt[1]+pnt[2]);
                pnt_[2] = 0.5*(pnt[2]+pnt[3]);
                pnt_[3] = 0.5*(pnt[3]+pnt[0]);
                
                scalarList distPnt(4);
                for(int a=0;a<4;a++)
                {
                    bool foundFlag;
                    distPnt[a] = sideToNurbs(pnt_[a],foundFlag);
                    if(!foundFlag)
                        FatalErrorInFunction<<"No Nurbs near point!"<<exit(FatalError);
                }
                
                scalar edgPnt_0_2 = distPnt[0]*distPnt[2];
                scalar edgPnt_1_3 = distPnt[1]*distPnt[3];
                
                edgeList edgesToAdd(2);

                if(edgPnt_0_2>0 && edgPnt_1_3<=0)
                //the face cut edges is going through pnt 1-2 and 3-0 because the midpoints 1 and 3 multiplied is <= 0
                {
                    edgesToAdd[0] = edge(pt[1],pt[2]);
                    edgesToAdd[1] = edge(pt[3],pt[0]);
                }
                else if(edgPnt_0_2<=0 && edgPnt_1_3>0)
                //the face cut edges is going through pnt 0-1 and 2-3 because the midpoints 0 and 2 multiplied is <= 0                
                {
                    edgesToAdd[0] = edge(pt[0],pt[1]);
                    edgesToAdd[1] = edge(pt[2],pt[3]);
                }
                else if(edgPnt_0_2<=0 && edgPnt_1_3<=0)
                    FatalErrorInFunction<<"Both twin midpoints are near a cut. This can not happen!"<<exit(FatalError);
                else if(edgPnt_0_2>0 && edgPnt_1_3>0)
                {
                    if((pointsToSide_[facePoints[0]]*distPnt[0] <=0) &&
                       (pointsToSide_[facePoints[1]]*distPnt[1] <=0) &&
                       (pointsToSide_[facePoints[2]]*distPnt[2] <=0) &&
                       (pointsToSide_[facePoints[3]]*distPnt[3] <=0))
                    {
                        FatalErrorInFunction<<"There seem to be four new cut edges in the face. This can not happen!"<<exit(FatalError);
                    }
                    else 
                    if((pointsToSide_[facePoints[0]]*distPnt[0] > 0) &&
                       (pointsToSide_[facePoints[1]]*distPnt[1] <=0) &&
                       (pointsToSide_[facePoints[2]]*distPnt[2] > 0) &&
                       (pointsToSide_[facePoints[3]]*distPnt[3] <=0))
                    {
                        edgesToAdd[0] = edge(pt[1],pt[2]);
                        edgesToAdd[1] = edge(pt[3],pt[0]);
                    }
                    else
                    if((pointsToSide_[facePoints[0]]*distPnt[0] <=0) &&
                       (pointsToSide_[facePoints[1]]*distPnt[1] > 0) &&
                       (pointsToSide_[facePoints[2]]*distPnt[2] <=0) &&
                       (pointsToSide_[facePoints[3]]*distPnt[3] > 0))
                    {
                        edgesToAdd[0] = edge(pt[0],pt[1]);
                        edgesToAdd[1] = edge(pt[2],pt[3]);
                    }
                    else
                    if((pointsToSide_[facePoints[0]]*distPnt[0] > 0) &&
                       (pointsToSide_[facePoints[1]]*distPnt[1] > 0) &&
                       (pointsToSide_[facePoints[2]]*distPnt[2] > 0) &&
                       (pointsToSide_[facePoints[3]]*distPnt[3] > 0))
                    {
                        scalarList edgeLen(4);
                        edgeLen[0] = norm2(pnt[0]-pnt[1]);
                        edgeLen[1] = norm2(pnt[1]-pnt[2]);
                        edgeLen[2] = norm2(pnt[2]-pnt[3]);
                        edgeLen[3] = norm2(pnt[3]-pnt[0]);
                        
                        scalarList normedZeroDist(4);
                        for(int a=0;a<4;a++)
                            normedZeroDist[a] = abs(distPnt[a])/edgeLen[a];
                        
                        if((normedZeroDist[0]>normedZeroDist[1] && normedZeroDist[0]>normedZeroDist[3]) &&
                           (normedZeroDist[2]>normedZeroDist[1] && normedZeroDist[2]>normedZeroDist[3]))
                        {
                            edgesToAdd[0] = edge(pt[0],pt[1]);
                            edgesToAdd[1] = edge(pt[2],pt[3]);
                        }
                        else
                        if((normedZeroDist[1]>normedZeroDist[0] && normedZeroDist[1]>normedZeroDist[2]) &&
                           (normedZeroDist[3]>normedZeroDist[0] && normedZeroDist[3]>normedZeroDist[2]))
                        {
                            edgesToAdd[0] = edge(pt[1],pt[2]);
                            edgesToAdd[1] = edge(pt[3],pt[0]);
                        }
                        else
                            FatalErrorInFunction<<"Invalid distance. Can not determine the correct two edges !"<<normedZeroDist<<exit(FatalError);                            
                    }
                    else
                        FatalErrorInFunction<<"Whatever happens here should have been treated above. This can not happen!"<<exit(FatalError);
                }
                
                if(edgesToAdd[0][0] == edgesToAdd[0][1] || edgesToAdd[1][0] == edgesToAdd[1][1])
                {
                    Info<<endl;
                    Info<<"edge 1: "<<edgesToAdd[0]<<endl;
                    Info<<"edge 2: "<<edgesToAdd[1]<<endl;
                    FatalErrorInFunction<< "Edge has duplicate vertices"<< exit(FatalError);
                }
                
                edgeToFaces_.append(DynamicList<label>(0));                
                newMeshEdges_.append(edgesToAdd[0]);
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
                
                edgeToFaces_.append(DynamicList<label>(0));
                newMeshEdges_.append(edgesToAdd[1]);
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
                */
                
                edgeToFaces_.append(DynamicList<label>(0));                
                newMeshEdges_.append(edge(pt[0],pt[1]));
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
                
                edgeToFaces_.append(DynamicList<label>(0));
                newMeshEdges_.append(edge(pt[1],pt[2]));
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
                
                edgeToFaces_.append(DynamicList<label>(0));                
                newMeshEdges_.append(edge(pt[2],pt[3]));
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
                
                edgeToFaces_.append(DynamicList<label>(0));
                newMeshEdges_.append(edge(pt[3],pt[0]));
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
                
                problematicFace[i] = true;
                problematicFacePoints[i] = 4;
                problematicFaceNewPoints[i] = 4;
            }
        }
        else if(thisFacePoints.size() == 3)
        {
            if(i==589240)
            {
                Info<<"nbrOfPrevPoints: "<<nbrOfPrevPoints<<endl;
                Info<<"thisFacePoints: "<<thisFacePoints<<endl;
                //FatalErrorInFunction<<"Temp Stop!"<< exit(FatalError);
            }
            if(basisFaces[i].size() != 4)
            {
                FatalErrorInFunction<< "Faces with "<<basisFaces[i].size()<<" vertices not implemented!"<<exit(FatalError);
            }
            bool allPointsOld = true;
            label newPointsNum = 0;
            for(int k=0;k<thisFacePoints.size();k++)
            {
                if(thisFacePoints[k] >= nbrOfPrevPoints)
                {
                    allPointsOld = false;
                    newPointsNum++;
                }
            }
            if(allPointsOld)
            {
                if(basisFaces[i].size() != 4)
                {
                    FatalErrorInFunction<< "Faces with "<<basisFaces[i].size()<<" vertices not implemented!"<<exit(FatalError);
                }
                label fourthVertice = -1;
                std::unordered_set<label> cutPoints;
                for(int k=0;k<thisFacePoints.size();k++)
                {
                    cutPoints.insert(thisFacePoints[k]);
                }
                for(int k=0;k<basisFaces[i].size();k++)
                {
                    if(cutPoints.count(basisFaces[i][k])==0)
                    {
                        fourthVertice = basisFaces[i][k];
                        break;
                    }
                }
                if(fourthVertice==-1)
                {
                    FatalErrorInFunction<<"No fourth vertice assigned"<<exit(FatalError);
                }
                edge addingEdge = edge(-1,-1);
                bool edgeFound = false;
                for(int k=0;k<thisFacePoints.size();k++)
                {
                    label localIndx = basisFaces[i].which(thisFacePoints[k]);
                    label nextGlobalIndx = basisFaces[i].nextLabel(localIndx);
                    label prevGlobalIndx = basisFaces[i].prevLabel(localIndx);
                    if(thisFacePoints[k]!=fourthVertice && nextGlobalIndx!=fourthVertice && prevGlobalIndx!=fourthVertice)
                    {
                        addingEdge = edge(prevGlobalIndx,nextGlobalIndx);
                        edgeFound = true;
                    }
                }
                if(!edgeFound)
                {
                    Info<<endl;
                    Info<<"basisFaces[i]: "<<basisFaces[i]<<endl;
                    Info<<"thisFacePoints: "<<thisFacePoints<<endl;
                    Info<<"fourthVertice: "<<fourthVertice<<endl;
                    Info<<"addingEdge: "<<addingEdge<<endl;
                    FatalErrorInFunction<<"No added edge found"<<exit(FatalError);
                }
                edgeToFaces_.append(DynamicList<label>(0));                
                newMeshEdges_.append(addingEdge);
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);

                problematicFace[i] = true;
                problematicFacePoints[i] = 3;
                problematicFaceNewPoints[i] = 0;
            }
            else if(newPointsNum == 1)
            {
                if(i==589240)
                {
                    Info<<"nbrOfPrevPoints: "<<nbrOfPrevPoints<<endl;
                    Info<<"thisFacePoints: "<<thisFacePoints<<endl;
                    //FatalErrorInFunction<<"Temp Stop!"<< exit(FatalError);
                }
                if(basisFaces[i].size() != 4)
                {
                    FatalErrorInFunction<< "Faces with "<<basisFaces[i].size()<<" vertices not implemented!"<<exit(FatalError);
                }
                label newPoint = -1;
                for(int k=0;k<thisFacePoints.size();k++)
                {
                    if(thisFacePoints[k] >= nbrOfPrevPoints)
                    {
                        newPoint = thisFacePoints[k];
                        break;
                    }
                }
                edgeList addedEdges(2);
                label ind = 0;
                for(int k=0;k<thisFacePoints.size();k++)
                {
                    if(thisFacePoints[k] != newPoint)
                    {
                        addedEdges[ind] = edge(newPoint,thisFacePoints[k]);
                        ind++;
                    }
                }
                edgeToFaces_.append(DynamicList<label>(0));                
                newMeshEdges_.append(addedEdges[0]);
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
                
                edgeToFaces_.append(DynamicList<label>(0));                
                newMeshEdges_.append(addedEdges[1]);
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
            
                problematicFace[i] = true;
                problematicFacePoints[i] = 3;
                problematicFaceNewPoints[i] = 1;
            }
            else if(newPointsNum == 2)
            {
                edge addedEdge1(thisFacePoints[0],thisFacePoints[1]);
                edge addedEdge2(thisFacePoints[1],thisFacePoints[2]);
                edge addedEdge3(thisFacePoints[2],thisFacePoints[0]);
                
                edgeToFaces_.append(DynamicList<label>(0));                
                newMeshEdges_.append(addedEdge1);
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
                
                edgeToFaces_.append(DynamicList<label>(0));                
                newMeshEdges_.append(addedEdge2);
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
                
                edgeToFaces_.append(DynamicList<label>(0));                
                newMeshEdges_.append(addedEdge3);
                edgesToSide_.append(0);
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
                
                problematicFace[i] = true;
                problematicFacePoints[i] = 3;
                problematicFaceNewPoints[i] = 2;
            }
            else if(newPointsNum == 3)
            {
                FatalErrorInFunction
                << "Face with three cut point and three are new! This can not happen "<<exit(FatalError);
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
                    edgeToFaces_.append(DynamicList<label>(0));
                    edgeToFaces_[edgeToFaces_.size()-1].append(i);
                    
                    problematicFace[i] = true;
                    problematicFacePoints[i] = 2;
                    problematicFaceNewPoints[i] = 0;
                }
            }
            //one or more new point
            else
            {
                newMeshEdges_.append(edge(pt0,pt1));
                edgesToSide_.append(0);
                edgeToFaces_.append(DynamicList<label>(0));
                edgeToFaces_[edgeToFaces_.size()-1].append(i);
            }
        }
        else if(thisFacePoints.size() == 1)
        {
            if(thisFacePoints[0]>=nbrOfPrevPoints)
                FatalErrorInFunction
                << "Face with one but three of them are old ones and one is new! This can not happen "<<exit(FatalError);
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
    Info<<"edge to face done"<<endl;
    */
    
    faceToEdges_.setSize(basisFaces.size());
    //Info<<faceToEdges_.size()<<endl;
    //Info<<newMeshEdges_.size()<<endl;
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        //Info<<"Edge "<<i<<endl;
        if(edgeToFaces_[i].size() != 0)
        {
            //Info<<"\t"<<"edgeToFaces.size():"<<edgeToFaces_[i].size()<<endl;
            for(int k=0;k<edgeToFaces_[i].size();k++)
            {
                faceToEdges_[edgeToFaces_[i][k]].append(i);
            }
        }
    }
    //Info<<"face to edge done"<<endl;
    
    //Info<<faceToEdges_.size()<<endl;
    
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

    
    edgeToCells_.setSize(newMeshEdges_.size());
    const labelList& owner = this->faceOwner();
    const labelList& neighbour = this->faceNeighbour();
    const labelListList& edgeCells = this->edgeCells();
    
    //Info<<"newMeshEdges_: "<<newMeshEdges_.size()<<endl;
    //Info<<"edgeToCells_: "<<edgeToCells_.size()<<endl;
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        //Info<<"Edge: "<<i<<endl;
        if(edgeToFaces_[i].size() != 0)
        {
            //Info<<"edgeToFaces_["<<i<<"]: "<<edgeToFaces_[i].size()<<endl;
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
                    << exit(FatalError);
                }
                label thisFace = edgeToFaces_[i][0];
                //Info<<"thisFace: "<<thisFace<<endl;
                //Info<<"edgeToCells_[i].size(): "<<edgeToCells_[i].size()<<endl;
                //Info<<"Owner of thisFace: "<<owner[thisFace]<<endl;
                edgeToCells_[i].append(owner[thisFace]);
                //Info<<"Peter Pan"<<endl;
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
    
    //Info<<"edge to cell done"<<endl;

    
    cellToEdges_.setCapacity(meshCells.size());
    //Info<<"nCells: "<<meshCells.size()<<endl;
    //Info<<"newMeshEdges_: "<<newMeshEdges_.size()<<endl;
    //Info<<"edgeToCells_: "<<edgeToCells_.size()<<endl;
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        if(edgeToCells_[i].size() != 0)
        {
            //Info<<"edgeToCells_["<<i<<"]: "<<edgeToCells_[i].size()<<endl;
            for(int k=0;k<edgeToCells_[i].size();k++)
            {
                cellToEdges_[edgeToCells_[i][k]].append(i);
            }
        }
    }
    //Info<<"cell to edge done"<<endl;
    
    newMeshEdges_.setCapacity(newMeshEdges_.size());
    edgesToSide_.setCapacity(edgesToSide_.size());
    edgeToFaces_.setCapacity(edgeToFaces_.size());
    faceToEdges_.setCapacity(faceToEdges_.size());
    edgeToCells_.setCapacity(edgeToCells_.size());
    
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        edge oneEdge = newMeshEdges_[i];
        if(oneEdge[0] == oneEdge[1])
        {
            Info<<endl;
            Info<<"i:"<<i<<endl;
            Info<<"edge: "<<oneEdge<<endl;
            FatalErrorInFunction<< "Edge has duplicate vertices"<< exit(FatalError);
        }
    }
    
    List<DynamicList<DynamicList<label>>> cellNonconnectedEdges(meshCells.size());
    List<DynamicList<DynamicList<DynamicList<label>>>> cellNonConnectedMultiPoints(meshCells.size());
    List<DynamicList<DynamicList<DynamicList<label>>>> cellNonConnectedMultiEdges(meshCells.size());
    List<DynamicList<DynamicList<face>>> cellNonConnectedMultiFaces(meshCells.size());
    List<DynamicList<DynamicList<std::unordered_set<label>>>> cellNonConnectedMultiPointMap(meshCells.size());
    List<DynamicList<DynamicList<std::unordered_set<label>>>> cellNonConnectedMultiEdgeMap(meshCells.size());
    for(int i=0;i<meshCells.size();i++)
    // iterate across all cells
    {
        Info<<"cell:"<<i<<" - 1"<<endl;;
        labelList oneCellEdges = cellToEdges_[i];
        if(oneCellEdges.size()>0)
            Info<<"oneCellEdges:"<<oneCellEdges<<endl;
        
    //Seperate edges into chunks non connected by points
    //Start 
        std::unordered_set<label> usedEdges;
        for(int k=0;k<oneCellEdges.size();k++)
        // iterate until all edges are blocked
        {
            if(usedEdges.count(oneCellEdges[k])==0)
            {                        
                cellNonconnectedEdges[i].append(DynamicList<label>());
                DynamicList<label> nextEdgesInd;
                nextEdgesInd.append(oneCellEdges[k]);
                usedEdges.insert(oneCellEdges[k]);
                cellNonconnectedEdges[i].last().append(oneCellEdges[k]);
                for(int l=0;l<oneCellEdges.size();l++)
                {
                    DynamicList<label> frontEdges;
                    for(int m=0;m<oneCellEdges.size();m++)
                    // iteration to cover each edge |edges| times
                    {
                        if(usedEdges.count(oneCellEdges[m])==0)
                        {
                            label tryEdgeInd = oneCellEdges[m];
                            edge tryEdge = newMeshEdges_[tryEdgeInd];
                            for(int n=0;n<nextEdgesInd.size();n++)
                            {
                                edge currEdge = newMeshEdges_[nextEdgesInd[n]];
                                if(currEdge.connected(tryEdge))
                                {
                                    frontEdges.append(tryEdgeInd);
                                    break;
                                }
                            }
                        }
                    }
                    cellNonconnectedEdges[i].last().append(frontEdges);
                    nextEdgesInd = frontEdges;
                    for(int n=0;n<frontEdges.size();n++)
                    {
                        usedEdges.insert(frontEdges[n]);
                    }
                }
            }
        }
        for(int k=0;k<oneCellEdges.size();k++)
        {
            if(usedEdges.count(oneCellEdges[k])==0)
            {
                if(problematicFacePoints[i]==-1 && problematicFaceNewPoints[i] == -1)
                {}
                else
                    FatalErrorInFunction<< "Non used edge remains!"<< exit(FatalError);
            }
        }
        if(cellNonconnectedEdges[i].size()>0)
            Info<<"cellNonconnectedEdges["<<i<<"]:"<<cellNonconnectedEdges[i]<<endl;
    //End

        Info<<"cell:"<<i<<" - 2"<<endl;
    //Test chunks of edges
    //Start
        if(cellNonconnectedEdges[i].size() == 0)
        {
            //FatalErrorInFunction<< "Problematic face in cell but no edges colllected!"<< exit(FatalError);
        }
        else if(cellNonconnectedEdges[i].size() > 1)
        {
            /*
            Info<<"---------------------------------------------"<<endl;
            Info<<"Face: "<<i<<endl;
            Info<<"Owner Cell: "<<faceCells[0]<<" edges: "<<cellToEdges_[faceCells[0]]<<endl;
            Info<<"Neighbour Cell: "<<faceCells[1]<<" edges: "<<cellToEdges_[faceCells[1]]<<endl;
            Info<<"edgeList:"<<cellNonconnectedEdges<<endl;
            Info<<"---------------------------------------------"<<endl;
            */
        }
        else
        {
            std::list<label> edges;
            for(int k=0;k<cellNonconnectedEdges[i][0].size();k++)
            {
                edges.push_back(cellNonconnectedEdges[i][0][k]);
            }
            DynamicList<label> oneFace;
            std::unordered_set<label> oneFaceBlckd;
            oneFace.append(edges.front());
            oneFaceBlckd.insert(edges.front());
            edges.remove(edges.front());
            while(edges.size()!=0)
            {
                bool isConnected = false;
                DynamicList<label> connectedToFace;
                for(int k=0;k<oneFace.size();k++)
                {
                    label faceEdgeInd = oneFace[k];
                    edge faceEdge = newMeshEdges_[faceEdgeInd];
                    for(auto l=edges.cbegin();l!=edges.cend();l++)
                    {
                        label testEdgeInd = *l;
                        edge testEdge = newMeshEdges_[testEdgeInd];
                        if(faceEdgeInd==testEdgeInd)
                            FatalErrorInFunction<< "Invalid"<< exit(FatalError);
                        if(testEdge.connected(faceEdge))
                        {
                            connectedToFace.append(testEdgeInd);
                            isConnected = true;
                            break;
                        }
                    }
                    if(isConnected)
                        break;
                }
                if(!isConnected)
                {
                    FatalErrorInFunction<< "Edges do not build a connected face"<< exit(FatalError);
                }
                else
                {
                    oneFace.append(connectedToFace);
                    for(int k=0;k<connectedToFace.size();k++)
                    {
                        edges.remove(connectedToFace[k]);
                    }
                }
            }
        }
    //End
        
        Info<<"cell:"<<i<<" - 3"<<endl;
    //Seperate chunks of edges into chunks non connected by edges
    //Seperate chunks of edges into faces
    //Start
        for(int k=0;k<cellNonconnectedEdges[i].size();k++)
        // seperate edge connections
        {
            Info<<"cellNonconnectedEdges["<<i<<"]:"<<cellNonconnectedEdges[i]<<endl;
            Info<<"cell:"<<i<<" - 3  "<<k<<" - 1"<<endl;

        //Seperate edge chunks into closed faces
        //Start
            DynamicList<label>& edgeConnection = cellNonconnectedEdges[i][k];
            std::unordered_set<label> pointMap;
            std::unordered_map<label,std::unordered_set<label>> pointGraphData;
            DynamicList<std::pair<label,label>> pointEdgeComb;
            DynamicList<label> pointList;
            for(int l=0;l<edgeConnection.size();l++)
            {
                edge oneEdge = newMeshEdges_[edgeConnection[l]];
                for(int m=0;m<oneEdge.size();m++)
                {
                    label point = oneEdge[m];
                    if(pointMap.count(point)==0)
                    {
                        pointMap.insert(point);
                        pointList.append(point);
                    }
                    label connectedPoint = oneEdge.otherVertex(point);
                    if(connectedPoint==-1)
                        FatalErrorInFunction<< "Can not happen"<< exit(FatalError);
                    std::pair<label,label> pointEdge(connectedPoint,edgeConnection[l]);
                    pointEdgeComb.append(pointEdge);
                    pointGraphData[point].insert(pointEdgeComb.size()-1);
                }
            }
            Info<<"cell:"<<i<<" - 3  "<<k<<" - 2"<<endl;
            DynamicList<DynamicList<label>> closedCyclePoints;
            DynamicList<std::unordered_set<label>> closedCycleEdges;
            DynamicList<DynamicList<label>> closedCycleEdgesList;
            for(const std::pair<label,std::unordered_set<label>>& n : pointGraphData ) 
            {
                DynamicList<label> cyclePath;
                DynamicList<label> cycleEdgePath;
                std::unordered_set<label> coveredPoints;
                std::unordered_set<label> usedEdges;
                findCycles(i,n.first,-1,-1,cyclePath,cycleEdgePath,coveredPoints,usedEdges,pointGraphData,pointEdgeComb,closedCyclePoints,closedCycleEdges,closedCycleEdgesList);
            }
            if(closedCyclePoints.size() > 3)
                FatalErrorInFunction<< "Temp stop"<< exit(FatalError);
        //End
            
            Info<<"cell:"<<i<<" - 3  "<<k<<" - 3"<<endl;
            Info<<"closedCycleEdgesList:"<<closedCycleEdgesList<<endl;
            
        //Speperate chunks of faces only connected by one point
        //Start
            List<DynamicList<label>> edgeConnectedFaces(closedCycleEdges.size());
            //Build graph of connected closed cycles
            for(int l=0;l<closedCycleEdges.size();l++)
            {
                for(int m=0;m<closedCycleEdges.size();m++)
                {
                    if(l==m)
                    {
                        edgeConnectedFaces[l].append(l);
                    }
                    else
                    {
                        bool connected = false;
                        for(int n=0;n<closedCycleEdgesList[m].size();n++)
                        {
                            if(closedCycleEdges[l].count(closedCycleEdgesList[m][n])==0)
                            {
                                if(!connected)
                                {
                                    edgeConnectedFaces[l].append(m);
                                    connected = true;
                                }
                                else
                                {
                                    //Info<<"closedCycleEdgesList[m][n]:"<<closedCycleEdgesList[m][n]<<endl;
                                    //FatalErrorInFunction<< "Ducplicate share edge!"<< exit(FatalError);
                                }
                            }
                        }
                    }
                }
            }
            Info<<"cell:"<<i<<" - 3  "<<k<<" - 4"<<endl;
            Info<<"edgeConnectedFaces:"<<edgeConnectedFaces<<endl;

            DynamicList<DynamicList<label>> connectedFaces;
            std::queue<label> nextCycles;
            std::unordered_set<label> takenCycles;
            for(int l=0;l<edgeConnectedFaces.size();l++)
            {
                Info<<"connectedFaces:"<<connectedFaces<<endl;
                Info<<"cell:"<<i<<" - 3  "<<k<<" - 4  "<<l<<endl;
                if(takenCycles.count(l)==0)
                {
                    connectedFaces.append(DynamicList<label>());
                    connectedFaces.last().append(l);
                    takenCycles.insert(l);
                    Info<<"cell:"<<i<<" - 3  "<<k<<" - 4  "<<l<<" - 5 1"<<endl;
                    for(int m=0;m<edgeConnectedFaces[l].size();m++)
                    {
                        if(takenCycles.count(edgeConnectedFaces[l][m])==0)
                        {
                            nextCycles.push(edgeConnectedFaces[l][m]);
                            Info<<"push: "<<edgeConnectedFaces[l][m]<<endl;
                        }
                    }
                    Info<<"cell:"<<i<<" - 3  "<<k<<" - 4  "<<l<<" - 5 2"<<" : empty:"<<nextCycles.empty()<<endl;
                    Info<<"connectedFaces:"<<connectedFaces<<endl;
                    while(!nextCycles.empty())
                    {
                        label faceInd = nextCycles.front();
                        for(int m=0;m<edgeConnectedFaces[faceInd].size();m++)
                        {
                            if(takenCycles.count(edgeConnectedFaces[faceInd][m])==0)
                            {
                                nextCycles.push(edgeConnectedFaces[faceInd][m]);
                                Info<<"push: "<<edgeConnectedFaces[faceInd][m]<<endl;
                            }
                        }
                        if(takenCycles.count(faceInd)==0)
                        {
                            connectedFaces.last().append(faceInd);
                            takenCycles.insert(faceInd);
                        }
                        nextCycles.pop();
                        Info<<"cell:"<<i<<" - 3  "<<k<<" - 4  "<<l<<" - 5 3 "<<faceInd<<" : empty:"<<nextCycles.empty()<<endl;
                    }
                }
            }
            Info<<"cell:"<<i<<" - 3  "<<k<<" - 5"<<endl;
            DynamicList<DynamicList<DynamicList<label>>> edgeFaces;
            DynamicList<DynamicList<DynamicList<label>>> pointFaces;
            
            Info<<"connectedFaces:"<<connectedFaces<<endl;

            for(int l=0;l<connectedFaces.size();l++)
            {
                edgeFaces.append(DynamicList<DynamicList<label>>());
                pointFaces.append(DynamicList<DynamicList<label>>());
                for(int m=0;m<connectedFaces[l].size();m++)
                {
                    edgeFaces.last().append(closedCycleEdgesList[connectedFaces[l][m]]);
                    pointFaces.last().append(closedCyclePoints[connectedFaces[l][m]]);
                }
            }
            Info<<"edgeFaces"<<edgeFaces<<endl;
            Info<<"pointFaces"<<pointFaces<<endl;
            Info<<"cell:"<<i<<" - 3  "<<k<<" - 6"<<endl;
            for(int l=0;l<edgeFaces.size();l++)
            {
                cellNonConnectedMultiEdges[i].append(edgeFaces[l]);
                cellNonConnectedMultiPoints[i].append(pointFaces[l]);
            }
            if(i==198219)
            {
                Info<<"cellNonConnectedMultiEdges[i]:"<<cellNonConnectedMultiEdges[i]<<endl;
                Info<<"cellNonConnectedMultiPoints[i]:"<<cellNonConnectedMultiPoints[i]<<endl;
                //FatalErrorInFunction<< "Temp stop"<< exit(FatalError);
            }
        //End
        }
    //End
    }
    Info<<"cellNonConnectedMultiEdges[198219]:"<<cellNonConnectedMultiEdges[198219]<<endl;
    Info<<"cellNonConnectedMultiPoints[198219]:"<<cellNonConnectedMultiPoints[198219]<<endl;
    for(int i=0;i<meshCells.size();i++)
    {
        for(int j=0;j<cellNonConnectedMultiPoints[i].size();j++)
        {
            cellNonConnectedMultiFaces[i].append(DynamicList<face>());
            cellNonConnectedMultiPointMap[i].append(DynamicList<std::unordered_set<label>>());
            cellNonConnectedMultiEdgeMap[i].append(DynamicList<std::unordered_set<label>>());
            for(int k=0;k<cellNonConnectedMultiPoints[i][j].size();k++)
            {
                bool faceInFace = false;
                DynamicList<label> usedPoints;
                for(int l=0;l<cellNonConnectedMultiPoints[i][j][k].size();l++)
                {
                    label cutFacePnt = cellNonConnectedMultiPoints[i][j][k][l];
                    Info<<"cutFacePnt:"<<cutFacePnt<<endl;
                    Info<<"nbrOfPrevPoints:"<<nbrOfPrevPoints<<endl;
                    Info<<"basisPoints.size():"<<basisPoints.size()<<endl;
                    Info<<"face:"<<face(cellNonConnectedMultiPoints[i][j][k])<<endl;
                    for(int m=0;m<cellNonConnectedMultiPoints[i][j][k].size();m++)
                    {
                        Info<<cellNonConnectedMultiPoints[i][j][k][m]<<":"<<pointToEgde_[cellNonConnectedMultiPoints[i][j][k][m]]<<endl;
                    }                    
                    
                    if(cutFacePnt<nbrOfPrevPoints)
                    {
                        usedPoints.append(cutFacePnt);
                    }
                    else
                    {
                        label pntEdgeInd = pointToEgde_[cutFacePnt];
                        if(pntEdgeInd==-1)
                            FatalErrorInFunction<<"Can not happen!"<<exit(FatalError);
                        edge pntEdge = newMeshEdges_[pntEdgeInd];
                        usedPoints.append(pntEdge.start());
                        usedPoints.append(pntEdge.end());
                    }
                }
                cell thisCell = meshCells[i];
                for(int l=0;l<thisCell.size();l++)
                {
                    face thisFace = basisFaces[thisCell[l]];
                    bool allInFace = true;
                    for(int m=0;m<usedPoints.size();m++)
                    {
                        if(thisFace.which(usedPoints[m])==-1)
                            allInFace = false;
                    }
                    if(allInFace)
                        faceInFace = true;
                }
                if(!faceInFace)
                {
                    cellNonConnectedMultiFaces[i].last().append(face(cellNonConnectedMultiPoints[i][j][k]));
                    cellNonConnectedMultiPointMap[i].last().append(std::unordered_set<label>());
                    cellNonConnectedMultiEdgeMap[i].last().append(std::unordered_set<label>());
                    for(int l=0;l<cellNonConnectedMultiPoints[i][j][k].size();l++)
                    {
                        cellNonConnectedMultiPointMap[i].last().last().insert(cellNonConnectedMultiPoints[i][j][k][l]);
                        cellNonConnectedMultiEdgeMap[i].last().last().insert(cellNonConnectedMultiEdges[i][j][k][l]);
                    }
                }
            }
        }
    }
    //FatalErrorInFunction<< "Temp stop"<< exit(FatalError);
    
    const labelListList& pointToCells = this->pointCells();
    const labelListList& edgeToCells = this->edgeCells();
    DynamicList<label> addedEdgeToDelete;
    for(int i=0;i<faceToEdges_.size();i++)
    {
        if(problematicFace[i])
        {
            labelList edgeInd = faceToEdges_[i];
            label faceOwner = cellOwner[i];
            label faceNeighbour = -1;
            if(i<cellNeighbor.size())
                faceNeighbour = cellNeighbor[i];        
            labelList faceCells(0);
            faceCells.append(faceOwner);
            if(i<cellNeighbor.size())
            {
                faceCells.append(cellNeighbor[i]);
            }
            
        /*  Compute information if vertical edges of faces vertives are cut
         *  Start
         */
            List<List<bool>> verticalEdgesAreCut(4,List<bool>(faceCells.size(),false));
            List<List<label>> verticalEdges(4,List<label>(faceCells.size()));
            std::unordered_set<label> edgesMapOfFace;
            List<std::unordered_set<label>> edgesOfCellNeighbor(faceCells.size());
            for(int j=0;j<faceToEdge[i].size();j++)
                edgesMapOfFace.insert(faceToEdge[i][j]);
            for(int j=0;j<edgesOfCellNeighbor.size();j++)
                for(int k=0;k<cellToEdge[faceCells[j]].size();k++)
                    edgesOfCellNeighbor[j].insert(cellToEdge[faceCells[j]][k]);
            if(basisFaces[i].size()!=4)
                FatalErrorInFunction<<"Faces with other than 4 vertices not allowed!"<< exit(FatalError);
            for(int j=0;j<verticalEdges.size();j++)
            {
                label pnt = basisFaces[i][j];
                const labelList& onePointEdges = pointToEdge[pnt];
                for(int k=0;k<faceCells.size();k++)
                {
                    bool verticalEdgeFound = false;
                    label verticalEdge = -1;
                    for(int l=0;l<onePointEdges.size();l++)
                    {
                        if(edgesMapOfFace.count(onePointEdges[l])==0)
                        {
                            if(edgesOfCellNeighbor[k].count(onePointEdges[l])!=0)
                            {
                                if(verticalEdgeFound)
                                    FatalErrorInFunction<<"Double edge. Can not happen!"<< exit(FatalError);
                                verticalEdge = onePointEdges[l];
                                verticalEdgeFound=true;
                            }
                        }
                    }
                    if(!verticalEdgeFound)
                    {
                        Info<<endl;
                        Info<<"Face "<<i<<" face:"<<basisFaces[i]<<endl;
                        Info<<"Point:"<<pnt<<endl;
                        Info<<"onePointEdges:"<<onePointEdges<<endl;
                        Info<<"faceToEdge:"<<faceToEdge[i]<<endl;
                        FatalErrorInFunction<<"No edge found. Can not happen!"<< exit(FatalError);
                    }
                    verticalEdges[j][k]=verticalEdge;
                }
            }
            for(int j=0;j<verticalEdgesAreCut.size();j++)
            {
                for(int k=0;k<verticalEdgesAreCut[j].size();k++)
                {
                    if(edgeToPoint_[verticalEdges[j][k]]!=-1)
                        verticalEdgesAreCut[j][k] = true;
                }
            }
        /*  End
         */
            
            if(i==589240)
            {
                Info<<"Face: "<<i<<endl;
                Info<<"Owner Cell: "<<faceCells[0]<<" edges: "<<cellToEdges_[faceCells[0]]<<endl;
                for(int j=0;j<cellToEdges_[faceCells[0]].size();j++)
                {
                    Info<<" "<<cellToEdges_[faceCells[0]][j]<<" "<<newMeshEdges_[cellToEdges_[faceCells[0]][j]]<<endl;
                }
                Info<<"Neighbour Cell: "<<faceCells[1]<<" edges: "<<cellToEdges_[faceCells[1]]<<endl;
                for(int j=0;j<cellToEdges_[faceCells[1]].size();j++)
                {
                    Info<<" "<<cellToEdges_[faceCells[1]][j]<<" "<<newMeshEdges_[cellToEdges_[faceCells[1]][j]]<<endl;
                }
                Info<<"edgeList:"<<cellNonconnectedEdges[i]<<endl;
            }
            if(i==589240)
            {
                Info<<"problematicFacePoints:"<<problematicFacePoints[i]<<endl;
                Info<<"problematicNewFacePoints:"<<problematicFaceNewPoints[i]<<endl;
                Info<<"closedCycles:";
                Info<<cellNonConnectedMultiPoints[faceOwner]<<endl;
                Info<<"closedEdgeCycle:";
                Info<<cellNonConnectedMultiEdges[faceOwner]<<endl;
                Info<<"Face To Edges:"<<faceToEdges_[i]<<endl;
                for(int m=0;m<faceToEdges_[i].size();m++)
                {
                    bool foundDist;
                    edge oneEdge = newMeshEdges_[faceToEdges_[i][m]];
                    point edgeCenter = 0.5*(newMeshPoints_[oneEdge.start()]+newMeshPoints_[oneEdge.end()]);
                    Info<<"edge: "<<faceToEdges_[i][m]<<oneEdge<<" ("<<newMeshPoints_[oneEdge.start()]<<"  "<<newMeshPoints_[oneEdge.end()]<<")"<<endl;
                    Info<<"center:"<<edgeCenter<<distToNurbs(edgeCenter,foundDist)<<endl;
                }
                Info<<endl;
            }

        /*  Compute face connection for cut edges of faces
         *  Start
         */
            List<DynamicList<std::unordered_set<label>>> edgesOfEdgeConnectedFacesPerCell(faceCells.size());    // HashMap of the edge index of all edges in a face cluster
            for(int nC = 0;nC<faceCells.size();nC++)
            {
                edgesOfEdgeConnectedFacesPerCell[nC].setSize(cellNonConnectedMultiEdges[faceCells[nC]].size());
                for(int j=0;j<cellNonConnectedMultiEdges[faceCells[nC]].size();j++)
                //Iterate across clusters of faces
                {
                    for(int k=0;k<cellNonConnectedMultiEdges[faceCells[nC]][j].size();k++)
                    //Iterate across faces in clusters
                    {
                        for(int l=0;l<cellNonConnectedMultiEdges[faceCells[nC]][j][k].size();l++)
                        //Iterate across edges in faces
                        {
                            edgesOfEdgeConnectedFacesPerCell[nC][j].insert(cellNonConnectedMultiEdges[faceCells[nC]][j][k][l]);
                        }
                    }
                }
            }
            label nbrEdges;
            if(problematicFacePoints[i]==4 || problematicFacePoints[i]==3)
                nbrEdges = problematicFacePoints[i];
            else
                nbrEdges = 1;
            List<List<DynamicList<bool>>> connectedToFacePerEdge(nbrEdges); // Bool Map zeroedges - cell - clusterofFaces: true if connected
            for(int j=0;j<connectedToFacePerEdge.size();j++)
            {
                connectedToFacePerEdge[j].setSize(faceCells.size());
                for(int k=0;k<connectedToFacePerEdge[j].size();k++)
                {
                    for(int l=0;l<edgesOfEdgeConnectedFacesPerCell[k].size();l++)
                    {
                        if(edgesOfEdgeConnectedFacesPerCell[k][l].count(edgeInd[j])==0)
                            connectedToFacePerEdge[j][k].append(false);
                        else
                            connectedToFacePerEdge[j][k].append(true);
                    }
                    if(connectedToFacePerEdge[j][k].size()>1)
                    {
                        bool oneTrue = false;
                        for(int l=0;l<connectedToFacePerEdge[j][k].size();l++)
                        {
                            if(connectedToFacePerEdge[j][k][l])
                            {
                                if(!oneTrue)
                                    oneTrue = true;
                                else
                                    FatalErrorInFunction<< "Can not happen!"<< exit(FatalError);

                            }
                        }
                    }
                }
            }
            List<List<bool>> connectedToCellPerEdge(nbrEdges); // Bool Map zeroedges - cell: true if connected
            for(int j=0;j<connectedToCellPerEdge.size();j++)
            {
                connectedToCellPerEdge[j].setSize(faceCells.size());
                for(int k=0;k<connectedToCellPerEdge[j].size();k++)
                {
                    bool connected = false;
                    for(int l=0;l<connectedToFacePerEdge[j][k].size();l++)
                    {
                        connected = connected || connectedToFacePerEdge[j][k][l];
                    }
                    connectedToCellPerEdge[j][k] = connected;
                }
            }
        /*  End
         */
            
            if(i==589240)
            {
                for(int c=0;c<edgesOfEdgeConnectedFacesPerCell.size();c++)
                {
                    Info<<"Cell:"<<faceCells[c]<<" "<<endl;
                    for(int patch=0;patch<edgesOfEdgeConnectedFacesPerCell[c].size();patch++)
                    {
                        Info<<" - - - "<<"patch:"<<patch;
                        Info<<"{ ";
                        for(label ind : edgesOfEdgeConnectedFacesPerCell[c][patch]) {Info<<" "<<ind<<" ";}
                        Info<<" }"<<endl;
                    }
                }
                for(int e=0;e<connectedToFacePerEdge.size();e++)
                {
                    Info<<"edge:"<<edgeInd[e]<<" "<<endl;
                    for(int c=0;c<connectedToFacePerEdge[e].size();c++)
                    {
                        Info<<" -- -- -- Cell:"<<faceCells[c]<<" "<<endl;
                        for(int patch=0;patch<connectedToFacePerEdge[e][c].size();patch++)
                        {
                            Info<<" -- -- -- -- -- -- "<<"patch:"<<connectedToFacePerEdge[e][c][patch]<<endl;
                        }
                    }
                }
                for(int e=0;e<connectedToCellPerEdge.size();e++)
                {
                    Info<<"edge:"<<edgeInd[e]<<" "<<endl;
                    for(int c=0;c<connectedToCellPerEdge[e].size();c++)
                    {
                        Info<<" -- -- -- Cell:"<<faceCells[c]<<" "<<"patch:"<<connectedToCellPerEdge[e][c]<<endl;
                    }
                }
            }
            
         /*  Correct Edge computation by using distance vectors
          *  Start
          */
         /*  Remove-------------------------
            List<bool> correctEdge(edgeInd.size());
            edgeList edgesOfFace(edgeInd.size());
            List<List<vector>> vectorsOfEdges(edgeInd.size());
            List<List<scalar>> distOfEdges(edgeInd.size());
            for(int j=0;j<edgeInd.size();j++)
            {
                edgesOfFace[j] = newMeshEdges_[edgeInd[j]];
                vectorsOfEdges[j] = vectorsToNurbsOfEdge(newMeshPoints_[edgesOfFace[j].start()],newMeshPoints_[edgesOfFace[j].end()],11);
                distOfEdges[j] = distToNursOfEdge(newMeshPoints_[edgesOfFace[j].start()],newMeshPoints_[edgesOfFace[j].end()],11);
                scalar p19 = vectorsOfEdges[j][1] && vectorsOfEdges[j][9];
                scalar p28 = vectorsOfEdges[j][2] && vectorsOfEdges[j][8];
                scalar p37 = vectorsOfEdges[j][3] && vectorsOfEdges[j][7];
                scalar p46 = vectorsOfEdges[j][4] && vectorsOfEdges[j][6];

                if(p19<=0 && p28<=0 && p37<=0 && p46<=0)
                //Vectors point in different directions
                {
                    correctEdge[j] = false;
                }
                else if(p19>0 && p28>0 && p37>0 && p46>0)
                //Vectors point in same direction
                {
                    correctEdge[j] = true;
                }
                else
                //Inconsistent edge vectors directions
                {
                    Info<<"Edge: "<<edgesOfFace[j]<<endl;
                    Info<<"Startpoint: "<<newMeshPoints_[edgesOfFace[j].start()]<<endl;
                    Info<<"Endpoint: "<<newMeshPoints_[edgesOfFace[j].end()]<<endl;
                    Info<<"vectorsOfEdges["<<j<<"]: "<<vectorsOfEdges[j]<<endl;
                    //Info<<"vectorsOfEdges: "<<vectorsOfEdges<<endl;
                    Info<<"distOfEdge["<<j<<"]:"<<distOfEdges[j]<<endl;
                    Info<<"p19: "<<p19<<endl;
                    Info<<"p28: "<<p28<<endl;
                    Info<<"p37: "<<p37<<endl;
                    Info<<"p46: "<<p46<<endl;
                    bool found;
                    label nurbsInd;
                    scalar nurbsPara;
                    scalar dis;
                    nearestNurbsIndexPara(newMeshPoints_[edgesOfFace[j].start()],found,nurbsInd,nurbsPara);
                    Info<<"Start ind:"<<nurbsInd<<" para:"<<nurbsPara<<endl;
                    nurbsPara = 0;
                    nearestNurbsIndexPara(newMeshPoints_[edgesOfFace[j].end()],found,nurbsInd,nurbsPara);
                    Info<<"End ind:"<<nurbsInd<<" para:"<<nurbsPara<<endl;
                    FatalErrorInFunction<< "Inconsistent edge vectors directions"<< exit(FatalError);                            
                }
            }
            Remove -------------------------
        */  
        /*  End
         */  

        /*  Data separation of edges and connection to vertices
         *  Start
         */
            label oldEdge = 0;
            label newEdge = 0;
            DynamicList<label> newEdgeLocalInd;
            DynamicList<label> oldEdgeLocalInd;
            labelList egdeLocalIndToNewEdgeListInd(edgeInd.size(),-1);
            labelList edgeLocalIndToOldEdgeListInd(edgeInd.size(),-1);
            for(int j=0;j<edgeInd.size();j++)
            {
                if(edgeInd[j]>=nbrOfPrevEdges)
                {
                    newEdge++;
                    newEdgeLocalInd = j;
                }
                else
                {
                    oldEdge++;
                    oldEdgeLocalInd.append(j);
                }
            }
            for(int j=0;j<newEdgeLocalInd.size();j++)
                egdeLocalIndToNewEdgeListInd[newEdgeLocalInd[j]]=j;
            for(int j=0;j<oldEdgeLocalInd.size();j++)
                edgeLocalIndToOldEdgeListInd[oldEdgeLocalInd[j]]=j;
            
            std::unordered_map<label,label> cutEdgeLocalIndToNonZeroPntLocalInd;
            std::unordered_map<label,label> nonZeroPntLocalIndTocutEdgeLocalInd;
            for(int j=0;j<newEdgeLocalInd.size();j++)
            {
                edge zeroEdge = newMeshEdges_[edgeInd[newEdgeLocalInd[j]]];
                List<DynamicList<label>> edgePntsEdgeGlobInds(2);
                labelList edgePnts(2);
                edgePnts[0] = zeroEdge.start();
                edgePnts[1] = zeroEdge.end();
                for(int k=0;k<2;k++)
                {
                    if(edgePnts[k]<nbrOfPrevPoints)
                        edgePntsEdgeGlobInds[k].append(pointToEgde_[edgePnts[k]]);
                    else
                        edgePntsEdgeGlobInds[k].append(pointToEdge[edgePnts[k]]);
                }
                label nonZeroPntLocalIndForCutEdge=-1;
                for(int k=0;k<edgePntsEdgeGlobInds[0].size();k++)
                {
                    edge startEdge = newMeshEdges_[edgePntsEdgeGlobInds[0][k]];
                    for(int l=0;l<edgePntsEdgeGlobInds[1].size();l++)
                    {
                        edge endEdge = newMeshEdges_[edgePntsEdgeGlobInds[1][k]];
                        label comVertLocalInd = basisFaces[i].which(startEdge.commonVertex(endEdge));
                        if(comVertLocalInd!=-1)
                        {
                            if(nonZeroPntLocalIndForCutEdge==-1)
                                nonZeroPntLocalIndForCutEdge=comVertLocalInd;
                            else
                                FatalErrorInFunction<<"Double point. Can not happen!"<< exit(FatalError);
                        }
                    }
                }
                if(nonZeroPntLocalIndForCutEdge==-1)
                    FatalErrorInFunction<<"No point. Can not happen!"<< exit(FatalError);
                cutEdgeLocalIndToNonZeroPntLocalInd[j]=nonZeroPntLocalIndForCutEdge;
                nonZeroPntLocalIndTocutEdgeLocalInd[nonZeroPntLocalIndForCutEdge]=j;
            }
        /* End
         */
            
            if(problematicFacePoints[i]==4 && problematicFaceNewPoints[i]==4)
            {
                Info<<"--------------------------Move in 4-4-----------------------------------"<<endl;
                if(edgeInd.size()!=4)
                    FatalErrorInFunction<< "No four added edges"<< exit(FatalError);
                
                bool alternatingFacesToCellPerEdge = false;
                for(int j=0;j<faceCells.size();j++)
                {
                    bool oneOrderAlternatingFaces = true;
                    for(int k=0;k<connectedToCellPerEdge.size();k++)
                    {
                        oneOrderAlternatingFaces = oneOrderAlternatingFaces && connectedToCellPerEdge[k][(j+k)%faceCells.size()];
                    }
                    alternatingFacesToCellPerEdge = alternatingFacesToCellPerEdge || oneOrderAlternatingFaces;
                }
                bool alternatingFacesToCellPerEdgeStrict = false;
                for(int j=0;j<faceCells.size();j++)
                {
                    bool oneOrderAlternatingFacesStrict = true;
                    for(int k=0;k<connectedToCellPerEdge.size();k++)
                    {
                        oneOrderAlternatingFacesStrict = oneOrderAlternatingFacesStrict && (connectedToCellPerEdge[k][(j+k)%faceCells.size()]&&!connectedToCellPerEdge[k][(j+k+1)%faceCells.size()]);
                    }
                    alternatingFacesToCellPerEdgeStrict = alternatingFacesToCellPerEdgeStrict || oneOrderAlternatingFacesStrict;
                }
                
                bool mustBeAllFourEdges = false;
                /* all four edges if the face is no boundary face and
                 * there are alternating cut connections to cells
                 */
                if(faceCells.size()==2)
                {
                    if(alternatingFacesToCellPerEdge)
                    // Must be
                    {
                        if(alternatingFacesToCellPerEdgeStrict)
                        // Probable but wrong faces
                        {
                            mustBeAllFourEdges = true;
                        }
                    }
                }
                
                bool mustBeOnlyTwoEdges = false;
                labelList edgeLocalIndsToRemove(2,-1);
                /* only two edges if face has only two connections to the other cell
                 */
                if(!mustBeAllFourEdges)
                {
                    bool correctFaces13 = false;
                    bool correctNonFaces13 = false;
                    bool correctFaces02 = false;
                    bool correctNonFaces02 = false;
                    if(faceCells.size()==2)
                    {
                        if(connectedToCellPerEdge[0][0] && connectedToCellPerEdge[0][1] && connectedToCellPerEdge[2][0] && connectedToCellPerEdge[2][1])
                            correctFaces02 = true;
                        if(connectedToCellPerEdge[1][0] && connectedToCellPerEdge[1][1] && connectedToCellPerEdge[3][0] && connectedToCellPerEdge[3][1])
                            correctFaces13 = true;
                        if(!(connectedToCellPerEdge[0][0] && connectedToCellPerEdge[0][1] && connectedToCellPerEdge[2][0] && connectedToCellPerEdge[2][1]))
                            correctNonFaces13 = true;
                        if(!(connectedToCellPerEdge[1][0] && connectedToCellPerEdge[1][1] && connectedToCellPerEdge[3][0] && connectedToCellPerEdge[3][1]))
                            correctNonFaces02 = true;
                    }
                    else if(faceCells.size()==1)
                    {
                        if(connectedToCellPerEdge[0][0] && connectedToCellPerEdge[2][0])
                            correctFaces02 = true;
                        if(connectedToCellPerEdge[1][0] && connectedToCellPerEdge[3][0])
                            correctFaces13 = true;
                        if(!(connectedToCellPerEdge[0][0] && connectedToCellPerEdge[2][0]))
                            correctNonFaces13 = true;
                        if(!(connectedToCellPerEdge[1][0] && connectedToCellPerEdge[3][0]))
                            correctNonFaces02 = true;
                    }
                    if(correctFaces02 && correctNonFaces02)
                    {
                        mustBeOnlyTwoEdges = true;
                        edgeLocalIndsToRemove[0] = 1;
                        edgeLocalIndsToRemove[1] = 3;
                    }
                    else if(correctFaces13 && correctNonFaces13)
                    {
                        mustBeOnlyTwoEdges = true;
                        edgeLocalIndsToRemove[0] = 0;
                        edgeLocalIndsToRemove[1] = 2;
                    }
                    else if(!correctFaces02 && !correctFaces13)
                        FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                    else
                    {
                        DynamicList<label> zeroPoints;
                        for(int j=0;j<edgeInd.size();j++)
                        {
                            label firstEdge = j;
                            label secondEdge = (j+1)%edgeInd.size();
                            label cutPoint = newMeshEdges_[edgeInd[firstEdge]].commonVertex(newMeshEdges_[edgeInd[secondEdge]]);
                            if(cutPoint==-1 || cutPoint<nbrOfPrevPoints)
                                FatalErrorInFunction<<"Point between two in a four new edge face. Can not happen!"<< exit(FatalError);
                            zeroPoints.append(cutPoint);
                        }
                        if(zeroPoints.size()!=4)
                            FatalErrorInFunction<<"Must be four points. Can not happen!"<< exit(FatalError);
                        // four zero Points to list
                        List<DynamicList<DynamicList<face>>> zeroPointsClosedFaces(zeroPoints.size());
                        List<DynamicList<DynamicList<std::unordered_set<label>>>> zeroPointsClosedFaceMap(zeroPoints.size());
                        for(int j=0;j<zeroPoints.size();j++)
                        {
                            //collect all faces in neighboring cells connected to the point in posClosedFacesAroundPoint
                            label edgeInd = pointToEgde_[zeroPoints[j]];
                            const labelList& cellsAtPoint = edgeToCells[edgeInd];
                            label nbrRelCells = cellsAtPoint.size();
                            List<DynamicList<face>> posClosedFacesAroundPoint(nbrRelCells);
                            List<DynamicList<std::unordered_set<label>>> posClosedPointMapAroundPoint(nbrRelCells);
                            for(int k=0;k<cellsAtPoint.size();k++)
                            {
                                DynamicList<DynamicList<std::unordered_set<label>>> thisCellFaceGroupsMap = cellNonConnectedMultiPointMap[cellsAtPoint[k]];
                                DynamicList<DynamicList<face>> thisCellFaceGroups = cellNonConnectedMultiFaces[cellsAtPoint[k]];
                                bool oneFaceContains = false;
                                for(int l=0;l<thisCellFaceGroupsMap.size();l++)
                                {
                                    bool faceGroupContainsPoint = false;
                                    for(int m=0;m<thisCellFaceGroupsMap[l].size();m++)
                                    {
                                        if(thisCellFaceGroupsMap[l][m].count(zeroPoints[j])!=0)
                                        {
                                            if(oneFaceContains)
                                                FatalErrorInFunction<<"Seperate cut face groups share a point!"<< exit(FatalError);
                                            faceGroupContainsPoint = true;
                                        }
                                    }
                                    if(faceGroupContainsPoint)
                                    {
                                        posClosedPointMapAroundPoint[k] = thisCellFaceGroupsMap[l];
                                        posClosedFacesAroundPoint[k] = thisCellFaceGroups[l];
                                        oneFaceContains = true;
                                    }
                                }
                            }
                            computeClosedFaceFront(zeroPoints[j],posClosedFacesAroundPoint,posClosedPointMapAroundPoint,
                                                    zeroPointsClosedFaces[j],zeroPointsClosedFaceMap[j]);
                        }

                        List<List<scalar>> edgeToOuterAngle(edgeInd.size());
                        for(int j=0;j<edgeInd.size();j++)
                        {
                            edgeToOuterAngle[j] = List<scalar>(2);
                            for(int k=0;k<edgeToOuterAngle.size();k++)
                                edgeToOuterAngle[j][k] = -1;
                        }
                            
                        for(int j=0;j<edgeInd.size();j++)
                        {
                            label firstEdge = j;
                            label secondEdge = (j+1)%edgeInd.size();
                            label cutPoint = newMeshEdges_[edgeInd[firstEdge]].commonVertex(newMeshEdges_[edgeInd[secondEdge]]);
                            if(cutPoint==-1 || cutPoint<nbrOfPrevPoints)
                                FatalErrorInFunction<<"Point between two in a four new edge face. Can not happen!"<< exit(FatalError);
                            
                            edgeToOuterAngle[firstEdge][0]  = computeFaceFrontAngle(zeroPointsClosedFaces[j],zeroPointsClosedFaceMap[j],newMeshEdges_[edgeInd[firstEdge]],cutPoint,newMeshPoints_);
                            edgeToOuterAngle[secondEdge][1] = computeFaceFrontAngle(zeroPointsClosedFaces[j],zeroPointsClosedFaceMap[j],newMeshEdges_[edgeInd[secondEdge]],cutPoint,newMeshPoints_);
                        }
                        // chose edge by angle
                        List<label> connectPntToEdgeWithLowerAngle(4);
                        for(int j=0;j<connectPntToEdgeWithLowerAngle.size();j++)
                            connectPntToEdgeWithLowerAngle[j] = -1;
                        List<scalar> connectPntToMarginWithLowerAngle(4);
                        for(int j=0;j<connectPntToMarginWithLowerAngle.size();j++)
                            connectPntToMarginWithLowerAngle[j] = 0;
                        for(int j=0;j<edgeInd.size();j++)
                        {
                            label firstEdge = j;
                            label secondEdge = (j+1)%edgeInd.size();
                            scalar firstEdgeAngle = edgeToOuterAngle[firstEdge][0];
                            scalar secondEdgeAngle = edgeToOuterAngle[secondEdge][1];
                            
                            if(firstEdgeAngle>secondEdgeAngle)
                            {
                                connectPntToEdgeWithLowerAngle[j] = firstEdge;
                                connectPntToMarginWithLowerAngle[j] = firstEdgeAngle-secondEdgeAngle;
                            }
                            else if(firstEdgeAngle<secondEdgeAngle)
                            {
                                connectPntToEdgeWithLowerAngle[j] = secondEdge;
                                connectPntToMarginWithLowerAngle[j] = secondEdgeAngle-firstEdgeAngle;
                            }
                            else
                                connectPntToEdgeWithLowerAngle[j] = -1;
                        }
                        bool edges02Valid=true;
                        bool edges13Valid=true;
                        bool onlyEqualAngles=true;
                        for(int j=0;j<connectPntToEdgeWithLowerAngle.size();j++)
                        {
                            if(connectPntToEdgeWithLowerAngle[j]==1 || connectPntToEdgeWithLowerAngle[j]==3)
                                edges02Valid=false;
                            if(connectPntToEdgeWithLowerAngle[j]==0 || connectPntToEdgeWithLowerAngle[j]==2)
                                edges13Valid=false;
                            if(connectPntToEdgeWithLowerAngle[j]!=-1)
                                onlyEqualAngles=false;
                        }
                        if(onlyEqualAngles)
                            FatalErrorInFunction<<"All edges have neigboring edges with equal angles! Can not happen."<< exit(FatalError);
                        if(edges02Valid && edges13Valid)
                            FatalErrorInFunction<<"All edges valid ones! Can not happen."<< exit(FatalError);
                        if(edges02Valid)
                        {
                            mustBeOnlyTwoEdges = true;
                            edgeLocalIndsToRemove[0] = 1;
                            edgeLocalIndsToRemove[1] = 3;
                        }
                        else if(edges13Valid)
                        {
                            mustBeOnlyTwoEdges = true;
                            edgeLocalIndsToRemove[0] = 0;
                            edgeLocalIndsToRemove[1] = 2;
                        }
                        else
                        {
                            scalar totalAngle02 = 0;
                            scalar totalAngle13 = 0;
                            for(int j=0;j<edgeInd.size();j++)
                            {
                                if(j==0 || j==2)
                                    totalAngle02 += (edgeToOuterAngle[j][0] + edgeToOuterAngle[j][1]);
                                if(j==1 || j==3)
                                    totalAngle13 += (edgeToOuterAngle[j][0] + edgeToOuterAngle[j][1]);
                            }
                            if(totalAngle02>totalAngle13)
                            {
                                mustBeOnlyTwoEdges = true;
                                edgeLocalIndsToRemove[0] = 1;
                                edgeLocalIndsToRemove[1] = 3;
                            }
                            else if(totalAngle13>totalAngle02)
                            {
                                mustBeOnlyTwoEdges = true;
                                edgeLocalIndsToRemove[0] = 0;
                                edgeLocalIndsToRemove[1] = 2;
                            }
                            else
                            {
                                scalar smallestEdgeToOuterEdgeAngle = 1;
                                label smallestEdgeToOuterEdgeAngleInd = -1;
                                for(int j=0;j<edgeInd.size();j++)
                                {
                                    if(smallestEdgeToOuterEdgeAngle > edgeToOuterAngle[j][0])
                                    {
                                        smallestEdgeToOuterEdgeAngle = edgeToOuterAngle[j][0];
                                        smallestEdgeToOuterEdgeAngleInd = j;
                                    }
                                    if(smallestEdgeToOuterEdgeAngle > edgeToOuterAngle[j][1])
                                    {
                                        smallestEdgeToOuterEdgeAngle = edgeToOuterAngle[j][1];
                                        smallestEdgeToOuterEdgeAngleInd = j;
                                    }
                                }
                                if(smallestEdgeToOuterEdgeAngleInd==-1)
                                    FatalErrorInFunction<<"No edge with smallest angle! Can not happen."<< exit(FatalError);
                                if(smallestEdgeToOuterEdgeAngleInd==0 || smallestEdgeToOuterEdgeAngleInd==2)
                                {
                                    mustBeOnlyTwoEdges = true;
                                    edgeLocalIndsToRemove[0] = 0;
                                    edgeLocalIndsToRemove[1] = 2;
                                }
                                else if(smallestEdgeToOuterEdgeAngleInd==1 || smallestEdgeToOuterEdgeAngleInd==3)
                                {
                                    mustBeOnlyTwoEdges = true;
                                    edgeLocalIndsToRemove[0] = 1;
                                    edgeLocalIndsToRemove[1] = 3;
                                }
                                else
                                    FatalErrorInFunction<<"No edge with smallest angle! Can not happen."<< exit(FatalError);
                            }
                        }
                    }
                }
                if(mustBeAllFourEdges)
                {}
                else if(mustBeOnlyTwoEdges)
                {
                    addedEdgeToDelete.append(edgeInd[edgeLocalIndsToRemove[0]]);
                    addedEdgeToDelete.append(edgeInd[edgeLocalIndsToRemove[1]]);
                }
                else
                {
                    FatalErrorInFunction<< "Face with four cut points but no true cut edges!"<< exit(FatalError);
                    addedEdgeToDelete.append(edgeInd[0]);
                    addedEdgeToDelete.append(edgeInd[1]);
                    addedEdgeToDelete.append(edgeInd[2]);
                    addedEdgeToDelete.append(edgeInd[3]);
                }
            }
            else if(problematicFacePoints[i]==3)
            {
                Info<<"--------------------------Move in 3-----------------------------------"<<endl;
                if(edgeInd.size()!=3)
                    FatalErrorInFunction<< "No three added edges"<< exit(FatalError);

                DynamicList<edge> oldEdges;
                DynamicList<edge> newEdges;
                for(int j=0;j<newEdgeLocalInd.size();j++)
                    newEdges.append(newMeshEdges_[edgeInd[newEdgeLocalInd[j]]]);
                for(int j=0;j<oldEdgeLocalInd.size();j++)
                    oldEdges.append(newMeshEdges_[edgeInd[oldEdgeLocalInd[j]]]);
                
                labelList oldEdgesWithFacesTotalNbr(oldEdgeLocalInd.size(),0);
                labelList oldEdgesWithFacesInFaceCellsNbr(oldEdgeLocalInd.size(),0);
                // Counting as one for each cell with face
                for(int j=0;j<oldEdgeLocalInd.size();j++)
                {
                    //collect connected faces at old Edges in all connected cells
                    const labelList& cellsOfEdge = edgeToCells[oldEdgeLocalInd[j]];
                    label nbrRelCells = cellsOfEdge.size();
                    List<DynamicList<face>> posClosedFacesAroundEdge(nbrRelCells);
                    List<DynamicList<std::unordered_set<label>>> posClosedPointMapAroundEdge(nbrRelCells);
                    for(int k=0;k<cellsOfEdge.size();k++)
                    {
                        DynamicList<DynamicList<std::unordered_set<label>>> thisCellFaceGroupsMap = cellNonConnectedMultiEdgeMap[cellsOfEdge[k]];
                        DynamicList<DynamicList<face>> thisCellFaceGroups = cellNonConnectedMultiFaces[cellsOfEdge[k]];
                        bool oneFaceContains = false;
                        for(int l=0;l<thisCellFaceGroupsMap.size();l++)
                        {
                            bool faceGroupContainsPoint = false;
                            for(int m=0;m<thisCellFaceGroupsMap[l].size();m++)
                            {
                                if(thisCellFaceGroupsMap[l][m].count(oldEdgeLocalInd[j])!=0)
                                {
                                    if(oneFaceContains)
                                        FatalErrorInFunction<<"Seperate cut face groups share a point!"<< exit(FatalError);
                                    faceGroupContainsPoint = true;
                                }
                            }
                            if(faceGroupContainsPoint)
                            {
                                posClosedPointMapAroundEdge[k] = thisCellFaceGroupsMap[l];
                                posClosedFacesAroundEdge[k] = thisCellFaceGroups[l];
                                oneFaceContains = true;
                            }
                        }
                    }
                    label faceCounter = 0;
                    for(int k=0;k<nbrRelCells;k++)
                    {
                        if(posClosedFacesAroundEdge[k].size()!=0 && posClosedPointMapAroundEdge[k].size()!=0)
                        {
                            oldEdgesWithFacesTotalNbr[j]++;
                            bool cellIsOwnerOrNeighbor = false;
                            for(int l=0;l<faceCells.size();l++)
                            {
                                if(faceCells[l]==cellsOfEdge[j])
                                    cellIsOwnerOrNeighbor = true;
                            }
                            if(cellIsOwnerOrNeighbor)
                                oldEdgesWithFacesInFaceCellsNbr[j]++;
                        }
                        else if(posClosedFacesAroundEdge[k].size()==0 && posClosedPointMapAroundEdge[k].size()==0)
                        {}
                        else
                            FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                    }
                }
                
                label centralZeroPoint;
                labelList otherZeroPoints(2);
                if(problematicFaceNewPoints[i]==0)
                {
                    if(oldEdge!=2 || newEdge!=1)
                        FatalErrorInFunction<< "Invalid old and new Edges number!"<< exit(FatalError);
                    centralZeroPoint = oldEdges[0].commonVertex(oldEdges[1]);
                    otherZeroPoints[0]=(oldEdges[0].otherVertex(centralZeroPoint));
                    otherZeroPoints[1]=(oldEdges[1].otherVertex(centralZeroPoint));
                }
                else if(problematicFaceNewPoints[i]==1)
                {
                    if(oldEdge!=1 || newEdge!=2)
                        FatalErrorInFunction<< "Invalid old and new Edges number!"<< exit(FatalError);
                    centralZeroPoint = newEdges[0].commonVertex(newEdges[1]);
                    otherZeroPoints[0]=(newEdges[0].otherVertex(centralZeroPoint));
                    otherZeroPoints[1]=(newEdges[1].otherVertex(centralZeroPoint));
                }
                else if(problematicFaceNewPoints[i]==2)
                {
                    if(oldEdge!=0 || newEdge!=3)
                        FatalErrorInFunction<< "Invalid old and new Edges number!"<< exit(FatalError);
                    bool oldZeroPointFound = false;
                    for(int j=0;j<newEdges.size();j++)
                    {
                        label nextInd = (j+1)%newEdges.size();
                        label connectVertex = newEdges[j].commonVertex(newEdges[nextInd]);
                        if(connectVertex < basisPoints.size())
                        {
                            if(oldZeroPointFound)
                                FatalErrorInFunction<<"Cut face seems to have old vertices. Can not happen!"<< exit(FatalError);
                            centralZeroPoint = connectVertex;
                            otherZeroPoints[0]=(newEdges[j].otherVertex(connectVertex));
                            otherZeroPoints[1]=(newEdges[nextInd].otherVertex(connectVertex));
                            oldZeroPointFound = true;
                        }
                    }
                    if(oldZeroPointFound)
                        FatalErrorInFunction<<"Cut face seems to have old vertices. Can not happen!"<< exit(FatalError);
                }
                else
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                if(centralZeroPoint==-1)
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                
                std::unordered_map<label,label> otherZeroPntLocalIndToEdgeLocalInd;
                std::unordered_map<label,label> edgeLocalIndTootherZeroPntLocalInd;
                label centerEdgeLocalInd=-1;
                for(int j=0;j<edgeInd.size();j++)
                {
                    edge nEdge = newMeshEdges_[edgeInd[j]];
                    if(nEdge.otherVertex(otherZeroPoints[0])==otherZeroPoints[1])
                    {
                        if(centerEdgeLocalInd==-1)
                            centerEdgeLocalInd=j;
                        else
                            FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                    }
                    else if(nEdge.otherVertex(otherZeroPoints[0])==centralZeroPoint)
                    {
                        otherZeroPntLocalIndToEdgeLocalInd[0]=j;
                        edgeLocalIndTootherZeroPntLocalInd[j]=0;
                    }
                    else if(nEdge.otherVertex(otherZeroPoints[1])==centralZeroPoint)
                    {
                        otherZeroPntLocalIndToEdgeLocalInd[1]=j;
                        edgeLocalIndTootherZeroPntLocalInd[j]=1;
                    }
                    else
                        FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                }                
                DynamicList<label> zeroPoints({otherZeroPoints[0],centralZeroPoint,otherZeroPoints[1]});
                
                List<DynamicList<DynamicList<face>>> zeroPointsClosedFaces(zeroPoints.size());
                List<DynamicList<DynamicList<std::unordered_set<label>>>> zeroPointsClosedFaceMap(zeroPoints.size());
                
                for(int j=0;j<zeroPoints.size();j++)
                {
                    //collect all faces in neighboring cells connected to the point in posClosedFacesAroundPoint
                    const labelList& cellsAtPoint = pointToCells[zeroPoints[j]];
                    if(zeroPoints[j]>=nbrOfPrevPoints)
                    {
                        Info<<"zeroPoints[j]:"<<zeroPoints[j]<<endl;
                        Info<<"zeroPoints:"<<zeroPoints<<endl;
                        Info<<"cellsAtPoint:"<<cellsAtPoint<<endl;
                        FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);
                    }
                    label nbrRelCells = cellsAtPoint.size();
                    List<DynamicList<face>> posClosedFacesAroundPoint(nbrRelCells);
                    List<DynamicList<std::unordered_set<label>>> posClosedPointMapAroundPoint(nbrRelCells);
                    for(int k=0;k<cellsAtPoint.size();k++)
                    {
                        DynamicList<DynamicList<std::unordered_set<label>>> thisCellFaceGroupsMap = cellNonConnectedMultiPointMap[cellsAtPoint[k]];
                        DynamicList<DynamicList<face>> thisCellFaceGroups = cellNonConnectedMultiFaces[cellsAtPoint[k]];
                        bool oneFaceContains = false;
                        for(int l=0;l<thisCellFaceGroupsMap.size();l++)
                        {
                            bool faceGroupContainsPoint = false;
                            for(int m=0;m<thisCellFaceGroupsMap[l].size();m++)
                            {
                                if(thisCellFaceGroupsMap[l][m].count(zeroPoints[j])!=0)
                                {
                                    if(oneFaceContains)
                                        FatalErrorInFunction<<"Seperate cut face groups share a point!"<< exit(FatalError);
                                    faceGroupContainsPoint = true;
                                }
                            }
                            if(faceGroupContainsPoint)
                            {
                                posClosedPointMapAroundPoint[k] = thisCellFaceGroupsMap[l];
                                posClosedFacesAroundPoint[k] = thisCellFaceGroups[l];
                                oneFaceContains = true;
                            }
                        }
                        Info<<"k:"<<k<<endl;
                        Info<<"zeroPoints:"<<zeroPoints<<endl;
                        Info<<"thisCellFaceGroups:"<<thisCellFaceGroups<<endl;
                    }
                    Info<<"----------------------In here-------------------____"<<endl;
                    Info<<"j:"<<j<<"/"<<zeroPoints.size()<<endl;
                    Info<<"i:"<<i<<endl;
                    Info<<"cellsAtPoint:"<<cellsAtPoint<<endl;
                    Info<<"posClosedFacesAroundPoint:"<<posClosedFacesAroundPoint<<endl;
                    for(int a=0;a<posClosedFacesAroundPoint.size();a++)
                    {
                        Info<<"-------------------------------"<<endl;
                        Info<<"posClosedFacesAroundPoint["<<a<<"]:"<<posClosedFacesAroundPoint[a]<<endl;
                        for(int b=0;b<posClosedFacesAroundPoint[a].size();b++)
                        {
                            for(int c=0;c<posClosedFacesAroundPoint[a][b].size();c++)
                            {
                                if(posClosedFacesAroundPoint[a][b][c]>=nbrOfPrevPoints)
                                {
                                    if(pointToEgde_[posClosedFacesAroundPoint[a][b][c]]==-1)
                                        FatalErrorInFunction<<"Wrong!"<< exit(FatalError);
                                    edge pointEdge = newMeshEdges_[pointToEgde_[posClosedFacesAroundPoint[a][b][c]]];
                                    Info<<"  "<<posClosedFacesAroundPoint[a][b][c]<<":"<<pointEdge;
                                }
                                else
                                {
                                    if(pointToEgde_[posClosedFacesAroundPoint[a][b][c]]!=-1)
                                        FatalErrorInFunction<<"Wrong!"<< exit(FatalError);
                                    Info<<"  "<<posClosedFacesAroundPoint[a][b][c]<<":"<<"-1";
                                }
                            }
                            Info<<endl;
                        }
                    }
                    Info<<"zeroPoints[j]:"<<zeroPoints[j]<<endl;
                    Info<<"nbrOfPrevPoints:"<<nbrOfPrevPoints<<endl;
                    Info<<"problematicFaceNewPoints[i]:"<<problematicFaceNewPoints[i]<<endl;
                    Info<<"basisFaces[i]:"<<basisFaces[i]<<endl;
                    for(int a=0;a<cellsAtPoint.size();a++)
                    {
                        cell oneCell = meshCells[cellsAtPoint[a]];
                        Info<<"Cell:"<<cellsAtPoint[a]<<"  "<<oneCell.labels(basisFaces)<<endl<<"\t\t"<<oneCell.points(basisFaces,basisPoints)<<endl;
                    }

                    computeClosedFaceFront(zeroPoints[j],posClosedFacesAroundPoint,posClosedPointMapAroundPoint,
                                           zeroPointsClosedFaces[j],zeroPointsClosedFaceMap[j]);
                    Info<<"zeroPointsClosedFaces[j]:"<<zeroPointsClosedFaces[j]<<endl;
                    if(i==585222 && j==1)
                        FatalErrorInFunction<<"Temp Stop!"<< exit(FatalError);
                }
                
                List<List<bool>> vertexConnectedToVertex0(3);
                List<List<bool>> vertexConnectedToVertexCenter(3);
                List<List<bool>> vertexConnectedToVertex1(3);                    
                for(int j=0;j<3;j++)
                {
                    vertexConnectedToVertex0[j] = pointInFaceFront(zeroPointsClosedFaces[j],zeroPointsClosedFaceMap[j],otherZeroPoints[0]);
                    vertexConnectedToVertexCenter[j] = pointInFaceFront(zeroPointsClosedFaces[j],zeroPointsClosedFaceMap[j],centralZeroPoint);
                    vertexConnectedToVertex1[j] = pointInFaceFront(zeroPointsClosedFaces[j],zeroPointsClosedFaceMap[j],otherZeroPoints[1]);
                }
                for(int j=0;j<3;j++)
                {
                    if((vertexConnectedToVertex0[j].size()!=vertexConnectedToVertexCenter[j].size()) ||
                       (vertexConnectedToVertex1[j].size()!=vertexConnectedToVertexCenter[j].size()))
                        FatalErrorInFunction<<"Inequal number of face front on vertex "<<j<<" !"<< exit(FatalError);
                }
                List<bool> centerVertexHasFaceConnectionOnlyToOtherVertexX(2,false);
                List<bool> otherVertexXHasFaceConnectionOnlyToCenter(2,false);
                List<bool> otherVertexXHasFaceConnectionOnlyWithOppositeVertex(2,false);

                List<bool> otherVertexXHasNoFreeFaceConnection(2,true);
                List<bool> otherVertexXHasOnlyFreeFaceConnection(2,true);

                List<bool> otherVertexXHasNoFaceConnectionFullyConnected(2,true);
                List<bool> otherVertexXHasOnlyFaceConnectionFullyConnected(2,true);
                
                bool centerVertexHasNoFreeFaceConnection = true;
                bool centerVertexHasOnlyFreeFaceConnection = true;
                
                for(int j=0;j<vertexConnectedToVertex0[0].size();j++)
                {
                    if(!vertexConnectedToVertex0[0][j])
                        FatalErrorInFunction<< "Face front must have its own central Point!"<< exit(FatalError);                        
                    if(vertexConnectedToVertexCenter[0][j] && !vertexConnectedToVertex1[0][j])
                        otherVertexXHasFaceConnectionOnlyToCenter[0] = true;
                    if(!vertexConnectedToVertexCenter[0][j] && vertexConnectedToVertex1[0][j])
                        otherVertexXHasFaceConnectionOnlyWithOppositeVertex[0] = true;
                    
                    if(!vertexConnectedToVertexCenter[0][j] && !vertexConnectedToVertex1[0][j])
                        otherVertexXHasNoFreeFaceConnection[0] = false;
                    if(vertexConnectedToVertexCenter[0][j] || vertexConnectedToVertex1[0][j])
                        otherVertexXHasOnlyFreeFaceConnection[0] = false;
                    
                    if(vertexConnectedToVertexCenter[0][j] && vertexConnectedToVertex1[0][j])
                        otherVertexXHasNoFaceConnectionFullyConnected[0]=false;
                    if(!vertexConnectedToVertexCenter[0][j] || !vertexConnectedToVertex1[0][j])
                        otherVertexXHasOnlyFaceConnectionFullyConnected[0]=false;
                    
                }
                for(int j=0;j<vertexConnectedToVertex1[2].size();j++)
                {
                    if(!vertexConnectedToVertex1[2][j])
                        FatalErrorInFunction<< "Face front must have its own central Point!"<< exit(FatalError);                        
                    if(vertexConnectedToVertexCenter[2][j] && !vertexConnectedToVertex0[2][j])
                        otherVertexXHasFaceConnectionOnlyToCenter[1] = true;
                    if(!vertexConnectedToVertexCenter[2][j] && vertexConnectedToVertex0[2][j])
                        otherVertexXHasFaceConnectionOnlyWithOppositeVertex[1] = true;
                    
                    if(!vertexConnectedToVertexCenter[2][j] && !vertexConnectedToVertex0[2][j])
                        otherVertexXHasNoFreeFaceConnection[1] = false;
                    if(vertexConnectedToVertexCenter[2][j] || vertexConnectedToVertex0[2][j])
                        otherVertexXHasOnlyFreeFaceConnection[1] = false;
                    
                    if(vertexConnectedToVertexCenter[2][j] && vertexConnectedToVertex0[2][j])
                        otherVertexXHasNoFaceConnectionFullyConnected[1]=false;
                    if(!vertexConnectedToVertexCenter[2][j] || !vertexConnectedToVertex0[2][j])
                        otherVertexXHasOnlyFaceConnectionFullyConnected[1]=false;
                }
                for(int j=0;j<vertexConnectedToVertexCenter[1].size();j++)
                {        
                    if(!vertexConnectedToVertexCenter[1][j])
                        FatalErrorInFunction<< "Face front must have its own central Point!"<< exit(FatalError);
                    if(vertexConnectedToVertex0[1][j] && !vertexConnectedToVertex1[1][j])
                        centerVertexHasFaceConnectionOnlyToOtherVertexX[0] = true;
                    if(!vertexConnectedToVertex0[1][j] && vertexConnectedToVertex1[1][j])
                        centerVertexHasFaceConnectionOnlyToOtherVertexX[1] = true;
                    if(!vertexConnectedToVertex0[1][j] && !vertexConnectedToVertex1[1][j])
                        centerVertexHasNoFreeFaceConnection = false;
                    if(vertexConnectedToVertex0[1][j] || vertexConnectedToVertex1[1][j])
                        centerVertexHasOnlyFreeFaceConnection = false;
                }
                
                if(otherVertexXHasFaceConnectionOnlyToCenter[0] && otherVertexXHasOnlyFreeFaceConnection[0])
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                if(otherVertexXHasFaceConnectionOnlyToCenter[1] && otherVertexXHasOnlyFreeFaceConnection[1])
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                if(otherVertexXHasOnlyFaceConnectionFullyConnected[0] && otherVertexXHasOnlyFreeFaceConnection[0])
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                if(otherVertexXHasOnlyFaceConnectionFullyConnected[1] && otherVertexXHasOnlyFreeFaceConnection[1])
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                
                if((otherVertexXHasNoFreeFaceConnection[0] && otherVertexXHasOnlyFreeFaceConnection[0]) ||
                   (otherVertexXHasNoFreeFaceConnection[1] && otherVertexXHasOnlyFreeFaceConnection[1]))
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                if((otherVertexXHasNoFaceConnectionFullyConnected[0] && otherVertexXHasOnlyFaceConnectionFullyConnected[0]) ||
                   (otherVertexXHasNoFaceConnectionFullyConnected[1] && otherVertexXHasOnlyFaceConnectionFullyConnected[1]))
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                if(centerVertexHasNoFreeFaceConnection && centerVertexHasOnlyFreeFaceConnection)
                {
                    Info<<"i:"<<i<<endl;
                    Info<<"basisFaces[i]:"<<basisFaces[i]<<endl;
                    Info<<"problematicFaceNewPoints[i]:"<<problematicFaceNewPoints[i]<<endl;
                    Info<<"nbrOfPrevPoints:"<<nbrOfPrevPoints<<endl;
                    Info<<"zeroPoints:"<<zeroPoints<<endl;
                    Info<<"zeroPointsClosedFaces:"<<zeroPointsClosedFaces<<endl;
                    Info<<"vertexConnectedToVertex0:"<<vertexConnectedToVertex0<<endl;
                    Info<<"vertexConnectedToVertexCenter:"<<vertexConnectedToVertexCenter<<endl;
                    Info<<"vertexConnectedToVertex1:"<<vertexConnectedToVertex1<<endl;
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                }
                
                if(problematicFaceNewPoints[i]==0)
                {
                    if(newEdgeLocalInd.size()!=1 || oldEdgeLocalInd.size()!=2)
                        FatalErrorInFunction<<"Incorrect number of cut edges!"<< exit(FatalError);
                        
                    label nonZeroPntLocalVerticeInd = -1;
                    for(int j=0;j<basisFaces[i].size();j++)
                    {
                        if(pointDist[basisFaces[i][j]]!=0)
                        {
                            if(nonZeroPntLocalVerticeInd!=-1)
                                nonZeroPntLocalVerticeInd = j;
                            else
                                FatalErrorInFunction<<"Double non zero point. Can not happen!"<< exit(FatalError);
                        }
                    }
                    if(nonZeroPntLocalVerticeInd==-1)
                        FatalErrorInFunction<<"No non zero point. Can not happen!"<< exit(FatalError);
                    
                    bool mustBeThreeEdges = false;
                    // three edge face is when the new edge has a connected face to one cell but not to the other
                    if(faceCells.size()==2)
                    {
                        if(oldEdgesWithFacesTotalNbr[0]>0 && oldEdgesWithFacesTotalNbr[1]>0 && connectedToCellPerEdge[newEdgeLocalInd[0]][0])
                        //Must
                        {
                            if(!connectedToCellPerEdge[newEdgeLocalInd[0]][1])
                            //Would be invalid face
                            {
                                if(verticalEdgesAreCut[nonZeroPntLocalVerticeInd][0])                             
                                //Improbable but valid face
                                {
                                    mustBeThreeEdges = true;
                                }
                            }
                        }
                        if(oldEdgesWithFacesTotalNbr[0]>0 && oldEdgesWithFacesTotalNbr[1]>0 && connectedToCellPerEdge[newEdgeLocalInd[0]][0])
                        //Must
                        {
                            if(!connectedToCellPerEdge[newEdgeLocalInd[0]][1])
                            //Would be invalid face
                            {
                                if(verticalEdgesAreCut[nonZeroPntLocalVerticeInd][0])                             
                                //Improbable but valid face
                                {
                                    mustBeThreeEdges = true;
                                }
                            }
                        }
                    }
                    else if(faceCells.size()==1)
                    {                    
                        if(oldEdgesWithFacesTotalNbr[0]>0 && oldEdgesWithFacesTotalNbr[1]>0 && connectedToCellPerEdge[newEdgeLocalInd[0]][0])
                        //Must
                        {
                            if(verticalEdgesAreCut[nonZeroPntLocalVerticeInd][0])                             
                            //Improbable but valid face
                            {
                                mustBeThreeEdges = true;
                            }
                        }
                    }
                    else
                        FatalErrorInFunction<< "Face must have either one or two neighbor cells!"<< exit(FatalError);
                    
                    bool mustBeOnlyNewEdge = false;
                    /* only new edge face is when there are faces connected to new edge in both cells and 
                    * there are four faces connected to centralZeroPoint that are not connected to the face at all
                    * other than by this point
                    */
                    if(!mustBeThreeEdges)
                    {
                        if(faceCells.size()==2)
                        {
                            if(connectedToCellPerEdge[newEdgeLocalInd[0]][0] && connectedToCellPerEdge[newEdgeLocalInd[0]][1])
                            {                            
                                if(!centerVertexHasNoFreeFaceConnection)
                                    mustBeOnlyNewEdge = true;
                            }
                        }
                        else if(faceCells.size()==1)
                        {
                            if(connectedToCellPerEdge[newEdgeLocalInd[0]][0])
                            {
                                if(!centerVertexHasNoFreeFaceConnection)
                                    mustBeOnlyNewEdge = true;
                            }
                        }
                        else
                            FatalErrorInFunction<< "Face must have either one or two neighbor cells!"<< exit(FatalError);
                    }
                    
                    bool mustBeOnlyOneOldEdge = false;
                    label otherOldEdgeToRemoveLocalInd = -1;
                    /* only one old edge face is when there are faces connected to one old edge in four cells and 
                    * there are four faces connected to one outer old zero point that are not connected to the one
                    * single point
                    */
                    if(!mustBeThreeEdges && !mustBeOnlyNewEdge)
                    {
                        if(!otherVertexXHasNoFreeFaceConnection[0] && centerVertexHasFaceConnectionOnlyToOtherVertexX[1] && otherVertexXHasFaceConnectionOnlyToCenter[1])
                        {
                            label oldEdgeListIndOfOtherVertex = edgeLocalIndToOldEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[1]];
                            if((faceCells.size()==2 && oldEdgesWithFacesTotalNbr[oldEdgeListIndOfOtherVertex]>=2)||(faceCells.size()==1 && oldEdgesWithFacesTotalNbr[oldEdgeListIndOfOtherVertex]>=1))
                            {
                                otherOldEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[0];
                                mustBeOnlyOneOldEdge = true;
                            }
                        }
                        else if(!otherVertexXHasNoFreeFaceConnection[1] && centerVertexHasFaceConnectionOnlyToOtherVertexX[0] && otherVertexXHasFaceConnectionOnlyToCenter[0])
                        {
                            label oldEdgeListIndOfOtherVertex = edgeLocalIndToOldEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[1]];
                            if((faceCells.size()==2 && oldEdgesWithFacesTotalNbr[oldEdgeListIndOfOtherVertex]>=2)||(faceCells.size()==1 && oldEdgesWithFacesTotalNbr[oldEdgeListIndOfOtherVertex]>=1))
                            {
                                otherOldEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[1];
                                mustBeOnlyOneOldEdge = true;
                            }
                        }
                    }
                    
                    bool mustBeOTheNewAndOneOldEdge = false;
                    label singleOldEdgeToRemoveLocalInd = -1;
                    /* only one old edge and the new edge is when there is a face connection to one otherVertex and
                    * this face connection does only connect to the other "otherVertex".  
                    * in addition the old edge needs two faces
                    */
                    if(!mustBeThreeEdges && !mustBeOnlyNewEdge && !mustBeOnlyOneOldEdge)
                    {
                        if(centerVertexHasFaceConnectionOnlyToOtherVertexX[0] && otherVertexXHasNoFreeFaceConnection[0] && otherVertexXHasFaceConnectionOnlyWithOppositeVertex[1])
                        {
                            label oldEdgeListIndOfOtherVertex = edgeLocalIndToOldEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[0]];
                            if((faceCells.size()==2 && oldEdgesWithFacesTotalNbr[oldEdgeListIndOfOtherVertex]>=2 && connectedToCellPerEdge[newEdgeLocalInd[0]][0] && connectedToCellPerEdge[newEdgeLocalInd[0]][1])||
                            (faceCells.size()==1 && oldEdgesWithFacesTotalNbr[oldEdgeListIndOfOtherVertex]>=1 && connectedToCellPerEdge[newEdgeLocalInd[0]][0]))
                            {
                                singleOldEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[1];
                                mustBeOnlyOneOldEdge = true;
                            }
                        } 
                        else if(centerVertexHasFaceConnectionOnlyToOtherVertexX[1] && otherVertexXHasNoFreeFaceConnection[1] && otherVertexXHasFaceConnectionOnlyWithOppositeVertex[0])
                        {
                            label oldEdgeListIndOfOtherVertex = edgeLocalIndToOldEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[1]];
                            if((faceCells.size()==2 && oldEdgesWithFacesTotalNbr[oldEdgeListIndOfOtherVertex]>=2 && connectedToCellPerEdge[newEdgeLocalInd[0]][0] && connectedToCellPerEdge[newEdgeLocalInd[0]][1])||
                            (faceCells.size()==1 && oldEdgesWithFacesTotalNbr[oldEdgeListIndOfOtherVertex]>=1 && connectedToCellPerEdge[newEdgeLocalInd[0]][0]))
                            {
                                singleOldEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[0];
                                mustBeOnlyOneOldEdge = true;
                            }
                        }
                    }
                    
                    bool mustBeOnlyBothOldEdges = false;
                    /* only the old edges are correct if no other situation is correct
                    * and there are two faces for both old edges
                    */
                    if(!mustBeThreeEdges && !mustBeOnlyNewEdge && !mustBeOnlyOneOldEdge && !mustBeOTheNewAndOneOldEdge)
                    {
                        if(otherVertexXHasFaceConnectionOnlyToCenter[0] && otherVertexXHasFaceConnectionOnlyToCenter[1] && centerVertexHasNoFreeFaceConnection)
                        {
                        if((faceCells.size()==2 && oldEdgesWithFacesTotalNbr[0]>=2 && oldEdgesWithFacesTotalNbr[1]>=2)||
                            (faceCells.size()==1 && oldEdgesWithFacesTotalNbr[0]>=1 && oldEdgesWithFacesTotalNbr[1]>=1))
                            {
                                mustBeOnlyBothOldEdges = true;
                            }
                        }
                    }
                    
                    if(mustBeThreeEdges)
                    {}
                    else if(mustBeOnlyNewEdge)
                    {
                        addedEdgeToDelete.append(edgeInd[oldEdgeLocalInd[0]]);
                        addedEdgeToDelete.append(edgeInd[oldEdgeLocalInd[1]]);
                    }
                    else if(mustBeOnlyOneOldEdge)
                    {
                        addedEdgeToDelete.append(edgeInd[otherOldEdgeToRemoveLocalInd]);
                        addedEdgeToDelete.append(edgeInd[newEdgeLocalInd[0]]);
                    }
                    else if(mustBeOTheNewAndOneOldEdge)
                    {
                        addedEdgeToDelete.append(edgeInd[singleOldEdgeToRemoveLocalInd]);
                    }
                    else if(mustBeOnlyBothOldEdges)
                    {
                        addedEdgeToDelete.append(edgeInd[newEdgeLocalInd[0]]);
                    }
                    else
                    {
                        if(faceCells.size()==2)
                            FatalErrorInFunction<< "Completely false face can not happen for inner face"<< exit(FatalError);
                        
                        addedEdgeToDelete.append(edgeInd[oldEdgeLocalInd[0]]);
                        addedEdgeToDelete.append(edgeInd[oldEdgeLocalInd[1]]);
                        addedEdgeToDelete.append(edgeInd[newEdgeLocalInd[0]]);
                    }
                }
                else if(problematicFaceNewPoints[i]==1)
                {
                    if(newEdgeLocalInd.size()!=2 || oldEdgeLocalInd.size()!=1)
                        FatalErrorInFunction<<"Incorrect number of cut edges!"<< exit(FatalError);
                    
                    DynamicList<label> nonZeroPntLocalVerticeInd;
                    for(int j=0;j<basisFaces[i].size();j++)
                    {
                        if(pointDist[basisFaces[i][j]]!=0)
                        {
                            if((nonZeroPntLocalVerticeInd.size()==0)||(nonZeroPntLocalVerticeInd.size()==1))
                                nonZeroPntLocalVerticeInd.append(j);
                            else
                                FatalErrorInFunction<<"Other than two zero point. Can not happen!"<< exit(FatalError);
                        }
                    }
                    if(nonZeroPntLocalVerticeInd.size()!=2)
                        FatalErrorInFunction<<"Other than two zero point. Can not happen!"<< exit(FatalError);
                    
                    bool mustBeThreeEdges = false;
                    // three edge face is when all edges have exactly  a connected face to one cell but not to the other                
                    if(faceCells.size()==2)
                    {
                        bool validVerticalEdgeCuts = false;
                        if((connectedToCellPerEdge[newEdgeLocalInd[0]][0] && !connectedToCellPerEdge[newEdgeLocalInd[0]][1]) && 
                           (connectedToCellPerEdge[newEdgeLocalInd[1]][1] && !connectedToCellPerEdge[newEdgeLocalInd[1]][0]))
                        // Must be
                        {
                            if(verticalEdgesAreCut[cutEdgeLocalIndToNonZeroPntLocalInd[0]][0] && verticalEdgesAreCut[cutEdgeLocalIndToNonZeroPntLocalInd[1]][1])
                            // Valid if not true
                                validVerticalEdgeCuts=true;
                        }
                        if((!connectedToCellPerEdge[newEdgeLocalInd[0]][0] && connectedToCellPerEdge[newEdgeLocalInd[0]][1]) && 
                           (!connectedToCellPerEdge[newEdgeLocalInd[1]][1] && connectedToCellPerEdge[newEdgeLocalInd[1]][0]))
                        // Must be
                        {
                            if(verticalEdgesAreCut[cutEdgeLocalIndToNonZeroPntLocalInd[0]][1] && verticalEdgesAreCut[cutEdgeLocalIndToNonZeroPntLocalInd[1]][0])
                            // Valid if not true
                                validVerticalEdgeCuts=true;
                        }                            
                        if(oldEdgesWithFacesTotalNbr[0]>=1 && validVerticalEdgeCuts)
                        // Must be
                            mustBeThreeEdges = true;
                    }
                    else if(faceCells.size()==1)
                    {                    
                        bool validVerticalEdgeCuts = false;
                        if((connectedToCellPerEdge[newEdgeLocalInd[0]][0] && !connectedToCellPerEdge[newEdgeLocalInd[1]][0]))
                        // Must be
                        {
                            if(verticalEdgesAreCut[cutEdgeLocalIndToNonZeroPntLocalInd[0]][0])
                            // Valid if not true
                                validVerticalEdgeCuts=true;
                        }
                        if((!connectedToCellPerEdge[newEdgeLocalInd[0]][0] && connectedToCellPerEdge[newEdgeLocalInd[1]][0]))
                        // Must be
                        {
                            if(verticalEdgesAreCut[cutEdgeLocalIndToNonZeroPntLocalInd[1]][0])
                            // Valid if not true
                                validVerticalEdgeCuts=true;
                        }
                        if(validVerticalEdgeCuts)
                            mustBeThreeEdges = true;
                    }
                    else
                        FatalErrorInFunction<< "Face must have either one or two neighbor cells!"<< exit(FatalError);
                    
                    bool mustBeOnlyOneNewEdge = false;
                    label newEdgeToRemoveLocalInd = -1;
                    /* only new edge face is if one Vertex has free face connection and the other and the center
                     * vertex have a connection
                    */
                    if(!mustBeThreeEdges)
                    {
                        if(faceCells.size()==2)
                        {
                            if(!otherVertexXHasNoFreeFaceConnection[0] && otherVertexXHasFaceConnectionOnlyToCenter[1] && centerVertexHasFaceConnectionOnlyToOtherVertexX[1])
                            // Must be
                            {
                                label newEdgeListIndOfOtherVertex = egdeLocalIndToNewEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[1]];
                                if(connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][0] && connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][1])
                                // Must be
                                {
                                    mustBeOnlyOneNewEdge = true;
                                    newEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[0];
                                }
                            }
                            else if(!otherVertexXHasNoFreeFaceConnection[1] && otherVertexXHasFaceConnectionOnlyToCenter[0] && centerVertexHasFaceConnectionOnlyToOtherVertexX[0])
                            // Must be
                            {
                                label newEdgeListIndOfOtherVertex = egdeLocalIndToNewEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[0]];
                                if(connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][0] && connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][1])
                                // Must be
                                {
                                    mustBeOnlyOneNewEdge = true;
                                    newEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[1];
                                }
                            }
                        }
                        else if(faceCells.size()==1)
                        {
                            if(!otherVertexXHasNoFreeFaceConnection[0] && otherVertexXHasFaceConnectionOnlyToCenter[1] && centerVertexHasFaceConnectionOnlyToOtherVertexX[1])
                            // Must be
                            {
                                label newEdgeListIndOfOtherVertex = egdeLocalIndToNewEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[1]];
                                if(connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][0])
                                // Must be
                                {
                                    mustBeOnlyOneNewEdge = true;
                                    newEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[0];
                                }
                            }
                            else if(!otherVertexXHasNoFreeFaceConnection[1] && otherVertexXHasFaceConnectionOnlyToCenter[0] && centerVertexHasFaceConnectionOnlyToOtherVertexX[0])
                            // Must be
                            {
                                label newEdgeListIndOfOtherVertex = egdeLocalIndToNewEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[0]];
                                if(connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][0])
                                // Must be
                                {
                                    mustBeOnlyOneNewEdge = true;
                                    newEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[1];
                                }
                            }
                        }
                        else
                            FatalErrorInFunction<< "Face must have either one or two neighbor cells!"<< exit(FatalError);
                    }
                    
                    bool mustBeOneNewAndOneOldEdge = false;
                    label newEdgeToSingleRemoveLocalInd = -1;
                    /* only one new and one old edge if one vertex has no free face connection and other vertex 
                    */
                    if(!mustBeThreeEdges && !mustBeOnlyOneNewEdge)
                    {
                        if(faceCells.size()==2)
                        {
                            if(!otherVertexXHasNoFaceConnectionFullyConnected[0] && centerVertexHasFaceConnectionOnlyToOtherVertexX[1] && otherVertexXHasFaceConnectionOnlyWithOppositeVertex[0])
                            // Must be
                            {
                                label newEdgeListIndOfOtherVertex = egdeLocalIndToNewEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[0]];
                                if(connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][0] && connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][1] && oldEdgesWithFacesTotalNbr[0]>=2)
                                // Must be
                                {
                                    if(!otherVertexXHasNoFreeFaceConnection[0])
                                    // Probable but faces would be invalid
                                    {
                                        mustBeOnlyOneNewEdge = true;
                                        newEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[1];
                                    }
                                }
                            }
                            else if(!otherVertexXHasNoFaceConnectionFullyConnected[1] && centerVertexHasFaceConnectionOnlyToOtherVertexX[0] && otherVertexXHasFaceConnectionOnlyWithOppositeVertex[1])
                            // Must be
                            {
                                label newEdgeListIndOfOtherVertex = egdeLocalIndToNewEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[1]];
                                if(connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][0] && connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][1] && oldEdgesWithFacesTotalNbr[0]>=2)
                                // Must be
                                {
                                    if(!otherVertexXHasNoFreeFaceConnection[1])
                                    // Probable but faces would be invalid
                                    {
                                        mustBeOnlyOneNewEdge = true;
                                        newEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[0];
                                    }
                                }
                            }
                        }
                        else if(faceCells.size()==1)
                        {
                            if(!otherVertexXHasNoFaceConnectionFullyConnected[0] && centerVertexHasFaceConnectionOnlyToOtherVertexX[1] && otherVertexXHasFaceConnectionOnlyWithOppositeVertex[0])
                            // Must be
                            {
                                label newEdgeListIndOfOtherVertex = egdeLocalIndToNewEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[0]];
                                if(connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][0] && oldEdgesWithFacesTotalNbr[0]>=1)
                                // Must be
                                {
                                    if(!otherVertexXHasNoFreeFaceConnection[0])
                                    // Probable but faces would be invalid
                                    {
                                        mustBeOnlyOneNewEdge = true;
                                        newEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[1];
                                    }
                                }
                            }
                            else if(!otherVertexXHasNoFaceConnectionFullyConnected[1] && centerVertexHasFaceConnectionOnlyToOtherVertexX[0] && otherVertexXHasFaceConnectionOnlyWithOppositeVertex[1])
                            // Must be
                            {
                                label newEdgeListIndOfOtherVertex = egdeLocalIndToNewEdgeListInd[otherZeroPntLocalIndToEdgeLocalInd[1]];
                                if(connectedToCellPerEdge[newEdgeLocalInd[newEdgeListIndOfOtherVertex]][0] && oldEdgesWithFacesTotalNbr[0]>=1)
                                // Must be
                                {
                                    if(!otherVertexXHasNoFreeFaceConnection[1])
                                    // Probable but faces would be invalid
                                    {
                                        mustBeOnlyOneNewEdge = true;
                                        newEdgeToRemoveLocalInd = otherZeroPntLocalIndToEdgeLocalInd[0];
                                    }
                                }
                            }
                        }
                        else
                            FatalErrorInFunction<< "Face must have either one or two neighbor cells!"<< exit(FatalError);
                    }
                    
                    if(mustBeThreeEdges)
                    {}
                    else if(mustBeOnlyOneNewEdge)
                    {
                        addedEdgeToDelete.append(edgeInd[newEdgeToRemoveLocalInd]);
                        addedEdgeToDelete.append(edgeInd[oldEdgeLocalInd[0]]);
                    }
                    else if(mustBeOneNewAndOneOldEdge)
                    {
                        addedEdgeToDelete.append(edgeInd[newEdgeToSingleRemoveLocalInd]);
                    }
                    else
                    {
                        FatalErrorInFunction<< "Completely false face can not happen for this zero point face"<< exit(FatalError);
                    }                
                }
                else if(problematicFaceNewPoints[i]==2)
                {
                    if(newEdgeLocalInd.size()!=1 || oldEdgeLocalInd.size()!=2)
                        FatalErrorInFunction<<"Incorrect number of cut edges!"<< exit(FatalError);
                    
                    label nonZeroOppositePntLocalVerticeInd = -1;
                    for(int j=0;j<basisFaces[i].size();j++)
                    {
                        if(pointDist[basisFaces[i][j]]==0)
                        {
                            if(nonZeroOppositePntLocalVerticeInd!=-1)
                                nonZeroOppositePntLocalVerticeInd = (j+2)%basisFaces[i].size();
                            else
                                FatalErrorInFunction<<"Double opposite non zero point. Can not happen!"<< exit(FatalError);
                        }
                    }
                    if(nonZeroOppositePntLocalVerticeInd==-1)
                        FatalErrorInFunction<<"No opposite non zero point. Can not happen!"<< exit(FatalError);
                    
                    DynamicList<label> nonZeroOtherPntLocalVerticeInd;
                    for(int j=0;j<basisFaces[i].size();j++)
                    {
                        if(j==nonZeroOppositePntLocalVerticeInd)
                            continue;
                        if(pointDist[basisFaces[i][j]]!=0)
                        {
                            if((nonZeroOtherPntLocalVerticeInd.size()==0)||(nonZeroOtherPntLocalVerticeInd.size()==1))
                                nonZeroOtherPntLocalVerticeInd.append(j);
                            else
                                FatalErrorInFunction<<"Other than two other zero point. Can not happen!"<< exit(FatalError);
                        }
                    }
                    if(nonZeroOtherPntLocalVerticeInd.size()!=2)
                        FatalErrorInFunction<<"Other than two other zero point. Can not happen!"<< exit(FatalError);
                    
                    label centralEdgeLocalInd = newEdgeLocalInd[nonZeroPntLocalIndTocutEdgeLocalInd[nonZeroOppositePntLocalVerticeInd]];
                    labelList otherEdgesLocalInd(2);
                    otherEdgesLocalInd[0] = newEdgeLocalInd[nonZeroPntLocalIndTocutEdgeLocalInd[nonZeroOtherPntLocalVerticeInd[0]]];
                    otherEdgesLocalInd[1] = newEdgeLocalInd[nonZeroPntLocalIndTocutEdgeLocalInd[nonZeroOtherPntLocalVerticeInd[1]]];
                    if(centralEdgeLocalInd<0||centralEdgeLocalInd>2||otherEdgesLocalInd[0]<0||otherEdgesLocalInd[0]>2||otherEdgesLocalInd[1]<0||otherEdgesLocalInd[1]>2)
                        FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                    if((centralEdgeLocalInd+otherEdgesLocalInd[0]+otherEdgesLocalInd[1])!=3)
                        FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                    
                    //two other vertex can not have free face connection
                    if(!otherVertexXHasNoFreeFaceConnection[0] || !otherVertexXHasNoFreeFaceConnection[1])
                        FatalErrorInFunction<<"Added cut points can not have free face connections!"<< exit(FatalError);
                    
                    if(newEdges.size()!=3 || oldEdges.size()!=0)
                        FatalErrorInFunction<<"Invalid number of new and old Edges!"<< exit(FatalError);
                    label newPointEdgeInd=-1;
                    DynamicList<label> oldPointEdgesInd;
                    for(int j=0;j<newEdges.size();j++)
                    {
                        if(newEdges[j].start()<nbrOfPrevPoints && newEdges[j].end()<nbrOfPrevPoints)
                            newPointEdgeInd = j;
                        else
                            oldPointEdgesInd.append(j);
                    }
                    if(newPointEdgeInd==-1 || oldPointEdgesInd.size()!=2)
                        FatalErrorInFunction<<"Invalid number of new and old Edges!"<< exit(FatalError);
                    
                    bool mustBeThreeEdges = false;
                    // three edge face is when all edges have exactly  a connected face to one cell but not to the other                
                    if(faceCells.size()==2)
                    {
                        bool validVerticalEdgeCuts = false;
                        if((connectedToCellPerEdge[centralEdgeLocalInd][0]  && connectedToCellPerEdge[otherEdgesLocalInd[0]][1]  && connectedToCellPerEdge[otherEdgesLocalInd[1]][1]) &&
                           (!connectedToCellPerEdge[centralEdgeLocalInd][1] && !connectedToCellPerEdge[otherEdgesLocalInd[0]][0] && !connectedToCellPerEdge[otherEdgesLocalInd[1]][0]))
                        {
                            if(verticalEdgesAreCut[nonZeroOppositePntLocalVerticeInd][0] && verticalEdgesAreCut[nonZeroOtherPntLocalVerticeInd[0]][1] && verticalEdgesAreCut[nonZeroOtherPntLocalVerticeInd[1]][1])
                                validVerticalEdgeCuts=true;
                        }
                        if((!connectedToCellPerEdge[centralEdgeLocalInd][0] && !connectedToCellPerEdge[otherEdgesLocalInd[0]][1] && !connectedToCellPerEdge[otherEdgesLocalInd[1]][1]) &&
                           (connectedToCellPerEdge[centralEdgeLocalInd][1]  && connectedToCellPerEdge[otherEdgesLocalInd[0]][0]  && connectedToCellPerEdge[otherEdgesLocalInd[1]][0]))
                        {
                            if(verticalEdgesAreCut[nonZeroOppositePntLocalVerticeInd][1] && verticalEdgesAreCut[nonZeroOtherPntLocalVerticeInd[0]][0] && verticalEdgesAreCut[nonZeroOtherPntLocalVerticeInd[1]][0])
                                validVerticalEdgeCuts=true;
                        }                      
                        if(validVerticalEdgeCuts)
                            mustBeThreeEdges = true;
                    }
                    else if(faceCells.size()==1)
                    {                    
                        bool validVerticalEdgeCuts = false;
                        if((connectedToCellPerEdge[centralEdgeLocalInd][0]) && !connectedToCellPerEdge[otherEdgesLocalInd[0]][0] && !connectedToCellPerEdge[otherEdgesLocalInd[1]][0])
                        {
                            if(verticalEdgesAreCut[nonZeroOppositePntLocalVerticeInd][0])
                                validVerticalEdgeCuts=true;
                        }
                        if((!connectedToCellPerEdge[centralEdgeLocalInd][0]) && (connectedToCellPerEdge[otherEdgesLocalInd[0]][0]  && connectedToCellPerEdge[otherEdgesLocalInd[1]][0]))
                        {
                            if(verticalEdgesAreCut[nonZeroOtherPntLocalVerticeInd[0]][0] && verticalEdgesAreCut[nonZeroOtherPntLocalVerticeInd[1]][0])
                                validVerticalEdgeCuts=true;
                        }                      
                        if(validVerticalEdgeCuts)
                            mustBeThreeEdges = true;
                    }
                    else
                        FatalErrorInFunction<< "Face must have either one or two neighbor cells!"<< exit(FatalError);
                    
                    bool mustBeOnlyOneEdge = false;
                    /* The edge connecting new point to new point
                    */
                    if(!mustBeThreeEdges)
                    {
                        if(!centerVertexHasNoFreeFaceConnection && otherVertexXHasFaceConnectionOnlyWithOppositeVertex[0] && otherVertexXHasFaceConnectionOnlyWithOppositeVertex[1])
                        {
                            if((faceCells.size()==2 && connectedToCellPerEdge[centralEdgeLocalInd][0] && connectedToCellPerEdge[centralEdgeLocalInd][0]) ||
                               (faceCells.size()==1 && connectedToCellPerEdge[centralEdgeLocalInd][0]))
                            {
                                mustBeOnlyOneEdge = true;
                            }
                        }
                    }
                    
                    bool mustBeTwoEdges = false;
                    /* The edges connected to the center vertex
                     */
                    if(!mustBeOnlyOneEdge)
                    {
                        if(otherVertexXHasFaceConnectionOnlyToCenter[0] && otherVertexXHasFaceConnectionOnlyToCenter[1] && centerVertexHasNoFreeFaceConnection)
                        {
                            if(faceCells.size()==2 && connectedToCellPerEdge[otherEdgesLocalInd[0]][0] && connectedToCellPerEdge[otherEdgesLocalInd[0]][1] && connectedToCellPerEdge[otherEdgesLocalInd[1]][0] && connectedToCellPerEdge[otherEdgesLocalInd[1]][1])
                            {
                                mustBeTwoEdges = true;
                            }
                            else if(faceCells.size()==1 && connectedToCellPerEdge[otherEdgesLocalInd[0]][0] && connectedToCellPerEdge[otherEdgesLocalInd[1]][0])
                            {
                                mustBeTwoEdges = true;
                            }
                        }
                    }
                    
                    if(mustBeThreeEdges)
                    {}
                    if(mustBeOnlyOneEdge)
                    {
                        addedEdgeToDelete.append(edgeInd[oldPointEdgesInd[0]]);
                        addedEdgeToDelete.append(edgeInd[oldPointEdgesInd[1]]);
                    }
                    else if(mustBeTwoEdges)
                    {
                        addedEdgeToDelete.append(edgeInd[newPointEdgeInd]);
                    }
                    else
                    {
                        FatalErrorInFunction<< "Completely false face can not happen for this zero point face"<< exit(FatalError);
                    }
                }
                else
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
            }
            else if(problematicFacePoints[i]==2 && problematicFaceNewPoints[i]==0)
            {
               if(edgeInd.size()!=1)
                    FatalErrorInFunction<< "No one added edges"<< exit(FatalError);
               
               edge oneEdge = newMeshEdges_[edgeInd[0]];
               
                DynamicList<label> zeroPoints{oneEdge.start(),oneEdge.end()};
                
                List<DynamicList<DynamicList<face>>> zeroPointsClosedFaces(zeroPoints.size());
                List<DynamicList<DynamicList<std::unordered_set<label>>>> zeroPointsClosedFaceMap(zeroPoints.size());
                for(int j=0;j<zeroPoints.size();j++)
                {
                    //collect all faces in neighboring cells connected to the point in posClosedFacesAroundPoint
                    const labelList& cellsAtPoint = pointToCells[zeroPoints[j]];
                    label nbrRelCells = cellsAtPoint.size();
                    List<DynamicList<face>> posClosedFacesAroundPoint(nbrRelCells);
                    List<DynamicList<std::unordered_set<label>>> posClosedPointMapAroundPoint(nbrRelCells);
                    for(int k=0;k<cellsAtPoint.size();k++)
                    {
                        DynamicList<DynamicList<std::unordered_set<label>>> thisCellFaceGroupsMap = cellNonConnectedMultiPointMap[cellsAtPoint[k]];
                        DynamicList<DynamicList<face>> thisCellFaceGroups = cellNonConnectedMultiFaces[cellsAtPoint[k]];
                        bool oneFaceContains = false;
                        for(int l=0;l<thisCellFaceGroupsMap.size();l++)
                        {
                            bool faceGroupContainsPoint = false;
                            for(int m=0;m<thisCellFaceGroupsMap[l].size();m++)
                            {
                                if(thisCellFaceGroupsMap[l][m].count(zeroPoints[j])!=0)
                                {
                                    if(oneFaceContains)
                                        FatalErrorInFunction<<"Seperate cut face groups share a point!"<< exit(FatalError);
                                    faceGroupContainsPoint = true;
                                }
                            }
                            if(faceGroupContainsPoint)
                            {
                                posClosedPointMapAroundPoint[k] = thisCellFaceGroupsMap[l];
                                posClosedFacesAroundPoint[k] = thisCellFaceGroups[l];
                                oneFaceContains = true;
                            }
                        }
                    }
                    computeClosedFaceFront(zeroPoints[j],posClosedFacesAroundPoint,posClosedPointMapAroundPoint,
                                           zeroPointsClosedFaces[j],zeroPointsClosedFaceMap[j]);
                }                
                List<List<bool>> vertexConnectedToVertex0(2);
                List<List<bool>> vertexConnectedToVertex1(2);                    
                for(int j=0;j<2;j++)
                {
                    vertexConnectedToVertex0[j] = pointInFaceFront(zeroPointsClosedFaces[j],zeroPointsClosedFaceMap[j],zeroPoints[0]);
                    vertexConnectedToVertex1[j] = pointInFaceFront(zeroPointsClosedFaces[j],zeroPointsClosedFaceMap[j],zeroPoints[1]);
                }
                List<bool> vertexXHasNoFreeFaceConnection(2,true);
                for(int j=0;j<vertexConnectedToVertex0[0].size();j++)
                {
                    if(!vertexConnectedToVertex0[0][j])
                        FatalErrorInFunction<< "Face front must have its own central Point!"<< exit(FatalError);
                }
                for(int j=0;j<vertexConnectedToVertex0[1].size();j++)
                {
                    if(!vertexConnectedToVertex0[1][j])
                        vertexXHasNoFreeFaceConnection[1] = false;
                }
                for(int j=0;j<vertexConnectedToVertex1[0].size();j++)
                {
                    if(!vertexConnectedToVertex1[0][j])
                        vertexXHasNoFreeFaceConnection[0] = false;
                }
                for(int j=0;j<vertexConnectedToVertex1[1].size();j++)
                {
                    if(!vertexConnectedToVertex1[1][j])
                        FatalErrorInFunction<< "Face front must have its own central Point!"<< exit(FatalError);
                }
                if(!vertexXHasNoFreeFaceConnection[0] && !vertexXHasNoFreeFaceConnection[1])
                    addedEdgeToDelete.append(edgeInd[0]);
            }
        }
    }
    
    std::unordered_set<label> mapEdgesToDelete;
    for(int i=0;i<addedEdgeToDelete.size();i++)
        mapEdgesToDelete.insert(addedEdgeToDelete[i]);
    
    DynamicList<edge> newMeshEdges_Cp;
    DynamicList<label> edgesToSide_Cp;
    DynamicList<DynamicList<label>> edgeToFaces_Cp;
    DynamicList<DynamicList<label>> edgeToCells_Cp;

    DynamicList<DynamicList<label>> faceToEdges_Cp;    
    DynamicList<DynamicList<label>> cellToEdges_Cp;
            
    for(int i=0;i<nbrOfPrevEdges;i++)
    {
        newMeshEdges_Cp.append(newMeshEdges_[i]);
        edgesToSide_Cp.append(edgesToSide_[i]);
        if(mapEdgesToDelete.count(i)==0)
        {
            edgeToFaces_Cp.append(edgeToFaces_[i]);
            edgeToCells_Cp.append(edgeToCells_[i]);
        }
        else
        {
            for(int j=0;j<edgeToFaces_[i].size();j++)
            {
                label faceInd = edgeToFaces_[i][j];
                DynamicList<label> fEdges = faceToEdges_[faceInd];
                DynamicList<label> fEdgesN;
                for(int k=0;k<fEdges.size();k++)
                    if(fEdges[k]!=i)
                       fEdgesN.append(fEdges[k]);
                faceToEdges_[faceInd] = fEdgesN;
            }
            for(int j=0;j<edgeToCells_[i].size();j++)
            {
                label cellInd = edgeToCells_[i][j];
                DynamicList<label> fEdges = cellToEdges_[cellInd];
                DynamicList<label> fEdgesN;
                for(int k=0;k<fEdges.size();k++)
                    if(fEdges[k]!=i)
                       fEdgesN.append(fEdges[k]);
                cellToEdges_[cellInd] = fEdgesN;
            }
        }
    }
    for(int i=nbrOfPrevEdges;i<newMeshEdges_Cp.size();i++)
    {
        if(mapEdgesToDelete.count(i)==0)
        {
            edgeToFaces_Cp.append(edgeToFaces_[i]);
            edgeToCells_Cp.append(edgeToCells_[i]);
            newMeshEdges_Cp.append(newMeshEdges_[i]);
            edgesToSide_Cp.append(edgesToSide_[i]);
        }
        else
        {
            for(int j=0;j<edgeToFaces_[i].size();j++)
            {
                label faceInd = edgeToFaces_[i][j];
                DynamicList<label> fEdges = faceToEdges_[faceInd];
                DynamicList<label> fEdgesN;
                for(int k=0;k<fEdges.size();k++)
                    if(fEdges[k]!=i)
                       fEdgesN.append(fEdges[k]);
                faceToEdges_[faceInd] = fEdgesN;
            }
            for(int j=0;j<edgeToCells_[i].size();j++)
            {
                label cellInd = edgeToCells_[i][j];
                DynamicList<label> fEdges = cellToEdges_[cellInd];
                DynamicList<label> fEdgesN;
                for(int k=0;k<fEdges.size();k++)
                    if(fEdges[k]!=i)
                       fEdgesN.append(fEdges[k]);
                cellToEdges_[cellInd] = fEdgesN;
            }
        }
    }        
    
    Info<<"Prev Nbr Edges:"<<newMeshEdges_.size()<<endl;
    Info<<"After Nbr Edges:"<<newMeshEdges_Cp.size()<<endl;
    
    newMeshEdges_ = newMeshEdges_Cp;
    edgesToSide_ = edgesToSide_Cp;
    edgeToFaces_ = edgeToFaces_Cp;
    edgeToCells_ = edgeToCells_Cp;
    faceToEdges_ = faceToEdges_Cp;
    cellToEdges_ = cellToEdges_Cp;
    
    
    //FatalErrorInFunction<< "Temp stop end"<< exit(FatalError);
}

void Foam::cutCellFvMesh::findCycles
(
    int i,
    label startPoint,
    label nextPoint,
    label prevPoint,
    DynamicList<label> cyclePath,
    DynamicList<label> cycleEdgePath,
    std::unordered_set<label> coveredPoints,
    std::unordered_set<label> usedEdges,
    std::unordered_map<label,std::unordered_set<label>>& pointGraphData,
    DynamicList<std::pair<label,label>>& pointEdgeComb,
    DynamicList<DynamicList<label>>& closedCyclePoints,
    DynamicList<std::unordered_set<label>>& closedCycleEdges,
    DynamicList<DynamicList<label>>& closedCycleEdgesList
)
{
    if(i==589240)
        Info<<"startPoint:"<<startPoint<<" nextPoint:"<<nextPoint<<" prevPoint:"<<prevPoint<<endl;
    if(nextPoint==-1 && prevPoint==-1)
    {
        for(const label val : pointGraphData[startPoint]) 
        {
            std::pair<label,label> posNext = pointEdgeComb[val];
            cyclePath.append(startPoint);
            cycleEdgePath.append(posNext.second);
            if(i==589240)
                Info<<"val:"<<val<<" nextPoint: "<<posNext.first<<" nextEdge: "<<posNext.second<<"  append:"<<startPoint<<endl;
            coveredPoints.insert(startPoint);
            usedEdges.insert(posNext.second);
            if(coveredPoints.count(posNext.first)!=0)
                FatalErrorInFunction<< "Can not happen!"<< exit(FatalError);
            
            if(i==589240)
                Info<<"Start: "<<startPoint<<"->"<<posNext.first<<endl;
            findCycles(i,startPoint,posNext.first,startPoint,cyclePath,cycleEdgePath,coveredPoints,usedEdges,
                       pointGraphData,pointEdgeComb,closedCyclePoints,closedCycleEdges,closedCycleEdgesList);
            
            cyclePath.setSize(cyclePath.size()-1);
            cycleEdgePath.setSize(cycleEdgePath.size()-1);
            coveredPoints.erase(startPoint);
            usedEdges.erase(posNext.second);
        }
    }
    else
    {
        if(nextPoint == startPoint)
        {
            if(i==589240)
                Info<<"Possible end?"<<endl;

            bool faceOld = false;
            for(int i=0;i<closedCycleEdges.size();i++)
            {
                if(equalEdges(closedCycleEdges[i],usedEdges))
                    faceOld = true;
            }
            if(!faceOld)
            {
                closedCyclePoints.append(cyclePath);
                closedCycleEdgesList.append(cycleEdgePath);
                closedCycleEdges.append(usedEdges);
                if(i==589240)
                    Info<<"Closed face: "<<cyclePath<<endl;
            }
        }
        else
        {
            for(const label val : pointGraphData[nextPoint]) 
            {
                std::pair<label,label> posNext = pointEdgeComb[val];
                if((coveredPoints.count(posNext.first)==0 && posNext.first != prevPoint) ||
                   (posNext.first == startPoint && posNext.first != prevPoint))
                {
                    cyclePath.append(nextPoint);
                    if(i==589240)
                        Info<<"val:"<<val<<" posNext.first: "<<posNext.first<<" posNext.second: "<<posNext.second<<"  append:"<<nextPoint<<endl;
                    cycleEdgePath.append(posNext.second);
                    coveredPoints.insert(nextPoint);
                    usedEdges.insert(posNext.second);

                    if(i==589240)
                        Info<<nextPoint<<"->"<<posNext.first<<endl;
                    findCycles(i,startPoint,posNext.first,nextPoint,cyclePath,cycleEdgePath,coveredPoints,usedEdges,
                               pointGraphData,pointEdgeComb,closedCyclePoints,closedCycleEdges,closedCycleEdgesList);
                    if(i==589240)
                        Info<<"-----split-------"<<cyclePath<<endl;
                    cyclePath.setSize(cyclePath.size()-1);
                    cycleEdgePath.setSize(cycleEdgePath.size()-1);
                    coveredPoints.erase(nextPoint);
                    usedEdges.erase(posNext.second);
                    if(i==589240)
                        Info<<"-----split-------"<<cyclePath<<"-----------------------"<<endl;
                }
            }
        }
    }
    if(i==589240)
        Info<<"Back"<<"   startPoint:"<<startPoint<<" nextPoint:"<<nextPoint<<" prevPoint:"<<prevPoint<<endl;
}

bool Foam::cutCellFvMesh::equalEdges
(
    std::unordered_set<label> setA,
    std::unordered_set<label> setB
)
{
    for(const label k : setA)
    {
        if(setB.count(k)==0)
            return false;
    }
    for(const label k : setB)
    {
        if(setA.count(k)==0)
            return false;
    }
    return true;
}


void Foam::cutCellFvMesh::printAddedEdges
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

void Foam::cutCellFvMesh::newMeshFaces
(
)
{
    //Info<<"Starting adding Faces"<<endl;
    const cellList& meshCells = this->cells();
    //const pointField& basisPoints = this->points();
    const faceList& basisFaces = this->faces();
    const edgeList& basisEdges = this->edges();
    
    nbrOfPrevFaces = basisFaces.size();
    
    newMeshFaces_.setCapacity(basisFaces.size()*2);
    newMeshFaces_.append(basisFaces);

    facesToSide(newMeshFaces_);
    
    facesToSide_.setCapacity(basisFaces.size()*2);
    
    /*
     * A existing face is a cut face if all its points are cut points.
     * Only if thats true the respective face is connected to its connected
     * cells via the faceToCells_ list. A non empty faceToCells_ for a face 
     * is the sign that the face is a cut face
     */
    faceToCells_.setSize(basisFaces.size());
    const labelList& owner = this->faceOwner();
    const labelList& neighbour = this->faceNeighbour();
    for(int i=0;i<basisFaces.size();i++)
    {
        labelList facePoints = basisFaces[i];

        bool isCutFace = true;
        for(int k=0;k<facePoints.size();k++)
        {
            if(pointsToSide_[facePoints[k]] != 0)
                isCutFace = false;
        }
        //Info<<"\t"<<"is CutFace:"<<isCutFace<<endl;
        if(isCutFace)
        {
            faceToCells_[i].setCapacity(2);
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
    cellToFaces_.setSize(meshCells.size());
    for(int i=0;i<meshCells.size();i++)
    {   
        DynamicList<label> facePoints;
        facePoints.setCapacity(10);
        //labelList facePoints;
        labelList cellCutEdgeList = cellToEdges_[i];
        /* 1)
         * If a cell has no cut edges it is not cut by a face
         */
        if(cellToEdges_[i].size() == 0)
        {
            continue;
        }
        
        /* 2)
         * if cell has one cut edge thats an old edge that cell is not cut by a face
         */
        if(cellToEdges_[i].size() == 1)
        {
            if(cellToEdges_[i][0] < basisEdges.size())
            {
                continue;
            }
            else
            /* 3)
            * If a cell has one cut edge it is not cut. It should be impossible that a
            * cell has two cut edges because the minimum cut face thinkable is a three edge
            * face. Because of that a failure exit is called if these states appear 
            */
            {
            FatalErrorInFunction
            << "A cell cannot be cut by "<< cellCutEdgeList.size()
            << " edges! "<<"Cell is "<<i
            << exit(FatalError);
            }
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
                << exit(FatalError);
            }
            if(equalFace.size() == 0)
            {
                Info<<"\nCell "<<i<<" is cut by ";
                for(int i=0;i<cellCutEdgeList.size();i++)
                    Info<<cellCutEdgeList[i]<<" ";
                Info<<endl;
                FatalErrorInFunction
                << "Cell is cut by two old edges that do not"
                << " share a face! This must not happen! "
                << exit(FatalError);
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
                << exit(FatalError);
            }
            cellToFaces_[i].setSize(equalFace.size());
            for(int n=0;n<equalFace.size();n++) cellToFaces_[i][n] = equalFace[n];
        }
        
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
        facePoints.setCapacity(facePoints.size());
        //add to addedFaces
        face newFace(facePoints);
        //newFace = newFace.reverseFace();
        
        newMeshFaces_.append(newFace);
        facesToSide_.append(0);
        
        DynamicList<label> newFaceCell(1);
        newFaceCell[0] = i;
        faceToCells_.append(newFaceCell);
        
        //input to oldCellsToAddedFace
        cellToFaces_[i].append(newMeshFaces_.size()-1);
    }
    //Info<<"End adding faces"<<endl;
    newMeshFaces_.setCapacity(newMeshFaces_.size());
    faceToCells_.setCapacity(faceToCells_.size());
    cellToFaces_.setCapacity(cellToFaces_.size());
}

void Foam::cutCellFvMesh::newMeshFaces_plus
(
)
{
    //Info<<"Starting adding Faces"<<endl;
    const cellList& meshCells = this->cells();
    const pointField& basisPoints = this->points();
    const faceList& basisFaces = this->faces();
    const labelList& cellOwner = this->owner();
    const labelList& cellNeighbor = this->neighbour();
    const edgeList& basisEdges = this->edges();
    const labelListList& edgeToFaces = this->edgeFaces();
    const labelListList& faceToEdge = this->faceEdges();
    const labelListList& pointToFace = this->pointFaces();
    
    nbrOfPrevFaces = basisFaces.size();
    
    newMeshFaces_.setCapacity(basisFaces.size()*2);
    newMeshFaces_.append(basisFaces);

    facesToSide(newMeshFaces_);
    
    facesToSide_.setCapacity(basisFaces.size()*2);
    
    /*
     * A existing face is a cut face if all its points are cut points.
     * Only if thats true the respective face is connected to its connected
     * cells via the faceToCells_ list. A non empty faceToCells_ for a face 
     * is the sign that the face is a cut face
     */
    faceToCells_.setSize(basisFaces.size());
    const labelList& owner = this->faceOwner();
    const labelList& neighbour = this->faceNeighbour();
    for(int i=0;i<basisFaces.size();i++)
    {
        labelList facePoints = basisFaces[i];

        bool isCutFace = true;
        for(int k=0;k<facePoints.size();k++)
        {
            if(pointsToSide_[facePoints[k]] != 0)
                isCutFace = false;
        }
        //Info<<"\t"<<"is CutFace:"<<isCutFace<<endl;
        if(isCutFace)
        {
            faceToCells_[i].setCapacity(2);
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
    cellToFaces_.setSize(meshCells.size());
    for(int i=0;i<meshCells.size();i++)
    {
        labelList edgeOfCellList = cellToEdges_[i];
        bool facesToAdd = false;
        label countNewEdges = 0;
        for(int j=0;j<edgeOfCellList.size();j++)
        {
            if(edgeOfCellList[j]>=nbrOfPrevEdges)
            {
                countNewEdges++;
            }
        }
        if(countNewEdges>1)
        {
            facesToAdd = true;
        }

        
        if(i==196879)
        {
            Info<<endl;
            Info<<"edgeOfCellList: "<<edgeOfCellList<<endl;
            Info<<"nbrOfPrevEdges: "<<nbrOfPrevEdges<<endl;
            Info<<"facesToAdd: "<<facesToAdd<<endl;
            Info<<"edgeToFace 0: "<<edgeToFaces_[edgeOfCellList[0]]<<endl;
            Info<<"edgeToFace 1: "<<edgeToFaces_[edgeOfCellList[1]]<<endl;
            Info<<"edgeToFace 2: "<<edgeToFaces_[edgeOfCellList[2]]<<endl;
            Info<<"face: "<<basisFaces[edgeToFaces_[edgeOfCellList[2]][0]]<<endl;
            label faceInd = edgeToFaces_[edgeOfCellList[2]][0];
            for(int s=0;s<basisFaces[edgeToFaces_[edgeOfCellList[2]][0]].size();s++)
            {
                Info<<" "<<basisFaces[edgeToFaces_[edgeOfCellList[2]][0]][s]<<"|"<<
                pointsToSide_[basisFaces[edgeToFaces_[edgeOfCellList[2]][0]][s]]<<"|"<<pointDist[basisFaces[edgeToFaces_[edgeOfCellList[2]][0]][s]];
            }
            Info<<endl;
            Info<<"Face Owner: "<<cellOwner[faceInd]<<endl;
            Info<<"Face Neighbor: "<<cellNeighbor[faceInd]<<endl;
            
            Info<<"Owner cellTo Edges: "<<cellToEdges_[cellOwner[faceInd]]<<endl;
            Info<<"Neighbor cellTo Edges: "<<cellToEdges_[cellNeighbor[faceInd]]<<endl;
            
            Info<<"Owner points:";
            labelList ownerPoints = meshCells[cellOwner[faceInd]].labels(basisFaces);
            for(int s=0;s<ownerPoints.size();s++)
            {
                Info<<"-----"<<ownerPoints[s]<<basisPoints[ownerPoints[s]]<<"|"<<pointsToSide_[ownerPoints[s]]<<"|"<<pointDist[ownerPoints[s]];
            }
            Info<<endl;
            
            Info<<"Neighbor points:";
            labelList neighborPoints = meshCells[cellNeighbor[faceInd]].labels(basisFaces);
            for(int s=0;s<neighborPoints.size();s++)
            {
                Info<<"-----"<<neighborPoints[s]<<basisPoints[neighborPoints[s]]<<"|"<<pointsToSide_[neighborPoints[s]]<<"|"<<pointDist[neighborPoints[s]];
            }
            Info<<endl;
            bool flag;
            for(int s=0;s<edgeOfCellList.size();s++)
            {
                Info<<"Edge Center: "<<edgeOfCellList[s]<<"  :"<<distToNurbs(newMeshEdges_[edgeOfCellList[s]].centre(basisPoints),flag)<<endl;
            }
            
            //FatalErrorInFunction<< " Temp stop."<< exit(FatalError);
        }
        
        if(!facesToAdd)
            continue;
        
        if(i==196879)
        {
            Info<<endl;
            Info<<"edgeOfCellList: "<<edgeOfCellList<<endl;
            Info<<"nbrOfPrevEdges: "<<nbrOfPrevEdges<<endl;
            Info<<"facesToAdd: "<<facesToAdd<<endl;
            Info<<"edgeToFace 0: "<<edgeToFaces_[edgeOfCellList[0]]<<endl;
            Info<<"edgeToFace 1: "<<edgeToFaces_[edgeOfCellList[1]]<<endl;
            Info<<"edgeToFace 2: "<<edgeToFaces_[edgeOfCellList[2]]<<endl;
            FatalErrorInFunction<< " Temp stop."<< exit(FatalError);
        }
            
        bool allEdgesDone = false;
        std::unordered_set<label> edgeTaken;
        DynamicList<DynamicList<label>> closedEdgeLoops;
        label nbrFace=0;
        while(!allEdgesDone)
        {
            closedEdgeLoops.append(DynamicList<label>(0));
            
            //Append starting edge
            for(int j=0;j<edgeOfCellList.size();j++)
            {
                if(edgeTaken.count(edgeOfCellList[j]) == 0)
                {
                    closedEdgeLoops[nbrFace].append(edgeOfCellList[j]);
                    edgeTaken.insert(edgeOfCellList[j]);
                    break;
                }
            }
            if(closedEdgeLoops[nbrFace].size() != 1)
            {
                Info<<endl<<endl;
                Info<<closedEdgeLoops<<endl;
                Info<<edgeOfCellList<<endl;
                FatalErrorInFunction<< "Starting Point not found"<< exit(FatalError);
            }
            
            //Append all subsequent edges
            label faceEdgeCounter = 0;
            edge firstEdge = newMeshEdges_[closedEdgeLoops[nbrFace][0]];
            edge lastAddedEdge;
            bool faceClosed = false;
            do
            {
                lastAddedEdge = newMeshEdges_[closedEdgeLoops[nbrFace][faceEdgeCounter]];
                DynamicList<label> nextEdges;
                for(int j=0;j<edgeOfCellList.size();j++)
                {
                    if(edgeTaken.count(edgeOfCellList[j]) == 0 && newMeshEdges_[edgeOfCellList[j]].commonVertex(lastAddedEdge)!=-1)
                    {
                        nextEdges.append(edgeOfCellList[j]);
                    }
                }
                if(faceEdgeCounter==0)
                {
                    if(nextEdges.size() != 2)
                    {
                        Info<<"nbrFace: "<<nbrFace<<endl;
                        Info<<"edgeOfCellList: "<<edgeOfCellList<<endl;
                        Info<<"nextEdges: "<<nextEdges<<endl;
                        Info<<"firstEdge: "<<firstEdge<<endl;
                        Info<<"lastAddedEdge: "<<lastAddedEdge<<endl;
                        Info<<endl;
                        Info<<"nbrOfPrevEdges: "<<nbrOfPrevEdges<<endl;
                        Info<<"nbrOfPrevPoints: "<<nbrOfPrevPoints<<endl;
                        Info<<"cellFaces: "<<meshCells[i]<<endl;
                        for(int j=0;j<meshCells[i].size();j++)
                        {
                            Info<<"face:"<<meshCells[i][j]<<"  "<<basisFaces[meshCells[i][j]]<<endl;
                        }
                        for(int j=0;j<edgeOfCellList.size();j++)
                        {
                            Info<<"edgeNum: "<<edgeOfCellList[j]<<" edge: "<<newMeshEdges_[edgeOfCellList[j]]<<"on face: "<<edgeToFaces_[edgeOfCellList[j]]<<endl;
                        }
                        FatalErrorInFunction<< "Second edge selection in a face must have two options but has not"<< exit(FatalError);
                    }
                }
                else
                {
                    if(nextEdges.size() != 1)
                    {
                        Info<<endl;
                        Info<<"newMeshEdges_["<<1594182<<"]: "<<newMeshEdges_[1594182]<<endl;
                        Info<<"edgeOfCellList: "<<edgeOfCellList<<endl;
                        Info<<"lastAddedEdge: "<<lastAddedEdge<<endl;
                        Info<<"faceEdgeCounter: "<<faceEdgeCounter<<endl;
                        Info<<"closedEdgeLoops: "<<closedEdgeLoops<<endl;
                        Info<<"nextEdges: "<<nextEdges<<endl;
                        FatalErrorInFunction<< "Edge selection for third of further edge in a face must have one options but has not"<< exit(FatalError);
                    }
                }
                label nextEdge = nextEdges[0];
                closedEdgeLoops[nbrFace].append(nextEdge);
                edgeTaken.insert(nextEdge);
                faceEdgeCounter++;
                if(closedEdgeLoops[nbrFace].size()>2)
                {
                    edge lastAddedEdge = newMeshEdges_[nextEdge];
                    if(lastAddedEdge.commonVertex(firstEdge)!=-1)
                    {
                        faceClosed = true;
                    }
                }
            }
            while(!faceClosed);
            
            if(closedEdgeLoops[nbrFace].size()<3)
            {
                Info<<endl<<endl;
                Info<<closedEdgeLoops[nbrFace]<<endl;
                Info<<edgeOfCellList<<endl;
                FatalErrorInFunction<< "Closed edge with less than 3 edges"<<exit(FatalError);
            }

            //Check if edges remain untreated
            bool allDone = true;
            label numNotDone = 0;
            for(int j=0;j<edgeOfCellList.size();j++)
            {
                if(edgeTaken.count(edgeOfCellList[j]) == 0)
                {
                    allDone = false;
                    numNotDone++;
                }
            }
            if(numNotDone < 3 && !allDone)
            {
                Info<<endl<<endl;
                for (auto const& i: edgeTaken)
                    Info<< i << " ";
                for (auto const& i: edgeTaken)
                    Info<< edgeTaken.count(i) << " ";
                Info<<endl;
                Info<<closedEdgeLoops[nbrFace]<<endl;
                Info<<edgeOfCellList<<endl;
                FatalErrorInFunction<< "Less than three edges remaining!"<< exit(FatalError);
            }
            if(allDone)
            {
                allEdgesDone = true;
            }
            
            //Check face
            for(int j=0;j<closedEdgeLoops[nbrFace].size();j++)
            {
                edge nextEdge = newMeshEdges_[closedEdgeLoops[nbrFace][(j+1)%closedEdgeLoops[nbrFace].size()]];
                edge currEdge = newMeshEdges_[closedEdgeLoops[nbrFace][(j)%closedEdgeLoops[nbrFace].size()]];
                if(currEdge.commonVertex(nextEdge)==-1)
                {
                    FatalErrorInFunction<< "Wrong edge path for face!"<< exit(FatalError);
                }
            }
                      
            //Next step preparation
            nbrFace++;
        }
        
        DynamicList<DynamicList<label>> newFaceEdges;
        for(int j=0;j<closedEdgeLoops.size();j++)
        {
            bool newFace = false;
            for(int k=0;k<closedEdgeLoops[j].size();k++)
            {
                if(closedEdgeLoops[j][k]>=nbrOfPrevEdges)
                    newFace = true;
            }
            if(newFace)
            {
                newFaceEdges.append(closedEdgeLoops[j]);
            }
        }
        
        //Tranfer edge description of face to point description of face
        DynamicList<DynamicList<label>> newFacesPoints;
        for(int j=0;j<newFaceEdges.size();j++)
        {
            newFacesPoints.append(DynamicList<label>(0));
            for(int k=0;k<newFaceEdges[j].size();k++)
            {
                edge nextEdge = newMeshEdges_[newFaceEdges[j][(k+1)%newFaceEdges[j].size()]];
                edge currEdge = newMeshEdges_[newFaceEdges[j][k%newFaceEdges[j].size()]];
                label comVertex = currEdge.commonVertex(nextEdge);
                if(comVertex == -1)
                {
                    FatalErrorInFunction<< "Wrong edge path for face!"<< exit(FatalError);
                }
                newFacesPoints[j].append(comVertex);
            }
        }
        
        for(int j=0;j<newFacesPoints.size();j++)
        {
            face newFace(newFacesPoints[j]);
            
            newMeshFaces_.append(newFace);
            facesToSide_.append(0);
        
            DynamicList<label> newFaceCell(1);
            newFaceCell[0] = i;
            faceToCells_.append(newFaceCell);
            cellToFaces_[i].append(newMeshFaces_.size()-1);
        }
    }
    
    newMeshFaces_.setCapacity(newMeshFaces_.size());
    faceToCells_.setCapacity(faceToCells_.size());
    cellToFaces_.setCapacity(cellToFaces_.size());
}

void Foam::cutCellFvMesh::printAddedFaces
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

void Foam::cutCellFvMesh::cutOldFaces
(
)
{
    //const cellList& meshCells = this->cells();
    const faceList& meshFaces = this->faces();
    //const edgeList& meshEdges = this->edges();
    //const pointField& meshPoints = this->points();
    
    cellsToSide();
    
    cutFaces_.setCapacity(meshFaces.size());
    oldFacesToCutFaces_.setCapacity(meshFaces.size());
    cutFacesToSide_.setCapacity(meshFaces.size());
    
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
            
            /*
            Info<<endl<<"OrigFace "<<i<<" size:"<<currFace.size()<<"| ";
            for(int l=0;l<currFace.size();l++)
            {                
                label point = currFace[l];
                Info<<point<<newMeshPoints_[point]<<" ";
            }
            Info<<"Added Edge: "<<newMeshPoints_[addedEdge.start()]<<"->"<<newMeshPoints_[addedEdge.end()];
            Info<<endl;
            */
            
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
            
            /*
                Info<<"1 Face 1 size:"<<newFace1.size()<<"| ";
                for(int l=0;l<newFace1.size();l++)
                {                
                    label point = newFace1[l];
                    Info<<point<<newMeshPoints_[point]<<" ";
                }
                Info<<relPointIndex<<endl;
            */
            
            //Info<<"First part of face 1"<<endl;
            while(pointsToSide_[currPoint] == pointsToSide_[currFace.nextLabel(currPointIndex)])
            {
                currPoint = currFace.nextLabel(currPointIndex);
                currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                relPointIndex++;
                newFace1.append(currPoint);
    
                /*
                    Info<<"2 Face 1 size:"<<newFace1.size()<<"| ";
                    for(int l=0;l<newFace1.size();l++)
                    {                
                        label point = newFace1[l];
                        Info<<point<<newMeshPoints_[point]<<" ";
                    }
                    Info<<relPointIndex<<endl;
                */
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
                    << exit(FatalError);
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
                    << exit(FatalError);
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
                    << exit(FatalError);
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
                    << exit(FatalError);
                }
                //Info<<firstCutPoint<<endl;
            }
            else
            {
                FatalErrorInFunction
                << "Next point is neither a cut point nor a point at the " 
                << "other side of the cut."
                << exit(FatalError);
            }
            
            //Info<<firstCutPoint<<endl;

            //Info<<"Add cut points"<<endl;
            newFace1.append(firstCutPoint);
            newFace2.append(firstCutPoint);
            
            /*
            Info<<"3 Face 1 size:"<<newFace1.size()<<"| ";
            for(int l=0;l<newFace1.size();l++)
            {                
                label point = newFace1[l];
                Info<<point<<newMeshPoints_[point]<<" ";
            }
            Info<<relPointIndex<<endl;
            */
            
            currPoint = nextPoint;
            currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
            relPointIndex++;
            
            newFace2.append(currPoint);
            newFace2Sign = pointsToSide_[currPoint];
            
            /*
            Info<<"1 Face 2 size:"<<newFace2.size()<<"| ";
            for(int l=0;l<newFace2.size();l++)
            {                
                label point = newFace2[l];
                Info<<point<<newMeshPoints_[point]<<" ";
            }
            Info<<relPointIndex<<endl;
            */
            
            //Info<<"Create face 2"<<endl;
            while(pointsToSide_[currPoint] == pointsToSide_[currFace.nextLabel(currPointIndex)])
            {
                //Info<<newMeshPoints_[currPoint]<<pointsToSide_[currPoint]<<endl;
                //Info<<newMeshPoints_[currFace.nextLabel(currPointIndex)]<<pointsToSide_[currFace.nextLabel(currPointIndex)]<<endl;
                //Info<<"Index:"<<faceInd<<endl;
                currPoint = currFace.nextLabel(currPointIndex);
                currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                relPointIndex++;
                newFace2.append(currPoint);
                
                /*
                    Info<<"2 Face 2 size:"<<newFace2.size()<<"| ";
                    for(int l=0;l<newFace2.size();l++)
                    {                
                        label point = newFace2[l];
                        Info<<point<<newMeshPoints_[point]<<" ";
                    }
                    Info<<relPointIndex<<endl;
                */
                
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
            
                /*
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
                */

            
            //Info<<"Finish face 1"<<endl;
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
                    << exit(FatalError);
                }                    
            }
            
            /*
            Info<<"5 Face 1 size:"<<newFace1.size()<<"| ";
            for(int l=0;l<newFace1.size();l++)
            {                
                label point = newFace1[l];
                Info<<point<<newMeshPoints_[point]<<" ";
            }
            Info<<relPointIndex<<endl;
            */
            
            oldFacesToCutFaces_[i].setCapacity(2);
            cutFaces_.append(face(newFace1));
            cutFacesToSide_.append(newFace1Sign);
            oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
            cutFaces_.append(face(newFace2));
            cutFacesToSide_.append(newFace2Sign);
            oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
            
            /*
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
            */
        }
    }
    cutFaces_.setCapacity(cutFaces_.size());
    oldFacesToCutFaces_.setCapacity(oldFacesToCutFaces_.size());
    cutFacesToSide_.setCapacity(cutFacesToSide_.size());    
}

void Foam::cutCellFvMesh::cutOldFaces_plus
(
)
{
    //const cellList& meshCells = this->cells();
    const faceList& meshFaces = this->faces();
    //const edgeList& meshEdges = this->edges();
    //const pointField& meshPoints = this->points();
    const labelList& cellOwner = this->owner();
    const labelList& cellNeighbor = this->neighbour();
    
    cellsToSide();
    
    cutFaces_.setCapacity(meshFaces.size());
    oldFacesToCutFaces_.setCapacity(meshFaces.size());
    cutFacesToSide_.setCapacity(meshFaces.size());
    
    for(int i=0;i<meshFaces.size();i++)
    {
        label ownerCell = cellOwner[i];        
        DynamicList<face> ownerNewFaces;
        DynamicList<std::unordered_set<label>> ownerFacesPoints;
        ownerFacesPoints.setSize(cellToFaces_[ownerCell].size());
        for(int j=0;j<cellToFaces_[ownerCell].size();j++)
        {
            ownerNewFaces.append(newMeshFaces_[cellToFaces_[ownerCell][j]]);
            for(int k=0;k<newMeshFaces_[cellToFaces_[ownerCell][j]].size();k++)
            {
                ownerFacesPoints[j].insert(newMeshFaces_[cellToFaces_[ownerCell][j]][k]);
            }
        }

        label neighborCell = -1;
        DynamicList<face> neighborNewFaces;
        DynamicList<std::unordered_set<label>> neighborFacesPoints;
        if(i<cellNeighbor.size())
        {
            neighborCell = cellNeighbor[i];
            neighborFacesPoints.setSize(cellToFaces_[neighborCell].size());
            for(int j=0;j<cellToFaces_[neighborCell].size();j++)
            {
                ownerNewFaces.append(newMeshFaces_[cellToFaces_[neighborCell][j]]);
                for(int k=0;k<newMeshFaces_[cellToFaces_[neighborCell][j]].size();k++)
                {
                    neighborFacesPoints[j].insert(newMeshFaces_[cellToFaces_[neighborCell][j]][k]);
                }
            }
        }
        
        DynamicList<label> zeroEdgesOfThisFace = faceToEdges_[i];
        DynamicList<label> newZeroEdgesOfThisFace;
        for(int j=0;j<zeroEdgesOfThisFace.size();j++)
        {
            if(zeroEdgesOfThisFace[j] >= nbrOfPrevEdges)
            {
                newZeroEdgesOfThisFace.append(zeroEdgesOfThisFace[j]);
            }
        }
        /*
        if(i==585222)
        {
            Info<<endl;
            Info<<"zeroEdgesOfThisFace: "<<zeroEdgesOfThisFace<<endl;
            Info<<"newZeroEdgesOfThisFace: "<<newZeroEdgesOfThisFace<<endl;
        }
        */
        
        Info<<"i:"<<i<<"   zeroEdgesOfThisFace.size():"<<zeroEdgesOfThisFace.size()<<
        "  newZeroEdgesOfThisFace.size():"<<newZeroEdgesOfThisFace.size()<<endl;
        
        if(zeroEdgesOfThisFace.size()>4)
        {
            FatalErrorInFunction<< "Face has more than four zero edges! That can not happen"<< exit(FatalError);
        }
        else if(zeroEdgesOfThisFace.size()==4)
        {
            if(newZeroEdgesOfThisFace.size()!=0)
            {
                FatalErrorInFunction<< "Face has four zero edges but not all are old ones! "<< exit(FatalError);
            }
        }
        else if(zeroEdgesOfThisFace.size()==3)
        {
            if(newZeroEdgesOfThisFace.size()==3)
            {
                // Face is seperated into four triangles. One is a zero triangle.
                if(meshFaces[i].size()!=4)
                    FatalErrorInFunction<<"Face is not a quadrat. This must not happen!"<< exit(FatalError);
                labelList zeroPointsOfFace = faceToPoints_[i];
                if(zeroPointsOfFace.size()!=3)
                    FatalErrorInFunction<<"Face should have 3 zero Points but has not!"<< exit(FatalError);
                DynamicList<label> newZeroPoints;
                label oldZeroPoint = -1;
                for(int j=0;j<zeroPointsOfFace.size();j++)
                {
                    if(zeroPointsOfFace[j]<nbrOfPrevPoints)
                    {
                        if(oldZeroPoint==-1)
                            oldZeroPoint = zeroPointsOfFace[j];
                        else
                            FatalErrorInFunction<< "Face has more than one old zero Point. This can not happen! "<< exit(FatalError);
                    }
                    else
                    {
                        newZeroPoints.append(zeroPointsOfFace[j]);
                    }
                }
                if(oldZeroPoint==-1)
                    FatalErrorInFunction<<"Face has no old zero Point. This can not happen!"<< exit(FatalError);
                if(newZeroPoints.size()!=2)
                    FatalErrorInFunction<<"Face has not two new zero Point. This can not happen!"<< exit(FatalError);
                label oldZeroPointInd = meshFaces[i].which(oldZeroPoint);
                label oppositeFacePoint = meshFaces[i][(oldZeroPointInd+2)%meshFaces[i].size()];
                if(pointsToSide_[oppositeFacePoint] == 0)
                    FatalErrorInFunction<<"Non Zero Point is marked as zero!"<< exit(FatalError);
                
                DynamicList<label> edgesIndNeighboringOppositeFacePoint;
                DynamicList<edge> edgesNeighboringOppositeFacePoint;
                DynamicList<label> otherVertexOfEdgesNeighboringOppositeFacePoint;
                for(int j=0;j<newZeroPoints.size();j++)
                {
                    edgesIndNeighboringOppositeFacePoint.append(pointToEgde_[newZeroPoints[j]]);
                    edgesNeighboringOppositeFacePoint.append(newMeshEdges_[pointToEgde_[newZeroPoints[j]]]);
                    otherVertexOfEdgesNeighboringOppositeFacePoint.append(edgesNeighboringOppositeFacePoint.last().otherVertex(oppositeFacePoint));
                    if(otherVertexOfEdgesNeighboringOppositeFacePoint.last()==-1)
                        FatalErrorInFunction<<"Edges not connected to opposite face!"<< exit(FatalError);
                }
                if(edgesNeighboringOppositeFacePoint.size()!=2 &&
                   edgesNeighboringOppositeFacePoint[0].connected(edgesNeighboringOppositeFacePoint[1]))
                    FatalErrorInFunction<<"Edges not connected!"<< exit(FatalError);
                
                DynamicList<label> zeroPointList;
                zeroPointList.append(oldZeroPoint);
                zeroPointList.append(newZeroPoints);
                
                label matchingAddedOwnerFaces=0;
                std::unordered_set<label> takenOwnerFaces;
                List<bool>edgeWithNewOwnerFace(zeroEdgesOfThisFace.size(),false);
                for(int j=0;j<zeroPointList.size();j++)
                {
                    for(int k=0;k<ownerFacesPoints.size();k++)
                    {
                        if(ownerFacesPoints[k].count(zeroPointList[j])!=0 && ownerFacesPoints[k].count(zeroPointList[(j+1)%zeroPointList.size()])!=0)
                        {
                            if(takenOwnerFaces.count(k)!=0)
                                FatalErrorInFunction<<"Added face to more than one edge."<< exit(FatalError);
                            if(edgeWithNewOwnerFace[j])
                                FatalErrorInFunction<<"Added edge to more than one face."<< exit(FatalError);
                            takenOwnerFaces.insert(k);
                            matchingAddedOwnerFaces++;
                            edgeWithNewOwnerFace[j] = true;
                        }
                    }
                }
                label matchingAddedNeighborFaces=0;
                std::unordered_set<label> takenNeighborFaces;
                List<bool>edgeWithNewNeighborFace(zeroEdgesOfThisFace.size(),false);
                if(i<cellNeighbor.size())
                {
                    for(int j=0;j<zeroPointList.size();j++)
                    {
                        for(int k=0;k<neighborFacesPoints.size();k++)
                        {
                            if(neighborFacesPoints[k].count(zeroPointList[j])!=0 && neighborFacesPoints[k].count(zeroPointList[(j+1)%zeroPointList.size()])!=0)
                            {
                                if(takenNeighborFaces.count(k)!=0)
                                    FatalErrorInFunction<<"Added face to more than one edge."<< exit(FatalError);
                                if(edgeWithNewNeighborFace[j])
                                    FatalErrorInFunction<<"Added edge to more than one face."<< exit(FatalError);
                                takenNeighborFaces.insert(k);
                                matchingAddedNeighborFaces++;
                                edgeWithNewNeighborFace[j] = true;
                            }
                        }
                    }
                    for(int j=0;j<edgeWithNewNeighborFace.size();j++)
                    {
                        if(edgeWithNewNeighborFace[j] != edgeWithNewOwnerFace[j])
                            FatalErrorInFunction<<"Inconsistent face to edge for owner & neighbor."<< exit(FatalError);
                    }
                    if(matchingAddedNeighborFaces!=matchingAddedOwnerFaces)
                        FatalErrorInFunction<<"Inconsistent face to edge count."<< exit(FatalError);
                }
                
                if(matchingAddedOwnerFaces==3)
                {
                    oldFacesToCutFaces_[i].setCapacity(4);
                
                    labelList face1(3);
                    face1[0] = oldZeroPoint;
                    face1[1] = otherVertexOfEdgesNeighboringOppositeFacePoint[0];
                    face1[2] = newZeroPoints[0];
                    cutFaces_.append(face(face1));
                    cutFacesToSide_.append(pointsToSide_[face1[1]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                
                    labelList face2(3);
                    face2[0] = oldZeroPoint;
                    face2[1] = newZeroPoints[0];
                    face2[2] = newZeroPoints[1];
                    cutFaces_.append(face(face2));
                    cutFacesToSide_.append(pointsToSide_[face2[1]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                
                    labelList face3(3);
                    face3[0] = oldZeroPoint;
                    face3[1] = otherVertexOfEdgesNeighboringOppositeFacePoint[1];
                    face3[2] = newZeroPoints[1];
                    cutFaces_.append(face(face3));
                    cutFacesToSide_.append(pointsToSide_[face3[1]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                
                    labelList face4(3);
                    face4[0] = oppositeFacePoint;
                    face4[1] = newZeroPoints[0];
                    face4[2] = newZeroPoints[1];
                    cutFaces_.append(face(face4));
                    cutFacesToSide_.append(pointsToSide_[face4[0]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);             
                }
                else if(matchingAddedOwnerFaces==2)
                {
                    oldFacesToCutFaces_[i].setCapacity(3);
                    if(!edgeWithNewNeighborFace[1])
                        FatalErrorInFunction<<"New Point edge must have face"<< exit(FatalError);
                    
                    labelList face4(3);
                    face4[0] = oppositeFacePoint;
                    face4[1] = newZeroPoints[0];
                    face4[2] = newZeroPoints[1];
                    cutFaces_.append(face(face4));
                    cutFacesToSide_.append(pointsToSide_[face4[0]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                    
                    if(edgeWithNewNeighborFace[0])
                    {                
                        labelList face1(3);
                        face1[0] = oldZeroPoint;
                        face1[1] = otherVertexOfEdgesNeighboringOppositeFacePoint[0];
                        face1[2] = newZeroPoints[0];
                        cutFaces_.append(face(face1));
                        cutFacesToSide_.append(pointsToSide_[face1[1]]);
                        oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                        
                        labelList face2(4);
                        face2[0] = oldZeroPoint;
                        face2[1] = newZeroPoints[0];
                        face2[2] = newZeroPoints[1];
                        face2[3] = otherVertexOfEdgesNeighboringOppositeFacePoint[1];
                        cutFaces_.append(face(face2));
                        cutFacesToSide_.append(pointsToSide_[face2[3]]);
                        oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                    }
                    else if(edgeWithNewNeighborFace[2])
                    {
                        labelList face3(3);
                        face3[0] = oldZeroPoint;
                        face3[1] = otherVertexOfEdgesNeighboringOppositeFacePoint[1];
                        face3[2] = newZeroPoints[1];
                        cutFaces_.append(face(face3));
                        cutFacesToSide_.append(pointsToSide_[face3[1]]);
                        oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                        
                        labelList face2(4);
                        face2[0] = oldZeroPoint;
                        face2[1] = newZeroPoints[1];
                        face2[2] = newZeroPoints[0];
                        face2[3] = otherVertexOfEdgesNeighboringOppositeFacePoint[0];
                        cutFaces_.append(face(face2));
                        cutFacesToSide_.append(pointsToSide_[face2[3]]);
                        oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                    }
                    else
                        FatalErrorInFunction<<"Can not happen."<< exit(FatalError);
                }
                else if(matchingAddedOwnerFaces==1)
                {
                    oldFacesToCutFaces_[i].setCapacity(2);
                    if(edgeWithNewNeighborFace[1] == false)
                        FatalErrorInFunction<<"New Point edge must have face"<< exit(FatalError);

                    labelList face4(3);
                    face4[0] = oppositeFacePoint;
                    face4[1] = newZeroPoints[0];
                    face4[2] = newZeroPoints[1];
                    cutFaces_.append(face(face4));
                    cutFacesToSide_.append(pointsToSide_[face4[0]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                    
                    labelList face2(5);
                    face2[0] = oldZeroPoint;
                    face2[1] = otherVertexOfEdgesNeighboringOppositeFacePoint[0];
                    face2[2] = newZeroPoints[0];
                    face2[3] = newZeroPoints[1];
                    face2[4] = otherVertexOfEdgesNeighboringOppositeFacePoint[1];
                    cutFaces_.append(face(face2));
                    cutFacesToSide_.append(pointsToSide_[face2[1]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                }
                else
                    FatalErrorInFunction<<"Can not happen."<< exit(FatalError);
            }
            else if(newZeroEdgesOfThisFace.size()==2)
            {
                // Face is seperated into three triangles. One is a zero triangle.
                if(meshFaces[i].size()!=4)
                    FatalErrorInFunction<<"Face is not a quadrat. This must not happen!"<< exit(FatalError);
                labelList zeroPointsOfFace = faceToPoints_[i];
                if(zeroPointsOfFace.size()!=3)
                    FatalErrorInFunction<<"Face should have 3 zero Points but has not!"<< exit(FatalError);
                label newZeroPoints = -1;
                DynamicList<label> oldZeroPoint;
                for(int j=0;j<zeroPointsOfFace.size();j++)
                {
                    if(zeroPointsOfFace[j]>=nbrOfPrevPoints)
                    {
                        if(newZeroPoints==-1)
                            newZeroPoints = zeroPointsOfFace[j];
                        else
                            FatalErrorInFunction<< "Face has more than one new zero Point. This can not happen! "<< exit(FatalError);
                    }
                    else
                    {
                        oldZeroPoint.append(zeroPointsOfFace[j]);
                    }
                }
                if(newZeroPoints==-1)
                    FatalErrorInFunction<<"Face has no old zero Point. This can not happen!"<< exit(FatalError);
                if(oldZeroPoint.size()!=2)
                    FatalErrorInFunction<<"Face has not two new zero Point. This can not happen!"<< exit(FatalError);

                DynamicList<label> otherVertices(oldZeroPoint.size());
                for(int j=0;j<oldZeroPoint.size();j++)
                {
                    label otherInd = (j+1)%oldZeroPoint.size();
                    label localInd = meshFaces[i].which(oldZeroPoint[j]);
                    label nextInd = (localInd+1)%meshFaces[i].size();
                    label prevInd = (localInd-1)%meshFaces[i].size();
                    if(meshFaces[i][nextInd]==oldZeroPoint[otherInd])
                    {
                        otherVertices[j] = meshFaces[i][prevInd];
                    }
                    else if(meshFaces[i][prevInd]==oldZeroPoint[otherInd])
                    {
                        otherVertices[j] = meshFaces[i][nextInd];
                    }
                }
                
                DynamicList<label> zeroPointList;
                zeroPointList.append(oldZeroPoint);
                zeroPointList.append(newZeroPoints);
                
                label matchingAddedOwnerFaces=0;
                std::unordered_set<label> takenOwnerFaces;
                List<bool>edgeWithNewOwnerFace(zeroEdgesOfThisFace.size(),false);
                for(int j=0;j<zeroPointList.size();j++)
                {
                    for(int k=0;k<ownerFacesPoints.size();k++)
                    {
                        if(ownerFacesPoints[k].count(zeroPointList[j])!=0 && ownerFacesPoints[k].count(zeroPointList[(j+1)%zeroPointList.size()])!=0)
                        {
                            if(takenOwnerFaces.count(k)!=0)
                                FatalErrorInFunction<<"Added face to more than one edge."<< exit(FatalError);
                            if(edgeWithNewOwnerFace[j])
                                FatalErrorInFunction<<"Added edge to more than one face."<< exit(FatalError);
                            takenOwnerFaces.insert(k);
                            matchingAddedOwnerFaces++;
                            edgeWithNewOwnerFace[j] = true;
                        }
                    }
                }
                label matchingAddedNeighborFaces=0;
                std::unordered_set<label> takenNeighborFaces;
                List<bool>edgeWithNewNeighborFace(zeroEdgesOfThisFace.size(),false);
                if(i<cellNeighbor.size())
                {
                    for(int j=0;j<zeroPointList.size();j++)
                    {
                        for(int k=0;k<neighborFacesPoints.size();k++)
                        {
                            if(neighborFacesPoints[k].count(zeroPointList[j])!=0 && neighborFacesPoints[k].count(zeroPointList[(j+1)%zeroPointList.size()])!=0)
                            {
                                if(takenNeighborFaces.count(k)!=0)
                                    FatalErrorInFunction<<"Added face to more than one edge."<< exit(FatalError);
                                if(edgeWithNewNeighborFace[j])
                                    FatalErrorInFunction<<"Added edge to more than one face."<< exit(FatalError);
                                takenNeighborFaces.insert(k);
                                matchingAddedNeighborFaces++;
                                edgeWithNewNeighborFace[j] = true;
                            }
                        }
                    }
                    for(int j=0;j<edgeWithNewNeighborFace.size();j++)
                    {
                        if(edgeWithNewNeighborFace[j] != edgeWithNewOwnerFace[j])
                            FatalErrorInFunction<<"Inconsistent face to edge for owner & neighbor."<< exit(FatalError);
                    }
                    if(matchingAddedNeighborFaces!=matchingAddedOwnerFaces)
                        FatalErrorInFunction<<"Inconsistent face to edge count."<< exit(FatalError);
                }
                
                if(!edgeWithNewOwnerFace[0] && !edgeWithNewOwnerFace[0])
                    FatalErrorInFunction<<"One or more edges at new point without face!"<< exit(FatalError);

                oldFacesToCutFaces_[i].setCapacity(3);

                if(edgeWithNewOwnerFace[0])
                {               
                    labelList face1(3);
                    face1[0] = oldZeroPoint[0];
                    face1[1] = otherVertices[0];
                    face1[2] = newZeroPoints;
                    cutFaces_.append(face(face1));
                    cutFacesToSide_.append(pointsToSide_[face1[1]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                
                    labelList face2(4);
                    face2[0] = oldZeroPoint[0];;
                    face2[1] = oldZeroPoint[1];
                    face2[2] = otherVertices[1];
                    face2[3] = newZeroPoints;
                    cutFaces_.append(face(face2));
                    cutFacesToSide_.append(pointsToSide_[face2[2]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                }
                else if(edgeWithNewOwnerFace[2])
                {                
                    labelList face2(3);
                    face2[0] = oldZeroPoint[1];;
                    face2[1] = otherVertices[1];
                    face2[2] = newZeroPoints;
                    cutFaces_.append(face(face2));
                    cutFacesToSide_.append(pointsToSide_[face2[1]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                
                    labelList face3(4);
                    face3[0] = oldZeroPoint[0];
                    face3[1] = oldZeroPoint[1];
                    face3[2] = newZeroPoints;
                    face3[3] = otherVertices[0];
                    cutFaces_.append(face(face3));
                    cutFacesToSide_.append(pointsToSide_[face3[3]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                }
                else if(edgeWithNewOwnerFace[0] && edgeWithNewOwnerFace[2])
                {
                    labelList face1(3);
                    face1[0] = oldZeroPoint[0];
                    face1[1] = otherVertices[0];
                    face1[2] = newZeroPoints;
                    cutFaces_.append(face(face1));
                    cutFacesToSide_.append(pointsToSide_[face1[1]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                
                    labelList face2(3);
                    face2[0] = oldZeroPoint[1];;
                    face2[1] = otherVertices[1];
                    face2[2] = newZeroPoints;
                    cutFaces_.append(face(face2));
                    cutFacesToSide_.append(pointsToSide_[face2[1]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                
                    labelList face3(3);
                    face3[0] = oldZeroPoint[0];
                    face3[1] = oldZeroPoint[1];
                    face3[2] = newZeroPoints;
                    cutFaces_.append(face(face3));
                    cutFacesToSide_.append(pointsToSide_[face3[1]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                }
            }
            else if(newZeroEdgesOfThisFace.size()==1)
            {
                // Face is seperated into two triangles. One is a zero triangle.
                if(meshFaces[i].size()!=4)
                    FatalErrorInFunction<<"Face is not a quadrat. This must not happen!"<< exit(FatalError);
                labelList zeroPointsOfFace = faceToPoints_[i];
                if(zeroPointsOfFace.size()!=3)
                    FatalErrorInFunction<<"Face should have 3 zero Points but has not!"<< exit(FatalError);
                
                label nonZeroFacePoint = -1;
                DynamicList<label> zeroBoundaryPoints;
                label centralZeroPoint = -1;
                Info<<"-------------------------------------"<<endl;
                for(int j=1;j<meshFaces[i].size()+1;j++)
                {
                    label thisInd = (j)%meshFaces[i].size();
                    if(pointsToSide_[meshFaces[i][thisInd]]!=0)
                    {
                        if(nonZeroFacePoint!=-1)
                            FatalErrorInFunction<< "Face more than one nonZeroFacePoint!"<< exit(FatalError);
                        nonZeroFacePoint = meshFaces[i][thisInd];
                    }
                    else
                    {
                        label nextInd = (j+1)%meshFaces[i].size();
                        label prevInd = (j-1)%meshFaces[i].size();
                        if(pointsToSide_[meshFaces[i][nextInd]]==0 && pointsToSide_[meshFaces[i][prevInd]]==0 
                           && pointsToSide_[meshFaces[i][thisInd]]==0)
                        {
                            if(centralZeroPoint!=-1)
                                FatalErrorInFunction<< "Face more than one centralZeroPoint!"<< exit(FatalError);
                            centralZeroPoint = meshFaces[i][thisInd];
                            zeroBoundaryPoints.append(meshFaces[i][nextInd]);
                            zeroBoundaryPoints.append(meshFaces[i][prevInd]);
                        }
                        Info<<endl;
                        Info<<"j: "<<j<<endl;
                        Info<<"prevInd: "<<prevInd<<" = "<<pointsToSide_[meshFaces[i][prevInd]]<<endl;
                        Info<<"thisInd: "<<thisInd<<" = "<<pointsToSide_[meshFaces[i][thisInd]]<<endl;
                        Info<<"nextInd: "<<nextInd<<" = "<<pointsToSide_[meshFaces[i][nextInd]]<<endl;
                    }
                }
                if(nonZeroFacePoint==-1)
                    FatalErrorInFunction<<"No nonZeroFacePoint found!"<< exit(FatalError);
                if(centralZeroPoint==-1)
                {
                    Info<<endl;
                    Info<<"nonZeroFacePoint: "<<nonZeroFacePoint<<endl;
                    Info<<"meshFaces[i]:"<<meshFaces[i]<<endl;
                    Info<<"Side: ";
                    for(int j=0;j<meshFaces[i].size();j++)
                        Info<<pointsToSide_[meshFaces[i][j]]<<"   ";
                    Info<<endl;
                    FatalErrorInFunction<<"No centralZeroPoint found!"<< exit(FatalError);
                }                    

                DynamicList<label> zeroPointList;
                zeroPointList.append(zeroBoundaryPoints);
                zeroPointList.append(centralZeroPoint);
                
                label matchingAddedOwnerFaces=0;
                std::unordered_set<label> takenOwnerFaces;
                List<bool>edgeWithNewOwnerFace(zeroEdgesOfThisFace.size(),false);
                for(int j=0;j<zeroPointList.size();j++)
                {
                    for(int k=0;k<ownerFacesPoints.size();k++)
                    {
                        if(ownerFacesPoints[k].count(zeroPointList[j])!=0 && ownerFacesPoints[k].count(zeroPointList[(j+1)%zeroPointList.size()])!=0)
                        {
                            if(takenOwnerFaces.count(k)!=0)
                                FatalErrorInFunction<<"Added face to more than one edge."<< exit(FatalError);
                            if(edgeWithNewOwnerFace[j])
                                FatalErrorInFunction<<"Added edge to more than one face."<< exit(FatalError);
                            takenOwnerFaces.insert(k);
                            matchingAddedOwnerFaces++;
                            edgeWithNewOwnerFace[j] = true;
                        }
                    }
                }
                label matchingAddedNeighborFaces=0;
                std::unordered_set<label> takenNeighborFaces;
                List<bool>edgeWithNewNeighborFace(zeroEdgesOfThisFace.size(),false);
                if(i<cellNeighbor.size())
                {
                    for(int j=0;j<zeroPointList.size();j++)
                    {
                        for(int k=0;k<neighborFacesPoints.size();k++)
                        {
                            if(neighborFacesPoints[k].count(zeroPointList[j])!=0 && neighborFacesPoints[k].count(zeroPointList[(j+1)%zeroPointList.size()])!=0)
                            {
                                if(takenNeighborFaces.count(k)!=0)
                                    FatalErrorInFunction<<"Added face to more than one edge."<< exit(FatalError);
                                if(edgeWithNewNeighborFace[j])
                                    FatalErrorInFunction<<"Added edge to more than one face."<< exit(FatalError);
                                takenNeighborFaces.insert(k);
                                matchingAddedNeighborFaces++;
                                edgeWithNewNeighborFace[j] = true;
                            }
                        }
                    }
                    for(int j=0;j<edgeWithNewNeighborFace.size();j++)
                    {
                        if(edgeWithNewNeighborFace[j] != edgeWithNewOwnerFace[j])
                            FatalErrorInFunction<<"Inconsistent face to edge for owner & neighbor."<< exit(FatalError);
                    }
                    if(matchingAddedNeighborFaces!=matchingAddedOwnerFaces)
                        FatalErrorInFunction<<"Inconsistent face to edge count."<< exit(FatalError);
                }
                
                if(edgeWithNewOwnerFace[0])
                {
                    oldFacesToCutFaces_[i].setCapacity(2);
                
                    labelList face1(3);
                    face1[0] = centralZeroPoint;
                    face1[1] = zeroBoundaryPoints[0];
                    face1[2] = zeroBoundaryPoints[1];
                    cutFaces_.append(face(face1));
                    cutFacesToSide_.append(pointsToSide_[face1[1]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                
                    labelList face2(3);
                    face2[0] = zeroBoundaryPoints[0];
                    face2[1] = nonZeroFacePoint;
                    face2[2] = zeroBoundaryPoints[1];
                    cutFaces_.append(face(face2));
                    cutFacesToSide_.append(pointsToSide_[face2[1]]);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                }
                
                /*
                if(i==585222)
                {
                    Info<<endl;
                    Info<<"face1: "<<face1<<endl;
                    Info<<"face2: "<<face2<<endl;
                    Info<<"meshFaces[i]: "<<meshFaces[i]<<endl;
                    Info<<"edgePoints: "<<newMeshEdges_[1594182]<<endl;
                    Info<<"zeroEdgesOfThisFace: "<<zeroEdgesOfThisFace<<endl;
                    Info<<"newZeroEdgesOfThisFace: "<<newZeroEdgesOfThisFace<<endl;
                    FatalErrorInFunction<< "Temp Stop"<< exit(FatalError);
                }
                */
            }
            else
            {
                FatalErrorInFunction<< "Face has three zero edges but all are old ones! "<< exit(FatalError);
            }
        }
        else if(zeroEdgesOfThisFace.size()==2)
        {
            if(newZeroEdgesOfThisFace.size()==2)
            {
                edge newEdge1 = newMeshEdges_[newZeroEdgesOfThisFace[0]]; 
                edge faceEdge1_OfNewEdge1 = newMeshEdges_[pointToEgde_[newEdge1.start()]];
                edge faceEdge2_OfNewEdge1 = newMeshEdges_[pointToEgde_[newEdge1.end()]];
                
                edge newEdge2 = newMeshEdges_[newZeroEdgesOfThisFace[1]]; 
                edge faceEdge1_OfNewEdge2 = newMeshEdges_[pointToEgde_[newEdge2.start()]];
                edge faceEdge2_OfNewEdge2 = newMeshEdges_[pointToEgde_[newEdge2.end()]];
                
                label thirdPointByNewEdge1 = -1;
                thirdPointByNewEdge1 = faceEdge1_OfNewEdge1.commonVertex(faceEdge2_OfNewEdge1);
                if(thirdPointByNewEdge1==-1)
                    FatalErrorInFunction<< "thirdPointByNewEdge1 not found! "<< exit(FatalError);
                
                label thirdPointByNewEdge2 = -1;
                thirdPointByNewEdge2 = faceEdge1_OfNewEdge2.commonVertex(faceEdge2_OfNewEdge2);
                if(thirdPointByNewEdge2==-1)
                    FatalErrorInFunction<< "thirdPointByNewEdge2 not found! "<< exit(FatalError);                
                
                labelList face1(3);
                face1[0] = newEdge1.start();
                face1[1] = thirdPointByNewEdge1;
                face1[2] = newEdge1.end();
                cutFaces_.append(face(face1));
                cutFacesToSide_.append(pointsToSide_[face1[1]]);
                oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                
                labelList face2(3);
                face2[0] = newEdge2.start();
                face2[1] = thirdPointByNewEdge2;
                face2[2] = newEdge2.end();
                cutFaces_.append(face(face2));
                cutFacesToSide_.append(pointsToSide_[face2[1]]);
                oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                
                labelList face3(6);
                face3[0] = newEdge1.start();
                face3[1] = newEdge1.end();
                face3[2] = faceEdge2_OfNewEdge1.otherVertex(thirdPointByNewEdge1);
                face3[3] = -1;
                face3[4] = -1;
                if(faceEdge2_OfNewEdge1.connected(faceEdge1_OfNewEdge2))
                {
                    face3[3] = newEdge2.start();
                    face3[4] = newEdge2.end();
                }
                else if(faceEdge2_OfNewEdge1.connected(faceEdge2_OfNewEdge2))
                {
                    face3[3] = newEdge2.end();
                    face3[4] = newEdge2.start();
                }
                else
                    FatalErrorInFunction<< "Fourth and fifth point not found! "<< exit(FatalError);
                face3[5] = faceEdge1_OfNewEdge1.otherVertex(thirdPointByNewEdge1);
                cutFaces_.append(face(face3));
                cutFacesToSide_.append(pointsToSide_[face3[2]]);
                oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                
                DynamicList<label> zeroPointList;
                zeroPointList.append(face3[0]);
                zeroPointList.append(face3[1]);
                zeroPointList.append(face3[3]);
                zeroPointList.append(face3[4]);
                
                label matchingAddedOwnerFaces=0;
                std::unordered_set<label> takenOwnerFaces;
                List<bool>edgeWithNewOwnerFace(zeroEdgesOfThisFace.size(),false);
                for(int j=0;j<zeroPointList.size();j++)
                {
                    for(int k=0;k<ownerFacesPoints.size();k++)
                    {
                        if(ownerFacesPoints[k].count(zeroPointList[j])!=0 && ownerFacesPoints[k].count(zeroPointList[(j+1)%zeroPointList.size()])!=0)
                        {
                            if(takenOwnerFaces.count(k)!=0)
                                FatalErrorInFunction<<"Added face to more than one edge."<< exit(FatalError);
                            if(edgeWithNewOwnerFace[j])
                                FatalErrorInFunction<<"Added edge to more than one face."<< exit(FatalError);
                            takenOwnerFaces.insert(k);
                            matchingAddedOwnerFaces++;
                            edgeWithNewOwnerFace[j] = true;
                        }
                    }
                }
                label matchingAddedNeighborFaces=0;
                std::unordered_set<label> takenNeighborFaces;
                List<bool>edgeWithNewNeighborFace(zeroEdgesOfThisFace.size(),false);
                if(i<cellNeighbor.size())
                {
                    for(int j=0;j<zeroPointList.size();j++)
                    {
                        for(int k=0;k<neighborFacesPoints.size();k++)
                        {
                            if(neighborFacesPoints[k].count(zeroPointList[j])!=0 && neighborFacesPoints[k].count(zeroPointList[(j+1)%zeroPointList.size()])!=0)
                            {
                                if(takenNeighborFaces.count(k)!=0)
                                    FatalErrorInFunction<<"Added face to more than one edge."<< exit(FatalError);
                                if(edgeWithNewNeighborFace[j])
                                    FatalErrorInFunction<<"Added edge to more than one face."<< exit(FatalError);
                                takenNeighborFaces.insert(k);
                                matchingAddedNeighborFaces++;
                                edgeWithNewNeighborFace[j] = true;
                            }
                        }
                    }
                    for(int j=0;j<edgeWithNewNeighborFace.size();j++)
                    {
                        if(edgeWithNewNeighborFace[j] != edgeWithNewOwnerFace[j])
                            FatalErrorInFunction<<"Inconsistent face to edge for owner & neighbor."<< exit(FatalError);
                    }
                    if(matchingAddedNeighborFaces!=matchingAddedOwnerFaces)
                        FatalErrorInFunction<<"Inconsistent face to edge count."<< exit(FatalError);
                }
                
                if(!edgeWithNewOwnerFace[0] || edgeWithNewOwnerFace[1] || 
                    !edgeWithNewOwnerFace[2] || edgeWithNewOwnerFace[3])
                {
                    FatalErrorInFunction<< "Invalid cell for new edge! "<< exit(FatalError);
                }
            }
            else
            {
                FatalErrorInFunction<< "Face has two zero edges but not all are new ones! "<< exit(FatalError);
            }
        }
        else if(zeroEdgesOfThisFace.size()==1)
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
            
            
                //Info<<"First part of face 1"<<endl;
                while(pointsToSide_[currPoint] == pointsToSide_[currFace.nextLabel(currPointIndex)])
                {
                    currPoint = currFace.nextLabel(currPointIndex);
                    currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                    relPointIndex++;
                    newFace1.append(currPoint);
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
                        << exit(FatalError);
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
                        << exit(FatalError);
                    }                
                }
                // Cut point is added point
                else if(pointsToSide_[currPoint] * pointsToSide_[currFace.nextLabel(currPointIndex)] < 0)
                {
                    nextPoint = currFace.nextLabel(currPointIndex);
                    labelList currPointEdges = this->pointEdges()[currPoint];
                    labelList nextPointEdges = this->pointEdges()[nextPoint];
                
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
                        << exit(FatalError);
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
                        << exit(FatalError);
                    }
                    //Info<<firstCutPoint<<endl;
                }
                else
                {
                    FatalErrorInFunction
                    << "Next point is neither a cut point nor a point at the " 
                    << "other side of the cut."
                    << exit(FatalError);
                }
            
                //Info<<firstCutPoint<<endl;

                //Info<<"Add cut points"<<endl;
                newFace1.append(firstCutPoint);
                newFace2.append(firstCutPoint);
            
                currPoint = nextPoint;
                currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                relPointIndex++;
            
                newFace2.append(currPoint);
                newFace2Sign = pointsToSide_[currPoint];
            
                //Info<<"Create face 2"<<endl;
                while(pointsToSide_[currPoint] == pointsToSide_[currFace.nextLabel(currPointIndex)])
                {
                    currPoint = currFace.nextLabel(currPointIndex);
                    currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                    relPointIndex++;
                    newFace2.append(currPoint);
                }
                //Info<<"Finished loop"<<endl;
                newFace2.append(secondCutPoint);
                newFace1.append(secondCutPoint);
            
                if(secondCutPoint < nbrOfPrevPoints)
                {
                    currPointIndex = (currPointIndex >= currFace.size()-1) ? 0 : currPointIndex+1;
                    relPointIndex++;
                }
            
                //Info<<"Finish face 1"<<endl;
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
                        << exit(FatalError);
                    }                    
                }
                
                Info<<"Test for split face"<<endl;
                DynamicList<label> zeroPointList;
                zeroPointList.append(addedEdge.start());
                zeroPointList.append(addedEdge.end());
                
                label matchingAddedOwnerFaces=0;
                std::unordered_set<label> takenOwnerFaces;
                List<bool>edgeWithNewOwnerFace(zeroEdgesOfThisFace.size(),false);
                for(int j=0;j<zeroPointList.size()-1;j++)
                {
                    Info<<"j:"<<j<<endl;
                    for(int k=0;k<ownerFacesPoints.size();k++)
                    {
                        Info<<"k:"<<k<<" /"<<ownerFacesPoints.size()<<endl;
                        if(ownerFacesPoints[k].count(zeroPointList[j])!=0 && ownerFacesPoints[k].count(zeroPointList[(j+1)%zeroPointList.size()])!=0)
                        {
                            Info<<"Has face"<<endl;
                            if(takenOwnerFaces.count(k)!=0)
                                FatalErrorInFunction<<"Added face to more than one edge."<< exit(FatalError);
                            Info<<"1"<<endl;
                            if(edgeWithNewOwnerFace[j])
                                FatalErrorInFunction<<"Added edge to more than one face."<< exit(FatalError);
                            Info<<"2"<<endl;
                            takenOwnerFaces.insert(k);
                            Info<<"3"<<endl;
                            matchingAddedOwnerFaces++;
                            Info<<"4"<<endl;
                            edgeWithNewOwnerFace[j] = true;
                        }
                        Info<<"No face"<<endl;
                    }
                }
                Info<<"Got owner info"<<endl;
                label matchingAddedNeighborFaces=0;
                std::unordered_set<label> takenNeighborFaces;
                List<bool>edgeWithNewNeighborFace(zeroEdgesOfThisFace.size(),false);
                if(i<cellNeighbor.size())
                {
                    for(int j=0;j<zeroPointList.size()-1;j++)
                    {
                        for(int k=0;k<neighborFacesPoints.size();k++)
                        {
                            if(neighborFacesPoints[k].count(zeroPointList[j])!=0 && neighborFacesPoints[k].count(zeroPointList[(j+1)%zeroPointList.size()])!=0)
                            {
                                if(takenNeighborFaces.count(k)!=0)
                                    FatalErrorInFunction<<"Added face to more than one edge."<< exit(FatalError);
                                if(edgeWithNewNeighborFace[j])
                                    FatalErrorInFunction<<"Added edge to more than one face."<< exit(FatalError);
                                takenNeighborFaces.insert(k);
                                matchingAddedNeighborFaces++;
                                edgeWithNewNeighborFace[j] = true;
                            }
                        }
                    }
                    for(int j=0;j<edgeWithNewNeighborFace.size();j++)
                    {
                        if(edgeWithNewNeighborFace[j] != edgeWithNewOwnerFace[j])
                            FatalErrorInFunction<<"Inconsistent face to edge for owner & neighbor."<< exit(FatalError);
                    }
                    if(matchingAddedNeighborFaces!=matchingAddedOwnerFaces)
                        FatalErrorInFunction<<"Inconsistent face to edge count."<< exit(FatalError);
                }
                Info<<"Get face info"<<endl;
                bool splitFace = false;
                if(addedEdge.start()>=nbrOfPrevPoints || addedEdge.end()>=nbrOfPrevPoints)
                {
                    if(!edgeWithNewOwnerFace[0])
                        FatalErrorInFunction<< "New point edge but no new face! "<< exit(FatalError);
                    else
                        splitFace = true;
                }
                else
                {
                    if(edgeWithNewOwnerFace[0])
                        splitFace = true;
                }
                if(splitFace)
                {
                    oldFacesToCutFaces_[i].setCapacity(2);
                    cutFaces_.append(face(newFace1));
                    cutFacesToSide_.append(newFace1Sign);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                    cutFaces_.append(face(newFace2));
                    cutFacesToSide_.append(newFace2Sign);
                    oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                }
            }
        }
    }
    cutFaces_.setCapacity(cutFaces_.size());
    oldFacesToCutFaces_.setCapacity(oldFacesToCutFaces_.size());
    cutFacesToSide_.setCapacity(cutFacesToSide_.size());    
}

void Foam::cutCellFvMesh::printCutFaces
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

void Foam::cutCellFvMesh::createNewMeshData
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
    
    oldSplittedCellToNewPlusCell = List<DynamicList<label>>(meshCells.size());
    oldSplittedCellToNewMinusCell = List<DynamicList<label>>(meshCells.size());
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
    cutCellsMinusAndPlus = labelListList(0);
    oldCellToMinusCutCell = labelList(meshCells.size());
    oldCellToPlusCutCell = labelList(meshCells.size());
    
    Info<<"Insert Split cell faces"<<endl;
    // Compute List of new faces splitting old cells
    label addedCutFacesNbr = 0;
    addedCutFaces = faceList(0);
    addedCutFaceNeighbor = labelList(0);
    addedCutFaceOwner = labelList(0);
    for(int i=0;i<cellToFaces_.size();i++)
    {
        oldCellToMinusCutCell[i] = -1;
        oldCellToPlusCutCell[i] = -1;
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
            
            
            labelList addedIndex = {addedCutFaces.size()-1};
            cutCellsMinusAndPlus.append(addedIndex);
            oldCellToMinusCutCell[i] = cutCellsMinusAndPlus.size()-1;
            cutCellsMinusAndPlus.append(addedIndex);
            oldCellToPlusCutCell[i] = cutCellsMinusAndPlus.size()-1;

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
    
    for(int i=0;i<oldCellToPlusCutCell.size();i++)
    {
        Info<<"Cell "<<i<<" has ";
        if(oldCellToMinusCutCell[i] != -1 && oldCellToPlusCutCell[i] != -1)
        {
            Info<<
            cutCellsMinusAndPlus[oldCellToMinusCutCell[i]].size()<<" minus cut faces and "<<
            cutCellsMinusAndPlus[oldCellToPlusCutCell[i]].size()<<" plus cut faces"<<endl;
        }
        else
        {
            Info<<
            oldCellToMinusCutCell[i]<<" and "<<
            oldCellToPlusCutCell[i]<<endl;
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
        
        Info<<"Face "<<" size: "<<faceToEdges_[i].size()<<" on side: "<<facesToSide_[i]<<" owner:"<<owner[i]<<
        " on side "<<cellsToSide_[owner[i]]<<" and neighbour:"<<neighbour[i]<<" on side "<<cellsToSide_[neighbour[i]]<<endl;
        for(int k=0;k<meshFaces[i].size();k++)
        {
            Info<<"point: "<<meshPoints[meshFaces[i][k]]<<" at dist: "<<pointDist[meshFaces[i][k]]<<endl;
        }
        
        if(faceToEdges_[i].size() == 1 && faceToEdges_[i][0] >= nbrOfPrevEdges)
        {
            if(oldFacesToCutFaces_[i].size() != 2)
            {
                FatalErrorInFunction
                << " Splitted interior cell is cut into"<<oldFacesToCutFaces_[i].size()
                << " faces instead of the expected 2."
                << exit(FatalError);
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
            
            if( cutCellsMinusAndPlus[oldCellToPlusCutCell[neighbour[i]]].size() == 0 ||
                cutCellsMinusAndPlus[oldCellToPlusCutCell[owner[i]]].size() == 0
                )
            {
                FatalErrorInFunction
                << " Split face interior inserted but cell has no cut face."
                << exit(FatalError);
            }
            cutCellsMinusAndPlus[oldCellToPlusCutCell[neighbour[i]]].append
            (
                splitAndUnsplitFacesInterior.size()-1+addedCutFaces.size()
            );
            cutCellsMinusAndPlus[oldCellToPlusCutCell[owner[i]]].append
            (
                splitAndUnsplitFacesInterior.size()-1+addedCutFaces.size()
            );

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

            if( cutCellsMinusAndPlus[oldCellToMinusCutCell[neighbour[i]]].size() == 0 ||
                cutCellsMinusAndPlus[oldCellToMinusCutCell[owner[i]]].size() == 0
                )
            {
                FatalErrorInFunction
                << " Split face interior inserted but cell has no cut face."
                << exit(FatalError);
            }
            cutCellsMinusAndPlus[oldCellToMinusCutCell[neighbour[i]]].append
            (
                splitAndUnsplitFacesInterior.size()-1+addedCutFaces.size()
            );
            cutCellsMinusAndPlus[oldCellToMinusCutCell[owner[i]]].append
            (
                splitAndUnsplitFacesInterior.size()-1+addedCutFaces.size()
            );
            
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
                
                if(cellsToSide_[neighbour[i]] == 0)
                {
                    if(cutCellsMinusAndPlus[oldCellToPlusCutCell[neighbour[i]]].size() == 0)
                    {
                        FatalErrorInFunction
                        << " Unsplit face interior inserted but cell has no cut face."
                        << exit(FatalError);
                    }
                    cutCellsMinusAndPlus[oldCellToPlusCutCell[neighbour[i]]].append
                    (
                        splitAndUnsplitFacesInterior.size()-1+addedCutFaces.size()
                    );
                }
                else if(cellsToSide_[owner[i]] == 0)
                {
                    if(cutCellsMinusAndPlus[oldCellToPlusCutCell[owner[i]]].size() == 0)
                    {
                        FatalErrorInFunction
                        << " Unsplit face interior inserted but cell has no cut face."
                        << exit(FatalError);
                    }
                    cutCellsMinusAndPlus[oldCellToPlusCutCell[owner[i]]].append
                    (
                        splitAndUnsplitFacesInterior.size()-1+addedCutFaces.size()
                    );
                }
                
                /*
                 * Error is wrong
                if( (oldCellToPlusCutCell[neighbour[i]] != -1 &&
                     oldCellToPlusCutCell[owner[i]] != -1)
                    ||
                    (cellsToSide_[neighbour[i]] == 0 && cellsToSide_[owner[i]] == 0)
                    )
                {
                    Info<<"Neighbor "<<neighbour[i]<<" Index: "<<oldCellToPlusCutCell[neighbour[i]]<<
                        " Size: "<<cutCellsMinusAndPlus[oldCellToPlusCutCell[neighbour[i]]].size()<<endl;
                    Info<<"Owner "<<owner[i]<<" Index: "<<oldCellToPlusCutCell[owner[i]]<<
                        " Size: "<<cutCellsMinusAndPlus[oldCellToPlusCutCell[owner[i]]].size()<<endl;
                    FatalErrorInFunction
                    << " Unsplit face interior can not be owner and neighboring a cut cell."
                    << exit(FatalError);
                }
                */
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

                if(cellsToSide_[neighbour[i]] == 0)
                {
                    if(cutCellsMinusAndPlus[oldCellToMinusCutCell[neighbour[i]]].size() == 0)
                    {
                        FatalErrorInFunction
                        << " Unsplit face interior inserted but cell has no cut face."
                        << exit(FatalError);
                    }
                    cutCellsMinusAndPlus[oldCellToMinusCutCell[neighbour[i]]].append
                    (
                        splitAndUnsplitFacesInterior.size()-1+addedCutFaces.size()
                    );
                }
                else if(cellsToSide_[owner[i]] == 0)
                {
                    if(cutCellsMinusAndPlus[oldCellToMinusCutCell[owner[i]]].size() == 0)
                    {
                        FatalErrorInFunction
                        << " Unsplit face interior inserted but cell has no cut face."
                        << exit(FatalError);
                    }
                    cutCellsMinusAndPlus[oldCellToMinusCutCell[owner[i]]].append
                    (
                        splitAndUnsplitFacesInterior.size()-1+addedCutFaces.size()
                    );
                }
                
                /*
                 * Error is wrong
                if( (oldCellToMinusCutCell[neighbour[i]] != -1 &&
                     oldCellToMinusCutCell[owner[i]] != -1)
                    ||
                    (cellsToSide_[neighbour[i]] == 0 && cellsToSide_[owner[i]] == 0)
                    )                    
                {
                    Info<<"Neighbor "<<neighbour[i]<<" Index: "<<oldCellToMinusCutCell[neighbour[i]]<<
                        " Size: "<<cutCellsMinusAndPlus[oldCellToMinusCutCell[neighbour[i]]].size()<<endl;
                    Info<<"Owner "<<owner[i]<<" Index: "<<oldCellToMinusCutCell[owner[i]]<<
                        " Size: "<<cutCellsMinusAndPlus[oldCellToMinusCutCell[owner[i]]].size()<<endl;
                    FatalErrorInFunction
                    << " Unsplit face interior can not be owner and neighboring a cut cell."
                    << exit(FatalError);
                }
                */
            }
            else
            {
                FatalErrorInFunction
                << "A face with the side: "<<facesToSide_[i]<<" was not treated."
                << " This must not happen."
                << exit(FatalError);
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
            << exit(FatalError);
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
            
            if(cutCellsMinusAndPlus[oldCellToPlusCutCell[owner[i]]].size() == 0)
            {
                FatalErrorInFunction
                << " Split face interior inserted but cell has no cut face."
                << exit(FatalError);
            }
            cutCellsMinusAndPlus[oldCellToPlusCutCell[owner[i]]].append
            (
                splitAndUnsplitFacesBoundary.size()-1+
                addedCutFaces.size()+splitAndUnsplitFacesInterior.size()
            );

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
            
            if(cutCellsMinusAndPlus[oldCellToMinusCutCell[owner[i]]].size() == 0)
            {
                FatalErrorInFunction
                << " Split face interior inserted but cell has no cut face."
                << exit(FatalError);
            }
            cutCellsMinusAndPlus[oldCellToMinusCutCell[owner[i]]].append
            (
                splitAndUnsplitFacesBoundary.size()-1+
                addedCutFaces.size()+splitAndUnsplitFacesInterior.size()
            );
            
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
                
                if(cellsToSide_[owner[i]] == 0)
                {
                    if(cutCellsMinusAndPlus[oldCellToPlusCutCell[owner[i]]].size() == 0)
                    {
                        FatalErrorInFunction
                        << " Unsplit face interior inserted but cell has no cut face."
                        << exit(FatalError);
                    }
                    cutCellsMinusAndPlus[oldCellToPlusCutCell[owner[i]]].append
                    (
                        splitAndUnsplitFacesBoundary.size()-1+
                        addedCutFaces.size()+splitAndUnsplitFacesInterior.size()
                    );
                }
            }
            else if(facesToSide_[i] == -1)
            {
                splitAndUnsplitFaceBoundaryNeighbor.append(-1);
                
                if(cellsToSide_[owner[i]] == 0)
                    splitAndUnsplitFaceBoundaryOwner.append(oldCellsToAddedMinusSideCellIndex[owner[i]]);
                else
                    splitAndUnsplitFaceBoundaryOwner.append(owner[i]);
                
                if(cellsToSide_[owner[i]] == 0)
                {
                    if(cutCellsMinusAndPlus[oldCellToMinusCutCell[owner[i]]].size() == 0)
                    {
                        FatalErrorInFunction
                        << " Unsplit face interior inserted but cell has no cut face."
                        << exit(FatalError);
                    }
                    cutCellsMinusAndPlus[oldCellToMinusCutCell[owner[i]]].append
                    (
                        splitAndUnsplitFacesBoundary.size()-1+
                        addedCutFaces.size()+splitAndUnsplitFacesInterior.size()
                    );
                }
            }
            else
            {
                FatalErrorInFunction
                << "A face with the side: "<<facesToSide_[i]<<" was not treated."
                << " This must not happen."
                << exit(FatalError);
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

void Foam::cutCellFvMesh::createNewMeshData_cutNeg
(
)
{
    const cellList& meshCells = this->cells();
    const faceList& meshFaces = this->faces();
    //const edgeList& meshEdges = this->edges();
    //const pointField& meshPoints = this->points();
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
        Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
    }
    
    oldSplittedCellToNewPlusCell = List<DynamicList<label>>(meshCells.size());
    oldSplittedCellToNewMinusCell = List<DynamicList<label>>(meshCells.size());
    for(int i=0;i<meshCells.size();i++)
    {
        oldSplittedCellToNewPlusCell[i] = -1;
        oldSplittedCellToNewMinusCell[i] = -1;
    }
    
    // Compute new cellIndexes for added cells
    labelList oldCellsToAddedMinusSideCellIndex(meshCells.size());
    deletedCellsList = labelList(meshCells.size());
    label addedCellIndex = 0;
    //label deletedCellNumber = 0;
    for(int i=0;i<meshCells.size();i++)
    {
        deletedCellsList[i] = 0;
        if(cellToFaces_[i].size() == 1 && cellToFaces_[i][0] >= nbrOfPrevFaces)
        {
            oldCellsToAddedMinusSideCellIndex[i] = addedCellIndex+oldCellsToAddedMinusSideCellIndex.size();
            oldSplittedCellToNewPlusCell[i] = i;
            oldSplittedCellToNewMinusCell[i]= oldCellsToAddedMinusSideCellIndex[i];
            //Info<<i<<"->"<<oldCellsToAddedMinusSideCellIndex[i]<<endl;
            addedCellIndex++;
        }
        if(cellsToSide_[i] == -1)
        {
            deletedCellsList[i] = 1;
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
    
    //Info<<"Insert Split cell faces"<<endl;
    // Compute List of new faces splitting old cells
    //label addedCutFacesNbr = 0;
    //addedCutFaces = DynamicList<face>();
    addedCutFaces.setCapacity(cellToFaces_.size());
    //addedCutFaceNeighbor = DynamicList<label>();
    addedCutFaceNeighbor.setCapacity(cellToFaces_.size());
    //addedCutFaceOwner = DynamicList<label>();
    addedCutFaceOwner.setCapacity(cellToFaces_.size());
    for(int i=0;i<cellToFaces_.size();i++)
    {
        if(cellToFaces_[i].size() == 1 && cellToFaces_[i][0] >= nbrOfPrevFaces)
        {
            face addedFace = newMeshFaces_[cellToFaces_[i][0]];
            
            labelList thisCellPointLabels = meshCells[i].labels(meshFaces);
            cell thisCell = meshCells[i];
            vector thisNormal = addedFace.normal(newMeshPoints_);
            //Info<<"This Normal: "<<thisNormal<<endl;
            point thisCentre = addedFace.centre(newMeshPoints_);
            //Info<<"This Centre: "<<thisCentre<<endl;
            
            label testInd = -1;
            for(int i=0;i<thisCellPointLabels.size();i++)
            {
                if(pointsToSide_[thisCellPointLabels[i]] == -1)
                {
                    testInd = thisCellPointLabels[i];
                    break;
                }
            }
            //Info<<"test Point:"<<newMeshPoints_[testInd]<<endl;
            vector centreToPointInd = newMeshPoints_[testInd] - thisCentre;
            //centreToPointInd -= thisCentre;
            //Info<<"centreToPointInd: "<<centreToPointInd<<endl;
            scalar dir = centreToPointInd && thisNormal;
            //Info<<"dir: "<<dir<<endl;
            if(dir < 0)
                addedFace = addedFace.reverseFace();
            
            //Info<<centreToPointInd<<endl;
            
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

    
    //Info<<"Insert split faces interior"<<endl;
    // Compute the List of new faces resulting from the splitting of old faces
    label addedSplitCellsInteriorNbr = 0;
    splitAndUnsplitFacesInterior.setCapacity(neighbour.size());
    splitAndUnsplitFaceInteriorNeighbor.setCapacity(neighbour.size());
    splitAndUnsplitFaceInteriorOwner.setCapacity(neighbour.size());
    bool addedOneFace;
    for(int i=0;i<neighbour.size();i++)
    {
        //Info<<"Face "<<i<<" size: "<<faceToEdges_[i].size()<<" on side: "<<facesToSide_[i]<<endl;
        
        addedOneFace = false;
        if(faceToEdges_[i].size() == 1 && faceToEdges_[i][0] >= nbrOfPrevEdges)
        {
            if(oldFacesToCutFaces_[i].size() != 2)
            {
                FatalErrorInFunction
                << " Splitted interior cell is cut into"<<oldFacesToCutFaces_[i].size()
                << " faces instead of the expected 2."
                << exit(FatalError);
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
            addedOneFace = true;
            
            addedSplitCellsInteriorNbr++;
        }
        else
        {
            //Info<<"GonetoElse"<<endl;
            // Interior uncut face on positive side is appended  without change
            if(facesToSide_[i] == 1)
            {
                splitAndUnsplitFacesInterior.append(meshFaces[i]);
                splitAndUnsplitFaceInteriorNeighbor.append(neighbour[i]);
                splitAndUnsplitFaceInteriorOwner.append(owner[i]);
                addedOneFace = true;
                
                addedSplitCellsInteriorNbr++;
                //Info<<"Inserted Split face"<<endl;
            }
            // Interior cell on that is neither +1 nor -1 must be 0 and be treated in the first if part
            else if(facesToSide_[i] != -1)
            {
                FatalErrorInFunction
                << "A face with the side: "<<facesToSide_[i]<<" was not treated."
                << " This must not happen."
                << exit(FatalError);
            }
            //Info<<"Jumped"<<endl;
        }
        //Info<<splitAndUnsplitFaceInteriorOwner.size()<<endl;
        if(addedOneFace && splitAndUnsplitFaceInteriorOwner[splitAndUnsplitFaceInteriorOwner.size()-1] == -1)
        {
            FatalErrorInFunction
            << " Owner of face must not be -1 as happend in face "<<i
            << exit(FatalError);
        }
    }

    for(int i=0;i<boundMesh.size();i++)
    {
        patchStarts[i] += addedSplitCellsInteriorNbr-neighbour.size();
    }
    
    //Info<<"Insert split faces boundary"<<endl;
    label currBoundaryPatch = 0;
    label countOldBoundaryFaces = 0;
    label countNewBoundaryFaces = 0;
    splitAndUnsplitFacesBoundary.setCapacity(meshFaces.size()-neighbour.size());
    splitAndUnsplitFaceBoundaryNeighbor.setCapacity(meshFaces.size()-neighbour.size());
    splitAndUnsplitFaceBoundaryOwner.setCapacity(meshFaces.size()-neighbour.size());
    for(int i=neighbour.size();i<meshFaces.size();i++)
    {
        //Info<<"Boundary face "<<i;
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
                << exit(FatalError);
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
        //Info<<i<<":"<<patchStarts[i-1]<<"+"<<patchSizes[i-1]<<"="<<patchStarts[i-1] + patchSizes[i-1]<<endl;
        patchStarts[i] = patchStarts[i-1] + patchSizes[i-1]; 
    }
    
    
    //reduce for empty cells
    labelList cellReductionNumb(meshCells.size());
    mapOldCellsToNewCells = List<DynamicList<label>>(meshCells.size());
    label count = 0;
    label newCells = 0;
    for(int i=0;i<cellReductionNumb.size();i++)
    {
        if(deletedCellsList[i] == 1)
        {
            count++;
            cellReductionNumb[i] = -1;
            mapOldCellsToNewCells[i] = -1;
        }
        else
        {
            cellReductionNumb[i] = count;
            mapOldCellsToNewCells[i] = newCells;
            newCells++;
        }
    }
    mapNewCellsToOldCells = labelList(meshCells.size()-cellReductionNumb.last());
    for(int i=0;i<mapOldCellsToNewCells.size();i++)
    {
        if(mapOldCellsToNewCells[i].size()!=0)
        {
            //Correction just to allow compiling
            mapNewCellsToOldCells[mapOldCellsToNewCells[i][0]] = i;
        }
    }
    
    
    for(int i=0;i<addedCutFaces.size();i++)
    {
        if(cellReductionNumb[addedCutFaceOwner[i]] != -1 &&
           cellReductionNumb[addedCutFaceNeighbor[i]] != -1)
        {
            addedCutFaceOwner[i] -= cellReductionNumb[addedCutFaceOwner[i]];
            addedCutFaceNeighbor[i] -= cellReductionNumb[addedCutFaceNeighbor[i]];
        }
        else
        {
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
    for(int i=0;i<splitAndUnsplitFacesInterior.size();i++)
    {
        if(cellReductionNumb[splitAndUnsplitFaceInteriorOwner[i]] != -1 &&
           cellReductionNumb[splitAndUnsplitFaceInteriorNeighbor[i]] != -1)
        {
            splitAndUnsplitFaceInteriorOwner[i] -= cellReductionNumb[splitAndUnsplitFaceInteriorOwner[i]];
            splitAndUnsplitFaceInteriorNeighbor[i] -= cellReductionNumb[splitAndUnsplitFaceInteriorNeighbor[i]];
        }
        else
        {
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
    for(int i=0;i<splitAndUnsplitFacesBoundary.size();i++)
    {
        if(cellReductionNumb[splitAndUnsplitFaceBoundaryOwner[i]] != -1)
        {
            splitAndUnsplitFaceBoundaryOwner[i] -= cellReductionNumb[splitAndUnsplitFaceBoundaryOwner[i]];
        }
        else
        {
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
    
    addedCutFaces.setCapacity(addedCutFaces.size());
    addedCutFaceOwner.setCapacity(addedCutFaceOwner.size());
    addedCutFaceNeighbor.setCapacity(addedCutFaceNeighbor.size());
    
    splitAndUnsplitFacesInterior.setCapacity(splitAndUnsplitFacesInterior.size());
    splitAndUnsplitFaceInteriorOwner.setCapacity(splitAndUnsplitFaceInteriorOwner.size());
    splitAndUnsplitFaceInteriorNeighbor.setCapacity(splitAndUnsplitFaceInteriorNeighbor.size());

    splitAndUnsplitFacesBoundary.setCapacity(splitAndUnsplitFacesBoundary.size());
    splitAndUnsplitFaceBoundaryOwner.setCapacity(splitAndUnsplitFaceBoundaryOwner.size());
    splitAndUnsplitFaceBoundaryNeighbor.setCapacity(splitAndUnsplitFaceBoundaryNeighbor.size());
}

void Foam::cutCellFvMesh::createNewMeshData_cutNeg_plus
(
)
{
    const cellList& meshCells = this->cells();
    const faceList& meshFaces = this->faces();
    //const edgeList& meshEdges = this->edges();
    //const pointField& meshPoints = this->points();
    const labelList owner   = this->faceOwner();
    const labelList neighbour = this->faceNeighbour();    
    const polyBoundaryMesh& boundMesh = this->boundaryMesh();
    const labelListList& pointToEdges = this->pointEdges();
    const labelListList& cellToEdges = this->cellEdges();
    const labelListList& cellToPoints = this->cellPoints();
    
    // Store old boundary patches
    patchStarts = labelList(boundMesh.size());
    patchSizes = labelList(boundMesh.size());
    Info<<endl;
    for(int i=0;i<boundMesh.size();i++)
    {
        patchStarts[i] = boundMesh[i].start();
        patchSizes[i] = boundMesh[i].faceCentres().size();
        Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
    }
    
    oldSplittedCellToNewPlusCell = List<DynamicList<label>>(meshCells.size());
    oldSplittedCellToNewMinusCell = List<DynamicList<label>>(meshCells.size());
    for(int i=0;i<meshCells.size();i++)
    {
        oldSplittedCellToNewPlusCell[i].setSize(0);
        oldSplittedCellToNewMinusCell[i].setSize(0);
    }
    
    // Compute new cellIndexes for added cells
    List<DynamicList<label>> oldCellsToAddedMinusSideCellIndex(meshCells.size());
    deletedCellsList = labelList(meshCells.size());
    label addedCellIndex = 0;
    DynamicList<DynamicList<DynamicList<label>>> cellToNewMinusCellsPointLabels;
    cellToNewMinusCellsPointLabels.setSize(meshCells.size());
    DynamicList<DynamicList<DynamicList<label>>> cellToNewPlusCellsPointLabels;
    cellToNewPlusCellsPointLabels.setSize(meshCells.size());
    DynamicList<DynamicList<label>> cellToNewMinusCellsIndexes;
    cellToNewMinusCellsIndexes.setSize(meshCells.size());
    DynamicList<DynamicList<label>> cellToNewPlusCellsIndexes;
    cellToNewPlusCellsIndexes.setSize(meshCells.size());
    
    //Info<<endl<<"cellToFaces_:"<<cellToFaces_<<endl;
    
    for(int i=0;i<meshCells.size();i++)
    {
        Info<<"---------------------------------------------------"<<endl;
        if(i==196879)
        {
            Info<<endl;
            Info<<"cellToFaces_[196879].size(): "<<cellToFaces_[i].size()<<endl;
            //FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
        }
        deletedCellsList[i] = 0;
        if(cellToFaces_[i].size() > 0)
        {
            std::unordered_set<label> pointsTreated;
            std::unordered_set<label> edgesTreated;
            DynamicList<DynamicList<label>> minusCells;
            DynamicList<DynamicList<label>> plusCells;
            const labelList& thisCellEdges = cellToEdges[i];
            const labelList& thisCellPoints = cellToPoints[i];
            oldSplittedCellToNewPlusCell[i] = i;
            for(int j=0;j<cellToFaces_[i].size();j++)
            {
                if(cellToFaces_[i][j] >= nbrOfPrevFaces)
                {
                    DynamicList<label> minusCell;
                    DynamicList<label> plusCell;
                    face cutFace = newMeshFaces_[cellToFaces_[i][j]];
                    DynamicList<DynamicList<edge>> facePointEdges;
                    facePointEdges.setSize(cutFace.size());
                    //Info<<"facePointEdges.size():"<<facePointEdges.size()<<endl;
                    //store the points around each cutFace cutted edge via facePointEdges
                    for(int k=0;k<cutFace.size();k++)
                    {
                        label pointInd = cutFace[k];
                        //Info<<"point: "<<pointInd<<"  ";
                        DynamicList<label> edgeInds;
                        if((pointInd<nbrOfPrevPoints && pointToEgde_[pointInd]!=-1)||
                           (pointToEgde_[pointInd]==-1 && pointInd>=nbrOfPrevPoints))
                            FatalErrorInFunction<<"Cut Point to edge assignment is wrong."<<endl;
                            
                        if(pointInd>=nbrOfPrevPoints)
                            edgeInds.append(pointToEgde_[pointInd]);
                        else
                        {
                            labelList pointEdges = pointToEdges[pointInd];
                            for(int l=0;l<pointEdges.size();l++)
                            {
                                labelList edgeCells = this->edgeCells(pointEdges[l]);
                                for(int m=0;m<edgeCells.size();m++)
                                {
                                    if(edgeCells[m] == i)
                                    {
                                        edgeInds.append(pointEdges[l]);
                                    }
                                }
                            }
                        }

                        //Info<<"edge: "<<edgeInds<<"  ";
                        
                        for(int l=0;l<edgeInds.size();l++)
                        {
                            if(edgesTreated.count(edgeInds[l]) == 0)
                            {
                                facePointEdges[k].append(newMeshEdges_[edgeInds[l]]);
                                edgesTreated.insert(edgeInds[l]);
                            }
                        }
                        //Info<<"facePointEdges[k]: "<<facePointEdges[k]<<endl;
                        pointsTreated.insert(pointInd);
                    }
                    
                    Info<<"pointsTreated: ";
                    for(auto& i : pointsTreated)
                        Info<<i<<" ";
                    Info<<endl;
                    
                    
                    //store the front points from the cut edges of each cut face
                    DynamicList<label> plusCellFrontPoints;
                    DynamicList<label> minusCellFrontPoints;
                    //Info<<"facePointEdges.size():"<<facePointEdges.size()<<endl;
                    for(int l=0;l<facePointEdges.size();l++)
                    {
                        for(int m=0;m<facePointEdges[l].size();m++)
                        {
                            labelList edgePoints(2);
                            edge oneEdge = facePointEdges[l][m];
                            edgePoints[0] = oneEdge.start();
                            edgePoints[1] = oneEdge.end();
                            for(int n=0;n<edgePoints.size();n++)
                            {
                                //Info<<"Try point: "<<edgePoints[n]<<endl;
                                if(pointsToSide_[edgePoints[n]] == 0)
                                {
                                    if(facePointEdges[l].size()==1)
                                        FatalErrorInFunction<<"An edge point can not be a zero Point if there is only one ."<<exit(FatalError);                                    
                                }
                                else if(pointsToSide_[edgePoints[n]] == 1)
                                {
                                    if(pointsTreated.count(edgePoints[n]) == 0)
                                    {
                                        plusCell.append(edgePoints[n]);
                                        pointsTreated.insert(edgePoints[n]);
                                        plusCellFrontPoints.append(edgePoints[n]);
                                    }
                                }
                                else if(pointsToSide_[edgePoints[n]] == -1)
                                {
                                    if(pointsTreated.count(edgePoints[n]) == 0)
                                    {
                                        minusCell.append(edgePoints[n]);
                                        pointsTreated.insert(edgePoints[n]);
                                        minusCellFrontPoints.append(edgePoints[n]);
                                    }
                                }
                                else
                                    FatalErrorInFunction<<"Point side must bei -1,0,1."<<exit(FatalError);                                    
                            }
                        }
                    }
                    
                    Info<<endl;
                    for(auto& i : pointsTreated)
                        Info<<i<<" ";
                    Info<<endl;
                    Info<<"plusCell  "<<plusCell<<endl;
                    Info<<"minusCell "<<minusCell<<endl;
                    Info<<"plusCellFrontPoints  "<<plusCellFrontPoints<<endl;
                    Info<<"minusCellFrontPoints "<<minusCellFrontPoints<<endl;
                    
                    
                    while(plusCellFrontPoints.size() > 0)
                    {
                        //Info<<"Step"<<endl;
                        DynamicList<label> temp = plusCellFrontPoints;
                        plusCellFrontPoints.setSize(0);
                        for(int l=0;l<temp.size();l++)
                        {
                            label frontPointInd = temp[l];
                            for(int m=0;m<thisCellEdges.size();m++)
                            {
                                if(edgesTreated.count(thisCellEdges[m])==0 && edgesToSide_[thisCellEdges[m]] == 1)
                                {
                                    label otherPoint = newMeshEdges_[thisCellEdges[m]].otherVertex(frontPointInd);
                                    //Info<<"startpoint:"<<frontPointInd<<" edge:"<<thisCellEdges[m]<<" otherPoint: "<<otherPoint;
                                    if(otherPoint!=-1 && pointsTreated.count(otherPoint)==0)
                                    {
                                        //Info<<" added";
                                        if(pointsToSide_[otherPoint] != 1)
                                            FatalErrorInFunction<<"Edge with 1 side has non 1 point"<<exit(FatalError);;                                        
                                        plusCell.append(otherPoint);
                                        plusCellFrontPoints.append(otherPoint);
                                        pointsTreated.insert(otherPoint);
                                        edgesTreated.insert(thisCellEdges[m]);
                                    }
                                    //Info<<endl;
                                }                                
                            }
                        }
                    }
                    while(minusCellFrontPoints.size() > 0)
                    {
                        //Info<<"Step"<<endl;
                        DynamicList<label> temp = minusCellFrontPoints;
                        minusCellFrontPoints.setSize(0);
                        for(int l=0;l<temp.size();l++)
                        {
                            label frontPointInd = temp[l];
                            for(int m=0;m<thisCellEdges.size();m++)
                            {
                                if(edgesTreated.count(thisCellEdges[m])==0 && edgesToSide_[thisCellEdges[m]] == -1)
                                {
                                    label otherPoint = newMeshEdges_[thisCellEdges[m]].otherVertex(frontPointInd);
                                    //Info<<"startpoint:"<<frontPointInd<<" edge:"<<thisCellEdges[m]<<" otherPoint: "<<otherPoint;
                                    if(otherPoint!=-1  && pointsTreated.count(otherPoint)==0)
                                    {
                                        //Info<<" added";
                                        if(pointsToSide_[otherPoint] != -1)
                                            FatalErrorInFunction<<"Edge with -1 side has non -1 point"<<exit(FatalError);                                     
                                        minusCell.append(otherPoint);
                                        minusCellFrontPoints.append(otherPoint);
                                        pointsTreated.insert(otherPoint);
                                        edgesTreated.insert(thisCellEdges[m]);
                                    }
                                    //Info<<endl;
                                }                                
                            }
                        }
                    }

                    if(plusCell.size()==0 || minusCell.size()==0)
                    {
                        Info<<endl<<endl;
                        Info<<"i:"<<i<<endl;
                        Info<<"j:"<<j<<endl;
                        Info<<"facePointEdges:"<<facePointEdges<<endl;
                        Info<<"cutFace:"<<cutFace<<endl;
                        Info<<"cell:"<<meshCells[i]<<endl;
                        Info<<"cell Points: "<<meshCells[i].labels(meshFaces)<<endl;
                        Info<<"cellToFaces_:"<<cellToFaces_[i]<<endl;
                        Info<<"nbrNewFaces:"<<newMeshFaces_.size()<<endl;
                        Info<<"newMeshFace:"<<newMeshFaces_[cellToFaces_[i][0]]<<endl;
                        Info<<plusCellFrontPoints<<endl;
                        Info<<minusCellFrontPoints<<endl;
                        Info<<plusCell<<endl;
                        Info<<minusCell<<endl;
                        FatalErrorInFunction<<"Face does not create at least one cell."<<exit(FatalError);
                    }
                    minusCells.append(minusCell);
                    plusCells.append(plusCell);                                      
                }
            }
            for(int j=0;j<thisCellPoints.size();j++)
            {
                if(pointsTreated.count(thisCellPoints[j])==0)
                    FatalErrorInFunction<<"Untreated point remains"<<exit(FatalError);
            }

            label countOldPoints = 0;
            for(int j=0;j<plusCells.size();j++)
                countOldPoints+=plusCells[j].size();
            for(int j=0;j<minusCells.size();j++)
                countOldPoints+=minusCells[j].size();
            labelList cellVerticeLabels = meshCells[i].labels(meshFaces);
            for(int j=0;j<cellVerticeLabels.size();j++)
                if(pointsToSide_[cellVerticeLabels[j]] == 0)
                    countOldPoints++;         

            if(countOldPoints!=8)
            {
                Info<<endl<<endl;
                for(int j=0;j<cellToFaces_[i].size();j++)
                    Info<<"cutface: "<<newMeshFaces_[cellToFaces_[i][0]]<<endl;
                Info<<"cellVerticeLabels: "<<cellVerticeLabels<<endl;
                Info<<"plusCells: "<<plusCells<<endl;
                Info<<"minusCells: "<<minusCells<<endl;
                FatalErrorInFunction<<"Old points unequal 8 but original cells have to be a cube"<<exit(FatalError);
            }

            if(minusCells.size()==0 || plusCells.size()==0)
                FatalErrorInFunction<<"Zero plus or minus cells"<<exit(FatalError);
            if(minusCells.size()>1 && plusCells.size()>1)
                FatalErrorInFunction<<"More than one plus and minus cells"<<exit(FatalError);
            if(!(minusCells.size()==1 || plusCells.size()==1))
                FatalErrorInFunction<<"Not one plus or minus cells"<<exit(FatalError);

            label countNewSplitFaces=0;
            for(int j=0;j<cellToFaces_[i].size();j++)
                if(cellToFaces_[i][j] >= nbrOfPrevFaces)
                    countNewSplitFaces++;
            if((countNewSplitFaces+1)!=(minusCells.size()+plusCells.size()))
                FatalErrorInFunction<<"Number of plus and minus cells does not match the number of cut faces"<<exit(FatalError);
            
            cellToNewMinusCellsPointLabels[i] = minusCells;
            cellToNewPlusCellsPointLabels[i] = plusCells;
            
            if(i==196879)
            {
                Info<<endl;
                Info<<"minusCells: "<<minusCells<<endl;
                Info<<"plusCells: "<<plusCells<<endl;
                FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
            }
            
            if(minusCells.size()==1 && plusCells.size()>1)
            {
                for(int k=0;k<plusCells.size();k++)
                {
                    oldSplittedCellToNewPlusCell[i].append(addedCellIndex+meshCells.size());
                    cellToNewPlusCellsIndexes[i].append(addedCellIndex+meshCells.size());
                    addedCellIndex++;                    
                }
                oldSplittedCellToNewMinusCell[i].append(i);
                cellToNewMinusCellsIndexes[i].append(i);
                deletedCellsList[i] = 1;
            }
            else if(minusCells.size()>1 && plusCells.size()==1)
            {
                for(int k=0;k<minusCells.size();k++)
                {
                    oldSplittedCellToNewMinusCell[i].append(addedCellIndex+meshCells.size());
                    cellToNewMinusCellsIndexes[i].append(addedCellIndex+meshCells.size());
                    addedCellIndex++;                    
                }
                oldSplittedCellToNewPlusCell[i].append(i);
                cellToNewPlusCellsIndexes[i].append(i);
            }
            else if(minusCells.size()==1 && plusCells.size()==1)
            {
                oldSplittedCellToNewMinusCell[i].append(addedCellIndex+meshCells.size());
                cellToNewMinusCellsIndexes[i].append(addedCellIndex+meshCells.size());
                addedCellIndex++;
                oldSplittedCellToNewPlusCell[i].append(i);
                cellToNewPlusCellsIndexes[i].append(i);
            }
            else
                FatalErrorInFunction<<"This combination is not possible"<<exit(FatalError);            
        }
        if(cellsToSide_[i] == -1)
        {
            deletedCellsList[i] = 1;
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
    
    //Info<<"Insert Split cell faces"<<endl;
    // Compute List of new faces splitting old cells
    addedCutFaces.setCapacity(cellToFaces_.size());
    addedCutFaceNeighbor.setCapacity(cellToFaces_.size());
    addedCutFaceOwner.setCapacity(cellToFaces_.size());
    for(int i=0;i<cellToFaces_.size();i++)
    {
        for(int j=0;j<cellToFaces_[i].size();j++)
        {
            if(cellToFaces_[i][j] >= nbrOfPrevFaces)
            {
                face addedFace = newMeshFaces_[cellToFaces_[i][j]];
            
                labelList thisCellPointLabels = meshCells[i].labels(meshFaces);
                cell thisCell = meshCells[i];
                vector thisNormal = addedFace.normal(newMeshPoints_);
                //Info<<"This Normal: "<<thisNormal<<endl;
                point thisCentre = addedFace.centre(newMeshPoints_);
                //Info<<"This Centre: "<<thisCentre<<endl;
            
                label testInd = -1;
                for(int i=0;i<thisCellPointLabels.size();i++)
                {
                    if(pointsToSide_[thisCellPointLabels[i]] == -1)
                    {
                        testInd = thisCellPointLabels[i];
                        break;
                    }
                }
                //Info<<"test Point:"<<newMeshPoints_[testInd]<<endl;
                vector centreToPointInd = newMeshPoints_[testInd] - thisCentre;
                //centreToPointInd -= thisCentre;
                //Info<<"centreToPointInd: "<<centreToPointInd<<endl;
                scalar dir = centreToPointInd && thisNormal;
                //Info<<"dir: "<<dir<<endl;
                if(dir < 0)
                    addedFace = addedFace.reverseFace();
            
                //Info<<centreToPointInd<<endl;
                if(oldSplittedCellToNewPlusCell[i].size()==1 && oldSplittedCellToNewMinusCell[i].size()==1)
                {
                    addedCutFaces.append(addedFace);
                    addedCutFaceNeighbor.append(-1);
                    addedCutFaceOwner.append(i);
                }
                else if(oldSplittedCellToNewPlusCell[i].size()==1 && oldSplittedCellToNewMinusCell[i].size()>1)
                {
                    addedCutFaces.append(addedFace);
                    addedCutFaceNeighbor.append(-1);
                    addedCutFaceOwner.append(i);
                }
                else if(oldSplittedCellToNewPlusCell[i].size()>1 && oldSplittedCellToNewMinusCell[i].size()==1)
                {
                    addedCutFaces.append(addedFace);
                    addedCutFaceNeighbor.append(-1);
                    addedCutFaceOwner.append(oldSplittedCellToNewPlusCell[i][j]);
                }
                else
                    FatalErrorInFunction<<"This combination is not possible"<<exit(FatalError);
            }
        }
    }

    //Info<<endl<<"cellToNewPlusCellsIndexes:"<<cellToNewPlusCellsIndexes<<endl;
    //Info<<endl<<"cellToNewMinusCellsIndexes:"<<cellToNewMinusCellsIndexes<<endl;
    
    //Info<<"Insert split faces interior"<<endl;
    // Compute the List of new faces resulting from the splitting of old faces
    label addedSplitCellsInteriorNbr = 0;
    splitAndUnsplitFacesInterior.setCapacity(neighbour.size());
    splitAndUnsplitFaceInteriorNeighbor.setCapacity(neighbour.size());
    splitAndUnsplitFaceInteriorOwner.setCapacity(neighbour.size());
    bool addedOneFace;
    for(int i=0;i<neighbour.size();i++)
    {
        if(oldFacesToCutFaces_[i].size()==1 || oldFacesToCutFaces_[i].size()>4)
        {
            FatalErrorInFunction<< "Face cut in one or more than four cut faces."<< exit(FatalError);
        }
        if(i==585222)
        {
            Info<<endl;
            Info<<"i:"<<i<<endl;
            Info<<"oldFacesToCutFaces_["<<i<<"]:"<<oldFacesToCutFaces_[i]<<endl;
        }
        //Info<<"Face "<<i<<" size: "<<faceToEdges_[i].size()<<" on side: "<<facesToSide_[i]<<endl;

        for(int j=0;j<oldFacesToCutFaces_[i].size();j++)
        {
            face face1      = cutFaces_[oldFacesToCutFaces_[i][j]];
            label signFace1 = cutFacesToSide_[oldFacesToCutFaces_[i][j]];
            
            if(signFace1>=0)
            {
                splitAndUnsplitFacesInterior.append(face1);
                label oldOwnerCell = owner[i];
                label oldNeighbourCell = neighbour[i];
                
                DynamicList<DynamicList<label>> ownerNewMinusCellsPointLabels = cellToNewMinusCellsPointLabels[oldOwnerCell];
                DynamicList<DynamicList<label>> ownerNewPlusCellsPointLabels = cellToNewPlusCellsPointLabels[oldOwnerCell];
                DynamicList<DynamicList<label>> neighbourNewMinusCellsPointLabels = cellToNewMinusCellsPointLabels[oldNeighbourCell];
                DynamicList<DynamicList<label>> neighbourNewPlusCellsPointLabels = cellToNewPlusCellsPointLabels[oldNeighbourCell];
                   
                DynamicList<label> ownerNewMinusCellsIndex = cellToNewMinusCellsIndexes[oldOwnerCell];
                DynamicList<label> ownerNewPlusCellsIndex = cellToNewPlusCellsIndexes[oldOwnerCell];
                DynamicList<label> neighbourNewMinusCellsIndex = cellToNewMinusCellsIndexes[oldNeighbourCell];
                DynamicList<label> neighbourNewPlusCellsIndex = cellToNewPlusCellsIndexes[oldNeighbourCell];
                
                DynamicList<std::unordered_set<label>> ownerNewMinusCellsPointLabelsMap;
                ownerNewMinusCellsPointLabelsMap.setSize(ownerNewMinusCellsPointLabels.size());
                for(int k=0;k<ownerNewMinusCellsPointLabels.size();k++)
                    for(int l=0;l<ownerNewMinusCellsPointLabels[k].size();l++)
                        ownerNewMinusCellsPointLabelsMap[k].insert(ownerNewMinusCellsPointLabels[k][l]);
                    
                DynamicList<std::unordered_set<label>> ownerNewPlusCellsPointLabelsMap;
                ownerNewPlusCellsPointLabelsMap.setSize(ownerNewPlusCellsPointLabels.size());
                for(int k=0;k<ownerNewPlusCellsPointLabels.size();k++)
                    for(int l=0;l<ownerNewPlusCellsPointLabels[k].size();l++)
                        ownerNewPlusCellsPointLabelsMap[k].insert(ownerNewPlusCellsPointLabels[k][l]);

                DynamicList<std::unordered_set<label>> neighbourNewMinusCellsPointLabelsMap;
                neighbourNewMinusCellsPointLabelsMap.setSize(neighbourNewMinusCellsPointLabels.size());
                for(int k=0;k<neighbourNewMinusCellsPointLabels.size();k++)
                    for(int l=0;l<neighbourNewMinusCellsPointLabels[k].size();l++)
                        neighbourNewMinusCellsPointLabelsMap[k].insert(neighbourNewMinusCellsPointLabels[k][l]);
                    
                DynamicList<std::unordered_set<label>> neighbourNewPlusCellsPointLabelsMap;
                neighbourNewPlusCellsPointLabelsMap.setSize(neighbourNewPlusCellsPointLabels.size());
                for(int k=0;k<neighbourNewPlusCellsPointLabels.size();k++)
                    for(int l=0;l<neighbourNewPlusCellsPointLabels[k].size();l++)
                        neighbourNewPlusCellsPointLabelsMap[k].insert(neighbourNewPlusCellsPointLabels[k][l]);
                    
                label newOwnerCell = -1;
                for(int k=0;k<ownerNewMinusCellsPointLabelsMap.size();k++)
                {
                    bool thisFaceNeighbors = false;
                    for(int l=0;l<face1.size();l++)
                    {
                        if(ownerNewMinusCellsPointLabelsMap[k].count(face1[l]))
                        {
                            if(newOwnerCell!=-1 && !thisFaceNeighbors)
                                FatalErrorInFunction<< "Old splitted face neigbours two newSplit Cells." << exit(FatalError);
                            newOwnerCell = ownerNewMinusCellsIndex[k];
                            thisFaceNeighbors = true;
                            Info<<"Assign: "<<i<<endl;
                        }
                    }
                }
                for(int k=0;k<ownerNewPlusCellsPointLabelsMap.size();k++)
                {
                    bool thisFaceNeighbors = false;
                    for(int l=0;l<face1.size();l++)
                    {
                        if(ownerNewPlusCellsPointLabelsMap[k].count(face1[l]))
                        {
                            if(newOwnerCell!=-1 && !thisFaceNeighbors)
                                FatalErrorInFunction<< "Old splitted face neigbours two newSplit Cells." << exit(FatalError);
                            newOwnerCell = ownerNewPlusCellsIndex[k];
                            thisFaceNeighbors = true;
                            Info<<"Assign: "<<i<<endl;
                        }
                    }
                }
                if(newOwnerCell==-1)
                {
                    Info<<endl;
                    Info<<"oldOwnerCell: "<<oldOwnerCell<<endl;
                    Info<<"oldNeighbourCell: "<<oldNeighbourCell<<endl;
                    Info<<"face: "<<i<<endl;
                    Info<<"oldFaceToCutFace: "<<j<<endl;
                    Info<<"ownerNewMinusCellsIndex: "<<ownerNewMinusCellsIndex<<endl;
                    Info<<"ownerNewPlusCellsIndex: "<<ownerNewPlusCellsIndex<<endl;
                    Info<<"ownerNewMinusCellsPointLabels: "<<ownerNewMinusCellsPointLabels<<endl;
                    Info<<"ownerNewPlusCellsPointLabels: "<<ownerNewPlusCellsPointLabels<<endl;
                    Info<<"faceInd: "<<i<<endl;
                    Info<<"face1: "<<face1<<endl;
                    FatalErrorInFunction<< "Old splitted face does not neigbour a newSplit Cell." << exit(FatalError);
                }
                label newNeighborCell = -1;
                for(int k=0;k<neighbourNewMinusCellsPointLabelsMap.size();k++)
                {
                    bool thisFaceNeighbors = false;
                    for(int l=0;l<face1.size();l++)
                    {
                        if(neighbourNewMinusCellsPointLabelsMap[k].count(face1[l]))
                        {
                            if(newNeighborCell!=-1 && !thisFaceNeighbors)
                                FatalErrorInFunction<< "Old splitted face neigbours two newSplit Cells." << exit(FatalError);
                            newNeighborCell = neighbourNewMinusCellsIndex[k];
                            thisFaceNeighbors = true;
                        }
                    }
                }
                for(int k=0;k<neighbourNewPlusCellsPointLabelsMap.size();k++)
                {
                    bool thisFaceNeighbors = false;
                    for(int l=0;l<face1.size();l++)
                    {
                        if(neighbourNewPlusCellsPointLabelsMap[k].count(face1[l]))
                        {
                            if(newNeighborCell!=-1 && !thisFaceNeighbors)
                                FatalErrorInFunction<< "Old splitted face neigbours two newSplit Cells." << exit(FatalError);
                            newNeighborCell = neighbourNewPlusCellsIndex[k];
                            thisFaceNeighbors = true;
                        }
                    }
                }
                if((newNeighborCell==-1 && signFace1 != 0)||(signFace1==0 && newNeighborCell!=-1))
                {
                    FatalErrorInFunction<< "Zero face must have no neighbor" << exit(FatalError);
                }
                if(signFace1==0)
                {
                    splitAndUnsplitFaceInteriorNeighbor.append(-1);
                }
                else
                {
                    splitAndUnsplitFaceInteriorNeighbor.append(newNeighborCell);
                    addedSplitCellsInteriorNbr++;
                }

                splitAndUnsplitFaceInteriorOwner.append(newOwnerCell);

                addedOneFace = true;
            }
        }
        if(oldFacesToCutFaces_[i].size()==0)
        {
            //Info<<"GonetoElse"<<endl;
            // Interior uncut face on positive side is appended  without change
            if(facesToSide_[i] == 1)
            {
                splitAndUnsplitFacesInterior.append(meshFaces[i]);
                splitAndUnsplitFaceInteriorNeighbor.append(neighbour[i]);
                splitAndUnsplitFaceInteriorOwner.append(owner[i]);
                addedOneFace = true;
                
                addedSplitCellsInteriorNbr++;
                //Info<<"Inserted Split face"<<endl;
            }
            else if(facesToSide_[i] == 0)
            {
                splitAndUnsplitFacesInterior.append(meshFaces[i]);
                splitAndUnsplitFaceInteriorNeighbor.append(-1);
                splitAndUnsplitFaceInteriorOwner.append(owner[i]);
                addedOneFace = true;
                
            }
            // Interior cell on that is neither +1 nor -1 must be 0 and be treated in the first if part
            else if(facesToSide_[i] != -1)
            {
                FatalErrorInFunction
                << "A face with the side: "<<facesToSide_[i]<<" was not treated."
                << " This must not happen."
                << exit(FatalError);
            }
        }
        //Info<<splitAndUnsplitFaceInteriorOwner.size()<<endl;
        if(addedOneFace && splitAndUnsplitFaceInteriorOwner[splitAndUnsplitFaceInteriorOwner.size()-1] == -1)
        {
            FatalErrorInFunction
            << " Owner of face must not be -1 as happend in face "<<i
            << exit(FatalError);
        }
    }
    
    //Info<<endl<<"splitAndUnsplitFaceInteriorOwner"<<splitAndUnsplitFaceInteriorOwner<<endl;
    //Info<<endl<<"splitAndUnsplitFaceInteriorNeighbor"<<splitAndUnsplitFaceInteriorNeighbor<<endl;
    
    Info<<"New Boundary"<<endl;
    for(int i=0;i<boundMesh.size();i++)
    {
        patchStarts[i] += addedSplitCellsInteriorNbr-neighbour.size();
        Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
    }
    
    //Info<<"Insert split faces boundary"<<endl;
    label currBoundaryPatch = 0;
    label countOldBoundaryFaces = 0;
    label countNewBoundaryFaces = 0;
    splitAndUnsplitFacesBoundary.setCapacity(meshFaces.size()-neighbour.size());
    splitAndUnsplitFaceBoundaryNeighbor.setCapacity(meshFaces.size()-neighbour.size());
    splitAndUnsplitFaceBoundaryOwner.setCapacity(meshFaces.size()-neighbour.size());
    for(int i=neighbour.size();i<meshFaces.size();i++)
    {
        if(oldFacesToCutFaces_[i].size()==1 || oldFacesToCutFaces_[i].size()>4)
        {
            FatalErrorInFunction<< "Face cut in one or more than four cut faces."<< exit(FatalError);
        }
        
        for(int j=0;j<oldFacesToCutFaces_[i].size();j++)
        {
            face face1      = cutFaces_[oldFacesToCutFaces_[i][j]];
            label signFace1 = cutFacesToSide_[oldFacesToCutFaces_[i][j]];
            
            if(signFace1>=0)
            {
                splitAndUnsplitFacesBoundary.append(face1);
                label oldOwnerCell = owner[i];
                
                DynamicList<DynamicList<label>> ownerNewMinusCellsPointLabels = cellToNewMinusCellsPointLabels[oldOwnerCell];
                DynamicList<DynamicList<label>> ownerNewPlusCellsPointLabels = cellToNewPlusCellsPointLabels[oldOwnerCell];
                   
                DynamicList<label> ownerNewMinusCellsIndex = cellToNewMinusCellsIndexes[oldOwnerCell];
                DynamicList<label> ownerNewPlusCellsIndex = cellToNewPlusCellsIndexes[oldOwnerCell];
                
                DynamicList<std::unordered_set<label>> ownerNewMinusCellsPointLabelsMap;
                ownerNewMinusCellsPointLabelsMap.setSize(ownerNewMinusCellsPointLabels.size());
                for(int k=0;k<ownerNewMinusCellsPointLabels.size();k++)
                    for(int l=0;l<ownerNewMinusCellsPointLabels[k].size();l++)
                        ownerNewMinusCellsPointLabelsMap[k].insert(ownerNewMinusCellsPointLabels[k][l]);
                    
                DynamicList<std::unordered_set<label>> ownerNewPlusCellsPointLabelsMap;
                ownerNewPlusCellsPointLabelsMap.setSize(ownerNewPlusCellsPointLabels.size());
                for(int k=0;k<ownerNewPlusCellsPointLabels.size();k++)
                    for(int l=0;l<ownerNewPlusCellsPointLabels[k].size();l++)
                        ownerNewPlusCellsPointLabelsMap[k].insert(ownerNewPlusCellsPointLabels[k][l]);
                  
                label newOwnerCell = -1;
                for(int k=0;k<ownerNewMinusCellsPointLabelsMap.size();k++)
                {
                    bool thisFaceNeighbors = false;
                    for(int l=0;l<face1.size();l++)
                    {
                        if(ownerNewMinusCellsPointLabelsMap[k].count(face1[l]))
                        {
                            if(newOwnerCell!=-1 && !thisFaceNeighbors)
                                FatalErrorInFunction<< "Old splitted face neigbours two newSplit Cells." << exit(FatalError);
                            newOwnerCell = ownerNewMinusCellsIndex[k];
                            thisFaceNeighbors = true;
                        }
                    }
                }
                for(int k=0;k<ownerNewPlusCellsPointLabelsMap.size();k++)
                {
                    bool thisFaceNeighbors = false;
                    for(int l=0;l<face1.size();l++)
                    {
                        if(ownerNewPlusCellsPointLabelsMap[k].count(face1[l]))
                        {
                            if(newOwnerCell!=-1 && !thisFaceNeighbors)
                                FatalErrorInFunction<< "Old splitted face neigbours two newSplit Cells." << exit(FatalError);
                            newOwnerCell = ownerNewPlusCellsIndex[k];
                            thisFaceNeighbors = true;
                        }
                    }
                }
                if(newOwnerCell==-1)
                    FatalErrorInFunction<< "Old splitted face does not neigbour a newSplit Cell." << exit(FatalError);
                
                splitAndUnsplitFaceBoundaryNeighbor.append(-1);
                splitAndUnsplitFaceBoundaryOwner.append(newOwnerCell);

                addedOneFace = true;
                countNewBoundaryFaces++;
            }
        }
        //Info<<"Boundary face "<<i;
        if(oldFacesToCutFaces_[i].size()==0)
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
                countNewBoundaryFaces++;
                splitAndUnsplitFacesBoundary.append(meshFaces[i]);
                splitAndUnsplitFaceBoundaryNeighbor.append(-1);
                splitAndUnsplitFaceBoundaryOwner.append(owner[i]);
            }
            else if(facesToSide_[i] != -1)
            {
                FatalErrorInFunction
                << "A face with the side: "<<facesToSide_[i]<<" was not treated."
                << " This must not happen."
                << exit(FatalError);
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
            Info<<"NewPatchSize["<<currBoundaryPatch<<"]: "<<patchSizes[currBoundaryPatch]<<endl;
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
    
    
    //reduce for empty cells
    DynamicList<label> cellReductionNumb;
    cellReductionNumb.setSize(meshCells.size());
    mapOldCellsToNewCells.setSize(meshCells.size());
    label count = 0;
    label countDel = 0;
    //Info<<endl<<"cellReductionNumb.size(): "<<cellReductionNumb.size()<<endl;
    //Info<<endl<<"deletedCellsList: "<<deletedCellsList<<endl;
    for(int i=0;i<cellReductionNumb.size();i++)
    {
        if(deletedCellsList[i] == 1)
        {
            if(!(cellToNewPlusCellsIndexes[i].size()>1 && cellToNewMinusCellsIndexes[i].size()==1) &&
               !(cellToNewPlusCellsIndexes[i].size()==0 && cellToNewMinusCellsIndexes[i].size()==0))
            {
                Info<<endl<<"cellToNewPlusCellsIndexes["<<i<<"].size():"<<cellToNewPlusCellsIndexes[i].size()<<endl;
                Info<<endl<<"cellToNewMinusCellsIndexes["<<i<<"].size():"<<cellToNewMinusCellsIndexes[i].size()<<endl;
                FatalErrorInFunction<<"Deleted Cell. Invalid cell splitting."<<exit(FatalError);
            }

            if(cellToNewPlusCellsIndexes[i].size()==0)
            {            
                mapOldCellsToNewCells[i].setSize(0);
            }
            else
            {
                if(!(cellToNewPlusCellsIndexes[i].size()>1))
                    FatalErrorInFunction<<"Face neighbors or ownes deleted cell. This can not happen."<<exit(FatalError);
                    
                for(int j=0;j<cellToNewPlusCellsIndexes[i].size();j++)
                {
                    mapOldCellsToNewCells[i].append(cellToNewPlusCellsIndexes[i][j]);
                }
            }
            cellReductionNumb[i] = -1;
            countDel++;
        }
        else
        {
            if(!(cellToNewPlusCellsIndexes[i].size()==1 && cellToNewMinusCellsIndexes[i].size()==1) &&
               !(cellToNewPlusCellsIndexes[i].size()==1 && cellToNewMinusCellsIndexes[i].size()>1) &&
               !(cellToNewPlusCellsIndexes[i].size()==0 && cellToNewMinusCellsIndexes[i].size()==0))
            {
                Info<<endl<<"cellToNewPlusCellsIndexes["<<i<<"].size():"<<cellToNewPlusCellsIndexes[i].size()<<endl;
                Info<<endl<<"cellToNewMinusCellsIndexes["<<i<<"].size():"<<cellToNewMinusCellsIndexes[i].size()<<endl;
                FatalErrorInFunction<<"Nondeleted Cell. Invalid cell splitting."<<exit(FatalError);
            }
                
            cellReductionNumb[i] = countDel;
            mapOldCellsToNewCells[i].append(count);
        }
        count++;
    }
    //Info<<endl<<"cellReductionNumb: "<<cellReductionNumb<<endl;
    //Info<<endl<<"mapOldCellsToNewCells: "<<mapOldCellsToNewCells<<endl;


    label maxNumOfNewCells = 0;
    for(int i=0;i<mapOldCellsToNewCells.size();i++)
    {
        if((mapOldCellsToNewCells[i].size()==0 && cellReductionNumb[i]!=-1)||(mapOldCellsToNewCells[i].size()!=0 && cellReductionNumb[i]==-1))
            FatalErrorInFunction<<"Face neighbors or ownes deleted cell. This can not happen."<<exit(FatalError);
        for(int j=0;j<mapOldCellsToNewCells[i].size();j++)
        {
            if(mapOldCellsToNewCells[i][j]<i)
                FatalErrorInFunction<<"Old cells can only be mapped to new cells of higher index."<<exit(FatalError);
            mapOldCellsToNewCells[i][j] -= cellReductionNumb[i];
            if(maxNumOfNewCells<mapOldCellsToNewCells[i][j])
                maxNumOfNewCells = mapOldCellsToNewCells[i][j];
            
        }
    }
    maxNumOfNewCells++;
    //Info<<endl<<"mapOldCellsToNewCells: "<<mapOldCellsToNewCells<<endl;
    //Info<<endl<<"maxNumOfNewCells: "<<maxNumOfNewCells<<endl;
    
    for(int i=cellReductionNumb.size();i<maxNumOfNewCells;i++)
    {
        cellReductionNumb.append(countDel);
    }
        
    mapNewCellsToOldCells = labelList(maxNumOfNewCells,-1);
    for(int i=0;i<mapOldCellsToNewCells.size();i++)
    {
        for(int j=0;j<mapOldCellsToNewCells[i].size();j++)
        {
            if(mapNewCellsToOldCells[mapOldCellsToNewCells[i][j]]!=-1)
            {
                Info<<endl<<"mapOldCellsToNewCells["<<i<<"]["<<j<<"]:"<<mapOldCellsToNewCells[i][j]<<endl;
                FatalErrorInFunction<<"Multiple assign to mapNewCellsToOldCells."<<exit(FatalError);
            }
            mapNewCellsToOldCells[mapOldCellsToNewCells[i][j]] = i;                
        }
    }
    //Info<<endl<<"mapNewCellsToOldCells: "<<mapNewCellsToOldCells<<endl;
        
    for(int i=0;i<addedCutFaces.size();i++)
    {
        if(cellReductionNumb[addedCutFaceOwner[i]] != -1 &&
           cellReductionNumb[addedCutFaceNeighbor[i]] != -1)
        {
            addedCutFaceOwner[i] -= cellReductionNumb[addedCutFaceOwner[i]];
        }
        else
        {
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }

    //Info<<endl<<"splitAndUnsplitFaceInteriorNeighbor 1"<<splitAndUnsplitFaceInteriorNeighbor<<endl;

    
    for(int i=0;i<splitAndUnsplitFacesInterior.size();i++)
    {
        if(cellReductionNumb[splitAndUnsplitFaceInteriorOwner[i]] != -1 &&
           cellReductionNumb[splitAndUnsplitFaceInteriorNeighbor[i]] != -1)
        {
            splitAndUnsplitFaceInteriorOwner[i] -= cellReductionNumb[splitAndUnsplitFaceInteriorOwner[i]];
            if(splitAndUnsplitFaceInteriorNeighbor[i]!=-1)
                splitAndUnsplitFaceInteriorNeighbor[i] -= cellReductionNumb[splitAndUnsplitFaceInteriorNeighbor[i]];
        }
        else
        {
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
    
    //Info<<endl<<"splitAndUnsplitFaceInteriorNeighbor 2"<<splitAndUnsplitFaceInteriorNeighbor<<endl;
    
    
    DynamicList<face> tempFaceNonNeighb;
    DynamicList<label> tempOwnerNonNeighb;
    DynamicList<label> tempNeighbourNonNeighb;
    DynamicList<face> tempFaceWithNeighb;
    DynamicList<label> tempOwnerWithNeighb;
    DynamicList<label> tempNeighbourWithNeighb;
    for(int i=0;i<splitAndUnsplitFacesInterior.size();i++)
    {
        if(splitAndUnsplitFaceInteriorNeighbor[i] == -1)
        {
            tempFaceNonNeighb.append(splitAndUnsplitFacesInterior[i]);
            tempOwnerNonNeighb.append(splitAndUnsplitFaceInteriorOwner[i]);
            tempNeighbourNonNeighb.append(splitAndUnsplitFaceInteriorNeighbor[i]);
        }
        else
        {
            tempFaceWithNeighb.append(splitAndUnsplitFacesInterior[i]);
            tempOwnerWithNeighb.append(splitAndUnsplitFaceInteriorOwner[i]);
            tempNeighbourWithNeighb.append(splitAndUnsplitFaceInteriorNeighbor[i]);
        }
    }
    
    //Info<<endl;
    //Info<<"tempNeighbourWithNeighb:"<<tempNeighbourWithNeighb<<endl;
    //Info<<"tempNeighbourNonNeighb:"<<tempNeighbourNonNeighb<<endl;
    
    splitAndUnsplitFacesInterior.setSize(0);
    splitAndUnsplitFaceInteriorOwner.setSize(0);
    splitAndUnsplitFaceInteriorNeighbor.setSize(0);
    splitAndUnsplitFacesInterior.append(tempFaceWithNeighb);
    splitAndUnsplitFaceInteriorOwner.append(tempOwnerWithNeighb);
    splitAndUnsplitFaceInteriorNeighbor.append(tempNeighbourWithNeighb);
    splitAndUnsplitFacesInteriorToBoundary.append(tempFaceNonNeighb);
    splitAndUnsplitFaceInteriorToBoundaryOwner.append(tempOwnerNonNeighb);
    splitAndUnsplitFaceInteriorToBoundaryNeighbor.append(tempNeighbourNonNeighb);
    
    for(int i=0;i<patchStarts.size()-1;i++)
        patchStarts[i] -= splitAndUnsplitFacesInteriorToBoundary.size();
    
    for(int i=0;i<splitAndUnsplitFacesBoundary.size();i++)
    {
        if(cellReductionNumb[splitAndUnsplitFaceBoundaryOwner[i]] != -1)
        {
            splitAndUnsplitFaceBoundaryOwner[i] -= cellReductionNumb[splitAndUnsplitFaceBoundaryOwner[i]];
        }
        else
        {
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
    
    addedCutFaces.setCapacity(addedCutFaces.size());
    addedCutFaceOwner.setCapacity(addedCutFaceOwner.size());
    addedCutFaceNeighbor.setCapacity(addedCutFaceNeighbor.size());
    
    splitAndUnsplitFacesInterior.setCapacity(splitAndUnsplitFacesInterior.size());
    splitAndUnsplitFaceInteriorOwner.setCapacity(splitAndUnsplitFaceInteriorOwner.size());
    splitAndUnsplitFaceInteriorNeighbor.setCapacity(splitAndUnsplitFaceInteriorNeighbor.size());
    
    //Info<<endl<<"splitAndUnsplitFaceInteriorNeighbor"<<splitAndUnsplitFaceInteriorNeighbor<<endl;

    splitAndUnsplitFacesBoundary.setCapacity(splitAndUnsplitFacesBoundary.size());
    splitAndUnsplitFaceBoundaryOwner.setCapacity(splitAndUnsplitFaceBoundaryOwner.size());
    splitAndUnsplitFaceBoundaryNeighbor.setCapacity(splitAndUnsplitFaceBoundaryNeighbor.size());
}

void Foam::cutCellFvMesh::printNewMeshData
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

void Foam::cutCellFvMesh::printMesh
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
        Info<<"Cell:"<<i<<" with "<<meshCells[i].nFaces()<<" faces |";
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


void Foam::cutCellFvMesh::selfTestMesh()
{  
    Info<<"START MESH SELF TEST"<<endl;
    const cellList& meshCells = this->cells();
    const faceList& meshFaces = this->faces();
    //const edgeList& meshEdges = this->edges();
    const pointField& meshPoints = this->points();
    const labelList owner   = this->faceOwner();
    const labelList neighbour = this->faceNeighbour();    
    //const polyBoundaryMesh& boundMesh = this->boundaryMesh();
    
    scalar maxFaceSize = -1,avgFaceSize=0,minFaceSize;
    if(meshFaces.size()>0) minFaceSize = meshFaces[0].mag(meshPoints);
    for(int i=0;i<meshFaces.size();i++)
    {
        scalar thisFaceSize = meshFaces[i].mag(meshPoints);
        avgFaceSize+=thisFaceSize;
        if(maxFaceSize<thisFaceSize)
            maxFaceSize=thisFaceSize;
        if(minFaceSize>thisFaceSize)
            minFaceSize=thisFaceSize;
    }
    avgFaceSize /= meshFaces.size();
    
    for(int i=0;i<meshCells.size();i++)
    {
        if(meshCells[i].size() < 4)
            FatalErrorInFunction<<"Cell with less than four faces!"<< exit(FatalError);
        for(int k=0;k<meshCells[i].size();k++)
        {
            if(meshCells[i][k]<0 || meshCells[i][k]>=meshFaces.size())
                FatalErrorInFunction<<"Cell with faces out of bound!"<< exit(FatalError);
        }

    }
    
    for(int i=0;i<meshFaces.size();i++)
    {
        if(meshFaces[i].size() < 3)
            FatalErrorInFunction<<"Face with less than three vertices!"<< exit(FatalError);
        for(int k=0;k<meshFaces[i].size();k++)
        {
            if(meshFaces[i][k]<0 || meshFaces[i][k]>=meshPoints.size())
                FatalErrorInFunction<<"Face with points out of bound!"<< exit(FatalError);
        }
    }
    
    scalar maxCellSize = -1,avgCellSize=0,minCellSize;
    if(meshCells.size()>0) minCellSize = meshCells[0].mag(meshPoints,meshFaces);
    for(int i=0;i<meshCells.size();i++)
    {
        scalar thisCellSize = meshCells[i].mag(meshPoints,meshFaces);
        avgCellSize+=thisCellSize;
        if(maxCellSize<thisCellSize)
            maxCellSize=thisCellSize;
        if(minCellSize>thisCellSize)
            minCellSize=thisCellSize;
    }
    avgCellSize /= meshCells.size();
    
    label countExtrSmall = 0;
    for(int i=0;i<meshCells.size();i++)
    {
        scalar thisCellSize = meshCells[i].mag(meshPoints,meshFaces);
        if(thisCellSize < (maxCellSize*partialThreeshold*(1.f/1e10)))
        {
            countExtrSmall++;
            Info<<"i:"<<i<<"  size:"<<thisCellSize<<endl;
        }
    }
    Info<<"countExtrSmall: "<<countExtrSmall<<endl;

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
        List<scalar> crossProdsArea(meshFaces[i].size()+1);
        for(int k=0;k<meshFaces[i].size();k++)
        {
            prev = meshPoints[meshFaces[i].prevLabel(k)];
            curr = meshPoints[meshFaces[i][k]];
            next = meshPoints[meshFaces[i].nextLabel(k)];
            
            edge1 = curr-prev;
            edge2 = next-curr;
            crossProds[k] = crossProd(edge1,edge2);
            crossProdsArea[k] = norm2(crossProds[k]);
        }
        crossProds[crossProds.size()-1] = crossProds[0];
        crossProdsArea[crossProds.size()-1] = norm2(crossProds[0]);
        scalar res;
        for(int k=0;k<meshFaces[i].size();k++)
        {
            res = crossProds[k] && crossProds[k+1];
            if(res<0 && area>=partialThreeshold*maxFaceSize*(1.f/1e10) && 
               crossProdsArea[k] >=partialThreeshold*maxFaceSize*(1.f/1e10) &&
               crossProdsArea[k+1] >=partialThreeshold*maxFaceSize*(1.f/1e10)
              )
            {
                Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
                if(i < neighbour.size())
                    Info<<" Neighbor:"<<neighbour[i]<<" ";
                Info<<endl;
                for(int k=0;k<meshFaces[i].size();k++)
                {
                    Info<<meshFaces[i][k]<<meshPoints[meshFaces[i][k]]<<"->";
                }
                Info<<endl;
                Info<<" with centre:"<<centreFace;
                Info<<" and normal vector:"<<normalFace;
                Info<<" and area:"<<area<<endl;
                Info<<"crossProds["<<k<<"]:"<<crossProds[k]<<endl;
                Info<<"crossProds["<<k+1<<"]:"<<crossProds[k+1]<<endl;
                Info<<"crossProdsArea["<<k<<"]:"<<crossProdsArea[k]<<endl;
                Info<<"crossProdsArea["<<k+1<<"]:"<<crossProdsArea[k+1]<<endl;
                Info<<"maxArea: "<<partialThreeshold*maxFaceSize*(1.f/1e10)<<endl;
                
                FatalErrorInFunction
                << "Face must not have a concave shape!"
                << exit(FatalError);
            }
        }
        
        //Test for centre point internal
        vector toCentre1,toCentre2;
        crossProds = List<vector>(meshFaces[i].size()+1);
        crossProdsArea = List<scalar>(meshFaces[i].size()+1);
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
            crossProdsArea[k] = norm2(crossProds[k]);
            //Info<<"crossProds: "<<crossProds[k]<<endl;

        }
        crossProds[crossProds.size()-1] = crossProds[0];
        crossProdsArea[crossProds.size()-1] = norm2(crossProds[0]);
        for(int k=0;k<meshFaces[i].size();k++)
        {
            res = crossProds[k] && crossProds[k+1];
            if(res<0 && area>=partialThreeshold*maxFaceSize*(1.f/1e10) && 
               crossProdsArea[k] >=partialThreeshold*maxFaceSize*(1.f/1e10) &&
               crossProdsArea[k+1] >=partialThreeshold*maxFaceSize*(1.f/1e10)
              )
            {
                Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
                if(i < neighbour.size())
                    Info<<" Neighbor:"<<neighbour[i]<<" ";
                Info<<endl;
                for(int k=0;k<meshFaces[i].size();k++)
                {
                    Info<<meshFaces[i][k]<<meshPoints[meshFaces[i][k]]<<"->";
                }
                Info<<endl;
                Info<<" with centre:"<<centreFace;
                Info<<" and normal vector:"<<normalFace;
                Info<<" and area:"<<area<<endl;
                
                Info<<"CrossProds k: "<<crossProds[k]<<endl;
                Info<<"CrossProds k+1: "<<crossProds[k+1]<<endl;
                Info<<"res: "<<res<<endl;
                
                
                FatalErrorInFunction
                << "Face  has a centre thats not strictly inside!"
                << exit(FatalError);
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
                    Info<<endl;
                    for(int s=0;s<meshFaces[i].size();s++)
                    {
                        Info<<meshFaces[i][s]<<meshPoints[meshFaces[i][s]]<<"->";
                    }
                    Info<<endl;
                    Info<<" with centre:"<<centreFace;
                    Info<<" and normal vector:"<<normalFace;
                    Info<<" and area:"<<area<<endl;
            
                    FatalErrorInFunction
                    << "Face had a double point!"
                    << exit(FatalError);
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
            << exit(FatalError); 
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
            << exit(FatalError); 
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
                << exit(FatalError); 
            }
        }
        
        
        vector faceCentreToOwnerCentre = meshCells[ownerCell].centre(meshPoints,meshFaces)-centreFace;
        if((faceCentreToOwnerCentre && normalFace)>=0)
        {
            cell owneCell = meshCells[ownerCell];
            scalar thisCellSize = owneCell.mag(meshPoints,meshFaces);
            if(thisCellSize >= (maxCellSize*partialThreeshold*(1.f/1e10)) && norm2(normalFace)!=0 &&
               area>=partialThreeshold*maxFaceSize*(1.f/1e10))
            // too small cells might fail in this condition because of rounding error.
            // as a result they are exempt here.
            // the same is for zero faces with a zero normal vector
            {
                Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
                if(i < neighbour.size())
                    Info<<"Neighbor:"<<neighbour[i]<<" ";
                Info<<endl;
                Info<<"faceCentreToOwnerCentre: "<<faceCentreToOwnerCentre<<endl;
                Info<<"ownerCentre: "<<meshCells[ownerCell].centre(meshPoints,meshFaces)<<endl;
                Info<<"normalFace: "<<normalFace<<endl;
                Info<<endl;
                for(int k=0;k<meshFaces[i].size();k++)
                {
                    Info<<meshFaces[i][k]<<meshPoints[meshFaces[i][k]]<<"->";
                }
                Info<<endl;
                for(int k=0;k<meshFaces[i].size();k++)
                {
                    Info<<meshPoints[meshFaces[i][(k+1)%meshFaces[i].size()]]-meshPoints[meshFaces[i][k%meshFaces[i].size()]]<<"->";
                }
                Info<<endl;
                Info<<" with centre:"<<centreFace;
                Info<<" and normal vector:"<<normalFace;
                Info<<" and area:"<<area<<endl<<endl;
                
                Info<<"Owner Cell"<<endl;
                cell oneCell = meshCells[ownerCell];
                Info<<"Cell centre is: "<<oneCell.centre(meshPoints,meshFaces)<<endl;
                Info<<"Cell: "<<oneCell<<endl;
                Info<<"Cell Size: "<<oneCell.size()<<endl;
                Info<<"Cell volume: "<<oneCell.mag(meshPoints,meshFaces)<<endl;
            
                for(int k=0;k<oneCell.size();k++)
                {
                    label oneFaceInd = oneCell[k];
                    face oneFace = meshFaces[oneFaceInd];
                    Info<<"Face "<<k<<" area: "<<oneFace.mag(meshPoints)<<"::";
                    for(int kk=0;kk<oneFace.size();kk++)
                    {
                        Info<<oneFace[kk]<<meshPoints[oneFace[kk]]<<"->";
                    }
                    Info<<endl;
                }
                
                Info<<endl;
                Info<<"maxCellSize:"<<maxCellSize<<endl;
                Info<<"minCellSize:"<<minCellSize<<endl;
                Info<<"avgCellSize:"<<avgCellSize<<endl;
            
                FatalErrorInFunction
                <<"Normal vector is "<<normalFace<<" while faceCentreToNeighbourCentre is "<<faceCentreToOwnerCentre<<"!"
                <<" They must have the same direction"
                << exit(FatalError);
            }
        }
        
        if(i<neighbour.size())
        {
            label neighbourCell = neighbour[i];
            vector centreToNeighbourCentre = meshCells[neighbourCell].centre(meshPoints,meshFaces)-centreFace;
            if((centreToNeighbourCentre && normalFace)<=0)
            {
                cell neighCell = meshCells[neighbourCell];
                scalar thisCellSize = neighCell.mag(meshPoints,meshFaces);
                if(thisCellSize >= (maxCellSize*partialThreeshold*(1.f/1e10)) && norm2(normalFace)!=0 &&
                   area>=partialThreeshold*maxFaceSize*(1.f/1e10)
                )
                // too small cells might fail in this condition because of rounding error.
                // as a result they are exempt here.
                // the same is for zero faces with a zero normal vector
                {
                    Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
                    if(i < neighbour.size())
                        Info<<"Neighbor:"<<neighbour[i]<<" ";
                    Info<<endl;
                    Info<<"centreToNeighbourCentre:"<<centreToNeighbourCentre;
                    Info<<"neighbourCentre: "<<meshCells[neighbourCell].centre(meshPoints,meshFaces);
                    Info<<endl;
                    for(int k=0;k<meshFaces[i].size();k++)
                    {
                        Info<<meshFaces[i][k]<<meshPoints[meshFaces[i][k]]<<"->";
                    }
                    Info<<endl;
                    for(int k=0;k<meshFaces[i].size();k++)
                    {
                        Info<<meshPoints[meshFaces[i][(k+1)%meshFaces[i].size()]]-meshPoints[meshFaces[i][k%meshFaces[i].size()]]<<"->";
                    }
                    Info<<endl;
                    Info<<" with centre:"<<centreFace;
                    Info<<" and normal vector:"<<normalFace;
                    Info<<" and area:"<<area<<endl<<endl;
                
                    Info<<"Neighbour Cell"<<endl;
                    cell oneCell = meshCells[neighbourCell];
                    Info<<"Cell centre is: "<<oneCell.centre(meshPoints,meshFaces)<<endl;
                    Info<<"Cell: "<<oneCell<<endl;
                    Info<<"Cell Size: "<<oneCell.size()<<endl;
                    Info<<"Cell volume: "<<oneCell.mag(meshPoints,meshFaces)<<endl;
            
                    for(int k=0;k<oneCell.size();k++)
                    {
                        label oneFaceInd = oneCell[k];
                        face oneFace = meshFaces[oneFaceInd];
                        Info<<"Face "<<k<<" area: "<<oneFace.mag(meshPoints)<<"::";
                        for(int kk=0;kk<oneFace.size();kk++)
                        {
                            Info<<oneFace[kk]<<meshPoints[oneFace[kk]]<<"->";
                        }
                        Info<<endl;
                    }
                
                    Info<<endl;
                    Info<<"maxCellSize:"<<maxCellSize<<endl;
                    Info<<"minCellSize:"<<minCellSize<<endl;
                    Info<<"avgCellSize:"<<avgCellSize<<endl;
            
                    FatalErrorInFunction
                    <<"Normal vector is "<<normalFace<<" while faceCentreToNeighbourCentre is "<<centreToNeighbourCentre<<"!"
                    <<" They must have the same direction"
                    << exit(FatalError);
                }
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
            << exit(FatalError);
        }
        if(mag == 0 && cellsToSide_[i] != -1 && norm2(cellCentre)!=0)
        {
            Info<<"Cell:"<<i<<" with "<<meshCells[i].nFaces()<<" faces |"<<endl;
            Info<<"meshPoints.size(): "<<meshPoints.size()<<endl;
            for(int k=0;k<meshCells[i].size();k++)
            {
                Info<<meshCells[i][k]<<"->";
            }
            Info<<endl;
            Info<<" with centre:"<<cellCentre;
            Info<<" and volume:"<<mag<<endl;
            
            Info<<"Cell"<<endl;
            cell oneCell = meshCells[i];
            Info<<"Cell centre is: "<<oneCell.centre(meshPoints,meshFaces)<<endl;
            Info<<"Cell: "<<oneCell<<endl;
            Info<<"Cell Size: "<<oneCell.size()<<endl;
            Info<<"Cell volume: "<<oneCell.mag(meshPoints,meshFaces)<<endl;
            
            for(int k=0;k<oneCell.size();k++)
            {
                label oneFaceInd = oneCell[k];
                face oneFace = meshFaces[oneFaceInd];
                Info<<"Face "<<k<<" area: "<<oneFace.mag(meshPoints)<<"::";
                for(int kk=0;kk<oneFace.size();kk++)
                {
                    Info<<oneFace[kk]<<meshPoints[oneFace[kk]]<<"->";
                }
                Info<<endl;
            }
            
            FatalErrorInFunction
            << "Cell cannot have Volume equal zero while being on side:"<<cellsToSide_[i]
            << exit(FatalError);
        }
        
        
        //Test if centre is really inside cell
        for(int a=0;a<meshCells[i].size();a++)
        {
            label oneFaceInd = meshCells[i][a];
            scalar oneFaceArea = meshFaces[oneFaceInd].mag(meshPoints);

            vector thisFaceNormal = meshFaces[oneFaceInd].normal(meshPoints);
            if(owner[oneFaceInd] != i)
                thisFaceNormal = -1*thisFaceNormal;
            vector thisFaceCentre = meshFaces[oneFaceInd].centre(meshPoints);
            if((((thisFaceCentre-cellCentre) && thisFaceNormal)<=0))
            {
                if(mag >= (maxCellSize*partialThreeshold*(1.f/1e10)) && norm2(thisFaceNormal) && 
                   oneFaceArea >= (maxFaceSize*partialThreeshold*(1.f/1e10))
                )
                {
                    Info<<"Cell:"<<i<<" with "<<meshCells[i].nFaces()<<"faces |";
                    for(int k=0;k<meshCells[i].size();k++)
                    {
                        Info<<meshCells[i][k]<<"->";
                    }
                    Info<<" with centre:"<<cellCentre;
                    Info<<" and volume:"<<mag<<endl;
                    
                    cell oneCell = meshCells[i];
                    Info<<"Cell centre is: "<<oneCell.centre(meshPoints,meshFaces)<<endl;
                    Info<<"Cell: "<<oneCell<<endl;
                    Info<<"Cell Size: "<<oneCell.size()<<endl;
                    Info<<"Cell volume: "<<oneCell.mag(meshPoints,meshFaces)<<endl;
            
                    for(int k=0;k<oneCell.size();k++)
                    {
                        label oneFaceInd = oneCell[k];
                        face oneFace = meshFaces[oneFaceInd];
                        Info<<"Face "<<k<<" area: "<<oneFace.mag(meshPoints)<<"::";
                        for(int kk=0;kk<oneFace.size();kk++)
                        {
                            Info<<oneFace[kk]<<meshPoints[oneFace[kk]]<<"->";
                        }
                        Info<<endl;
                    }
                
                    FatalErrorInFunction
                    << "Cell Face "<<a<<" has a normal "<<thisFaceNormal<<" but cellCentreToFaceCentre is "<<thisFaceCentre-cellCentre
                    << exit(FatalError);
                }                
            }
        }
    }
    Info<<": MESH IS CORRECT"<<endl;

}

void Foam::cutCellFvMesh::agglomerateSmallCells_cutNeg
(
    scalarList& newCellVolume,
    scalarList& oldCellVolume,
    scalar partialThreeshold
)
{
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span;
    
    Info<<endl;
    Info<<"Preprocessing of small cells ";
    t1 = std::chrono::high_resolution_clock::now();
    /*
    for(int i=0;i<oldCellVolume.size();i++)
    {
        Info<<"Cell "<<i<<" volume: "<<oldCellVolume[i]<<" splitted to minus side cell: "<<oldSplittedCellToNewMinusCell[i]<<endl;
        Info<<"Cell "<<i<<" splitted to plus side cell: "<<oldSplittedCellToNewPlusCell[i]<<endl;
    }
    */
    scalarList partialVolumeScale = scalarList(newCellVolume.size());
    //Info<<"new Cell Size: "<<newCellVolume.size()<<endl;

    label deletedCellsCount = 0;
    for(int i=0;i<deletedCellsList.size();i++)
    {
        if(deletedCellsList[i] == 1)
            deletedCellsCount++;
    }
    Info<<"1"<<endl;
    if(newCellVolume.size()+deletedCellsCount != oldCellVolume.size())
    {
        FatalErrorInFunction
        << "Must not happen!"
        << exit(FatalError); 
    }
    Info<<"2"<<endl;
    
    
    Info<<endl<<"deletedCellsList: "<<deletedCellsList<<endl;
    Info<<endl<<"oldSplittedCellToNewPlusCell: "<<oldSplittedCellToNewPlusCell<<endl;


    Info<<"3"<<endl;
    
    if(newCellVolume.size() != mapNewCellsToOldCells.size())
        FatalErrorInFunction<< "Must not happen!"<< exit(FatalError);
    
    for(int i=0;i<newCellVolume.size();i++)
    {
        if(mapNewCellsToOldCells[i] == -1)
            partialVolumeScale[i] = 1;
        else
            partialVolumeScale[i] = newCellVolume[i]/oldCellVolume[mapNewCellsToOldCells[i]];
    }
    Info<<"4"<<endl;
    /*
    for(int i=0;i<newCellVolume.size();i++)
    {
        Info<<"new Cell "<<i<<" has volume scale: "<<partialVolumeScale[i]<<endl;
    }
    */
    
    const cellList& newCells = this->cells();
    const faceList& faces = this->faces();
    const labelList& owner   = this->faceOwner();
    const labelList& neighbour = this->faceNeighbour();
    const pointField& points = this->points();
    
    DynamicList<DynamicList<DynamicList<label>>> possibleMergeFaces;
    possibleMergeFaces.setSize(newCellVolume.size());
    DynamicList<DynamicList<DynamicList<label>>> possibleMergeCells;
    possibleMergeCells.setSize(newCellVolume.size());
    DynamicList<DynamicList<scalar>> possibleMergeFaceArea;
    possibleMergeFaceArea.setSize(newCellVolume.size());
    DynamicList<DynamicList<bool>> possibleMergeFaceSufficient;
    possibleMergeFaceSufficient.setSize(newCellVolume.size());
    DynamicList<bool> oneMergeFaceSufficient;
    oneMergeFaceSufficient.setSize(newCellVolume.size());
    DynamicList<bool> mergeNecessary;
    mergeNecessary.setSize(newCellVolume.size());
    DynamicList<label> MergeCell;
    label neighbourCell = -1;
    scalar neighbourCellPartialVolume;
    Info<<"Create merge Cells"<<endl;
    for(int i=0;i<newCellVolume.size();i++)
    {
        Info<<i<<endl;
        mergeNecessary[i] = false;
        if((partialVolumeScale[i] < 1) && (partialVolumeScale[i] < partialThreeshold))
        {
            mergeNecessary[i] = true;
            /*
            Info<<"new Cell "<<i<<" is signed for merge via faces ";
            for(int k=0;k<newCells[i].size();k++)
            {
                Info<<newCells[i][k]<<",";
            }
            Info<<" with cells:";
            */
            
            /*
             * Find merging cells for one to one merging
             */
            Info<<"1-1"<<endl;
            for(int k=0;k<newCells[i].size();k++)
            {
                if(newCells[i][k] < neighbour.size())
                {
                    if(owner[newCells[i][k]] == i)
                        neighbourCell = neighbour[newCells[i][k]];
                    else if(neighbour[newCells[i][k]] == i)
                        neighbourCell = owner[newCells[i][k]];
                    else
                        FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<< exit(FatalError);
                    
                    neighbourCellPartialVolume = partialVolumeScale[neighbourCell];
                    
                    if(neighbourCellPartialVolume + partialVolumeScale[i] >= partialThreeshold)
                    {
                        DynamicList<label> temp;
                        temp.append(newCells[i][k]);
                        possibleMergeFaces[i].append(temp);
                        
                        temp.setSize(0);
                        temp.append(neighbourCell);
                        possibleMergeCells[i].append(temp);

                        possibleMergeFaceArea[i].append(faces[possibleMergeFaces[i][possibleMergeCells[i].size()-1][0]].mag(points));
                        
                        possibleMergeFaceSufficient[i].append(true);
                    }
                }
            }
            
            /*
             * Find merging cells for one to four merging
             */
            Info<<"1-4"<<endl;
            for(int a=0;a<newCells[i].size();a++)
            {
                if(newCells[i][a] < neighbour.size())
                {
                    for(int b=a+1;b<newCells[i].size();b++)
                    {
                        if(newCells[i][b] < neighbour.size())
                        {
                            Info<<"Try at :"<<a<<" "<<b<<endl;
                            label face_a,face_b;
                            DynamicList<label> mergeFace;
                            label neighbour_a,neighbour_b,fourthCell_a,fourthCell_b,fourthCell_F;
                            
                            neighbour_a=-1;
                            neighbour_b=-1;
                            fourthCell_a=-1;
                            fourthCell_b=-1;
                            fourthCell_F=-1;
                            face_a = newCells[i][a];
                            face_b = newCells[i][b];
                            
                            if(owner[face_a] == i)
                                neighbour_a = neighbour[face_a];
                            else if(neighbour[face_a] == i)
                                neighbour_a = owner[face_a];
                            else
                                FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<< exit(FatalError);  

                            if(owner[face_b] == i)
                                neighbour_b = neighbour[face_b];
                            else if(neighbour[face_b] == i)
                                neighbour_b = owner[face_b];
                            else
                                FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<< exit(FatalError);  

                            Info<<"Cell "<<i<<" has neighbours:"<<neighbour_a<<" "<<neighbour_b<<endl;
                            for(int x=0;x<newCells[neighbour_a].size();x++)
                            {
                                Info<<"x:"<<x<<endl;
                                if(newCells[neighbour_a][x]==face_a ||
                                   newCells[neighbour_a][x] >= neighbour.size())
                                    continue;
                                    
                                if(owner[newCells[neighbour_a][x]]==neighbour_a)
                                    fourthCell_a = neighbour[newCells[neighbour_a][x]];
                                else if(neighbour[newCells[neighbour_a][x]]==neighbour_a)
                                    fourthCell_a = owner[newCells[neighbour_a][x]];
                                else
                                    FatalErrorInFunction<< "Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<< exit(FatalError);
                                
                                for(int y=0;y<newCells[neighbour_b].size();y++)
                                {
                                    Info<<"y:"<<y<<endl;
                                    if(newCells[neighbour_b][y]==face_b ||
                                       newCells[neighbour_b][y] >= neighbour.size())
                                        continue;
                                        
                                    if(owner[newCells[neighbour_b][y]]==neighbour_b)
                                        fourthCell_b = neighbour[newCells[neighbour_b][y]];
                                    else if(neighbour[newCells[neighbour_b][y]]==neighbour_b)
                                        fourthCell_b = owner[newCells[neighbour_b][y]];
                                    else
                                        FatalErrorInFunction<< "Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<< exit(FatalError);  

                                    Info<<"Fourth cell candidates: "<<fourthCell_a<<" "<<fourthCell_b<<endl;
                                    if(fourthCell_a==fourthCell_b)
                                    {
                                        mergeFace.append(newCells[neighbour_a][x]);
                                        mergeFace.append(newCells[neighbour_b][y]);
                                        fourthCell_F = fourthCell_a;
                                    }
                                }
                            }
                            Info<<mergeFace<<endl;
                            if(mergeFace.size()==0)
                                continue;
                            if(mergeFace.size()!=2)
                            {
                                FatalErrorInFunction
                                << "Agglomeration cell not found for fourthCell_F!"
                                << exit(FatalError);
                            }
                            if(fourthCell_F==-1)
                            {
                                FatalErrorInFunction
                                << "Agglomeration cell not found for all cells!"
                                << exit(FatalError);
                            }
                            DynamicList<label> mergeCells;
                            mergeCells.append(neighbour_a);
                            mergeCells.append(neighbour_b);
                            mergeCells.append(fourthCell_F);
                            DynamicList<label> temp;
                            temp.append(face_a);
                            temp.append(face_b);
                            temp.append(mergeFace);
                            
                            neighbourCellPartialVolume = 0;
                            for(int p=0;p<mergeCells.size();p++)
                            {
                                neighbourCellPartialVolume += partialVolumeScale[mergeCells[p]];
                            }
                    
                            if(neighbourCellPartialVolume + partialVolumeScale[i] >= partialThreeshold)
                            {                         
                                possibleMergeFaces[i].append(temp);
                                possibleMergeCells[i].append(mergeCells);
                                possibleMergeFaceArea[i].append(0.0); // Set to zero to make the merge last priority
                                possibleMergeFaceSufficient[i].append(true);
                            }
                        }
                    }
                }
            }

            /*
             * Find merging cells for one to eight merging
             */
            Info<<"1-8"<<endl;
            for(int a=0;a<newCells[i].size();a++)
            {
                if(newCells[i][a] < neighbour.size())
                {
                    for(int b=a+1;b<newCells[i].size();b++)
                    {
                        if(newCells[i][b] < neighbour.size())
                        {
                            for(int c=b+1;c<newCells[i].size();c++)
                            {
                                if(newCells[i][c] < neighbour.size())
                                {
                                    Info<<"Try at :"<<a<<" "<<b<<" "<<c<<endl;
                                    label face_a,face_b,face_c;
                                    DynamicList<label> mergeFace;
                                    label neighbour_a,neighbour_b,fourthCell_a,fourthCell_b,
                                          fourthCell_F,neighbour_c,sixthCell_a,sixthCell_b,sixthCell_F,
                                          seventhCell_a,seventhCell_b,seventhCell_F,
                                          eightCell_a,eightCell_b,eightCell_c,eightCell_F;
                                    fourthCell_a=-1;
                                    fourthCell_b=-1;
                                    sixthCell_a=-1;
                                    sixthCell_b=-1;
                                    seventhCell_a=-1;
                                    seventhCell_b=-1;
                                    seventhCell_F=-1;
                                    sixthCell_F=-1;
                                    neighbour_c=-1;
                                    fourthCell_F=-1;
                                    neighbour_a=-1;
                                    neighbour_b=-1;
                                    eightCell_a=-1;
                                    eightCell_b=-1;
                                    eightCell_c=-1;
                                    eightCell_F=-1;
                                    
                                    face_a = newCells[i][a];
                                    face_b = newCells[i][b];
                                    face_c = newCells[i][c];                               
                            
                                    if(owner[face_a] == i)
                                        neighbour_a = neighbour[face_a];
                                    else if(neighbour[face_a] == i)
                                        neighbour_a = owner[face_a];
                                    else
                                        FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<<exit(FatalError);  

                                    if(owner[face_b] == i)
                                        neighbour_b = neighbour[face_b];
                                    else if(neighbour[face_b] == i) 
                                        neighbour_b = owner[face_b];
                                    else
                                        FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<<exit(FatalError);
                                    
                                    if(neighbour_a==-1||neighbour_b==-1)
                                    {
                                        FatalErrorInFunction
                                        << "Agglomeration cell not found for neighbour_a and neighbour_b!"
                                        << exit(FatalError);
                                    }
                                    
                                    for(int x=0;x<newCells[neighbour_a].size();x++)
                                    {
                                        if(newCells[neighbour_a][x]==face_a ||
                                           newCells[neighbour_a][x] >= neighbour.size())
                                            continue;
                                        if(owner[newCells[neighbour_a][x]]==neighbour_a)
                                        {
                                            fourthCell_a = neighbour[newCells[neighbour_a][x]];
                                        }
                                        else if(neighbour[newCells[neighbour_a][x]]==neighbour_a)
                                        {   
                                            fourthCell_a = owner[newCells[neighbour_a][x]];
                                        }
                                        else
                                        {
                                            FatalErrorInFunction
                                            << "Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"
                                            << exit(FatalError);  
                                        }                                
                                        for(int y=0;y<newCells[neighbour_b].size();y++)
                                        {                                
                                            if(newCells[neighbour_b][y]==face_b ||
                                               newCells[neighbour_b][y] >= neighbour.size())
                                                continue;
                                            if(owner[newCells[neighbour_b][y]]==neighbour_b)
                                            {
                                                fourthCell_b = neighbour[newCells[neighbour_b][y]];
                                            }
                                            else if(neighbour[newCells[neighbour_b][y]]==neighbour_b)
                                            {   
                                                fourthCell_b = owner[newCells[neighbour_b][y]];
                                            }
                                            else
                                            {
                                                FatalErrorInFunction
                                                << "Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"
                                                << exit(FatalError);  
                                            }
                                            if(fourthCell_a==fourthCell_b)
                                            {
                                                mergeFace.append(newCells[neighbour_a][x]);
                                                mergeFace.append(newCells[neighbour_b][y]);
                                                fourthCell_F = fourthCell_a;
                                            }
                                        }
                                    }
                                    if(mergeFace.size()==0)
                                        continue;
                                    if(mergeFace.size()!=2)
                                    {
                                        FatalErrorInFunction
                                        << "Agglomeration cell not found for all cells!"
                                        << exit(FatalError);
                                    }
                                    if(fourthCell_F==-1)
                                    {
                                        FatalErrorInFunction
                                        << "Agglomeration cell not found for fourthCell_F!"
                                        << exit(FatalError);
                                    }   
                                    
                                    DynamicList<label> mergeCells;
                                    mergeCells.append(neighbour_a);
                                    mergeCells.append(neighbour_b);
                                    mergeCells.append(fourthCell_F);
                                    DynamicList<label> mergeFaces;
                                    mergeFaces.append(face_a);
                                    mergeFaces.append(face_b);
                                    mergeFaces.append(mergeFace);
                                    
                                    Info<<"newCells["<<i<<"]: "<<newCells[i]<<endl;
                                    Info<<"i:"<<i<<endl;
                                    Info<<"owner["<<face_c<<"]:"<<owner[face_c]<<endl;
                                    Info<<"neighbour["<<face_c<<"]:"<<neighbour[face_c]<<endl;
                                    if(owner[face_c] == i)
                                        neighbour_c = neighbour[face_c];
                                    else if(neighbour[face_c] == i)
                                        neighbour_c = owner[face_c];
                                    else
                                        FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<<exit(FatalError);
                                    
                                    if(neighbour_c==-1)
                                    {
                                        FatalErrorInFunction
                                        << "Agglomeration cell not found for neighbour_c!"
                                        << exit(FatalError);
                                    }   
                                    
                                    mergeCells.append(neighbour_c);
                                    mergeFaces.append(face_c);
                                    
                                    /*
                                     * Find connecting cell neighbour_a -> neighbour_c
                                     */
                                    mergeFace.setSize(0);
                                    for(int x=0;x<newCells[neighbour_a].size();x++)
                                    {
                                        if(newCells[neighbour_a][x]==face_a ||
                                           newCells[neighbour_a][x] >= neighbour.size())
                                            continue;
                                            
                                        if(owner[newCells[neighbour_a][x]]==neighbour_a)
                                            sixthCell_a = neighbour[newCells[neighbour_a][x]];
                                        else if(neighbour[newCells[neighbour_a][x]]==neighbour_a) 
                                            sixthCell_a = owner[newCells[neighbour_a][x]];
                                        else
                                            FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<<exit(FatalError);  
                            
                                        for(int y=0;y<newCells[neighbour_c].size();y++)
                                        {                                
                                            if(newCells[neighbour_c][y]==face_b ||
                                               newCells[neighbour_c][y] >= neighbour.size())
                                                continue;
                                                
                                            if(owner[newCells[neighbour_c][y]]==neighbour_c)
                                                sixthCell_b = neighbour[newCells[neighbour_c][y]];
                                            else if(neighbour[newCells[neighbour_c][y]]==neighbour_c) 
                                                sixthCell_b = owner[newCells[neighbour_c][y]];
                                            else
                                                FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<<exit(FatalError);  
                                            
                                            if(sixthCell_a==sixthCell_b)
                                            {
                                                mergeFace.append(newCells[neighbour_a][x]);
                                                mergeFace.append(newCells[neighbour_c][y]);
                                                sixthCell_F = sixthCell_a;
                                            }
                                        }
                                    }
                                    if(mergeFace.size()==0)
                                        continue;
                                    if(mergeFace.size()!=2)
                                    {
                                        FatalErrorInFunction
                                        << "Agglomeration cell not found for all cells!"
                                        << exit(FatalError);
                                    }
                                    if(sixthCell_F==-1)
                                    {
                                        FatalErrorInFunction
                                        << "Agglomeration cell not found for sixthCell_F!"
                                        << exit(FatalError);
                                    }
                                    mergeCells.append(sixthCell_F);
                                    mergeFaces.append(mergeFace);

                                    /*
                                     * Find connecting cell neighbour_b -> neighbour_c
                                     */
                                    mergeFace.setSize(0);
                                    for(int x=0;x<newCells[neighbour_b].size();x++)
                                    {
                                        if(newCells[neighbour_b][x]==face_b ||
                                           newCells[neighbour_b][x] >= neighbour.size())
                                            continue;
                                            
                                        if(owner[newCells[neighbour_b][x]]==neighbour_b)
                                            seventhCell_a = neighbour[newCells[neighbour_b][x]];
                                        else if(neighbour[newCells[neighbour_b][x]]==neighbour_b) 
                                            seventhCell_a = owner[newCells[neighbour_b][x]];
                                        else
                                            FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<<exit(FatalError);  
                            
                                        for(int y=0;y<newCells[neighbour_c].size();y++)
                                        {                                
                                            if(newCells[neighbour_c][y]==face_b ||
                                               newCells[neighbour_c][y] >= neighbour.size())
                                                continue;
                                                
                                            if(owner[newCells[neighbour_c][y]]==neighbour_c)
                                                seventhCell_b = neighbour[newCells[neighbour_c][y]];
                                            else if(neighbour[newCells[neighbour_c][y]]==neighbour_c) 
                                                seventhCell_b = owner[newCells[neighbour_c][y]];
                                            else
                                                FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<<exit(FatalError);  
                                            
                                            if(seventhCell_a==seventhCell_b)
                                            {
                                                mergeFace.append(newCells[neighbour_b][x]);
                                                mergeFace.append(newCells[neighbour_c][y]);
                                                seventhCell_F = seventhCell_a;
                                            }
                                        }
                                    }
                                    if(mergeFace.size()==0)
                                        continue;
                                    if(mergeFace.size()!=2)
                                    {
                                        FatalErrorInFunction
                                        << "Agglomeration cell not found for all cells!"
                                        << exit(FatalError);
                                    }
                                    if(seventhCell_F==-1)
                                    {
                                        FatalErrorInFunction
                                        << "Agglomeration cell not found for seventhCell_F!"
                                        << exit(FatalError);
                                    }                                    
                                    mergeCells.append(seventhCell_F);
                                    mergeFaces.append(mergeFace);
                                    
                                    /*
                                     * Find connecting cell fourthCell -> sixthCell -> seventhCell
                                     */
                                    mergeFace.setSize(0);
                                    for(int x=0;x<newCells[fourthCell_F].size();x++)
                                    {
                                        if(newCells[fourthCell_F][x] >= neighbour.size())
                                            continue;
                                        
                                        Info<<"newCells["<<fourthCell_F<<"]: "<<newCells[fourthCell_F]<<endl;
                                        Info<<"i:"<<fourthCell_F<<endl;
                                        Info<<"owner["<<newCells[fourthCell_F][x]<<"]:"<<owner[newCells[fourthCell_F][x]]<<endl;
                                        Info<<"neighbour["<<newCells[fourthCell_F][x]<<"]:"<<neighbour[newCells[fourthCell_F][x]]<<endl;
                                        
                                        if(owner[newCells[fourthCell_F][x]]==fourthCell_F)
                                            eightCell_a = neighbour[newCells[fourthCell_F][x]];
                                        else if(neighbour[newCells[fourthCell_F][x]]==fourthCell_F) 
                                            eightCell_a = owner[newCells[fourthCell_F][x]];
                                        else
                                            FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<<exit(FatalError);  
                            
                                        for(int y=0;y<newCells[sixthCell_F].size();y++)
                                        {                                
                                            if(newCells[sixthCell_F][y] >= neighbour.size())
                                                continue;
                                                
                                            if(owner[newCells[sixthCell_F][y]]==sixthCell_F)
                                                eightCell_b = neighbour[newCells[sixthCell_F][y]];
                                            else if(neighbour[newCells[sixthCell_F][y]]==sixthCell_F) 
                                                eightCell_b = owner[newCells[sixthCell_F][y]];
                                            else
                                                FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<<exit(FatalError);
                                            
                                            for(int z=0;z<newCells[seventhCell_F].size();z++)
                                            {                                
                                                if(newCells[seventhCell_F][y] >= neighbour.size())
                                                    continue;
                                                
                                                if(owner[newCells[seventhCell_F][z]]==seventhCell_F)
                                                    eightCell_c = neighbour[newCells[seventhCell_F][z]];
                                                else if(neighbour[newCells[seventhCell_F][z]]==seventhCell_F) 
                                                    eightCell_c = owner[newCells[seventhCell_F][z]];
                                                else
                                                    FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<<exit(FatalError);  
                                            
                                                
                                            
                                                if(eightCell_a==eightCell_b && eightCell_b==eightCell_c)
                                                {
                                                    mergeFace.append(newCells[fourthCell_F][x]);
                                                    mergeFace.append(newCells[sixthCell_F][y]);
                                                    mergeFace.append(newCells[seventhCell_F][z]);
                                                    eightCell_F = eightCell_a;
                                                }
                                            }
                                        }
                                    }
                                    if(mergeFace.size()==0)
                                        continue;
                            //Problem here
                                    if(mergeFace.size()!=3)
                                    {
                                        Info<<"mergeFace.size():"<<mergeFace.size()<<endl;
                                        FatalErrorInFunction
                                        << "Agglomeration cell not found for all cells!"
                                        << exit(FatalError);
                                    }
                                    mergeCells.append(eightCell_F);
                                    mergeFaces.append(mergeFace);
                                    
                                    neighbourCellPartialVolume = 0;
                                    for(int p=0;p<mergeCells.size();p++)
                                    {
                                        neighbourCellPartialVolume += partialVolumeScale[mergeCells[p]];
                                    }
                    
                                    if(neighbourCellPartialVolume + partialVolumeScale[i] >= partialThreeshold)
                                    { 
                                        Info<<"Merge Cells at :"<<a<<" "<<b<<" "<<c<<endl;
                                        possibleMergeFaces[i].append(mergeFaces);
                                        possibleMergeCells[i].append(mergeCells);
                                        possibleMergeFaceArea[i].append(0.0); // Set to zero to make the merge last priority
                                        possibleMergeFaceSufficient[i].append(true);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            //Info<<endl;
            /*
            oneMergeFaceSufficient[i] = false;
            for(int k=0;k<possibleMergeFaceSufficient.size();k++)
            {
                if(possibleMergeFaceSufficient[i][k])
                {
                    oneMergeFaceSufficient[i] = true;
                    break;
                }
            }
            */
            if(possibleMergeCells[i].size() == 0)
            {
                Info<<endl<<"Problem in cell "<<i<<" with partial volume "<< partialVolumeScale[i]<<endl;
                Info<<"\tNeighbours are:"<<endl;
                for(int k=0;k<newCells[i].size();k++)
                {
                    if(newCells[i][k] < neighbour.size())
                    {
                        if(owner[newCells[i][k]] == i)
                        {
                            neighbourCell = neighbour[newCells[i][k]];
                        }
                        else if(neighbour[newCells[i][k]] == i)
                        {   
                            neighbourCell = owner[newCells[i][k]];
                        }
                        neighbourCellPartialVolume = partialVolumeScale[neighbourCell];
                        
                        Info<<"\tCell:"<<neighbourCell<<" partialVol:"<<neighbourCellPartialVolume<<
                            " combinedPartialVol:"<<(neighbourCellPartialVolume + partialVolumeScale[i])<<" with threshold:"<<partialThreeshold
                            <<endl;
                    }
                }
                FatalErrorInFunction
                << "Agglomeration cell not found for all cells!"
                << exit(FatalError);
            }
        }
        /*
        Info<<"cell:"<<i<<" mergeNecessary:"<<mergeNecessary[i]<<" partialVolumeScale:"<<partialVolumeScale[i]<<endl;
        */
    }
    
    Info<<"Created merge Cells"<<endl;
    for(int i=0;i<newCellVolume.size();i++)
    {
        for(int j=0;j<possibleMergeFaces[i].size();j++)
        {
            label numCellMerge=-1;
            if(possibleMergeFaces[i][j].size() == 1)
                numCellMerge = 2;
            else if(possibleMergeFaces[i][j].size() == 4)
                numCellMerge = 4;
            else if(possibleMergeFaces[i][j].size() == 12)
                numCellMerge = 8;
            else
                FatalErrorInFunction<<"Wrong number of merge faces"<< exit(FatalError);
            
            if(possibleMergeCells[i][j].size()+1 != numCellMerge)
            {
                Info<<"numCellMerge: "<<numCellMerge<<endl;
                Info<<"possibleMergeCells[i][j].size(): "<<possibleMergeCells[i][j].size()<<endl;
                Info<<"possibleMergeFaces[i][j].size(): "<<possibleMergeFaces[i][j].size()<<endl;
                Info<<"possibleMergeCells[i][j]: "<<possibleMergeCells[i][j]<<endl;
                Info<<"possibleMergeFaces[i][j]: "<<possibleMergeFaces[i][j]<<endl;
                FatalErrorInFunction<<"Number of merge cells does not match"<< exit(FatalError);
            }


            std::unordered_multiset<label> cellSet;
            for(int k=0;k<possibleMergeFaces[i][j].size();k++)
            {
                cellSet.insert(owner[possibleMergeFaces[i][j][k]]);
                cellSet.insert(neighbour[possibleMergeFaces[i][j][k]]);
            }
            
            label cellDuplicationNum=-1;
            if(numCellMerge==2) cellDuplicationNum = 1;
            else if(numCellMerge==4) cellDuplicationNum = 2;
            else if(numCellMerge==8) cellDuplicationNum = 3;
            else FatalErrorInFunction<<"Somethings wrong here"<< exit(FatalError);
            
            for(int k=0;k<possibleMergeCells[i][j].size();k++)
            {
                if(cellSet.count(possibleMergeCells[i][j][k])!=static_cast<long unsigned int>(cellDuplicationNum))
                    FatalErrorInFunction<<"Merge Cells and Faces do not match!"<< exit(FatalError);
            }
        }
    }
    
    // Sort possible merging cell by respect to face area biggest to smallest
    for(int i=0;i<possibleMergeFaceArea.size();i++)
    {
        int j;
        scalar keyArea;
        DynamicList<label> keyFaces,keyCells;
        bool keySuff;
        for(int k=1;k<possibleMergeFaceArea[i].size();k++)
        {
            keyArea = possibleMergeFaceArea[i][k];
            keyFaces = possibleMergeFaces[i][k];
            keyCells = possibleMergeCells[i][k];
            keySuff = possibleMergeFaceSufficient[i][k];
            j = k-1;
            while(j>=0 && possibleMergeFaceArea[i][j] < keyArea)
            {
                possibleMergeFaceArea[i][j+1] = possibleMergeFaceArea[i][j];
                possibleMergeFaces[i][j+1] = possibleMergeFaces[i][j];
                possibleMergeCells[i][j+1] = possibleMergeCells[i][j];
                possibleMergeFaceSufficient[i][j+1] = possibleMergeFaceSufficient[i][j];
                j--;
            }
            possibleMergeFaceArea[i][j+1] = keyArea;
            possibleMergeFaces[i][j+1] = keyFaces;
            possibleMergeCells[i][j+1] = keyCells;
            possibleMergeFaceSufficient[i][j+1] = keySuff;
        }
    }
    Info<<"Sorted merge Cells"<<endl;
    
    for(int i=0;i<possibleMergeFaces.size();i++)
    {
        for(int j=1;j<possibleMergeFaces[i].size();j++)
        {
            if(possibleMergeFaces[i][j-1].size() > possibleMergeFaces[i][j].size())
                FatalErrorInFunction<<"Invalid sorting"<<exit(FatalError);
        }
    }
    
//Test for correct merge candidates
    scalar factor = 1/partialThreeshold;
    scalar minCellVol = newCells[0].mag(points,faces);
    label minCellInd = 0;
    scalar maxCellVol = newCells[0].mag(points,faces);
    label maxCellInd = 0;
    scalar CellVolAvg = 0;
    
    scalar vol,neighborVol;
    for(int i=0;i<newCells.size();i++)
    {
        vol = newCells[i].mag(points,faces);
        CellVolAvg += vol;
        if(vol > maxCellVol)
        {
            maxCellVol = vol;
            maxCellInd = i;
        }
        if(vol < minCellVol)
        {
            minCellVol = vol;
            minCellInd = i;
        }
    }
    CellVolAvg /= newCells.size();

    
    Info<<endl<<"Minimum cell "<<minCellInd<<" vol:"<<minCellVol
    <<endl<<"Maximum cell "<<maxCellInd<<" vol:"<<maxCellVol<<endl;
    Info<<" Average vol was:"<<CellVolAvg<<endl;    
    
    for(int i=0;i<newCells.size();i++)
    {
        vol = newCells[i].mag(points,faces); 
        if((vol*factor) < maxCellVol)
        {
            if(mergeNecessary[i] == false)
            {
                FatalErrorInFunction
                << "Not to merge but necessary"
                << exit(FatalError);
            }
            if(possibleMergeCells[i].size() == 0)
            {
                FatalErrorInFunction
                << "No merge data available"
                << exit(FatalError);
            }            
            for(int k=0;k<possibleMergeCells[i].size();k++)
            {
                neighborVol = 0;
                for(int s=0;s<possibleMergeCells[i][k].size();s++)
                {
                    neighborVol += newCells[possibleMergeCells[i][k][s]].mag(points,faces);
                }
                if((neighborVol+vol)*factor < maxCellVol)
                {
                    FatalErrorInFunction
                    << "Not sufficient merge data"
                    << exit(FatalError);
                }
            }
        }
    }
    
    scalar partialVol;
    for(int i=0;i<newCells.size();i++)
    {
        if((partialVolumeScale[i] < 1) && (partialVolumeScale[i] < partialThreeshold))
        {
            if(possibleMergeCells[i].size()==0)
            {
                FatalErrorInFunction
                << "Merge Face "<<i<<" with no merge partners!"
                << exit(FatalError);  
            }
            if( (possibleMergeFaces[i].size()!=possibleMergeCells[i].size())&&
                (possibleMergeCells[i].size()!=possibleMergeFaceArea[i].size()))
            {
                FatalErrorInFunction
                << "Data error"
                << exit(FatalError);  
            }
            for(int k=0;k<possibleMergeCells[i].size();k++)
            {
                partialVol = partialVolumeScale[i];
                for(int s=0;s<possibleMergeCells[i][k].size();s++)
                {
                    partialVol += partialVolumeScale[possibleMergeCells[i][k][s]];
                }
                if(partialVol < partialThreeshold)
                {
                    FatalErrorInFunction
                    << "Data 2 error"
                    << exit(FatalError);  
                }
            }
        }
    }
    
    for(int i=0;i<possibleMergeFaces.size();i++)
    {
        for(int j=0;j<possibleMergeFaces[i].size();j++)
        {
            // Test one option
            std::unordered_multiset<label> optionMergeCellsSet;
            DynamicList<label> optionMergeCellsList;
            for(int k=0;k<possibleMergeFaces[i][j].size();k++)
            {
                if(possibleMergeFaces[i][j][k] >= neighbour.size())
                {
                    FatalErrorInFunction<<"Boundary face is merge face"<<exit(FatalError);
                }
                label ownerCell = owner[possibleMergeFaces[i][j][k]];
                label neighborCell = neighbour[possibleMergeFaces[i][j][k]];
                if(optionMergeCellsSet.find(ownerCell) == optionMergeCellsSet.end())
                    optionMergeCellsList.append(ownerCell);
                if(optionMergeCellsSet.find(neighborCell) == optionMergeCellsSet.end())
                    optionMergeCellsList.append(neighborCell);
                optionMergeCellsSet.insert(ownerCell);
                optionMergeCellsSet.insert(neighborCell);
            }
            if(optionMergeCellsList.size()<=0)
                FatalErrorInFunction<<"Merging option with no cells!"<<exit(FatalError);
            
            label mergeCellMult = optionMergeCellsSet.count(optionMergeCellsList[0]);
            for(int k=1;k<optionMergeCellsList.size();k++)
            {
                if(optionMergeCellsSet.count(optionMergeCellsList[k]) != static_cast<long unsigned int>(mergeCellMult))
                {
                    FatalErrorInFunction<<"Different multiplcity of cells in "<<
                    optionMergeCellsList.size()<<" cell merging "<<exit(FatalError);
                }
            }
            
            if(mergeCellMult==1 && optionMergeCellsList.size() == 2)
            {
                //2 Cell merging
                if(possibleMergeFaces[i][j].size() != 1)
                    FatalErrorInFunction<<"Error in 2 Cell merging! "<<exit(FatalError);
            }
            else if(mergeCellMult==2 && optionMergeCellsList.size() == 4)
            {
                //4 Cell merging
                if(possibleMergeFaces[i][j].size() != 4)
                    FatalErrorInFunction<<"Error in 4 Cell merging! "<<exit(FatalError);
            }
            else if(mergeCellMult==3 && optionMergeCellsList.size() == 8)
            {
                //4 Cell merging
                if(possibleMergeFaces[i][j].size() != 12)
                    FatalErrorInFunction<<"Error in 8 Cell merging! "<<exit(FatalError);
            }
            else
            {
                Info<<"mergeCellMult: "<<mergeCellMult<<endl;
                Info<<"optionMergeCellsList.size() == "<<optionMergeCellsList.size()<<endl;
                Info<<"possibleMergeFaces[i][j].size() == "<<possibleMergeFaces[i][j].size()<<endl;
                FatalErrorInFunction<<"Inconsistent merge Option!"<<exit(FatalError);
            }
            std::unordered_set<label> avoidCellDuplicates;
            for(int k=0;k<possibleMergeCells[i][j].size();k++)
            {
                if(avoidCellDuplicates.find(possibleMergeCells[i][j][k]) == avoidCellDuplicates.end())
                    avoidCellDuplicates.insert(possibleMergeCells[i][j][k]);
                else
                    FatalErrorInFunction<<"Duplicate cell in Merge option!"<<exit(FatalError);
            }
            for(int k=0;k<possibleMergeCells[i][j].size();k++)
            {
                if(optionMergeCellsSet.count(possibleMergeCells[i][j][k]) != static_cast<long unsigned int>(mergeCellMult))
                    FatalErrorInFunction<<"Non matching mergeCell in Merge option!"<<exit(FatalError);
            }
        }
    }
//End: Test for correct merge candidates

    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Processing of small cells ";
    t1 = std::chrono::high_resolution_clock::now(); 

    /*
    std::unordered_set<label> cellReserved;
    DynamicList<DynamicList<label>> blockedCells;
    blockedCells.setSize(possibleMergeCells.size());
    
    
    labelList mergeFaceOfCell = searchDown_rec(possibleMergeFaceArea,possibleMergeFaces,possibleMergeCells,oneMergeFaceSufficient,mergeNecessary,0,blockedCells,cellReserved);
    */

    
    List<DynamicList<label>> mergeFaceOfCell = searchDown_iter_preBlock(owner,neighbour,
                                                possibleMergeFaceArea, possibleMergeFaces,
                                                possibleMergeCells,oneMergeFaceSufficient,
                                                mergeNecessary
                                               );
    
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Test for duplicate merging selection ";
    t1 = std::chrono::high_resolution_clock::now();
    
    std::unordered_set<label> usedFace;
    for(int i=0;i<mergeFaceOfCell.size();i++)
    {
        if(mergeFaceOfCell[i].size() == 0)
        {
            FatalErrorInFunction
            << "Cell with unwanted results"
            << exit(FatalError);  
        }
        if(mergeFaceOfCell[i][0] == -3)
        {
            FatalErrorInFunction
            << "Too small cell was not treated by backtracking algorithm"
            << exit(FatalError);  
        }
        if(mergeFaceOfCell[i][0] == -4)
        {
            FatalErrorInFunction
            << "Cell with unwanted results"
            << exit(FatalError);  
        }
        if(mergeFaceOfCell[i][0] < 0)
            continue;
        
        for(int s=0;s<mergeFaceOfCell[i].size();s++)
        {
            if(usedFace.find(mergeFaceOfCell[i][s]) == usedFace.end())
                usedFace.insert(mergeFaceOfCell[i][s]);
            else
            {
                Info<<endl<<endl;
                label fc = mergeFaceOfCell[i][s];
                label wnr = owner[fc];
                label nghbr = neighbour[fc];
            
                Info<<"Face "<<fc<<" merging of cell:"<<wnr<<" with Vol:"<<newCellVolume[wnr]<<
                "cell:"<<nghbr<<" with Vol:"<<newCellVolume[nghbr];
                FatalErrorInFunction
                << "Merge Face used twice!"
                << exit(FatalError);  
            }
        }
    }

    label selectedFace;
    std::unordered_set<label> usedCells;
    for(int i=0;i<mergeFaceOfCell.size();i++)
    {
        for(int s=0;s<mergeFaceOfCell[i].size();s++)
        {
            selectedFace = mergeFaceOfCell[i][s];
            if(selectedFace == -3)
            {
                Info<<"Face is: "<<selectedFace<<endl;
            
                FatalErrorInFunction
                << "One cell not treated by algorithm"
                << exit(FatalError);  
            }
        
            // Test that each cell to merge is merged
            if(selectedFace == -1 && mergeNecessary[i])
            {
                Info<<"Selected Face is: "<<selectedFace<<" but mergeNecessary["<<i<<"]"<<mergeNecessary[i]<<endl;
            
                FatalErrorInFunction
                << "Agglomeration cell not found but necessary!"
                << exit(FatalError);  
            }
            
            if((selectedFace == -1 || selectedFace == -2) && mergeFaceOfCell[i].size() != 1)
            {
                FatalErrorInFunction
                << "Not merged or merged by other cell but size of mergeList is not one!"
                << exit(FatalError);  
            }
        }
        
        if(mergeFaceOfCell[i][0] == -1 || mergeFaceOfCell[i][0] == -2)
        {
            continue;
        }
        
        for(int s=0;s<mergeFaceOfCell[i].size();s++)
        {
            selectedFace = mergeFaceOfCell[i][s];
            //Test that merge faces are legit faces
            if(selectedFace > neighbour.size() || selectedFace < 0)
            {
                Info<<endl<<"Merging Face: "<<selectedFace<<endl;
                FatalErrorInFunction
                << "Merging face is not an existing face!"
                << exit(FatalError);
            }
        }
        
        bool selectedFaceExists = false;
        for(int s=0;s<possibleMergeFaces[i].size();s++)
        {
            if(possibleMergeFaces[i][s] == mergeFaceOfCell[i])
                selectedFaceExists = true;
        }
        if(!selectedFaceExists)
        {
            FatalErrorInFunction
            << "Merging face is not in the mergingFace list!"
            << exit(FatalError);
        }
        
        label numMergFaces = mergeFaceOfCell[i].size();
        if(numMergFaces!=1 && numMergFaces!=4 && numMergFaces!=12)
        {
            Info<<endl<<"numMergFaces["<<i<<"]:"<<numMergFaces<<endl;
            FatalErrorInFunction
            << "Number of merging faces does not match!"
            << exit(FatalError);
        }
        
        //Test that merging faces are the correct ones
        std::unordered_multiset<label> cellSet;
        DynamicList<label> allCells;
        for(int s=0;s<mergeFaceOfCell[i].size();s++)
        {
            selectedFace = mergeFaceOfCell[i][s];
            if(neighbour.size() <= selectedFace)
            {
                Info<<"mergeFaceOfCell["<<i<<"]["<<s<<"]="<<mergeFaceOfCell[i][s]<<endl;
                Info<<"neighbour.size():"<<neighbour.size()<<endl;
                FatalErrorInFunction
                << "Merge via boundary face not possible "
                << exit(FatalError);
            }
            if(cellSet.count(neighbour[selectedFace]) == 0)
                allCells.append(neighbour[selectedFace]);
            if(cellSet.count(owner[selectedFace]) == 0)
                allCells.append(owner[selectedFace]);
            cellSet.insert(neighbour[selectedFace]);
            cellSet.insert(owner[selectedFace]);
        }
        
        if(allCells.size()!=2 && allCells.size()!=4 && allCells.size()!=8)
        {
            Info<<endl<<allCells<<endl;
            FatalErrorInFunction
            << "Number of merging cells does not match!"
            << exit(FatalError);
        }
        
        label cellMult=-1;
        if(numMergFaces==1)
        {
            cellMult = 1; //two cell merge
        }
        else if(numMergFaces==4)
        {
            cellMult = 2; //four cell merge
        }
        else if(numMergFaces==12)
        {
            cellMult = 3; //eight cell merge
        }
        else
            FatalErrorInFunction<<"Number of merge faces is wrong: "<<numMergFaces<<"! "<< exit(FatalError);
        
        for(int s=0;s<allCells.size();s++)
        {
            if(cellSet.count(allCells[s])!=static_cast<long unsigned int>(cellMult))
            {
                FatalErrorInFunction
                << "Wrong number of cell count "
                << exit(FatalError);
            }
        }
        
        Info<<endl;
        bool cellsMatch = false;
        bool partMatch;
        for(int s=0;s<possibleMergeCells[i].size();s++)
        {
            partMatch = false;
            if(possibleMergeCells[i][s].size() == allCells.size()-1)
            {
                partMatch = true;
                bool match=false;
                for(int w=0;w<allCells.size();w++)
                {
                    Info<<allCells[w]<<"|"<<match<<"|"<<partMatch<<endl;
                    match = false;
                    if(allCells[w] == i)
                    {
                        match = true;
                        continue;
                    }
                    for(int x=0;x<possibleMergeCells[i][s].size();x++)
                    {
                        if(possibleMergeCells[i][s][x] == allCells[w])
                            match = true;                            
                    }
                    if(match == false)
                        partMatch = false;
                }
            }
            if(partMatch)
                cellsMatch = true;
        }
        if(!cellsMatch)
        {
            Info<<endl<<"i:"<<i<<endl;
            Info<<"possibleMergeCells["<<i<<"]:"<<possibleMergeCells[i]<<endl;
            Info<<"allCells["<<i<<"]:"<<allCells<<endl;
            FatalErrorInFunction
            << "Cells do not match "
            << exit(FatalError);
        }
        
        for(int s=0;s<allCells.size();s++)
        {
            if(usedCells.find(allCells[s]) == usedCells.end())
                usedCells.insert(allCells[s]);
            else
            {
                Info<<"All Cells:"<<allCells<<endl;
                Info<<"Cell: "<<allCells[s]<<" used twice"<<endl;
                FatalErrorInFunction
                << "Merge Cell used twice!"
                << exit(FatalError);
            }
        }
    }
    
//TestSection
    {
        scalar factor = 1/partialThreeshold;
        const cellList& cell = this->cells();
        const faceList& face = this->faces();
        const pointField& point = this->points();
        const labelList& owner   = this->faceOwner();
        const labelList& neighbour = this->faceNeighbour();
    
        scalar minCellVol = cell[0].mag(point,face);
        label minCellInd = 0;
        scalar maxCellVol = cell[0].mag(point,face);
        label maxCellInd = 0;
        scalar CellVolAvg = 0;
    
        scalar vol,neighbourVol;
        for(int i=0;i<cell.size();i++)
        {
            vol = cell[i].mag(point,face);
            CellVolAvg += vol;
            if(vol > maxCellVol)
            {
                maxCellVol = vol;
                maxCellInd = i;
            }
            if(vol < minCellVol)
            {
                minCellVol = vol;
                minCellInd = i;
            }
        }
        CellVolAvg /= cell.size();
        
        Info<<endl;
        Info<<"MinCell "<<minCellInd<<" vol: "<<minCellVol<<endl;
        Info<<"MaxCell "<<maxCellInd<<" vol: "<<maxCellVol<<endl;
        Info<<"Average cell vol "<<CellVolAvg<<endl;
        scalar MAXCELLVOL = maxCellVol;
        minCellVol = cell[0].mag(point,face);
        minCellInd = 0;
        maxCellVol = cell[0].mag(point,face);
        maxCellInd = 0;
        CellVolAvg = 0;
        Info<<"di dumm"<<endl;
        //label mergeFace,mergeCell;
        for(int i=0;i<cell.size();i++)
        {
            if(mergeFaceOfCell[i][0] < 0)
                continue;

            Info<<"i:"<<i<<endl;
            vol = cell[i].mag(points,faces);
            if((vol*factor) < maxCellVol)
            {
                std::unordered_multiset<label> cellSet;
                DynamicList<label> allCells;
                for(int s=0;s<mergeFaceOfCell[i].size();s++)
                {
                    selectedFace = mergeFaceOfCell[i][s];
                    if(neighbour.size() <= selectedFace)
                    {
                        FatalErrorInFunction
                        << "Merge via boundary face not possible "
                        << exit(FatalError);
                    }
                    if(cellSet.count(neighbour[selectedFace]) == 0)
                        allCells.append(neighbour[selectedFace]);
                    if(cellSet.count(owner[selectedFace]) == 0)
                        allCells.append(owner[selectedFace]);
                    cellSet.insert(neighbour[selectedFace]);
                    cellSet.insert(owner[selectedFace]);
                }
                Info<<">->"<<endl;
                neighbourVol = 0;
                Info<<mergeFaceOfCell[i]<<endl;
                Info<<allCells<<endl;
                Info<<cell.size()<<endl;
                for(int s=0;s<allCells.size();s++)
                {
                    neighbourVol += cell[allCells[s]].mag(points,faces);
                }
                Info<<"<ö.a"<<endl;
                
                if((neighbourVol+vol)*factor < MAXCELLVOL)
                {
                    Info<<endl;
                    Info<<"Cell "<<i<<" merged with "<<allCells<<" but not sufficient"<<
                    " because vol cell:"<<vol<<" and vol neighbor:"<<neighbourVol<<endl;
                    FatalErrorInFunction
                    << "Not sufficient merge data"
                    << exit(FatalError);
                }
                CellVolAvg += neighbourVol+vol;
                if(neighbourVol+vol > maxCellVol)
                {
                    maxCellVol = neighbourVol+vol;
                    maxCellInd = i;
                }
                if(neighbourVol+vol < minCellVol)
                {
                    minCellVol = neighbourVol+vol;
                    minCellInd = i;
                }
            }
            else
            {
                CellVolAvg += vol;
                if(vol > maxCellVol)
                {
                    maxCellVol = vol;
                    maxCellInd = i;
                }
                if(vol < minCellVol)
                {
                    minCellVol = vol;
                    minCellInd = i;
                }
            }
        }
        CellVolAvg /= cell.size();
        Info<<"MinCell "<<minCellInd<<" vol: "<<minCellVol<<endl;
        Info<<"MaxCell "<<maxCellInd<<" vol: "<<maxCellVol<<endl;
        Info<<"Average cell vol "<<CellVolAvg<<endl;
    }
//End:TestSection
    Info<<"--"<<endl;
    t2 = std::chrono::high_resolution_clock::now();
    Info<<"---"<<endl;
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<<"-"<<endl;
    Info<< "took \t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Test for -1 merging selection ";
    t1 = std::chrono::high_resolution_clock::now();
    
    // Remove function
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
    
    /*
    FatalErrorInFunction
    << "Temporary stop"
    << exit(FatalError);
    */
    
    if(mergeFaceOfCell.size() == 0)
    {
        FatalErrorInFunction
        << "Agglomeration cell not found for all cells!"
        << abort(FatalError);  
    }
    if(mergeFaceOfCell.size() != newCells.size())
    {
        FatalErrorInFunction
        << "Agglomeration cell list size unequal to cell list size!"
        << abort(FatalError);  
    }    
    
    Info<<"Remove merged cell from list ";
    t1 = std::chrono::high_resolution_clock::now();
    // Remove agglomerated cell with too low volume for merging
    faceList newFaces_ = faces;
    labelList newOwner_ = owner;
    labelList newNeighbour_ = neighbour;
    Info<<"-----------------------_"<<endl;
    /*
    for(int i=0;i<newFaces_.size();i++)
    {
        Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
        if(i < neighbour.size())
            Info<<" Neighbor:"<<neighbour[i]<<" ";
        for(int k=0;k<newFaces_[i].size();k++)
        {
            Info<<points[newFaces_[i][k]]<<"->";
        }
        Info<<" with centre:"<<newFaces_[i].centre(points);
        Info<<" and normal vector:"<<newFaces_[i].normal(points);
        Info<<" and area:"<<newFaces_[i].mag(points)<<endl;
        
        if(owner[i] == -1)
        {
            FatalErrorInFunction
            << "Cell owner is -1"
            << exit(FatalError);
        }
    }
    */
    
    DynamicList<label> oldCellNumToNewCellNum;
    oldCellNumToNewCellNum.setSize(newCells.size(),-1);
    
    labelList deletedCells(newCells.size());
    for(int i=0;i<deletedCells.size();i++)
    {
        deletedCells[i] = 0;
    }
    int countDeleteFaces = 0;
    for(int i=0;i<newCells.size();i++)
    {
        //Info<<"newCell:"<<i<<endl;
        //if(mergeFaceOfCell[i] != -1 && mergeFaceOfCell[i] != -2)
        if(mergeFaceOfCell[i].size() > 1)
        {
            std::unordered_multiset<label> cellSet;
            DynamicList<label> mergeCells;
            for(int s=0;s<mergeFaceOfCell[i].size();s++)
            {
                selectedFace = mergeFaceOfCell[i][s];
                if(selectedFace >= neighbour.size())
                {
                    Info<<"selectedFace:"<<selectedFace<<endl;
                    Info<<"Merging faces of cell: "<<i<<" :"<<mergeFaceOfCell[i]<<endl;
                    Info<<"neighbour.size()="<<neighbour.size()<<endl;
                    Info<<"Bool:"<<(neighbour.size() >= selectedFace)<<endl;
                    FatalErrorInFunction
                    << "Merge via boundary face not possible "
                    << exit(FatalError);
                }
                if(cellSet.count(neighbour[selectedFace]) == 0 && neighbour[selectedFace] != i)
                    mergeCells.append(neighbour[selectedFace]);
                if(cellSet.count(owner[selectedFace]) == 0 && owner[selectedFace] != i)
                    mergeCells.append(owner[selectedFace]);
                cellSet.insert(neighbour[selectedFace]);
                cellSet.insert(owner[selectedFace]);
            }
            
            countDeleteFaces += mergeFaceOfCell[i].size();
            label myCell = i;
            
            //label viaFace = mergeFaceOfCell[i];
            /*
            label myCell = i;
            label mergedWithCell = -1;
            if(owner[mergeFaceOfCell[i]] == i)
                mergedWithCell = neighbour[mergeFaceOfCell[i]];
            else if(neighbour[mergeFaceOfCell[i]] == i)
                mergedWithCell = owner[mergeFaceOfCell[i]];
            else
            {
                FatalErrorInFunction
                << "No Merge Cell found. That can not happen!"
                << exit(FatalError);
            }
            */
            if(oldCellNumToNewCellNum[i] != -1)
            {
                FatalErrorInFunction
                << "oldCell to new Cell already taken: own"
                << exit(FatalError);
            }
            oldCellNumToNewCellNum[i] = i;
            for(int s=0;s<mergeCells.size();s++)
            {
                if(oldCellNumToNewCellNum[mergeCells[s]] != -1 && 
                   oldCellNumToNewCellNum[mergeCells[s]] != mergeCells[s])
                {
                    FatalErrorInFunction
                    << "oldCell already used to new Cell already taken: own"
                    << exit(FatalError);
                }
                oldCellNumToNewCellNum[mergeCells[s]] = i;
            
                if(deletedCells[mergeCells[s]] == 1)
                {
                    FatalErrorInFunction
                    << "Cell is multiple times deleted!"
                    << exit(FatalError);
                }
                else
                {
                    deletedCells[mergeCells[s]] = 1;
                }
            }

            // set informations of faces to delete to -1
            for(int s=0;s<mergeFaceOfCell[i].size();s++)
            {
                newFaces_[mergeFaceOfCell[i][s]] = face(); // Only viable if empty face is doable
                if(newOwner_[mergeFaceOfCell[i][s]] == -1)
                {
                    FatalErrorInFunction
                    << "Deletion Face has already owner -1!"
                    << exit(FatalError);
                }
                newOwner_[mergeFaceOfCell[i][s]] = -1;
                if(mergeFaceOfCell[i][s] >= newNeighbour_.size())
                {
                    FatalErrorInFunction
                    << "Merge Face is has no neighbour that can not happen!"
                    << exit(FatalError);
                }
                newNeighbour_[mergeFaceOfCell[i][s]] = -1;
            }
            
            std::unordered_set<label> mergeFaces;
            for(int s=0;s<mergeFaceOfCell[i].size();s++)
            {
                mergeFaces.insert(mergeFaceOfCell[i][s]);
            }
            // Renumber faces of main Cell
            labelList facesMyCell = newCells[myCell];
            for(int k=0;k<facesMyCell.size();k++)
            {
                if(mergeFaces.count(facesMyCell[k]) != 0)
                    continue;

                if(owner[facesMyCell[k]] == myCell)
                    newOwner_[facesMyCell[k]] = myCell;
                else if(neighbour[facesMyCell[k]] == myCell)
                    newNeighbour_[facesMyCell[k]] = myCell;
                else
                {
                    FatalErrorInFunction
                    << "Face of merging cells is neither in owner nor in neighbour cell!"
                    << exit(FatalError);
                }
            }
            
            for(int s=0;s<mergeCells.size();s++)
            {
                labelList facesMergeCell = newCells[mergeCells[s]];
                for(int k=0;k<facesMergeCell.size();k++)
                {
                    if(mergeFaces.count(facesMergeCell[k]) != 0)
                        continue;
                
                    if(owner[facesMergeCell[k]] == mergeCells[s])
                        newOwner_[facesMergeCell[k]] = myCell;
                    else if(neighbour[facesMergeCell[k]] == mergeCells[s])
                        newNeighbour_[facesMergeCell[k]] = myCell;
                    else
                    {
                        FatalErrorInFunction
                        << "Face of merging cells is neither in owner nor in neighbour cell!"
                        << exit(FatalError);
                    }                    
                }
            }
        }
        else
        {
            if(oldCellNumToNewCellNum[i] == -1)
                oldCellNumToNewCellNum[i] = i;
        }
    }
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
    /*
    for(int i=0;i<newFaces_.size();i++)
    {
        Info<<"Face:"<<i<<" Owner:"<<newOwner_[i]<<" ";
        if(i < newNeighbour_.size())
            Info<<" Neighbor:"<<newNeighbour_[i]<<" ";
        for(int k=0;k<newFaces_[i].size();k++)
        {
            Info<<points[newFaces_[i][k]]<<"->";
        }
        if(newFaces_[i].size() > 0)
        {
            Info<<" with centre:"<<newFaces_[i].centre(points);
            Info<<" and normal vector:"<<newFaces_[i].normal(points);
            Info<<" and area:"<<newFaces_[i].mag(points)<<endl;
        }
        else
            Info<<endl;
    }
    */
    
    Info<<"Recompute faces, owner, neighbor ";
    t1 = std::chrono::high_resolution_clock::now();
    
    faceList newFaces__(newFaces_.size()-countDeleteFaces);
    labelList newOwner__(newOwner_.size()-countDeleteFaces);
    labelList newNeighbour__(newNeighbour_.size()-countDeleteFaces);
    //Info<<"Created Data Struc"<<endl;

    
    int countDel = 0;
    for(int i=0;i<newFaces_.size();i++)
        if(newOwner_[i] == -1)
            countDel++;
        
    if(countDel != countDeleteFaces)
    {
        FatalErrorInFunction
        << countDel<<"!="<<countDeleteFaces
        << exit(FatalError);
    }
    
    int insertCounter = 0;
    for(int i = 0;i<newFaces_.size();i++)
    {
        /*
        Info<<"Move "<<i<<" of "<<newFaces_.size()<<" to "<<insertCounter
        <<"/"<<newFaces_.size()-countDeleteFaces<<
        "-"<<newOwner_.size()-countDeleteFaces<<endl;
        */
        if(newOwner_[i] != -1)
        {
            //Info<<"0"<<endl;
            newFaces__[insertCounter] = newFaces_[i];
            //Info<<"1"<<endl;
            newOwner__[insertCounter] = newOwner_[i];
            //Info<<"2"<<endl;
            if(newOwner__[insertCounter] == -1)
            {            
                FatalErrorInFunction
                << "newOwner["<<insertCounter<<"]: "<<newOwner__[insertCounter]
                <<" from "<<newOwner_[i]
                << exit(FatalError);
            }
            //Info<<"3"<<endl;
            if(i < newNeighbour_.size())
                newNeighbour__[insertCounter] = newNeighbour_[i];
            insertCounter++;
        }
    }
    
    /*
    Info<<"End"<<endl;
    for(int i=0;i<newFaces__.size();i++)
    {
        Info<<"Face:"<<i<<" Owner:"<<newOwner__[i]<<" ";
        if(i < newNeighbour__.size())
            Info<<" Neighbor:"<<newNeighbour__[i]<<" ";
        for(int k=0;k<newFaces__[i].size();k++)
        {
            Info<<points[newFaces__[i][k]]<<"->";
        }
        if(newFaces__[i].size() > 0)
        {
            Info<<" with centre:"<<newFaces__[i].centre(points);
            Info<<" and normal vector:"<<newFaces__[i].normal(points);
            Info<<" and area:"<<newFaces__[i].mag(points)<<endl;
        }
        else
            Info<<endl;
        if(newOwner__[i] == -1)
        {            
            FatalErrorInFunction
            << "Stop."
            << exit(FatalError);
        }
    }
    */
    
    const polyBoundaryMesh& boundMesh = this->boundaryMesh();
    patchStarts = labelList(boundMesh.size());
    patchSizes = labelList(boundMesh.size());
    for(int i=0;i<boundMesh.size();i++)
    {
        patchStarts[i] = boundMesh[i].start()-countDeleteFaces;
        patchSizes[i] = boundMesh[i].faceCentres().size();
    }
    
    
    labelList cellReductionNumb(newCells.size());
    label count = 0;
    for(int i=0;i<cellReductionNumb.size();i++)
    {
        if(deletedCells[i] == 1)
        {
            count++;
            cellReductionNumb[i] = -1;
        }
        else
        {
            //Info<<"One none deleted"<<endl;
            cellReductionNumb[i] = count;
        }
    }
    label delNum = 0;
    for(int i=0;i<oldCellNumToNewCellNum.size();i++)
    {
        if(cellReductionNumb[i] != -1)
        {
            oldCellNumToNewCellNum[i] -= cellReductionNumb[i];
        }
        else
        {
            oldCellNumToNewCellNum[i] -= cellReductionNumb[oldCellNumToNewCellNum[i]];
            delNum++;
        }
    }
    Info<<"delNum:"<<delNum<<endl;
    label numDeletedCells = 0;
    for(int i=0;i<deletedCells.size();i++)
    {
        if(deletedCells[i] == 1)
            numDeletedCells++;
    }
    Info<<"numDeletedCells: "<<numDeletedCells<<endl;
    
    for(int i=0;i<newNeighbour__.size();i++)
    {
        if(cellReductionNumb[newOwner__[i]] != -1 &&
           cellReductionNumb[newNeighbour__[i]] != -1)
        {
            int temp = newOwner__[i];
            newOwner__[i] -= cellReductionNumb[newOwner__[i]];
            if(newOwner__[i] == -1)
            {
                FatalErrorInFunction
                << "Owner original "<<temp<<endl
                << "newOwner "<<newOwner__[i]<<endl
                << "because of reduction "<<cellReductionNumb[temp]
                << exit(FatalError);
            }
            newNeighbour__[i] -= cellReductionNumb[newNeighbour__[i]];
        }
        else
        {
            /*
            Info<<"newOwner__["<<i<<"]:"<<cellReductionNumb[newOwner__[i]]<<endl
            <<"newNeighbour__["<<i<<"]:"<<cellReductionNumb[newNeighbour__[i]]<<endl;
            */
            
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
    for(int i=newNeighbour__.size();i<newFaces__.size();i++)
    {
        if(cellReductionNumb[newOwner__[i]] != -1)
        {
            int temp = newOwner__[i];
            newOwner__[i] -= cellReductionNumb[newOwner__[i]];
            if(newOwner__[i] == -1)
            {
                FatalErrorInFunction
                << "Owner original "<<temp<<endl
                << "newOwner "<<newOwner__[i]<<endl
                << "because of reduction "<<cellReductionNumb[temp]
                << exit(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
    
    Info<<"cellReductionNumb.size() ="<<cellReductionNumb.size()<<endl;
    DynamicList<DynamicList<label>> newCellNumToOldCellNum;
    label numberDeletedCells = -1;
    for(int i=cellReductionNumb.size()-1;i>=0;i--)
    {
        if(cellReductionNumb[i] != -1)
        {
            Info<<"Non deleted:"<<i<<endl;
            numberDeletedCells = cellReductionNumb[i];
            Info<<"numberDeletedCells:"<<numberDeletedCells<<endl;
            break;
        }
    }
    if(numberDeletedCells == -1)
    {
        FatalErrorInFunction
        << "All cells deleted!"
        << exit(FatalError);  
    }
    newCellNumToOldCellNum.setSize(newCells.size()-numberDeletedCells);
    for(int i=0;i<oldCellNumToNewCellNum.size();i++)
    {
        if(oldCellNumToNewCellNum[i] >= newCellNumToOldCellNum.size())
        {
            Info<<"oldCellNumToNewCellNum["<<i<<"]:"<<oldCellNumToNewCellNum[i]<<" newCellNumToOldCellNum.size():"<<newCellNumToOldCellNum.size()<<endl;
            FatalErrorInFunction
            << "Wrong assignment!"
            << exit(FatalError); 
        }
        else
        {
            newCellNumToOldCellNum[oldCellNumToNewCellNum[i]].append(i);
        }
    }
    
    /*
    Info<<"newCell: "<<50491<<" was "<<newCellNumToOldCellNum[50941]<<endl;
    Info<<"mergeFaceOfCell:"<<mergeFaceOfCell[newCellNumToOldCellNum[50941][0]]
    <<" mergeNecessary: "<<mergeNecessary[newCellNumToOldCellNum[50941][0]]<<" with "<<
    "partialVolumeScale: "<<partialVolumeScale[newCellNumToOldCellNum[50941][0]]<<endl;
    */
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Test reset and test ";
    t1 = std::chrono::high_resolution_clock::now();
    
    testNewMeshData(newFaces__,newOwner__,newNeighbour__,patchStarts,patchSizes);
    
    resetPrimitives(Foam::clone(points),
                    Foam::clone(newFaces__),
                    Foam::clone(newOwner__),
                    Foam::clone(newNeighbour__),
                    patchSizes,
                    patchStarts,
                    true);
    
    /*
     * Uncomment later!!!
    testForCellSize
    (
        possibleMergeFaceArea,possibleMergeFaces, possibleMergeCells,
        oneMergeFaceSufficient, mergeNecessary, mergeFaceOfCell,partialThreeshold
    );
    */
    
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t" << time_span.count() << " seconds."<<endl;
}

void Foam::cutCellFvMesh::agglomerateSmallCells_cutNeg_plus
(
    scalarList& newCellVolume,
    scalarList& oldCellVolume,
    scalar partialThreeshold
)
{
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span;
    
    Info<<endl;
    Info<<"Preprocessing of small cells ";
    t1 = std::chrono::high_resolution_clock::now();
    /*
    for(int i=0;i<oldCellVolume.size();i++)
    {
        Info<<"Cell "<<i<<" volume: "<<oldCellVolume[i]<<" splitted to minus side cell: "<<oldSplittedCellToNewMinusCell[i]<<endl;
        Info<<"Cell "<<i<<" splitted to plus side cell: "<<oldSplittedCellToNewPlusCell[i]<<endl;
    }
    */
    scalarList partialVolumeScale = scalarList(newCellVolume.size());
    //Info<<"new Cell Size: "<<newCellVolume.size()<<endl;

    label deletedCellsCount = 0;
    for(int i=0;i<deletedCellsList.size();i++)
    {
        if(deletedCellsList[i] == 1)
            deletedCellsCount++;
    }
    Info<<"1"<<endl;
    if(newCellVolume.size()+deletedCellsCount != oldCellVolume.size())
    {
        FatalErrorInFunction
        << "Must not happen!"
        << exit(FatalError); 
    }
    Info<<"2"<<endl;
    
    
    Info<<endl<<"deletedCellsList: "<<deletedCellsList<<endl;
    Info<<endl<<"oldSplittedCellToNewPlusCell: "<<oldSplittedCellToNewPlusCell<<endl;


    Info<<"3"<<endl;
    
    if(newCellVolume.size() != mapNewCellsToOldCells.size())
        FatalErrorInFunction<< "Must not happen!"<< exit(FatalError);
    
    for(int i=0;i<newCellVolume.size();i++)
    {
        if(mapNewCellsToOldCells[i] == -1)
            partialVolumeScale[i] = 1;
        else
            partialVolumeScale[i] = newCellVolume[i]/oldCellVolume[mapNewCellsToOldCells[i]];
    }
    Info<<"4"<<endl;
    /*
    for(int i=0;i<newCellVolume.size();i++)
    {
        Info<<"new Cell "<<i<<" has volume scale: "<<partialVolumeScale[i]<<endl;
    }
    */
    
    const cellList& newCells = this->cells();
    const faceList& faces = this->faces();
    const labelList& owner   = this->faceOwner();
    const labelList& neighbour = this->faceNeighbour();
    const pointField& points = this->points();
    
    DynamicList<DynamicList<DynamicList<label>>> possibleMergeFaces;
    possibleMergeFaces.setSize(newCellVolume.size());
    DynamicList<DynamicList<DynamicList<label>>> possibleMergeCells;
    possibleMergeCells.setSize(newCellVolume.size());
    DynamicList<DynamicList<scalar>> possibleMergeCellsPartialSize;
    possibleMergeCellsPartialSize.setSize(newCellVolume.size());
    DynamicList<DynamicList<scalar>> possibleMergeFaceArea;
    possibleMergeFaceArea.setSize(newCellVolume.size());
    DynamicList<DynamicList<bool>> possibleMergeFaceSufficient;
    possibleMergeFaceSufficient.setSize(newCellVolume.size());
    DynamicList<bool> oneMergeFaceSufficient;
    oneMergeFaceSufficient.setSize(newCellVolume.size());
    DynamicList<bool> mergeNecessary;
    mergeNecessary.setSize(newCellVolume.size());
    DynamicList<label> MergeCell;
    label neighbourCell = -1;
    scalar neighbourCellPartialVolume;
    Info<<"Create merge Cells"<<endl;
    for(int i=0;i<newCellVolume.size();i++)
    {
        Info<<i<<endl;
        mergeNecessary[i] = false;
        if((partialVolumeScale[i] < 1) && (partialVolumeScale[i] < partialThreeshold))
        {
            mergeNecessary[i] = true;
            /*
            Info<<"new Cell "<<i<<" is signed for merge via faces ";
            for(int k=0;k<newCells[i].size();k++)
            {
                Info<<newCells[i][k]<<",";
            }
            Info<<" with cells:";
            */
            
            /*
             * Find merging cells for one to one merging
             */
            Info<<"1-1"<<endl;
            for(int k=0;k<newCells[i].size();k++)
            {
                if(newCells[i][k] < neighbour.size())
                {
                    if(owner[newCells[i][k]] == i)
                        neighbourCell = neighbour[newCells[i][k]];
                    else if(neighbour[newCells[i][k]] == i)
                        neighbourCell = owner[newCells[i][k]];
                    else
                        FatalErrorInFunction<<"Agglomeration face does not belong to the agglomerated cell. Something is wrong here!"<< exit(FatalError);
                    
                    neighbourCellPartialVolume = partialVolumeScale[neighbourCell];
                    
                    DynamicList<label> temp;
                    temp.append(newCells[i][k]);
                    possibleMergeFaces[i].append(temp);
                        
                    temp.setSize(0);
                    temp.append(neighbourCell);
                    possibleMergeCells[i].append(temp);

                    possibleMergeFaceArea[i].append(faces[possibleMergeFaces[i][possibleMergeCells[i].size()-1][0]].mag(points));
                    possibleMergeFaceSufficient[i].append(true);
                    possibleMergeCellsPartialSize[i].append(neighbourCellPartialVolume+partialVolumeScale[i]);
                }
            }
            
            //Info<<endl;
            /*
            oneMergeFaceSufficient[i] = false;
            for(int k=0;k<possibleMergeFaceSufficient.size();k++)
            {
                if(possibleMergeFaceSufficient[i][k])
                {
                    oneMergeFaceSufficient[i] = true;
                    break;
                }
            }
            */
            if(possibleMergeCells[i].size() == 0)
            {
                Info<<endl<<"Problem in cell "<<i<<" with partial volume "<< partialVolumeScale[i]<<endl;
                Info<<"\tNeighbours are:"<<endl;
                for(int k=0;k<newCells[i].size();k++)
                {
                    if(newCells[i][k] < neighbour.size())
                    {
                        if(owner[newCells[i][k]] == i)
                        {
                            neighbourCell = neighbour[newCells[i][k]];
                        }
                        else if(neighbour[newCells[i][k]] == i)
                        {   
                            neighbourCell = owner[newCells[i][k]];
                        }
                        neighbourCellPartialVolume = partialVolumeScale[neighbourCell];
                        
                        Info<<"\tCell:"<<neighbourCell<<" partialVol:"<<neighbourCellPartialVolume<<
                            " combinedPartialVol:"<<(neighbourCellPartialVolume + partialVolumeScale[i])<<" with threshold:"<<partialThreeshold
                            <<endl;
                    }
                }
                FatalErrorInFunction
                << "Agglomeration cell not found for all cells!"
                << exit(FatalError);
            }
        }
        /*
        Info<<"cell:"<<i<<" mergeNecessary:"<<mergeNecessary[i]<<" partialVolumeScale:"<<partialVolumeScale[i]<<endl;
        */
    }
    
    Info<<"Created merge Cells --"<<endl;
    for(int i=0;i<newCellVolume.size();i++)
    {
        Info<<"i:"<<i<<endl;
        Info<<"possibleMergeFaces[i]:"<<possibleMergeFaces[i]<<endl;
        for(int j=0;j<possibleMergeFaces[i].size();j++)
        {
            Info<<"j:"<<j<<endl;
            
            label numCellMerge=-1;
            if(possibleMergeFaces[i][j].size() == 1)
                numCellMerge = 2;
            else if(possibleMergeFaces[i][j].size() == 4)
                numCellMerge = 4;
            else if(possibleMergeFaces[i][j].size() == 12)
                numCellMerge = 8;
            else
                FatalErrorInFunction<<"Wrong number of merge faces"<< exit(FatalError);
            
            if(possibleMergeCells[i][j].size()+1 != numCellMerge)
            {
                Info<<"numCellMerge: "<<numCellMerge<<endl;
                Info<<"possibleMergeCells[i][j].size(): "<<possibleMergeCells[i][j].size()<<endl;
                Info<<"possibleMergeFaces[i][j].size(): "<<possibleMergeFaces[i][j].size()<<endl;
                Info<<"possibleMergeCells[i][j]: "<<possibleMergeCells[i][j]<<endl;
                Info<<"possibleMergeFaces[i][j]: "<<possibleMergeFaces[i][j]<<endl;
                FatalErrorInFunction<<"Number of merge cells does not match"<< exit(FatalError);
            }


            std::unordered_multiset<label> cellSet;
            for(int k=0;k<possibleMergeFaces[i][j].size();k++)
            {
                cellSet.insert(owner[possibleMergeFaces[i][j][k]]);
                cellSet.insert(neighbour[possibleMergeFaces[i][j][k]]);
            }
            
            label cellDuplicationNum=-1;
            if(numCellMerge==2) cellDuplicationNum = 1;
            else if(numCellMerge==4) cellDuplicationNum = 2;
            else if(numCellMerge==8) cellDuplicationNum = 3;
            else FatalErrorInFunction<<"Somethings wrong here"<< exit(FatalError);
            
            for(int k=0;k<possibleMergeCells[i][j].size();k++)
            {
                if(cellSet.count(possibleMergeCells[i][j][k])!=static_cast<long unsigned int>(cellDuplicationNum))
                    FatalErrorInFunction<<"Merge Cells and Faces do not match!"<< exit(FatalError);
            }
        }
        Info<<"End"<<endl;
    }
    
    // Sort possible merging cell by respect to face area biggest to smallest
    for(int i=0;i<possibleMergeFaceArea.size();i++)
    {
        Info<<"Sort: "<<i<<endl;
        //Info<<"possibleMergeCellsPartialSize: "<<possibleMergeCellsPartialSize<<endl;
        int j;
        scalar keyArea,keySize;
        DynamicList<label> keyFaces,keyCells;
        bool keySuff;
        for(int k=1;k<possibleMergeFaceArea[i].size();k++)
        {
            keyArea = possibleMergeFaceArea[i][k];
            keyFaces = possibleMergeFaces[i][k];
            keyCells = possibleMergeCells[i][k];
            keySuff = possibleMergeFaceSufficient[i][k];
            keySize = possibleMergeCellsPartialSize[i][k];            
            j = k-1;
            while(j>=0 && possibleMergeFaceArea[i][j] < keyArea)
            {
                possibleMergeFaceArea[i][j+1] = possibleMergeFaceArea[i][j];
                possibleMergeFaces[i][j+1] = possibleMergeFaces[i][j];
                possibleMergeCells[i][j+1] = possibleMergeCells[i][j];
                possibleMergeFaceSufficient[i][j+1] = possibleMergeFaceSufficient[i][j];
                possibleMergeCellsPartialSize[i][j+1] = possibleMergeCellsPartialSize[i][j];
                j--;
            }
            possibleMergeFaceArea[i][j+1] = keyArea;
            possibleMergeFaces[i][j+1] = keyFaces;
            possibleMergeCells[i][j+1] = keyCells;
            possibleMergeFaceSufficient[i][j+1] = keySuff;
            possibleMergeCellsPartialSize[i][j+1] = keySize;
        }
    }
    Info<<"Sorted merge Cells"<<endl;
    
    for(int i=0;i<possibleMergeFaces.size();i++)
    {
        for(int j=1;j<possibleMergeFaces[i].size();j++)
        {
            if(possibleMergeFaces[i][j-1].size() > possibleMergeFaces[i][j].size())
                FatalErrorInFunction<<"Invalid sorting"<<exit(FatalError);
        }
    }
    
//Test for correct merge candidates
    scalar factor = 1/partialThreeshold;
    scalar minCellVol = newCells[0].mag(points,faces);
    label minCellInd = 0;
    scalar maxCellVol = newCells[0].mag(points,faces);
    label maxCellInd = 0;
    scalar CellVolAvg = 0;
    
    scalar vol,neighborVol;
    for(int i=0;i<newCells.size();i++)
    {
        vol = newCells[i].mag(points,faces);
        CellVolAvg += vol;
        if(vol > maxCellVol)
        {
            maxCellVol = vol;
            maxCellInd = i;
        }
        if(vol < minCellVol)
        {
            minCellVol = vol;
            minCellInd = i;
        }
    }
    CellVolAvg /= newCells.size();

    
    Info<<endl<<"Minimum cell "<<minCellInd<<" vol:"<<minCellVol
    <<endl<<"Maximum cell "<<maxCellInd<<" vol:"<<maxCellVol<<endl;
    Info<<" Average vol was:"<<CellVolAvg<<endl;    
    
    for(int i=0;i<newCells.size();i++)
    {
        vol = newCells[i].mag(points,faces); 
        if((vol*factor) < maxCellVol)
        {
            if(mergeNecessary[i] == false)
            {
                FatalErrorInFunction
                << "Not to merge but necessary"
                << exit(FatalError);
            }
            if(possibleMergeCells[i].size() == 0)
            {
                FatalErrorInFunction
                << "No merge data available"
                << exit(FatalError);
            }            
        }
    }
    
    scalar partialVol;
    for(int i=0;i<newCells.size();i++)
    {
        if((partialVolumeScale[i] < 1) && (partialVolumeScale[i] < partialThreeshold))
        {
            if(possibleMergeCells[i].size()==0)
            {
                FatalErrorInFunction
                << "Merge Face "<<i<<" with no merge partners!"
                << exit(FatalError);  
            }
            if( (possibleMergeFaces[i].size()!=possibleMergeCells[i].size())&&
                (possibleMergeCells[i].size()!=possibleMergeFaceArea[i].size()))
            {
                FatalErrorInFunction
                << "Data error"
                << exit(FatalError);  
            }
        }
    }
    
    for(int i=0;i<possibleMergeFaces.size();i++)
    {
        for(int j=0;j<possibleMergeFaces[i].size();j++)
        {
            // Test one option
            std::unordered_multiset<label> optionMergeCellsSet;
            DynamicList<label> optionMergeCellsList;
            for(int k=0;k<possibleMergeFaces[i][j].size();k++)
            {
                if(possibleMergeFaces[i][j][k] >= neighbour.size())
                {
                    FatalErrorInFunction<<"Boundary face is merge face"<<exit(FatalError);
                }
                label ownerCell = owner[possibleMergeFaces[i][j][k]];
                label neighborCell = neighbour[possibleMergeFaces[i][j][k]];
                if(optionMergeCellsSet.find(ownerCell) == optionMergeCellsSet.end())
                    optionMergeCellsList.append(ownerCell);
                if(optionMergeCellsSet.find(neighborCell) == optionMergeCellsSet.end())
                    optionMergeCellsList.append(neighborCell);
                optionMergeCellsSet.insert(ownerCell);
                optionMergeCellsSet.insert(neighborCell);
            }
            if(optionMergeCellsList.size()<=0)
                FatalErrorInFunction<<"Merging option with no cells!"<<exit(FatalError);
            
            label mergeCellMult = optionMergeCellsSet.count(optionMergeCellsList[0]);
            for(int k=1;k<optionMergeCellsList.size();k++)
            {
                if(optionMergeCellsSet.count(optionMergeCellsList[k]) != static_cast<long unsigned int>(mergeCellMult))
                {
                    FatalErrorInFunction<<"Different multiplcity of cells in "<<
                    optionMergeCellsList.size()<<" cell merging "<<exit(FatalError);
                }
            }
            
            if(mergeCellMult==1 && optionMergeCellsList.size() == 2)
            {
                //2 Cell merging
                if(possibleMergeFaces[i][j].size() != 1)
                    FatalErrorInFunction<<"Error in 2 Cell merging! "<<exit(FatalError);
            }
            else if(mergeCellMult==2 && optionMergeCellsList.size() == 4)
            {
                //4 Cell merging
                if(possibleMergeFaces[i][j].size() != 4)
                    FatalErrorInFunction<<"Error in 4 Cell merging! "<<exit(FatalError);
            }
            else if(mergeCellMult==3 && optionMergeCellsList.size() == 8)
            {
                //4 Cell merging
                if(possibleMergeFaces[i][j].size() != 12)
                    FatalErrorInFunction<<"Error in 8 Cell merging! "<<exit(FatalError);
            }
            else
            {
                Info<<"mergeCellMult: "<<mergeCellMult<<endl;
                Info<<"optionMergeCellsList.size() == "<<optionMergeCellsList.size()<<endl;
                Info<<"possibleMergeFaces[i][j].size() == "<<possibleMergeFaces[i][j].size()<<endl;
                FatalErrorInFunction<<"Inconsistent merge Option!"<<exit(FatalError);
            }
            std::unordered_set<label> avoidCellDuplicates;
            for(int k=0;k<possibleMergeCells[i][j].size();k++)
            {
                if(avoidCellDuplicates.find(possibleMergeCells[i][j][k]) == avoidCellDuplicates.end())
                    avoidCellDuplicates.insert(possibleMergeCells[i][j][k]);
                else
                    FatalErrorInFunction<<"Duplicate cell in Merge option!"<<exit(FatalError);
            }
            for(int k=0;k<possibleMergeCells[i][j].size();k++)
            {
                if(optionMergeCellsSet.count(possibleMergeCells[i][j][k]) != static_cast<long unsigned int>(mergeCellMult))
                    FatalErrorInFunction<<"Non matching mergeCell in Merge option!"<<exit(FatalError);
            }
        }
    }
//End: Test for correct merge candidates

    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Processing of small cells ";
    t1 = std::chrono::high_resolution_clock::now(); 

    /*
    std::unordered_set<label> cellReserved;
    DynamicList<DynamicList<label>> blockedCells;
    blockedCells.setSize(possibleMergeCells.size());
    
    
    labelList mergeFaceOfCell = searchDown_rec(possibleMergeFaceArea,possibleMergeFaces,possibleMergeCells,oneMergeFaceSufficient,mergeNecessary,0,blockedCells,cellReserved);
    */

    
    List<DynamicList<label>> mergeFaceOfCell = assignMergeFaces(owner,neighbour,
                                                possibleMergeFaceArea, possibleMergeFaces,
                                                possibleMergeCells,oneMergeFaceSufficient,
                                                mergeNecessary,possibleMergeCellsPartialSize,
                                                partialThreeshold);
    
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Test for duplicate merging selection ";
    t1 = std::chrono::high_resolution_clock::now();
    
    std::unordered_set<label> usedFace;
    for(int i=0;i<mergeFaceOfCell.size();i++)
    {
        if(mergeFaceOfCell[i].size() == 0)
        {
            FatalErrorInFunction
            << "Cell with unwanted results"
            << exit(FatalError);  
        }
        if(mergeFaceOfCell[i][0] == -3)
        {
            FatalErrorInFunction
            << "Too small cell was not treated by backtracking algorithm"
            << exit(FatalError);  
        }
        if(mergeFaceOfCell[i][0] == -4)
        {
            FatalErrorInFunction
            << "Cell with unwanted results"
            << exit(FatalError);  
        }
        if(mergeFaceOfCell[i][0] < 0)
            continue;
        
        for(int s=0;s<mergeFaceOfCell[i].size();s++)
        {
            if(usedFace.find(mergeFaceOfCell[i][s]) == usedFace.end())
                usedFace.insert(mergeFaceOfCell[i][s]);
            else
            {
                Info<<endl<<endl;
                label fc = mergeFaceOfCell[i][s];
                label wnr = owner[fc];
                label nghbr = neighbour[fc];
            
                Info<<"Face "<<fc<<" merging of cell:"<<wnr<<" with Vol:"<<newCellVolume[wnr]<<
                " cell:"<<nghbr<<" with Vol:"<<newCellVolume[nghbr];
                FatalErrorInFunction
                << "Merge Face used twice!"
                << exit(FatalError);  
            }
        }
    }

    label selectedFace;
    std::unordered_set<label> usedCells;
    for(int i=0;i<mergeFaceOfCell.size();i++)
    {
        for(int s=0;s<mergeFaceOfCell[i].size();s++)
        {
            selectedFace = mergeFaceOfCell[i][s];
            if(selectedFace == -3)
            {
                Info<<"Face is: "<<selectedFace<<endl;
            
                FatalErrorInFunction
                << "One cell not treated by algorithm"
                << exit(FatalError);  
            }
        
            // Test that each cell to merge is merged
            if(selectedFace == -1 && mergeNecessary[i])
            {
                Info<<"Selected Face is: "<<selectedFace<<" but mergeNecessary["<<i<<"]"<<mergeNecessary[i]<<endl;
            
                FatalErrorInFunction
                << "Agglomeration cell not found but necessary!"
                << exit(FatalError);  
            }
            
            if((selectedFace == -1 || selectedFace == -2) && mergeFaceOfCell[i].size() != 1)
            {
                FatalErrorInFunction
                << "Not merged or merged by other cell but size of mergeList is not one!"
                << exit(FatalError);  
            }
        }
        
        if(mergeFaceOfCell[i][0] == -1 || mergeFaceOfCell[i][0] == -2)
        {
            continue;
        }
        
        for(int s=0;s<mergeFaceOfCell[i].size();s++)
        {
            selectedFace = mergeFaceOfCell[i][s];
            //Test that merge faces are legit faces
            if(selectedFace > neighbour.size() || selectedFace < 0)
            {
                Info<<endl<<"Merging Face: "<<selectedFace<<endl;
                FatalErrorInFunction
                << "Merging face is not an existing face!"
                << exit(FatalError);
            }
        }
        
        bool selectedFaceExists = false;
        for(int s=0;s<possibleMergeFaces[i].size();s++)
        {
            if(possibleMergeFaces[i][s] == mergeFaceOfCell[i])
                selectedFaceExists = true;
        }
        if(!selectedFaceExists)
        {
            FatalErrorInFunction
            << "Merging face is not in the mergingFace list!"
            << exit(FatalError);
        }
        
        label numMergFaces = mergeFaceOfCell[i].size();
        if(numMergFaces!=1 && numMergFaces!=4 && numMergFaces!=12)
        {
            Info<<endl<<"numMergFaces["<<i<<"]:"<<numMergFaces<<endl;
            FatalErrorInFunction
            << "Number of merging faces does not match!"
            << exit(FatalError);
        }
        
        //Test that merging faces are the correct ones
        std::unordered_multiset<label> cellSet;
        DynamicList<label> allCells;
        for(int s=0;s<mergeFaceOfCell[i].size();s++)
        {
            selectedFace = mergeFaceOfCell[i][s];
            if(neighbour.size() <= selectedFace)
            {
                Info<<"mergeFaceOfCell["<<i<<"]["<<s<<"]="<<mergeFaceOfCell[i][s]<<endl;
                Info<<"neighbour.size():"<<neighbour.size()<<endl;
                FatalErrorInFunction
                << "Merge via boundary face not possible "
                << exit(FatalError);
            }
            if(cellSet.count(neighbour[selectedFace]) == 0)
                allCells.append(neighbour[selectedFace]);
            if(cellSet.count(owner[selectedFace]) == 0)
                allCells.append(owner[selectedFace]);
            cellSet.insert(neighbour[selectedFace]);
            cellSet.insert(owner[selectedFace]);
        }
        
        if(allCells.size()!=2 && allCells.size()!=4 && allCells.size()!=8)
        {
            Info<<endl<<allCells<<endl;
            FatalErrorInFunction
            << "Number of merging cells does not match!"
            << exit(FatalError);
        }
        
        label cellMult=-1;
        if(numMergFaces==1)
        {
            cellMult = 1; //two cell merge
        }
        else if(numMergFaces==4)
        {
            cellMult = 2; //four cell merge
        }
        else if(numMergFaces==12)
        {
            cellMult = 3; //eight cell merge
        }
        else
            FatalErrorInFunction<<"Number of merge faces is wrong: "<<numMergFaces<<"! "<< exit(FatalError);
        
        for(int s=0;s<allCells.size();s++)
        {
            if(cellSet.count(allCells[s])!=static_cast<long unsigned int>(cellMult))
            {
                FatalErrorInFunction
                << "Wrong number of cell count "
                << exit(FatalError);
            }
        }
        
        Info<<endl;
        bool cellsMatch = false;
        bool partMatch;
        for(int s=0;s<possibleMergeCells[i].size();s++)
        {
            partMatch = false;
            if(possibleMergeCells[i][s].size() == allCells.size()-1)
            {
                partMatch = true;
                bool match=false;
                for(int w=0;w<allCells.size();w++)
                {
                    Info<<allCells[w]<<"|"<<match<<"|"<<partMatch<<endl;
                    match = false;
                    if(allCells[w] == i)
                    {
                        match = true;
                        continue;
                    }
                    for(int x=0;x<possibleMergeCells[i][s].size();x++)
                    {
                        if(possibleMergeCells[i][s][x] == allCells[w])
                            match = true;                            
                    }
                    if(match == false)
                        partMatch = false;
                }
            }
            if(partMatch)
                cellsMatch = true;
        }
        if(!cellsMatch)
        {
            Info<<endl<<"i:"<<i<<endl;
            Info<<"possibleMergeCells["<<i<<"]:"<<possibleMergeCells[i]<<endl;
            Info<<"allCells["<<i<<"]:"<<allCells<<endl;
            FatalErrorInFunction
            << "Cells do not match "
            << exit(FatalError);
        }
    }
    
//TestSection
    {
        scalar factor = 1/partialThreeshold;
        const cellList& cell = this->cells();
        const faceList& face = this->faces();
        const pointField& point = this->points();
        const labelList& owner   = this->faceOwner();
        const labelList& neighbour = this->faceNeighbour();
    
        scalar minCellVol = cell[0].mag(point,face);
        label minCellInd = 0;
        scalar maxCellVol = cell[0].mag(point,face);
        label maxCellInd = 0;
        scalar CellVolAvg = 0;
    
        scalar vol,neighbourVol;
        for(int i=0;i<cell.size();i++)
        {
            vol = cell[i].mag(point,face);
            CellVolAvg += vol;
            if(vol > maxCellVol)
            {
                maxCellVol = vol;
                maxCellInd = i;
            }
            if(vol < minCellVol)
            {
                minCellVol = vol;
                minCellInd = i;
            }
        }
        CellVolAvg /= cell.size();
        
        Info<<endl;
        Info<<"MinCell "<<minCellInd<<" vol: "<<minCellVol<<endl;
        Info<<"MaxCell "<<maxCellInd<<" vol: "<<maxCellVol<<endl;
        Info<<"Average cell vol "<<CellVolAvg<<endl;
        scalar MAXCELLVOL = maxCellVol;
        minCellVol = cell[0].mag(point,face);
        minCellInd = 0;
        maxCellVol = cell[0].mag(point,face);
        maxCellInd = 0;
        CellVolAvg = 0;
        Info<<"di dumm"<<endl;
        //label mergeFace,mergeCell;
        for(int i=0;i<cell.size();i++)
        {
            if(mergeFaceOfCell[i][0] < 0)
                continue;

            Info<<"i:"<<i<<endl;
            vol = cell[i].mag(points,faces);
            if((vol*factor) < maxCellVol)
            {
                std::unordered_multiset<label> cellSet;
                DynamicList<label> allCells;
                for(int s=0;s<mergeFaceOfCell[i].size();s++)
                {
                    selectedFace = mergeFaceOfCell[i][s];
                    if(neighbour.size() <= selectedFace)
                    {
                        FatalErrorInFunction
                        << "Merge via boundary face not possible "
                        << exit(FatalError);
                    }
                    if(cellSet.count(neighbour[selectedFace]) == 0)
                        allCells.append(neighbour[selectedFace]);
                    if(cellSet.count(owner[selectedFace]) == 0)
                        allCells.append(owner[selectedFace]);
                    cellSet.insert(neighbour[selectedFace]);
                    cellSet.insert(owner[selectedFace]);
                }
                Info<<">->"<<endl;
                neighbourVol = 0;
                Info<<mergeFaceOfCell[i]<<endl;
                Info<<allCells<<endl;
                Info<<cell.size()<<endl;
                for(int s=0;s<allCells.size();s++)
                {
                    neighbourVol += cell[allCells[s]].mag(points,faces);
                }

                CellVolAvg += neighbourVol+vol;
                if(neighbourVol+vol > maxCellVol)
                {
                    maxCellVol = neighbourVol+vol;
                    maxCellInd = i;
                }
                if(neighbourVol+vol < minCellVol)
                {
                    minCellVol = neighbourVol+vol;
                    minCellInd = i;
                }
            }
            else
            {
                CellVolAvg += vol;
                if(vol > maxCellVol)
                {
                    maxCellVol = vol;
                    maxCellInd = i;
                }
                if(vol < minCellVol)
                {
                    minCellVol = vol;
                    minCellInd = i;
                }
            }
        }
        CellVolAvg /= cell.size();
        Info<<"MinCell "<<minCellInd<<" vol: "<<minCellVol<<endl;
        Info<<"MaxCell "<<maxCellInd<<" vol: "<<maxCellVol<<endl;
        Info<<"Average cell vol "<<CellVolAvg<<endl;
    }
//End:TestSection
    Info<<"--"<<endl;
    t2 = std::chrono::high_resolution_clock::now();
    Info<<"---"<<endl;
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<<"-"<<endl;
    Info<< "took \t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Test for -1 merging selection ";
    t1 = std::chrono::high_resolution_clock::now();
    
    // Remove function
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
    
    /*
    FatalErrorInFunction
    << "Temporary stop"
    << exit(FatalError);
    */
    
    if(mergeFaceOfCell.size() == 0)
    {
        FatalErrorInFunction
        << "Agglomeration cell not found for all cells!"
        << abort(FatalError);  
    }
    if(mergeFaceOfCell.size() != newCells.size())
    {
        FatalErrorInFunction
        << "Agglomeration cell list size unequal to cell list size!"
        << abort(FatalError);  
    }    
    
    Info<<"Remove merged cell from list ";
    t1 = std::chrono::high_resolution_clock::now();
    // Remove agglomerated cell with too low volume for merging
    faceList newFaces_ = faces;
    labelList newOwner_ = owner;
    labelList newNeighbour_ = neighbour;
    Info<<"-----------------------_"<<endl;
    /*
    for(int i=0;i<newFaces_.size();i++)
    {
        Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" ";
        if(i < neighbour.size())
            Info<<" Neighbor:"<<neighbour[i]<<" ";
        for(int k=0;k<newFaces_[i].size();k++)
        {
            Info<<points[newFaces_[i][k]]<<"->";
        }
        Info<<" with centre:"<<newFaces_[i].centre(points);
        Info<<" and normal vector:"<<newFaces_[i].normal(points);
        Info<<" and area:"<<newFaces_[i].mag(points)<<endl;
        
        if(owner[i] == -1)
        {
            FatalErrorInFunction
            << "Cell owner is -1"
            << exit(FatalError);
        }
    }
    */
    
    DynamicList<label> oldCellNumToNewCellNum;
    oldCellNumToNewCellNum.setSize(newCells.size(),-1);
    
    labelList deletedCells(newCells.size());
    for(int i=0;i<deletedCells.size();i++)
    {
        deletedCells[i] = 0;
    }
    int countDeleteFaces = 0;
    for(int i=0;i<newCells.size();i++)
    {
        //Info<<"newCell:"<<i<<endl;
        //if(mergeFaceOfCell[i] != -1 && mergeFaceOfCell[i] != -2)
        if(mergeFaceOfCell[i].size() > 1)
        {
            std::unordered_multiset<label> cellSet;
            DynamicList<label> mergeCells;
            for(int s=0;s<mergeFaceOfCell[i].size();s++)
            {
                selectedFace = mergeFaceOfCell[i][s];
                if(selectedFace >= neighbour.size())
                {
                    Info<<"selectedFace:"<<selectedFace<<endl;
                    Info<<"Merging faces of cell: "<<i<<" :"<<mergeFaceOfCell[i]<<endl;
                    Info<<"neighbour.size()="<<neighbour.size()<<endl;
                    Info<<"Bool:"<<(neighbour.size() >= selectedFace)<<endl;
                    FatalErrorInFunction
                    << "Merge via boundary face not possible "
                    << exit(FatalError);
                }
                if(cellSet.count(neighbour[selectedFace]) == 0 && neighbour[selectedFace] != i)
                    mergeCells.append(neighbour[selectedFace]);
                if(cellSet.count(owner[selectedFace]) == 0 && owner[selectedFace] != i)
                    mergeCells.append(owner[selectedFace]);
                cellSet.insert(neighbour[selectedFace]);
                cellSet.insert(owner[selectedFace]);
            }
            
            countDeleteFaces += mergeFaceOfCell[i].size();
            label myCell = i;
            
            //label viaFace = mergeFaceOfCell[i];
            /*
            label myCell = i;
            label mergedWithCell = -1;
            if(owner[mergeFaceOfCell[i]] == i)
                mergedWithCell = neighbour[mergeFaceOfCell[i]];
            else if(neighbour[mergeFaceOfCell[i]] == i)
                mergedWithCell = owner[mergeFaceOfCell[i]];
            else
            {
                FatalErrorInFunction
                << "No Merge Cell found. That can not happen!"
                << exit(FatalError);
            }
            */
            if(oldCellNumToNewCellNum[i] != -1)
            {
                FatalErrorInFunction
                << "oldCell to new Cell already taken: own"
                << exit(FatalError);
            }
            oldCellNumToNewCellNum[i] = i;
            for(int s=0;s<mergeCells.size();s++)
            {
                if(oldCellNumToNewCellNum[mergeCells[s]] != -1 && 
                   oldCellNumToNewCellNum[mergeCells[s]] != mergeCells[s])
                {
                    FatalErrorInFunction
                    << "oldCell already used to new Cell already taken: own"
                    << exit(FatalError);
                }
                oldCellNumToNewCellNum[mergeCells[s]] = i;
            
                if(deletedCells[mergeCells[s]] == 1)
                {
                    FatalErrorInFunction
                    << "Cell is multiple times deleted!"
                    << exit(FatalError);
                }
                else
                {
                    deletedCells[mergeCells[s]] = 1;
                }
            }

            // set informations of faces to delete to -1
            for(int s=0;s<mergeFaceOfCell[i].size();s++)
            {
                newFaces_[mergeFaceOfCell[i][s]] = face(); // Only viable if empty face is doable
                if(newOwner_[mergeFaceOfCell[i][s]] == -1)
                {
                    FatalErrorInFunction
                    << "Deletion Face has already owner -1!"
                    << exit(FatalError);
                }
                newOwner_[mergeFaceOfCell[i][s]] = -1;
                if(mergeFaceOfCell[i][s] >= newNeighbour_.size())
                {
                    FatalErrorInFunction
                    << "Merge Face is has no neighbour that can not happen!"
                    << exit(FatalError);
                }
                newNeighbour_[mergeFaceOfCell[i][s]] = -1;
            }
            
            std::unordered_set<label> mergeFaces;
            for(int s=0;s<mergeFaceOfCell[i].size();s++)
            {
                mergeFaces.insert(mergeFaceOfCell[i][s]);
            }
            // Renumber faces of main Cell
            labelList facesMyCell = newCells[myCell];
            for(int k=0;k<facesMyCell.size();k++)
            {
                if(mergeFaces.count(facesMyCell[k]) != 0)
                    continue;

                if(owner[facesMyCell[k]] == myCell)
                    newOwner_[facesMyCell[k]] = myCell;
                else if(neighbour[facesMyCell[k]] == myCell)
                    newNeighbour_[facesMyCell[k]] = myCell;
                else
                {
                    FatalErrorInFunction
                    << "Face of merging cells is neither in owner nor in neighbour cell!"
                    << exit(FatalError);
                }
            }
            
            for(int s=0;s<mergeCells.size();s++)
            {
                labelList facesMergeCell = newCells[mergeCells[s]];
                for(int k=0;k<facesMergeCell.size();k++)
                {
                    if(mergeFaces.count(facesMergeCell[k]) != 0)
                        continue;
                
                    if(owner[facesMergeCell[k]] == mergeCells[s])
                        newOwner_[facesMergeCell[k]] = myCell;
                    else if(neighbour[facesMergeCell[k]] == mergeCells[s])
                        newNeighbour_[facesMergeCell[k]] = myCell;
                    else
                    {
                        FatalErrorInFunction
                        << "Face of merging cells is neither in owner nor in neighbour cell!"
                        << exit(FatalError);
                    }                    
                }
            }
        }
        else
        {
            if(oldCellNumToNewCellNum[i] == -1)
                oldCellNumToNewCellNum[i] = i;
        }
    }
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
    /*
    for(int i=0;i<newFaces_.size();i++)
    {
        Info<<"Face:"<<i<<" Owner:"<<newOwner_[i]<<" ";
        if(i < newNeighbour_.size())
            Info<<" Neighbor:"<<newNeighbour_[i]<<" ";
        for(int k=0;k<newFaces_[i].size();k++)
        {
            Info<<points[newFaces_[i][k]]<<"->";
        }
        if(newFaces_[i].size() > 0)
        {
            Info<<" with centre:"<<newFaces_[i].centre(points);
            Info<<" and normal vector:"<<newFaces_[i].normal(points);
            Info<<" and area:"<<newFaces_[i].mag(points)<<endl;
        }
        else
            Info<<endl;
    }
    */
    
    Info<<"Recompute faces, owner, neighbor ";
    t1 = std::chrono::high_resolution_clock::now();
    
    faceList newFaces__(newFaces_.size()-countDeleteFaces);
    labelList newOwner__(newOwner_.size()-countDeleteFaces);
    labelList newNeighbour__(newNeighbour_.size()-countDeleteFaces);
    //Info<<"Created Data Struc"<<endl;

    
    int countDel = 0;
    for(int i=0;i<newFaces_.size();i++)
        if(newOwner_[i] == -1)
            countDel++;
        
    if(countDel != countDeleteFaces)
    {
        FatalErrorInFunction
        << countDel<<"!="<<countDeleteFaces
        << exit(FatalError);
    }
    
    int insertCounter = 0;
    for(int i = 0;i<newFaces_.size();i++)
    {
        /*
        Info<<"Move "<<i<<" of "<<newFaces_.size()<<" to "<<insertCounter
        <<"/"<<newFaces_.size()-countDeleteFaces<<
        "-"<<newOwner_.size()-countDeleteFaces<<endl;
        */
        if(newOwner_[i] != -1)
        {
            //Info<<"0"<<endl;
            newFaces__[insertCounter] = newFaces_[i];
            //Info<<"1"<<endl;
            newOwner__[insertCounter] = newOwner_[i];
            //Info<<"2"<<endl;
            if(newOwner__[insertCounter] == -1)
            {            
                FatalErrorInFunction
                << "newOwner["<<insertCounter<<"]: "<<newOwner__[insertCounter]
                <<" from "<<newOwner_[i]
                << exit(FatalError);
            }
            //Info<<"3"<<endl;
            if(i < newNeighbour_.size())
                newNeighbour__[insertCounter] = newNeighbour_[i];
            insertCounter++;
        }
    }
    
    /*
    Info<<"End"<<endl;
    for(int i=0;i<newFaces__.size();i++)
    {
        Info<<"Face:"<<i<<" Owner:"<<newOwner__[i]<<" ";
        if(i < newNeighbour__.size())
            Info<<" Neighbor:"<<newNeighbour__[i]<<" ";
        for(int k=0;k<newFaces__[i].size();k++)
        {
            Info<<points[newFaces__[i][k]]<<"->";
        }
        if(newFaces__[i].size() > 0)
        {
            Info<<" with centre:"<<newFaces__[i].centre(points);
            Info<<" and normal vector:"<<newFaces__[i].normal(points);
            Info<<" and area:"<<newFaces__[i].mag(points)<<endl;
        }
        else
            Info<<endl;
        if(newOwner__[i] == -1)
        {            
            FatalErrorInFunction
            << "Stop."
            << exit(FatalError);
        }
    }
    */
    
    const polyBoundaryMesh& boundMesh = this->boundaryMesh();
    patchStarts = labelList(boundMesh.size());
    patchSizes = labelList(boundMesh.size());
    for(int i=0;i<boundMesh.size();i++)
    {
        patchStarts[i] = boundMesh[i].start()-countDeleteFaces;
        patchSizes[i] = boundMesh[i].faceCentres().size();
    }
    
    
    labelList cellReductionNumb(newCells.size());
    label count = 0;
    for(int i=0;i<cellReductionNumb.size();i++)
    {
        if(deletedCells[i] == 1)
        {
            count++;
            cellReductionNumb[i] = -1;
        }
        else
        {
            //Info<<"One none deleted"<<endl;
            cellReductionNumb[i] = count;
        }
    }
    label delNum = 0;
    for(int i=0;i<oldCellNumToNewCellNum.size();i++)
    {
        if(cellReductionNumb[i] != -1)
        {
            oldCellNumToNewCellNum[i] -= cellReductionNumb[i];
        }
        else
        {
            oldCellNumToNewCellNum[i] -= cellReductionNumb[oldCellNumToNewCellNum[i]];
            delNum++;
        }
    }
    Info<<"delNum:"<<delNum<<endl;
    label numDeletedCells = 0;
    for(int i=0;i<deletedCells.size();i++)
    {
        if(deletedCells[i] == 1)
            numDeletedCells++;
    }
    Info<<"numDeletedCells: "<<numDeletedCells<<endl;
    
    for(int i=0;i<newNeighbour__.size();i++)
    {
        if(cellReductionNumb[newOwner__[i]] != -1 &&
           cellReductionNumb[newNeighbour__[i]] != -1)
        {
            int temp = newOwner__[i];
            newOwner__[i] -= cellReductionNumb[newOwner__[i]];
            if(newOwner__[i] == -1)
            {
                FatalErrorInFunction
                << "Owner original "<<temp<<endl
                << "newOwner "<<newOwner__[i]<<endl
                << "because of reduction "<<cellReductionNumb[temp]
                << exit(FatalError);
            }
            newNeighbour__[i] -= cellReductionNumb[newNeighbour__[i]];
        }
        else
        {
            /*
            Info<<"newOwner__["<<i<<"]:"<<cellReductionNumb[newOwner__[i]]<<endl
            <<"newNeighbour__["<<i<<"]:"<<cellReductionNumb[newNeighbour__[i]]<<endl;
            */
            
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
    for(int i=newNeighbour__.size();i<newFaces__.size();i++)
    {
        if(cellReductionNumb[newOwner__[i]] != -1)
        {
            int temp = newOwner__[i];
            newOwner__[i] -= cellReductionNumb[newOwner__[i]];
            if(newOwner__[i] == -1)
            {
                FatalErrorInFunction
                << "Owner original "<<temp<<endl
                << "newOwner "<<newOwner__[i]<<endl
                << "because of reduction "<<cellReductionNumb[temp]
                << exit(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
    
    Info<<"cellReductionNumb.size() ="<<cellReductionNumb.size()<<endl;
    DynamicList<DynamicList<label>> newCellNumToOldCellNum;
    label numberDeletedCells = -1;
    for(int i=cellReductionNumb.size()-1;i>=0;i--)
    {
        if(cellReductionNumb[i] != -1)
        {
            Info<<"Non deleted:"<<i<<endl;
            numberDeletedCells = cellReductionNumb[i];
            Info<<"numberDeletedCells:"<<numberDeletedCells<<endl;
            break;
        }
    }
    if(numberDeletedCells == -1)
    {
        FatalErrorInFunction
        << "All cells deleted!"
        << exit(FatalError);  
    }
    newCellNumToOldCellNum.setSize(newCells.size()-numberDeletedCells);
    for(int i=0;i<oldCellNumToNewCellNum.size();i++)
    {
        if(oldCellNumToNewCellNum[i] >= newCellNumToOldCellNum.size())
        {
            Info<<"oldCellNumToNewCellNum["<<i<<"]:"<<oldCellNumToNewCellNum[i]<<" newCellNumToOldCellNum.size():"<<newCellNumToOldCellNum.size()<<endl;
            FatalErrorInFunction
            << "Wrong assignment!"
            << exit(FatalError); 
        }
        else
        {
            newCellNumToOldCellNum[oldCellNumToNewCellNum[i]].append(i);
        }
    }
    
    /*
    Info<<"newCell: "<<50491<<" was "<<newCellNumToOldCellNum[50941]<<endl;
    Info<<"mergeFaceOfCell:"<<mergeFaceOfCell[newCellNumToOldCellNum[50941][0]]
    <<" mergeNecessary: "<<mergeNecessary[newCellNumToOldCellNum[50941][0]]<<" with "<<
    "partialVolumeScale: "<<partialVolumeScale[newCellNumToOldCellNum[50941][0]]<<endl;
    */
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Test reset and test ";
    t1 = std::chrono::high_resolution_clock::now();
    
    testNewMeshData(newFaces__,newOwner__,newNeighbour__,patchStarts,patchSizes);
    
    resetPrimitives(Foam::clone(points),
                    Foam::clone(newFaces__),
                    Foam::clone(newOwner__),
                    Foam::clone(newNeighbour__),
                    patchSizes,
                    patchStarts,
                    true);
    
    /*
     * Uncomment later!!!
    testForCellSize
    (
        possibleMergeFaceArea,possibleMergeFaces, possibleMergeCells,
        oneMergeFaceSufficient, mergeNecessary, mergeFaceOfCell,partialThreeshold
    );
    */
    
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t" << time_span.count() << " seconds."<<endl;
}

labelList Foam::cutCellFvMesh::searchDown
(
    scalarListList& possibleMergeFaceArea,
    labelListList& possibleMergeFaces,
    labelListList& possibleMergeCells,
    boolList& oneMergeFaceSufficient,
    boolList& mergeNecessary,
    label count,
    std::unordered_set<label> cellReserved
)
{
    /*
    Info<<"possibleMergeFaceArea:"<<possibleMergeFaceArea.size()<<endl;
    Info<<"possibleMergeFaces:"<<possibleMergeFaces.size()<<endl;
    Info<<"possibleMergeCells:"<<possibleMergeCells.size()<<endl;
    Info<<"oneMergeFaceSufficient:"<<oneMergeFaceSufficient.size()<<endl;
    Info<<"mergeNecessary:"<<mergeNecessary.size()<<endl;
    
    Info<<"Starts"<<endl;
    Info<<"Height:"<<count<<"/"<<possibleMergeCells.size();
    */
    
    if(count < possibleMergeCells.size()-1)
    {
        //Info<<" in first";
        if(mergeNecessary[count] && cellReserved.find(count) == cellReserved.end())
        {
            cellReserved.insert(count);
            //Info<<" merge"<<endl;
            for(int i=0;i<possibleMergeCells[count].size();i++)
            {
                label oneCell = possibleMergeCells[count][i];
                label oneFace = possibleMergeFaces[count][i];
                std::unordered_set<label> cellReservedCpy;
                if(cellReserved.find(oneCell) == cellReserved.end())
                {
                    labelList retList;
                    cellReservedCpy = cellReserved;
                    cellReservedCpy.insert(oneCell);
                    retList = searchDown
                    (possibleMergeFaceArea,possibleMergeFaces,possibleMergeCells,
                     oneMergeFaceSufficient,mergeNecessary,count+1,cellReservedCpy);
                    if(retList.size() != 0)
                    {
                        labelList returnList = {oneFace};
                        returnList.append(retList);
                        //Info<<"Return first merge: "<<returnList.size()<<" from "<<count<<endl;
                        return returnList;
                    }
                }
            }
            labelList returnList(0);
            //Info<<"Return first: "<<returnList.size()<<" from "<<count<<endl;
            return returnList;
        }
        else
        {
            //Info<<" empty"<<endl;
            labelList returnList = {-1};
            //Info<<" 1"<<endl;
            labelList retList;
            retList = searchDown
            (possibleMergeFaceArea,possibleMergeFaces,possibleMergeCells,
             oneMergeFaceSufficient,mergeNecessary,count+1,cellReserved);
            //Info<<"Recursion"<<endl;
            if(retList.size() != 0)
            {
                returnList.append(retList);
                //Info<<"Return first empty: "<<returnList.size()<<" from "<<count<<endl;
                return returnList;
            }
            else
            {
                returnList = labelList(0);
                //Info<<"Return first empty: "<<returnList.size()<<" from "<<count<<endl;
                return returnList;
            }
        }
    }
    else
    {
        //Info<<" in second";
        if(mergeNecessary[count] && cellReserved.find(count) == cellReserved.end())
        {
            cellReserved.insert(count);
            //Info<<" merge"<<endl;
            for(int i=0;i<possibleMergeCells[count].size();i++)
            {
                label oneCell = possibleMergeCells[count][i];
                label oneFace = possibleMergeFaces[count][i];
                if(cellReserved.find(oneCell) == cellReserved.end())
                {
                    labelList returnList = {oneFace};
                    //Info<<"Return second merge: "<<returnList.size()<<" from "<<count<<endl;
                    return returnList;
                }
            }
            labelList returnList(0);
            //Info<<"Return second: "<<returnList.size()<<" from "<<count<<endl;
            return returnList;
        }
        else
        {
            labelList returnList = {-1};
            //Info<<"Return second empty: "<<returnList.size()<<" from "<<count<<endl;
            return returnList;
        }
    }
}

labelList Foam::cutCellFvMesh::searchDown_rec
(
    DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
    DynamicList<DynamicList<label>>& possibleMergeFaces,
    DynamicList<DynamicList<label>>& possibleMergeCells,
    DynamicList<bool>& oneMergeFaceSufficient,
    DynamicList<bool>& mergeNecessary,
    label count,
    DynamicList<DynamicList<label>>& blockedCells,
    std::unordered_set<label>& cellReserved
)
{
    //Info<<count;
    /*
    Info<<"possibleMergeFaceArea:"<<possibleMergeFaceArea.size()<<endl;
    Info<<"possibleMergeFaces:"<<possibleMergeFaces.size()<<endl;
    Info<<"possibleMergeCells:"<<possibleMergeCells.size()<<endl;
    Info<<"oneMergeFaceSufficient:"<<oneMergeFaceSufficient.size()<<endl;
    Info<<"mergeNecessary:"<<mergeNecessary.size()<<endl;
    
    Info<<"Starts"<<endl;
    Info<<"Height:"<<count<<"/"<<possibleMergeCells.size();
    */
    
    if(count < possibleMergeCells.size()-1)
    {
        //Info<<" mergeNecessary:"<<mergeNecessary[count];
        if(mergeNecessary[count] && cellReserved.find(count) == cellReserved.end())
        {
            cellReserved.insert(count);
            blockedCells[count].append(count);
            //Info<<" merge"<<endl;
            for(int i=0;i<possibleMergeCells[count].size();i++)
            {
                if(cellReserved.find(possibleMergeCells[count][i])
                    == cellReserved.end())
                {
                    //Info<<" -> "<<possibleMergeCells[count][i]<<endl;
                    cellReserved.insert(possibleMergeCells[count][i]);
                    blockedCells[count].append(possibleMergeCells[count][i]);

                    labelList retList = searchDown_rec
                    (possibleMergeFaceArea,possibleMergeFaces,possibleMergeCells,oneMergeFaceSufficient,mergeNecessary,count+1,blockedCells,cellReserved);
                    if(retList.size() != 0)
                    {
                        labelList returnList = {possibleMergeFaces[count][i]};
                        returnList.append(retList);
                        //Info<<"Return first merge: "<<returnList.size()<<" from "<<count<<endl;
                        return returnList;
                    }
                }
            }
            labelList returnList(0);
            for(int k=count; k<blockedCells.size();k++)
            {
                for(int l=0;l<blockedCells[k].size();l++)
                {
                    cellReserved.erase(blockedCells[k][l]);
                }
            }
            //Info<<"-> "<<-1<<endl;
            blockedCells.setSize(count);
            //Info<<"Return first: "<<returnList.size()<<" from "<<count<<endl;
            return returnList;
        }
        else
        {
            //Info<<"-> "<<-1<<endl;
            //Info<<" empty"<<endl;
            labelList returnList = {-1};
            //Info<<" 1"<<endl;
            labelList retList = searchDown_rec
            (possibleMergeFaceArea,possibleMergeFaces,possibleMergeCells,oneMergeFaceSufficient,mergeNecessary,count+1,blockedCells,cellReserved);
            //Info<<"Recursion"<<endl;
            if(retList.size() != 0)
            {
                returnList.append(retList);
                //Info<<"Return first empty: "<<returnList.size()<<" from "<<count<<endl;
                return returnList;
            }
            else
            {
                returnList = labelList(0);
                for(int k=count; k<blockedCells.size();k++)
                {
                    for(int l=0;l<blockedCells[k].size();l++)
                    {
                        cellReserved.erase(blockedCells[k][l]);
                    }
                }
                blockedCells.setSize(count);
                //Info<<"Return first empty: "<<returnList.size()<<" from "<<count<<endl;
                return returnList;
            }
        }
    }
    else
    {
        //Info<<" mergeNecessary:"<<mergeNecessary[count];
        if(mergeNecessary[count] && cellReserved.find(count) == cellReserved.end())
        {
            cellReserved.insert(count);
            blockedCells[count].append(count);

            //Info<<" merge"<<endl;
            for(int i=0;i<possibleMergeCells[count].size();i++)
            {
                if(cellReserved.find(possibleMergeCells[count][i]) == cellReserved.end())
                {
                    //Info<<"-> "<<possibleMergeCells[count][i]<<endl;
                    labelList returnList = {possibleMergeFaces[count][i]};
                    //Info<<"Return second merge: "<<returnList.size()<<" from "<<count<<endl;
                    return returnList;
                }
            }
            labelList returnList(0);
            for(int k=count; k<blockedCells.size();k++)
            {
                for(int l=0;l<blockedCells[k].size();l++)
                {
                    cellReserved.erase(blockedCells[k][l]);
                }
            }
            blockedCells.setSize(count);
            //Info<<"Return second: "<<returnList.size()<<" from "<<count<<endl;
            return returnList;
        }
        else
        {
            //Info<<"-> "<<-1<<endl;
            labelList returnList = {-1};
            //Info<<"Return second empty: "<<returnList.size()<<" from "<<count<<endl;
            return returnList;
        }
    }
}

bool mergeCellSelectionBlocks
(
    const label cellInd_rec,
    const label redCellInd_rec,
    const DynamicList<label>& mergeCells_rec,
    const label cellInd_block,
    const label redCellInd_block,
    const DynamicList<DynamicList<label>> mergeCells_block,
    std::unordered_map<label,label>& cellReserved
)
{
    for(int i=0;i<mergeCells_rec.size();i++)
    {
        label ins = mergeCells_rec[i];
        auto keyIt = cellReserved.find(ins);
        if(keyIt != cellReserved.end() && keyIt->second < cellInd_rec)
        {
            return false;
        }
    }
    
    bool blockCellIncluded = false;
    for(int i=0;i<mergeCells_rec.size();i++)
    {
        if(mergeCells_rec[i]==cellInd_block)
            blockCellIncluded = true;
    }
    if(blockCellIncluded)
        return true;
    
    bool oneMergePossAtBlock = true;
    std::unordered_set<label> blockedCells;
    blockedCells.insert(cellInd_rec);
    for(int i=0;i<mergeCells_rec.size();i++)
    {
        label ins = mergeCells_rec[i];
        blockedCells.insert(ins);
    }
//Info<<"backPoint:"<<cellInd_rec<<" trackBackArray:"<<mergeCells_rec<<endl;
    
    bool nonBlocking;
    for(int i=0;i<mergeCells_block.size();i++)
    {
//Info<<"blockPoint:"<<cellInd_block<<" blockArray:"<<mergeCells_block[i]<<endl;
        nonBlocking = true;
        for(int ii=0;ii<mergeCells_block[i].size();ii++)
        {
            if(blockedCells.count(mergeCells_block[i][ii]) == 1)
                nonBlocking = false;
        }
        if(nonBlocking)
        {
//Info<<"does not block"<<endl;
            oneMergePossAtBlock = true;
        }
    }
//Info<<"Block:"<<oneMergePossAtBlock<<endl;
    
    return oneMergePossAtBlock;
}

List<DynamicList<label>> Foam::cutCellFvMesh::searchDown_iter
(
    DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
    DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
    DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
    DynamicList<bool>& oneMergeFaceSufficient,
    DynamicList<bool>& mergeNecessary
)
{
    label mergeCounter = 0;
    for(int i=0;i<possibleMergeCells.size();i++)
    {
        if(mergeNecessary[i])
            mergeCounter++;
        
        if(possibleMergeFaces[i].size() != possibleMergeCells[i].size())
        {
            FatalErrorInFunction
            << " Invalid Data 1!"<<endl
            << exit(FatalError);
        }
        
        if((mergeNecessary[i] && possibleMergeFaces[i].size() == 0) || (!mergeNecessary[i] && possibleMergeFaces[i].size() != 0))
        {
            FatalErrorInFunction
            << " Invalid Data 2!"<<endl
            << exit(FatalError);
        }
    }
    Info<<"mergeCounter:"<<mergeCounter<<endl;
    DynamicList<DynamicList<DynamicList<label>>> possibleMergeFaces_red;
    possibleMergeFaces_red.setSize(mergeCounter);
    DynamicList<DynamicList<DynamicList<label>>> possibleMergeCells_red;
    possibleMergeCells_red.setSize(mergeCounter);
    DynamicList<bool> mergeNecessary_red;
    mergeNecessary_red.setSize(mergeCounter);
    labelList cellToRedInd(mergeNecessary.size());
    labelList redIndToCell(mergeCounter);
    label k = 0;
    
    for(int i=0;i<mergeNecessary.size();i++)
    {
        if(mergeNecessary[i])
        {
            possibleMergeFaces_red[k] = possibleMergeFaces[i];
            possibleMergeCells_red[k] = possibleMergeCells[i];
            mergeNecessary_red[k] = mergeNecessary[i];
            cellToRedInd[i] = k;
            redIndToCell[k] = i;
            k++;
        }
        else
        {
            cellToRedInd[i] = -1;
        }
    }
    
    for(int i=0;i<cellToRedInd.size();i++)
    {
        Info<<"i:"<<i<<endl;
        if(cellToRedInd[i]>=possibleMergeFaces_red.size())
        {
            Info<<"cellToRedInd["<<i<<"]:"<<cellToRedInd[i]<<" / "<<possibleMergeFaces_red.size()<<endl;
            Info<<"Bool:"<<(cellToRedInd[i]>=possibleMergeFaces_red.size())<<endl;
            FatalErrorInFunction
            << " Invalid Data 1!"<<endl
            << exit(FatalError);
        }
        if(cellToRedInd[i] != -1 && redIndToCell[cellToRedInd[i]] != i)
        {
            Info<<"cellToRedInd["<<i<<"]:"<<cellToRedInd[i]<<" / "<<possibleMergeFaces_red.size()<<endl;
            Info<<"redIndToCell["<<cellToRedInd[i]<<"]:"<<redIndToCell[cellToRedInd[i]]<<endl;
            FatalErrorInFunction
            << " Invalid Data 2!"<<endl
            << exit(FatalError);
        }
    }
    
    Info<<mergeCounter<<"/"<<possibleMergeCells.size()<<endl;
    
    Info<<endl;
    label count = 0;
    DynamicList<DynamicList<label>> blockedCells;
    blockedCells.setSize(possibleMergeCells_red.size());
    std::unordered_map<label,label> cellReserved;
    bool MergeFaceFound;
    DynamicList<label> mergeFace;
    DynamicList<label> mergeCell;
    DynamicList<label> temp;
    temp.append(-3);
    List<DynamicList<label>> assignList(possibleMergeCells.size(),temp);
    labelList tryedCells(possibleMergeCells_red.size(),0);
    
    label maxDepth = 0;
    label minDepth = 0;
    
    /*
    Info<<"possibleMergeCells["<<13016<<"].size():"<<possibleMergeCells[13016].size()<<endl;
    Info<<"mergeNecessary["<<13016<<"][0]:"<<mergeNecessary[13016]<<endl;
    Info<<"possibleMergeCells["<<13016<<"][0]:"<<possibleMergeCells[13016][0]<<endl;

    
    Info<<"possibleMergeCells["<<3753<<"].size():"<<possibleMergeCells[3753].size()<<endl;
    Info<<"mergeNecessary["<<3753<<"][0]:"<<mergeNecessary[3753]<<endl;
    */
    /*
    for(int i=0;i<possibleMergeCells.size();i++)
    {
        bool mayMerge = false;
        for(int k=0;k<possibleMergeCells[i].size();k++)
        {
            if(possibleMergeCells[i][k] == 3753)
                mayMerge = true;
        }
        if(mayMerge)
        {
            Info<<"-----------------------------------"<<endl;
            Info<<"possibleMergeCells["<<i<<"].size():"<<possibleMergeCells[i].size()<<endl;
            Info<<"mergeNecessary["<<i<<"]:"<<mergeNecessary[i]<<endl;
            Info<<"possibleMergeCells["<<i<<"][0]:"<<possibleMergeCells[i][0]<<endl;
        }
    }

    FatalErrorInFunction
    << " Temporary stop!"<<endl
    << exit(FatalError);
    */
    
    bool Change = false;
    
    Info<<"index:5404:"<<possibleMergeCells_red[5404]<<endl;
    
    DynamicList<label> problemCells;
    problemCells.setSize(possibleMergeCells_red.size());
    for(int i=0;i<possibleMergeCells_red.size();i++)
    {
        
    }
    
    Info<<"Start iterate count:"<<count<<"possibleMergeCells_red.size()"<<possibleMergeCells_red.size()<<endl;
    for(;count<possibleMergeCells_red.size();)
    {
        for( const auto& n : cellReserved ) 
        {
            if(cellToRedInd[n.second] >= possibleMergeCells_red.size())
            {
                Info<< "Key:[" << n.first << "] Value:[" << n.second << "]"<<endl;
                Info<<"possibleMergeCells_red.size():"<<possibleMergeCells_red.size()<<endl;
                Info<<"TrackBackPoint:"<<n.second<<endl;
                Info<<"TrackBackPointRed:"<<cellToRedInd[n.second]<<endl;
                Info<<"count:"<<count<<endl;
                FatalErrorInFunction
                << " Invalid Data 2!"<<endl
                << exit(FatalError);
            }
        }
        
        /*
        if(redIndToCell[count] < minDepth)
        {
            minDepth = redIndToCell[count];
            Change = true;
        }
        if(redIndToCell[count] > maxDepth)
        {
            maxDepth = redIndToCell[count];
            minDepth = redIndToCell[count];
            Change = true;
        }
        if(Change)
        {
            Info<<"minDepth:"<<minDepth<<" maxDepth:"<<maxDepth<<endl;
            Change = false;
        }
        */
        
        
        if(count < minDepth)
        {
            minDepth = count;
            Change = true;
        }
        if(count > maxDepth)
        {
            maxDepth = count;
            minDepth = count;
            Change = true;
        }
        if(Change)
        {
            Info<<"minDepth:"<<minDepth<<" maxDepth:"<<maxDepth<<endl;
            /*
            Info<<" index:"<<maxDepth<<":"<<possibleMergeCells_red[maxDepth]<<endl;
            for(int a=0;a<possibleMergeCells_red[maxDepth].size();a++)
            {
                for(int b=0;b<possibleMergeCells_red[maxDepth][a].size();b++)
                {
                    auto keyIt = cellReserved.find(possibleMergeCells_red[count][a][b]);
                    if(keyIt != cellReserved.end())
                    {
                        Info<<"Cell: "<<possibleMergeCells_red[count][a][b]<<
                        " blocked by: "<<cellToRedInd[keyIt->second]<<endl;
                    }
                }
            }
            */
            Change = false;
        }
    


//Info<<"------"<<count<<"------"<<redIndToCell[count]<<"-------"<<"tryedCells["<<count<<"] = "<<tryedCells[count]<<"/"<<"possibleMergeCells_red["<<count<<"] = "<<possibleMergeCells_red[count].size()<<endl;

        if(mergeNecessary_red[count])
        /* Decision A: Enters if block if merge is necessary and the cell is not already used for
         * a merge with another cell
         */
        {
            if(cellReserved.count(redIndToCell[count]) == 0)
            /* Decision B: Enters block if the cell is not already used for
            * a merge with another cell
            */
            {
                DynamicList<label> temp;
                temp.append(-4);
                MergeFaceFound = false;
                mergeFace = temp;
                mergeCell = temp;
            
//Info<<"tryedCells["<<count<<"] = "<<tryedCells[count]<<"/"<<"possibleMergeCells_red["<<count<<"] = "<<possibleMergeCells_red[count].size()<<endl;

                for(int i=tryedCells[count];i<possibleMergeCells_red[count].size();i++,tryedCells[count]++)
                {
//Info<<"Merge Cell:"<<possibleMergeCells_red[count][i];

                    bool cellsNotBlocked = true;
                    DynamicList<label> blockedBecausOf;
                    for(int s=0;s<possibleMergeCells_red[count][i].size();s++)
                    {
                        if(cellReserved.count(possibleMergeCells_red[count][i][s]) != 0)
                        {
                            auto keyIt = cellReserved.find(possibleMergeCells_red[count][i][s]);
                            if(keyIt == cellReserved.end())
                            {
                                FatalErrorInFunction
                                << " Error1 !"<<endl
                                << exit(FatalError); 
                            }
                            blockedBecausOf.append(keyIt->second);
                            if(keyIt->second>=mergeNecessary.size())
                            {
                                Info<<"cell: "<<keyIt->first<<endl;
                                FatalErrorInFunction
                                << " Error!"<<endl
                                << exit(FatalError);
                            }
                            cellsNotBlocked = false;
                        }
                    }
                    if(cellsNotBlocked)
                    {
//Info<<" not blocked"<<endl;
                        MergeFaceFound = true;

//Info<<"Found face for cell: "<<count<<"  ";
//Info<<"tryedCells["<<count<<"] = "<<tryedCells[count]<<"/"<<"possibleMergeCells["<<count<<"] = "<<possibleMergeCells[count].size()<<endl;

                        mergeFace = possibleMergeFaces_red[count][i];
                        mergeCell = possibleMergeCells_red[count][i];

                        tryedCells[count]++;
                        break;
                    }
                    else
                    {
//Info<<" blocked because ";
                        for(int x=0;x<blockedBecausOf.size();x++)
                        {
//Info<<"("<<blockedBecausOf[x]<<"|"<<cellToRedInd[blockedBecausOf[x]]<<") ";
                        }
//Info<<endl;
                    }
                }
                if(MergeFaceFound == false)
                /* Decision B: The iteration across all possibleMergeCells failed. This happens if
                * all the possible merge cells between tryed[count] und possibleMergeCells[count].size
                * are already blocked.
                * The result is a backtracking to the first previous cell where another selection was still
                * possible
                */
                {
//Info<<"Not found"<<endl;
/*
Info<<"Merge face not found for "<<count<<endl;
Info<<"tryedCells["<<count<<"] = "<<tryedCells[count]<<"/"<<"possibleMergeCells["<<count<<"] = "<<possibleMergeCells[count].size()<<" assignList["<<count<<"]:"<<assignList[count]<<"  mergeNecessary["<<count<<"]:"<<mergeNecessary[count]<<endl;
Info<<"tryedCells["<<count-1<<"] = "<<tryedCells[count-1]<<"/"<<"possibleMergeCells["<<count-1<<"] = "<<possibleMergeCells[count-1].size()<<" assignList["<<count-1<<"]:"<<assignList[count-1]<<"  mergeNecessary["<<count-1<<"]:"<<mergeNecessary[count-1]<<endl;
*/
                    DynamicList<label> trackBackPoints;
                    DynamicList<label> trackBackPointsMergeInd;
//Info<<"Size:"<<possibleMergeCells_red[count].size()<<endl;
//Info<<"Size2:"<<trackBackPoints.size()<<endl;
                    for(int s=0;s<possibleMergeCells_red[count].size();s++)
                    {
//Info<<"s1:"<<s<<endl;
                        if(possibleMergeCells_red[count][s].size() == 0)
                        {
                            FatalErrorInFunction
                            << " Zero length merge section!"<<endl
                            << exit(FatalError);
                        }
                        bool blockedBecausThis = false;
                        for(int z=0;z<possibleMergeCells_red[count][s].size();z++)
                        {
                            auto keyIt = cellReserved.find(possibleMergeCells_red[count][s][z]);
                            if(keyIt != cellReserved.end())
                            {
if(count == 90){
//Info<<"-"<<s<<"-"<<"MergeCell: "<<possibleMergeCells_red[count][s][z]<<" blocked by "<<keyIt->second<<" from merge with "<<redIndToCell[count]<<endl;
}
                                blockedBecausThis = true;
                            }
                            else
                            {
if(count == 90){
//Info<<"-"<<s<<"-"<<"MergeCell: "<<possibleMergeCells_red[count][s][z]<<" not blocked from merge with "<<redIndToCell[count]<<endl;
}
                            }
                        }
//Info<<"s2:"<<s<<endl;
                        if(!blockedBecausThis)
                            continue;
//Info<<"s3:"<<s<<endl;
                        DynamicList<label> trackPointsForMergeOption;
                        for(int z=0;z<possibleMergeCells_red[count][s].size();z++)
                        {
//Info<<z<<"/"<<possibleMergeCells_red[count][s].size()<<":"<<possibleMergeCells_red[count][s][z]<<endl;
                            auto keyIt = cellReserved.find(possibleMergeCells_red[count][s][z]);
                            if(keyIt != cellReserved.end())
                            {
                                trackPointsForMergeOption.append(keyIt->second);
                            }
                            else
                                continue;
                            if(cellToRedInd[keyIt->second] >= count || cellToRedInd[keyIt->second] < 0)
                            {
                                Info<<"possibleMergeCells_red.size():"<<possibleMergeCells_red.size()<<endl;
                                Info<<"TrackBackPoint:"<<keyIt->second<<endl;
                                Info<<"TrackBackPointRed:"<<cellToRedInd[keyIt->second]<<endl;
                                Info<<"count:"<<count<<endl;
                                FatalErrorInFunction
                                << " More general possible Track Back Point is wrong!"<<endl
                                << exit(FatalError);
                            }           
                        }
                        
//Info<<"s4:"<<s<<endl;
                        /*
                        if(keyIt == cellReserved.end())
                        {
                            Info<<"s:"<<s<<endl;
                            Info<<"tryedCells["<<count<<"] = "<<tryedCells[count]<<"/"<<"possibleMergeCells_red["<<count<<"] = "<<possibleMergeCells_red[count].size()<<endl;
                            FatalErrorInFunction
                            << " Merge possiblity not taken but cell is not blocked!"<<endl
                            << exit(FatalError);
                        }
                        */
                        for(int z=1;z<possibleMergeCells_red[count][s].size();z++)
                        {
                            auto keyIt = cellReserved.find(possibleMergeCells_red[count][s][z]);
                            if(keyIt != cellReserved.end() && keyIt->first != possibleMergeCells_red[count][s][z])
                            {
                                FatalErrorInFunction
                                << " Something is wrong here!"<<endl
                                << exit(FatalError);
                            }
                        }
//Info<<"s6:"<<s<<endl;
                        if(trackPointsForMergeOption.size()==0)
                        {
                            FatalErrorInFunction
                            << " Something is wrong here!"<<endl
                            << exit(FatalError);
                        }
//Info<<"s7:"<<s<<endl;
                        label min = trackPointsForMergeOption[0];
                        for(int z=0;z<trackPointsForMergeOption.size();z++)
                        {
                            if(min > trackPointsForMergeOption[z])
                                min = trackPointsForMergeOption[z];
                        }
//Info<<"s8:"<<s<<endl;
                        if(cellToRedInd[min] >= count || cellToRedInd[min] < 0)
                        {
                            Info<<"possibleMergeCells_red.size():"<<possibleMergeCells_red.size()<<endl;
                            Info<<"TrackBackPoint:"<<min<<endl;
                            Info<<"TrackBackPointRed:"<<cellToRedInd[min]<<endl;
                            Info<<"count:"<<count<<endl;
                            FatalErrorInFunction
                            << " Possible Track Back Point is wrong!"<<endl
                            << exit(FatalError);
                        }
//Info<<"Append"<<endl;
                        trackBackPoints.append(min);
                        trackBackPointsMergeInd.append(s);
                    }                    
/*                    
Info<<"trackBackPoints:";
for(int s=0;s<trackBackPoints.size();s++) Info<<"  "<<trackBackPoints[s]<<"/"<<cellToRedInd[trackBackPoints[s]];
Info<<endl;
*/
                    if(trackBackPoints.size() < 1)
                    {
                        FatalErrorInFunction
                        << " Something wrong!"<<endl
                        << exit(FatalError);
                    }

                    int js;
                    label keyMergeInd,keyBackPoint;
                    for(int s=1;s<trackBackPoints.size();s++)
                    {
                        keyBackPoint = trackBackPoints[s];
                        keyMergeInd = trackBackPointsMergeInd[s];
                        js = s-1;
                        while(js>=0 && trackBackPoints[js] < keyBackPoint)
                        {
                            trackBackPoints[js+1] = trackBackPoints[js];
                            trackBackPointsMergeInd[js+1] = trackBackPointsMergeInd[js];
                            js--;
                        }
                        trackBackPoints[js+1] = keyBackPoint;
                        trackBackPointsMergeInd[js+1] = keyMergeInd;
                    }
/*
for(int ijk=0;ijk<=count;ijk++)
{
    Info<<"redInd:"<<ijk<<" cellInd:"<<redIndToCell[ijk]<<" "<<tryedCells[ijk]<<"/"<<possibleMergeCells_red[ijk].size()<<endl;
}
*/                    
//Info<<"Sort: tBP:"<<trackBackPoints<<" tBPMI:"<<trackBackPointsMergeInd<<endl;

                    label bestTrackBackPoint = trackBackPoints[trackBackPoints.size()-1];
/*
if(count == 90){
for(int s=0;s<trackBackPoints.size();s++)
{
Info<<"----------"<<s<<"----------------------"<<endl;
Info<<"Red Track back: "<<cellToRedInd[trackBackPoints[s]]<<" Merge:"<<possibleMergeCells_red[cellToRedInd[trackBackPoints[s]]]<<" tryed:"<<tryedCells[cellToRedInd[trackBackPoints[s]]]<<endl<<"Cell:"<<trackBackPoints[s]<<endl;
Info<<"----------"<<s<<"----------------------"<<endl;
}
Info<<"--------------------------------"<<endl;
Info<<"Red Track back: "<<count<<" Merge:"<<possibleMergeCells_red[count]<<"Cell:"<<redIndToCell[count]<<endl;
Info<<"--------------------------------"<<endl;
}
*/
                    bool oneTrackBackPointUnblocks = false;
                    for(int s=0;s<trackBackPoints.size();s++)
                    {
                        bool trackBackCanFreeCount = false;
                        label redIntBackPoint = cellToRedInd[trackBackPoints[s]];
                        label cellBackPoint = trackBackPoints[s];
                        label tryContPoint = tryedCells[redIntBackPoint];
//Info<<"redInt"<<redIntBackPoint<<" tryContPoint:"<<tryContPoint<<" possibleMergeCells_red["<<redIntBackPoint<<"].size()="<<possibleMergeCells_red[redIntBackPoint].size()<<endl;
                        for(int ss=tryContPoint;ss<possibleMergeCells_red[redIntBackPoint].size();ss++)
                        {
                            bool oneMergePossAtBlock = mergeCellSelectionBlocks
                                (cellBackPoint,redIntBackPoint,possibleMergeCells_red[redIntBackPoint][ss],                                                                   redIndToCell[count],count,possibleMergeCells_red[count],cellReserved);
                                
//Info<<ss<<" Blocks: "<<trackBackPoints[s]<<" Bool:"<<oneMergePossAtBlock<<endl;
                            
                            if(oneMergePossAtBlock)
                            {
                                trackBackCanFreeCount = true;
                                tryedCells[redIntBackPoint] = ss;
                                oneTrackBackPointUnblocks = true;
//Info<<"Non block: "<<trackBackPoints[s]<<endl;
                                break;
                            }
                        }
                        if(trackBackCanFreeCount)
                        {
                            bestTrackBackPoint = trackBackPoints[s];
                            break;
                        }
                    }
//Info<<"btbp: "<<bestTrackBackPoint<<endl;
//if(count == 90){ Info<<"bestTrackBackPoint:"<<bestTrackBackPoint<<endl;}
                    if(!oneTrackBackPointUnblocks && false)
                    {
                        for(;; bestTrackBackPoint--)
                        {
                            if(bestTrackBackPoint == -1)
                            {
                                Info<<"Backtracking from cell: "<<bestTrackBackPoint<<endl;
                                FatalErrorInFunction
                                << "Failed in Merging Selection! There is no found combination for all merging cells."<<endl
                                << exit(FatalError);
                            }
                            label redIndCntBck = cellToRedInd[bestTrackBackPoint];
                            if((redIndCntBck==-1 && assignList[bestTrackBackPoint][0]!=-3) || (assignList[bestTrackBackPoint][0]==-3 && redIndCntBck!=-1))
                            {
                                FatalErrorInFunction
                                << " Can not happen!"<<endl
                                << exit(FatalError);
                            }
                            if(redIndCntBck==-1)
                                continue;
                            
                            bool trackBackCanFreeCount = false;
                            label tryContPoint = tryedCells[redIndCntBck];
                            for(int ss=tryContPoint;ss<possibleMergeCells_red[redIndCntBck].size();ss++)
                            {
                                bool oneMergePossAtBlock = mergeCellSelectionBlocks
                                (bestTrackBackPoint,redIndCntBck,possibleMergeCells_red[redIndCntBck][ss],                                                                   redIndToCell[count],count,possibleMergeCells_red[count],cellReserved);
                                                            
                                if(oneMergePossAtBlock)
                                {
                                    trackBackCanFreeCount = true;
                                    tryedCells[redIndCntBck] = ss;
                                    break;
                                }
                            }
                            if(trackBackCanFreeCount)
                            {
                                break;
                            }
                        }                        
                    }
                        
                    bestTrackBackPoint = cellToRedInd[bestTrackBackPoint];
/*
if(count == 90){ 
    Info<<"bestTrackBackPointRed:"<<bestTrackBackPoint<<endl;
    FatalErrorInFunction<<"Temporary stop!"<<endl<<exit(FatalError);
}
*/

//Info<<"bestTrackBackPoint:"<<bestTrackBackPoint<<endl;
                    if(bestTrackBackPoint >= count || bestTrackBackPoint < 0)
                    {
                        Info<<"bestTrackBackPoint:"<<bestTrackBackPoint<<endl;
                        Info<<"count:"<<count<<endl;
                        FatalErrorInFunction
                        << " Track Back Point is wrong!"<<endl
                        << exit(FatalError);
                    }
                    //////////////////////////////////////////----------------------------

                    int backtrackingIndex = -1;                    
                    for(int cntBck = bestTrackBackPoint; cntBck>=0; cntBck--)
                    {
                        if(assignList[redIndToCell[cntBck]].size() == 1 &&
                           assignList[redIndToCell[cntBck]][0] == -3)
                        {
                            FatalErrorInFunction
                            << " Backtracking to untreated cell!"<<endl
                            << exit(FatalError);
                        }
                        if((assignList[redIndToCell[cntBck]].size() == 1 && 
                            assignList[redIndToCell[cntBck]][0] == -1)   ||
                            !mergeNecessary_red[cntBck])
                        {
                            FatalErrorInFunction
                            << " Backtracking to nonmerge Cell or -1 assign cell!"<<endl
                            << exit(FatalError);
                        }
                        if(assignList[redIndToCell[cntBck]].size() == 0)
                        {
                            FatalErrorInFunction
                            << " Backtracking empty assign list!"<<endl
                            << exit(FatalError);
                        }
                        
                        if((assignList[redIndToCell[cntBck]].size() > 0 && 
                            assignList[redIndToCell[cntBck]][0] != -2) && tryedCells[cntBck]<possibleMergeCells_red[cntBck].size())
                        {
                            backtrackingIndex = cntBck;
                            break;
                        }
                    }
                    if(backtrackingIndex != -1)
                    /* Decision C: If the backtracking resulted in a cell the backtracking was succesful.
                    * The iteration continues with this specific cell. 
                    */
                    {
                        if(!mergeNecessary_red[count])
                        {
                            FatalErrorInFunction
                            << " Backtracking to none merge cell!"<<endl
                            << exit(FatalError);
                        }
                        if(assignList[redIndToCell[backtrackingIndex]].size()==0)
                        {
                            FatalErrorInFunction
                            << " Backtracking to zero length assign List can not happen!"<<endl
                            << exit(FatalError);
                        }
                        if(assignList[redIndToCell[backtrackingIndex]][0] < 0)
                        {
                            Info<<endl<<"assignList["<<backtrackingIndex<<"][0]="<<assignList[backtrackingIndex][0]<<endl;
                            FatalErrorInFunction
                            << " Backtracking to cell with assign list < 0!"<<endl
                            << exit(FatalError);
                        }
                        if(tryedCells[backtrackingIndex] <= 0)
                        {
                            Info<<"tryedCells["<<backtrackingIndex<<"]="<<tryedCells[backtrackingIndex]<<endl;
                            FatalErrorInFunction
                            << " Backtracking to cell tryedCells <= 0!"<<endl
                            << exit(FatalError);
                        }
                        
                        for(int k=backtrackingIndex; k<=count;k++)
                        {
/*
Info<<"<<<<<";
*/
                            for(int l=0;l<blockedCells[k].size();l++)
                            {
                                cellReserved.erase(blockedCells[k][l]);
/*
Info<<">>>>>>>>>>>Cleared["<<k<<"]:"<<blockedCells[k][l]<<endl;
*/
                            }
                        }
                        //Clean list of blockedCells
                        for(int k=backtrackingIndex+1;k<=count;k++)
                        {
                            tryedCells[k] = 0;
                            blockedCells[k].setSize(0);
/*
Info<<"Clear tryedCells and blockedCells ["<<k<<"]"<<endl;
*/
                                    
                        }
                        for(int o=backtrackingIndex+1;o<=count;o++)
                        {
                            if(tryedCells[o] != 0)
                            {
                                Info<<endl<<endl<<"Backtracking from "<<redIndToCell[count]<<" to "<<redIndToCell[backtrackingIndex]<<endl;
                                Info<<"tryedCells["<<o<<"]="<<tryedCells[o]<<endl;
                                FatalErrorInFunction
                                << " Backtracking to cell while tryedCells is unequal 0!"<<endl
                                << exit(FatalError);
                            }
                        }
                        
                        blockedCells[backtrackingIndex].setSize(0);
                        
                        count = backtrackingIndex;
/*
Info<<"Go to "<<count<<endl;
*/
                    }
                    else
                    /* Decision C: The backtracking did not found a backtracking cell. The result is an abort.
                    */
                    {
                        Info<<"Backtracking from cell: "<<redIndToCell[count]<<endl;
                        FatalErrorInFunction
                        << " Failed in Merging Selection! There is no found combination for all merging cells."<<endl
                        << exit(FatalError);
                    }
                }
                else
                /* Decision B: A possible merge cell was found. As a result the two cells to merge are
                * inserted both in the cellReserved map as well as in the blockedCells list             * 
                */
                {
                    if(redIndToCell[count]>=mergeNecessary.size() || count>=possibleMergeCells_red.size())
                    {
                        Info<<"Inserted Backtracking number from count: "<<count<<endl;
                        Info<<"redIndToCell[count]: "<<redIndToCell[count]<<endl;
                        Info<<"Inversion: "<<cellToRedInd[redIndToCell[count]]<<endl;
                        FatalErrorInFunction
                        << " Temporary stop."<<endl
                        << exit(FatalError);
                    }
                    std::pair<label,label> ins1(redIndToCell[count],redIndToCell[count]);
                    cellReserved.insert(ins1);
                    blockedCells[count].append(redIndToCell[count]);

                    for(int u=0;u<mergeCell.size();u++)
                    {
                        if(redIndToCell[count]>=mergeNecessary.size() || count>=possibleMergeCells_red.size())
                        {
                            Info<<"Inserted Backtracking number from count: "<<count<<endl;
                            Info<<"redIndToCell[count]: "<<redIndToCell[count]<<endl;
                            Info<<"Inversion: "<<cellToRedInd[redIndToCell[count]]<<endl;
                            FatalErrorInFunction
                            << " Temporary stop."<<endl
                            << exit(FatalError);
                        }
                        std::pair<label,label> ins2(mergeCell[u],redIndToCell[count]);
                        cellReserved.insert(ins2);
                        blockedCells[count].append(mergeCell[u]);
                    }
/*
Info<<">>>> Added "<<count<<" and "<<mergeCell<<" to blockedCells["<<count<<"]"<<endl;
*/
                    assignList[redIndToCell[count]] = mergeFace;
                    count++;
                }
            }
            else
            /* Decision B: Enters else block for cells that are already used for merge.
            */
            {
                DynamicList<label> temp;
                temp.append(-2);
                assignList[redIndToCell[count]] = temp;
                tryedCells[count] = 0;
                count++;
            }
        }
        else
        /* Decision A: Enters else block for cells that are not too small. The assignList is
         * filled with -1 for these cells.
         */
        {
            DynamicList<label> temp;
            temp.append(-1);
            assignList[redIndToCell[count]] = temp;
            tryedCells[count] = 0;
            count++;
        }
    }
    
    Info<<"Fin"<<endl;
    
    for(int i=0;i<assignList.size();i++)
    {
        if(mergeNecessary[i])
        {
            if(assignList[i].size()==0)
            {
                FatalErrorInFunction
                << " Error empty assign List at i:"<<i<<endl
                << exit(FatalError);
            }
            if(assignList[i][0] == -3)
            {
                FatalErrorInFunction
                << " Too small cell was not treated by backtracking algorithm!"<<endl
                << exit(FatalError);
            }
            
            if(assignList[i][0] == -1)
            {
                FatalErrorInFunction
                << " Too small cell was listed as large enough!"<<endl
                << exit(FatalError);
            }
        }
        else
        {
            if(assignList[i].size()==0)
            {
                FatalErrorInFunction
                << " Error empty assign List at i:"<<i<<endl
                << exit(FatalError);
            }
            if(assignList[i].size()!=1)
            {
                FatalErrorInFunction
                << " Error overfull assign List at i:"<<i<<endl
                << exit(FatalError);
            }
            if(assignList[i][0] != -3)
            {
                FatalErrorInFunction
                << " Backtracking algorithm wrote inside large enough cell. Something is wrong here!"<<endl
                << exit(FatalError);
            }
            assignList[i][0] = -1;
        }
    }
    
    /*
    labelList assList(assignList.size());
    for(int i=0;i<assignList.size();i++)
    {
        if(assignList[i].size() != 1)
        {
            FatalErrorInFunction
            << " Error!"<<endl
            << exit(FatalError);
        }
        assList[i] = assignList[i][0];
    }
    return assList;
    */
    return assignList;
    /*
    Label index
    -1 : Cell is not too small
    -2 : cell is already merged by other cell
    -3 : No assignment
    */
}

List<DynamicList<label>> Foam::cutCellFvMesh::searchDown_iter_preBlock
(
    const labelList& owner,
    const labelList& neighbour,
    DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
    DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
    DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
    DynamicList<bool>& oneMergeFaceSufficient,
    DynamicList<bool>& mergeNecessary
)
{
    label mergeCounter = 0;
    for(int i=0;i<possibleMergeCells.size();i++)
    {
        if(mergeNecessary[i])
            mergeCounter++;
        
        if(possibleMergeFaces[i].size() != possibleMergeCells[i].size())
        {
            FatalErrorInFunction
            << " Invalid Data 1!"<<endl
            << exit(FatalError);
        }
        
        if((mergeNecessary[i] && possibleMergeFaces[i].size() == 0) || (!mergeNecessary[i] && possibleMergeFaces[i].size() != 0))
        {
            FatalErrorInFunction
            << " Invalid Data 2!"<<endl
            << exit(FatalError);
        }
    }
    Info<<"mergeCounter:"<<mergeCounter<<endl;
    DynamicList<DynamicList<DynamicList<label>>> possibleMergeFaces_red;
    possibleMergeFaces_red.setSize(mergeCounter);
    DynamicList<DynamicList<DynamicList<label>>> possibleMergeCells_red;
    possibleMergeCells_red.setSize(mergeCounter);
    DynamicList<bool> mergeNecessary_red;
    mergeNecessary_red.setSize(mergeCounter);
    labelList cellToRedInd(mergeNecessary.size());
    labelList redIndToCell(mergeCounter);
    label k = 0;
    
    for(int i=0;i<mergeNecessary.size();i++)
    {
        if(mergeNecessary[i])
        {
            possibleMergeFaces_red[k] = possibleMergeFaces[i];
            possibleMergeCells_red[k] = possibleMergeCells[i];
            mergeNecessary_red[k] = mergeNecessary[i];
            cellToRedInd[i] = k;
            redIndToCell[k] = i;
            k++;
        }
        else
        {
            cellToRedInd[i] = -1;
        }
    }
    
    for(int i=0;i<cellToRedInd.size();i++)
    {
        //Info<<"i:"<<i<<endl;
        if(cellToRedInd[i]>=possibleMergeFaces_red.size())
        {
            Info<<"cellToRedInd["<<i<<"]:"<<cellToRedInd[i]<<" / "<<possibleMergeFaces_red.size()<<endl;
            Info<<"Bool:"<<(cellToRedInd[i]>=possibleMergeFaces_red.size())<<endl;
            FatalErrorInFunction
            << " Invalid Data 1!"<<endl
            << exit(FatalError);
        }
        if(cellToRedInd[i] != -1 && redIndToCell[cellToRedInd[i]] != i)
        {
            Info<<"cellToRedInd["<<i<<"]:"<<cellToRedInd[i]<<" / "<<possibleMergeFaces_red.size()<<endl;
            Info<<"redIndToCell["<<cellToRedInd[i]<<"]:"<<redIndToCell[cellToRedInd[i]]<<endl;
            FatalErrorInFunction
            << " Invalid Data 2!"<<endl
            << exit(FatalError);
        }
    }
    
    Info<<mergeCounter<<"/"<<possibleMergeCells.size()<<endl;
    
    Info<<endl;
    label count = 0;
    DynamicList<DynamicList<label>> blockedCells;
    blockedCells.setSize(possibleMergeCells_red.size());
    std::unordered_map<label,label> cellReserved;
    bool MergeFaceFound;
    DynamicList<label> mergeFace;
    DynamicList<label> mergeCell;
    DynamicList<label> temp;
    temp.append(-3);
    List<DynamicList<label>> assignList(possibleMergeCells.size(),temp);
    labelList tryedCells(possibleMergeCells_red.size(),0);
    
    std::unordered_multimap<label,std::pair<label,label>> cellPreBlock;
    DynamicList<DynamicList<bool>> cellMergPosBlocked_red;
    cellMergPosBlocked_red.setSize(possibleMergeCells_red.size());
    
    DynamicList<DynamicList<label>> cellMergPosBlocked_red_Reason;
    cellMergPosBlocked_red_Reason.setSize(possibleMergeCells_red.size());
    
    DynamicList<DynamicList<label>> cellMergPosBlockedMulti_red;
    cellMergPosBlockedMulti_red.setSize(possibleMergeCells_red.size());
    
    List<bool> cellMergDone_red(possibleMergeCells_red.size(),false);
    List<label> cellMergDoneMult_red(possibleMergeCells_red.size(),0);
    List<DynamicList<std::pair<label,label>>> preBlockedOptions_red(possibleMergeCells_red.size());
    List<DynamicList<label>> preCellMergDone_red(possibleMergeCells_red.size());
    for(int a=0;a<possibleMergeCells_red.size();a++)
    {
        cellMergPosBlocked_red[a].setSize(possibleMergeCells_red[a].size(),false);
        cellMergPosBlocked_red_Reason[a].setSize(possibleMergeCells_red[a].size(),possibleMergeCells.size());
        cellMergPosBlockedMulti_red[a].setSize(possibleMergeCells_red[a].size(),0);
        for(int b=0;b<possibleMergeCells_red[a].size();b++)
        {
            for(int c=0;c<possibleMergeCells_red[a][b].size();c++)
            {
                std::pair<label,label> index(a,b);
                std::pair<label,std::pair<label,label>> data(possibleMergeCells_red[a][b][c],index);                
                cellPreBlock.insert(data);
            }
        }
    }
    
    label maxDepth = 0;
    label minDepth = 0;
    
    bool Change = false;
        
    DynamicList<label> problemCells;
    problemCells.setSize(possibleMergeCells_red.size());
    
    Info<<"Start iterate count:"<<count<<" possibleMergeCells_red.size() "<<possibleMergeCells_red.size()<<endl;
    for(;count<possibleMergeCells_red.size();)
    {   
        for( const auto n : cellReserved ) 
        {
            if(cellToRedInd[n.second] >= possibleMergeCells_red.size())
            {
                Info<< "Key:[" << n.first << "] Value:[" << n.second << "]"<<endl;
                Info<<"possibleMergeCells_red.size():"<<possibleMergeCells_red.size()<<endl;
                Info<<"TrackBackPoint:"<<n.second<<endl;
                Info<<"TrackBackPointRed:"<<cellToRedInd[n.second]<<endl;
                Info<<"count:"<<count<<endl;
                FatalErrorInFunction
                << " Invalid Data 2!"<<endl
                << exit(FatalError);
            }
        }
        
        if(count < minDepth)
        {
            minDepth = count;
            Change = true;
        }
        if(count > maxDepth)
        {
            maxDepth = count;
            minDepth = count;
            Change = true;
        }
        if(Change)
        {
            Info<<"minDepth:"<<minDepth<<"/"<<redIndToCell[minDepth]
                <<" maxDepth:"<<maxDepth<<"/"<<redIndToCell[maxDepth]<<endl;
            Change = false;
        }

//Info<<"------"<<count<<"------"<<redIndToCell[count]<<"-------"<<"tryedCells["<<count<<"] = "<<tryedCells[count]<<"/"<<"possibleMergeCells_red["<<count<<"] = "<<possibleMergeCells_red[count].size()<<endl;

        if(mergeNecessary_red[count])
        /* Decision A: Enters if block if merge is necessary and the cell is not already used for
         * a merge with another cell
         */
        {
            if(cellReserved.count(redIndToCell[count]) == 0 && !cellMergDone_red[count])
            /* Decision B: Enters block if the cell is not already used for
            * a merge with another cell
            */
            {
                DynamicList<label> temp;
                temp.append(-4);
                MergeFaceFound = false;
                mergeFace = temp;
                mergeCell = temp;
                label tryedCellsStartPoint = tryedCells[count];
                bool allCellsAreBlocked = true;
                bool allCellsWillBlock  = true;
                DynamicList<bool> cellsAreBlocked;
                cellsAreBlocked.setSize(possibleMergeCells_red[count].size()-tryedCellsStartPoint);
                DynamicList<bool> cellsWillBlock;
                cellsWillBlock.setSize(possibleMergeCells_red[count].size()-tryedCellsStartPoint);
                DynamicList<DynamicList<DynamicList<label>>> cellsThatWillBeBlocked;
                cellsThatWillBeBlocked.setSize(possibleMergeCells_red[count].size()-tryedCellsStartPoint);
                
//Info<<"tryedCells["<<count<<"] = "<<tryedCells[count]<<"/"<<"possibleMergeCells_red["<<count<<"] = "<<possibleMergeCells_red[count].size()<<endl;
Info<<"---------------------------------------------_________________-------"<<endl;
                for(int i=tryedCells[count];i<possibleMergeCells_red[count].size();i++,tryedCells[count]++)
                {
Info<<"i: "<<i<<"  Merge Cell:"<<possibleMergeCells_red[count][i]<<endl;

                    bool cellsNotBlocked = true;
                    DynamicList<label> blockedBecausOf;
//Info<<" test Non Blocking "<<endl;
                    for(int s=0;s<possibleMergeCells_red[count][i].size();s++)
                    {
                        if(cellReserved.count(possibleMergeCells_red[count][i][s]) != 0)
                        {
                            auto keyIt = cellReserved.find(possibleMergeCells_red[count][i][s]);
                            if(keyIt == cellReserved.end())
                            {
                                FatalErrorInFunction
                                << " Error1 !"<<endl
                                << exit(FatalError); 
                            }
                            blockedBecausOf.append(keyIt->second);
                            if(keyIt->second>=mergeNecessary.size())
                            {
                                Info<<"cell: "<<keyIt->first<<endl;
                                FatalErrorInFunction
                                << " Error!"<<endl
                                << exit(FatalError);
                            }
                            cellsNotBlocked = false;
                            Info<<"possibleMergeCells_red[count][i][s]:"<<possibleMergeCells_red[count][i][s]<<endl;
                            Info<<"blocked:"<<cellReserved.count(possibleMergeCells_red[count][i][s])<<endl;
                            Info<<"cellsNotBlocked:"<<cellsNotBlocked<<endl;
                        }
                        else
                        {
                            Info<<"possibleMergeCells_red[count][i][s]:"<<possibleMergeCells_red[count][i][s]<<endl;
                            Info<<"not blocked:"<<cellReserved.count(possibleMergeCells_red[count][i][s])<<endl;
                        }
                    }
Info<<"Hello"<<endl;
Info<<"1 cellsNotBlocked:"<<cellsNotBlocked<<endl;
                    cellsAreBlocked[i-tryedCellsStartPoint] = !cellsNotBlocked;
                    if(!cellsNotBlocked) Info<<" are blocked "<<endl;
Info<<"Not hello"<<endl;
                    
                    bool cellsWillNotBlock = true;
                    DynamicList<label> willBlock;
                    cellsThatWillBeBlocked[i-tryedCellsStartPoint].setSize(possibleMergeCells_red[count][i].size());
//Info<<" test will Block "<<endl;
                    // Iterate across all selected mergeCells in one mergeOptions
                    for(int s=0;s<possibleMergeCells_red[count][i].size();s++)
                    {
                        Info<<"s:"<<s<<endl;
//Info<<"count:"<<count<<" Test if "<<possibleMergeCells_red[count][i][s]<<" will block for merge with "<<redIndToCell[count]<<endl;
                        std::unordered_multimap<label,label> mergeOptionLocation; //Cells and mergeOptions that are blocked 
                        DynamicList<label> cellToBeBlocked; //List of cells that might be preBlocked
                    
                        for(auto keyIt = cellPreBlock.find(possibleMergeCells_red[count][i][s]);
                            keyIt != cellPreBlock.end() && keyIt->first == possibleMergeCells_red[count][i][s];
                            keyIt++)
                        {
                            if(((keyIt->second).first > count) && 
                               (mergeOptionLocation.find((keyIt->second).first)==mergeOptionLocation.end()))
                                cellToBeBlocked.append((keyIt->second).first);
                            mergeOptionLocation.insert(keyIt->second);
//Info<<" Insert Cell:"<<(keyIt->second).first<<" Option:"<<(keyIt->second).second<<" from "<<keyIt->first<<endl;
//Info<<"   "<<redIndToCell[(keyIt->second).first]<<"    "<<possibleMergeCells_red[(keyIt->second).first][(keyIt->second).second]<<endl;
                        }
                        
Info<<" cellToBeBlocked "<<cellToBeBlocked<<endl;
                        
                        //Iterate across all cells that might be totally blocked
                        for(int ss=0;ss<cellToBeBlocked.size();ss++)
                        {
Info<<" Is cell "<<cellToBeBlocked[ss]<<"/"<<redIndToCell[cellToBeBlocked[ss]]<<" blocked? ";
                            // select the specific merge option that might be blocked
                            DynamicList<bool> posBlocked = cellMergPosBlocked_red[cellToBeBlocked[ss]];
//Info<<posBlocked<<endl;
                            bool allBlocked = true;
                            std::unordered_set<label> optionsBlocking; //mergeOptions that are blocked
                            for(auto iter = mergeOptionLocation.find(cellToBeBlocked[ss]);
                                iter != mergeOptionLocation.end() && iter->first == cellToBeBlocked[ss];
                                iter++)
                            {
                                optionsBlocking.insert(iter->second);
//Info<<"   Option added: "<<iter->second<<endl;
                            }
Info<<"End"<<endl;
                            //Iterate across all mergeOptions and test if blocked
                            for(int sss=0;sss<posBlocked.size();sss++)
                            {
//Info<<"sss:"<<sss<<" posBlocked[sss]:"<<posBlocked[sss]<<" optionsBlocking:"<<!(optionsBlocking.find(sss)==optionsBlocking.end())<<" blocked:"<<allBlocked<<endl;
                                if(!posBlocked[sss] && (optionsBlocking.find(sss)==optionsBlocking.end()))
                                {
                                    bool localcellsNotBlocked = true;
                                    for(int ssss=0;ssss<possibleMergeCells_red[cellToBeBlocked[ss]][sss].size();ssss++)
                                    {
                                        if(cellReserved.count(possibleMergeCells_red[cellToBeBlocked[ss]][sss][ssss]) != 0)
                                        {
                                            localcellsNotBlocked = false;
                                        }
                                    }
                                    if(localcellsNotBlocked)
                                    {
                                        // Test for cellReserved blocking in the option!!!
                                        allBlocked = false;
                                    }
                                }
                            }
                            if(allBlocked)
                            {
Info<<"Set to Block"<<endl;
Info<<"i: "<<i<<" s: "<<s<<endl<<"cellsThatWillBeBlocked: "<<cellsThatWillBeBlocked<<endl;
Info<<"ss: "<<ss<<endl<<"cellToBeBlocked: "<<cellToBeBlocked<<endl;
                                cellsThatWillBeBlocked[i-tryedCellsStartPoint][s].append(cellToBeBlocked[ss]);
                                cellsWillNotBlock = false;
Info<<"End"<<endl;
                            }
                        }
Info<<"cellsWillNotBlock: "<<cellsWillNotBlock<<endl;
                    }
                    cellsWillBlock[i-tryedCellsStartPoint] = !cellsWillNotBlock;
                    if(!cellsWillNotBlock) Info<<" will block "<<endl;
Info<<"2 cellsNotBlocked:"<<cellsNotBlocked<<endl;
                    Info<<endl;
//Info<<" Ende"<<endl;
                    if(cellsNotBlocked)
                    {
                        allCellsAreBlocked = false;
                    }
                    if(cellsWillNotBlock)
                    {
                        allCellsWillBlock = false;
                    }
                    if(cellsNotBlocked && cellsWillNotBlock)
                    {
//Info<<" not blocked"<<endl;
                        MergeFaceFound = true;

//Info<<"Found face for cell: "<<count<<"  ";
//Info<<"tryedCells["<<count<<"] = "<<tryedCells[count]<<"/"<<"possibleMergeCells["<<count<<"] = "<<possibleMergeCells[count].size()<<endl;

                        mergeFace = possibleMergeFaces_red[count][i];
                        mergeCell = possibleMergeCells_red[count][i];
                        Info<<"Set new "<<count<<endl;
                        Info<<"Merge Face:"<<mergeFace<<endl;
                        Info<<"Merge Cell:"<<mergeCell<<endl;
                        
                        Info<<"3 cellsNotBlocked:"<<cellsNotBlocked<<endl;
                        Info<<"cellsWillNotBlock:"<<cellsWillNotBlock<<endl;

                        tryedCells[count]++;
                        break;
                    }
                    else
                    {
//Info<<"  "<<redIndToCell[count]<<"  "<<possibleMergeCells_red[count][i];
if(!cellsNotBlocked){
//Info<<" blocked because ";
                        for(int x=0;x<blockedBecausOf.size();x++)
                        {
//Info<<"("<<blockedBecausOf[x]<<"|"<<cellToRedInd[blockedBecausOf[x]]<<") ";
                        }
}
                        if(!cellsWillNotBlock)
                        {
//Info<<" will block-----------------------------------------------------------------------"<<endl;
                        }
                        else
                        {
//Info<<endl;
                        }
                    }
                }
                if(MergeFaceFound == false)
                /* Decision B: The iteration across all possibleMergeCells failed. This happens if
                * all the possible merge cells between tryed[count] und possibleMergeCells[count].size
                * are already blocked.
                * The result is a backtracking to the first previous cell where another selection was still
                * possible
                */
                {
Info<<"Not found"<<endl;
/*
Info<<"Merge face not found for "<<count<<endl;
Info<<"tryedCells["<<count<<"] = "<<tryedCells[count]<<"/"<<"possibleMergeCells["<<count<<"] = "<<possibleMergeCells[count].size()<<" assignList["<<count<<"]:"<<assignList[count]<<"  mergeNecessary["<<count<<"]:"<<mergeNecessary[count]<<endl;
Info<<"tryedCells["<<count-1<<"] = "<<tryedCells[count-1]<<"/"<<"possibleMergeCells["<<count-1<<"] = "<<possibleMergeCells[count-1].size()<<" assignList["<<count-1<<"]:"<<assignList[count-1]<<"  mergeNecessary["<<count-1<<"]:"<<mergeNecessary[count-1]<<endl;
*/
/*
 * Wrong code
                    if(allCellsAreBlocked)
                    {
                        FatalErrorInFunction
                        << " Must not happen!"<<endl
                        << exit(FatalError);
                    }
                    if(allCellsWillBlock)
                    {
                    }
*/

                    DynamicList<label> trackBackPoints;
                    DynamicList<label> trackBackPointsMergeInd;
//Info<<"Size:"<<possibleMergeCells_red[count].size()<<endl;
//Info<<"Size2:"<<trackBackPoints.size()<<endl;
                    int z=0;
Info<<"tryedCellsStartPoint:"<<tryedCellsStartPoint<<endl;
Info<<"cellsAreBlocked "<<cellsAreBlocked<<endl;
Info<<"cellsWillBlock "<<cellsWillBlock<<endl;
Info<<"cellsThatWillBeBlocked: "<<cellsThatWillBeBlocked<<endl;
                    for(int s=tryedCellsStartPoint;s<possibleMergeCells_red[count].size();s++,z++)
                    {
Info<<"s:"<<s<<" possibleMergeCells_red:"<<possibleMergeCells_red[count][s]<<endl;
                        if(possibleMergeCells_red[count][s].size() == 0)
                        {
                            FatalErrorInFunction
                            << " Zero length merge section!"<<endl
                            << exit(FatalError);
                        }
                        for(int zz=1;zz<possibleMergeCells_red[count][s].size();zz++)
                        {
                            auto keyIt = cellReserved.find(possibleMergeCells_red[count][s][zz]);
                            if(keyIt != cellReserved.end() && keyIt->first != possibleMergeCells_red[count][s][zz])
                            {
                                FatalErrorInFunction
                                << " Something is wrong here!"<<endl
                                << exit(FatalError);
                            }
                        }
                        
                        label minTrackBackPointFromBlocked = -1;
                        DynamicList<label> trackPointsForMergeOption;
                        if(cellsAreBlocked[z])
                        {
                            for(int zz=0;zz<possibleMergeCells_red[count][s].size();zz++)
                            {
//Info<<z<<"/"<<possibleMergeCells_red[count][s].size()<<":"<<possibleMergeCells_red[count][s][z]<<endl;
                                auto keyIt = cellReserved.find(possibleMergeCells_red[count][s][zz]);
                                if(keyIt != cellReserved.end())
                                {
                                    trackPointsForMergeOption.append(keyIt->second);
                                }
                                else
                                    continue;
                                if(cellToRedInd[keyIt->second] >= count || cellToRedInd[keyIt->second] < 0)
                                {
                                    Info<<"possibleMergeCells_red.size():"<<possibleMergeCells_red.size()<<endl;
                                    Info<<"TrackBackPoint:"<<keyIt->second<<endl;
                                    Info<<"TrackBackPointRed:"<<cellToRedInd[keyIt->second]<<endl;
                                    Info<<"count:"<<count<<endl;
                                    FatalErrorInFunction
                                    << " More general possible Track Back Point is wrong!"<<endl
                                    << exit(FatalError);
                                }           
                            }
                            if(trackPointsForMergeOption.size()==0)
                            {
                                FatalErrorInFunction
                                << " Something is wrong here!"<<endl
                                << exit(FatalError);
                            }                            
                            minTrackBackPointFromBlocked = trackPointsForMergeOption[0];
                            for(int zz=0;zz<trackPointsForMergeOption.size();zz++)
                            {
                                if(minTrackBackPointFromBlocked > trackPointsForMergeOption[zz])
                                    minTrackBackPointFromBlocked = trackPointsForMergeOption[zz];
                            }                        
                            if(cellToRedInd[minTrackBackPointFromBlocked] >= count ||
                               cellToRedInd[minTrackBackPointFromBlocked] < 0)
                            {
                                Info<<"possibleMergeCells_red.size():"<<possibleMergeCells_red.size()<<endl;
                                Info<<"TrackBackPoint:"<<minTrackBackPointFromBlocked<<endl;
                                Info<<"TrackBackPointRed:"<<cellToRedInd[minTrackBackPointFromBlocked]<<endl;
                                Info<<"count:"<<count<<endl;
                                FatalErrorInFunction
                                << " Possible Track Back Point is wrong!"<<endl
                                << exit(FatalError);                                
                            }
                        }
                        //minTrackBackPointFromBlocked is guarenteed to be !=-1 at the end of this path
                        label minTrackBackPointFromWillBlock = -1;
                        if(cellsWillBlock[z])
                        {
                            if(cellsThatWillBeBlocked[z].size()==0)
                            {                           
                                FatalErrorInFunction<<"Option will block but there are no cells that will be blocked! This can not happen!"<<exit(FatalError);
                            }
                            //Iterate across all mergeCells of the mergeOption
                            label minBackTrackingForMergeOption = possibleMergeCells.size();
                            bool minBackTrackingForMergeOptionSET = false;
                            for(int zz=0;zz<cellsThatWillBeBlocked[z].size();zz++)
                            {
                                bool noBlockedCells = true;
                                for(int inter=0;inter<cellsThatWillBeBlocked[z].size();inter++)
                                {
                                    if(cellsThatWillBeBlocked[z][inter].size()!=0)
                                        noBlockedCells = false;
                                }
                                if(noBlockedCells)
                                {                           
                                    FatalErrorInFunction<<"Option will block but there are no cells in the inner lists that will be blocked! This can not happen!"<<exit(FatalError);
                                }
                                
                                label minBackTrackingForOneMergeCell = possibleMergeCells.size();
                                bool minBackTrackingForOneMergeCellSET = false;
                                //Iterate across the cells that will be blocked by each mergeCell
                                for(int zzz=0;zzz<cellsThatWillBeBlocked[z][zz].size();zzz++)
                                {
                                    Info<<"zzz: "<<zzz<<endl;
                                //Begin testing if not all are blocked and the current object blocks; Testing block is self contained!!!
                                    DynamicList<label> nonBlockedMergeOption;
                                    label blockedCell = cellsThatWillBeBlocked[z][zz][zzz];
                                    bool allOptionsBlocked = true;
                                    //Iterate across all options of cells that will be blocked
                                    for(int g=0;g<cellMergPosBlocked_red[blockedCell].size();g++)
                                    {
                                        if(!cellMergPosBlocked_red[blockedCell][g])
                                        {
                                            nonBlockedMergeOption.append(g);
                                            allOptionsBlocked = false;
                                        }
                                    }
                                    if(allOptionsBlocked)
                                    {
                                        FatalErrorInFunction
                                        << " It can not happen that all options are already blocked! That should be prevented previously!"<<endl
                                        << exit(FatalError);  
                                    }
                                    bool allNonBlockedOptionsWillBeBlocked = true;
                                    for(int g=0;g<nonBlockedMergeOption.size();g++)
                                    {
                                        std::unordered_multimap<label,label> mergeOptionLocation; //Cells and mergeOptions that are blocked 
                                        DynamicList<label> cellToBeBlocked; //List of cells that might be preBlocked
                    
                                        for(auto keyIt = cellPreBlock.find(possibleMergeCells_red[count][s][z]);
                                            keyIt != cellPreBlock.end() && keyIt->first == possibleMergeCells_red[count][s][z];
                                            keyIt++)
                                        {
                                                if(((keyIt->second).first > count) && 
                                                    (mergeOptionLocation.find((keyIt->second).first)==mergeOptionLocation.end()))
                                                    cellToBeBlocked.append((keyIt->second).first);
                                                mergeOptionLocation.insert(keyIt->second);
                                        }
                                        bool oneBlocks = false;
                                        for(auto keyIt = mergeOptionLocation.find(blockedCell);
                                            keyIt != mergeOptionLocation.end() && keyIt->first == blockedCell;
                                            keyIt++)
                                            {
                                                if(keyIt->second == nonBlockedMergeOption[g])
                                                    oneBlocks = true;
                                            }
                                        if(!oneBlocks)
                                        {
                                            for(int gg=0;gg<possibleMergeCells_red[blockedCell][nonBlockedMergeOption[g]].size();gg++)
                                            {
                                                if(cellReserved.count(possibleMergeCells_red[blockedCell][nonBlockedMergeOption[g]][gg]) != 0)
                                                {
                                                    allNonBlockedOptionsWillBeBlocked = false;
                                                }
                                            }
                                        }
                                        // All options in nonBlockedMergeOption that are not already blocked
                                        // have to be blocked by the current mergeOption.
                                        // If not: failure
                                    }
                                    if(!allNonBlockedOptionsWillBeBlocked)
                                    {
                                        FatalErrorInFunction
                                        << " It can not happen that the non blocked options will not be blocked by the current cell!"<<endl
                                        << exit(FatalError);  
                                    }
                                //End testing

                                //Start finding trackBackPoints
                                    //Iterate across all merging options of the willBeBlocked cell
                                    label willBeBlockedCellTrackBack = -1;
                                    bool willBeBlockedCellTrackBackSET = false;
                                    for(int g=0;g<possibleMergeCells_red[blockedCell].size();g++)
                                    {
                                        DynamicList<label> reservedBackPoints;
                                        for(int gg=0;gg<possibleMergeCells_red[blockedCell][g].size();gg++)
                                        {
                                            auto keyIt = cellReserved.find(possibleMergeCells_red[blockedCell][g][gg]);
                                            if(keyIt != cellReserved.end())
                                            {
                                                reservedBackPoints.append(keyIt->second);
                                            }
                                        }
                                        label mergeOptionTrackBack = -1;
                                        label reservedBackPoint = possibleMergeCells.size();
                                        for(int gg=0;gg<reservedBackPoints.size();gg++)
                                        {
                                            if(reservedBackPoints[gg]<reservedBackPoint)
                                                reservedBackPoint = reservedBackPoints[gg];
                                        }
                                        label blockedBackPoint = cellMergPosBlocked_red_Reason[blockedCell][g];
                                        mergeOptionTrackBack = (reservedBackPoint<blockedBackPoint)?reservedBackPoint:blockedBackPoint;
                                        if(blockedBackPoint!=possibleMergeCells.size() || reservedBackPoint!=possibleMergeCells.size())
                                        {
                                            if(mergeOptionTrackBack>=possibleMergeCells.size() || mergeOptionTrackBack<0)
                                            {
                                                Info<<endl<<"mergeOptionTrackBack: "<<mergeOptionTrackBack<<endl;
                                                Info<<"blockedBackPoint: "<<blockedBackPoint<<endl;
                                                Info<<"reservedBackPoint: "<<reservedBackPoint<<endl;
                                                Info<<"possibleMergeCells.size(): "<<possibleMergeCells.size()<<endl;
                                                Info<<"cellsThatWillBeBlocked[z]:"<<cellsThatWillBeBlocked[z]<<endl;
                                                Info<<"cellsThatWillBeBlocked[z][zz]:"<<cellsThatWillBeBlocked[z][zz]<<endl;
                                                Info<<"blockedCell: "<<blockedCell<<endl;
                                                Info<<"cellMergPosBlocked_red_Reason[blockedCell]:"<<cellMergPosBlocked_red_Reason[blockedCell]<<endl;
                                                Info<<"reservedBackPoints: "<<reservedBackPoints<<endl;
                                                Info<<"The initial value can not still be there!"<<endl;
                                                Info<<"mergeOptionTrackBack: "<<endl;
                                                FatalErrorInFunction
                                                << " Track Back Point across merging Option of willBeBlocked cell is out of range! "<<endl<<exit(FatalError);
                                            }
                                            willBeBlockedCellTrackBack = (willBeBlockedCellTrackBack<mergeOptionTrackBack)?mergeOptionTrackBack:willBeBlockedCellTrackBack;
                                            if(willBeBlockedCellTrackBack<mergeOptionTrackBack)
                                                willBeBlockedCellTrackBackSET=true;
                                        }
                                    }
                                    if(willBeBlockedCellTrackBackSET)
                                    {
                                        if(willBeBlockedCellTrackBack>=possibleMergeCells.size() || willBeBlockedCellTrackBack<0)
                                        {
                                            Info<<endl<<"willBeBlockedCellTrackBack: "<<willBeBlockedCellTrackBack<<endl;
                                            Info<<"possibleMergeCells.size(): "<<possibleMergeCells.size()<<endl;
                                            Info<<"Has to be != -1 because mergeOptionTrackBack is a valid value!"<<endl;
                                            FatalErrorInFunction
                                            << " Track Back Point of willBeBlocked cell is out of range! "<<endl
                                            << exit(FatalError);
                                        }
                                        Info<<"-----------------------"<<willBeBlockedCellTrackBack<<endl;
                                        minBackTrackingForOneMergeCell = (minBackTrackingForOneMergeCell>willBeBlockedCellTrackBack)?  willBeBlockedCellTrackBack:minBackTrackingForOneMergeCell;
                                        if(minBackTrackingForOneMergeCell>willBeBlockedCellTrackBack)
                                            minBackTrackingForOneMergeCellSET = true;
                                    
                                        Info<<endl<<"willBeBlockedCellTrackBack: "<<willBeBlockedCellTrackBack<<" | minBackTrackingForOneMergeCell: "<<minBackTrackingForOneMergeCell<<endl;
                                    }
                                }
                                if(minBackTrackingForOneMergeCellSET)
                                {
                                    if(minBackTrackingForOneMergeCell>=possibleMergeCells.size() || minBackTrackingForOneMergeCell<0)
                                    {
                                        Info<<endl<<"cellsThatWillBeBlocked.size(): "<<cellsThatWillBeBlocked.size()<<endl;
                                        Info<<"cellsThatWillBeBlocked[z].size(): "<<cellsThatWillBeBlocked[z].size()<<endl;
                                        Info<<"cellsThatWillBeBlocked[z][zz].size(): "<<cellsThatWillBeBlocked[z][zz].size()<<endl;
                                        Info<<"minBackTrackingForOneMergeCellSET: "<<minBackTrackingForOneMergeCellSET<<endl;
                                        Info<<"minBackTrackingForOneMergeCell: "<<minBackTrackingForOneMergeCell<<endl;
                                        Info<<"possibleMergeCells.size(): "<<possibleMergeCells.size()<<endl;
                                        FatalErrorInFunction
                                        << " Track Back Point of all willBeBlocked cell is out of range! "<<endl
                                        << exit(FatalError);
                                    }
                                    minBackTrackingForMergeOption = (minBackTrackingForMergeOption<minBackTrackingForOneMergeCell)?minBackTrackingForMergeOption:minBackTrackingForOneMergeCell;
                                    if(!(minBackTrackingForMergeOption<minBackTrackingForOneMergeCell))
                                        minBackTrackingForMergeOptionSET = true;
                                }
                            }
                            if(minBackTrackingForMergeOptionSET)
                            {
                                if(minBackTrackingForMergeOption>=possibleMergeCells.size() || minBackTrackingForMergeOption<0)
                                {
                                    Info<<endl<<"minBackTrackingForMergeOption: "<<minBackTrackingForMergeOption<<endl;
                                    Info<<"possibleMergeCells.size(): "<<possibleMergeCells.size()<<endl;
                                    FatalErrorInFunction
                                    << " Track Back Point of all mergeOption is out of range! "<<endl
                                    << exit(FatalError);
                                }
                                minTrackBackPointFromWillBlock = minBackTrackingForMergeOption;
                            }
                            else
                            {
                                FatalErrorInFunction
                                << "No assignment to minBackTrackingForMergeOption! This can not happen!"<<endl
                                << exit(FatalError);
                            }
                        }
                        if(minTrackBackPointFromBlocked == -1 && minTrackBackPointFromWillBlock == -1)
                        {
                            FatalErrorInFunction
                            << " Neither blocked nor will block. Can not happen if nor merge face found!"<<endl
                            << exit(FatalError);
                        }
                        label min = -1;
                        if(minTrackBackPointFromBlocked != -1 && minTrackBackPointFromWillBlock != -1)
                        {
                            min = (minTrackBackPointFromBlocked
                                   < minTrackBackPointFromWillBlock)?
                                   minTrackBackPointFromBlocked :
                                   minTrackBackPointFromWillBlock;
                        }
                        else if(minTrackBackPointFromBlocked != -1)
                        {
                            min = minTrackBackPointFromBlocked;
                        }
                        else if(minTrackBackPointFromWillBlock != -1)
                        {
                            min = minTrackBackPointFromWillBlock;
                        }
                        else
                            FatalErrorInFunction<< " No track back point found! "<< exit(FatalError);
                        Info<<"minTrackBackPointFromWillBlock:"<<minTrackBackPointFromWillBlock<<"  minTrackBackPointFromBlocked:"<<minTrackBackPointFromBlocked<<endl;
                        if(min==-1)
                        {
                            FatalErrorInFunction
                            << " To append set back point is -1. Can not happen!"<<endl
                            << exit(FatalError);
                        }
                        trackBackPoints.append(min);
                        trackBackPointsMergeInd.append(s);
                    }

Info<<"cellToRedInd.size:"<<cellToRedInd.size()<<endl;
Info<<"trackBackPoints:";
for(int s=0;s<trackBackPoints.size();s++) Info<<"  "<<trackBackPoints[s]<<"/"<<cellToRedInd[trackBackPoints[s]];
Info<<endl;

                    if(trackBackPoints.size() < 1)
                    {
                        Info<<"trackBackPoints: "<<trackBackPoints<<endl;
                        FatalErrorInFunction
                        << " Something wrong!"<<endl
                        << exit(FatalError);
                    }

                    int js;
                    label keyMergeInd,keyBackPoint;
                    for(int s=1;s<trackBackPoints.size();s++)
                    {
                        keyBackPoint = trackBackPoints[s];
                        keyMergeInd = trackBackPointsMergeInd[s];
                        js = s-1;
                        while(js>=0 && trackBackPoints[js] < keyBackPoint)
                        {
                            trackBackPoints[js+1] = trackBackPoints[js];
                            trackBackPointsMergeInd[js+1] = trackBackPointsMergeInd[js];
                            js--;
                        }
                        trackBackPoints[js+1] = keyBackPoint;
                        trackBackPointsMergeInd[js+1] = keyMergeInd;
                    }
/*
for(int ijk=0;ijk<=count;ijk++)
{
    Info<<"redInd:"<<ijk<<" cellInd:"<<redIndToCell[ijk]<<" "<<tryedCells[ijk]<<"/"<<possibleMergeCells_red[ijk].size()<<endl;
}
*/                    
Info<<"Sort: tBP:"<<trackBackPoints<<" tBPMI:"<<trackBackPointsMergeInd<<endl;

                    label bestTrackBackPoint = trackBackPoints[trackBackPoints.size()-1];
/*
if(count == 90){
for(int s=0;s<trackBackPoints.size();s++)
{
Info<<"----------"<<s<<"----------------------"<<endl;
Info<<"Red Track back: "<<cellToRedInd[trackBackPoints[s]]<<" Merge:"<<possibleMergeCells_red[cellToRedInd[trackBackPoints[s]]]<<" tryed:"<<tryedCells[cellToRedInd[trackBackPoints[s]]]<<endl<<"Cell:"<<trackBackPoints[s]<<endl;
Info<<"----------"<<s<<"----------------------"<<endl;
}
Info<<"--------------------------------"<<endl;
Info<<"Red Track back: "<<count<<" Merge:"<<possibleMergeCells_red[count]<<"Cell:"<<redIndToCell[count]<<endl;
Info<<"--------------------------------"<<endl;
}
*/
                    bool oneTrackBackPointUnblocks = false;
                    for(int s=0;s<trackBackPoints.size();s++)
                    {
                        bool trackBackCanFreeCount = false;
                        label redIntBackPoint = cellToRedInd[trackBackPoints[s]];
                        label cellBackPoint = trackBackPoints[s];
                        label tryContPoint = tryedCells[redIntBackPoint];
//Info<<"redInt"<<redIntBackPoint<<" tryContPoint:"<<tryContPoint<<" possibleMergeCells_red["<<redIntBackPoint<<"].size()="<<possibleMergeCells_red[redIntBackPoint].size()<<endl;
                        for(int ss=tryContPoint;ss<possibleMergeCells_red[redIntBackPoint].size();ss++)
                        {
                            bool oneMergePossAtBlock = mergeCellSelectionBlocks
                                (cellBackPoint,redIntBackPoint,possibleMergeCells_red[redIntBackPoint][ss],                                                                   redIndToCell[count],count,possibleMergeCells_red[count],cellReserved);
                                
//Info<<ss<<" Blocks: "<<trackBackPoints[s]<<" Bool:"<<oneMergePossAtBlock<<endl;
                            
                            if(oneMergePossAtBlock)
                            {
                                trackBackCanFreeCount = true;
                                tryedCells[redIntBackPoint] = ss;
                                oneTrackBackPointUnblocks = true;
//Info<<"Non block: "<<trackBackPoints[s]<<endl;
                                break;
                            }
                        }
                        if(trackBackCanFreeCount)
                        {
                            bestTrackBackPoint = trackBackPoints[s];
                            break;
                        }
                    }
//Info<<"btbp: "<<bestTrackBackPoint<<endl;
//if(count == 90){ Info<<"bestTrackBackPoint:"<<bestTrackBackPoint<<endl;}
                    if(!oneTrackBackPointUnblocks && false)
                    {
                        for(;; bestTrackBackPoint--)
                        {
                            if(bestTrackBackPoint == -1)
                            {
                                Info<<"Backtracking from cell: "<<bestTrackBackPoint<<endl;
                                FatalErrorInFunction
                                << "Failed in Merging Selection! There is no found combination for all merging cells."<<endl
                                << exit(FatalError);
                            }
                            label redIndCntBck = cellToRedInd[bestTrackBackPoint];
                            if((redIndCntBck==-1 && assignList[bestTrackBackPoint][0]!=-3) || (assignList[bestTrackBackPoint][0]==-3 && redIndCntBck!=-1))
                            {
                                FatalErrorInFunction
                                << " Can not happen!"<<endl
                                << exit(FatalError);
                            }
                            if(redIndCntBck==-1)
                                continue;
                            
                            bool trackBackCanFreeCount = false;
                            label tryContPoint = tryedCells[redIndCntBck];
                            for(int ss=tryContPoint;ss<possibleMergeCells_red[redIndCntBck].size();ss++)
                            {
                                bool oneMergePossAtBlock = mergeCellSelectionBlocks
                                (bestTrackBackPoint,redIndCntBck,possibleMergeCells_red[redIndCntBck][ss],                                                                   redIndToCell[count],count,possibleMergeCells_red[count],cellReserved);
                                                            
                                if(oneMergePossAtBlock)
                                {
                                    trackBackCanFreeCount = true;
                                    tryedCells[redIndCntBck] = ss;
                                    break;
                                }
                            }
                            if(trackBackCanFreeCount)
                            {
                                break;
                            }
                        }                        
                    }
                        
                    bestTrackBackPoint = cellToRedInd[bestTrackBackPoint];
/*
if(count == 90){ 
    Info<<"bestTrackBackPointRed:"<<bestTrackBackPoint<<endl;
    FatalErrorInFunction<<"Temporary stop!"<<endl<<exit(FatalError);
}
*/

//Info<<"bestTrackBackPoint:"<<bestTrackBackPoint<<endl;
                    if(bestTrackBackPoint >= count || bestTrackBackPoint < 0)
                    {
                        Info<<"bestTrackBackPoint:"<<bestTrackBackPoint<<endl;
                        Info<<"count:"<<count<<endl;
                        FatalErrorInFunction
                        << " Track Back Point is wrong!"<<endl
                        << exit(FatalError);
                    }
                    //////////////////////////////////////////----------------------------

                    int backtrackingIndex = -1;                    
                    for(int cntBck = bestTrackBackPoint; cntBck>=0; cntBck--)
                    {
                        if(assignList[redIndToCell[cntBck]].size() == 1 &&
                           assignList[redIndToCell[cntBck]][0] == -3)
                        {
                            FatalErrorInFunction
                            << " Backtracking to untreated cell!"<<endl
                            << exit(FatalError);
                        }
                        if((assignList[redIndToCell[cntBck]].size() == 1 && 
                            assignList[redIndToCell[cntBck]][0] == -1)   ||
                            !mergeNecessary_red[cntBck])
                        {
                            FatalErrorInFunction
                            << " Backtracking to nonmerge Cell or -1 assign cell!"<<endl
                            << exit(FatalError);
                        }
                        if(assignList[redIndToCell[cntBck]].size() == 0)
                        {
                            FatalErrorInFunction
                            << " Backtracking empty assign list!"<<endl
                            << exit(FatalError);
                        }
                        
                        if((assignList[redIndToCell[cntBck]].size() > 0 && 
                            assignList[redIndToCell[cntBck]][0] != -2) && tryedCells[cntBck]<possibleMergeCells_red[cntBck].size())
                        {
                            backtrackingIndex = cntBck;
                            break;
                        }
                    }
                    if(backtrackingIndex != -1)
                    /* Decision C: If the backtracking resulted in a cell the backtracking was succesful.
                    * The iteration continues with this specific cell. 
                    */
                    {
                        if(!mergeNecessary_red[count])
                        {
                            FatalErrorInFunction
                            << " Backtracking to none merge cell!"<<endl
                            << exit(FatalError);
                        }
                        if(assignList[redIndToCell[backtrackingIndex]].size()==0)
                        {
                            FatalErrorInFunction
                            << " Backtracking to zero length assign List can not happen!"<<endl
                            << exit(FatalError);
                        }
                        if(assignList[redIndToCell[backtrackingIndex]][0] < 0)
                        {
                            Info<<endl<<"assignList["<<backtrackingIndex<<"][0]="<<assignList[backtrackingIndex][0]<<endl;
                            FatalErrorInFunction
                            << " Backtracking to cell with assign list < 0!"<<endl
                            << exit(FatalError);
                        }
                        if(tryedCells[backtrackingIndex] <= 0)
                        {
                            Info<<"tryedCells["<<backtrackingIndex<<"]="<<tryedCells[backtrackingIndex]<<endl;
                            FatalErrorInFunction
                            << " Backtracking to cell tryedCells <= 0!"<<endl
                            << exit(FatalError);
                        }
                        
                        for(int k=backtrackingIndex; k<=count;k++)
                        {
/*
Info<<"<<<<<";
*/
                            for(int l=0;l<blockedCells[k].size();l++)
                            {
                                cellReserved.erase(blockedCells[k][l]);
/*
Info<<">>>>>>>>>>>Cleared["<<k<<"]:"<<blockedCells[k][l]<<endl;
*/
                            }
                            for(int l=0;l<preBlockedOptions_red[k].size();l++)
                            {
                                label cellCount_red = preBlockedOptions_red[k][l].first;
                                label mergeSelectionCount = preBlockedOptions_red[k][l].second;
                                
                                if(!cellMergPosBlocked_red[cellCount_red][mergeSelectionCount] ||
                                    (cellMergPosBlockedMulti_red[cellCount_red][mergeSelectionCount] <= 0))
                                {
                                    FatalErrorInFunction
                                    << "Cell is either not blocked or the mergePos multiplcity is smaller than one! Can not happen"<<endl
                                    << exit(FatalError);
                                }
                                if(cellMergPosBlocked_red_Reason[cellCount_red][mergeSelectionCount] == 
                                    possibleMergeCells.size())
                                {
                                    FatalErrorInFunction
                                    << "mergePos multiplcity is set to initial value! Can not happend."<<endl
                                    << exit(FatalError);
                                }
                                
                                
                                
                                cellMergPosBlockedMulti_red[cellCount_red][mergeSelectionCount]--;
                                if(cellMergPosBlockedMulti_red[cellCount_red][mergeSelectionCount] == 0)
                                {
                                    cellMergPosBlocked_red[cellCount_red][mergeSelectionCount] = false;
                                    cellMergPosBlocked_red_Reason[cellCount_red][mergeSelectionCount] = possibleMergeCells.size();
                                }
                            }
                            preBlockedOptions_red[k].setSize(0);
                            for(int l=0;l<preCellMergDone_red[k].size();l++)
                            {
                                if(!cellMergDone_red[cellToRedInd[preCellMergDone_red[k][l]]])
                                {
                                    FatalErrorInFunction
                                    << "Can not happen!"<<endl
                                    << exit(FatalError);
                                }
                                cellMergDoneMult_red[cellToRedInd[preCellMergDone_red[k][l]]]--;
                                if(cellMergDoneMult_red[cellToRedInd[preCellMergDone_red[k][l]]] == 0)
                                    cellMergDone_red[cellToRedInd[preCellMergDone_red[k][l]]] = false;                                
                            }
                            preCellMergDone_red[k].setSize(0);
                        }
                        //Clean list of blockedCells
                        for(int k=backtrackingIndex+1;k<=count;k++)
                        {
                            tryedCells[k] = 0;
                            blockedCells[k].setSize(0);
/*
Info<<"Clear tryedCells and blockedCells ["<<k<<"]"<<endl;
*/
                                    
                        }
                        for(int o=backtrackingIndex+1;o<=count;o++)
                        {
                            if(tryedCells[o] != 0)
                            {
                                Info<<endl<<endl<<"Backtracking from "<<redIndToCell[count]<<" to "<<redIndToCell[backtrackingIndex]<<endl;
                                Info<<"tryedCells["<<o<<"]="<<tryedCells[o]<<endl;
                                FatalErrorInFunction
                                << " Backtracking to cell while tryedCells is unequal 0!"<<endl
                                << exit(FatalError);
                            }
                        }                        
                        
                        
                        blockedCells[backtrackingIndex].setSize(0);
                        
                        count = backtrackingIndex;
/*
Info<<"Go to "<<count<<endl;
*/
                    }
                    else
                    /* Decision C: The backtracking did not found a backtracking cell. The result is an abort.
                    */
                    {
                        Info<<"Backtracking from cell: "<<redIndToCell[count]<<endl;
                        FatalErrorInFunction
                        << " Failed in Merging Selection! There is no found combination for all merging cells."<<endl
                        << exit(FatalError);
                    }
                }
                else
                /* Decision B: A possible merge cell was found. As a result the two cells to merge are
                * inserted both in the cellReserved map as well as in the blockedCells list             * 
                */
                {
                    if(redIndToCell[count]>=mergeNecessary.size() || count>=possibleMergeCells_red.size())
                    {
                        Info<<"Inserted Backtracking number from count: "<<count<<endl;
                        Info<<"redIndToCell[count]: "<<redIndToCell[count]<<endl;
                        Info<<"Inversion: "<<cellToRedInd[redIndToCell[count]]<<endl;
                        FatalErrorInFunction
                        << " Temporary stop."<<endl
                        << exit(FatalError);
                    }
                    
                    /* Set the reservation table cellReserved
                     * and the list blockedCells for the too small cell
                     */
                    std::pair<label,label> ins1(redIndToCell[count],redIndToCell[count]);
                    cellReserved.insert(ins1);
                    blockedCells[count].append(redIndToCell[count]);

                    /* Set the reservation table cellReserved
                     * and the list blockedCells for the merging cells
                     */
                    for(int u=0;u<mergeCell.size();u++)
                    {
                        if(redIndToCell[count]>=mergeNecessary.size() || count>=possibleMergeCells_red.size())
                        {
                            Info<<"Inserted Backtracking number from count: "<<count<<endl;
                            Info<<"redIndToCell[count]: "<<redIndToCell[count]<<endl;
                            Info<<"Inversion: "<<cellToRedInd[redIndToCell[count]]<<endl;
                            FatalErrorInFunction
                            << " Temporary stop."<<endl
                            << exit(FatalError);
                        }
                        std::pair<label,label> ins2(mergeCell[u],redIndToCell[count]);
                        cellReserved.insert(ins2);
                        blockedCells[count].append(mergeCell[u]);
                    }
                    
                    /* Set the Lists cellMergPosBlocked_red, cellMergPosBlockedMulti_red,
                     * preBlockedOptions_red, cellMergDone_red, cellMergDoneMult_red
                     * and preCellMergDone_red for the merging cells
                     */                    
                    for(int u=0;u<mergeCell.size();u++)
                    {
                        for(auto keyIter=cellPreBlock.find(mergeCell[u]);
                            keyIter!=cellPreBlock.end() && keyIter->first == mergeCell[u];
                            keyIter++)
                        {
                            label cellCount_red = (keyIter->second).first;
                            label mergeSelectionCount = (keyIter->second).second;
                            
                            if((!cellMergPosBlocked_red[cellCount_red][mergeSelectionCount] && 
                               cellMergPosBlockedMulti_red[cellCount_red][mergeSelectionCount]>0) ||
                               (cellMergPosBlocked_red[cellCount_red][mergeSelectionCount] && 
                               cellMergPosBlockedMulti_red[cellCount_red][mergeSelectionCount]==0)
                            )
                            {
                                FatalErrorInFunction
                                << "Can not happen!"<<endl
                                << exit(FatalError);
                            }
                            
                            cellMergPosBlocked_red[cellCount_red][mergeSelectionCount] = true;
                            cellMergPosBlockedMulti_red[cellCount_red][mergeSelectionCount]++;
                            preBlockedOptions_red[count].append(keyIter->second);
                            
                            if(cellMergPosBlocked_red_Reason[cellCount_red][mergeSelectionCount]>redIndToCell[count])
                            {
                                cellMergPosBlocked_red_Reason[cellCount_red][mergeSelectionCount] = redIndToCell[count];
                            }
                        }
                        if(cellToRedInd[mergeCell[u]] != -1)
                        {
                            if((!cellMergDone_red[cellToRedInd[mergeCell[u]]] && 
                               cellMergDoneMult_red[cellToRedInd[mergeCell[u]]]>0) ||
                               (cellMergDone_red[cellToRedInd[mergeCell[u]]] && 
                               cellMergDoneMult_red[cellToRedInd[mergeCell[u]]]==0)
                            )
                            {
                                FatalErrorInFunction
                                << "Can not happen!"<<endl
                                << exit(FatalError);
                            }
                            cellMergDone_red[cellToRedInd[mergeCell[u]]] = true;
                            cellMergDoneMult_red[cellToRedInd[mergeCell[u]]]++;
                            preCellMergDone_red[count].append(mergeCell[u]);
                        }
                    }
                    
                    /* Set the Lists cellMergPosBlocked_red, cellMergPosBlockedMulti_red,
                     * preBlockedOptions_red, cellMergDone_red, cellMergDoneMult_red
                     * and preCellMergDone_red for the too small cell
                     */     
                    for(auto keyIter=cellPreBlock.find(redIndToCell[count]);
                        keyIter!=cellPreBlock.end() && keyIter->first == redIndToCell[count];
                        keyIter++)
                    {
                        label cellCount_red = (keyIter->second).first;
                        label mergeSelectionCount = (keyIter->second).second;
                        
                        if((!cellMergPosBlocked_red[cellCount_red][mergeSelectionCount] && 
                            cellMergPosBlockedMulti_red[cellCount_red][mergeSelectionCount]>0) ||
                            (cellMergPosBlocked_red[cellCount_red][mergeSelectionCount] && 
                            cellMergPosBlockedMulti_red[cellCount_red][mergeSelectionCount]==0)
                        )
                        {
                            FatalErrorInFunction
                            << "Can not happen!"<<endl
                            << exit(FatalError);
                        }
                            
                        cellMergPosBlocked_red[cellCount_red][mergeSelectionCount] = true;
                        cellMergPosBlockedMulti_red[cellCount_red][mergeSelectionCount]++;
                        preBlockedOptions_red[count].append(keyIter->second);
                        
                        if(cellMergPosBlocked_red_Reason[cellCount_red][mergeSelectionCount]>redIndToCell[count])
                        {
                            cellMergPosBlocked_red_Reason[cellCount_red][mergeSelectionCount] = redIndToCell[count];
                        }
                    }
                    if(cellToRedInd[redIndToCell[count]] != -1)
                    {
                        if((!cellMergDone_red[cellToRedInd[redIndToCell[count]]] && 
                            cellMergDoneMult_red[cellToRedInd[redIndToCell[count]]]>0) ||
                            (cellMergDone_red[cellToRedInd[redIndToCell[count]]] && 
                            cellMergDoneMult_red[cellToRedInd[redIndToCell[count]]]==0)
                        )
                        {
                            FatalErrorInFunction
                            << "Can not happen!"<<endl
                            << exit(FatalError);
                        }
                        cellMergDone_red[cellToRedInd[redIndToCell[count]]] = true;
                        cellMergDoneMult_red[cellToRedInd[redIndToCell[count]]]++;
                        preCellMergDone_red[count].append(redIndToCell[count]);
                    }

//Info<<"count:"<<count<<" Cell :"<<redIndToCell[count]<<" merged with: "<<mergeCell<<endl;
//Info<<">>>> Added "<<count<<" and "<<mergeCell<<" to blockedCells["<<count<<"]"<<endl;

                    assignList[redIndToCell[count]] = mergeFace;
                    
/*
if(redIndToCell[count] == 838)
{
    Info<<endl;
    Info<<"Merge Face:"<<mergeFace<<endl;
    //Info<<"possibleMergeFace_red["<<count<<"]:"<<possibleMergeFaces_red[count]<<endl;
                        
    Info<<endl;
    Info<<"Merge Cell:"<<mergeCell<<endl;
    //Info<<"possibleMergeCells_red["<<count<<"]:"<<possibleMergeCells_red[count]<<endl;
                        
    FatalErrorInFunction<< "Temporary stop!"<<exit(FatalError);
}
*/
                    
                    count++;
                }
            }
            else
            /* Decision B: Enters else block for cells that are already used for merge.
            */
            {
                DynamicList<label> temp;
                temp.append(-2);
                assignList[redIndToCell[count]] = temp;
                tryedCells[count] = 0;
                count++;
            }
        }
        else
        /* Decision A: Enters else block for cells that are not too small. The assignList is
         * filled with -1 for these cells.
         */
        {
            DynamicList<label> temp;
            temp.append(-1);
            assignList[redIndToCell[count]] = temp;
            tryedCells[count] = 0;
            count++;
        }
    }
    
    Info<<"Fin"<<endl;
    
//Test
    std::unordered_set<label> allUsedCells;
    for(int i=0;i<assignList.size();i++)
    {
        if(mergeNecessary[i])
        {
            if(assignList[i].size()==0)
            {
                FatalErrorInFunction
                << " Error empty assign List at i:"<<i<<endl
                << exit(FatalError);
            }
            if(assignList[i][0] == -3)
            {
                FatalErrorInFunction
                << " Too small cell was not treated by backtracking algorithm!"<<endl
                << exit(FatalError);
            }            
            if(assignList[i][0] == -1)
            {
                FatalErrorInFunction
                << " Too small cell was listed as large enough!"<<endl
                << exit(FatalError);
            }
            
            std::unordered_multiset<label> mergeCellsSet;
            DynamicList<label> mergeCellsList;
            for(int k=0;k<assignList[i].size();k++)
            {
                if(assignList[i][k] >= neighbour.size())
                {
                    FatalErrorInFunction<<"Boundary face is merge face"<<exit(FatalError);
                }
                label ownerCell = owner[assignList[i][k]];
                label neighborCell = neighbour[assignList[i][k]];
                if(mergeCellsSet.find(ownerCell) == mergeCellsSet.end())
                    mergeCellsList.append(ownerCell);
                if(mergeCellsSet.find(neighborCell) == mergeCellsSet.end())
                    mergeCellsList.append(neighborCell);
                mergeCellsSet.insert(ownerCell);
                mergeCellsSet.insert(neighborCell);
            }
            if(mergeCellsList.size()<=0)
                FatalErrorInFunction<<"Merging option with no cells!"<<exit(FatalError);
            
            label mergeCellMult = mergeCellsSet.count(mergeCellsList[0]);
            for(int k=1;k<mergeCellsList.size();k++)
            {
                if(mergeCellsSet.count(mergeCellsList[k]) != static_cast<long unsigned int>(mergeCellMult))
                {
                    FatalErrorInFunction<<"Different multiplcity of cells in "<<
                    mergeCellsList.size()<<" cell merging "<<exit(FatalError);
                }
            }
            
            if(mergeCellMult==1 && mergeCellsList.size() == 2)
            {
                //2 Cell merging
                if(assignList[i].size() != 1)
                    FatalErrorInFunction<<"Error in 2 Cell merging! "<<exit(FatalError);
            }
            else if(mergeCellMult==2 && mergeCellsList.size() == 4)
            {
                //4 Cell merging
                if(assignList[i].size() != 4)
                    FatalErrorInFunction<<"Error in 4 Cell merging! "<<exit(FatalError);
            }
            else if(mergeCellMult==3 && mergeCellsList.size() == 8)
            {
                //4 Cell merging
                if(assignList[i].size() != 12)
                    FatalErrorInFunction<<"Error in 8 Cell merging! "<<exit(FatalError);
            }
            else
            {
                Info<<"mergeCellMult: "<<mergeCellMult<<endl;
                Info<<"optionMergeCellsList.size() == "<<mergeCellsList.size()<<endl;
                Info<<"possibleMergeFaces[i][j].size() == "<<assignList[i].size()<<endl;
                FatalErrorInFunction<<"Inconsistent merge Option!"<<exit(FatalError);
            }
            for(int k=0;k<mergeCellsList.size();k++)
            {
                if(allUsedCells.find(mergeCellsList[k]) != allUsedCells.end())
                {
                    Info<<"All Cells: "<<mergeCellsList<<endl;
                    Info<<"Cell: "<<mergeCellsList[k]<<" used twice!"<<endl;
                    Info<<"Error at Indx:"<<i<<endl;
                    FatalErrorInFunction<<"Merging cell already taken!"<<exit(FatalError);
                }
                else
                {
                    allUsedCells.insert(mergeCellsList[k]);
                }
            }
        }
        else
        {
            if(assignList[i].size()==0)
            {
                FatalErrorInFunction
                << " Error empty assign List at i:"<<i<<endl
                << exit(FatalError);
            }
            if(assignList[i].size()!=1)
            {
                FatalErrorInFunction
                << " Error overfull assign List at i:"<<i<<endl
                << exit(FatalError);
            }
            if(assignList[i][0] != -3)
            {
                FatalErrorInFunction
                << " Backtracking algorithm wrote inside large enough cell. Something is wrong here!"<<endl
                << exit(FatalError);
            }
            assignList[i][0] = -1;
        }
    }
//Ende Test
    
    /*
    labelList assList(assignList.size());
    for(int i=0;i<assignList.size();i++)
    {
        if(assignList[i].size() != 1)
        {
            FatalErrorInFunction
            << " Error!"<<endl
            << exit(FatalError);
        }
        assList[i] = assignList[i][0];
    }
    return assList;
    */
    return assignList;
    /*
    Label index
    -1 : Cell is not too small
    -2 : cell is already merged by other cell
    -3 : No assignment
    */
}

List<DynamicList<label>> Foam::cutCellFvMesh::assignMergeFaces
(
    const labelList& owner,
    const labelList& neighbour,
    DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
    DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
    DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
    DynamicList<bool>& oneMergeFaceSufficient,
    DynamicList<bool>& mergeNecessary,
    DynamicList<DynamicList<scalar>>& possibleMergeCellsPartialSize,
    scalar partialThreeshold
)
{
    if(possibleMergeFaceArea.size() != possibleMergeFaces.size() || possibleMergeCells.size() != possibleMergeFaces.size())
        FatalErrorInFunction<<"Invalid input parameters!"<<exit(FatalError);
        
    List<DynamicList<label>> assignList(possibleMergeCells.size());
    List<label> usedCellsMultiplicity(possibleMergeCells.size(),0);
    List<label> cellRequestedMultiplicity(possibleMergeCells.size(),0);
    for(int i=0;i<cellRequestedMultiplicity.size();i++)
    {
        if((mergeNecessary[i] && possibleMergeCells[i].size()==0) || (!mergeNecessary[i] && possibleMergeCells[i].size()!=0))
        {
            Info<<endl;
            Info<<"mergeNecessary["<<i<<"]:"<<mergeNecessary[i]<<endl;
            Info<<"possibleMergeCells["<<i<<"]:"<<possibleMergeCells[i]<<endl;
            FatalErrorInFunction<<"Invalid input parameters!"<<exit(FatalError);
        }

        for(int j=0;j<possibleMergeCells[i].size();j++)
        {
            if(possibleMergeCells[i][j].size()!=1)
                FatalErrorInFunction<<"Invalid input parameters!"<<exit(FatalError);
            
            cellRequestedMultiplicity[possibleMergeCells[i][j][0]]++;
        }
    }
    
    labelList faceUsedByCell(owner.size(),-1);
    for(int count=0;count<possibleMergeCells.size();count++)
    {   
        if(mergeNecessary[count])
        /* Decision A: Enters if block if merge is necessary and the cell is not already used for
         * a merge with another cell
         */
        {
            if(usedCellsMultiplicity[count]==0)
            {
                DynamicList<label> largeEnoughMergeCells;
                for(int j=0;j<possibleMergeCells[count].size();j++)
                {
                    if(possibleMergeCells[count][j].size()!=1)
                        FatalErrorInFunction<<"More than one merge cell in one option!"<<exit(FatalError);
                
                    if(possibleMergeCellsPartialSize[count][j]>=partialThreeshold)
                        largeEnoughMergeCells.append(j);
                }
                if(largeEnoughMergeCells.size()>0)
                {
                    DynamicList<label> nonTakenMergeCells;
                    for(int j=0;j<largeEnoughMergeCells.size();j++)
                    {
                        int ind = largeEnoughMergeCells[j];
                        if(usedCellsMultiplicity[possibleMergeCells[count][ind][0]]==0)
                            nonTakenMergeCells.append(ind);
                    }
                    if(nonTakenMergeCells.size()>0)
                    {
                        label minimumRequestedCell = possibleMergeCells.size()+1;
                        label minimumRequestedCellInd = -1;
                        for(int j=0;j<nonTakenMergeCells.size();j++)
                        {
                            if(minimumRequestedCell>cellRequestedMultiplicity[possibleMergeCells[count][nonTakenMergeCells[j]][0]])
                            {
                                minimumRequestedCell=cellRequestedMultiplicity[possibleMergeCells[count][nonTakenMergeCells[j]][0]];
                                minimumRequestedCellInd = nonTakenMergeCells[j];
                            }
                        }
                        if(minimumRequestedCellInd==-1 || minimumRequestedCell==possibleMergeCells.size()+1)
                            FatalErrorInFunction<<"Can not happen!"<<exit(FatalError);
                
                        assignList[count].append(possibleMergeFaces[count][minimumRequestedCellInd][0]);
                        usedCellsMultiplicity[possibleMergeCells[count][minimumRequestedCellInd][0]]++;
                    
                        if(faceUsedByCell[possibleMergeFaces[count][minimumRequestedCellInd][0]] != -1)
                        {
                            Info<<endl;
                            Info<<"cell "<<faceUsedByCell[possibleMergeFaces[count][minimumRequestedCellInd][0]];
                            label face = assignList[faceUsedByCell[possibleMergeFaces[count][minimumRequestedCellInd][0]]][0];
                            Info<<" to merge: "<<owner[face]<<" and "<<neighbour[face]<<" via "<<face<<endl;

                            Info<<"cell "<<count<<" merged with "<<possibleMergeCells[count][minimumRequestedCellInd][0]<<" via face: "<<possibleMergeFaces[count][minimumRequestedCellInd][0]<<endl;
                            Info<<"usedCellsMultiplicity["<<count<<"]:"<<usedCellsMultiplicity[count]<<endl;
                        
                            FatalErrorInFunction<<"Face used twice!"<<exit(FatalError);
                        }
                        faceUsedByCell[possibleMergeFaces[count][minimumRequestedCellInd][0]] = count;
                    }
                    else
                    {
                        label minimumUsedCell = possibleMergeCells.size()+1;
                        label minimumUsedCellInd = -1;
                        for(int j=0;j<largeEnoughMergeCells.size();j++)
                        {                
                            if(minimumUsedCell>usedCellsMultiplicity[possibleMergeCells[count][largeEnoughMergeCells[j]][0]])
                            {
                                minimumUsedCell=usedCellsMultiplicity[possibleMergeCells[count][largeEnoughMergeCells[j]][0]];
                                minimumUsedCellInd=largeEnoughMergeCells[j];
                            }
                        }
                        assignList[count].append(possibleMergeFaces[count][minimumUsedCellInd][0]);
                        usedCellsMultiplicity[possibleMergeCells[count][minimumUsedCellInd][0]]++;
                    
                        if(faceUsedByCell[possibleMergeFaces[count][minimumUsedCellInd][0]] != -1)
                        {
                            Info<<endl;
                            Info<<"cell "<<faceUsedByCell[possibleMergeFaces[count][minimumUsedCellInd][0]];
                            label face = assignList[faceUsedByCell[possibleMergeFaces[count][minimumUsedCellInd][0]]][0];
                            Info<<" to merge: "<<owner[face]<<" and "<<neighbour[face]<<" via "<<face<<endl;

                            Info<<"cell "<<count<<" merged with "<<possibleMergeCells[count][minimumUsedCellInd][0]<<" via face: "<<possibleMergeFaces[count][minimumUsedCellInd][0]<<endl;
                            FatalErrorInFunction<<"Face used twice!"<<exit(FatalError);
                        }                    
                        faceUsedByCell[possibleMergeFaces[count][minimumUsedCellInd][0]] = count;
                    }
                }
                else
                {
                    DynamicList<label> largerThanOneTenthPartialThreeshold;
                    for(int j=0;j<possibleMergeCells[count].size();j++)
                    {
                        if(possibleMergeCellsPartialSize[count][j]>=0.1*partialThreeshold)
                            largerThanOneTenthPartialThreeshold.append(j);
                    }
                    if(largerThanOneTenthPartialThreeshold.size()>0)
                    {
                        scalar largestCell = 0;
                        label largestCellInd = -1;
                        for(int j=0;j<largerThanOneTenthPartialThreeshold.size();j++)
                        {
                            if(largestCell<possibleMergeCellsPartialSize[count][largerThanOneTenthPartialThreeshold[j]])
                            {
                                largestCell=possibleMergeCellsPartialSize[count][largerThanOneTenthPartialThreeshold[j]];
                                largestCellInd = largerThanOneTenthPartialThreeshold[j];
                            }
                        }
                        if(largestCell==0 || largestCellInd==-1)
                            FatalErrorInFunction<<"Can not happen!"<<exit(FatalError);
                
                        assignList[count].append(possibleMergeFaces[count][largestCellInd][0]);
                        usedCellsMultiplicity[possibleMergeCells[count][largestCellInd][0]]++;
                    
                        if(faceUsedByCell[possibleMergeFaces[count][largestCellInd][0]] != -1)
                        {
                            Info<<endl;
                            Info<<"cell "<<faceUsedByCell[possibleMergeFaces[count][largestCellInd][0]];
                            label face = assignList[faceUsedByCell[possibleMergeFaces[count][largestCellInd][0]]][0];
                            Info<<" to merge: "<<owner[face]<<" and "<<neighbour[face]<<" via "<<face<<endl;

                            Info<<"cell "<<count<<" merged with "<<possibleMergeCells[count][largestCellInd][0]<<" via face: "<<possibleMergeFaces[count][largestCellInd][0]<<endl;
                            Info<<"usedCellsMultiplicity["<<count<<"]:"<<usedCellsMultiplicity[count]<<endl;
                        
                            FatalErrorInFunction<<"Face used twice!"<<exit(FatalError);
                        }
                        faceUsedByCell[possibleMergeFaces[count][largestCellInd][0]] = count;
                    }
                    else
                    {
                        label maximumRequestedCell = -1;
                        label maximumRequestedCellInd = -1;
                        for(int j=0;j<possibleMergeCells[count].size();j++)
                        {
                            if(maximumRequestedCell<cellRequestedMultiplicity[possibleMergeCells[count][j][0]])
                            {
                                maximumRequestedCell=cellRequestedMultiplicity[possibleMergeCells[count][j][0]];
                                maximumRequestedCellInd = j;
                            }
                        }
                        if(maximumRequestedCellInd==-1 || maximumRequestedCell==-1)
                            FatalErrorInFunction<<"Can not happen!"<<exit(FatalError);
                
                        assignList[count].append(possibleMergeFaces[count][maximumRequestedCellInd][0]);
                        usedCellsMultiplicity[possibleMergeCells[count][maximumRequestedCellInd][0]]++;
                    
                        if(faceUsedByCell[possibleMergeFaces[count][maximumRequestedCellInd][0]] != -1)
                        {
                            Info<<endl;
                            Info<<"cell "<<faceUsedByCell[possibleMergeFaces[count][maximumRequestedCellInd][0]];
                            label face = assignList[faceUsedByCell[possibleMergeFaces[count][maximumRequestedCellInd][0]]][0];
                            Info<<" to merge: "<<owner[face]<<" and "<<neighbour[face]<<" via "<<face<<endl;

                            Info<<"cell "<<count<<" merged with "<<possibleMergeCells[count][maximumRequestedCellInd][0]<<" via face: "<<possibleMergeFaces[count][maximumRequestedCellInd][0]<<endl;
                            Info<<"usedCellsMultiplicity["<<count<<"]:"<<usedCellsMultiplicity[count]<<endl;
                        
                            FatalErrorInFunction<<"Face used twice!"<<exit(FatalError);
                        }
                        faceUsedByCell[possibleMergeFaces[count][maximumRequestedCellInd][0]] = count;
                    }
                }
            }
            else
            {
                assignList[count].append(-2);
            }
        }
        else
        /* Decision A: Enters else block for cells that are not too small. The assignList is
         * filled with -1 for these cells.
         */
        {
            assignList[count].append(-1);
        }
    }
    
    Info<<"Fin"<<endl;
    
//Test
    for(int i=0;i<assignList.size();i++)
    {
        if(mergeNecessary[i])
        {
            if(assignList[i].size()==0)
            {
                FatalErrorInFunction
                << " Error empty assign List at i:"<<i<<endl
                << exit(FatalError);
            }
            if(assignList[i][0] == -1)
            {
                FatalErrorInFunction
                << " Too small cell was not treated by backtracking algorithm!"<<endl
                << exit(FatalError);
            }            
            
            std::unordered_multiset<label> mergeCellsSet;
            DynamicList<label> mergeCellsList;
            for(int k=0;k<assignList[i].size();k++)
            {
                if(assignList[i][k] >= neighbour.size())
                {
                    FatalErrorInFunction<<"Boundary face is merge face"<<exit(FatalError);
                }
                label ownerCell = owner[assignList[i][k]];
                label neighborCell = neighbour[assignList[i][k]];
                if(mergeCellsSet.find(ownerCell) == mergeCellsSet.end())
                    mergeCellsList.append(ownerCell);
                if(mergeCellsSet.find(neighborCell) == mergeCellsSet.end())
                    mergeCellsList.append(neighborCell);
                mergeCellsSet.insert(ownerCell);
                mergeCellsSet.insert(neighborCell);
            }
            if(mergeCellsList.size()<=0)
                FatalErrorInFunction<<"Merging option with no cells!"<<exit(FatalError);
            
            label mergeCellMult = mergeCellsSet.count(mergeCellsList[0]);
            for(int k=1;k<mergeCellsList.size();k++)
            {
                if(mergeCellsSet.count(mergeCellsList[k]) != static_cast<long unsigned int>(mergeCellMult))
                {
                    FatalErrorInFunction<<"Different multiplcity of cells in "<<
                    mergeCellsList.size()<<" cell merging "<<exit(FatalError);
                }
            }
            
            if(mergeCellMult==1 && mergeCellsList.size() == 2)
            {
                //2 Cell merging
                if(assignList[i].size() != 1)
                    FatalErrorInFunction<<"Error in 2 Cell merging! "<<exit(FatalError);
            }
            else
            {
                Info<<"mergeCellMult: "<<mergeCellMult<<endl;
                Info<<"optionMergeCellsList.size() == "<<mergeCellsList.size()<<endl;
                Info<<"possibleMergeFaces[i][j].size() == "<<assignList[i].size()<<endl;
                FatalErrorInFunction<<"Inconsistent merge Option!"<<exit(FatalError);
            }
        }
        else
        {
            if(assignList[i].size()==0)
            {
                FatalErrorInFunction
                << " Error empty assign List at i:"<<i<<endl
                << exit(FatalError);
            }
            if(assignList[i][0] != -1)
            {
                FatalErrorInFunction
                << " Backtracking algorithm wrote inside large enough cell. Something is wrong here!"<<endl
                << exit(FatalError);
            }
            assignList[i][0] = -1;
        }
    }
//Ende Test
    
    /*
    labelList assList(assignList.size());
    for(int i=0;i<assignList.size();i++)
    {
        if(assignList[i].size() != 1)
        {
            FatalErrorInFunction
            << " Error!"<<endl
            << exit(FatalError);
        }
        assList[i] = assignList[i][0];
    }
    return assList;
    */
    return assignList;
    /*
    Label index
    -1 : Cell is not too small
    -2 : cell is already merged by other cell
    -3 : No assignment
    */
}

void Foam::cutCellFvMesh::testNewMeshData
(
    const faceList& newFaces,
    const labelList& newFaceOwner,
    const labelList& newFaceNeighbor,
    const labelList& patchStarts,
    const labelList& patchSizes
)
{
    /*
    Info<<"Face number:"<<newFaces.size()
    <<" Owner number:"<<newFaceOwner.size()
    <<" Neighbor number:"<<newFaceNeighbor.size()<<endl;
    */
    
    if(newFaces.size() != newFaceOwner.size())
    {
        FatalErrorInFunction
        <<"Face number:"<<newFaces.size()
        <<" Owner number:"<<newFaceOwner.size()
        <<" Neighbor number:"<<newFaceNeighbor.size()
        << exit(FatalError);
    }
    bool boundaryFacesReached = false;
    for(int i=0;i<newFaceOwner.size();i++)
    {
        if(newFaceOwner[i] == -1)
        {       
            FatalErrorInFunction
            <<"Face "<<i<<" is owned by: -1 "
            <<"while neighbour is "<<newFaceNeighbor[i]
            << exit(FatalError);
        }
        if(i<newFaceNeighbor.size())
        {
            if(newFaceOwner[i] == newFaceNeighbor[i])
            {
                FatalErrorInFunction
                <<"Face "<<i<<" is owned by: "<<newFaceOwner[i]
                <<"while neighbour is "<<newFaceNeighbor[i]
                << exit(FatalError);
            }
            
            if(boundaryFacesReached && newFaceNeighbor[i] != -1)
            {
                if(i-1>0 && i+1<newFaceOwner.size())
                {
                    Info<<"Neighbour["<<i-1<<"]:"<<newFaceNeighbor[i-1]
                    <<" Neighbour["<<i<<"]:"<<newFaceNeighbor[i]
                    <<" Neighbour["<<i+1<<"]:"<<newFaceNeighbor[i+1]<<endl;
                }            
                FatalErrorInFunction
                <<"Boundary face before inner face."
                << exit(FatalError);
            }
            else if(!boundaryFacesReached && newFaceNeighbor[i] == -1)
            {
                boundaryFacesReached = true;
            }
        }
    }    
}

void Foam::cutCellFvMesh::testForCellSize
(
    DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
    DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
    DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
    DynamicList<bool>& oneMergeFaceSufficient,
    DynamicList<bool>& mergeNecessary,
    labelList& mergeFaceOfCell,
    scalar partialThreeshold
)
{
    const cellList& newCells = this->cells();
    const faceList& faces = this->faces();
    const labelList& owner   = this->faceOwner();
    const labelList& neighbour = this->faceNeighbour();
    const pointField& points = this->points();
    
    scalar factor = (1+partialThreeshold)/partialThreeshold;
    const cellList& cell = this->cells();
    const faceList& face = this->faces();
    const pointField& point = this->points();
    
    scalar minCellVol = cell[0].mag(point,face);
    label minCellInd = 0;
    scalar maxCellVol = cell[0].mag(point,face);
    label maxCellInd = 0;
    scalar CellVolAvg = 0;
    
    scalar vol;
    for(int i=0;i<cell.size();i++)
    {
        vol = cell[i].mag(point,face);
        CellVolAvg += vol;
        if(vol > maxCellVol)
        {
            maxCellVol = vol;
            maxCellInd = i;
        }
        if(vol < minCellVol)
        {
            minCellVol = vol;
            minCellInd = i;
        }
    }
    CellVolAvg /= cell.size();
    
    if((minCellVol*factor) < maxCellVol)
    {
        Info<<endl<<"Minimum cell "<<minCellInd<<" vol:"<<minCellVol
        <<" is more than four times smaller than Maximum "
        <<"cell "<<maxCellInd<<" vol:"<<maxCellVol<<endl;
        Info<<" Average vol was:"<<CellVolAvg<<endl;        
        
        label neighbourCell=-1;
        scalar neighbourCellVolume;
        scalar smallCellVolume = newCells[minCellInd].mag(points,faces);
        
        Info<<endl<<"Problem in cell "<<minCellInd<<" with volume "<<smallCellVolume<<endl;
        Info<<"\tNeighbours are:"<<endl;
        
        for(int k=0;k<newCells[minCellInd].size();k++)
        {
            if(newCells[minCellInd][k] < neighbour.size())
            {
                if(owner[newCells[minCellInd][k]] == minCellInd)
                {
                    neighbourCell = neighbour[newCells[minCellInd][k]];
                }
                else if(neighbour[newCells[minCellInd][k]] == minCellInd)
                {   
                    neighbourCell = owner[newCells[minCellInd][k]];
                }
                else
                    FatalErrorInFunction<<"Found no neighbour cell in testForCellSize"<< exit(FatalError);
                
                neighbourCellVolume = newCells[neighbourCell].mag(points,faces);
                
                Info<<"\tCell:"<<neighbourCell<<" partialVol:"<<neighbourCellVolume<<
                    " combinedPartialVol:"<<(neighbourCellVolume + smallCellVolume)<<endl;
            }
        }
        /*
        scalar vol;

        Info<<"The following cells are too small:"<<endl;
        for(int i=0;i<cell.size();i++)
        {
            vol = cell[i].mag(point,face);
            if(4*vol < maxCellVol)
            {
                Info<<"\tCell "<<i<<" with vol:"<<vol<<endl;
            }
        }
        
        Info<<"The following cells are too large:"<<endl;
        for(int i=0;i<cell.size();i++)
        {
            vol = cell[i].mag(point,face);
            if(0.25*vol > minCellVol)
            {
                Info<<"\tCell "<<i<<" with vol:"<<vol<<endl;
            }
        }
        */
        FatalErrorInFunction
        << "Cell size problem"
        << exit(FatalError);  
    }
}

void Foam::cutCellFvMesh::correctFaceNormalDir
(
    const pointField& points,
    faceList& faces,
    const labelList& owner,
    const labelList& neighbour
)
{
    scalar maxFaceSize = -1,avgFaceSize=0,minFaceSize;
    if(faces.size()>0) minFaceSize = faces[0].mag(points);
    for(int i=0;i<faces.size();i++)
    {
        scalar thisFaceSize = faces[i].mag(points);
        avgFaceSize+=thisFaceSize;
        if(maxFaceSize<thisFaceSize)
            maxFaceSize=thisFaceSize;
        if(minFaceSize>thisFaceSize)
            minFaceSize=thisFaceSize;
    }
    avgFaceSize /= faces.size();
    //Info<<"owner.size()="<<owner.size()<<endl;
    //Info<<"neighbour.size()="<<neighbour.size()<<endl;

    //Info<<"Correct1"<<endl;
    label numCells = 0;
    for(int i=0;i<owner.size();i++)
    {
        numCells = std::max(numCells,owner[i]);
    }
    numCells++;
    List<DynamicList<label>> cellsFaces(numCells);

    //Info<<"Correct2::"<<cellsFaces.size()<<endl;

    //Info<<"Correct2"<<endl;
    
    if(owner.size() != neighbour.size())
        FatalErrorInFunction<<"Can not happen"<< exit(FatalError); 

    /*
    for(int i=0;i<owner.size();i++)
        Info<<"i:"<<i<<" owner["<<i<<"]:"<<owner[i]<<" neighbour["<<i<<"]:"<<neighbour[i]<<" / "<<neighbour.size()<<" /-/ "<<cellsFaces.size()<<" numCells: "<<numCells<<endl;
    */
    
    //FatalErrorInFunction<<"Temporary stop"<< exit(FatalError); 
    
    for(int i=0;i<owner.size();i++)
    {
        //Info<<"i:"<<i<<" owner[i]:"<<owner[i]<<" neighbour[i]:"<<neighbour[i]<<" / "<<neighbour.size()<<" /-/ "<<cellsFaces.size()<<" numCells: "<<numCells<<endl;
        cellsFaces[owner[i]].append(i);
        if(neighbour[i]!=-1)
            cellsFaces[neighbour[i]].append(i);
    }

    DynamicList<cell> cellList;
    cellList.setSize(numCells);
    
    Info<<"Correct3"<<endl;
    
    for(int i=0;i<numCells;i++)
    {
        labelList cellFaces = cellsFaces[i];
        cell oneCell(cellFaces);
        cellList[i] = oneCell;
    }
    
    Info<<"Correct4"<<endl;

    for(int i=0;i<faces.size();i++)
    {
        //Info<<"Test face "<<i<<endl;
        point centreFace = faces[i].centre(points);
        vector normalFace = faces[i].normal(points);
    
        vector faceCentreToOwnerCentre = cellList[owner[i]].centre(points,faces)-centreFace;
        if((faceCentreToOwnerCentre && normalFace)>=0)
        {
            faces[i].flip();
            /*
            DynamicList<label> reversedFaceList;
            for(int k=faces[i].size()-1;k>=0;k++)
            {
                reversedFaceList.append(faces[i][k]);
            }
            face newFace(reversedFaceList);
            face[i] = newFace;
            */
        }
    }
    Info<<"Correct5"<<endl;
    
    for(int i=0;i<faces.size();i++)
    {
        //Info<<"Test face "<<i<<endl;
        point centreFace = faces[i].centre(points);
        vector normalFace = faces[i].normal(points);
        scalar area = faces[i].mag(points);
    
        vector faceCentreToOwnerCentre = cellList[owner[i]].centre(points,faces)-centreFace;
        if((faceCentreToOwnerCentre && normalFace)>=0)
        {
            if((area >= maxFaceSize*partialThreeshold*(1.f/1e10)) && norm2(normalFace)!=0 && norm2(faceCentreToOwnerCentre)!=0)
            {
                Info<<endl<<endl;
                Info<<"Face:"<<i<<" Owner:"<<owner[i]<<" "<<endl;;
                if(i < neighbour.size())
                    Info<<"Neighbor:"<<neighbour[i]<<" ";
                Info<<endl;
                for(int k=0;k<faces[i].size();k++)
                {
                    Info<<faces[i][k]<<points[faces[i][k]]<<"->";
                }
                Info<<endl<<" with face centre:"<<centreFace;
                Info<<endl<<" and face normal vector:"<<normalFace;
                Info<<endl<<" and area:"<<area<<endl;
                Info<<"nbrOfPrevPoints: "<<nbrOfPrevPoints<<endl;
                Info<<"nbrOfPrevFaces: "<<nbrOfPrevFaces<<endl;
                
                Info<<"------------------------------------"<<endl;
            
                cell oneCell = cellList[owner[i]];
                Info<<"Owner Cell centre is: "<<oneCell.centre(points,faces)<<endl;
                Info<<"Owner Cell: "<<oneCell<<endl;
                Info<<"Owner Cell Size: "<<oneCell.size()<<endl;
                Info<<"Owner Cell volume: "<<oneCell.mag(points,faces)<<endl;
            
                for(int k=0;k<oneCell.size();k++)
                {
                    label oneFaceInd = oneCell[k];
                    face oneFace = faces[oneFaceInd];
                    Info<<"Face "<<k<<" area: "<<oneFace.mag(points)<<"::";
                    for(int kk=0;kk<oneFace.size();kk++)
                    {
                        Info<<oneFace[kk]<<points[oneFace[kk]]<<"->";
                    }
                    Info<<endl;
                }
                
                Info<<"------------------------------------"<<endl;
            
                if(i < neighbour.size())
                {
                    cell oneCell = cellList[neighbour[i]];
                    Info<<"Neighbour Cell centre is: "<<oneCell.centre(points,faces)<<endl;
                    Info<<"Neighbour Cell: "<<oneCell<<endl;
                    Info<<"Neighbour Cell Size: "<<oneCell.size()<<endl;
                    Info<<"Neighbour Cell volume: "<<oneCell.mag(points,faces)<<endl;
            
                    for(int k=0;k<oneCell.size();k++)
                    {
                    label oneFaceInd = oneCell[k];
                    face oneFace = faces[oneFaceInd];
                    Info<<"Face "<<k<<" area: "<<oneFace.mag(points)<<"::";
                    for(int kk=0;kk<oneFace.size();kk++)
                    {
                        Info<<oneFace[kk]<<points[oneFace[kk]]<<"->";
                    }
                    Info<<endl;
                    }
                }
            
                FatalErrorInFunction
                <<"Normal vector is "<<normalFace<<" while faceCentreToOwnerCentre is "<<faceCentreToOwnerCentre<<"!"
                <<" They must have a opposite direction"
                << exit(FatalError);
            }
            else
            {
                for(int j=0;j<faces[i].size();j++)
                {
                    for(int k=0;k<faces[i].size();k++)
                    {
                        if(j!=k && faces[i][j]==faces[i][k])
                        {
                            Info<<"faces["<<i<<"]: "<<faces[i]<<endl;
                            FatalErrorInFunction<<"Face with zero area has duplicate points"<<exit(FatalError);
                        }
                    }
                }
            }
        }
    }
}
