#include "cutCellPolyMesh.H"

Foam::cutCellPolyMesh::cutCellPolyMesh
(
    const IOobject& io,
    List<std::shared_ptr<Nurbs>> Curves,
    Time& runTime,
    std::unique_ptr<volScalarField>& solidFraction
):
fvMesh(io),
Curves(std::move(Curves)),
MainTree(new KdTree(this->Curves)),
NurbsTrees(List<std::unique_ptr<BsTree>>(this->Curves.size()))
{
    for(int i=0;i<this->Curves.size();i++)
    {
        NurbsTrees[i] = std::move(std::unique_ptr<BsTree>(new BsTree(this->Curves[i])));
    }
    Info<<"Prepared all Data"<<endl;
    
    projectNurbsSurface();
    Info<<"Projected Nurbs Surface"<<endl;
    
    newMeshPoints();
    Info<<"Added Mesh Points"<<endl;
    
    printAddedPoints();
    newMeshEdges();
    edgesToSide();
    newMeshFaces();
    cutOldFaces();
    
    faceList faces(0);
    labelList owner(0);
    labelList neighbour(0);
        
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

    printMesh();
    
    cellList cutCells(0);
    for(int i=0;i<cutCellsMinusAndPlus.size();i++)
    {
        cell cutCell(Foam::clone(cutCellsMinusAndPlus[i]));
        cutCells.append(cutCell);
    }
    
    pointField oldPoints(Foam::clone(this->points()));
    faceList oldFaces(Foam::clone(this->faces()));
    labelList oldOwners(Foam::clone(this->faceOwner()));
    labelList oldNeighbours(Foam::clone(this->faceNeighbour()));
    const cellList& oldCells = this->cells();
    const polyBoundaryMesh& oldboundMesh = this->boundaryMesh();
    labelList oldPatchStarts(oldboundMesh.size());
    labelList oldPatchSizes(oldboundMesh.size());
    for(int i=0;i<oldboundMesh.size();i++)
    {
        oldPatchStarts[i] = oldboundMesh[i].start();
        oldPatchSizes[i] = oldboundMesh[i].faceCentres().size();
    }
    
    oldCellVolume = scalarList(oldCells.size());
    for(int i=0;i<oldCellVolume.size();i++)
    {
        oldCellVolume[i] = oldCells[i].mag(oldPoints,oldFaces);
        Info<<"Old Cell "<<i<<": "<<oldCellVolume[i]<<endl;
    }
    scalar wholeCellVol = 0;
    scalar plusCutCellVol = 0;
    scalar minusCutCellVol = 0;
    for(int i=0;i<oldCells.size();i++)
    {
        wholeCellVol = oldCellVolume[i];
        Info<<"Cell "<<i<<" Vol: "<<wholeCellVol<<endl;
        if((oldCellToMinusCutCell[i] != -1 && oldCellToPlusCutCell[i] == -1) ||
           (oldCellToMinusCutCell[i] == -1 && oldCellToPlusCutCell[i] != -1))
            Info<<"Can not happeen"<<endl;
        if(oldCellToPlusCutCell[i] != -1)
        {
            cell minusCell = cutCells[oldCellToMinusCutCell[i]];
            cell plusCell = cutCells[oldCellToPlusCutCell[i]];
            minusCutCellVol = minusCell.mag(oldPoints,faces);
            plusCutCellVol = plusCell.mag(oldPoints,faces);
            Info<<"\t MinusCell: "<<minusCutCellVol<<endl;
            Info<<"\t PlusCell: "<<plusCutCellVol<<endl;
        }
    }
    
      
    resetPrimitives(Foam::clone(newMeshPoints_),
                    Foam::clone(faces),
                    Foam::clone(owner),
                    Foam::clone(neighbour),
                    patchSizes,
                    patchStarts,
                    true);
    
    const cellList& newCells = this->cells();
    Info<<"Size: "<<newCells.size()<<endl; 
    newCellVolume = scalarList(newCells.size());
    for(int i=0;i<newCellVolume.size();i++)
    {
        newCellVolume[i] = newCells[i].mag(this->points(),this->faces());
        Info<<"Cell "<<i<<": "<<newCellVolume[i]<<endl;
    }
    
    scalarList solidFrac(oldCellVolume.size());
    for(int i=0;i<oldCellVolume.size();i++)
    {
        if(cellsToSide_[i] == 1)
            solidFrac[i] = 0;
        else if(cellsToSide_[i] == -1)
            solidFrac[i] = 1;
        else
            solidFrac[i] = newCellVolume[oldSplittedCellToNewMinusCell[i]] / oldCellVolume[i];
        
        Info<<"Cell "<<i<<" Partial: "<<solidFrac[i]<<endl;
    }
    
    resetPrimitives(Foam::clone(oldPoints),
                    Foam::clone(oldFaces),
                    Foam::clone(oldOwners),
                    Foam::clone(oldNeighbours),
                    oldPatchSizes,
                    oldPatchStarts,
                    true);
    
    //this->clearGeom();
    //this->clearOut();
    //this->clearPrimitives();
    
    const cellList& newCell2 = this->cells();
    Info<<"Size: "<<newCell2.size()<<endl; 
    
    //printMesh();
    this->write();
    //printMesh();
    //selfTestMesh();
    
    solidFraction = std::unique_ptr<volScalarField>
    (
        new volScalarField
        (
            Foam::IOobject
            (
                "solidFraction",
                runTime.timeName(),
                (*this),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            (*this)
        )
    );
    
    for(int i=0;i<oldCellVolume.size();i++)
    {
        (*solidFraction)[i] = solidFrac[i];
    }
    solidFraction->write();
}

Foam::cutCellPolyMesh::cutCellPolyMesh
(
    const IOobject& io,
    List<std::shared_ptr<Nurbs>> Curves,
    cutStatus state
):
fvMesh(io),
Curves(std::move(Curves)),
MainTree(new KdTree(this->Curves)),
NurbsTrees(List<std::unique_ptr<BsTree>>(this->Curves.size()))
{
    for(int i=0;i<this->Curves.size();i++)
    {
        NurbsTrees[i] = std::move(std::unique_ptr<BsTree>(new BsTree(this->Curves[i])));
    }
    Info<<"Prepared all Data"<<endl;
    
    projectNurbsSurface();
    Info<<"Projected Nurbs Surface"<<endl;
    
    newMeshPoints();
    Info<<"Added Mesh Points"<<endl;
    
    printAddedPoints();
    newMeshEdges();
    edgesToSide();
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
    for(int i=0;i<oldCellVolume.size();i++)
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
    
    const cellList& newCells = this->cells();
    newCellVolume = scalarList(newCells.size());
    for(int i=0;i<newCellVolume.size();i++)
    {
        newCellVolume[i] = newCells[i].mag(newMeshPoints_,this->faces());
    }
    
    agglomerateSmallCells_cutNeg(newCellVolume,oldCellVolume);

    //printMesh();
    this->write();
    //printMesh();
    //selfTestMesh();
}
 
void Foam::cutCellPolyMesh::projectNurbsSurface()
{
    const pointField& points = this->points();
    pointDist = scalarList(points.size());
    
    for(int i=0;i<points.size();i++)
    {
        //Info<<"Worked on Point: "<<points[i]<<endl;
        
        std::unique_ptr<labelList> firstOrderNearNurbs = MainTree->nearNurbsCurves(points[i]);
        if(firstOrderNearNurbs->size() == 0)
            continue;
        //Info<<"Got list size:"<<firstOrderNearNurbs->size()<<endl;
        
        scalarList distToNurbsSurface(0);
        scalarList paraToNurbsSurface(0);
        for(int k=0;k<firstOrderNearNurbs->size();k++)
        {
            label thisNurbs = (*firstOrderNearNurbs)[k];
            //Info<<"Index of nurbs:"<<thisNurbs<<endl;
            scalar thisNodePara = NurbsTrees[thisNurbs]->closestParaOnNurbsToPoint(points[i]);
            if(thisNodePara < this->Curves[thisNurbs]->min_U())
            {
                pointDist[i] = 1;
                continue;
            }
            paraToNurbsSurface.append(thisNodePara);
            distToNurbsSurface.append(this->Curves[thisNurbs]->distanceToNurbsSurface(thisNodePara,points[i]));
        }
        
        scalar minDistToNurbsSurface = std::numeric_limits<scalar>::max();
        for(int k=0;k<distToNurbsSurface.size();k++)
        {
            if(distToNurbsSurface[k] < minDistToNurbsSurface)
                minDistToNurbsSurface = distToNurbsSurface[k];
        }
        pointDist[i] = minDistToNurbsSurface;
        
        //Info<<"Finished working on Point: "<<points[i]<<" "<<pointDist[i]<<endl;
    }
    
    //Info<<"Final point writing"<<endl;
    for(int i=0;i<points.size();i++)
    {
        //Info<<points[i]<<" "<<pointDist[i]<<endl;
    }
}
