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

    //printMesh();
    
    Info<<"a";
    cellList cutCells(0);
    for(int i=0;i<cutCellsMinusAndPlus.size();i++)
    {
        cell cutCell(Foam::clone(cutCellsMinusAndPlus[i]));
        cutCells.append(cutCell);
    }
    Info<<cutCells.size();

    const cellList& cellsLis = this->cells();
    const faceList& faceList = this->faces();
    const pointField& pointLis = this->points();
    
    Info<<"c";

    oldCellVolume = scalarList(cellsLis.size());
    for(int i=0;i<oldCellVolume.size();i++)
    {
        oldCellVolume[i] = cellsLis[i].mag(pointLis,faceList);
        Info<<"Old Cell "<<i<<": "<<oldCellVolume[i]<<endl;
    }
    
    
    scalar wholeCellVol = 0;
    scalar plusCutCellVol = 0;
    scalar minusCutCellVol = 0;
    for(int i=0;i<cellsLis.size();i++)
    {
        wholeCellVol = oldCellVolume[i];
        Info<<"Cell "<<i<<" Vol: "<<wholeCellVol<<endl;
        if((oldCellToMinusCutCell[i] != -1 && oldCellToPlusCutCell[i] == -1) ||
           (oldCellToMinusCutCell[i] == -1 && oldCellToPlusCutCell[i] != -1))
            Info<<"Can not happen"<<endl;
        if(oldCellToPlusCutCell[i] != -1)
        {
            Info<<oldCellToMinusCutCell[i]<<endl;
            Info<<oldCellToPlusCutCell[i]<<endl;
            
            cell minusCell = cutCells[oldCellToMinusCutCell[i]];
            Info<<"\t MinusCell: ";
            for(int i=0;i<minusCell.size();i++)
                Info<<minusCell[i]<<" ";
            Info<<endl;
            
            cell plusCell = cutCells[oldCellToPlusCutCell[i]];
            Info<<"\t PlusCell: ";
            for(int i=0;i<plusCell.size();i++)
                Info<<plusCell[i]<<" ";
            Info<<endl;
            
            Info<<"Got cells"<<endl;
            minusCutCellVol = minusCell.mag(newMeshPoints_,faces);
            plusCutCellVol = plusCell.mag(newMeshPoints_,faces);
            Info<<"\t MinusCell: "<<minusCutCellVol<<endl;
            Info<<"\t PlusCell: "<<plusCutCellVol<<endl;
        }
    }
    
    scalarList cutCellVol(cutCells.size());
    for(int i=0;i<cutCells.size();i++)
    {
        cutCellVol[i] = cutCells[i].mag(newMeshPoints_,faces);
    }
    
    
    /*
    resetPrimitives(Foam::clone(newMeshPoints_),
                    Foam::clone(faces),
                    Foam::clone(owner),
                    Foam::clone(neighbour),
                    patchSizes,
                    patchStarts,
                    true);
    */
    
    printMesh();
    
    scalarList solidFrac(oldCellVolume.size());
    for(int i=0;i<oldCellVolume.size();i++)
    {
        Info<<i<<endl;
        if(cellsToSide_[i] == 1)
            solidFrac[i] = 0;
        else if(cellsToSide_[i] == -1)
            solidFrac[i] = 1;
        else
        {
            scalar plusSideVol = cutCellVol[oldCellToPlusCutCell[i]];
            Info<<plusSideVol<<endl;            
            scalar minusSideVol = cutCellVol[oldCellToMinusCutCell[i]];
            Info<<minusSideVol<<endl;
            solidFrac[i] = minusSideVol / (plusSideVol+minusSideVol);
        }
        
        Info<<"Cell "<<i<<" Partial: "<<solidFrac[i]<<endl;
    }
    
    //this->clearGeom();
    //this->clearOut();
    //this->clearPrimitives();
    
    //printMesh();
    //this->write();
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
    //printMesh();
    
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
    Info<<"Please write"<<endl;
    this->write();
    printMesh();
    //selfTestMesh();
}
 
void Foam::cutCellPolyMesh::projectNurbsSurface()
{
    const pointField& points = this->points();
    pointDist = scalarList(points.size());
    
    for(int i=0;i<points.size();i++)
    {
        Info<<"Working on Point: "<<i<<" "<<points[i]<<endl;
        
        std::unique_ptr<labelList> firstOrderNearNurbs = MainTree->nearNurbsCurves(points[i]);
        if(firstOrderNearNurbs->size() == 0)
        {
            pointDist[i] = 1;
            Info<<"\tSkipped because far away"<<endl;
            continue;
        }
        Info<<"\tGot list size:"<<firstOrderNearNurbs->size()<<endl;
        
        scalarList distToNurbsSurface(0);
        scalarList paraToNurbsSurface(0);
        bool outSideNurbsBox = false;
        for(int k=0;k<firstOrderNearNurbs->size();k++)
        {
            label thisNurbs = (*firstOrderNearNurbs)[k];
            scalar thisNodePara = NurbsTrees[thisNurbs]->closestParaOnNurbsToPoint(points[i]);
            Info<<"\tIndex of nurbs:"<<thisNurbs<<" with para: "<<thisNodePara<<endl;
            if(thisNodePara < this->Curves[thisNurbs]->min_U())
            {
                pointDist[i] = 1;
                outSideNurbsBox = true;
                continue;
            }
            paraToNurbsSurface.append(thisNodePara);
            distToNurbsSurface.append(this->Curves[thisNurbs]->distanceToNurbsSurface(thisNodePara,points[i]));
        }
        if(outSideNurbsBox)
            continue;
        
        scalar minDistToNurbsSurface = std::numeric_limits<scalar>::max();
        for(int k=0;k<distToNurbsSurface.size();k++)
        {
            if(distToNurbsSurface[k] < minDistToNurbsSurface)
                minDistToNurbsSurface = distToNurbsSurface[k];
        }
        pointDist[i] = minDistToNurbsSurface;
        
        //Info<<"Finished working on Point: "<<points[i]<<" "<<pointDist[i]<<endl;
    }
    
    Info<<"Final point writing"<<endl;
    for(int i=0;i<points.size();i++)
    {
        Info<<i<<" "<<points[i]<<" "<<pointDist[i]<<endl;
    }
}
