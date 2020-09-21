#include "cutCellPolyMesh.H"

Foam::cutCellPolyMesh::cutCellPolyMesh
(
    const IOobject& io,
    List<std::shared_ptr<Nurbs>> Curves,
    cutStatus state
):
polyMesh(io),
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
    
    /*
    if(state == delNegMesh)
    {
        const polyBoundaryMesh& boundMesh = this->boundaryMesh();
        for(int i=0;i<boundMesh.size();i++)
        {
            polyPatch thisPatch = boundMesh[i];
            Info<<"-------Patch "<<thisPatch.index()<<"----------"<<endl;
            Info<<"| Name:"<<thisPatch.name()<<endl;
            Info<<"| Size:"<<thisPatch.size()<<endl;
            Info<<"| Start:"<<thisPatch.start()<<endl;
            Info<<"| patchI:"<<thisPatch.index()<<endl;
            Info<<"| patchType:"<<thisPatch.physicalType()<<endl;
        }

        patchStarts = labelList(boundMesh.size());
        patchSizes = labelList(boundMesh.size());
        for(int i=0;i<boundMesh.size();i++)
        {
            patchStarts[i] = boundMesh[i].start();
            patchSizes[i] = boundMesh[i].faceCentres().size();
        }
        patchStarts.append(patchStarts.last()+patchSizes.last());
        
        for(int i=0;i<boundMesh.size();i++)
        {        
            Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
        }
        
        Info<<boundMesh.size()<<endl;

        label insertPatchi = boundMesh.size()+1;
        
        word name = "cutWall";
        word patchType = "wall";

        polyPatch* patch= new polyPatch
        (
            name,
            addedCutFaces.size(),
            patchStarts.last(),
            insertPatchi,
            boundMesh,
            patchType
        );
        List<polyPatch*> p(1,patch);
        
        fileName dictName = "wallDict";
        dictionary patchFieldDict(dictName);
        word defaultPatchFieldType = "wall";
        
        Info<<"-------Patch "<<patch->index()<<"----------"<<endl;
        Info<<"| Name:"<<patch->name()<<endl;
        Info<<"| Size:"<<patch->size()<<endl;
        Info<<"| Start:"<<patch->start()<<endl;
        Info<<"| patchI:"<<patch->index()<<endl;
        Info<<"| patchType:"<<patch->physicalType()<<endl;
        
        boundMesh.append(patch);
        Info<<boundary_.size()<<endl;
        
        
        //this->addPatches(p);
        
        
        this->addPatch(insertPatchi,patch,patchFieldDict,
                       defaultPatchFieldType,true);  
        
    }
    */
    
    
    this->write();
    //printMesh();
    selfTestMesh();
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
