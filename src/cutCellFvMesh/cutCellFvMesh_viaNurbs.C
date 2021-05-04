#include "cutCellFvMesh.H"

Foam::cutCellFvMesh::cutCellFvMesh
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

Foam::cutCellFvMesh::cutCellFvMesh
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
    /*
    const objectRegistry& reg = this->thisDb();
    fileName path = reg.time().timePath();
    Info<<"path():"<<this->fvMesh::polyMesh::objectRegistry::path()<<endl;
    Info<<"rootPath():"<<this->fvMesh::polyMesh::objectRegistry::rootPath()<<endl;
    Info<<"caseName():"<<this->fvMesh::polyMesh::objectRegistry::caseName()<<endl;
    Info<<"instance():"<<this->fvMesh::polyMesh::objectRegistry::instance()<<endl;
    Info<<"local():"<<this->fvMesh::polyMesh::objectRegistry::local()<<endl;
    
    const objectRegistry& regis = this->fvMesh::polyMesh::objectRegistry::db();
    Info<<"reg Name:"<<regis.names()<<endl;


    if(path[path.size()-1]=='0' && path[path.size()-2]=='/')
    {
        Info<<"Is old"<<endl;
        Info<<path<<endl;
    }
    else
    {
        Info<<"Is new"<<endl;
        Info<<path<<endl;
    }
    FatalErrorInFunction<< "Temp stop"<<endl<< exit(FatalError);
    */
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span;
    
    Info<<"Created Main Tree"<<endl;
    for(int i=0;i<this->Curves.size();i++)
    {
        NurbsTrees[i] = std::move(std::unique_ptr<BsTree>(new BsTree(this->Curves[i])));
    }
    Info<<"Prepared all Data"<<endl;
    
    Info<<"Projection of distance field ";
    t1 = std::chrono::high_resolution_clock::now();
    projectNurbsSurface();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Adding of cut points";
    t1 = std::chrono::high_resolution_clock::now();
    newMeshPoints();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Adding of cut edges";
    t1 = std::chrono::high_resolution_clock::now();
    newMeshEdges();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t\t" << time_span.count() << " seconds."<<endl;
        
    edgesToSide();
    
    Info<<"Adding of cut faces";
    t1 = std::chrono::high_resolution_clock::now();
    newMeshFaces_plus();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t\t"<< time_span.count() << " seconds."<<endl;
        
    Info<<"Cutting old faces";
    t1 = std::chrono::high_resolution_clock::now();
    cutOldFaces_plus();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t\t\t" << time_span.count() << " seconds."<<endl;
        
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
        Info<<"-------------------------------------------"<<endl;
        Info<<"Create new Mesh data and cut negative cells";
        t1 = std::chrono::high_resolution_clock::now();
        createNewMeshData_cutNeg_plus();
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        Info<< " took \t"<< time_span.count() << " seconds."<<endl;
        Info<<"-------------------------------------------"<<endl;

        Info<<"Combine resetPrimitives data";
        t1 = std::chrono::high_resolution_clock::now();
        faces.append(splitAndUnsplitFacesInterior);
        faces.append(splitAndUnsplitFacesBoundary);
        faces.append(addedCutFaces);
        faces.append(splitAndUnsplitFacesInteriorToBoundary);
    
        owner.append(splitAndUnsplitFaceInteriorOwner);
        owner.append(splitAndUnsplitFaceBoundaryOwner);
        owner.append(addedCutFaceOwner);
        owner.append(splitAndUnsplitFaceInteriorToBoundaryOwner);
        
        neighbour.append(splitAndUnsplitFaceInteriorNeighbor);
        neighbour.append(splitAndUnsplitFaceBoundaryNeighbor);
        neighbour.append(addedCutFaceNeighbor);
        neighbour.append(splitAndUnsplitFaceInteriorToBoundaryNeighbor);
        
        Info<<"3"<<endl;
        for(int i=0;i<patchStarts.size();i++)
        {
            Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
        }
        
        patchStarts[patchStarts.size()-1] = patchStarts.last()+patchSizes.last();
        patchSizes[patchSizes.size()-1] = addedCutFaces.size()+splitAndUnsplitFacesInteriorToBoundary.size();
        
        Info<<endl;
        //Info<<"faces: "<<faces<<endl;
        //Info<<"owner: "<<owner<<endl;
        //Info<<"neighbour: "<<neighbour<<endl;
        
        
        Info<<"--"<<endl;
        for(int i=0;i<patchStarts.size();i++)
        {
            Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
        }
        Info<<"InteriorFacesSize: "<<splitAndUnsplitFacesInterior.size()<<endl;
        Info<<"BoundaryFaceSize: "<<splitAndUnsplitFacesBoundary.size()+addedCutFaces.size()+splitAndUnsplitFacesInteriorToBoundary.size()<<endl;

        Info<<"splitAndUnsplitFacesInterior: "<<splitAndUnsplitFacesInterior.size()<<endl;
        Info<<"splitAndUnsplitFacesBoundary: "<<splitAndUnsplitFacesBoundary.size()<<endl;
        Info<<"addedCutFaces: "<<addedCutFaces.size()<<endl;
        Info<<"splitAndUnsplitFacesInteriorToBoundary: "<<splitAndUnsplitFacesInteriorToBoundary.size()<<endl;

        
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        Info<< " took \t\t\t"<< time_span.count() << " seconds."<<endl;
    }
    //printMesh();
    //FatalErrorInFunction<<"Temporary stop!"<<exit(FatalError);

    correctFaceNormalDir(newMeshPoints_,faces,owner,neighbour);

    const pointField& oldPoints = this->points();
    const faceList& oldFaceList = this->faces();
    const cellList& oldCells = this->cells();
    
    oldCellVolume = scalarList(oldCells.size());
    for(int i=0;i<oldCellVolume.size();i++)
    {
        oldCellVolume[i] = oldCells[i].mag(oldPoints,oldFaceList);
    }
    
    testNewMeshData(faces,owner,neighbour,patchStarts,patchSizes);
        
    resetPrimitives(Foam::clone(newMeshPoints_),
                    Foam::clone(faces),
                    Foam::clone(owner),
                    Foam::clone(neighbour),
                    patchSizes,
                    patchStarts,
                    true);
    selfTestMesh();
        
    const cellList& newCells = this->cells();
    newCellVolume = scalarList(newCells.size());
    for(int i=0;i<newCellVolume.size();i++)
    {
        newCellVolume[i] = newCells[i].mag(newMeshPoints_,this->faces());
    }
    
    Info<<"Agglomerate small cut-cells";
    t1 = std::chrono::high_resolution_clock::now();
    agglomerateSmallCells_cutNeg_plus(newCellVolume,oldCellVolume);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t" << time_span.count() << " seconds."<<endl;

    //printMesh();
    Info<<"Please write"<<endl;
    this->write();
    Info<<"Written"<<endl;
    //printMesh();
    selfTestMesh();
}
 
void Foam::cutCellFvMesh::projectNurbsSurface()
{
    const pointField& points = this->points();
    pointDist = scalarList(points.size());
    
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::high_resolution_clock::time_point t3;
    std::chrono::high_resolution_clock::time_point t4;
    std::chrono::high_resolution_clock::time_point t5;
    std::chrono::high_resolution_clock::time_point t6;
    std::chrono::high_resolution_clock::time_point t7;
    std::chrono::high_resolution_clock::time_point t8;
    std::chrono::duration<double> time_span1(0);
    std::chrono::duration<double> time_span2(0);
    std::chrono::duration<double> time_span3(0);
    std::chrono::duration<double> time_span4(0);
    
    //Info<<endl;
    for(int i=0;i<points.size();i++)
    {
        //Info<<"Working on Point: "<<i<<" "<<points[i]<<endl;
        
        
        t1 = std::chrono::high_resolution_clock::now();
        std::unique_ptr<labelList> firstOrderNearNurbs = MainTree->nearNurbsCurves(points[i]);
        if(firstOrderNearNurbs->size() == 0)
        {
            pointDist[i] = 1000;
            //Info<<"\tSkipped because far away: "<<points[i]<<endl;
            continue;
        }
        //Info<<"\tIs inside: "<<points[i]<<endl;

        //Info<<"\tGot list size:"<<firstOrderNearNurbs->size()<<endl;
        t2 = std::chrono::high_resolution_clock::now();
        time_span1 += t2-t1;
        
        t3 = std::chrono::high_resolution_clock::now();
        DynamicList<scalar> distToNurbsSurface;
        distToNurbsSurface.setCapacity(10);
        DynamicList<scalar> paraToNurbsSurface;
        paraToNurbsSurface.setCapacity(10);
        bool outSideNurbsBox = false;
        for(int k=0;k<firstOrderNearNurbs->size();k++)
        {
            t5 = std::chrono::high_resolution_clock::now();
            label thisNurbs = (*firstOrderNearNurbs)[k];
            scalar thisNodePara = NurbsTrees[thisNurbs]->closestParaOnNurbsToPoint(points[i]);
            //Info<<"\tIndex of nurbs:"<<thisNurbs<<" with para: "<<thisNodePara<<endl;
            if(thisNodePara < this->Curves[thisNurbs]->min_U())
            {
                pointDist[i] = 1;
                outSideNurbsBox = true;
                continue;
            }
            t6 = std::chrono::high_resolution_clock::now();
            time_span3 += t6-t5; 
            
            t7 = std::chrono::high_resolution_clock::now();
            paraToNurbsSurface.append(thisNodePara);
            distToNurbsSurface.append(this->Curves[thisNurbs]->distanceToNurbsSurface(thisNodePara,points[i]));
            t8 = std::chrono::high_resolution_clock::now();
            time_span4 += t8-t7; 
        }
        if(outSideNurbsBox)
            continue;
        
        scalar minDistToNurbsSurface = std::numeric_limits<scalar>::max();
        for(int k=0;k<distToNurbsSurface.size();k++)
        {
            if(distToNurbsSurface[k] < minDistToNurbsSurface)
                minDistToNurbsSurface = distToNurbsSurface[k];
        }
        t4 = std::chrono::high_resolution_clock::now();
        time_span2 += t4-t3;
        
        pointDist[i] = minDistToNurbsSurface;
        
        //Info<<"Finished working on Point: "<<points[i]<<" "<<pointDist[i]<<endl;
        /*
        if(points[i] == vector(1.6,-0.1,-0.1))
        {
            FatalErrorInFunction
            << "Temp stop"<<endl
            << exit(FatalError);
        }
        */
    }
    Info<<endl;
    Info<<"Kd-Tree took \t\t\t\t" << time_span1.count() << " seconds."<<endl;
    Info<<"BsTree took \t\t\t\t" << time_span3.count() << " seconds."<<endl;
    Info<<"Newton took \t\t\t\t" << time_span4.count() << " seconds."<<endl;
    Info<<"BsTree + Newton took \t\t\t" << time_span2.count() << " seconds."<<endl;

    /*
    Info<<"Final point writing"<<endl;
    for(int i=0;i<points.size();i++)
    {
        Info<<i<<" "<<points[i]<<" "<<pointDist[i]<<endl;
    }
    */
}
