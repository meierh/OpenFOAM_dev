#include "cutCellFvMesh.H"
void testForNonHexMesh(fvMesh& mesh);

Foam::cutCellFvMesh::cutCellFvMesh
(
    const IOobject& io,
    std::shared_ptr<std::vector<Nurbs>> Curves,
    Time& runTime,
    std::unique_ptr<volScalarField>& solidFraction
):
dynamicRefineFvMesh(io),
Curves(Curves),
MainTree(std::unique_ptr<KdTree>(new KdTree(this->Curves))),
NurbsTrees(List<std::unique_ptr<BsTree>>((*(this->Curves)).size()))
{
    for(unsigned long i=0;i<(*(this->Curves)).size();i++)
    {
        NurbsTrees[i] = std::move(std::unique_ptr<BsTree>(new BsTree((*(this->Curves))[i])));
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
    
    owner.append(addedCutFacesOwner);
    owner.append(splitAndUnsplitFacesInteriorOwner);
    owner.append(splitAndUnsplitFacesBoundaryOwner);
    
    neighbour.append(addedCutFacesNeighbor);
    neighbour.append(splitAndUnsplitFacesInteriorNeighbor);
    neighbour.append(splitAndUnsplitFacesBoundaryNeighbor);

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
    std::shared_ptr<std::vector<Nurbs>> Curves,
    cutStatus state
):
dynamicRefineFvMesh(io),
Curves(Curves),
MainTree(std::unique_ptr<KdTree>(new KdTree(this->Curves))),
NurbsTrees(List<std::unique_ptr<BsTree>>((*(this->Curves)).size())),
ibAlgorithm(state)
{
    //MainTree = std::unique_ptr<KdTree>(new KdTree(this->Curves));
    //NurbsTrees = List<std::unique_ptr<BsTree>>((*(this->Curves)).size());
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
    for(unsigned long i=0;i<(*(this->Curves)).size();i++)
    {
        NurbsTrees[i] = std::move(std::unique_ptr<BsTree>(new BsTree((*(this->Curves))[i])));
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
    if(ibAlgorithm == internalCut)
    {
        createNewMeshData();
    
        faces.append(addedCutFaces);
        faces.append(splitAndUnsplitFacesInterior);
        faces.append(splitAndUnsplitFacesBoundary);
    
        owner.append(addedCutFacesOwner);
        owner.append(splitAndUnsplitFacesInteriorOwner);
        owner.append(splitAndUnsplitFacesBoundaryOwner);
    
        neighbour.append(addedCutFacesNeighbor);
        neighbour.append(splitAndUnsplitFacesInteriorNeighbor);
        neighbour.append(splitAndUnsplitFacesBoundaryNeighbor);
    }
    else if(ibAlgorithm == delNegMesh)
    {
        Info<<"-------------------------------------------"<<endl;
        Info<<"Create new Mesh data and cut negative cells";
        t1 = std::chrono::high_resolution_clock::now();
        createNewMeshData_cutNeg_plus();
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        Info<< " took \t"<< time_span.count() << " seconds."<<endl;
        Info<<"-------------------------------------------"<<endl;

        Info<<"Combine resetPrimitives data"<<endl;
        t1 = std::chrono::high_resolution_clock::now();
        faces.append(splitAndUnsplitFacesInterior);
        faces.append(splitAndUnsplitFacesBoundary);
        faces.append(addedCutFaces);
        faces.append(splitAndUnsplitFacesInteriorToBoundary);
    
        owner.append(splitAndUnsplitFacesInteriorOwner);
        owner.append(splitAndUnsplitFacesBoundaryOwner);
        owner.append(addedCutFacesOwner);
        owner.append(splitAndUnsplitFacesInteriorToBoundaryOwner);
        
        neighbour.append(splitAndUnsplitFacesInteriorNeighbor);
        neighbour.append(splitAndUnsplitFacesBoundaryNeighbor);
        neighbour.append(addedCutFacesNeighbor);
        neighbour.append(splitAndUnsplitFacesInteriorToBoundaryNeighbor);
        
        for(int i=0;i<owner.size();i++)
        {
            if(owner[i]<0)
            {
                Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
                Info<<"faces["<<i<<"]:"<<faces[i]<<endl;
                Info<<"owner[i]:"<<owner[i]<<endl;
                Info<<"neighbour[i]:"<<neighbour[i]<<endl;            
                FatalErrorInFunction<<"Owner fail stop"<< exit(FatalError); 
            }
            if(neighbour[i]<-1)
            {
                Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
                Info<<"faces["<<i<<"]:"<<faces[i]<<endl;
                Info<<"owner[i]:"<<owner[i]<<endl;
                Info<<"neighbour[i]:"<<neighbour[i]<<endl;            
                FatalErrorInFunction<<"Neighbour fail stop"<< exit(FatalError); 
            }
        }

        /*
        Info<<"addedCutFaceOwner[4420]:"<<addedCutFacesOwner[4420]<<endl;
        Info<<"addedCutFaceNeighbor[4420]:"<<addedCutFacesNeighbor[4420]<<endl;
        Info<<"addedCutFaces[4420]:"<<addedCutFaces[4420]<<endl;
        
        Info<<"splitAndUnsplitFacesInterior: "<<splitAndUnsplitFacesInterior.size()<<endl;
        Info<<"splitAndUnsplitFacesBoundary: "<<splitAndUnsplitFacesBoundary.size()<<endl;
        Info<<"addedCutFaces: "<<addedCutFaces.size()<<endl;
        Info<<"splitAndUnsplitFacesInteriorToBoundary: "<<splitAndUnsplitFacesInteriorToBoundary.size()<<endl;
        
        FatalErrorInFunction<< "Temp stop"<<endl<< exit(FatalError);
        */
        
        Info<<"3"<<endl;
        for(int i=0;i<patchStarts.size();i++)
        {
            Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
        }

        /*
        patchStarts[patchStarts.size()-1] = patchStarts.last()+patchSizes.last();
        patchSizes[patchSizes.size()-1] = addedCutFaces.size()+splitAndUnsplitFacesInteriorToBoundary.size();
        */
        
        //Info<<endl;
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

    Info<<"Correcting face normal direction";
    t1 = std::chrono::high_resolution_clock::now();
    correctFaceNormalDir(newMeshPoints_,faces,owner,neighbour);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t\t\t" << time_span.count() << " seconds."<<endl;

    const pointField& oldPoints = this->points();
    const faceList& oldFaceList = this->faces();
    const cellList& oldCells = this->cells();
    
    oldCellVolume = scalarList(oldCells.size());
    for(int i=0;i<oldCellVolume.size();i++)
    {
        oldCellVolume[i] = oldCells[i].mag(oldPoints,oldFaceList);
    }
    
    testNewMeshData(faces,owner,neighbour,patchStarts,patchSizes);
        
    Info<<"Reset:"<<endl;
    resetPrimitives(Foam::clone(newMeshPoints_),
                    Foam::clone(faces),
                    Foam::clone(owner),
                    Foam::clone(neighbour),
                    patchSizes,
                    patchStarts,
                    true);
    Info<<"First self test"<<endl;
    selfTestMesh();
    
    const cellList& newCells = this->cells();
    newCellVolume = scalarList(newCells.size());
    for(int i=0;i<newCellVolume.size();i++)
    {
        Info<<"i:"<<i<<endl;
        newCellVolume[i] = newCells[i].mag(newMeshPoints_,this->faces());
    }
    
    Info<<"Agglomerate small cut-cells";
    t1 = std::chrono::high_resolution_clock::now();
    agglomerateSmallCells_cutNeg_plus(newCellVolume,oldCellVolume,partialThreeshold);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t" << time_span.count() << " seconds."<<endl;

    //printMesh();
    Info<<"Please write"<<endl;
    this->write();
    Info<<"Written"<<endl;
    //printMesh();
    selfTestMesh();
    Info<<"Ending"<<endl;
}
 
void Foam::cutCellFvMesh::projectNurbsSurface()
{
    const pointField& points = this->points();
    pointDist = scalarList(points.size());
    meshPointNurbsReference.setSize(points.size());
    
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
    std::unordered_set<label> twoCellInd;
    twoCellInd.insert(202310);
    twoCellInd.insert(202631);
    twoCellInd.insert(215792);
    twoCellInd.insert(215471);
    twoCellInd.insert(202630);
    twoCellInd.insert(215791);
    twoCellInd.insert(215470);
    twoCellInd.insert(202309);
    
    twoCellInd.insert(202311);
    twoCellInd.insert(202632);
    twoCellInd.insert(215793);
    twoCellInd.insert(215472);

    for(int i=0;i<points.size();i++)
    {
        //Info<<"Working on Point: "<<i<<" "<<points[i]<<endl;
        //label testInd = 202631;
        
        t1 = std::chrono::high_resolution_clock::now();
        std::unique_ptr<labelList> firstOrderNearNurbs = MainTree->nearNurbsCurves(points[i]);
        /*
        if(i==testInd)
        {
            Info<<"Near Nurbs: "<<firstOrderNearNurbs->size()<<endl;
            //FatalErrorInFunction<< " Temp stop."<< exit(FatalError);
        }
        */
        if(firstOrderNearNurbs->size() == 0)
        {
            pointDist[i] = std::numeric_limits<scalar>::max();
            //Info<<"\tSkipped because far away: "<<points[i]<<endl;
            continue;
        }
        //Info<<"\tIs inside: "<<points[i]<<endl;

        //Info<<"\tGot list size:"<<firstOrderNearNurbs->size()<<endl;
        t2 = std::chrono::high_resolution_clock::now();
        time_span1 += t2-t1;
        
        t3 = std::chrono::high_resolution_clock::now();
        DynamicList<scalar> distToNurbsSurface;
        DynamicList<scalar> paraToNurbsSurface;
        DynamicList<label> indToNurbsSurface;
        bool allOutSideNurbsBox = true;
        for(int k=0;k<firstOrderNearNurbs->size();k++)
        {
            t5 = std::chrono::high_resolution_clock::now();
            label thisNurbs = (*firstOrderNearNurbs)[k];
            scalar thisNodePara = NurbsTrees[thisNurbs]->closestParaOnNurbsToPoint(points[i]);
            //Info<<"\tIndex of nurbs:"<<thisNurbs<<" with para: "<<thisNodePara<<endl;
            if(thisNodePara < (*(this->Curves))[thisNurbs].min_U())
            {
                /*
                if(i==testInd)
                {
                    Info<<k<<" Nurbs dist smaller: "<<thisNodePara<<"/"<<this->Curves[thisNurbs]->min_U()<<endl;
                }
                */
                pointDist[i] = std::numeric_limits<scalar>::max();
                continue;
            }
            allOutSideNurbsBox = false;
            t6 = std::chrono::high_resolution_clock::now();
            time_span3 += t6-t5; 
            
            t7 = std::chrono::high_resolution_clock::now();
            /*
            if(i==testInd)
            {
                Info<<k<<" Nurbs dist smaller: "<<thisNodePara<<"/"<<this->Curves[thisNurbs]->min_U()<<"   NurbsPoint"<<
                this->Curves[thisNurbs]->Curve_Derivative(0,thisNodePara)<<"   point:"<<points[i]<<endl;
            }
            */
            paraToNurbsSurface.append(thisNodePara);
            distToNurbsSurface.append((*(this->Curves))[thisNurbs].distanceToNurbsSurface(thisNodePara,points[i]));
            indToNurbsSurface.append(thisNurbs);
            t8 = std::chrono::high_resolution_clock::now();
            time_span4 += t8-t7; 
        }
        if(allOutSideNurbsBox)
            continue;
        
        scalar minDistToNurbsSurface = std::numeric_limits<scalar>::max();
        scalar minDistparaToNurbsSurface;
        label minDistindToNurbsSurface;
        for(int k=0;k<distToNurbsSurface.size();k++)
        {
            if(distToNurbsSurface[k] < minDistToNurbsSurface)
            {
                minDistToNurbsSurface = distToNurbsSurface[k];
                minDistparaToNurbsSurface = paraToNurbsSurface[k];
                minDistindToNurbsSurface = indToNurbsSurface[k];
            }
        }
        t4 = std::chrono::high_resolution_clock::now();
        time_span2 += t4-t3;
        /*
        if(i==testInd)
        {
            Info<<"paraToNurbsSurface: "<<paraToNurbsSurface<<endl;
            Info<<"distToNurbsSurface: "<<distToNurbsSurface<<endl;
            Info<<"allOutSideNurbsBox: "<<allOutSideNurbsBox<<endl;
            Info<<"minDistToNurbsSurface: "<<minDistToNurbsSurface<<endl;
            //FatalErrorInFunction<< " Temp stop."<< exit(FatalError);
        }
        */
        
        pointDist[i] = minDistToNurbsSurface;
        nurbsReference temp;
        temp.nurbsInd = minDistindToNurbsSurface;
        temp.nurbsPara = minDistparaToNurbsSurface;
        meshPointNurbsReference[i] = temp;
        
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

Foam::cutCellFvMesh::cutCellFvMesh
(
    const IOobject& io,
    cutStatus state,
    scalar cellDimToStructureDimLimit
):
dynamicRefineFvMesh(io),
ibAlgorithm(state),
motionPtr_(motionSolver::New(*this,dynamicMeshDict())),
cellDimToStructureDimLimit(cellDimToStructureDimLimit)
{
    Info<<"cellDimToStructureDimLimit:"<<cellDimToStructureDimLimit<<endl;
    const objectRegistry& reg = this->thisDb();
    wordList namesIO = reg.names();
    const objectRegistry& parentReg = reg.parent();
    wordList parNamesIO = parentReg.names();
    const objectRegistry& parentparentReg = parentReg.parent();
    wordList parparNamesIO = parentparentReg.names();
    //Info<<"namesIO:"<<namesIO<<endl;
    //Info<<"parNamesIO:"<<parNamesIO<<endl;
    //Info<<"parparNamesIO:"<<parparNamesIO<<endl;
    
    bool refineIsHex = false;
    const dictionary& dynDict = this->dynamicMeshDict();
    if(dynDict.found("useHexTopology",true,true))
    {
        const entry& hexEntry = dynDict.lookupEntry("useHexTopology",true,true);
        if(hexEntry.isStream())
        {
            ITstream& hexTopStream = hexEntry.stream();
            token hexTopToken;
            hexTopStream.read(hexTopToken);
            if(hexTopToken.wordToken().match("true"))
                refineIsHex=true;
        }
    }
    if(!refineIsHex)
        FatalIOError<<"\"useHexTopology yes\" must be defined and set in dynamicMeshDict"<< exit(FatalIOError);     
    
    // Initialize mesh class
    Info<<"Read Nurbs"<<endl;
    fileName runDirectory = this->fvMesh::polyMesh::objectRegistry::rootPath();
    fileName caseName = this->fvMesh::polyMesh::objectRegistry::caseName();
    NurbsReader Reader(runDirectory,caseName);
    Curves = Reader.getNurbsCurves();
    MainTree = std::unique_ptr<KdTree>(new KdTree(this->Curves));
    NurbsTrees = List<std::unique_ptr<BsTree>>((*(this->Curves)).size());
    Info<<"Created Main Tree"<<endl;
    for(unsigned long i=0;i<(*(this->Curves)).size();i++)
    {
        NurbsTrees[i] = std::move(std::unique_ptr<BsTree>(new BsTree((*(this->Curves))[i])));
    }
    Info<<"Prepared all Data"<<endl;
    Info<<"Before refinement"<<endl;
    testForNonHexMesh(*this);    
    // Refine nurbs cut surface
    refineTheImmersedBoundary();
    Info<<"After refinement"<<endl;
    testForNonHexMesh(*this);
    
    //Apply the immersed boundary cut method
    cutTheImmersedBoundary();
    this->write();
}

void Foam::cutCellFvMesh::refineTheImmersedBoundary()
{
    scalar oldCellDimToStructureDim = std::numeric_limits<scalar>::max();
    label refinementIteration=0;
    while(true)
    {
        Info<<"-------------------------------------------------"<<endl;
        // Project distance to points
        std::chrono::high_resolution_clock::time_point t1;
        std::chrono::high_resolution_clock::time_point t2;
        std::chrono::duration<double> time_span;
        Info<<"Projection of distance field ";
        t1 = std::chrono::high_resolution_clock::now();
        projectNurbsSurface();
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
        
        // Collect refinement canidate cells
        Info<<"Collect all cutCells"<<endl;
        const pointField& points = this->points();
        const faceList& faces = this->faces();
        const cellList& cells = this->cells();
        const labelList& owner  = this->owner();
        const labelList& neighbour  = this->neighbour();
        const hexRef8& cutterEngine = this->meshCutter();
        const labelList& cellRefinementLevel = cutterEngine.cellLevel();
        DynamicList<label> iBCells;
        for(int i=0;i<cells.size();i++)
        {
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
            
            if(nullExist || (posExist && negExist))
            {
                iBCells.append(i);
            }
        }
        for(int i=0;i<iBCells.size();i++)
        {
            if(iBCells[i]<0 || iBCells[i]>=cells.size())
            {
                Info<<"iBCells[i]:"<<iBCells[i]<<endl;
                FatalErrorInFunction<<"Invalid cell"<< exit(FatalError);
            }
        }
        std::unordered_set<label> refineCellsMap;
        for(int i=0;i<iBCells.size();i++)
        {
            refineCellsMap.insert(iBCells[i]);
            for(int j=0;j<cells[iBCells[i]].size();j++)
            {
                label faceLabel = cells[iBCells[i]][j];
                label neighborCell = -1;
                if(owner[faceLabel]==iBCells[i] && faceLabel<neighbour.size())
                    neighborCell = neighbour[faceLabel];
                else if(neighbour[faceLabel]==iBCells[i])
                    neighborCell = owner[faceLabel];
                
                if(neighborCell<-1 || neighborCell>=cells.size())
                {
                    Info<<endl;
                    Info<<"iBCells[i]:"<<iBCells[i]<<endl;
                    Info<<"face:"<<faceLabel<<endl;
                    Info<<"owner[faceLabel]:"<<owner[faceLabel]<<endl;
                    Info<<"neighbour[faceLabel]:"<<neighbour[faceLabel]<<endl;
                    Info<<"neighbour[faceLabel]:"<<neighbour[faceLabel]<<endl;
                    Info<<"neighborCell:"<<neighborCell<<endl;
                    FatalErrorInFunction<<"False cell"<< exit(FatalError);
                }
                refineCellsMap.insert(neighborCell);
            }
        }
        refineCellsMap.erase(-1);
        DynamicList<label> refineCells;
        for(const label& cell:refineCellsMap)
            refineCells.append(cell);
        
        Info<<"refineCells.size():"<<refineCells.size()<<endl;
        //FatalErrorInFunction<<"Temp stop"<< exit(FatalError);
        
        // Collect reRefinement cells
        Info<<"Collect all reRefinment cells"<<endl;
        DynamicList<label> reRefinementCells;
        reRefinementCells.append(refineCells);
        std::unordered_set<label> reRefinementCellsMap;
        for(const label& cell: reRefinementCells)
            reRefinementCellsMap.insert(cell);
        bool reRefinementDone=false;
        label cellLoopStartInd=0;
        label cellLoopEndInd=reRefinementCells.size();
        while(true)
        {
            std::unordered_set<label> duplicateTest;
            for(const label& insTest: reRefinementCells)
            {
                if(duplicateTest.count(insTest)!=0)
                    FatalErrorInFunction<<"Duplicate"<< exit(FatalError);
                else
                    duplicateTest.insert(insTest);
            }
            
            Info<<"cellLoopStartInd:"<<cellLoopStartInd<<endl;
            Info<<"cellLoopEndInd:"<<cellLoopEndInd<<endl;
            bool cellsAppended = false;
            for(int i=cellLoopStartInd;i<cellLoopEndInd;i++)
            {
                label refineCellLevel = cellRefinementLevel[reRefinementCells[i]];
                for(int j=0;j<cells[reRefinementCells[i]].size();j++)
                {
                    label faceLabel = cells[reRefinementCells[i]][j];
                    label neighborCell = -1;
                    if(owner[faceLabel]==reRefinementCells[i] && faceLabel<neighbour.size())
                        neighborCell = neighbour[faceLabel];
                    else if(neighbour[faceLabel]==reRefinementCells[i])
                        neighborCell = owner[faceLabel];
                    
                    if(neighborCell!=-1 && reRefinementCellsMap.count(neighborCell)==0)
                    {
                        label neighbourRefineCellLevel = cellRefinementLevel[neighborCell];
                        label refineLvlDiff = refineCellLevel-neighbourRefineCellLevel;
                        if(refineLvlDiff>0 || refineLvlDiff<-1)
                        {
                            if(refineLvlDiff<-1 || refineLvlDiff>1)
                            {
                                Info<<"cell:"<<reRefinementCells[i]<<" refLevel:"<<refineCellLevel<<endl;
                                Info<<"neighborCell:"<<neighborCell<<" refLevel:"<<neighbourRefineCellLevel<<endl;
                                FatalErrorInFunction<<"Can not happen."<< exit(FatalError);
                            }
                            if(refineLvlDiff>0)
                            {
                                reRefinementCells.append(neighborCell);
                                reRefinementCellsMap.insert(neighborCell);
                                cellsAppended=true;
                            }
                        }
                    }
                }
            }
            if(cellsAppended)
            {
                cellLoopStartInd=cellLoopEndInd;
                cellLoopEndInd=reRefinementCells.size();
            }
            else
                break;
        }        
        DynamicList<label> tempCells;
        for(int i=refineCells.size();i<reRefinementCells.size();i++)
        {
            tempCells.append(reRefinementCells[i]);
            if(refineCellsMap.count(reRefinementCells[i])!=0)
                FatalErrorInFunction<<"reRefinement of refine cells"<< exit(FatalError);
        }
        reRefinementCells = tempCells;
        
        for(int i=0;i<reRefinementCells.size();i++)
        {
            if(reRefinementCells[i]<0 || reRefinementCells[i]>=cells.size())
            {
                Info<<"reRefinementCells[i]:"<<reRefinementCells[i]<<endl;
                FatalErrorInFunction<<"Invalid cell"<< exit(FatalError);
            }
        }
        
        Info<<"cells:"<<cells.size()<<endl;
        this->refine(reRefinementCells);
        const cellList& cells2 = this->cells();
        const faceList& faces2 = this->faces();
        Info<<"cells:"<<cells2.size()<<endl;        
        
        for(int i=0;i<refineCells.size();i++)
        {
            if(refineCells[i]<0 || refineCells[i]>=cells2.size())
            {
                Info<<"refineCellsMap:"<<refineCellsMap.count(refineCells[i])<<endl;
                Info<<"refineCells[i]:"<<refineCells[i]<<endl;
                FatalErrorInFunction<<"Invalid cell"<< exit(FatalError);
            }
        }
        
        Info<<"Test for cell size ratio"<<endl;
        // Test for cell size ratio
        scalar maxEdgeLen = std::numeric_limits<scalar>::min();
        for(int i=0;i<refineCells.size();i++)
        {
            //Info<<"Cell:"<<refineCells[i]<<endl;
            edgeList cellEdges = cells2[refineCells[i]].edges(faces2);
            for(const edge& oneEdge: cellEdges)
            {
                scalar oneEdgeLen = oneEdge.mag(points);
                maxEdgeLen = (maxEdgeLen<oneEdgeLen)?oneEdgeLen:maxEdgeLen;
            }        
        }
        Info<<"End"<<endl;
        scalar minRadius = minNurbsRadius();
        scalar cellDimToStructureDim = maxEdgeLen/minRadius;
        
        if(cellDimToStructureDim>=oldCellDimToStructureDim)
            FatalErrorInFunction<<"No progress in near nurbs refinement"<< exit(FatalError);

        Info<<"cellRefinementLevel[3874]:"<<cellRefinementLevel[3874]<<endl;
        Info<<"cellRefinementLevel[15225]:"<<cellRefinementLevel[15225]<<endl;
        Info<<"Refine"<<endl;
        Info<<cellDimToStructureDim<<"/"<<cellDimToStructureDimLimit<<endl;
        
        /*
        if(refinementIteration==2)
            FatalErrorInFunction<<"Temp stop"<< exit(FatalError);
        */
        
        if(cellDimToStructureDim > cellDimToStructureDimLimit)
        {
            this->refine(refineCells);
            refinementIteration++;
        }
        else
            break;
    }
    //this->write();    
}

void Foam::cutCellFvMesh::cutTheImmersedBoundary
(
)
{
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span;
    
    Info<<"Projection of distance field ";
    t1 = std::chrono::high_resolution_clock::now();
    projectNurbsSurface();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
    
    //checkForHexCellsInCutArea();
    
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
    if(ibAlgorithm == internalCut)
    {
        createNewMeshData();
    
        faces.append(addedCutFaces);
        faces.append(splitAndUnsplitFacesInterior);
        faces.append(splitAndUnsplitFacesBoundary);
    
        owner.append(addedCutFacesOwner);
        owner.append(splitAndUnsplitFacesInteriorOwner);
        owner.append(splitAndUnsplitFacesBoundaryOwner);
    
        neighbour.append(addedCutFacesNeighbor);
        neighbour.append(splitAndUnsplitFacesInteriorNeighbor);
        neighbour.append(splitAndUnsplitFacesBoundaryNeighbor);
    }
    else if(ibAlgorithm == delNegMesh)
    {
        Info<<"-------------------------------------------"<<endl;
        Info<<"Create new Mesh data and cut negative cells";
        t1 = std::chrono::high_resolution_clock::now();
        createNewMeshData_cutNeg_plus();
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        Info<< " took \t"<< time_span.count() << " seconds."<<endl;
        Info<<"-------------------------------------------"<<endl;

        Info<<"Combine resetPrimitives data"<<endl;
        t1 = std::chrono::high_resolution_clock::now();
        faces.append(splitAndUnsplitFacesInterior);
        faces.append(splitAndUnsplitFacesBoundary);
        faces.append(addedCutFaces);
        faces.append(splitAndUnsplitFacesInteriorToBoundary);
    
        owner.append(splitAndUnsplitFacesInteriorOwner);
        owner.append(splitAndUnsplitFacesBoundaryOwner);
        owner.append(addedCutFacesOwner);
        owner.append(splitAndUnsplitFacesInteriorToBoundaryOwner);
        
        neighbour.append(splitAndUnsplitFacesInteriorNeighbor);
        neighbour.append(splitAndUnsplitFacesBoundaryNeighbor);
        neighbour.append(addedCutFacesNeighbor);
        neighbour.append(splitAndUnsplitFacesInteriorToBoundaryNeighbor);
        
        for(int i=0;i<owner.size();i++)
        {
            if(owner[i]<0)
            {
                Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
                Info<<"faces["<<i<<"]:"<<faces[i]<<endl;
                Info<<"owner[i]:"<<owner[i]<<endl;
                Info<<"neighbour[i]:"<<neighbour[i]<<endl;            
                FatalErrorInFunction<<"Owner fail stop"<< exit(FatalError); 
            }
            if(neighbour[i]<-1)
            {
                Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
                Info<<"faces["<<i<<"]:"<<faces[i]<<endl;
                Info<<"owner[i]:"<<owner[i]<<endl;
                Info<<"neighbour[i]:"<<neighbour[i]<<endl;            
                FatalErrorInFunction<<"Neighbour fail stop"<< exit(FatalError); 
            }
        }

        /*
        Info<<"addedCutFaceOwner[4420]:"<<addedCutFacesOwner[4420]<<endl;
        Info<<"addedCutFaceNeighbor[4420]:"<<addedCutFacesNeighbor[4420]<<endl;
        Info<<"addedCutFaces[4420]:"<<addedCutFaces[4420]<<endl;
        
        Info<<"splitAndUnsplitFacesInterior: "<<splitAndUnsplitFacesInterior.size()<<endl;
        Info<<"splitAndUnsplitFacesBoundary: "<<splitAndUnsplitFacesBoundary.size()<<endl;
        Info<<"addedCutFaces: "<<addedCutFaces.size()<<endl;
        Info<<"splitAndUnsplitFacesInteriorToBoundary: "<<splitAndUnsplitFacesInteriorToBoundary.size()<<endl;
        
        FatalErrorInFunction<< "Temp stop"<<endl<< exit(FatalError);
        */

        
        Info<<"3"<<endl;
        for(int i=0;i<patchStarts.size();i++)
        {
            Info<<"BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
        }

        /*
        patchStarts[patchStarts.size()-1] = patchStarts.last()+patchSizes.last();
        patchSizes[patchSizes.size()-1] = addedCutFaces.size()+splitAndUnsplitFacesInteriorToBoundary.size();
        */
        
        //Info<<endl;
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

    Info<<"Correcting face normal direction";
    t1 = std::chrono::high_resolution_clock::now();
    correctFaceNormalDir(newMeshPoints_,faces,owner,neighbour);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t\t\t" << time_span.count() << " seconds."<<endl;

    const pointField& oldPoints = this->points();
    const faceList& oldFaceList = this->faces();
    const cellList& oldCells = this->cells();
    
    oldCellVolume = scalarList(oldCells.size());
    for(int i=0;i<oldCellVolume.size();i++)
    {
        oldCellVolume[i] = oldCells[i].mag(oldPoints,oldFaceList);
    }
    
    testNewMeshData(faces,owner,neighbour,patchStarts,patchSizes);
        
    Info<<"Reset:"<<endl;
    resetPrimitives(Foam::clone(newMeshPoints_),
                    Foam::clone(faces),
                    Foam::clone(owner),
                    Foam::clone(neighbour),
                    patchSizes,
                    patchStarts,
                    true);
    Info<<"First self test"<<endl;
    selfTestMesh();
    
    const cellList& newCells = this->cells();
    newCellVolume = scalarList(newCells.size());
    for(int i=0;i<newCellVolume.size();i++)
    {
        Info<<"i:"<<i<<endl;
        newCellVolume[i] = newCells[i].mag(newMeshPoints_,this->faces());
    }
    
    Info<<"Agglomerate small cut-cells";
    t1 = std::chrono::high_resolution_clock::now();
    agglomerateSmallCells_cutNeg_plus(newCellVolume,oldCellVolume,partialThreeshold);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t" << time_span.count() << " seconds."<<endl;

    //printMesh();
    Info<<"Please write"<<endl;
    this->write();
    Info<<"Written"<<endl;
    //printMesh();
    selfTestMesh();
    Info<<"Ending"<<endl;
}

bool Foam::cutCellFvMesh::update()
{
    dynamicFvMesh::movePoints(motionPtr_->newPoints());
    return true;
}

bool Foam::cutCellFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    motionPtr_->write();
    return fvMesh::writeObject(fmt, ver, cmp, write);
}

void Foam::cutCellFvMesh::moveTheMesh()
{
    point minVec = point(-1,-1,-1);
    scalar minLen = std::numeric_limits<scalar>::max();
    point avgVec = point(0,0,0);
    label cntAvg = 0;
    point maxVec = point(0,0,0);
    scalar maxLen = 0;
    
    const pointField& meshPoints = this->points();
    pointField newPoints(meshPoints.size());
    
    if(meshPoints.size()!=meshPointNurbsReference.size())
        FatalErrorInFunction<<"Unequal point number for reference ("<<meshPointNurbsReference.size()<<") and points ("<<meshPoints.size()<<"). Can not happen!"<< exit(FatalError);  
    
    for(int i=0;i<meshPoints.size();i++)
    {
        newPoints[i] = meshPoints[i];
        if(pointsToSide_[i]==0)
        {
            nurbsReference reference = meshPointNurbsReference[i];
            if(reference.nurbsInd==-1 || reference.nurbsPara==-1)
                FatalErrorInFunction<<"Zero point has no reference nurbs. Can not happen!"<< exit(FatalError);            
            vector pointMotion = (*(this->Curves))[reference.nurbsInd].movementVector(reference.nurbsPara);
            if(Foam::mag(pointMotion)>1)
            {
                Info<<"meshPoints["<<i<<"]:"<<meshPoints[i]<<endl;
                Info<<"reference.nurbsInd:"<<reference.nurbsInd<<endl;
                Info<<"reference.nurbsPara:"<<reference.nurbsPara<<endl;
                Info<<"pointMotion:"<<endl;
                FatalErrorInFunction<<"Temp Stop!"<< exit(FatalError);            
            }
            newPoints[i] +=  pointMotion;
            
            
            if(maxLen<Foam::mag(pointMotion))
            {
                maxLen = Foam::mag(pointMotion);
                maxVec = pointMotion;
            }
            if(minLen>Foam::mag(pointMotion))
            {
                minLen = Foam::mag(pointMotion);
                minVec = pointMotion;
            }
            avgVec+=pointMotion;
            cntAvg++;
        }
    }
    avgVec/=cntAvg;
    
    //Info<<"maxVec:"<<maxVec<<endl;
    //Info<<"minVec:"<<minVec<<endl;
    //Info<<"avgVec:"<<avgVec<<endl;
    
    motionPtr_->movePoints(newPoints);
    //Info<<"Moved "<<newPoints.size()<<" points"<<endl;
}

void Foam::cutCellFvMesh::moveNurbsCurves
(
    List<List<vector>> movedControlPoints
)
{
    if((*(this->Curves)).size()!=static_cast<unsigned long>(movedControlPoints.size()))
        FatalErrorInFunction<<"Given number of controlPoint update blocks must equal the number of Nurbs"<< exit(FatalError);
    for(int i=0;i<movedControlPoints.size();i++)
    {
        (*(this->Curves))[i].moveNurbs(movedControlPoints[i]);
    }
}

void testForNonHexMesh(fvMesh& mesh)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    label cubeCells=0;
    label nonCubeCells=0;
    label nonCubeBorderCell=0;
        
    for(int i=0;i<cells.size();i++)
    {
        bool nonCubeCell = false;
        labelList vertices = cells[i].labels(faces);
        label nVertices = vertices.size();
        
        edgeList edges = cells[i].edges(faces);
        label nEdges = edges.size();
        
        label nFaces = cells[i].size();
        
        if(nVertices!=8 || nEdges!=12 || nFaces!=6)
        {
            nonCubeCell==true;
        }
        
        bool notFourVerticeFace = false;
        for(int j=0;j<cells[i].size();j++)
        {
            face oneFace = faces[cells[i][j]];
            if(oneFace.size()!=4)
            {
                notFourVerticeFace=true;
            }
        }
        if(notFourVerticeFace||nonCubeCell)
        {
            nonCubeCells++;
            continue;
        }
        cubeCells++;
    }
    
    Info<<nonCubeCells<<" non cube cells"<<endl;
    Info<<cubeCells<<" cube cells"<<endl;
}

void Foam::cutCellFvMesh::checkForHexCellsInCutArea()
{
    const cellList& cells = this->cells();
    const faceList& faces = this->faces();
    label cubeCells=0;
    label nonCubeCells=0;
    label nonCubeBorderCell=0;
    
    for(int i=0;i<cells.size();i++)
    {
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
        
        if(nullExist || (posExist && negExist))
        {
            bool nonCubeCell = false;
            labelList vertices = cells[i].labels(faces);
            label nVertices = vertices.size();
            
            edgeList edges = cells[i].edges(faces);
            label nEdges = edges.size();
            
            label nFaces = cells[i].size();
            
            if(nVertices!=8 || nEdges!=12 || nFaces!=6)
            {
                nonCubeCell==true;
            }
            
            bool notFourVerticeFace = false;
            for(int j=0;j<cells[i].size();j++)
            {
                face oneFace = faces[cells[i][j]];
                if(oneFace.size()!=4)
                {
                    notFourVerticeFace=true;
                }
            }
            if(notFourVerticeFace||nonCubeCell)
            {
                FatalErrorInFunction<<"Non cube cell in immersed boundary cut area"<< exit(FatalError);
            }
        }
    }
}

scalar Foam::cutCellFvMesh::minNurbsRadius()
{
    scalar minRadius = std::numeric_limits<scalar>::max();
    for(const Nurbs& oneNurbs: *Curves)
    {
        scalar radius = oneNurbs.radius();
        minRadius = (radius<minRadius)?radius:minRadius;
    }
    return minRadius;
}
