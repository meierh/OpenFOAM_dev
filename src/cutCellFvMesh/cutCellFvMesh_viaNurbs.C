#include "cutCellFvMesh.H"
void testForNonHexMesh(fvMesh& mesh);

Foam::scalar normEuler(Foam::vector pnt)
{
    return std::sqrt(pnt.x()*pnt.x()+pnt.y()*pnt.y()+pnt.z()*pnt.z());
}

Foam::cutCellFvMesh::cutCellFvMesh
(
    const IOobject& io,
    std::shared_ptr<std::vector<Nurbs1D>> Curves,
    Time& runTime,
    std::unique_ptr<volScalarField>& solidFraction
):
dynamicRefineFvMesh(io),
marchingCubesAlgorithm(*this),
Curves(Curves),
MainTree(std::unique_ptr<KdTree>(new KdTree(this->Curves))),
NurbsTrees(List<std::unique_ptr<BsTree>>((*(this->Curves)).size()))
{
    FatalErrorInFunction<<"Depreceated!!!"<< exit(FatalError);
    
    for(unsigned long i=0;i<(*(this->Curves)).size();i++)
    {
        NurbsTrees[i] = std::move(std::unique_ptr<BsTree>(new BsTree((*(this->Curves))[i])));
    }
    Info<<"Prepared all Data"<<endl;
    
    intersectionRadius = 0.02;
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
    std::shared_ptr<std::vector<Nurbs1D>> Curves,
    cutStatus state
):
dynamicRefineFvMesh(io),
marchingCubesAlgorithm(*this),
Curves(Curves),
MainTree(std::unique_ptr<KdTree>(new KdTree(this->Curves))),
NurbsTrees(List<std::unique_ptr<BsTree>>((*(this->Curves)).size())),
ibAlgorithm(state),
motionPtr_(motionSolver::New(*this,dynamicMeshDict()))
{
    FatalErrorInFunction<<"Depreceated!!!"<< exit(FatalError); 

    
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span;
    
    intersectionRadius = 0;
    
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
 
void Foam::cutCellFvMesh::projectNurbsSurface(bool reset)
{
    const pointField& points = this->points();
    pointDist = scalarList(points.size());
    
    DynamicList<DynamicList<nurbsReference>> meshPointNurbsReferenceOld = meshPointNurbsReference;
    scalarList pointDistMapOld = pointDistMap;

    meshPointNurbsReference.clear();
    meshPointNurbsReference.setSize(points.size());
    pointDistMap.clear();
    pointDistMap.setSize(points.size());
    
    if(reset)
    {
        newToOldPointIndMap.clear();
    }
    
    std::chrono::high_resolution_clock::time_point t1,t2,t3,t4,t5,t6,t7,t8;
    std::chrono::duration<double> time_span1(0),time_span2(0),time_span3(0),time_span4(0);
    
    label newCalculatedPnts = 0;
    for(int i=0;i<points.size();i++)
    {
        auto oldInd = newToOldPointIndMap.find(i);
        if(oldInd!=newToOldPointIndMap.end())
        {            
            scalar distRes = pointDistMapOld[oldInd->second];
            meshPointNurbsReference[i] = meshPointNurbsReferenceOld[oldInd->second];
            pointDist[i] = distRes;
            pointDistMap[i] = pointDist[i];

            continue;
        }
        
        //Info<<"After"<<endl;
        
        newCalculatedPnts++;
        
        t1 = std::chrono::high_resolution_clock::now();
        std::unique_ptr<labelList> firstOrderNearNurbs = MainTree->nearNurbsCurves(points[i]);
        if(firstOrderNearNurbs->size() == 0)
        {
            pointDist[i] = std::numeric_limits<scalar>::max();
            meshPointNurbsReference[i].clear();
            pointDistMap[i] = pointDist[i];
            continue;
        }
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
            if(thisNodePara < (*(this->Curves))[thisNurbs].min_U())
            {
                pointDist[i] = std::numeric_limits<scalar>::max();
                continue;
            }
            allOutSideNurbsBox = false;
            t6 = std::chrono::high_resolution_clock::now();
            time_span3 += t6-t5; 
            
            t7 = std::chrono::high_resolution_clock::now();
            paraToNurbsSurface.append(thisNodePara);
            distToNurbsSurface.append((*(this->Curves))[thisNurbs].distanceToNurbsSurface(thisNodePara,points[i]));
            indToNurbsSurface.append(thisNurbs);
            t8 = std::chrono::high_resolution_clock::now();
            time_span4 += t8-t7; 
        }
        if(allOutSideNurbsBox)
        {
            pointDist[i] = std::numeric_limits<scalar>::max();
            meshPointNurbsReference[i].clear();
            pointDistMap[i] = pointDist[i];
            continue;
        }
        
        scalar minDistToNurbsSurface = std::numeric_limits<scalar>::max();
        scalar minDistparaToNurbsSurface = -1;
        label minDistindToNurbsSurface = -1;
        for(int k=0;k<distToNurbsSurface.size();k++)
        {
            if(distToNurbsSurface[k] < minDistToNurbsSurface)
            {
                minDistToNurbsSurface = distToNurbsSurface[k];
                minDistparaToNurbsSurface = paraToNurbsSurface[k];
                minDistindToNurbsSurface = indToNurbsSurface[k];
            }
        }

        bool nonSecondNurbs = true;
        scalar secondMinDistToNurbsSurface = std::numeric_limits<scalar>::max();
        scalar secondMinDistparaToNurbsSurface = -1;
        label secondMinDistindToNurbsSurface = -1;
        for(int k=0;k<distToNurbsSurface.size();k++)
        {
            if(indToNurbsSurface[k]!=minDistindToNurbsSurface && distToNurbsSurface[k]<secondMinDistToNurbsSurface)
            {
                nonSecondNurbs = false;
                secondMinDistToNurbsSurface = distToNurbsSurface[k];
                secondMinDistparaToNurbsSurface = paraToNurbsSurface[k];
                secondMinDistindToNurbsSurface = indToNurbsSurface[k];
            }
        }

        if(nonSecondNurbs)
        {
            pointDist[i] = minDistToNurbsSurface;
            
            meshPointNurbsReference[i].clear();
            nurbsReference temp;
            temp.nurbsInd = minDistindToNurbsSurface;
            temp.nurbsPara = minDistparaToNurbsSurface;
            meshPointNurbsReference[i].append(temp);
            pointDistMap[i] = pointDist[i];
        }
        else
        {
            if(secondMinDistToNurbsSurface < minDistToNurbsSurface)
                FatalErrorInFunction<<"Second smallest dist smaller than smallest one. Can not happen!"<< exit(FatalError);
            
            vector vecToMinDistNurbs = (*Curves)[minDistindToNurbsSurface].Curve_Derivative(0,minDistparaToNurbsSurface);
            vecToMinDistNurbs = vecToMinDistNurbs - points[i];
            vector vecToSecondMinDistNurbs = (*Curves)[secondMinDistindToNurbsSurface].Curve_Derivative(0,secondMinDistparaToNurbsSurface);
            vecToSecondMinDistNurbs = vecToSecondMinDistNurbs - points[i];
            
            scalar angle;
            bool vecToNurbsZeroOnce = false;
            if(normEuler(vecToMinDistNurbs) * normEuler(vecToSecondMinDistNurbs) != 0)
                angle  = (vecToMinDistNurbs && vecToSecondMinDistNurbs) / (normEuler(vecToMinDistNurbs) * normEuler(vecToSecondMinDistNurbs));
            else
            {
                angle = -1;
                vecToNurbsZeroOnce = true;
            }
            scalar radiusFactor = ((angle+1.)/2.);
            scalar smoothingRadius = radiusFactor * intersectionRadius;
            //Info<<"vecToMinDistNurbs:"<<vecToMinDistNurbs<<" vecToSecondMinDistNurbs:"<<vecToSecondMinDistNurbs;
            //Info<<"  agl:"<<angle<<" smoRad:"<<smoothingRadius<<"  ";
            
            if(vecToNurbsZeroOnce && std::abs(minDistToNurbsSurface)<smoothingRadius && std::abs(secondMinDistToNurbsSurface)<smoothingRadius)
            {
                // rounding gets reducing further away from zero surface
                scalar roundingFactor = 1-(std::abs(minDistindToNurbsSurface)/smoothingRadius);
                
                // rounding gets scaled in respect of similarity of distance measure
                scalar distFirstToSecondMin  = std::abs(minDistToNurbsSurface - secondMinDistToNurbsSurface);
                scalar distFirstToSecondMinFactor = distFirstToSecondMin / (2*smoothingRadius);

                meshPointNurbsReference[i].clear();
                nurbsReference temp1;
                temp1.nurbsInd = minDistindToNurbsSurface;
                temp1.nurbsPara = minDistparaToNurbsSurface;
                meshPointNurbsReference[i].append(temp1);
                nurbsReference temp2;
                temp2.nurbsInd = secondMinDistindToNurbsSurface;
                temp2.nurbsPara = secondMinDistparaToNurbsSurface;
                meshPointNurbsReference[i].append(temp2);
                
                pointDist[i] = minDistToNurbsSurface + smoothingRadius*(1-distFirstToSecondMinFactor)*roundingFactor;
                pointDistMap[i] = pointDist[i];
            }
            else
            {
                meshPointNurbsReference[i].clear();
                nurbsReference temp1;
                temp1.nurbsInd = minDistindToNurbsSurface;
                temp1.nurbsPara = minDistparaToNurbsSurface;
                meshPointNurbsReference[i].append(temp1);
                
                pointDist[i] = minDistToNurbsSurface;
                pointDistMap[i] = pointDist[i];
            }
        }
        
        t4 = std::chrono::high_resolution_clock::now();
        time_span2 += t4-t3;
    }
    Info<<endl;
    Info<<"Kd-Tree took \t\t\t\t" << time_span1.count() << " seconds."<<endl;
    Info<<"BsTree took \t\t\t\t" << time_span3.count() << " seconds."<<endl;
    Info<<"Newton took \t\t\t\t" << time_span4.count() << " seconds."<<endl;
    Info<<"BsTree + Newton took \t\t\t" << time_span2.count() << " seconds."<<endl;
    Info<<"newCalculatedPnts:"<<newCalculatedPnts<<" / "<<points.size()<<endl;
}

Foam::cutCellFvMesh::cutCellFvMesh
(
    const IOobject& io,
    cutStatus state,
    scalar cellDimToStructureDimLimit
):
dynamicRefineFvMesh(io),
ibAlgorithm(state),
marchingCubesAlgorithm(*this),
motionPtr_(motionSolver::New(*this,dynamicMeshDict())),
cellDimToStructureDimLimit(cellDimToStructureDimLimit)
{  
    intersectionRadius = 0;
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
    cutTheImmersedBoundary_MC33();
    
    
    Info<<"pointDist[4608]:"<<pointDist[4608]<<Foam::endl;    
    const fvBoundaryMesh& bound = this->boundary();
    const label IBPatchID = bound.findPatchID("cutCell");
    Info<<"IBPatchID:"<<IBPatchID<<Foam::endl;
    const fvBoundaryMesh& boundary = this->boundary();
    const fvPatch& nurbsBoundary = boundary[IBPatchID];
    const faceList& faces = this->faces();
    label i=0;
    label faceInd=nurbsBoundary.start();
    Info<<"Collect faces"<<Foam::endl;
    Info<<"this->points().size():"<<this->points().size()<<Foam::endl;
    Info<<"meshPointNurbsReference.size():"<<meshPointNurbsReference.size()<<Foam::endl;
    Info<<"pointDist.size():"<<pointDist.size()<<Foam::endl;

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
        
        newToOldPointIndMap.clear();
        
        // Collect refinement canidate cells
        Info<<"Collect all cutCells"<<endl;
        const pointField& points = this->points();
        const faceList& faces = this->faces();
        const cellList& cells = this->cells();
        const labelList& owner  = this->faceOwner();
        const labelList& neighbour  = this->faceNeighbour();
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
        
        autoPtr<mapPolyMesh> reRefineMap = this->refine(reRefinementCells);
        this->motionPtr_->updateMesh(reRefineMap);
        const labelList& pntMapNewToOld = reRefineMap->pointMap();
        for(int i=0;i<pntMapNewToOld.size();i++)
        {
            if(pntMapNewToOld[i]!=-1)
            {
                newToOldPointIndMap[i] = pntMapNewToOld[i];
            }
        }       
        for(auto iter = newToOldPointIndMap.begin();iter!=newToOldPointIndMap.end();iter++)
        {
            if(iter->second>=meshPointNurbsReference.size())
            {
                Info<<"oldPointNbr:"<<reRefineMap->nOldPoints()<<endl;
                Info<<"newPointNbr:"<<reRefineMap->pointMap().size()<<endl;
                Info<<"pntMapNewToOld.size():"<<pntMapNewToOld.size()<<endl;
                Info<<"iter-second:"<<iter->second<<endl;
                Info<<"meshPointNurbsReference.size():"<<meshPointNurbsReference.size()<<endl;
                FatalErrorInFunction<<"Can not happen"<< exit(FatalError);
            }
        }
        //Info<<"pntMapNewToOld.size():"<<pntMapNewToOld.size()<<endl;
        //Info<<"meshPointNurbsReference.size():"<<meshPointNurbsReference.size()<<endl;
        
        const cellList& cells2 = this->cells();
        const faceList& faces2 = this->faces();
        const pointField& points2 = this->points();
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
        scalar minEdgeLen = std::numeric_limits<scalar>::max();
        scalar avgEdgeLen = 0;
        for(int i=0;i<refineCells.size();i++)
        {
            //Info<<"Cell:"<<refineCells[i]<<endl;
            edgeList cellEdges = cells2[refineCells[i]].edges(faces2);
            for(const edge& oneEdge: cellEdges)
            {
                scalar oneEdgeLen = oneEdge.mag(points2);
                maxEdgeLen = (maxEdgeLen<oneEdgeLen)?oneEdgeLen:maxEdgeLen;
                minEdgeLen = (minEdgeLen>oneEdgeLen)?oneEdgeLen:minEdgeLen;
                avgEdgeLen += oneEdgeLen;
            }
        }
        avgEdgeLen /= refineCells.size();
        Info<<"maxEdgeLen:"<<maxEdgeLen<<endl;
        Info<<"minEdgeLen:"<<minEdgeLen<<endl;
        Info<<"avgEdgeLen:"<<avgEdgeLen<<endl;
        scalar minRadius = minNurbsRadius();
        Info<<"minRadius:"<<minRadius<<endl;
        scalar cellDimToStructureDim = maxEdgeLen/minRadius;
        
        if(cellDimToStructureDim>=oldCellDimToStructureDim)
            FatalErrorInFunction<<"No progress in near nurbs refinement"<< exit(FatalError);

        Info<<cellDimToStructureDim<<"/"<<cellDimToStructureDimLimit<<endl;
        
        if(cellDimToStructureDim > cellDimToStructureDimLimit)
        {
            Info<<"Pre mesh cell nbr:"<<this->cells().size()<<endl;
            autoPtr<mapPolyMesh> refineMap = this->refine(refineCells);
            const labelList& pntMapNewToOld = reRefineMap->pointMap();
            for(int i=0;i<pntMapNewToOld.size();i++)
            {
                if(pntMapNewToOld[i]!=-1)
                {
                    label newPnt = i;
                    label oldPnt = pntMapNewToOld[i];
                    label oldOldPnt;
                    auto connex = newToOldPointIndMap.find(oldPnt);
                    if(connex != newToOldPointIndMap.end())
                    {
                        oldOldPnt = connex->second;
                        label mapNew = connex->first;
                        newToOldPointIndMap.erase(mapNew);
                        newToOldPointIndMap[newPnt] = oldOldPnt;
                    }
                }
            }
            for(auto iter = newToOldPointIndMap.begin();iter!=newToOldPointIndMap.end();iter++)
            {
                if(iter->second>=meshPointNurbsReference.size())
                {
                    Info<<"oldPointNbr:"<<reRefineMap->nOldPoints()<<endl;
                    Info<<"newPointNbr:"<<reRefineMap->pointMap().size()<<endl;
                    Info<<"pntMapNewToOld.size():"<<pntMapNewToOld.size()<<endl;
                    Info<<"iter-second:"<<iter->second<<endl;
                    Info<<"meshPointNurbsReference.size():"<<meshPointNurbsReference.size()<<endl;
                    FatalErrorInFunction<<"Can not happen"<< exit(FatalError);
                }
            }
            //Info<<"pntMapNewToOld.size():"<<pntMapNewToOld.size()<<endl;
            //Info<<"meshPointNurbsReference.size():"<<meshPointNurbsReference.size()<<endl;
            //FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);

            this->motionPtr_->updateMesh(refineMap);
            refinementIteration++;
            Info<<"----------------------------------------"<<endl;
        }
        else
            break;
    }
    //this->write();    
}

void Foam::cutCellFvMesh::cutTheImmersedBoundary()
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
        
    pointField points(0);
    faceList faces(0);
    labelList owner(0);
    labelList neighbour(0);
    List<std::unordered_map<label,label>> oldPointIndToPatchInd;
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
        
        //pointsToSide_ = pointsToSide;
        Info<<"newMeshPoints_.size():"<<newMeshPoints_.size()<<endl;
        
        List<bool> pntDeleted(newMeshPoints_.size(),true);
        for(const face& oneFace: faces)
            for(const label& oneVertice: oneFace)
                pntDeleted[oneVertice] = false;
                
        labelList pntOldIndToNewInd(newMeshPoints_.size(),-1);
        label index=0;
        for(int i=0;i<pntDeleted.size();i++)
        {
            if(!pntDeleted[i])
            {
                pntOldIndToNewInd[i] = index;
                index++;
            }
        }
        points.setSize(index+1);
        index=0;
        for(int i=0;i<pntDeleted.size();i++)
        {
            if(!pntDeleted[i])
            {
                points[index] = newMeshPoints_[i];
                index++;
            }
        }
        
        for(face& oneFace: faces)
        {
            for(label& oneVertice: oneFace)
            {
                point oldPoint = newMeshPoints_[oneVertice];
                point newPoint = points[pntOldIndToNewInd[oneVertice]];
                if(oldPoint != newPoint)
                    FatalErrorInFunction<< "Can not happen"<<endl<< exit(FatalError);
                oneVertice = pntOldIndToNewInd[oneVertice];
            }
        }
        Info<<"Set Ref"<<endl;
        DynamicList<DynamicList<nurbsReference>> meshPointNurbsReference_new;
        meshPointNurbsReference_new.setSize(points.size());
        scalarList pointDist_new;
        pointDist_new.setSize(points.size());
        for(int i=0;i<meshPointNurbsReference.size();i++)
        {
            if(!pntDeleted[i])
            {                
                meshPointNurbsReference_new[pntOldIndToNewInd[i]] = meshPointNurbsReference[i];
                pointDist_new[pntOldIndToNewInd[i]] = pointDist[i];
            }
        }
        meshPointNurbsReference = meshPointNurbsReference_new;
        pointDist = pointDist_new;
                        
        const polyBoundaryMesh& boundaryMesh = this->boundaryMesh();
        oldPointIndToPatchInd.setSize(boundaryMesh.size());
        for(int i=0;i<boundaryMesh.size();i++)
        {
            const polyPatch& onePatch = boundaryMesh[i];
            const labelList& patchPoints = onePatch.boundaryPoints();
            for(int j=0;j<patchPoints.size();j++)
                oldPointIndToPatchInd[i].insert(std::pair<label,label>(patchPoints[j],j));
        }
                
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        Info<< " took \t\t\t"<< time_span.count() << " seconds."<<endl;
    }
    //printMesh();
    //FatalErrorInFunction<<"Temporary stop!"<<exit(FatalError);

    Info<<"Correcting face normal direction";
    t1 = std::chrono::high_resolution_clock::now();
    correctFaceNormalDir(points,faces,owner,neighbour);
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
        if(oldCellVolume[i]==0.0)
        {
            Info<<"oldCells["<<i<<"]:"<<oldCells[i]<<endl;
            FatalErrorInFunction<< "Temp stop"<<endl<< exit(FatalError);
        }   
    }
    
    testNewMeshData(faces,owner,neighbour,patchStarts,patchSizes);
    
    resetPrimitives(Foam::clone(points),
                    Foam::clone(faces),
                    Foam::clone(owner),
                    Foam::clone(neighbour),
                    patchSizes,
                    patchStarts,
                    false);
    Info<<"Reset"<<endl;
    
    std::unordered_set<label> activePnts;
    for(int i=0;i<this->faces().size();i++)
        for(int j=0;j<this->faces()[i].size();j++)
            activePnts.insert(this->faces()[i][j]);
    label cnt=0;
    for(int i=0;i<points.size();i++)
        if(activePnts.count(i)!=0)
            cnt++;
    cnt++;
    if(cnt!=this->points().size())
    {
        Info<<"cnt:"<<cnt<<endl;
        Info<<"this->points().size():"<<this->points().size()<<endl;
        FatalErrorInFunction<< "Invalid point number"<<endl<< exit(FatalError);
    }
    
    this->topoChanging(true);

    //Reset field size for motionSolver
    //Begin
    Foam::motionSolver* rawPtr = motionPtr_.ptr();  
    displacementLaplacianFvMotionSolver* dMS;
    try{
        dMS = dynamic_cast<displacementLaplacianFvMotionSolver*>(rawPtr);
        if(dMS==NULL)
            FatalErrorInFunction<< "Cast to displacementMotionSolver failed. Must use the one"<<endl<< exit(FatalError);
    }catch(...){
        FatalErrorInFunction<< "Cast to displacementMotionSolver failed. Must use the one"<<endl<< exit(FatalError);
    }    
    pointVectorField& pVF = dMS->pointDisplacement();
    if(this->points().size()!=pVF.size())
        pVF.setSize(this->points().size());
    pointVectorField::Boundary& boundFieldPointDispl = pVF.boundaryFieldRef();
    //Info<<"boundFieldPointDispl.size():"<<boundFieldPointDispl.size()<<endl;
    for(int i=0;i<boundFieldPointDispl.size();i++)
    {
        pointPatchField<vector>* boundaryPatchField = &boundFieldPointDispl[i];
        const pointPatch& boundaryPointPatch = boundaryPatchField->patch();
        const labelList& patchPoints = boundaryPointPatch.meshPoints();
        fixedValuePointPatchField<vector>* fVPPF;
        try{
            fVPPF = dynamic_cast<fixedValuePointPatchField<vector>*>(boundaryPatchField);
            if(fVPPF==NULL)
                FatalErrorInFunction<< "Cast to valuePointPatchField failed. Must use the one"<<endl<< exit(FatalError);
        }catch(...){
            FatalErrorInFunction<< "Cast to valuePointPatchField failed. Must use the one"<<endl<< exit(FatalError);
        }
        fVPPF->setSize(patchPoints.size());       
    }
    volVectorField& vVF = dMS->cellDisplacement();
    if(this->cells().size()!=vVF.size())
        vVF.setSize(this->cells().size());
    volVectorField::Boundary& vVF_Bound = vVF.boundaryFieldRef();
    for(int i=0;i<vVF_Bound.size();i++)
    {
        fvPatchField<vector>& boundField = vVF_Bound[i];
        const fvPatch& boundPatch = boundField.patch();
        if(boundPatch.size()!=boundField.size())
            boundField.setSize(boundPatch.size());
    }
    pointField& pointsNull = dMS->points0();
    pointsNull.setSize(this->points().size());
    for(int i=0;i<this->points().size();i++)
        pointsNull[i] = this->points()[i];
    motionPtr_.set(rawPtr);
    //End
   
    Info<<"First self test"<<endl;
    selfTestMesh();
    
    //Agglomeration for too small cells
    //Begin
    Info<<"Agglomerate small cut-cells";
    const cellList& newCells = this->cells();
    newCellVolume = scalarList(newCells.size());
    for(int i=0;i<newCellVolume.size();i++)
    {
        newCellVolume[i] = newCells[i].mag(points,this->faces());
    }
    t1 = std::chrono::high_resolution_clock::now();
    agglomerateSmallCells_cutNeg_plus(newCellVolume,oldCellVolume,partialThreeshold);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t" << time_span.count() << " seconds."<<endl;
    //End

    //printMesh();
    Info<<"Please write"<<endl;
    this->write();
    Info<<"Written"<<endl;
    //printMesh();
    selfTestMesh();
    Info<<"Ending"<<endl;
}

void Foam::cutCellFvMesh::cutTheImmersedBoundary_MC33()
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
    
    Info<<"Execution of Marching Cubes";
    t1 = std::chrono::high_resolution_clock::now();
    executeMarchingCubes();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t" << time_span.count() << " seconds."<<endl;
        
    Info<<"Adding of cut points";
    t1 = std::chrono::high_resolution_clock::now();
    newMeshPoints_MC33();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t\t" << time_span.count() << " seconds."<<endl;
    
    Info<<"Adding of cut edges";
    t1 = std::chrono::high_resolution_clock::now();
    newMeshEdges_MC33();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t\t" << time_span.count() << " seconds."<<endl;
        
    edgesToSide();
    
    Info<<"Adding of cut faces";
    t1 = std::chrono::high_resolution_clock::now();
    newMeshFaces_MC33();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t\t"<< time_span.count() << " seconds."<<endl;
        
    Info<<"Cutting old faces";
    t1 = std::chrono::high_resolution_clock::now();
    cutOldFaces_MC33();
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t\t\t" << time_span.count() << " seconds."<<endl;
        
    pointField points(0);
    faceList faces(0);
    labelList owner(0);
    labelList neighbour(0);
    List<std::unordered_map<label,label>> oldPointIndToPatchInd;
    
    Info<<"-------------------------------------------"<<endl;
    Info<<"Create new Mesh data and cut negative cells";
    t1 = std::chrono::high_resolution_clock::now();
    createNewMeshData_MC33();
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
    
    //pointsToSide_ = pointsToSide;
    Info<<"newMeshPoints_.size():"<<newMeshPoints_.size()<<endl;
    
    List<bool> pntDeleted(newMeshPoints_.size(),true);
    for(const face& oneFace: faces)
        for(const label& oneVertice: oneFace)
            pntDeleted[oneVertice] = false;
            
    labelList pntOldIndToNewInd(newMeshPoints_.size(),-1);
    label index=0;
    for(int i=0;i<pntDeleted.size();i++)
    {
        if(!pntDeleted[i])
        {
            pntOldIndToNewInd[i] = index;
            index++;
        }
    }
    points.setSize(index+1);
    index=0;
    for(int i=0;i<pntDeleted.size();i++)
    {
        if(!pntDeleted[i])
        {
            points[index] = newMeshPoints_[i];
            index++;
        }
    }
    
    for(face& oneFace: faces)
    {
        for(label& oneVertice: oneFace)
        {
            point oldPoint = newMeshPoints_[oneVertice];
            point newPoint = points[pntOldIndToNewInd[oneVertice]];
            if(oldPoint != newPoint)
                FatalErrorInFunction<< "Can not happen"<<endl<< exit(FatalError);
            oneVertice = pntOldIndToNewInd[oneVertice];
        }
    }
    Info<<"Set Ref"<<endl;
    DynamicList<DynamicList<nurbsReference>> meshPointNurbsReference_new;
    meshPointNurbsReference_new.setSize(points.size());
    scalarList pointDist_new;
    pointDist_new.setSize(points.size());
    for(int i=0;i<meshPointNurbsReference.size();i++)
    {
        if(!pntDeleted[i])
        {                
            meshPointNurbsReference_new[pntOldIndToNewInd[i]] = meshPointNurbsReference[i];
            pointDist_new[pntOldIndToNewInd[i]] = pointDist[i];
        }
    }
    meshPointNurbsReference = meshPointNurbsReference_new;
    pointDist = pointDist_new;
                    
    const polyBoundaryMesh& boundaryMesh = this->boundaryMesh();
    oldPointIndToPatchInd.setSize(boundaryMesh.size());
    for(int i=0;i<boundaryMesh.size();i++)
    {
        const polyPatch& onePatch = boundaryMesh[i];
        const labelList& patchPoints = onePatch.boundaryPoints();
        for(int j=0;j<patchPoints.size();j++)
            oldPointIndToPatchInd[i].insert(std::pair<label,label>(patchPoints[j],j));
    }
            
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t"<< time_span.count() << " seconds."<<endl;

    //printMesh();
    //FatalErrorInFunction<<"Temporary stop!"<<exit(FatalError);

    Info<<"Correcting face normal direction";
    t1 = std::chrono::high_resolution_clock::now();
    correctFaceNormalDir(points,faces,owner,neighbour);
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
        if(oldCellVolume[i]==0.0)
        {
            Info<<"oldCells["<<i<<"]:"<<oldCells[i]<<endl;
            FatalErrorInFunction<< "Temp stop"<<endl<< exit(FatalError);
        }   
    }
    
    testNewMeshData(faces,owner,neighbour,patchStarts,patchSizes);
    
    resetPrimitives(Foam::clone(points),
                    Foam::clone(faces),
                    Foam::clone(owner),
                    Foam::clone(neighbour),
                    patchSizes,
                    patchStarts,
                    false);
    Info<<"Reset"<<endl;
    
    std::unordered_set<label> activePnts;
    for(int i=0;i<this->faces().size();i++)
        for(int j=0;j<this->faces()[i].size();j++)
            activePnts.insert(this->faces()[i][j]);
    label cnt=0;
    for(int i=0;i<points.size();i++)
        if(activePnts.count(i)!=0)
            cnt++;
    cnt++;
    if(cnt!=this->points().size())
    {
        Info<<"cnt:"<<cnt<<endl;
        Info<<"this->points().size():"<<this->points().size()<<endl;
        FatalErrorInFunction<< "Invalid point number"<<endl<< exit(FatalError);
    }
    
    this->topoChanging(true);

    //Reset field size for motionSolver
    //Begin
    Foam::motionSolver* rawPtr = motionPtr_.ptr();  
    displacementLaplacianFvMotionSolver* dMS;
    try{
        dMS = dynamic_cast<displacementLaplacianFvMotionSolver*>(rawPtr);
        if(dMS==NULL)
            FatalErrorInFunction<< "Cast to displacementMotionSolver failed. Must use the one"<<endl<< exit(FatalError);
    }catch(...){
        FatalErrorInFunction<< "Cast to displacementMotionSolver failed. Must use the one"<<endl<< exit(FatalError);
    }    
    pointVectorField& pVF = dMS->pointDisplacement();
    if(this->points().size()!=pVF.size())
        pVF.setSize(this->points().size());
    pointVectorField::Boundary& boundFieldPointDispl = pVF.boundaryFieldRef();
    //Info<<"boundFieldPointDispl.size():"<<boundFieldPointDispl.size()<<endl;
    for(int i=0;i<boundFieldPointDispl.size();i++)
    {
        pointPatchField<vector>* boundaryPatchField = &boundFieldPointDispl[i];
        const pointPatch& boundaryPointPatch = boundaryPatchField->patch();
        const labelList& patchPoints = boundaryPointPatch.meshPoints();
        fixedValuePointPatchField<vector>* fVPPF;
        try{
            fVPPF = dynamic_cast<fixedValuePointPatchField<vector>*>(boundaryPatchField);
            if(fVPPF==NULL)
                FatalErrorInFunction<< "Cast to valuePointPatchField failed. Must use the one"<<endl<< exit(FatalError);
        }catch(...){
            FatalErrorInFunction<< "Cast to valuePointPatchField failed. Must use the one"<<endl<< exit(FatalError);
        }
        fVPPF->setSize(patchPoints.size());       
    }
    volVectorField& vVF = dMS->cellDisplacement();
    if(this->cells().size()!=vVF.size())
        vVF.setSize(this->cells().size());
    volVectorField::Boundary& vVF_Bound = vVF.boundaryFieldRef();
    for(int i=0;i<vVF_Bound.size();i++)
    {
        fvPatchField<vector>& boundField = vVF_Bound[i];
        const fvPatch& boundPatch = boundField.patch();
        if(boundPatch.size()!=boundField.size())
            boundField.setSize(boundPatch.size());
    }
    pointField& pointsNull = dMS->points0();
    pointsNull.setSize(this->points().size());
    for(int i=0;i<this->points().size();i++)
        pointsNull[i] = this->points()[i];
    motionPtr_.set(rawPtr);
    //End
   
    Info<<"First self test"<<endl;
    selfTestMesh();
    
    //Agglomeration for too small cells
    //Begin
    Info<<"Agglomerate small cut-cells";
    const cellList& newCells = this->cells();
    newCellVolume = scalarList(newCells.size());
    for(int i=0;i<newCellVolume.size();i++)
    {
        newCellVolume[i] = newCells[i].mag(points,this->faces());
    }
    t1 = std::chrono::high_resolution_clock::now();
    agglomerateSmallCells_cutNeg_plus(newCellVolume,oldCellVolume,partialThreeshold);
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< " took \t\t\t" << time_span.count() << " seconds."<<endl;
    //End

    //printMesh();
    Info<<"Please write"<<endl;
    //this->write();
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
    //Info<<"owner:"<<this->owner()<<endl;    
    //FatalErrorInFunction<< "Temp stop"<<endl<< exit(FatalError);
    
    point minVec = point(-1,-1,-1);
    scalar minLen = std::numeric_limits<scalar>::max();
    point avgVec = point(0,0,0);
    label cntAvg = 0;
    point maxVec = point(0,0,0);
    scalar maxLen = 0;
    
    const pointField& meshPoints = this->points();
    pointField newPoints(meshPoints.size());
    
    Info<<"meshPoints.size():"<<meshPoints.size()<<endl;
    Info<<"pointsToSide_.size():"<<pointsToSide_.size()<<endl;
    
    if(meshPoints.size()!=meshPointNurbsReference.size())
        FatalErrorInFunction<<"Unequal point number for reference ("<<meshPointNurbsReference.size()<<") and points ("<<meshPoints.size()<<"). Can not happen!"<< exit(FatalError);  
    
    std::unordered_map<label,vector> cutCellPointToMotion;
    for(int i=0;i<meshPoints.size();i++)
    {
        newPoints[i] = meshPoints[i];
        if(pointsToSide_[i]==0)
        {
            DynamicList<nurbsReference>& reference = meshPointNurbsReference[i];
            if(reference.size()==0)
                FatalErrorInFunction<<"Zero point has no reference nurbs. Can not happen!"<< exit(FatalError);
            vector pointMotion = vector(0,0,0);
            for(int j=0;j<reference.size();j++)
            {
                pointMotion += (*(this->Curves))[reference[j].nurbsInd].movementVector(reference[j].nurbsPara);
            }
            pointMotion /= reference.size();
            newPoints[i] +=  pointMotion;
            cutCellPointToMotion.insert(std::pair<label,vector>(i,pointMotion));
            
            
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
    
    Info<<"maxVec:"<<maxVec<<endl;
    Info<<"minVec:"<<minVec<<endl;
    Info<<"avgVec:"<<avgVec<<endl;
    Info<<"cntAvg:"<<cntAvg<<endl;
    
    Foam::motionSolver* rawPtr = motionPtr_.ptr();  
    Info<<"rawPtr null:"<<(rawPtr==0)<<endl;
    displacementLaplacianFvMotionSolver* dMS;
    try{
        dMS = dynamic_cast<displacementLaplacianFvMotionSolver*>(rawPtr);
        if(dMS==NULL)
            FatalErrorInFunction<< "Cast to displacementMotionSolver failed. Must use the one"<<endl<< exit(FatalError);
    }catch(...){
        FatalErrorInFunction<< "Cast to displacementMotionSolver failed. Must use the one"<<endl<< exit(FatalError);
    }
    
    pointVectorField& pVF = dMS->pointDisplacement();
    Info<<" nbrPoints:"<<this->points().size()<<endl;
    if(this->points().size()!=pVF.size())
        pVF.setSize(this->points().size());
    pointVectorField::Boundary& boundFieldPointDispl = pVF.boundaryFieldRef();
    Info<<"boundFieldPointDispl.size():"<<boundFieldPointDispl.size()<<endl;
    pointPatchField<vector>* boundaryPatchField = &boundFieldPointDispl[boundFieldPointDispl.size()-1];
    const pointPatch& boundaryPointPatch = boundaryPatchField->patch();
    const labelList& patchPoints = boundaryPointPatch.meshPoints();
    Info<<"   boundaryPatchfield size:"<<boundaryPatchField->size();
    Info<<"   meshPoints size:"<<patchPoints.size();
    fixedValuePointPatchField<vector>* fVPPF;
    try{
        fVPPF = dynamic_cast<fixedValuePointPatchField<vector>*>(boundaryPatchField);
        if(fVPPF==NULL)
            FatalErrorInFunction<< "Cast to valuePointPatchField failed. Must use the one"<<endl<< exit(FatalError);
    }catch(...){
        FatalErrorInFunction<< "Cast to valuePointPatchField failed. Must use the one"<<endl<< exit(FatalError);
    }
    Info<<"    fixedValuePointPatchField:"<<fVPPF->size()<<endl;
    Info<<"cutCellPointToMotion.size():"<<cutCellPointToMotion.size()<<endl;
    if(patchPoints.size()!=cutCellPointToMotion.size())
        FatalErrorInFunction<<"Wrong patch"<<endl<< exit(FatalError);
    
    Field<vector> movingPointsField(patchPoints.size());
    for(int j=0;j<patchPoints.size();j++)
    {
        movingPointsField[j] = cutCellPointToMotion[patchPoints[j]];
    }
    label intFieldSize = fVPPF->primitiveField().size();

    //Info<<"    movingPointsField:"<<movingPointsField.size()<<endl;
    //Info<<"movingPointsField:"<<movingPointsField<<endl;
    fVPPF->setSize(movingPointsField.size());
    for(int i=0;i<fVPPF->size();i++)
        (*fVPPF)[i] = movingPointsField[i];
    //valuePointPatchField<vector>::operator==(movingPointsField);
    Info<<"    fixedValuePointPatchField:"<<fVPPF->size()<<endl;
    Info<<"    fVPPF->updated():"<<fVPPF->updated()<<endl;
    if(fVPPF->updated())
        FatalErrorInFunction<< "Already updated. Can not happen!"<<endl<< exit(FatalError);
    
    /*
    volVectorField& vVF = dMS->cellDisplacement();
    Info<<"nbrCells:"<<this->cells().size()<<endl;
    Info<<"cellDisp Field size:"<<vVF.size()<<endl;
    Info<<"cellDisp Field:"<<vVF.primitiveFieldRef()<<endl;
    if(this->cells().size()!=vVF.size())
        vVF.setSize(this->cells().size());
        //vVF.Field<vector>::operator==(Field<vector>(this->cells().size()));
    Info<<"cellDisp Field size:"<<vVF.size()<<endl;
    Info<<"cellDisp Field:"<<vVF.primitiveFieldRef()<<endl;
    volVectorField::Boundary& vVF_Bound = vVF.boundaryFieldRef();
    for(int i=0;i<vVF_Bound.size();i++)
    {
        fvPatchField<vector>& boundField = vVF_Bound[i];
        const fvPatch& boundPatch = boundField.patch();
        Info<<"boundPatch size:"<<boundPatch.size()<<endl;
        //Info<<"boundField size:"<<boundField<<endl;
        Info<<"boundField size:"<<boundField.size()<<endl;
        if(boundPatch.size()!=boundField.size())
            boundField.setSize(boundPatch.size());
        Info<<"boundField size:"<<boundField.size()<<endl;
        //Info<<"boundField size:"<<boundField<<endl;
        Info<<"---"<<endl;
    }
    */
    
    
    //FatalErrorInFunction<< "Temp stop"<<endl<< exit(FatalError);
    
    
    //rawPtr->solve();
    
    motionPtr_.set(rawPtr);
    
    //Info<<"owner:"<<this->faceOwner()<<endl;

    
    //FatalErrorInFunction<< "Temp stop"<<endl<< exit(FatalError);
    
    Info<<"Moved Set movement to field"<<endl;
    
    /*
    surfaceScalarField* weights_ = new surfaceScalarField
     (
         IOobject
         (
             "weights",
             this->pointsInstance(),
             *this,
             IOobject::NO_READ,
             IOobject::NO_WRITE,
             false // Do not register
         ),
         *this,
         dimless
     );
     surfaceScalarField& weights = *weights_;
 
     // Set local references to mesh data
     // (note that we should not use fvMesh sliced fields at this point yet
     //  since this causes a loop when generating weighting factors in
     //  coupledFvPatchField evaluation phase)
     const labelUList& owner = this->faceOwner();
     const labelUList& neighbour = this->faceNeighbour();
 
     const vectorField& Cf = this->faceCentres();
     const vectorField& C = this->cellCentres();
     const vectorField& Sf = this->faceAreas();
 
     // ... and reference to the internal field of the weighting factors
     scalarField& w = weights.primitiveFieldRef();

         Info<<"1"<<endl;

    Info<<"owner:"<<owner<<endl;
     forAll(owner, facei)
     {
         Info<<endl<<"facei:"<<facei<<endl;
         Info<<"Sf[facei]:"<<Sf[facei]<<endl;
         Info<<"owner[facei]:"<<owner[facei]<<" / "<<C.size()<<endl;
         Info<<"C[owner[facei]]:"<<C[owner[facei]]<<endl;
         Info<<"C[neighbour[facei]]:"<<C[neighbour[facei]]<<endl;
         // Note: mag in the dot-product.
         // For all valid meshes, the non-orthogonality will be less that
         // 90 deg and the dot-product will be positive.  For invalid
         // meshes (d & s <= 0), this will stabilise the calculation
         // but the result will be poor.
         const scalar SfdOwn = mag(Sf[facei]&(Cf[facei] - C[owner[facei]]));
         const scalar SfdNei = mag(Sf[facei]&(C[neighbour[facei]] - Cf[facei]));
         const scalar SfdOwnNei = SfdOwn + SfdNei;
 
         Info<<facei<<"  SfdOwn:"<<SfdOwn<<"     SfdNei:"<<SfdNei<<"     SfdOwnNei:"<<SfdOwnNei<<"     vGreat:"<<vGreat<<endl;
         if (SfdNei/vGreat < SfdOwnNei)
         {
            Info<<"   "<<facei<<"  SfdOwn:"<<SfdOwn<<"     SfdNei:"<<SfdNei<<"     SfdOwnNei:"<<SfdOwnNei<<endl;
             w[facei] = SfdNei/SfdOwnNei;
         }
         else
         {
             const scalar dOwn = mag(Cf[facei] - C[owner[facei]]);
             const scalar dNei = mag(C[neighbour[facei]] - Cf[facei]);
             const scalar dOwnNei = dOwn + dNei;
             
            Info<<"   "<<facei<<"  dOwn:"<<dOwn<<"     dNei:"<<dNei<<"     dOwnNei:"<<dOwnNei<<endl;
             
             w[facei] = dNei/dOwnNei;
         }
     }
         Info<<"2"<<endl;

    
     surfaceScalarField::Boundary& wBf =
         weights.boundaryFieldRef();
 
     */
        /*
     forAll(this->boundary(), patchi)
     {
         this->boundary()[patchi].makeWeights(wBf[patchi]);
     }
        */

        
    /*
    Info<<"this->points.size():"<<this->points().size()<<endl;
    Info<<"this->nPoints():"<<this->nPoints()<<endl;
    Info<<"this->oldPoints().size():"<<this->oldPoints().size()<<endl;
    
    
    tmp<pointField> x = motionPtr_->newPoints();	
    Info<<"newPoints size:"<<x->size()<<endl;
    
    
    rawPtr = motionPtr_.ptr();  
    //Info<<"rawPtr null:"<<(rawPtr==0)<<endl;
    try{
        dMS = dynamic_cast<displacementLaplacianFvMotionSolver*>(rawPtr);
    }catch(...){
        FatalErrorInFunction<< "Cast to displacementMotionSolver failed. Must use the one"<<endl<< exit(FatalError);
    }    
    pointVectorField& pVF2 = dMS->pointDisplacement();
    Info<<"pointDisplacement.size():"<<pVF2.primitiveField().size()<<endl;
    Info<<"points0().size():"<<dMS->points0().size()<<endl;

    motionPtr_.set(rawPtr);        
    
    Info<<"this->points.size():"<<this->points().size()<<endl;
    Info<<"this->nPoints():"<<this->nPoints()<<endl;
    Info<<"this->oldPoints().size():"<<this->oldPoints().size()<<endl;
    
    FatalErrorInFunction<< "Temp stop"<<endl<< exit(FatalError);
    Info<<"Moved Set movement to field"<<endl;
    */
}

void Foam::cutCellFvMesh::moveNurbsCurves
(
    std::vector<List<List<vector>>>& movedDeformationControlPoints
)
{
    Info<<"Mesh assign Deformation curve in cutCellFvMesh"<<Foam::endl;
    if((*(this->Curves)).size()!=static_cast<unsigned long>(movedDeformationControlPoints.size()))
    {
        Info<<"(*(this->Curves)).size():"<<(*(this->Curves)).size()<<Foam::endl;
        Info<<"movedDeformationControlPoints.size():"<<movedDeformationControlPoints.size()<<Foam::endl;
        FatalErrorInFunction<<"Given number of controlPoint update blocks must equal the number of Nurbs"<< exit(FatalError);
    }
    for(int i=0;i<movedDeformationControlPoints.size();i++)
    {
        Info<<"Assign Deformation curve "<<i<<" with "<<movedDeformationControlPoints[i][0].size()<<" CPs"<<Foam::endl;
        (*(this->Curves))[i].moveNurbs(movedDeformationControlPoints[i]);
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
    for(const Nurbs1D& oneNurbs1D: *Curves)
    {
        scalar radius = oneNurbs1D.radius();
        minRadius = (radius<minRadius)?radius:minRadius;
    }
    return minRadius;
}

const DynamicList<DynamicList<Foam::cutCellFvMesh::nurbsReference>>& Foam::cutCellFvMesh::getMeshPointNurbsReference() const
{
    return meshPointNurbsReference;
}

const std::shared_ptr<std::vector<Nurbs1D>> Foam::cutCellFvMesh::getCurves()
const
{
    return Curves;
}

void Foam::cutCellFvMesh::setInitialDeformationCurve
(
    std::vector<scalarList>& nurbs_to_knots,
    std::vector<List<vector>>& nurbs_to_controlPoints,
    std::vector<scalarList>& nurbs_to_weights,
    std::vector<label>& nurbs_to_degree
)
{
    label nbrOfDefNurbsCurves = nurbs_to_knots.size();
    if(nbrOfDefNurbsCurves!=nurbs_to_controlPoints.size() || nbrOfDefNurbsCurves!=nurbs_to_weights.size() ||            
       nbrOfDefNurbsCurves!=nurbs_to_degree.size()        || nbrOfDefNurbsCurves!=Curves->size())
    {
        Info<<"nbrOfDefNurbsCurves:"<<nbrOfDefNurbsCurves<<endl;
        Info<<"nurbs_to_controlPoints.size():"<<nurbs_to_controlPoints.size()<<endl;
        Info<<"nurbs_to_weights.size():"<<nurbs_to_weights.size()<<endl;
        Info<<"Curves->size():"<<Curves->size()<<endl;
        Info<<"nurbs_to_degree.size()"<<nurbs_to_degree.size()<<endl;
        FatalErrorInFunction<<"Wrong Deformation curve number"<< exit(FatalError);
    }

    for(int i=0;i<nbrOfDefNurbsCurves;i++)
    {
        (*Curves)[i].createDeformationCurve(nurbs_to_knots[i],nurbs_to_controlPoints[i],
                                            nurbs_to_weights[i],nurbs_to_degree[i]);
    }
}
