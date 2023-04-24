#include "cutCellFvMesh.H"

void Foam::cutCellFvMesh::Barrier(bool stop)
{
    labelList test(Pstream::nProcs(),0);
    test[Pstream::myProcNo()] = Pstream::myProcNo();
    Pstream::gatherList(test);
    Pstream::scatterList(test);
    std::cout<<"XXXXXXXXXXXXX---------------Barrier:"<<Pstream::myProcNo()<<"----"<<test[0]<<","<<test[1]<<","<<test[2]<<","<<test[3]<<"---------------XXXXXXXXXXXX"<<std::endl;
    fflush(stdout);
    Pstream::gatherList(test);
    if(stop)
        FatalErrorInFunction<<Pstream::myProcNo()<<"Temp Stop!"<< exit(FatalError);

}

Foam::scalar norm2(Foam::vector pnt)
{
    return std::sqrt(pnt.x()*pnt.x()+pnt.y()*pnt.y()+pnt.z()*pnt.z());
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
        if(lvlSet >= 0)
            pointsToSide[i] = 1;
        else if(lvlSet < 0)
            pointsToSide[i] = -1;
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
        pointsToSide_[i] = (lvlSet >= 0) * 1 + (lvlSet < 0) * -1;
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
        
        if(lvlSetStart >= 0 && lvlSetEnd >= 0)
        {
            edgesToSide[i] = +1;
        }
        else if(lvlSetStart < 0 && lvlSetEnd < 0)
        {
            edgesToSide[i] = -1;
        }
        else
        {
            if(lvlSetStart < 0 && lvlSetEnd == 0)
            {
                edgesToSide[i] = -1;
            }
            else if(lvlSetStart == 0 && lvlSetEnd < 0)
            {
                edgesToSide[i] = -1;
            }
            else
            {
                edgesToSide[i] = 0;
            }
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
    
    for(int i=0;i<edges.size();i++)
    {
        scalar lvlSetStart = pointDist[edges[i].start()];
        scalar lvlSetEnd = pointDist[edges[i].end()];
        
        if(lvlSetStart >= 0 && lvlSetEnd >= 0)
        {
            edgesToSide_[i] = +1;
        }
        else if(lvlSetStart < 0 && lvlSetEnd < 0)
        {
            edgesToSide_[i] = -1;
        }
        else
        {
            if(lvlSetStart < 0 && lvlSetEnd == 0)
            {
                edgesToSide_[i] = -1;
            }
            else if(lvlSetStart == 0 && lvlSetEnd < 0)
            {
                edgesToSide_[i] = -1;
            }
            else
            {
                edgesToSide_[i] = 0;
            }
        }
    }
}

void Foam::cutCellFvMesh::facesToSide
(
)
{
    const faceList& faces = this->faces();
    facesToSide_.setSize(faces.size());
    
    for(int i=0;i<faces.size();i++)
    {
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<faces[i].size();k++)
        {
            if(pointDist[faces[i][k]] > 0)
                posExist = true;
            else if(pointDist[faces[i][k]] < 0)
                negExist = true;
        }
        
        if(!negExist)
        {
            if(posExist)
                facesToSide_[i] = +1;
            else
                facesToSide_[i] = +2;
        }
        else
        {
            if(posExist)
            {
                facesToSide_[i] = 0;
            }
            else
            {
                facesToSide_[i] = -1;
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
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<faces[i].size();k++)
        {
            if(pointDist[faces[i][k]] > 0)
                posExist = true;
            else if(pointDist[faces[i][k]] < 0)
                negExist = true;
        }
        
        if(!negExist)
        {
            if(posExist)
                facesToSide[i] = +1;
            else
                facesToSide[i] = +2;
        }
        else
        {
            if(posExist)
            {
                facesToSide[i] = 0;
            }
            else
            {
                facesToSide[i] = -1;
            }
        }
    }
    this->facesToSide_ = facesToSide;
}

void Foam::cutCellFvMesh::cellsToSide
(
)
{
    const faceList& faces = this->faces();
    const cellList& cells = this->cells();
    labelList cellsToSide(cells.size());
    
    for(int i=0;i<cells.size();i++)
    {
        labelList cellLabels = cells[i].labels(faces);
        
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<cellLabels.size();k++)
        {
            if(pointDist[cellLabels[k]] > 0)
                posExist = true;
            else if(pointDist[cellLabels[k]] < 0)
                negExist = true;
        }
        
        if(!negExist)
        {
            cellsToSide[i] = +1;
        }
        else
        {
            if(posExist)
            {
                cellsToSide[i] = 0;
            }
            else
            {
                cellsToSide[i] = -1;
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
    const faceList& faces = this->faces();
    labelList cellsToSide(cells.size());
    
    for(int i=0;i<cells.size();i++)
    {
        labelList cellLabels = cells[i].labels(faces);
        
        bool posExist = false;
        bool negExist = false;
        for(int k=0;k<cellLabels.size();k++)
        {
            if(pointDist[cellLabels[k]] > 0)
                posExist = true;
            else if(pointDist[cellLabels[k]] < 0)
                negExist = true;
        }
        
        if(!negExist)
        {
            cellsToSide[i] = +1;
        }
        else
        {
            if(posExist)
            {
                cellsToSide[i] = 0;
            }
            else
            {
                cellsToSide[i] = -1;
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

void Foam::cutCellFvMesh::executeMarchingCubes()
{
    const cellList& meshCells = this->cells();
    const edgeList& basisEdges = this->edges();
    const faceList& basisFaces = this->faces();
    
    mc33CutCellData.setCapacity(meshCells.size());
    for(int i=0;i<meshCells.size();i++)
    {
        const cell& thisCell = meshCells[i];
        const labelList& cellLabels = thisCell.labels(basisFaces);
        bool onePlus=false;
        bool oneNeg=false;
        
        for(label cellPointLabel : cellLabels)
        {
            if(pointDist[cellPointLabel] >= 0)
                onePlus = true;
            else 
                oneNeg = true;
        }
        
        MC33::MC33Cube oneCube;
        if(onePlus && oneNeg)
        {
            oneCube = marchingCubesAlgorithm.computeCutCell(i);
            mc33CutCellData.append(oneCube);
        }
        else if(!onePlus && !oneNeg)
        {
            FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
        }
        else
        {
            mc33CutCellData.append(oneCube);
        }
    }
    
    const labelListList& edgeCells = this->edgeCells();
    for(int edgInd=0;edgInd<edgeCells.size();edgInd++)
    {
        const edge& currentEdge = basisEdges[edgInd];
        for(int i=0;i<edgeCells[edgInd].size();i++)
        {
            label cellInd = edgeCells[edgInd][i];
            MC33::MC33Cube& thisCellMc33Cube = mc33CutCellData[cellInd];
            if(thisCellMc33Cube.cell!=-1)
            {
                bool assignOne = false;
                for(int j=0;j<thisCellMc33Cube.edges.size();j++)
                {
                    edge& cubeEdge = thisCellMc33Cube.edges[j];
                    if(cubeEdge.otherVertex(currentEdge.start())!=-1 &&
                       cubeEdge.otherVertex(currentEdge.end())!=-1)
                    {
                        if(thisCellMc33Cube.edgeGlobalInd[j]!=-1)
                            FatalErrorInFunction<<"Double assigment"<< exit(FatalError);                        
                        thisCellMc33Cube.edgeGlobalInd[j] = edgInd;
                        assignOne=true;
                    }
                }
                if(!assignOne)
                    FatalErrorInFunction<<"Missing assignment"<< exit(FatalError);
            }
        }
    }
    
    for(int i=0;i<meshCells.size();i++)
    {            
        if(mc33CutCellData[i].cell!=-1)
        {
            for(int j=0;j<mc33CutCellData[i].edgeGlobalInd.size();j++)
            {
                if(mc33CutCellData[i].edgeGlobalInd[j]==-1)
                    FatalErrorInFunction<<"Missing assignment"<< exit(FatalError);
            }
        }
    }  
}

void Foam::cutCellFvMesh::newMeshPoints_MC33()
{
    const cellList& meshCells = this->cells();
    const pointField& basisPoints = this->points();
    //const faceList& basisFaces = this->faces();
    const edgeList& basisEdges = this->edges();
    const labelListList& edgeCells = this->edgeCells();
    
    List<List<label>> edgeToNeigboringMC33Cubes_faceToCutEdge_localInd(edgeCells.size());
    List<bool> edgeIsCutEdge(edgeCells.size(),false);

    for(int edgInd=0;edgInd<edgeCells.size();edgInd++)
    {
        //Compute local edgeInd of treated edge for all neighborCells and store in neigboringMC33Cubes_sharedEdge_localInd
        List<label> neigboringMC33Cubes_sharedEdge_localInd(edgeCells[edgInd].size(),-1);
        for(int i=0;i<edgeCells[edgInd].size();i++)
        {
            label cellInd = edgeCells[edgInd][i];
            if(mc33CutCellData[cellInd].cell!=-1)
            {
                //neighboringMC33Cubes.append(&mc33CutCellData[cellInd]);
                label MC33Cube_cellLocal_Edge = -1;
                for(int j=0;j<mc33CutCellData[cellInd].edgeGlobalInd.size();j++)
                {
                    if(mc33CutCellData[cellInd].edgeGlobalInd[j]==edgInd)
                    {
                        if(MC33Cube_cellLocal_Edge==-1)
                            MC33Cube_cellLocal_Edge = j;
                        else
                            FatalErrorInFunction<<"Double assignment!"<< exit(FatalError);
                    }
                }
                if(MC33Cube_cellLocal_Edge==-1)
                {
                    Info<<Foam::endl<<"edgInd:"<<edgInd<<Foam::endl;
                    Info<<"mc33CutCellData[cellInd].edgeGlobalInd:"<<mc33CutCellData[cellInd].edgeGlobalInd<<Foam::endl;
                    FatalErrorInFunction<<"Missing edge assignment!"<< exit(FatalError);
                }
                neigboringMC33Cubes_sharedEdge_localInd[i] = MC33Cube_cellLocal_Edge;
            }
        }

        //Test if edge is cut and check for consistency
        bool oneAssgined = false;
        bool oneNotAssigned = false;
        for(int i=0;i<edgeCells[edgInd].size();i++)
        {
            edgeToNeigboringMC33Cubes_faceToCutEdge_localInd[edgInd].resize(edgeCells[edgInd].size(),-1);
            label cellInd = edgeCells[edgInd][i];
            label sharedEdge_localInd = neigboringMC33Cubes_sharedEdge_localInd[i];
            if(mc33CutCellData[cellInd].cell!=-1)
            {
                DynamicList<label> MC33Cube_cellLocal_cutfaceInd;
                DynamicList<label> neighborhood_count;
                for(unsigned int j=0;j<mc33CutCellData[cellInd].cutTriangles.size();j++)
                {
                    auto& oneTriangle = mc33CutCellData[cellInd].cutTriangles[j];
                    label cutEdge1 = std::get<0>(oneTriangle);
                    label cutEdge2 = std::get<1>(oneTriangle);
                    label cutEdge3 = std::get<2>(oneTriangle);

                    if(cutEdge1==sharedEdge_localInd ||
                       cutEdge2==sharedEdge_localInd ||
                       cutEdge3==sharedEdge_localInd)
                    {
                        MC33Cube_cellLocal_cutfaceInd.append(j);
                        neighborhood_count.append(0);
                    }
                }
                for(int k1=0;k1<MC33Cube_cellLocal_cutfaceInd.size();k1++)
                {
                    label cutFaceInd1 = MC33Cube_cellLocal_cutfaceInd[k1];
                    auto& oneTriangle = mc33CutCellData[cellInd].cutTriangles[cutFaceInd1];
                    label cutEdge1 = std::get<0>(oneTriangle);
                    label cutEdge2 = std::get<1>(oneTriangle);
                    label cutEdge3 = std::get<2>(oneTriangle);
                    bool oneFaceMatch = false;
                    for(int k2=0;k2<MC33Cube_cellLocal_cutfaceInd.size();k2++)
                    {
                        if(k1==k2)
                            continue;
                        
                        label cutFaceInd2 = MC33Cube_cellLocal_cutfaceInd[k2];
                        std::unordered_set<label> oneTriangleSet = {cutEdge1,cutEdge2,cutEdge3};
                        auto& prevTriangle = mc33CutCellData[cellInd].cutTriangles[cutFaceInd2];
                        label prevTri_cutEdge1 = std::get<0>(prevTriangle);
                        label prevTri_cutEdge2 = std::get<1>(prevTriangle);
                        label prevTri_cutEdge3 = std::get<2>(prevTriangle);
                        bool facesMatch = (oneTriangleSet.count(prevTri_cutEdge1)&&
                                            oneTriangleSet.count(prevTri_cutEdge2))||        
                                          (oneTriangleSet.count(prevTri_cutEdge2)&&
                                            oneTriangleSet.count(prevTri_cutEdge3))||
                                          (oneTriangleSet.count(prevTri_cutEdge3)&&
                                            oneTriangleSet.count(prevTri_cutEdge1));
                        if(facesMatch)
                        {
                            oneFaceMatch=true;
                            neighborhood_count[k2]++;
                        }
                    }
                    if(!oneFaceMatch && MC33Cube_cellLocal_cutfaceInd.size()>1)
                    {
                        Info<<"neighborhood_count:"<<neighborhood_count<<Foam::endl;
                        Info<<"MC33Cube_cellLocal_cutfaceInd:"<<MC33Cube_cellLocal_cutfaceInd<<Foam::endl;
                        FatalErrorInFunction<<"Triangle with missing alignment!"<< exit(FatalError);
                    }
                }   
                for(label count: neighborhood_count)
                {
                    if(count>2)
                    {
                        Info<<"neighborhood_count:"<<neighborhood_count<<Foam::endl;
                        FatalErrorInFunction<<"Triangle can only border maximaly two other trianlges!"<< exit(FatalError);
                    }
                }
                if(MC33Cube_cellLocal_cutfaceInd.size()==0)
                    oneNotAssigned=true;
                if(MC33Cube_cellLocal_cutfaceInd.size()>0)
                    oneAssgined=true;
                edgeToNeigboringMC33Cubes_faceToCutEdge_localInd[edgInd][i] = sharedEdge_localInd;
            }
            else
            {
                oneNotAssigned=true;
            }
        }
        if((oneNotAssigned && oneAssgined)||(!oneNotAssigned && !oneAssgined))
        {
            Info<<"edgInd:"<<edgInd<<Foam::endl;
            Info<<"oneNotAssigned:"<<oneNotAssigned<<Foam::endl;
            Info<<"oneAssgined:"<<oneAssgined<<Foam::endl;
            FatalErrorInFunction<<"Inconsistent edge assignment!"<< exit(FatalError);
        }
        if(oneAssgined && !oneNotAssigned)
        {
            edgeIsCutEdge[edgInd] = true;
        }
    }
    
    nOldPoints = nbrOfPrevPoints = basisPoints.size();
    preMotionPoints = basisPoints;
    
    DynamicList<label>provisional_pointMap;
    for(label oldInd=0; oldInd<nOldPoints; oldInd++)
    {
        if(pointsToSide_[oldInd]==+1)
        {
            provisional_pointMap.append(oldInd);
        }
    }
    reversePointMap.setSize(nOldPoints,-1);
    for(label newInd=0; newInd<provisional_pointMap.size(); newInd++)
    {
        reversePointMap[provisional_pointMap[newInd]] = newInd;
    }

    newMeshPointsInFunc.append(basisPoints);
    pointsToSide(basisPoints);
    
    pointToEgde_.setSize(basisPoints.size(),-1);
    edgeToPoint_.setSize(basisEdges.size(),-1);
        
    for(int edgInd=0;edgInd<basisEdges.size();edgInd++)
    {
        if(edgeIsCutEdge[edgInd])
        {
            label startLabel = basisEdges[edgInd].start();
            label endLabel = basisEdges[edgInd].end();
            scalar phiStart = pointDist[startLabel];
            scalar phiEnd = pointDist[endLabel];
            
            label cutPointInd=-1;
            if(phiStart==0)
            {
                cutPointInd = startLabel;
            }
            else if(phiEnd==0)
            {
                cutPointInd = endLabel;
            }
            else if(phiStart==0 && phiEnd==0)
            {
                FatalErrorInFunction<<"Edge can not be cut!"<< exit(FatalError);
            }
            else
            {
                scalar scale = phiStart / (phiStart-phiEnd);
                if(scale>=1 || scale<=0)
                    FatalErrorInFunction<<"New Point must be inside edge!"<< exit(FatalError);
                vector startPoint = basisPoints[startLabel];
                vector endPoint = basisPoints[endLabel];
                vector startEndVector = endPoint-startPoint;
                vector newPoint = startPoint + scale * startEndVector;
                
                bool found;
                DynamicList<nurbsReference> reference;
                distToNurbs(newPoint,found,reference);
                if(!found)
                    FatalErrorInFunction<<"New Point have a dist to Nurbs!"<< exit(FatalError);

                provisional_pointMap.append(-1);
                cutPointInd = newMeshPointsInFunc.size();
                newMeshPointsInFunc.append(newPoint);
                meshPointNurbsReference.append(reference);            
                pointsToSide_.append(1);            
                pointToEgde_.append(edgInd);
                edgeToPoint_[edgInd] = cutPointInd;
            }
            
            for(int i=0;i<edgeCells[edgInd].size();i++)
            {
                label cellInd = edgeCells[edgInd][i];
                if(mc33CutCellData[cellInd].cell!=-1)
                {
                    label mc33CubeLocalEdgeInd = edgeToNeigboringMC33Cubes_faceToCutEdge_localInd[edgInd][i];
                    if(mc33CutCellData[cellInd].cutEdgeVerticeIndex[mc33CubeLocalEdgeInd]!=-1)
                        FatalErrorInFunction<<"Inconsistent vertice assignment!"<< exit(FatalError);
                    mc33CutCellData[cellInd].cutEdgeVerticeIndex[mc33CubeLocalEdgeInd] = cutPointInd;
                }
            }
        }
    }
    
    for(int cellInd=0; cellInd<meshCells.size(); cellInd++)
    {
        if(mc33CutCellData[cellInd].cell!=-1)
        {
            std::unordered_set<label> cutEdges;
            for(unsigned int j=0;j<mc33CutCellData[cellInd].cutTriangles.size();j++)
            {
                auto triangle = mc33CutCellData[cellInd].cutTriangles[j];
                cutEdges.insert(std::get<0>(triangle));
                cutEdges.insert(std::get<1>(triangle));
                cutEdges.insert(std::get<2>(triangle));
            }
            for(int j=0;j<mc33CutCellData[cellInd].cutEdgeVerticeIndex.size();j++)
            {
                if(cutEdges.find(j)!=cutEdges.end())
                {
                    if(mc33CutCellData[cellInd].cutEdgeVerticeIndex[j]==-1)
                        FatalErrorInFunction<<"Missing Vertice assignment!"<< exit(FatalError);
                }
            }
        }
    }
    pointMap = provisional_pointMap;
    
    
    newMeshPoints_ = pointField(pointMap.size());
    label newInd=0;
    for(;newInd<pointMap.size() && pointMap[newInd]!=-1; newInd++)
    {
        newMeshPoints_[newInd] = newMeshPointsInFunc[pointMap[newInd]];
    }
    if((pointMap.size()-newInd) != (newMeshPointsInFunc.size()-nOldPoints))
    {
        FatalErrorInFunction<<"Invalid!"<< exit(FatalError);
    }
    for(label i=0; newInd<pointMap.size(); i++,newInd++)
    {
        newMeshPoints_[newInd] = newMeshPointsInFunc[nOldPoints+i];
    }
    
    pointToEgde_.setCapacity(pointToEgde_.size());   
}

void Foam::cutCellFvMesh::printAddedPoints
(
)
{
    Info<<"---------------------------------------AddedPoints-----------------------------------"<<endl;
    for(int k=nbrOfPrevPoints;k<newMeshPoints_.size();k++)
    {
        Info<<"Point added: "<<k<<"-"<<newMeshPoints_[k]<<"\t";
        Info<<" at edge: ";
        Info<<pointToEgde_[k]<<endl;
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
    DynamicList<nurbsReference> throwAway;
    return distToNurbs(pnt,foundFlag,throwAway);
}

scalar Foam::cutCellFvMesh::distToNurbs
(
    point pnt,
    bool& foundFlag,
    DynamicList<nurbsReference>& reference
)
{
    //Info<<"pnt:"<<pnt<<endl;
    foundFlag = true;
    scalar dist;
    std::unique_ptr<labelList> firstOrderNearNurbs = MainTree->nearNurbsCurves(pnt);
    if(firstOrderNearNurbs->size()==0)
    {
        Info<<"No firstOderNurbs found"<<endl;
        foundFlag = false;
        return -1;
    }
    //Info<<"fOnN:"<<*firstOrderNearNurbs<<" ";
    DynamicList<scalar> distToNurbsSurface;
    DynamicList<scalar> paraToNurbsSurface;
    DynamicList<label> indToNurbsSurface;
    bool allOutSideNurbsBox = true;
    for(int k=0;k<firstOrderNearNurbs->size();k++)
    {
        label thisNurbs = (*firstOrderNearNurbs)[k];
        scalar thisNodePara = NurbsTrees[thisNurbs]->closestParaOnNurbsToPoint(pnt);
        //Info<<"\tIndex of nurbs:"<<thisNurbs<<" with para: "<<thisNodePara<<endl;
        //Info<<"(*(this->Curves))[thisNurbs].min_U():"<<(*(this->Curves))[thisNurbs].min_U()<<endl;
        if(thisNodePara < (*(this->Curves))[thisNurbs].min_U())
        {
            dist = std::numeric_limits<scalar>::max();
            continue;
        }
        allOutSideNurbsBox = false;
        paraToNurbsSurface.append(thisNodePara);
        distToNurbsSurface.append((*(this->Curves))[thisNurbs].distanceToNurbsSurface(thisNodePara,pnt));
        indToNurbsSurface.append(thisNurbs);
        //Info<<"paraToNurbsSurface:"<<paraToNurbsSurface<<" distToNurbsSurface:"<<distToNurbsSurface<<" indToNurbsSurface:"<<indToNurbsSurface<<endl;
    }
    if(allOutSideNurbsBox)
    {
        foundFlag = false;
        return -1;
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
        dist = minDistToNurbsSurface;
        
        reference.clear();
        nurbsReference temp;
        temp.nurbsInd = minDistindToNurbsSurface;
        temp.nurbsPara = minDistparaToNurbsSurface;
        reference.append(temp);
    }
    else
    {
        if(secondMinDistToNurbsSurface < minDistToNurbsSurface)
            FatalErrorInFunction<<"Second smallest dist smaller than smallest one. Can not happen!"<< exit(FatalError);
        
        vector vecToMinDistNurbs = (*Curves)[minDistindToNurbsSurface].Curve_Derivative(0,minDistparaToNurbsSurface);
        vecToMinDistNurbs = vecToMinDistNurbs - pnt;
        vector vecToSecondMinDistNurbs = (*Curves)[secondMinDistindToNurbsSurface].Curve_Derivative(0,secondMinDistparaToNurbsSurface);
        vecToSecondMinDistNurbs = vecToSecondMinDistNurbs - pnt;
        
        scalar angle;
        bool vecToNurbsZeroOnce = false;
        if(norm2(vecToMinDistNurbs) * norm2(vecToSecondMinDistNurbs) != 0)
            angle  = (vecToMinDistNurbs && vecToSecondMinDistNurbs) / (norm2(vecToMinDistNurbs) * norm2(vecToSecondMinDistNurbs));
        else
        {
            angle = -1;
            vecToNurbsZeroOnce = true;
        }
        scalar radiusFactor = ((angle+1.)/2.);
        scalar smoothingRadius = radiusFactor * intersectionRadius;

        if(vecToNurbsZeroOnce && std::abs(minDistToNurbsSurface)<smoothingRadius && std::abs(secondMinDistToNurbsSurface)<smoothingRadius)
        {
            // rounding gets reducing further away from zero surface
            scalar roundingFactor = 1-(std::abs(minDistindToNurbsSurface)/smoothingRadius);
            
            // rounding gets scaled in respect of similarity of distance measure
            scalar distFirstToSecondMin  = std::abs(minDistToNurbsSurface - secondMinDistToNurbsSurface);
            scalar distFirstToSecondMinFactor = distFirstToSecondMin / (2*smoothingRadius);

            reference.clear();
            nurbsReference temp1;
            temp1.nurbsInd = minDistindToNurbsSurface;
            temp1.nurbsPara = minDistparaToNurbsSurface;
            reference.append(temp1);
            nurbsReference temp2;
            temp2.nurbsInd = secondMinDistindToNurbsSurface;
            temp2.nurbsPara = secondMinDistparaToNurbsSurface;
            reference.append(temp2);
            
            dist = minDistToNurbsSurface + smoothingRadius*(1-distFirstToSecondMinFactor)*roundingFactor;
        }
        else
        {
            reference.clear();
            nurbsReference temp1;
            temp1.nurbsInd = minDistindToNurbsSurface;
            temp1.nurbsPara = minDistparaToNurbsSurface;
            reference.append(temp1);
            
            dist = minDistToNurbsSurface;
        }
    }
    
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
        if(thisNodePara < (*(this->Curves))[thisNurbs].min_U())
        {
            continue;
        }
        allOutSideNurbsBox = false;
        indNurbs.append(thisNurbs);
        paraToNurbsSurface.append(thisNodePara);
        distToNurbsSurface.append((*(this->Curves))[thisNurbs].distanceToNurbsSurface(thisNodePara,pnt));
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

Foam::vector Foam::cutCellFvMesh::vectorToNurbs
(
    point pnt,
    bool& foundFlag
)
{
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
        if(thisNodePara < (*(this->Curves))[thisNurbs].min_U())
        {
            continue;
        }
        allOutSideNurbsBox = false;
        paraToNurbsSurface.append(thisNodePara);
        distToNurbsSurface.append((*(this->Curves))[thisNurbs].distanceToNurbsSurface(thisNodePara,pnt));
        nurbsInd.append(thisNurbs);
    }
    if(allOutSideNurbsBox)
    {
        foundFlag = false;
        return vector();
    }
    label minNurbs = -1;
    scalar minParaToNurbsSurface=-1;
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
    scalar para = minParaToNurbsSurface;
    vector nurbsPoint = (*(this->Curves))[minNurbs].Curve_Derivative(0,para);
    vector pointToNurbsVector = nurbsPoint-pnt;
    return pointToNurbsVector;
}

Foam::List<Foam::vector> Foam::cutCellFvMesh::vectorsToNurbsOfEdge
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

Foam::label Foam::cutCellFvMesh::sideToNurbs(point pnt,bool& foundFlag)
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

/*
void Foam::cutCellFvMesh::computeClosedFaceFront
(
    label centerPointInd,
    label faceInd,
    labelList& problematicFacePoints,
    const labelList& cellsInd,
    List<DynamicList<face>>& facesInCellsIn,
    List<DynamicList<std::unordered_set<label>>>& facesMapInCellsIn,
    DynamicList<bool>& closedFaceFrontOutFaceInFace,
    DynamicList<DynamicList<face>>& closedFaceFrontOut,
    DynamicList<DynamicList<std::unordered_set<label>>>& closedFaceFrontMapOut
)
{
    const pointField& basisPoints = this->points();
    const faceList& basisFaces = this->faces();
    const cellList& basisCells = this->cells();
    const edgeList& basisEdges = this->edges();
    Info<<"---------------------------In Func-------------------"<<endl;
    Info<<"centerPointInd:"<<centerPointInd<<endl;
    Info<<"facesInCellsIn:"<<facesInCellsIn<<endl;
    Info<<"cellsInd:"<<cellsInd<<endl;
    
    if(facesInCellsIn.size()!=facesMapInCellsIn.size())
        FatalErrorInFunction<<"Unequal range of parameters!"<< exit(FatalError);
    
    //Remove faces not touching the centerPointInd
    List<DynamicList<face>> facesOfCells(facesInCellsIn.size());
    List<DynamicList<std::unordered_set<label>>> facesOfCellsMap(facesInCellsIn.size());
    List<DynamicList<bool>> inFaceFaceOfCells(facesInCellsIn.size());
    List<DynamicList<DynamicList<std::pair<label,label>>>> equalFaces(facesInCellsIn.size());
    for(int i=0;i<facesInCellsIn.size();i++)
    {
        Info<<"--------------i:"<<i<<"-----------"<<endl;
        DynamicList<label> faceInFace;
        DynamicList<label> faceInCell;
        DynamicList<label> faceMixed;
        label removedFaces=0;
        cell currCell = basisCells[cellsInd[i]];
        for(int j=0;j<facesInCellsIn[i].size();j++)
        {
            Info<<"--------------j:"<<j<<"-----------:"<<facesInCellsIn[i][j]<<endl;
            Info<<"nbrOfPrevPoints:"<<nbrOfPrevPoints<<endl;
            //test if face contains center point
            bool faceContainsCenterPnt=false;
            if(facesMapInCellsIn[i][j].count(centerPointInd)!=0)
                faceContainsCenterPnt=true;
            else
                removedFaces++;
            
            bool isFaceInFace = false;
            bool isFaceMixed = false;
            bool isFaceInCell = false;
            
            //find all used points and count used points in face
            labelList nbrPntsInFace(currCell.size(),0);
            label maxNbrPntsInFace = 0;
            for(int k=0;k<facesInCellsIn[i][j].size();k++)
            {
                label cutFacePnt = facesInCellsIn[i][j][k];                    
                if(cutFacePnt<nbrOfPrevPoints)
                {
                    for(int l=0;l<currCell.size();l++)
                    {
                        face thisFace = basisFaces[currCell[l]];
                        if(thisFace.which(cutFacePnt)!=-1)
                        {
                            nbrPntsInFace[l]++;
                            Info<<cutFacePnt<<" is in"<<endl;
                        }
                    }
                }
                else
                {
                    label pntEdgeInd = pointToEgde_[cutFacePnt];
                    if(pntEdgeInd==-1)
                        FatalErrorInFunction<<"Can not happen!"<<exit(FatalError);
                    edge pntEdge = newMeshEdges_[pntEdgeInd];
                    Info<<"pntEdge:"<<pntEdge<<endl;
                    for(int l=0;l<currCell.size();l++)
                    {
                        face thisFace = basisFaces[currCell[l]];
                        if(thisFace.which(pntEdge.start())!=-1 && thisFace.which(pntEdge.end())!=-1)
                        {
                            nbrPntsInFace[l]++;
                            Info<<cutFacePnt<<" is in"<<endl;
                        }
                    }
                }
            }
            for(int l=0;l<currCell.size();l++)
            {
                if(maxNbrPntsInFace<nbrPntsInFace[l])
                    maxNbrPntsInFace = nbrPntsInFace[l];
            }
            
            //evaluate situation 
            if(problematicFacePoints[faceInd]==4)
            {
                if(maxNbrPntsInFace==2 || maxNbrPntsInFace==0)
                    isFaceInCell=true;
                else if(maxNbrPntsInFace==3)
                    FatalErrorInFunction<<"Can not happen"<< exit(FatalError);
                else if(maxNbrPntsInFace==4 && facesInCellsIn[i][j].size()==4)
                    isFaceInFace=true;
                else if(maxNbrPntsInFace==4 && facesInCellsIn[i][j].size()>4)
                    isFaceMixed=true;
                else
                    FatalErrorInFunction<<"Problematic four point face must have 3, 2 or 0 points in Face"<< exit(FatalError);
            }
            else if(problematicFacePoints[faceInd]==3)
            {
                if(maxNbrPntsInFace==2 || maxNbrPntsInFace==1 || maxNbrPntsInFace==0)
                    isFaceInCell=true;
                else if(maxNbrPntsInFace==3 && facesInCellsIn[i][j].size()==3)
                    isFaceInFace=true;
                else if(maxNbrPntsInFace==3 && facesInCellsIn[i][j].size()>3)
                    isFaceMixed=true;
                else
                {
                    Info<<"basisFaces["<<faceInd<<"]:"<<basisFaces[faceInd]<<endl;
                    Info<<"i:"<<i<<"  j:"<<j<<endl;
                    Info<<"maxNbrPntsInFace:"<<maxNbrPntsInFace<<endl;
                    Info<<"facesInCellsIn[i][j]:"<<facesInCellsIn[i][j]<<endl;
                    FatalErrorInFunction<<"Problematic four point face must have 3, 2 or 0 points in Face"<< exit(FatalError);      
                }
            }                
            else if(problematicFacePoints[faceInd]==2)
            {
                isFaceInCell=true;
            }
            else
            {
                Info<<"problematicFacePoints["<<faceInd<<"]:"<<problematicFacePoints[faceInd]<<endl;
                FatalErrorInFunction<<"ProblematicFacePoints must be 4,3 or 2!"<< exit(FatalError);
            }
            if(faceContainsCenterPnt)
            {
                if(isFaceInFace && !isFaceMixed && !isFaceInCell)
                    faceInFace.append(j);
                else if(!isFaceInFace && isFaceMixed && !isFaceInCell)
                    faceMixed.append(j);
                else if(!isFaceInFace && !isFaceMixed && isFaceInCell)
                    faceInCell.append(j);
                else
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
            }
            Info<<"nbrPntsInFace:"<<nbrPntsInFace<<endl;
            Info<<"isFaceInFace:"<<isFaceInFace<<endl;
            Info<<"isFaceMixed:"<<isFaceMixed<<endl;
            Info<<"isFaceInCell:"<<isFaceInCell<<endl;
            Info<<"faceToRemove:"<<!faceContainsCenterPnt<<endl;
        }
        
        if((removedFaces+faceInFace.size()+faceInCell.size()+faceMixed.size())!=facesInCellsIn[i].size())
        {
            Info<<"removedFaces:"<<removedFaces<<endl;
            Info<<"faceInFace:"<<faceInFace<<endl;
            Info<<"faceInCell:"<<faceInCell<<endl;
            Info<<"faceMixed:"<<faceMixed<<endl;
            Info<<"facesInCellsIn["<<i<<"]:"<<facesInCellsIn[i]<<endl;
            FatalErrorInFunction<<"Faces lost! Can not happen."<< exit(FatalError);
        }
        if((faceInFace.size()==0 && (facesInCellsIn[i].size()-removedFaces)>1))
        {
            Info<<endl<<endl;
            const cellList& cells = this->cells();
            const faceList& faces = this->faces();
            const pointField& points = this->points();
            for(const label& cellI: cellsInd)
            {
                const cell& oneCell = cells[cellI];
                Info<<"cellInd:"<<cellI<<endl;
                Info<<"cellLabels:"<<oneCell.labels(faces)<<endl;
                Info<<"cellPoints:"<<oneCell.points(faces,points)<<endl;
                Info<<"---"<<endl;
            }
            Info<<"faceInd:"<<faceInd<<endl;
            Info<<"face:"<<faces[faceInd]<<endl;
            Info<<"problematicFacePoints[faceInd]:"<<problematicFacePoints[faceInd]<<endl;
            Info<<"centerPointInd:"<<centerPointInd<<endl;
            Info<<"faceInFace:"<<faceInFace<<endl;
            Info<<"faceInCell:"<<faceInCell<<endl;
            Info<<"faceMixed:"<<faceMixed<<endl;
            Info<<"removedFaces:"<<removedFaces<<endl;
            Info<<"facesInCellsIn["<<i<<"]:"<<facesInCellsIn[i]<<endl;
            for(const labelList& face : facesInCellsIn[i])
            {
                for(const label& pntInd : face)
                {
                    Info<<pntInd<<" -> "<<newMeshPoints_[pntInd]<<endl;
                }
            }
            FatalErrorInFunction<<"One face in face must exist for multi face cells that call this function!"<< exit(FatalError);
        }
        if(faceMixed.size()>0 && (faceInFace.size()==0 || faceInCell.size()==0) && removedFaces==0)
            FatalErrorInFunction<<"Mixed face but not faceInFace and FaceInCell! Can not happen."<< exit(FatalError);
        
        Info<<"removedFaces:"<<removedFaces<<endl;
        Info<<"faceInFace:"<<faceInFace<<endl;
        Info<<"faceInCell:"<<faceInCell<<endl;
        Info<<"faceMixed:"<<faceMixed<<endl;
        for(int j=0;j<faceInFace.size();j++)
        {
            facesOfCells[i].append(facesInCellsIn[i][faceInFace[j]]);
            facesOfCellsMap[i].append(facesMapInCellsIn[i][faceInFace[j]]);
            inFaceFaceOfCells[i].append(true);
        }
        for(int j=0;j<faceInCell.size();j++)
        {
            facesOfCells[i].append(facesInCellsIn[i][faceInCell[j]]);
            facesOfCellsMap[i].append(facesMapInCellsIn[i][faceInCell[j]]);
            inFaceFaceOfCells[i].append(false);
        }
    }
    DynamicList<std::pair<label,label>> faceInFaceIndx;
    DynamicList<face> faceInFaceFace;
    for(int i=0;i<facesOfCells.size();i++)
    {
        for(int j=0;j<facesOfCells[i].size();j++)
        {
            if(inFaceFaceOfCells[i][j])
            {
                faceInFaceIndx.append(std::pair<label,label>(i,j));
                faceInFaceFace.append(facesOfCells[i][j]);
            }
        }
    }
    for(int i=0;i<facesOfCells.size();i++)
    {
        for(int j=0;j<facesOfCells[i].size();j++)
        {
            equalFaces[i].append(DynamicList<std::pair<label,label>>());
        }
    }
    for(int i=0;i<faceInFaceIndx.size();i++)
    {
        face oneFaceInFace = faceInFaceFace[i];
        std::unordered_set<label> oneFaceInFaceMap = facesOfCellsMap[faceInFaceIndx[i].first][faceInFaceIndx[i].second];
        for(int j=0;j<faceInFaceIndx.size();j++)
        {
            if(i==j)
                continue;
            
            face compFaceInFace = faceInFaceFace[j];
            std::unordered_set<label> compFaceInFaceMap = facesOfCellsMap[faceInFaceIndx[j].first][faceInFaceIndx[j].second];
            if(oneFaceInFace.size()==compFaceInFace.size())
            {
                bool isEqualFace = true;
                for(int k=0;k<oneFaceInFace.size();k++)
                {
                    if(compFaceInFaceMap.count(oneFaceInFace[k])==0)
                        isEqualFace=false;
                }
                if(isEqualFace)
                {
                    equalFaces[faceInFaceIndx[i].first][faceInFaceIndx[i].second].append(faceInFaceIndx[j]);
                }
            }
        }
    }
    
    Info<<"centerPointInd:"<<centerPointInd<<endl;
    Info<<"facesOfCells:"<<facesOfCells<<endl;
    for(int i=0;i<equalFaces.size();i++)
    {
        Info<<"Cell  [";
        for(int j=0;j<equalFaces[i].size();j++)
        {
            Info<<"{ face: ";
            for(int k=0;k<equalFaces[i][j].size();k++)
            {
                Info<<" ("<<equalFaces[i][j][k].first<<"|"<<equalFaces[i][j][k].second<<") ";
            }
            Info<<" }";
        }
        Info<<"]"<<endl;
    }
    
    DynamicList<DynamicList<std::pair<label,label>>> uniqueCycles;
    DynamicList<std::unordered_multimap<label,label>> uniqueCyclesMap; 
    DynamicList<bool> uniqueCyclesFaceInFace;
    //Combine facesInCells to faceFronts
    enum fcState {OPEN=0,CLOSED=1,FAIL=-1};    
    for(int i=0;i<facesOfCells.size();i++)
    {
        for(int j=0;j<facesOfCells[i].size();j++)
        {            
            std::pair<label,label> startFace(i,j);
            //Info<<"start: "<<"("<<i<<"|"<<j<<")---------------------------------------"<<endl;
            DynamicList<std::unordered_multimap<label,label>> usedFaces;
            usedFaces.setSize(1);
            usedFaces[0].insert(startFace);
            DynamicList<DynamicList<std::pair<label,label>>> faceCycle;
            faceCycle.setSize(1);
            faceCycle[0].append(startFace);
            DynamicList<fcState> state;
            state.setSize(1);
            state[0]=OPEN;
            DynamicList<DynamicList<std::pair<label,label>>> lastFaceEnterEdgeList;
            lastFaceEnterEdgeList.setSize(1);
            bool allDone = false;
            label loopCounter = 0;
            
            while(!allDone)
            {
                //Info<<"Loop:"<<loopCounter<<endl;
                
                label len = faceCycle.size();
                for(int k=0;k<len;k++)
                //Iterate over all cycles
                {                    
                    if(state[k] != OPEN)
                        continue;
                    
                    std::pair<label,label> currFace = faceCycle[k].last();
                    
                    //Info<<"Test is closed for k:"<<k<<endl;
                    // Test if cycle is closed
                    DynamicList<label> equalPntsList;
                    label equalPnts = 0;
                    for(int n=0;n<facesOfCells[startFace.first][startFace.second].size();n++)
                    {
                        if(facesOfCellsMap[currFace.first][currFace.second].count(facesOfCells[startFace.first][startFace.second][n])!=0)
                        {
                            equalPnts++;
                            equalPntsList.append(facesOfCells[startFace.first][startFace.second][n]);
                        }
                    }
                    if(faceCycle[k].size()>2 && equalPnts==2)
                    {                               
                        bool isEntryEdge = false;
                        std::pair<label,label> enterEdge = lastFaceEnterEdgeList[k].first();
                        if((enterEdge.first == equalPntsList[0] && enterEdge.second == equalPntsList[1]) ||
                           (enterEdge.first == equalPntsList[1] && enterEdge.second == equalPntsList[0]))
                            isEntryEdge = true;
                            
                        if(!isEntryEdge)
                        {                            
                            state[k] = CLOSED;
                            continue;
                        }
                    }
                    else if((faceCycle[k].size()==1 && equalPnts!=facesOfCells[faceCycle[k][0].first][faceCycle[k][0].second].size()) || (faceCycle[k].size()==2 && equalPnts!=2))
                    {
                        Info<<"faceCycle[k].size():"<<faceCycle[k].size()<<endl;
                        Info<<"equalPnts:"<<equalPnts<<endl;
                        FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                    }
                    
                    //Info<<"Find next faces for k:"<<k<<endl;
                    //Find the next faces
                    DynamicList<std::pair<label,label>> addedFaces;
                    for(int l=0;l<facesOfCells.size();l++)
                    {
                        auto facesRange = usedFaces[k].equal_range(l);
                        std::unordered_set<label> usedFacesInThisCell;
                        //Info<<"\t\t\t\t already taken in cell:"<<l;
                        for(auto it=facesRange.first; it!=facesRange.second; ++it)
                        {
                            //Info<<" face:"<<it->second;
                            usedFacesInThisCell.insert(it->second);
                        }
                        //Info<<endl;
                        
                        for(int m=0;m<facesOfCells[l].size();m++)
                        {
                            if(usedFacesInThisCell.count(m)!=0)
                                continue;
                            
                            //Info<<"("<<l<<"|"<<m<<")"<<endl;
                            DynamicList<label> equalPntsList;
                            label equalPnts = 0;
                            for(int n=0;n<facesOfCells[l][m].size();n++)
                            {
                                if(facesOfCellsMap[currFace.first][currFace.second].count(facesOfCells[l][m][n])!=0)
                                {
                                    equalPnts++;
                                    equalPntsList.append(facesOfCells[l][m][n]);
                                }
                            }
                            if(equalPnts==2)
                            {
                                if(equalPntsList.size()!=2)
                                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                                
                                bool isEntryEdge = false; 
                                if(lastFaceEnterEdgeList[k].size()>0)
                                {
                                    std::pair<label,label> enterEdge = lastFaceEnterEdgeList[k].last();
                                    if((enterEdge.first == equalPntsList[0] && enterEdge.second == equalPntsList[1]) ||
                                       (enterEdge.first == equalPntsList[1] && enterEdge.second == equalPntsList[0]))
                                        isEntryEdge = true;
                                }
                                    
                                if(!isEntryEdge)
                                {
                                    addedFaces.append(std::pair<label,label>(l,m));
                                    //Info<<"\t\tAdd "<<"cell:"<<l<<" face:"<<m<<endl;
                                }
                            }
                            else if(equalPnts!=1 && (facesOfCells[currFace.first][currFace.second].size()!=facesOfCells[l][m].size()))
                            {
                                Info<<"currFace: ("<<currFace.first<<"|"<<currFace.second<<")"<<endl;
                                Info<<"face: ("<<l<<"|"<<m<<") (";
                                for(int n=0;n<facesOfCells[l][m].size();n++)
                                {
                                    if(facesOfCellsMap[currFace.first][currFace.second].count(facesOfCells[l][m][n])!=0)
                                        Info<<" 1";
                                }
                                Info<<" )"<<endl;
                                FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                            }
                        }
                    }
                    
                    //Info<<"Find edges for k:"<<k<<endl;
                    //Find the edge points that connect 
                    DynamicList<std::pair<label,label>> addedFacesConnectEdge;
                    for(int n=0;n<addedFaces.size();n++)
                    {
                        face faceCurr = facesOfCells[currFace.first][currFace.second];
                        face faceAdd = facesOfCells[addedFaces[n].first][addedFaces[n].second];
                        DynamicList<label> connectingPnts;
                        for(int o=0;o<faceAdd.size();o++)
                        {
                            label faceInd;
                            if((faceInd=faceCurr.which(faceAdd[o]))!=-1)
                            {
                                connectingPnts.append(faceAdd[o]);
                                if(faceAdd[o]!=faceCurr[faceInd])
                                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                            }
                        }
                        if(connectingPnts.size()!=2)
                            FatalErrorInFunction<<"Connecting edge must have two points!"<< exit(FatalError);
                        addedFacesConnectEdge.append(std::pair<label,label>(connectingPnts[0],connectingPnts[1]));
                    }
                    
                    //Info<<"Add faces for k:"<<k<<endl;
                    //Add the next faces and duplicate if necessary
                    if(addedFaces.size()==0)
                    {
                        state[k] = FAIL;
                    }
                    else
                    {
                        //Create copy of datastructures
                        DynamicList<std::pair<label,label>> cpCycle = faceCycle[k];
                        fcState cpState = state[k];
                        std::unordered_multimap<label,label> cpUsedFaces = usedFaces[k];
                        DynamicList<std::pair<label,label>> cpLastFaceEnterEdgeList = lastFaceEnterEdgeList[k];

                        //append next items to datastructures
                        faceCycle[k].append(addedFaces[0]);                        
                        lastFaceEnterEdgeList[k].append(addedFacesConnectEdge[0]);
                        usedFaces[k].insert(addedFaces[0]);
                        //Info<<"-n blocked: ("<<addedFaces[0].first<<"|"<<addedFaces[0].second<<") ";
                        for(int n=0;n<equalFaces[addedFaces[0].first][addedFaces[0].second].size();n++)
                        {
                            usedFaces[k].insert(equalFaces[addedFaces[0].first][addedFaces[0].second][n]);
                            //Info<<" eq("<<equalFaces[addedFaces[0].first][addedFaces[0].second][n].first<<"|"<<equalFaces[addedFaces[0].first][addedFaces[0].second][n].second<<") ";
                        }
                        //Info<<endl;
                        
                        for(int o=1;o<addedFaces.size();o++)
                        {
                            faceCycle.append(cpCycle);
                            faceCycle.last().append(addedFaces[o]);
                            
                            state.append(cpState);
                            
                            usedFaces.append(cpUsedFaces);
                            usedFaces.last().insert(addedFaces[o]);
                            DynamicList<std::pair<label,label>> equalFaceList = equalFaces[addedFaces[o].first][addedFaces[o].second];
                            //Info<<"-n-o blocked: ("<<addedFaces[o].first<<"|"<<addedFaces[o].second<<") ";
                            for(int n=0;n<equalFaceList.size();n++)
                            {
                                usedFaces.last().insert(equalFaceList[n]);
                                //Info<<" eq("<<equalFaceList[n].first<<"|"<<equalFaceList[n].second<<") ";
                            }
                            //Info<<endl;
                            
                            lastFaceEnterEdgeList.append(cpLastFaceEnterEdgeList);
                            lastFaceEnterEdgeList.last().append(addedFacesConnectEdge[o]);
                        }
                    }
                    //Info<<"Done for k:"<<k<<endl;
                }
                bool noOpen=true;
                for(int k=0;k<state.size();k++)
                {
                    if(state[k]==OPEN)
                        noOpen=false;
                }
                if(noOpen)
                    allDone=true;
                Info<<endl;
                
                loopCounter++;
                if(loopCounter>32)
                {
                    Info<<"centerPointInd:"<<centerPointInd<<endl;
                    Info<<"facesOfCells:"<<facesOfCells<<endl;
                    Info<<"\t\t"<<"faceCycle"<<endl;
                    for(int k=0;k<faceCycle.size();k++)
                    {
                        Info<<"\t\t\t[";
                        for(int l=0;l<faceCycle[k].size();l++)
                        {
                            Info<<" ("<<faceCycle[k][l].first<<"|"<<faceCycle[k][l].second<<")";
                        }
                        Info<<" ]";
                        Info<<"   "<<"fcstate:"<<state[k];
                        Info<<endl;
                    }
                    FatalErrorInFunction<<"Loop iterates too often possibly infinite!"<< exit(FatalError);
                }
            }
            for(int k=0;k<faceCycle.size();k++)
            //iterate across all found closed face cycles
            {
                //Info<<"Insert in uniqueFace Cylce"<<endl;
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
                    std::unordered_multimap<label,label> oneCycleMap = uniqueCyclesMap[l];
                    for(int m=0;m<oneCylce.size();m++)
                    //test one face against one cycleMap    
                    {
                        std::pair<label,label> oneFace = oneCylce[m];
                        auto itemRange = oneCycleMap.equal_range(oneFace.first);
                        if(itemRange.first==itemRange.second)
                        //cell of face does not exist in map
                        {
                            allFaceMatch = false;
                        }
                        else
                        {
                            bool sameFaceInCellExists = false;
                            for(auto item=itemRange.first;item!=itemRange.second;++item)
                            {
                                if(item->second==oneFace.second)
                                //face does exist in map
                                {    
                                    sameFaceInCellExists = true;
                                }
                            }
                            if(!sameFaceInCellExists)
                            {
                                allFaceMatch = false;
                            }   
                        }
                    }
                    if(allFaceMatch)
                        noMatch = false;
                    //Info<<"Insert in uniqueFace Cylce DONE"<<endl;
                }
                if(noMatch)
                {
                    uniqueCycles.append(faceCycle[k]);
                    std::unordered_multimap<label,label> cycleMap;
                    for(int l=0;l<faceCycle[k].size();l++)
                    {
                        cycleMap.insert(faceCycle[k][l]);
                        DynamicList<std::pair<label,label>>& equalFacesToBlock =  equalFaces[faceCycle[k][l].first][faceCycle[k][l].second];
                        for(int m=0;m<equalFacesToBlock.size();m++)
                        {
                            cycleMap.insert(equalFacesToBlock[m]);
                        }
                    }
                    uniqueCyclesMap.append(cycleMap);
                    
                    bool oneFaceInFace = false;
                    for(int l=0;l<faceCycle[k].size();l++)
                    {
                        std::pair<label,label> thisFace = faceCycle[k][l];
                        oneFaceInFace |= inFaceFaceOfCells[thisFace.first][thisFace.second];
                    }
                    uniqueCyclesFaceInFace.append(oneFaceInFace);
                    
                    Info<<"Append Face  [";
                    for(int j=0;j<uniqueCycles.last().size();j++)
                    {
                        Info<<" ("<<uniqueCycles.last()[j].first<<"|"<<uniqueCycles.last()[j].second<<") ";
                    }
                    Info<<"]"<<endl;
                    Info<<"Map:";
                    for(auto item :uniqueCyclesMap.last())
                    {
                        Info<<" ("<<item.first<<"|"<<item.second<<")  ";
                    }
                    Info<<"End Map"<<endl;
                    Info<<endl;
                }
            }
        }
    }


    if(uniqueCycles.size()==0)
    {    
        Info<<"centerPointInd:"<<centerPointInd<<endl;
        Info<<"facesInCellsIn:"<<facesInCellsIn<<endl;
        Info<<"facesOfCells:"<<facesOfCells<<endl;
        Info<<"cellsInd:"<<cellsInd<<endl;
        FatalErrorInFunction<<"There must be at least one face front!"<< exit(FatalError);
    }
    for(int i=0;i<uniqueCycles.size();i++)
    {
        
        if(uniqueCycles[i].size()<3 || uniqueCycles[i].size()>8)
        {
            for(int j=0;j<cellsInd.size();j++)
            {
                cell oneCell = basisCells[cellsInd[j]];
                Info<<"cellInd:"<<cellsInd[j]<<" "<<oneCell.labels(basisFaces)<<endl<<oneCell.points(basisFaces,basisPoints)<<endl;
            }
            Info<<"cellsInd:"<<cellsInd<<endl;
            Info<<"facesInCellsIn:"<<facesInCellsIn<<endl;
            std::unordered_set<label> pntsLabel;
            for(int j=0;j<facesInCellsIn.size();j++)
            {
                for(int k=0;k<facesInCellsIn[j].size();k++)
                {
                    for(int l=0;l<facesInCellsIn[j][k].size();l++)
                    {
                        pntsLabel.insert(facesInCellsIn[j][k][l]);
                    }
                }
            }
            for(const label& pnt: pntsLabel)
            {
                if(pointToEgde_[pnt]!=-1)
                    Info<<pnt<<":"<<basisEdges[pointToEgde_[pnt]]<<endl;
            }
            for(int j=0;j<uniqueCycles.size();j++)
            {
                Info<<j<<" ";
                for(int k=0;k<uniqueCycles[j].size();k++)
                {
                    Info<<"("<<uniqueCycles[j][k].first<<" "<<uniqueCycles[j][k].second<<") ";
                }
                Info<<endl;
            }
            Info<<"i:"<<i<<endl;
            Info<<"facesOfCells:"<<facesOfCells<<endl;
            FatalErrorInFunction<<"There must be between 3 and 8 faces in a front!"<< exit(FatalError);
        }
        
    }
    Info<<"facesOfCells:"<<facesOfCells<<endl;
    for(int i=0;i<uniqueCycles.size();i++)
    {
        Info<<"Face "<<i<<"  [";
        for(int j=0;j<uniqueCycles[i].size();j++)
        {
            Info<<" ("<<uniqueCycles[i][j].first<<"|"<<uniqueCycles[i][j].second<<") ";
        }
        Info<<"]"<<endl;
        Info<<"Map:";
        for(auto item :uniqueCyclesMap[i])
        {
            Info<<" ("<<item.first<<"|"<<item.second<<")  ";
        }
        Info<<"End Map"<<endl;
        Info<<endl;
    }

    for(int i=0;i<uniqueCycles.size();i++)
    {
        closedFaceFrontOut.append(DynamicList<face>());
        closedFaceFrontMapOut.append(DynamicList<std::unordered_set<label>>());
        closedFaceFrontOutFaceInFace.append(uniqueCyclesFaceInFace[i]);
        for(int j=0;j<uniqueCycles[i].size();j++)
        {
            std::pair<label,label> face = uniqueCycles[i][j];
            closedFaceFrontOut.last().append(facesOfCells[face.first][face.second]);
            closedFaceFrontMapOut.last().append(facesOfCellsMap[face.first][face.second]);
        }
    }
    
    for(int i=0;i<closedFaceFrontOut.size();i++)
    {
        Info<<"nbrOfPrevPoints:"<<nbrOfPrevPoints<<endl;
        Info<<"Face "<<i<<"  [";
        for(int j=0;j<closedFaceFrontOut[i].size();j++)
        {
            Info<<closedFaceFrontOut[i][j]<<" ";
        }
        Info<<"]  faceInFace:"<<closedFaceFrontOutFaceInFace[i]<<endl;
        Info<<endl;
        
        if(closedFaceFrontOut[i].size()==3)
            FatalErrorInFunction<<"One face cycle around point consists of three faces! Can not happen! "<< exit(FatalError);
    }
    if(closedFaceFrontOutFaceInFace.size()!=closedFaceFrontOut.size())
        FatalErrorInFunction<<"Can not happen"<< exit(FatalError);
    if(closedFaceFrontOutFaceInFace.size()!=closedFaceFrontMapOut.size())
        FatalErrorInFunction<<"Can not happen"<< exit(FatalError);
}
*/

/*
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
*/

/*
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
    
    Foam::vector innerVector = connectVertex-otherVertexInner;
    Foam::vector outerVector = otherVertexOuter-connectVertex;
    
    scalar normInnerVector = norm2(outerVector);
    scalar normOuterVector = norm2(outerVector);
    scalar scalProd = innerVector & outerVector;
    
    scalar pseudoAngle = scalProd/(normInnerVector*normOuterVector);
    if(pseudoAngle>1 || pseudoAngle<-1)
        FatalErrorInFunction<<"Pseudo angle must be inside [-1,1]! Can not happen!"<< exit(FatalError);
    return pseudoAngle;
}
*/

void Foam::cutCellFvMesh::newMeshEdges_MC33
(
)
{
    const faceList& basisFaces = this->faces();
    const edgeList& basisEdges = this->edges();
    const labelList& faceOwner = this->faceOwner();
    const labelList& faceNeighbor = this->faceNeighbour();
    const labelListList& faceToEdge = this->faceEdges();
    const labelListList& pointToEdge = this->pointEdges();
    
    nbrOfPrevEdges = basisEdges.size();
    newMeshEdges_.append(basisEdges);
    edgesToSide(newMeshEdges_);
        
    faceToEdges_.setSize(basisFaces.size());
    edgeToFaces_.setSize(basisEdges.size());
    const labelListList& edgeFaces = this->edgeFaces();

    labelList thisFacePoints;
    for(int faceInd=0;faceInd<basisFaces.size();faceInd++)
    {
        const labelList thisFaceEdges = faceToEdge[faceInd];
        std::unordered_map<label,label> edgesOfFace;
        for(int j=0;j<thisFaceEdges.size();j++)
        {
            edgesOfFace.insert({thisFaceEdges[j],j});
        }
        std::vector<label> cellsOfThisFace = {faceOwner[faceInd]};
        if(faceInd < faceNeighbor.size())
            cellsOfThisFace.push_back(faceNeighbor[faceInd]);
        
        struct pairHash {
            std::size_t operator()(const std::pair<label,label>& pr) const{
                return std::hash<label>{}(pr.first)^std::hash<label>{}(pr.second);
            }
        };
        struct pairEqual {
            bool operator()(const std::pair<label,label>& a,const std::pair<label,label>& b)const{
                return (std::equal_to<label>{}(a.first,b.first) &&
                        std::equal_to<label>{}(a.second,b.second))||
                       (std::equal_to<label>{}(a.first,b.second) && 
                        std::equal_to<label>{}(a.second,b.first));
            }
        };
        std::vector<std::unordered_map<std::pair<label,label>,std::pair<label,label>,pairHash,pairEqual>>addedEdgesFromSide(cellsOfThisFace.size());
        std::vector<bool> areThereCutsInFace(cellsOfThisFace.size(),false);
        
        for(unsigned int j=0;j<cellsOfThisFace.size();j++)
        {
            label cellInd = cellsOfThisFace[j];
            MC33::MC33Cube& thisCellCube = mc33CutCellData[cellInd];
            if(thisCellCube.cell!=-1)
            {
                if(cellInd != thisCellCube.cell)
                    FatalErrorInFunction<<"More than two edges!"<< exit(FatalError);
                for(unsigned int k=0;k<thisCellCube.cutTriangles.size();k++)
                {
                    label localFaceEdgeInd1 = std::get<0>(mc33CutCellData[cellInd].cutTriangles[k]);
                    label faceCutEdge1 = mc33CutCellData[cellInd].edgeGlobalInd[localFaceEdgeInd1];
                    label localFaceEdgeInd2 = std::get<1>(mc33CutCellData[cellInd].cutTriangles[k]);
                    label faceCutEdge2 = mc33CutCellData[cellInd].edgeGlobalInd[localFaceEdgeInd2];
                    label localFaceEdgeInd3 = std::get<2>(mc33CutCellData[cellInd].cutTriangles[k]);
                    label faceCutEdge3 = mc33CutCellData[cellInd].edgeGlobalInd[localFaceEdgeInd3];
                    
                    bool foundCutFaceEdges=false;
                    label loc_edg_Ind1 = -1;
                    label loc_edg_Ind2 = -1;
                    label which_Triangle_Vertice = -1;
                    if(edgesOfFace.find(faceCutEdge1)!=edgesOfFace.end() &&
                       edgesOfFace.find(faceCutEdge2)!=edgesOfFace.end())
                    {
                        loc_edg_Ind1 = edgesOfFace.find(faceCutEdge1)->second;
                        loc_edg_Ind2 = edgesOfFace.find(faceCutEdge2)->second;
                        if(foundCutFaceEdges)
                            FatalErrorInFunction<<"More than two edges!"<< exit(FatalError);
                        foundCutFaceEdges=true;
                        which_Triangle_Vertice=0;
                    }
                    if(edgesOfFace.find(faceCutEdge2)!=edgesOfFace.end() &&
                       edgesOfFace.find(faceCutEdge3)!=edgesOfFace.end())
                    {
                        loc_edg_Ind1 = edgesOfFace.find(faceCutEdge2)->second;
                        loc_edg_Ind2 = edgesOfFace.find(faceCutEdge3)->second;
                        if(foundCutFaceEdges)
                            FatalErrorInFunction<<"More than two edges!"<< exit(FatalError);
                        foundCutFaceEdges=true;
                        which_Triangle_Vertice=1;
                    }
                    if(edgesOfFace.find(faceCutEdge3)!=edgesOfFace.end() &&
                       edgesOfFace.find(faceCutEdge1)!=edgesOfFace.end())
                    {
                        loc_edg_Ind1 = edgesOfFace.find(faceCutEdge3)->second;
                        loc_edg_Ind2 = edgesOfFace.find(faceCutEdge1)->second;
                        if(foundCutFaceEdges)
                            FatalErrorInFunction<<"More than two edges!"<< exit(FatalError);
                        foundCutFaceEdges=true;
                        which_Triangle_Vertice=2;
                    }
                    if(foundCutFaceEdges)
                    {
                        std::pair<std::pair<label,label>,std::pair<label,label>> keyValue =       
                                    {std::pair<label,label>(loc_edg_Ind1,loc_edg_Ind2),
                                     std::pair<label,label>(k,which_Triangle_Vertice)};
                        addedEdgesFromSide[j].insert(keyValue);
                                                      
                        areThereCutsInFace[j] = true;
                    }
                }
            }
        }
        if(cellsOfThisFace.size()==2)
        {
            for(auto cutEdgeIter0=addedEdgesFromSide[0].cbegin();
                cutEdgeIter0!=addedEdgesFromSide[0].cend();
                cutEdgeIter0++)
            {
                if(addedEdgesFromSide[1].find(cutEdgeIter0->first)==addedEdgesFromSide[1].end())
                    FatalErrorInFunction<<"Incompatible edges!"<< exit(FatalError);
            }
            for(auto cutEdgeIter1=addedEdgesFromSide[1].cbegin();
                cutEdgeIter1!=addedEdgesFromSide[1].cend();
                cutEdgeIter1++)
            {
                if(addedEdgesFromSide[0].find(cutEdgeIter1->first)==addedEdgesFromSide[0].end())
                    FatalErrorInFunction<<"Incompatible edges!"<< exit(FatalError);
            }
            if(addedEdgesFromSide[0].size()!=addedEdgesFromSide[1].size())
            {
                //Info<<"cellsOfThisFace.size():"<<cellsOfThisFace.size()<<Foam::endl;
                /*
                Info<<"cellsOfThisFace[0]:"<<cellsOfThisFace[0]<<Foam::endl;
                Info<<"cellsOfThisFace[1]:"<<cellsOfThisFace[1]<<Foam::endl;
                label cellInd0 = cellsOfThisFace[0];
                MC33::MC33Cube& thisCellCube0 = mc33CutCellData[cellInd0];
                label cellInd1 = cellsOfThisFace[1];
                MC33::MC33Cube& thisCellCube1 = mc33CutCellData[cellInd1];
                Info<<Foam::endl;
                Info<<"faceInd:"<<faceInd<<Foam::endl;
                Info<<"thisCellCube0.cell:"<<thisCellCube0.cell<<Foam::endl;
                Info<<"thisCellCube1.cell:"<<thisCellCube1.cell<<Foam::endl;
                */
                Info<<"addedEdgesFromSide[0].size():"<<addedEdgesFromSide[0].size()<<Foam::endl;
                Info<<"areThereCutsInFace[0]:"<<areThereCutsInFace[0]<<Foam::endl;
                Info<<"addedEdgesFromSide[1].size():"<<addedEdgesFromSide[1].size()<<Foam::endl;
                Info<<"areThereCutsInFace[1]:"<<areThereCutsInFace[1]<<Foam::endl;
                FatalErrorInFunction<<"Incompatible edge number!"<< exit(FatalError);
            }
            if(areThereCutsInFace[0] ^ areThereCutsInFace[1])
            {
                FatalErrorInFunction<<"Incompatible assignments!"<< exit(FatalError);
            }
        }
        else if(cellsOfThisFace.size()==1)
        {}
        else
            FatalErrorInFunction<<"Invalid!"<< exit(FatalError);
        
        if(areThereCutsInFace[0])
        {
            label cellInd = cellsOfThisFace[0];
            MC33::MC33Cube& thisCellCube = mc33CutCellData[cellInd];
            
            for(auto cutEdgeIter0=addedEdgesFromSide[0].cbegin();
                cutEdgeIter0!=addedEdgesFromSide[0].cend();
                cutEdgeIter0++)
            {
                label loc_edgInd1 = cutEdgeIter0->first.first;
                label loc_edgInd2 = cutEdgeIter0->first.second;
                if(loc_edgInd1>=thisFaceEdges.size() || loc_edgInd1<0 || 
                   loc_edgInd2>=thisFaceEdges.size() || loc_edgInd2<0)
                    FatalErrorInFunction<<"Invalid!"<< exit(FatalError);
                
                label triangleNbr = cutEdgeIter0->second.first;
                label verticeNbr = cutEdgeIter0->second.second;
                
                auto& cutTriangle = thisCellCube.cutTriangles[triangleNbr];
                std::vector<label> cutTriangleVec = 
                            {std::get<0>(cutTriangle),std::get<1>(cutTriangle), std::get<2>(cutTriangle)};

                label mc33_loc_edg1 = cutTriangleVec[verticeNbr%3];
                label mc33_loc_edg2 = cutTriangleVec[(verticeNbr+1)%3];
                
                label mc33_edg_vert1 = thisCellCube.cutEdgeVerticeIndex[mc33_loc_edg1];
                label mc33_edg_vert2 = thisCellCube.cutEdgeVerticeIndex[mc33_loc_edg2];
                
                if(mc33_edg_vert1 != mc33_edg_vert2)
                //Avoid null triangles
                {
                    if(mc33_edg_vert1<nbrOfPrevPoints && mc33_edg_vert2<nbrOfPrevPoints)
                    {
                        const labelList& mc33_edg_vert1_edges = pointToEdge[mc33_edg_vert1];
                        const labelList& mc33_edg_vert2_edges = pointToEdge[mc33_edg_vert2];
                        label sharedEdge = -1;
                        for(const label edgeInd1: mc33_edg_vert1_edges)
                        {
                            for(const label edgeInd2: mc33_edg_vert2_edges)
                            {
                                if(edgeInd1==edgeInd2)
                                {
                                    if(sharedEdge!=-1)
                                        FatalErrorInFunction<<"More than one shared edge"<< exit(FatalError);
                                    sharedEdge=edgeInd1;
                                }
                            }
                        }
                        if(sharedEdge!=-1)
                        {
                            edgeToFaces_[sharedEdge] = edgeFaces[sharedEdge];
                        }
                        else
                        {
                            edge addedEdge(mc33_edg_vert1,mc33_edg_vert2);
                            edgeToFaces_.append(DynamicList<label>(0));                
                            newMeshEdges_.append(addedEdge);
                            edgesToSide_.append(0);
                            edgeToFaces_[edgeToFaces_.size()-1].append(faceInd);
                        }
                    }
                    else
                    {
                        edge addedEdge(mc33_edg_vert1,mc33_edg_vert2);
                        edgeToFaces_.append(DynamicList<label>(0));                
                        newMeshEdges_.append(addedEdge);
                        edgesToSide_.append(0);
                        edgeToFaces_[edgeToFaces_.size()-1].append(faceInd);
                    }
                }
            }
        }
    }

    for(int i=0;i<newMeshEdges_.size();i++)
    {
        if(edgeToFaces_[i].size() != 0)
        {
            for(int k=0;k<edgeToFaces_[i].size();k++)
            {
                faceToEdges_[edgeToFaces_[i][k]].append(i);
            }
        }
    }
    
    /*
    edgeToCells_.setSize(newMeshEdges_.size());
    const labelList& owner = this->faceOwner();
    const labelList& neighbour = this->faceNeighbour();
    const labelListList& edgeCells = this->edgeCells();
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
                    << exit(FatalError);
                }
                label thisFace = edgeToFaces_[i][0];
                edgeToCells_[i].append(owner[thisFace]);
                if(thisFace < neighbour.size())
                    edgeToCells_[i].append(neighbour[thisFace]);
            }
        }
    }
    */
    
    /*
    cellToEdges_.setCapacity(meshCells.size());
    for(int i=0;i<newMeshEdges_.size();i++)
    {
        if(edgeToCells_[i].size() != 0)
        {
            for(int k=0;k<edgeToCells_[i].size();k++)
            {
                cellToEdges_[edgeToCells_[i][k]].append(i);
            }
        }
    }
    */
    
    newMeshEdges_.setCapacity(newMeshEdges_.size());
    edgesToSide_.setCapacity(edgesToSide_.size());
    edgeToFaces_.setCapacity(edgeToFaces_.size());
    faceToEdges_.setCapacity(faceToEdges_.size());
    //edgeToCells_.setCapacity(edgeToCells_.size());
}

/*
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
*/

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
        newMeshPoints_[newMeshEdges_[k].end()]<<"\t";
        Info<<"at faces: ";
        Info<<edgeToFaces_[k]<<endl;
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

/*
void Foam::cutCellFvMesh::groupEdgesToFaces
(
    List<DynamicList<DynamicList<label>>>& cellNonconnectedEdges,
    List<DynamicList<DynamicList<DynamicList<label>>>>& cellNonConnectedMultiPoints,
    List<DynamicList<DynamicList<DynamicList<label>>>>& cellNonConnectedMultiEdges,
    List<DynamicList<DynamicList<face>>>& cellNonConnectedMultiFaces,
    List<DynamicList<DynamicList<std::unordered_set<label>>>>& cellNonConnectedMultiPointMap,
    List<DynamicList<DynamicList<std::unordered_set<label>>>>& cellNonConnectedMultiEdgeMap
)
{
    const cellList& meshCells = this->cells();
    cellNonconnectedEdges           = List<DynamicList<DynamicList<label>>>(meshCells.size());
    cellNonConnectedMultiPoints     = List<DynamicList<DynamicList<DynamicList<label>>>>(meshCells.size());
    cellNonConnectedMultiEdges      = List<DynamicList<DynamicList<DynamicList<label>>>>(meshCells.size());
    cellNonConnectedMultiFaces      = List<DynamicList<DynamicList<face>>>(meshCells.size());
    cellNonConnectedMultiPointMap   = List<DynamicList<DynamicList<std::unordered_set<label>>>>(meshCells.size());
    cellNonConnectedMultiEdgeMap    = List<DynamicList<DynamicList<std::unordered_set<label>>>>(meshCells.size());
    
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
            {
                Info<<"closedCyclePoints:"<<closedCyclePoints<<endl;
                Info<<"closedCycleEdgesList:"<<closedCycleEdgesList<<endl;
                //FatalErrorInFunction<< "Temp stop"<< exit(FatalError);
            }
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
    
    for(int i=0;i<meshCells.size();i++)
    {
        //Info<<"cellNonConnectedMultiPoints[i]:"<<cellNonConnectedMultiPoints[i]<<endl;
        for(int j=0;j<cellNonConnectedMultiPoints[i].size();j++)
        {
            //Info<<"cellNonConnectedMultiPoints[i][j]:"<<cellNonConnectedMultiPoints[i][j]<<endl;
            bool facesWereRemoved = false;
            cellNonConnectedMultiFaces[i].append(DynamicList<face>());
            cellNonConnectedMultiPointMap[i].append(DynamicList<std::unordered_set<label>>());
            cellNonConnectedMultiEdgeMap[i].append(DynamicList<std::unordered_set<label>>());
            for(int k=0;k<cellNonConnectedMultiPoints[i][j].size();k++)
            {
                bool faceInFace = false;
                bool faceContainsFace = false;
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
                for(int l=0;l<thisCell.size();l++)
                {
                    if(problematicFace[thisCell[l]])
                    {
                        face thisFace = basisFaces[thisCell[l]];
                        label nbrOfPntsInFace=0;
                        for(int l=0;l<cellNonConnectedMultiPoints[i][j][k].size();l++)
                        {
                            label cutFacePnt = cellNonConnectedMultiPoints[i][j][k][l];
                            if(cutFacePnt<nbrOfPrevPoints)
                            {
                                if(thisFace.which(cutFacePnt)!=-1)
                                    nbrOfPntsInFace++;
                            }
                            else
                            {
                                label pntEdgeInd = pointToEgde_[cutFacePnt];
                                if(pntEdgeInd==-1)
                                    FatalErrorInFunction<<"Can not happen!"<<exit(FatalError);
                                edge pntEdge = newMeshEdges_[pntEdgeInd];
                                if(thisFace.which(pntEdge.start())!=-1 && thisFace.which(pntEdge.end())!=-1)
                                    nbrOfPntsInFace++;
                            }
                        }
                        if(problematicFacePoints[thisCell[l]]==4)
                        {
                            if(nbrOfPntsInFace==2 || nbrOfPntsInFace==0)
                            {}
                            else if(nbrOfPntsInFace==4)
                            {
                                faceContainsFace=true;
                            }
                            else
                                FatalErrorInFunction<<"Problematic four point face must have 4, 2 or 0 points in Face"<< exit(FatalError);                                
                        }
                        else if(problematicFacePoints[thisCell[l]]==3)
                        {
                            if(nbrOfPntsInFace==2 || nbrOfPntsInFace==0)
                            {}
                            else if(nbrOfPntsInFace==3)
                            {
                                faceContainsFace=true;
                            }
                            else
                                FatalErrorInFunction<<"Problematic four point face must have 3, 2 or 0 points in Face"<< exit(FatalError);                                
                        }
                    }
                }
                if(faceContainsFace && !faceInFace)
                    facesWereRemoved = true;

                if(!faceInFace && !faceContainsFace)
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
                //Info<<"faceInFace:"<<faceInFace<<endl;
                //Info<<"faceContainsFace:"<<faceContainsFace<<endl;
            }

            //Info<<"facesWereRemoved:"<<facesWereRemoved<<endl;
        }
        
        if(i==198218 || i==198219 || i==198538 || i==198539 || i==211018 || i==211019 || i==211338 || i==211339)
        {
            Info<<"------Cell: "<<i<<"--------------------------------------------"<<endl;
            Info<<"cellNonConnectedMultiPoints[i]:"<<cellNonConnectedMultiPoints[i]<<endl;
            Info<<"cellNonConnectedMultiFaces[i]:"<<cellNonConnectedMultiFaces[i]<<endl;
            std::unordered_set<label> newPnts;
            for(int a=0;a<cellNonConnectedMultiPoints[i].size();a++)
            {
                for(int b=0;b<cellNonConnectedMultiPoints[i][a].size();b++)
                {
                    for(int c=0;c<cellNonConnectedMultiPoints[i][a][b].size();c++)
                    {
                        if(cellNonConnectedMultiPoints[i][a][b][c]>=nbrOfPrevPoints)
                            newPnts.insert(cellNonConnectedMultiPoints[i][a][b][c]);
                    }
                }
            }
            for(label pnt:newPnts)
                Info<<"pnt:"<<pnt<<" -> "<<newMeshEdges_[pointToEgde_[pnt]]<<endl; 
        }
    }
}
*/

void Foam::cutCellFvMesh::newMeshFaces_MC33
(
)
{
    const cellList& meshCells = this->cells();
    const faceList& basisFaces = this->faces();
    
    nbrOfPrevFaces = basisFaces.size();
    //nOldFaces = nbrOfPrevFaces = basisFaces.size();
    
    newMeshFaces_.setCapacity(basisFaces.size()*2);
    newMeshFaces_.append(basisFaces);

    facesToSide(newMeshFaces_);
    
    facesToSide_.setCapacity(basisFaces.size()*2);

    faceToCells_.setSize(basisFaces.size());

    cellToFaces_.setSize(meshCells.size());
    for(int cellInd=0;cellInd<meshCells.size();cellInd++)
    {
        MC33::MC33Cube& thisCellCube = mc33CutCellData[cellInd];
        if(thisCellCube.cell==cellInd)
        {
            labelList& verticeIndex = thisCellCube.cutEdgeVerticeIndex;
            for(auto& triangle : thisCellCube.cutTriangles)
            {
                label vertice0 = verticeIndex[std::get<0>(triangle)];
                label vertice1 = verticeIndex[std::get<1>(triangle)];
                label vertice2 = verticeIndex[std::get<2>(triangle)];
                
                if(vertice0==vertice1 && vertice1==vertice2)
                {
                    if(vertice0<nbrOfPrevPoints)
                        continue;
                    else
                        FatalErrorInFunction<<"Multiple old pnts in one face"<<exit(FatalError);
                }
                else if(vertice0==vertice1)
                {
                    if(vertice0<nbrOfPrevPoints)
                        continue;
                    else
                        FatalErrorInFunction<<"Multiple old pnts in one face"<<exit(FatalError);
                }
                else if(vertice1==vertice2)
                {
                    if(vertice1<nbrOfPrevPoints)
                        continue;
                    else
                        FatalErrorInFunction<<"Multiple old pnts in one face"<<exit(FatalError);
                }
                else if(vertice2==vertice0)
                {
                    if(vertice2<nbrOfPrevPoints)
                        continue;
                    else
                        FatalErrorInFunction<<"Multiple old pnts in one face"<<exit(FatalError);
                }
                
                face addedFace(List<label>({vertice0,vertice1,vertice2}));
            
                newMeshFaces_.append(addedFace);
        
                facesToSide_.append(+2);
        
                DynamicList<label> newFaceCell;;
                newFaceCell.append(cellInd);
                faceToCells_.append(newFaceCell);
        
                cellToFaces_[cellInd].append(newMeshFaces_.size()-1);
            }
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

void Foam::cutCellFvMesh::cutOldFaces_MC33
(
)
{
    const faceList& meshFaces = this->faces();
    const labelList& cellOwner = this->faceOwner();
    const labelList& cellNeighbor = this->faceNeighbour();
    const labelListList& faceToEdge = this->faceEdges();
    
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
        
        
        DynamicList<label> addedEdges;
        for(const label cutEdge: faceToEdges_[i])
        {
            if(cutEdge >= nbrOfPrevEdges)
            {
                addedEdges.append(cutEdge);
            }
        }
        
        if(faceToEdges_[i].size()>2)
        {
            FatalErrorInFunction<<"Can not happen"<<exit(FatalError);
        }
        if(faceToEdges_[i].size()==2)
        {
            if(addedEdges.size()!=2)
                FatalErrorInFunction<<"Can not happen"<<exit(FatalError);
            
            edge edg0 = newMeshEdges_[faceToEdges_[i][0]];
            if(edg0.start()<nbrOfPrevPoints || edg0.end()<nbrOfPrevPoints)
                FatalErrorInFunction<<"Can not happen"<<exit(FatalError);
            
            edge edg1 = newMeshEdges_[faceToEdges_[i][1]];
            if(edg1.start()<nbrOfPrevPoints || edg1.end()<nbrOfPrevPoints)
                FatalErrorInFunction<<"Can not happen"<<exit(FatalError);
            
            edge faceEdge1_OfNewEdge0 = newMeshEdges_[pointToEgde_[edg0.start()]];
            edge faceEdge2_OfNewEdge0 = newMeshEdges_[pointToEgde_[edg0.end()]];
            
            edge faceEdge1_OfNewEdge1 = newMeshEdges_[pointToEgde_[edg1.start()]];
            edge faceEdge2_OfNewEdge1 = newMeshEdges_[pointToEgde_[edg1.end()]];
            
            label thirdPointByNewEdge0 = -1;
            thirdPointByNewEdge0 = faceEdge1_OfNewEdge0.commonVertex(faceEdge2_OfNewEdge0);
            if(thirdPointByNewEdge0==-1)
                FatalErrorInFunction<< "thirdPointByNewEdge0 not found! "<< exit(FatalError);
            
            label thirdPointByNewEdge1 = -1;
            thirdPointByNewEdge1 = faceEdge1_OfNewEdge1.commonVertex(faceEdge2_OfNewEdge1);
            if(thirdPointByNewEdge1==-1)
                FatalErrorInFunction<< "thirdPointByNewEdge1 not found! "<< exit(FatalError);                
            
            labelList face1(3);
            face1[0] = edg0.start();
            face1[1] = thirdPointByNewEdge0;
            face1[2] = edg0.end();
            cutFaces_.append(face(face1));
            cutFacesToSide_.append(pointsToSide_[face1[1]]);
            oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
            
            labelList face2(3);
            face2[0] = edg1.start();
            face2[1] = thirdPointByNewEdge1;
            face2[2] = edg1.end();
            cutFaces_.append(face(face2));
            cutFacesToSide_.append(pointsToSide_[face2[1]]);
            oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
            
            labelList face3(6);
            face3[0] = edg0.start();
            face3[1] = edg0.end();
            face3[2] = faceEdge2_OfNewEdge0.otherVertex(thirdPointByNewEdge0);
            face3[3] = -1;
            face3[4] = -1;
            if(faceEdge2_OfNewEdge0.connected(faceEdge1_OfNewEdge1))
            {
                face3[3] = edg1.start();
                face3[4] = edg1.end();
            }
            else if(faceEdge2_OfNewEdge0.connected(faceEdge2_OfNewEdge1))
            {
                face3[3] = edg1.end();
                face3[4] = edg1.start();
            }
            else
                FatalErrorInFunction<< "Fourth and fifth point not found! "<< exit(FatalError);
            face3[5] = faceEdge1_OfNewEdge0.otherVertex(thirdPointByNewEdge0);
            cutFaces_.append(face(face3));
            cutFacesToSide_.append(pointsToSide_[face3[2]]);
            oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
            
        }
        else if(faceToEdges_[i].size()==1)
        {
            if(faceToEdges_[i][0] >= nbrOfPrevEdges)
            {
                edge        addedEdge = newMeshEdges_[faceToEdges_[i][0]];
                std::vector<label> augmentedFace;
                std::vector<label> jumpVertices;
                for(int locEdgeNbr=0; locEdgeNbr<faceToEdge[i].size();locEdgeNbr++)
                {
                    label currEdgeInd = faceToEdge[i][locEdgeNbr];
                    label nextEdgeInd = faceToEdge[i][(locEdgeNbr+1)%faceToEdge[i].size()];
                    edge currEdge = newMeshEdges_[currEdgeInd];
                    edge nextEdge = newMeshEdges_[nextEdgeInd];
                    label interVert = currEdge.commonVertex(nextEdge);
                    if(interVert==-1)
                        FatalErrorInFunction<<"Can not happen. Edges must be in order"                       << exit(FatalError);
                    if(edgeToPoint_[currEdgeInd]!=-1)
                    {
                        augmentedFace.push_back(edgeToPoint_[currEdgeInd]);
                    }
                    augmentedFace.push_back(interVert);
                }
                jumpVertices.resize(augmentedFace.size(),-1);
                edge cutEdgeInds(-1,-1);
                for(unsigned int i=0;i<augmentedFace.size();i++)
                {
                    if(augmentedFace[i]==addedEdge[0])
                        cutEdgeInds[0] = i;
                    if(augmentedFace[i]==addedEdge[1])
                        cutEdgeInds[1] = i;
                }
                if(cutEdgeInds[0]==cutEdgeInds[1])
                    FatalErrorInFunction<<"Can not happen"<<exit(FatalError);
                jumpVertices[cutEdgeInds[0]] = cutEdgeInds[1];
                jumpVertices[cutEdgeInds[1]] = cutEdgeInds[0];
                
                labelList   newFace1(0);
                labelList   newFace2(0);              
                
                int face1Ind = cutEdgeInds[0];
                do
                {
                    newFace1.append(augmentedFace[face1Ind]);
                    face1Ind = (face1Ind+1)%augmentedFace.size();
                }
                while(face1Ind!=cutEdgeInds[1]);
                newFace1.append(augmentedFace[cutEdgeInds[1]]);
                label signFace1 = pointsToSide_[newFace1[1]];
                
                int face2Ind = cutEdgeInds[1];
                do
                {
                    newFace2.append(augmentedFace[face2Ind]);
                    face2Ind = (face2Ind+1)%augmentedFace.size();
                }
                while(face2Ind!=cutEdgeInds[0]);
                newFace2.append(augmentedFace[cutEdgeInds[0]]);
                label signFace2 = pointsToSide_[newFace2[1]];

                oldFacesToCutFaces_[i].setCapacity(2);
                cutFaces_.append(face(newFace1));
                cutFacesToSide_.append(signFace1);
                oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
                cutFaces_.append(face(newFace2));
                cutFacesToSide_.append(signFace2);
                oldFacesToCutFaces_[i].append(cutFaces_.size()-1);
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

bool faceBordersFace(const face& first, const face& second)
{
    label numberOfEqualEdges=0;
    for(edge oneEdgeFirst: first.edges())
    {
        for(edge oneEdgeSecond: second.edges())
        {
            if(oneEdgeFirst==oneEdgeSecond)
            {
                numberOfEqualEdges++;
                continue;
            }
            oneEdgeSecond.flip();
            if(oneEdgeFirst==oneEdgeSecond)
            {
                numberOfEqualEdges++;
                continue;
            }
        }
    }
    if(numberOfEqualEdges==0)
        return false;
    else if(numberOfEqualEdges==1)
        return true;
    else
    {
        Info<<"faceFirst:"<<first<<endl;
        Info<<"faceSecond:"<<second<<endl;
        FatalErrorInFunction<<"Can not happen"<<exit(FatalError);
        return false;
    }
}

void Foam::cutCellFvMesh::createNewMeshData_MC33
(
)
{
    const pointField& meshPoints = this->points();
    const cellList& meshCells = this->cells();
    const faceList& meshFaces = this->faces();
    const labelList owner   = this->faceOwner();
    const labelList neighbour = this->faceNeighbour();    
    const polyBoundaryMesh& boundMesh = this->boundaryMesh();
    const labelListList& pointToEdges = this->pointEdges();
    const labelListList& cellToEdges = this->cellEdges();
    
    nOldFaces = meshFaces.size();
    reverseFaceMap.setSize(nOldFaces,-1);
    DynamicList<label> provisional_faceMap;
    
    nOldCells = meshCells.size();
    reverseCellMap.setSize(nOldCells,-1);
    DynamicList<label> provisional_cellMap;
    oldCellVolumesPtr->setSize(nOldCells,0);
    for(int i=0;i<nOldCells;i++)
    {
        (*oldCellVolumesPtr)[i] = meshCells[i].mag(meshPoints,meshFaces);
    }
    

    // Store old boundary patches
    labelList oldFaceToPatchInd(meshFaces.size(),-1);
    patchStarts = labelList(boundMesh.size());
    patchSizes = labelList(boundMesh.size());
    label cutCellPatchIndex = -1;
    for(int i=0;i<boundMesh.size();i++)
    {
        patchStarts[i] = boundMesh[i].start();
        patchSizes[i] = boundMesh[i].faceCentres().size();
        word namePatch = boundMesh[i].name();
        if(namePatch=="cutCell")
        {
            cutCellPatchIndex = i;
        }
        Info<<namePatch<<"   BoundaryFaceStart:"<<patchStarts[i]<<" FacesSize:"<<patchSizes[i]<<endl;
    }
    for(int i=0;i<patchStarts.size();i++)
    {
        for(int j=0;j<patchSizes[i];j++)
        {
            oldFaceToPatchInd[patchStarts[i]+j] = i;
        }
    }
    oldPatchStarts = patchStarts;
    oldPatchNMeshPoints = patchSizes;
    if(cutCellPatchIndex==-1)
    {
        FatalErrorInFunction<<"Cut Cell patch does not exist"<<exit(FatalError);
    }

    //Prepare data storage for old to new cell index
    oldSplittedCellToNewPlusCell = List<DynamicList<label>>(meshCells.size());
    oldSplittedCellToNewMinusCell = List<DynamicList<label>>(meshCells.size());
    for(int i=0;i<meshCells.size();i++)
    {
        oldSplittedCellToNewPlusCell[i].setSize(0);
        oldSplittedCellToNewMinusCell[i].setSize(0);
    }
        
    // Compute new cellIndexes for added cells
    List<DynamicList<label>> oldCellsToAddedMinusSideCellIndex(meshCells.size());
    deletedCell = List<bool>(meshCells.size(),false);
    label addedCellIndex = 0;
    DynamicList<DynamicList<DynamicList<label>>> cellToNewMinusCellsPointLabels;
    cellToNewMinusCellsPointLabels.setSize(meshCells.size());
    DynamicList<DynamicList<DynamicList<label>>> cellToNewPlusCellsPointLabels;
    cellToNewPlusCellsPointLabels.setSize(meshCells.size());
    
    DynamicList<DynamicList<DynamicList<face>>> cellToNewMinusCellsCutFaces;
    cellToNewMinusCellsCutFaces.setSize(meshCells.size());
    DynamicList<DynamicList<DynamicList<face>>> cellToNewPlusCellsCutFaces;
    cellToNewPlusCellsCutFaces.setSize(meshCells.size());

    DynamicList<DynamicList<label>> cellToNewMinusCellsIndexes;
    cellToNewMinusCellsIndexes.setSize(meshCells.size());
    DynamicList<DynamicList<label>> cellToNewPlusCellsIndexes;
    cellToNewPlusCellsIndexes.setSize(meshCells.size());
    
    for(int i=0;i<meshCells.size();i++)
    {
        if(cellToFaces_[i].size() > 0)
        {
            if(cellsToSide_[i] != 0)
            {                
                FatalErrorInFunction<<"Non Zero cell has face inside!"<<exit(FatalError);
            }
            
            std::unordered_set<label> edgesTreated;
            DynamicList<DynamicList<label>> minusCells;
            DynamicList<DynamicList<label>> plusCells;
            DynamicList<DynamicList<face>> minusCellsCutFaces;
            DynamicList<DynamicList<face>> plusCellsCutFaces;
            DynamicList<std::unordered_set<label>> minusCellsCutFacesMap;
            DynamicList<std::unordered_set<label>> plusCellsCutFacesMap;
            const labelList& thisCellEdges = cellToEdges[i];
            oldSplittedCellToNewPlusCell[i] = i;
            
            std::vector<std::vector<label>> adjacentFaces(cellToFaces_[i].size());
            std::vector<bool> treatedFaces(cellToFaces_[i].size(),false);
            for(int j=0;j<cellToFaces_[i].size();j++)
            {
                face facej = newMeshFaces_[cellToFaces_[i][j]];
                for(int k=0;k<cellToFaces_[i].size();k++)
                {
                    face facek = newMeshFaces_[cellToFaces_[i][k]];
                    if(j!=k && faceBordersFace(facej,facek))
                    {
                        adjacentFaces[j].push_back(k);
                    }
                }
            }
            bool allFacesTreated=false;
            DynamicList<DynamicList<label>> connectedFacesLocInd;
            std::vector<label> front;
            while(!allFacesTreated)
            {
                if(front.size()==0)
                {
                    label nextFrontFace = -1;
                    for(int j=0;j<cellToFaces_[i].size();j++)
                    {
                        if(!treatedFaces[j])
                        {
                            nextFrontFace = j;
                            break;
                        }
                    }
                    if(nextFrontFace==-1)
                        allFacesTreated=true;
                    else
                    {
                        front.push_back(nextFrontFace);
                        connectedFacesLocInd.append(DynamicList<label>());
                    }
                }
                else
                {
                    std::vector<label> nextFront;
                    for(label faceInd: front)
                    {
                        treatedFaces[faceInd] = true;
                        for(label conFace : adjacentFaces[faceInd])
                        {
                            if(!treatedFaces[conFace])
                            {
                                nextFront.push_back(conFace);
                            }
                        }
                        connectedFacesLocInd.last().append(faceInd);
                    }
                    front = nextFront;
                }
            }
            
            if(connectedFacesLocInd.size()>1)
                FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
            
            DynamicList<DynamicList<label>> connectedFaces;
            connectedFaces.setSize(connectedFacesLocInd.size());
            for(int j=0;j<connectedFacesLocInd.size();j++)
            {
                for(label oneFaceInd: connectedFacesLocInd[j])
                {
                    connectedFaces[j].append(cellToFaces_[i][oneFaceInd]);
                }
            }            
            
            for(int j=0;j<connectedFaces.size();j++)
            {
                for(label oneFace: connectedFaces[j])
                {
                    if(oneFace < nbrOfPrevFaces)
                        FatalErrorInFunction<<"Can not happen"<<exit(FatalError);
                }
                
                std::unordered_set<label> uniConFacePnts;
                for(label oneFaceInd: connectedFaces[j])
                {
                    for(label onePnt: newMeshFaces_[oneFaceInd])
                    {
                        uniConFacePnts.insert(onePnt);
                    }
                }
                std::unordered_set<label> pointsTreated;                    
                DynamicList<label> minusCell;
                DynamicList<label> plusCell;
                
                DynamicList<face> minusCellCutFaces;
                DynamicList<face> plusCellCutFaces;
                for(label oneFaceInd: connectedFaces[j])
                {
                    face cutFace = newMeshFaces_[oneFaceInd];
                    minusCellCutFaces.append(cutFace);
                    plusCellCutFaces.append(cutFace);
                }
                
                DynamicList<DynamicList<edge>> facePointEdges;
                facePointEdges.setSize(uniConFacePnts.size());
                //store the points around each cutFace cutted edge via facePointEdges
                label k=0;
                for(auto iter=uniConFacePnts.begin();iter!=uniConFacePnts.end();iter++,k++)
                {
                    label pointInd = *iter;
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
                    for(int l=0;l<edgeInds.size();l++)
                    {
                        if(edgesTreated.count(edgeInds[l]) == 0)
                        {
                            facePointEdges[k].append(newMeshEdges_[edgeInds[l]]);
                            edgesTreated.insert(edgeInds[l]);
                        }
                    }
                    pointsTreated.insert(pointInd);
                }

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
                            if(pointsTreated.count(edgePoints[n])==0)
                            {
                                if(pointsToSide_[edgePoints[n]] == 0)
                                {
                                    FatalErrorInFunction<<"Can not happen."<<exit(FatalError);                                                                        
                                }
                                else if(pointsToSide_[edgePoints[n]] == 1)
                                {
                                    plusCell.append(edgePoints[n]);
                                    plusCellFrontPoints.append(edgePoints[n]);
                                    pointsTreated.insert(edgePoints[n]);
                                }
                                else if(pointsToSide_[edgePoints[n]] == -1)
                                {
                                    minusCell.append(edgePoints[n]);
                                    minusCellFrontPoints.append(edgePoints[n]);
                                    pointsTreated.insert(edgePoints[n]);
                                }
                                else
                                    FatalErrorInFunction<<"Point side must bei -1,0,1."<<exit(FatalError);
                            }
                        }
                    }
                }

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

                                if(otherPoint!=-1 && pointsTreated.count(otherPoint)==0)
                                {
                                    //Info<<" added";
                                    
                                    if(pointsToSide_[otherPoint] == -1)
                                    {
                                        Info<<"thisCellEdges:"<<thisCellEdges<<endl;
                                        Info<<"nbrOfPrevEdges:"<<nbrOfPrevEdges<<endl;
                                        for(int n=0;n<thisCellEdges.size();n++)
                                        {
                                            Info<<thisCellEdges[n]<<":"<<newMeshEdges_[thisCellEdges[n]]<<" -- side:"<<edgesToSide_[thisCellEdges[n]]<<" firstside:"<<pointsToSide_[newMeshEdges_[thisCellEdges[n]].start()]<<" secondside:"<<pointsToSide_[newMeshEdges_[thisCellEdges[n]].end()]<<endl;
                                        }
                                        Info<<"----------m:"<<m<<endl;
                                        Info<<thisCellEdges[m]<<":"<<newMeshEdges_[thisCellEdges[m]]<<" -- side:"<<edgesToSide_[thisCellEdges[m]]<<" firstside:"<<pointsToSide_[newMeshEdges_[thisCellEdges[m]].start()]<<" secondside:"<<pointsToSide_[newMeshEdges_[thisCellEdges[m]].end()]<<endl;                                            
                                        FatalErrorInFunction<<"Edge with 1 side has non 1 point"<<exit(FatalError);
                                    }
                                    
                                    plusCell.append(otherPoint);
                                    plusCellFrontPoints.append(otherPoint);
                                    edgesTreated.insert(thisCellEdges[m]);
                                    pointsTreated.insert(otherPoint);
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
                                if(otherPoint!=-1 && pointsTreated.count(otherPoint)==0)
                                {
                                    //Info<<" added";
                                    if(pointsToSide_[otherPoint] == 1)
                                    {
                                        Info<<"thisCellEdges:"<<thisCellEdges<<endl;
                                        Info<<"nbrOfPrevEdges:"<<nbrOfPrevEdges<<endl;
                                        for(int n=0;n<thisCellEdges.size();n++)
                                        {
                                            Info<<thisCellEdges[n]<<":"<<newMeshEdges_[thisCellEdges[n]]<<" -- side:"<<edgesToSide_[thisCellEdges[n]]<<" firstside:"<<pointsToSide_[newMeshEdges_[thisCellEdges[n]].start()]<<" secondside:"<<pointsToSide_[newMeshEdges_[thisCellEdges[n]].end()]<<endl;
                                        }
                                        Info<<"----------m:"<<m<<endl;
                                        Info<<thisCellEdges[m]<<":"<<newMeshEdges_[thisCellEdges[m]]<<" -- side:"<<edgesToSide_[thisCellEdges[m]]<<" firstside:"<<pointsToSide_[newMeshEdges_[thisCellEdges[m]].start()]<<" secondside:"<<pointsToSide_[newMeshEdges_[thisCellEdges[m]].end()]<<endl;                                            
                                        FatalErrorInFunction<<"Edge with 1 side has non 1 point"<<exit(FatalError);
                                    }                                    
                                    minusCell.append(otherPoint);
                                    minusCellFrontPoints.append(otherPoint);
                                    edgesTreated.insert(thisCellEdges[m]);
                                    pointsTreated.insert(otherPoint);
                                }
                                //Info<<endl;
                            }                                
                        }
                    }
                }
                if(plusCell.size()==0 && minusCell.size()==0)
                {
                    Info<<endl<<endl;
                    Info<<"cell:"<<i<<endl;
                    Info<<"j:"<<j<<endl;
                    Info<<"facePointEdges:"<<facePointEdges<<endl;
                    for(label oneFaceInd: connectedFaces[j])
                    {
                        face cutFace = newMeshFaces_[oneFaceInd];
                        Info<<"cutFace:"<<cutFace<<endl;
                    }
                    Info<<"cell faces:"<<meshCells[i]<<endl;
                    Info<<"cell Points: "<<meshCells[i].labels(meshFaces)<<endl;
                    Info<<"cellToFaces_:"<<cellToFaces_[i]<<endl;
                    for(int k=0;k<cellToFaces_[i].size();k++)
                    {
                        Info<<"   cellToFaces:"<<cellToFaces_[i][k]<<endl;
                        Info<<"   newMeshFaces_["<<cellToFaces_[i][k]<<"]:"<<newMeshFaces_[cellToFaces_[i][k]]<<endl;                            
                    }
                    Info<<"nbrNewFaces:"<<newMeshFaces_.size()<<endl;
                    Info<<"newMeshFace:"<<newMeshFaces_[cellToFaces_[i][j]]<<endl;
                    Info<<"plusCellFrontPoints:"<<plusCellFrontPoints<<endl;
                    Info<<"minusCellFrontPoints:"<<minusCellFrontPoints<<endl;
                    Info<<"plusCell:"<<plusCell<<endl;
                    Info<<"minusCell:"<<minusCell<<endl;
                    Info<<endl;
                    Info<<"plusCells:"<<plusCells<<endl;
                    Info<<"minusCells:"<<minusCells<<endl;
                    FatalErrorInFunction<<"Face does not create at least one cell."<<exit(FatalError);
                }
                
                label minusCellAlreadyDoneInd = -1;
                label plusCellAlreadyDoneInd = -1;
                for(int k=0;k<minusCellsCutFacesMap.size();k++)
                {
                    std::unordered_set<label>& oneCellMap = minusCellsCutFacesMap[k];
                    bool cellDone = true;
                    for(const label& vertice: minusCell)
                        if(oneCellMap.count(vertice)==0)
                            cellDone = false;
                    if(cellDone)
                    {
                        if(minusCellAlreadyDoneInd!=-1)
                            FatalErrorInFunction<<"Cell part already done multiple times. Can not happen"<<exit(FatalError);
                        minusCellAlreadyDoneInd = k;
                    }
                }
                for(int k=0;k<plusCellsCutFacesMap.size();k++)
                {
                    std::unordered_set<label>& oneCellMap = plusCellsCutFacesMap[k];
                    bool cellDone = true;
                    for(const label& vertice: plusCell)
                        if(oneCellMap.count(vertice)==0)
                            cellDone = false;
                    if(cellDone)
                    {
                        if(plusCellAlreadyDoneInd!=-1)
                            FatalErrorInFunction<<"Cell part already done multiple times. Can not happen"<<exit(FatalError);
                        plusCellAlreadyDoneInd = k;
                    }
                }
                
                if(minusCellAlreadyDoneInd!=-1 && plusCellAlreadyDoneInd!=-1)
                {
                    Info<<Foam::endl;
                    Info<<"j:"<<j<<Foam::endl;
                    Info<<"meshCell:"<<meshCells[i].labels(meshFaces)<<endl;
                    Info<<"pointToSide [";
                    for(label vertice: meshCells[i].labels(meshFaces))
                    {
                        Info<<vertice<<"  "<<pointsToSide_[vertice]<<Foam::endl;
                    }
                    Info<<"]"<<Foam::endl;                        Info<<"cellToFaces_["<<i<<"]:"<<cellToFaces_[i]<<Foam::endl;
                    Info<<"cut faces[";
                    for(label face: cellToFaces_[i])
                    {
                        Info<<newMeshFaces_[face]<<" ";
                    }
                    Info<<"]"<<Foam::endl;
                    Info<<"cellfaces [";
                    for(label face: meshCells[i])
                    {
                        Info<<newMeshFaces_[face]<<" ";
                    }
                    Info<<"]"<<Foam::endl;
                    Info<<"plusCells:"<<plusCells<<endl;
                    Info<<"minusCells:"<<minusCells<<endl;
                    FatalErrorInFunction<<"Cut Face neigbors two already done cell parts. Can not happen"<<exit(FatalError);
                }
                else if(minusCellAlreadyDoneInd!=-1)
                {
                    plusCells.append(plusCell);
                    
                    plusCellsCutFacesMap.append(std::unordered_set<label>());
                    for(const label& vertice : plusCell)
                        plusCellsCutFacesMap.last().insert(vertice);
                    
                    plusCellsCutFaces.append(plusCellCutFaces);
                    
                    for(label oneFaceInd: connectedFaces[j])
                    {
                        face cutFace = newMeshFaces_[oneFaceInd];
                        minusCellsCutFaces[minusCellAlreadyDoneInd].append(cutFace);
                    }
                }
                else if(plusCellAlreadyDoneInd!=-1)
                {
                    minusCells.append(minusCell);
                    
                    minusCellsCutFacesMap.append(std::unordered_set<label>());
                    for(const label& vertice : minusCell)
                        minusCellsCutFacesMap.last().insert(vertice);
                    
                    minusCellsCutFaces.append(minusCellCutFaces);
                    for(label oneFaceInd: connectedFaces[j])
                    {
                        face cutFace = newMeshFaces_[oneFaceInd];
                        plusCellsCutFaces[plusCellAlreadyDoneInd].append(cutFace);
                    }
                }
                else
                {
                    minusCells.append(minusCell);
                    minusCellsCutFacesMap.append(std::unordered_set<label>());
                    for(const label& vertice : minusCell)
                        minusCellsCutFacesMap.last().insert(vertice);
                    minusCellsCutFaces.append(minusCellCutFaces);
                    
                    plusCells.append(plusCell);
                    plusCellsCutFacesMap.append(std::unordered_set<label>());
                    for(const label& vertice : plusCell)
                        plusCellsCutFacesMap.last().insert(vertice);
                    plusCellsCutFaces.append(plusCellCutFaces);
                }
            }

            label countOldPoints = 0;
            for(int j=0;j<plusCells.size();j++)
                countOldPoints+=plusCells[j].size();
            for(int j=0;j<minusCells.size();j++)
                countOldPoints+=minusCells[j].size();
            DynamicList<label> cutFaces = cellToFaces_[i];
            std::unordered_set<label> cutPoints;
            for(int j=0;j<cutFaces.size();j++)
            {
                face cutFace = newMeshFaces_[cutFaces[j]];
                for(int k=0;k<cutFace.size();k++)
                {
                    if(cutFace[k]<nbrOfPrevPoints &&
                       pointsToSide_[cutFace[k]]==1 && 
                       cutPoints.count(cutFace[k])==0
                    )
                    {
                        countOldPoints++;
                        cutPoints.insert(cutFace[k]);
                    }
                }
            }

            if(countOldPoints!=8)
            {
                Info<<endl;
                Info<<"cellInd:"<<i<<endl;
                for(label vertice: meshCells[i].labels(meshFaces))
                {
                    Info<<"vertice:"<<vertice<<"  "<<meshPoints[vertice]<<" "<<pointDist[vertice]<<endl;
                }
                Info<<"countOldPoints:"<<countOldPoints<<endl;
                for(label point:cutPoints)
                    Info<<point<<" ";
                Info<<endl;
                for(int j=0;j<cellToFaces_[i].size();j++)
                    Info<<"cutface "<<j<<" "<<cellToFaces_[i][j]<<": "<<newMeshFaces_[cellToFaces_[i][j]]<<endl;
                Info<<"connectedFacesLocInd:"<<connectedFacesLocInd<<Foam::endl;
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

            if((connectedFacesLocInd.size()+1)!=(minusCells.size()+plusCells.size()))
                FatalErrorInFunction<<"Number of plus and minus cells does not match the number of cut faces"<<exit(FatalError);
            
            cellToNewMinusCellsPointLabels[i] = minusCells;
            cellToNewPlusCellsPointLabels[i] = plusCells;
            
            cellToNewMinusCellsCutFaces[i] = minusCellsCutFaces;
            cellToNewPlusCellsCutFaces[i] = plusCellsCutFaces;
            
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
                deletedCell[i] = true;
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
        else if(cellsToSide_[i] == -1)
        {
            deletedCell[i] = true;
        }
    }
    
    label maxCellInd = deletedCell.size()-1;
    for(int i=0;i<cellToNewMinusCellsIndexes.size();i++)
        for(int j=0;j<cellToNewMinusCellsIndexes[i].size();j++)
            maxCellInd = std::max(cellToNewMinusCellsIndexes[i][j],maxCellInd);
    for(int i=0;i<cellToNewPlusCellsIndexes.size();i++)
        for(int j=0;j<cellToNewPlusCellsIndexes[i].size();j++)
            maxCellInd = std::max(cellToNewPlusCellsIndexes[i][j],maxCellInd);
    maxCellInd++;
            
    List<bool> correctdeletedCell(maxCellInd,false);        
    for(int i=0;i<deletedCell.size();i++)
        correctdeletedCell[i] = deletedCell[i];
    
    /*
    labelList deletedCellSIZE(Pstream::nProcs(),0);
    deletedCellSIZE[Pstream::myProcNo()] = deletedCell.size();
    Pstream::gatherList(deletedCellSIZE);
    if(Pstream::master())
        std::cout<<"--------deletedCellSIZE:"<<Pstream::myProcNo()<<"----"<<deletedCellSIZE[0]<<","<<deletedCellSIZE[1]<<","<<deletedCellSIZE[2]<<","<<deletedCellSIZE[3]<<"----------"<<std::endl;
    
    labelList correctdeletedCellSIZE(Pstream::nProcs(),0);
    correctdeletedCellSIZE[Pstream::myProcNo()] = correctdeletedCell.size();
    Pstream::gatherList(correctdeletedCellSIZE);
    if(Pstream::master())
        std::cout<<"--------correctdeletedCellSIZE:"<<Pstream::myProcNo()<<"----"<<correctdeletedCellSIZE[0]<<","<<correctdeletedCellSIZE[1]<<","<<correctdeletedCellSIZE[2]<<","<<correctdeletedCellSIZE[3]<<"----------"<<std::endl;
    */

    for(int i=0;i<cellToNewMinusCellsIndexes.size();i++)
    {
        for(int j=0;j<cellToNewMinusCellsIndexes[i].size();j++)
        {
            label minusCellInd = cellToNewMinusCellsIndexes[i][j];
            if(minusCellInd<deletedCell.size())
            {
                if(!correctdeletedCell[minusCellInd])
                    FatalErrorInFunction<<"Minus Cell is not marked as deleted! Can not happen!"<<exit(FatalError);
            }
            else
                correctdeletedCell[minusCellInd] = true;
        }
    }    
    deletedCell = correctdeletedCell;
    
    /*
    deletedCellSIZE = labelList(Pstream::nProcs(),0);
    deletedCellSIZE[Pstream::myProcNo()] = deletedCell.size();
    Pstream::gatherList(deletedCellSIZE);
    if(Pstream::master())
        std::cout<<"--------deletedCellSIZE:"<<Pstream::myProcNo()<<"----"<<deletedCellSIZE[0]<<","<<deletedCellSIZE[1]<<","<<deletedCellSIZE[2]<<","<<deletedCellSIZE[3]<<"----------"<<std::endl;
    */
    
//Barrier(true);
        
    List<label> newCellToSide(maxCellInd+1,0);
    for(int i=0;i<cellToNewMinusCellsIndexes.size();i++)
        for(int j=0;j<cellToNewMinusCellsIndexes[i].size();j++)
            if(newCellToSide[cellToNewMinusCellsIndexes[i][j]]==0)
                newCellToSide[cellToNewMinusCellsIndexes[i][j]] = -1;
            else
                FatalErrorInFunction<<"Double assigned cell "<<i<<" "<<cellToNewMinusCellsIndexes[i][j]<<" from minus as "<<newCellToSide[cellToNewMinusCellsIndexes[i][j]]<<exit(FatalError);
    for(int i=0;i<cellToNewPlusCellsIndexes.size();i++)
        for(int j=0;j<cellToNewPlusCellsIndexes[i].size();j++)
            if(newCellToSide[cellToNewPlusCellsIndexes[i][j]]==0)
                newCellToSide[cellToNewPlusCellsIndexes[i][j]] = +1;
            else
                FatalErrorInFunction<<"Double assigned cell "<<i<<" "<<cellToNewPlusCellsIndexes[i][j]<<" from plus as "<<newCellToSide[cellToNewPlusCellsIndexes[i][j]]<<exit(FatalError);

    bool inZeroCell = false;
    for(int i=0;i<newCellToSide.size();i++)
    {
        if(newCellToSide[i]==0 && cellToFaces_[i].size()>0 && !inZeroCell)
        {
            inZeroCell = true;
        }
        else if(newCellToSide[i]!=0 && cellToFaces_[i].size()>0 && inZeroCell)
        {
            inZeroCell = false;
        }
    }
        
//Barrier(true);

    // Compute List of new faces splitting old cells
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
                    addedCutFacesNeighbor.append(-1);
                    addedCutFacesOwner.append(i);
                    addedCutFacesPatchInd.append(cutCellPatchIndex);
                }
                else if(oldSplittedCellToNewPlusCell[i].size()==1 && oldSplittedCellToNewMinusCell[i].size()>1)
                {
                    addedCutFaces.append(addedFace);
                    addedCutFacesNeighbor.append(-1);
                    addedCutFacesOwner.append(i);
                    addedCutFacesPatchInd.append(cutCellPatchIndex);
                }
                else if(oldSplittedCellToNewPlusCell[i].size()>1 && oldSplittedCellToNewMinusCell[i].size()==1)
                {
                    addedCutFaces.append(addedFace);
                    addedCutFacesNeighbor.append(-1);
                    addedCutFacesOwner.append(oldSplittedCellToNewPlusCell[i][j]);
                    addedCutFacesPatchInd.append(cutCellPatchIndex);                    
                }
                else
                    FatalErrorInFunction<<"This combination is not possible"<<exit(FatalError);
            }
        }
    }

//Barrier(true);
    
    // Compute the List of new faces resulting from the splitting of old faces
    for(int i=0;i<neighbour.size();i++)
    {
        if(oldFacesToCutFaces_[i].size()==1 || oldFacesToCutFaces_[i].size()>5)
        {
            Info<<"i:"<<i<<endl;
            Info<<"meshPoints.size():"<<meshPoints.size()<<endl;
            Info<<"newMeshPoints_.size():"<<newMeshPoints_.size()<<endl;
            Info<<"oldFacesToCutFaces_["<<i<<"]:"<<oldFacesToCutFaces_[i]<<endl;
            Info<<"meshFaces:"<<meshFaces[i]<<endl;
            Info<<"meshFaces:"<<meshFaces[i].points(meshPoints)<<endl;
            for(int j=0;j<oldFacesToCutFaces_[i].size();j++)
            {
                Info<<oldFacesToCutFaces_[i][j]<<" -- "<<cutFaces_[oldFacesToCutFaces_[i][j]]<<" -- ";
                for(int k=0;k<cutFaces_[oldFacesToCutFaces_[i][j]].size();k++)
                {
                    Info<<newMeshPoints_[cutFaces_[oldFacesToCutFaces_[i][j]][k]]<<" ";
                }
                Info<<endl;
            }
            FatalErrorInFunction<< "Face cut in one or more than five cut faces."<< exit(FatalError);
        }        
        for(int j=0;j<oldFacesToCutFaces_[i].size();j++)
        {
            face face1      = cutFaces_[oldFacesToCutFaces_[i][j]];
            label signFace1 = cutFacesToSide_[oldFacesToCutFaces_[i][j]];
                        
            if(signFace1>=0)
            {
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

                DynamicList<DynamicList<face>> ownerCellToNewMinusCellsCutFaces = cellToNewMinusCellsCutFaces[oldOwnerCell];
                DynamicList<DynamicList<face>> ownerCellToNewPlusCellsCutFaces = cellToNewPlusCellsCutFaces[oldOwnerCell];
                DynamicList<DynamicList<face>> neighbourCellToNewMinusCellsCutFaces = cellToNewMinusCellsCutFaces[oldNeighbourCell];
                DynamicList<DynamicList<face>> neighbourCellToNewPlusCellsCutFaces = cellToNewPlusCellsCutFaces[oldNeighbourCell];
                
                DynamicList<std::unordered_set<label>> ownerNewMinusCellsPointLabelsMap;
                ownerNewMinusCellsPointLabelsMap.setSize(ownerNewMinusCellsPointLabels.size());
                for(int k=0;k<ownerNewMinusCellsPointLabels.size();k++)
                {
                    for(int l=0;l<ownerNewMinusCellsPointLabels[k].size();l++)
                        ownerNewMinusCellsPointLabelsMap[k].insert(ownerNewMinusCellsPointLabels[k][l]);
                    for(int l=0;l<ownerCellToNewMinusCellsCutFaces[k].size();l++)
                        for(int m=0;m<ownerCellToNewMinusCellsCutFaces[k][l].size();m++)
                            ownerNewMinusCellsPointLabelsMap[k].insert(ownerCellToNewMinusCellsCutFaces[k][l][m]);
                }
                    
                DynamicList<std::unordered_set<label>> ownerNewPlusCellsPointLabelsMap;
                ownerNewPlusCellsPointLabelsMap.setSize(ownerNewPlusCellsPointLabels.size());
                for(int k=0;k<ownerNewPlusCellsPointLabels.size();k++)
                {
                    for(int l=0;l<ownerNewPlusCellsPointLabels[k].size();l++)
                        ownerNewPlusCellsPointLabelsMap[k].insert(ownerNewPlusCellsPointLabels[k][l]);
                    for(int l=0;l<ownerCellToNewPlusCellsCutFaces[k].size();l++)
                        for(int m=0;m<ownerCellToNewPlusCellsCutFaces[k][l].size();m++)
                            ownerNewPlusCellsPointLabelsMap[k].insert(ownerCellToNewPlusCellsCutFaces[k][l][m]);
                }

                DynamicList<std::unordered_set<label>> neighbourNewMinusCellsPointLabelsMap;
                neighbourNewMinusCellsPointLabelsMap.setSize(neighbourNewMinusCellsPointLabels.size());
                for(int k=0;k<neighbourNewMinusCellsPointLabels.size();k++)
                {
                    for(int l=0;l<neighbourNewMinusCellsPointLabels[k].size();l++)
                        neighbourNewMinusCellsPointLabelsMap[k].insert(neighbourNewMinusCellsPointLabels[k][l]);
                    for(int l=0;l<neighbourCellToNewMinusCellsCutFaces[k].size();l++)
                        for(int m=0;m<neighbourCellToNewMinusCellsCutFaces[k][l].size();m++)
                            neighbourNewMinusCellsPointLabelsMap[k].insert(neighbourCellToNewMinusCellsCutFaces[k][l][m]);
                }
                    
                DynamicList<std::unordered_set<label>> neighbourNewPlusCellsPointLabelsMap;
                neighbourNewPlusCellsPointLabelsMap.setSize(neighbourNewPlusCellsPointLabels.size());
                for(int k=0;k<neighbourNewPlusCellsPointLabels.size();k++)
                {
                    for(int l=0;l<neighbourNewPlusCellsPointLabels[k].size();l++)
                        neighbourNewPlusCellsPointLabelsMap[k].insert(neighbourNewPlusCellsPointLabels[k][l]);
                    for(int l=0;l<neighbourCellToNewPlusCellsCutFaces[k].size();l++)
                        for(int m=0;m<neighbourCellToNewPlusCellsCutFaces[k][l].size();m++)
                            neighbourNewPlusCellsPointLabelsMap[k].insert(neighbourCellToNewPlusCellsCutFaces[k][l][m]);
                }
                    
                label newOwnerCell = -1;
                label newOwnerCellPlusOrMinus = 0;
                for(int k=0;k<ownerNewMinusCellsPointLabelsMap.size();k++)
                {
                    bool allFacePntsInside = true;
                    for(int l=0;l<face1.size();l++)
                    {
                        if(ownerNewMinusCellsPointLabelsMap[k].count(face1[l])==0)
                        {
                            allFacePntsInside = false;
                        }
                    }
                    if(allFacePntsInside)
                    {
                        if(newOwnerCell==-1)
                        {
                            newOwnerCell = ownerNewMinusCellsIndex[k];
                            newOwnerCellPlusOrMinus = -1;
                        }
                        else
                            FatalErrorInFunction<< "Old splitted face neigbours two newSplit Cells." << exit(FatalError);
                    }
                }
                for(int k=0;k<ownerNewPlusCellsPointLabelsMap.size();k++)
                {
                    bool allFacePntsInside = true;
                    for(int l=0;l<face1.size();l++)
                    {
                        if(ownerNewPlusCellsPointLabelsMap[k].count(face1[l])==0)
                        {
                            allFacePntsInside = false;
                        }
                    }
                    if(allFacePntsInside)
                    {
                        if(newOwnerCell==-1)
                        {
                            newOwnerCell = ownerNewPlusCellsIndex[k];
                            newOwnerCellPlusOrMinus = +1;
                        }
                        else
                            FatalErrorInFunction<< "Old splitted face neigbours two newSplit Cells." << exit(FatalError);
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
                label newNeighborCellPlusOrMinus = 0;
                for(int k=0;k<neighbourNewMinusCellsPointLabelsMap.size();k++)
                {
                    bool allFacePntsInside = true;
                    for(int l=0;l<face1.size();l++)
                    {
                        if(neighbourNewMinusCellsPointLabelsMap[k].count(face1[l])==0)
                        {
                            allFacePntsInside = false;
                        }
                    }
                    if(allFacePntsInside)
                    {
                        if(newNeighborCell==-1)
                        {
                            newNeighborCell = neighbourNewMinusCellsIndex[k];
                            newNeighborCellPlusOrMinus = -1;
                        }
                        else
                            FatalErrorInFunction<< "Old splitted face neigbours two newSplit Cells." << exit(FatalError);
                    }
                }
                for(int k=0;k<neighbourNewPlusCellsPointLabelsMap.size();k++)
                {
                    bool allFacePntsInside = true;
                    for(int l=0;l<face1.size();l++)
                    {
                        if(neighbourNewPlusCellsPointLabelsMap[k].count(face1[l])==0)
                        {
                            allFacePntsInside = false;
                        }
                    }
                    if(allFacePntsInside)
                    {
                        if(newNeighborCell==-1)
                        {
                            newNeighborCell = neighbourNewPlusCellsIndex[k];
                            newNeighborCellPlusOrMinus = +1;
                        }
                        else
                            FatalErrorInFunction<< "Old splitted face neigbours two newSplit Cells." << exit(FatalError);
                    }
                }
                if(newNeighborCell==-1 && signFace1 != 0)
                {
                    /*
                    Info<<"-------------------------------------------------------"<<endl;
                    Info<<"oldOwnerCell:"<<oldOwnerCell<<endl;
                    Info<<"oldNeighbourCell:"<<oldNeighbourCell<<endl;
                    Info<<"labels meshCells["<<oldOwnerCell<<"]:"<<meshCells[oldOwnerCell].labels(meshFaces)<<endl;
                    Info<<"points meshCells["<<oldOwnerCell<<"]:"<<meshCells[oldOwnerCell].points(meshFaces,meshPoints)<<endl;
                    Info<<"labels meshCells["<<oldNeighbourCell<<"]:"<<meshCells[oldNeighbourCell].labels(meshFaces)<<endl;
                    Info<<"points meshCells["<<oldNeighbourCell<<"]:"<<meshCells[oldNeighbourCell].points(meshFaces,meshPoints)<<endl;
                    Info<<"signFace1:"<<signFace1<<endl;
                    Info<<"face1:"<<face1<<endl;
                    for(int x=0;x<face1.size();x++)
                    {
                        Info<<newMeshPoints_[face1[x]]<<"  ";
                    }
                    Info<<endl<<"pointsToSide_:  ";
                    for(int x=0;x<face1.size();x++)
                    {
                        Info<<pointsToSide_[face1[x]]<<"  ";
                    }
                    Info<<endl<<endl;
                    Info<<"oldFace:"<<meshFaces[i]<<endl;
                    Info<<"newMeshPoints_:  ";
                    for(int x=0;x<meshFaces[i].size();x++)
                    {
                        Info<<newMeshPoints_[meshFaces[i][x]]<<"  ";
                    }
                    Info<<endl<<"pointsToSide_:  ";
                    for(int x=0;x<meshFaces[i].size();x++)
                    {
                        Info<<pointsToSide_[meshFaces[i][x]]<<"  ";
                    }
                    Info<<endl<<endl;
                    Info<<"neighbourNewMinusCellsPointLabels:"<<neighbourNewMinusCellsPointLabels<<endl;
                    Info<<"neighbourNewPlusCellsPointLabels:"<<neighbourNewPlusCellsPointLabels<<endl;
                    Info<<endl<<endl;
                    */
                    
                    Info<<"---"<<endl;
                    Info<<"ownerNewMinusCellsIndex:"<<ownerNewMinusCellsIndex<<endl;
                    Info<<"ownerNewMinusCellsPointLabels:"<<ownerNewMinusCellsPointLabels<<endl;
                    Info<<"ownerCellToNewMinusCellsCutFaces:"<<ownerCellToNewMinusCellsCutFaces<<endl;
                    Info<<"---"<<endl;
                    Info<<"ownerNewPlusCellsIndex:"<<ownerNewPlusCellsIndex<<endl;
                    Info<<"ownerNewPlusCellsPointLabels:"<<ownerNewPlusCellsPointLabels<<endl;
                    Info<<"ownerCellToNewPlusCellsCutFaces:"<<ownerCellToNewPlusCellsCutFaces<<endl;
                    Info<<"---"<<endl;
                    Info<<"neighbourNewMinusCellsIndex:"<<neighbourNewMinusCellsIndex<<endl;
                    Info<<"neighbourNewMinusCellsPointLabels:"<<neighbourNewMinusCellsPointLabels<<endl;
                    Info<<"neighbourCellToNewMinusCellsCutFaces:"<<neighbourCellToNewMinusCellsCutFaces<<endl;
                    Info<<"---"<<endl;
                    Info<<"neighbourNewPlusCellsIndex:"<<neighbourNewPlusCellsIndex<<endl;
                    Info<<"neighbourNewPlusCellsPointLabels:"<<neighbourNewPlusCellsPointLabels<<endl;
                    Info<<"neighbourCellToNewPlusCellsCutFaces:"<<neighbourCellToNewPlusCellsCutFaces<<endl;
                    
                    Info<<"---"<<endl;
                    Info<<"oldOwnerCell:"<<oldOwnerCell<<endl;
                    Info<<"oldNeighbourCell:"<<oldNeighbourCell<<endl;
                    Info<<"newNeighborCell:"<<newNeighborCell<<endl;
                    Info<<"newOwnerCell:"<<newOwnerCell<<endl;
                    Info<<"signFace1:"<<signFace1<<endl;
                    Info<<"face1:"<<face1<<endl;
                    Info<<"newOwnerCellPlusOrMinus:"<<newOwnerCellPlusOrMinus<<endl;
                    Info<<"newNeighborCellPlusOrMinus:"<<newNeighborCellPlusOrMinus<<endl;
                    
                    FatalErrorInFunction<< "No neighbor face must be zero face" << exit(FatalError);
                }
                if(signFace1==0)
                {
                    FatalErrorInFunction<<"New cut face can not be a zero face" << exit(FatalError);
                    
                    /*
                    sAUFITB_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesInteriorToBoundary.size(),i));
                    sAUFITB_OldToNewMap.insert(std::pair<label,label>(i,-1));
                    splitAndUnsplitFacesInteriorToBoundary.append(face1);
                    splitAndUnsplitFacesInteriorToBoundaryPatchInd.append(cutCellPatchIndex);
                    splitAndUnsplitFacesInteriorToBoundaryOldInd.append(i);
                    if(newOwnerCellPlusOrMinus==+1 && newNeighborCellPlusOrMinus==-1)
                    {
                        splitAndUnsplitFacesInteriorToBoundaryNeighbor.append(-1);
                        splitAndUnsplitFacesInteriorToBoundaryOwner.append(newOwnerCell);
                    }
                    else if(newOwnerCellPlusOrMinus==-1 && newNeighborCellPlusOrMinus==+1)
                    {
                        splitAndUnsplitFacesInteriorToBoundaryNeighbor.append(-1);
                        splitAndUnsplitFacesInteriorToBoundaryOwner.append(newNeighborCell);
                    }
                    else
                        FatalErrorInFunction<<"Zero face in face does not border a negative and a positive cell" << exit(FatalError);
                    */
                }
                else if(signFace1==2)
                {
                    FatalErrorInFunction<<"Splitted face can not be an all zero face" << exit(FatalError);
                    
                    splitAndUnsplitFacesInteriorToBoundary.append(face1);
                    splitAndUnsplitFacesInteriorToBoundaryPatchInd.append(cutCellPatchIndex);
                    splitAndUnsplitFacesInteriorToBoundaryOldInd.append(-1);
                    
                    if(newOwnerCellPlusOrMinus==+1 && newNeighborCellPlusOrMinus==-1)
                    {
                        splitAndUnsplitFacesInteriorToBoundaryNeighbor.append(-1);
                        splitAndUnsplitFacesInteriorToBoundaryOwner.append(newOwnerCell);
                    }
                    else if(newOwnerCellPlusOrMinus==-1 && newNeighborCellPlusOrMinus==+1)
                    {
                        splitAndUnsplitFacesInteriorToBoundaryNeighbor.append(-1);
                        splitAndUnsplitFacesInteriorToBoundaryOwner.append(newNeighborCell);
                    }
                    else
                        FatalErrorInFunction<<"Zero face in face does not border a negative and a positive cell" << exit(FatalError);
                }
                else
                {
                    //sAUFI_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesInterior.size(),i));
                    //sAUFI_OldToNewMap.insert(std::pair<label,label>(i,-1));
                    splitAndUnsplitFacesInterior.append(face1);
                    splitAndUnsplitFacesInteriorNeighbor.append(newNeighborCell);
                    splitAndUnsplitFacesInteriorOwner.append(newOwnerCell);
                    splitAndUnsplitFacesInteriorPatchInd.append(-1);
                    splitAndUnsplitFacesInteriorOldInd.append(-1);
                }
            }
        }
        if(oldFacesToCutFaces_[i].size()==0)
        {
            //Info<<"GonetoElse"<<endl;
            // Interior uncut face on positive side is appended  without change
            if(facesToSide_[i] == 1)
            {
                //sAUFI_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesInterior.size(),i));
                //sAUFI_OldToNewMap.insert(std::pair<label,label>(i,splitAndUnsplitFacesInterior.size()));
                splitAndUnsplitFacesInterior.append(meshFaces[i]);
                splitAndUnsplitFacesInteriorNeighbor.append(neighbour[i]);
                splitAndUnsplitFacesInteriorOwner.append(owner[i]);
                splitAndUnsplitFacesInteriorPatchInd.append(-1);
                splitAndUnsplitFacesInteriorOldInd.append(i);
            }
            else if(facesToSide_[i] == 0)
            {
                FatalErrorInFunction<<"A face with the side 0 was not cut! "<< exit(FatalError);
                /*
                sAUFITB_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesInteriorToBoundary.size(),i));
                sAUFITB_OldToNewMap.insert(std::pair<label,label>(i,splitAndUnsplitFacesInteriorToBoundary.size()));
                splitAndUnsplitFacesInteriorToBoundary.append(meshFaces[i]);
                splitAndUnsplitFacesInteriorToBoundaryNeighbor.append(-1);
                splitAndUnsplitFacesInteriorToBoundaryPatchInd.append(cutCellPatchIndex);                
                
                label ownerCellSide = cellsToSide_[owner[i]];
                label neighbourCellSide = cellsToSide_[neighbour[i]];
                if(ownerCellSide==+1 && neighbourCellSide==-1)
                    splitAndUnsplitFacesInteriorToBoundaryOwner.append(owner[i]);
                else if(ownerCellSide==-1 && neighbourCellSide==+1)
                    splitAndUnsplitFacesInteriorToBoundaryOwner.append(neighbour[i]);
                else
                {
                    
                    Info<<"meshFaces["<<i<<"]:"<<meshFaces[i]<<Foam::endl;
                    for(label pnt : meshFaces[i])
                        Info<<"pnt:"<<pnt<<" side:"<<pointsToSide_[pnt]<<" dist:"<<pointDist[pnt]<<Foam::endl;
                    FatalErrorInFunction<<"Complete four point zero face must neighbor a negative and a positive cell"<< exit(FatalError);
                }
                */
            }
            else if(facesToSide_[i] == 2)
            {              
                splitAndUnsplitFacesInteriorToBoundary.append(meshFaces[i]);
                splitAndUnsplitFacesInteriorToBoundaryPatchInd.append(cutCellPatchIndex);
                splitAndUnsplitFacesInteriorToBoundaryOldInd.append(i);
                
                if(cellsToSide_[neighbour[i]]==1 && cellsToSide_[owner[i]]==0)
                {
                    splitAndUnsplitFacesInteriorToBoundaryNeighbor.append(-1);
                    splitAndUnsplitFacesInteriorToBoundaryOwner.append(neighbour[i]);
                }
                else if(cellsToSide_[neighbour[i]]==0 && cellsToSide_[owner[i]]==1)
                {
                    splitAndUnsplitFacesInteriorToBoundaryNeighbor.append(-1);
                    splitAndUnsplitFacesInteriorToBoundaryOwner.append(owner[i]);
                }
                else
                    FatalErrorInFunction<<"Cell side of owner or neighbor must be wrong" << exit(FatalError);
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
    }
    //Info<<"Insert split faces boundary"<<endl;
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
                
                if(signFace1==0 || signFace1==2)
                {
                    FatalErrorInFunction<< "New cut boundary face can not have a zero sign and must not be a zero face" << exit(FatalError);
                }

                if(oldFaceToPatchInd[i]<cutCellPatchIndex)
                {
                    //sAUFB_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesBoundary.size(),i));
                    //sAUFB_OldToNewMap.insert(std::pair<label,label>(i,-1));
                    splitAndUnsplitFacesBoundary.append(face1);
                    splitAndUnsplitFacesBoundaryNeighbor.append(-1);
                    splitAndUnsplitFacesBoundaryOwner.append(newOwnerCell);
                    splitAndUnsplitFacesBoundaryPatchInd.append(oldFaceToPatchInd[i]);
                    splitAndUnsplitFacesBoundaryOldInd.append(-1);
                }
                else if(oldFaceToPatchInd[i]>cutCellPatchIndex)
                {
                    //sAUFIf_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesBoundary.size(),i));
                    //sAUFIf_OldToNewMap.insert(std::pair<label,label>(i,-1));
                    splitAndUnsplitFacesInterface.append(face1);
                    splitAndUnsplitFacesInterfaceNeighbor.append(-1);
                    splitAndUnsplitFacesInterfaceOwner.append(newOwnerCell);
                    splitAndUnsplitFacesInterfacePatchInd.append(oldFaceToPatchInd[i]);
                    splitAndUnsplitFacesInterfaceOldInd.append(-1);
                }
                else
                    FatalErrorInFunction<< "Old splitted boundary face can not be a cutCellPatch face" << exit(FatalError);
            }
        }
        //Info<<"Boundary face "<<i;
        if(oldFacesToCutFaces_[i].size()==0)
        {            
            if(facesToSide_[i] == 1)
            {
                if(oldFaceToPatchInd[i]<cutCellPatchIndex)
                {
                    //sAUFB_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesBoundary.size(),i));
                    //sAUFB_OldToNewMap.insert(std::pair<label,label>(i,splitAndUnsplitFacesBoundary.size()));
                    splitAndUnsplitFacesBoundary.append(meshFaces[i]);
                    splitAndUnsplitFacesBoundaryNeighbor.append(-1);
                    splitAndUnsplitFacesBoundaryOwner.append(owner[i]);
                    splitAndUnsplitFacesBoundaryPatchInd.append(oldFaceToPatchInd[i]);
                    splitAndUnsplitFacesBoundaryOldInd.append(i);
                }
                else if(oldFaceToPatchInd[i]>cutCellPatchIndex)
                {
                    //sAUFIf_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesBoundary.size(),i));
                    //sAUFIf_OldToNewMap.insert(std::pair<label,label>(i,splitAndUnsplitFacesBoundary.size()));
                    splitAndUnsplitFacesInterface.append(meshFaces[i]);
                    splitAndUnsplitFacesInterfaceNeighbor.append(-1);
                    splitAndUnsplitFacesInterfaceOwner.append(owner[i]);
                    splitAndUnsplitFacesInterfacePatchInd.append(oldFaceToPatchInd[i]);
                    splitAndUnsplitFacesInterfaceOldInd.append(i);
                }
                else
                    FatalErrorInFunction<< "Old splitted boundary face can not be a cutCellPatch face" << exit(FatalError);
            }
            else if(facesToSide_[i] == 0)
            {
                FatalErrorInFunction
                << "A zero side face can not appear if the face is not cut!"
                << exit(FatalError);
                
                /*
                if(oldFaceToPatchInd[i]<cutCellPatchIndex)
                {
                    sAUFB_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesBoundary.size(),i));
                    sAUFB_OldToNewMap.insert(std::pair<label,label>(i,splitAndUnsplitFacesBoundary.size()));
                    splitAndUnsplitFacesBoundary.append(meshFaces[i]);
                    splitAndUnsplitFacesBoundaryNeighbor.append(-1);
                    splitAndUnsplitFacesBoundaryOwner.append(owner[i]);
                    splitAndUnsplitFacesBoundaryPatchInd.append(oldFaceToPatchInd[i]);
                }
                else if(oldFaceToPatchInd[i]>cutCellPatchIndex)
                {
                    sAUFIf_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesBoundary.size(),i));
                    sAUFIf_OldToNewMap.insert(std::pair<label,label>(i,splitAndUnsplitFacesBoundary.size()));
                    splitAndUnsplitFacesInterface.append(meshFaces[i]);
                    splitAndUnsplitFacesInterfaceNeighbor.append(-1);
                    splitAndUnsplitFacesInterfaceOwner.append(owner[i]);
                    splitAndUnsplitFacesInterfacePatchInd.append(oldFaceToPatchInd[i]);
                }
                else
                    FatalErrorInFunction<< "Old splitted boundary face can not be a cutCellPatch face" << exit(FatalError);
                */
            }
            else if(facesToSide_[i] == 2)
            {
                if(cellsToSide_[owner[i]]==0)
                {
                }
                else if(cellsToSide_[owner[i]]==1)
                {
                    if(oldFaceToPatchInd[i]<cutCellPatchIndex)
                    {
                        //sAUFB_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesBoundary.size(),i));
                        //sAUFB_OldToNewMap.insert(std::pair<label,label>(i,splitAndUnsplitFacesBoundary.size()));
                        splitAndUnsplitFacesBoundary.append(meshFaces[i]);
                        splitAndUnsplitFacesBoundaryNeighbor.append(-1);
                        splitAndUnsplitFacesBoundaryOwner.append(owner[i]);
                        splitAndUnsplitFacesBoundaryPatchInd.append(oldFaceToPatchInd[i]);
                        splitAndUnsplitFacesBoundaryOldInd.append(i);
                    }
                    else if(oldFaceToPatchInd[i]>cutCellPatchIndex)
                    {
                        //sAUFIf_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesBoundary.size(),i));
                        //sAUFIf_OldToNewMap.insert(std::pair<label,label>(i,splitAndUnsplitFacesBoundary.size()));
                        splitAndUnsplitFacesInterface.append(meshFaces[i]);
                        splitAndUnsplitFacesInterfaceNeighbor.append(-1);
                        splitAndUnsplitFacesInterfaceOwner.append(owner[i]);
                        splitAndUnsplitFacesInterfacePatchInd.append(oldFaceToPatchInd[i]);
                        splitAndUnsplitFacesInterfaceOldInd.append(i);
                    }
                    else
                        FatalErrorInFunction<< "Old splitted boundary face can not be a cutCellPatch face" << exit(FatalError);
                }
                else
                    FatalErrorInFunction<<"Cell side of owner or neighbor must be wrong" << exit(FatalError);
                
                /*
                if(oldFaceToPatchInd[i]<cutCellPatchIndex)
                {
                    sAUFB_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesBoundary.size(),i));
                    sAUFB_OldToNewMap.insert(std::pair<label,label>(i,splitAndUnsplitFacesBoundary.size()));
                    splitAndUnsplitFacesBoundary.append(meshFaces[i]);
                    splitAndUnsplitFacesBoundaryNeighbor.append(-1);
                    splitAndUnsplitFacesBoundaryOwner.append(owner[i]);
                    splitAndUnsplitFacesBoundaryPatchInd.append(oldFaceToPatchInd[i]);
                }
                else if(oldFaceToPatchInd[i]>cutCellPatchIndex)
                {
                    sAUFIf_NewToOldMap.insert(std::pair<label,label>(splitAndUnsplitFacesBoundary.size(),i));
                    sAUFIf_OldToNewMap.insert(std::pair<label,label>(i,splitAndUnsplitFacesBoundary.size()));
                    splitAndUnsplitFacesInterface.append(meshFaces[i]);
                    splitAndUnsplitFacesInterfaceNeighbor.append(-1);
                    splitAndUnsplitFacesInterfaceOwner.append(owner[i]);
                    splitAndUnsplitFacesInterfacePatchInd.append(oldFaceToPatchInd[i]);
                }
                else
                    FatalErrorInFunction<< "Old splitted boundary face can not be a cutCellPatch face" << exit(FatalError);
                */
            }
            else if(facesToSide_[i] != -1)
            {
                FatalErrorInFunction
                << "A face with the side: "<<facesToSide_[i]<<" was not treated."
                << " This must not happen."
                << exit(FatalError);
            }
        }    
    }    
    
//Barrier(true);
    
    //reduce for empty cells
    DynamicList<label> cellReductionNumb;
    cellReductionNumb.setSize(deletedCell.size());
    mapOldCellsToNewCells.setSize(meshCells.size());
    //label count = 0;
    label countDel = 0;
    for(int i=0;i<deletedCell.size();i++)
    {
        if(deletedCell[i])
        {
            if(i<meshCells.size())
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
            }
            cellReductionNumb[i] = -1;
            countDel++;
        }
        else
        {
            if(i<meshCells.size())
            {
                if(!(cellToNewPlusCellsIndexes[i].size()==1 && cellToNewMinusCellsIndexes[i].size()==1) &&
                !(cellToNewPlusCellsIndexes[i].size()==1 && cellToNewMinusCellsIndexes[i].size()>1) &&
                !(cellToNewPlusCellsIndexes[i].size()==0 && cellToNewMinusCellsIndexes[i].size()==0))
                {
                    Info<<endl<<"cellToNewPlusCellsIndexes["<<i<<"].size():"<<cellToNewPlusCellsIndexes[i].size()<<endl;
                    Info<<endl<<"cellToNewMinusCellsIndexes["<<i<<"].size():"<<cellToNewMinusCellsIndexes[i].size()<<endl;
                    FatalErrorInFunction<<"Nondeleted Cell. Invalid cell splitting."<<exit(FatalError);
                }
                mapOldCellsToNewCells[i].append(i);
                reverseCellMap[i] = i;
            }
            cellReductionNumb[i] = countDel;
        }
        //count++;
    }
    
    /*
    label maxNewCell = 0;
    for(int i=0;i<mapOldCellsToNewCells.size();i++)
    {
        for(int j=0;j<mapOldCellsToNewCells[i].size();j++)
            maxNewCell = (maxNewCell<mapOldCellsToNewCells[i][j])?mapOldCellsToNewCells[i][j]:maxNewCell;
    }
    */
    
    label maxNumOfNewCellsCorr = 0;
    label maxNumOfNewCellsNonCorr = 0; // added
    for(int i=0;i<mapOldCellsToNewCells.size();i++)
    {
        if((mapOldCellsToNewCells[i].size()==0 && cellReductionNumb[i]!=-1)) //||(mapOldCellsToNewCells[i].size()!=0 && cellReductionNumb[i]==-1))
        {
            Info<<"cellToNewPlusCellsIndexes["<<i<<"]:"<<cellToNewPlusCellsIndexes[i]<<endl;
            Info<<"cellToNewMinusCellsIndexes["<<i<<"]:"<<cellToNewMinusCellsIndexes[i]<<endl;
            Info<<"mapOldCellsToNewCells["<<i<<"]:"<<mapOldCellsToNewCells[i]<<endl;
            Info<<"cellReductionNumb["<<i<<"]:"<<cellReductionNumb[i]<<endl;
            FatalErrorInFunction<<"Face neighbors or ownes deleted cell. This can not happen."<<exit(FatalError);
        }
        for(int j=0;j<mapOldCellsToNewCells[i].size();j++)
        {
            if(mapOldCellsToNewCells[i][j]<i)
                FatalErrorInFunction<<"Old cells can only be mapped to new cells of higher index."<<exit(FatalError);
            if(maxNumOfNewCellsNonCorr<mapOldCellsToNewCells[i][j]) // added
                maxNumOfNewCellsNonCorr = mapOldCellsToNewCells[i][j]; // added
            if(cellReductionNumb[mapOldCellsToNewCells[i][j]]==-1)
                FatalErrorInFunction<<"Positive new cell is deleted. Can not happen"<<exit(FatalError);
            mapOldCellsToNewCells[i][j] -= cellReductionNumb[mapOldCellsToNewCells[i][j]];
            if(maxNumOfNewCellsCorr<mapOldCellsToNewCells[i][j])
                maxNumOfNewCellsCorr = mapOldCellsToNewCells[i][j];
        }
        if(reverseCellMap[i]!=-1)
        {
            reverseCellMap[i] -= cellReductionNumb[reverseCellMap[i]];
        }
    }
    cellMap.setSize(maxNumOfNewCellsCorr,-1);
    for(int oldCelli=0;oldCelli<reverseCellMap.size();oldCelli++)
    {
        if(reverseCellMap[oldCelli]!=-1)
        {
            cellMap[reverseCellMap[oldCelli]] = oldCelli;
        }
    }
    
    for(int i=cellReductionNumb.size();i<maxNumOfNewCellsNonCorr+1;i++)
    {
        cellReductionNumb.append(countDel);
    }
    
    /*
    mapNewCellsToOldCells = labelList(maxNumOfNewCellsCorr+1,-1);
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
    */
    
    /*
    label nbrDeletedCells = 0;
    for(int i=0;i<deletedCell.size();i++)
        if(deletedCell[i])
            nbrDeletedCells++;
    */
    
//Barrier(true);
    
    //correct face owner and neigbour 
    for(int i=0;i<addedCutFaces.size();i++)
    {
        if(cellReductionNumb[addedCutFacesOwner[i]] != -1 &&
           cellReductionNumb[addedCutFacesNeighbor[i]] != -1)
        {
            addedCutFacesOwner[i] -= cellReductionNumb[addedCutFacesOwner[i]];
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
        if(cellReductionNumb[splitAndUnsplitFacesInteriorOwner[i]] != -1 &&
           cellReductionNumb[splitAndUnsplitFacesInteriorNeighbor[i]] != -1)
        {
            if(splitAndUnsplitFacesInteriorNeighbor[i]==-1)
                FatalErrorInFunction<<"Interior face must have a neighbour. This can not happen."<<exit(FatalError);
            
            splitAndUnsplitFacesInteriorOwner[i] -= cellReductionNumb[splitAndUnsplitFacesInteriorOwner[i]];
            splitAndUnsplitFacesInteriorNeighbor[i] -= cellReductionNumb[splitAndUnsplitFacesInteriorNeighbor[i]];
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
        if(cellReductionNumb[splitAndUnsplitFacesBoundaryOwner[i]] != -1)
        {
            splitAndUnsplitFacesBoundaryOwner[i] -= cellReductionNumb[splitAndUnsplitFacesBoundaryOwner[i]];
        }
        else
        {
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
    for(int i=0;i<splitAndUnsplitFacesInteriorToBoundary.size();i++)
    {
        if(cellReductionNumb[splitAndUnsplitFacesInteriorToBoundaryOwner[i]] != -1)
        {
            splitAndUnsplitFacesInteriorToBoundaryOwner[i] -= cellReductionNumb[splitAndUnsplitFacesInteriorToBoundaryOwner[i]];
        }
        else
        {
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
    for(int i=0;i<splitAndUnsplitFacesInterface.size();i++)
    {
        if(cellReductionNumb[splitAndUnsplitFacesInterfaceOwner[i]] != -1)
        {
            splitAndUnsplitFacesInterfaceOwner[i] -= cellReductionNumb[splitAndUnsplitFacesInterfaceOwner[i]];
        }
        else
        {
            FatalErrorInFunction
            << "Face neighbors or ownes deleted cell. This can not happen."
            << exit(FatalError);
        }
    }
        
    if(patchStarts.size()==0)
        FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
    if(splitAndUnsplitFacesBoundaryPatchInd.size()==0)
        FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
    
    DynamicList<label> facesToBoundaryPatchInd;
    facesToBoundaryPatchInd.append(splitAndUnsplitFacesInteriorPatchInd);
    facesToBoundaryPatchInd.append(splitAndUnsplitFacesBoundaryPatchInd);
    facesToBoundaryPatchInd.append(addedCutFacesPatchInd);
    facesToBoundaryPatchInd.append(splitAndUnsplitFacesInteriorToBoundaryPatchInd);
    facesToBoundaryPatchInd.append(splitAndUnsplitFacesInterfacePatchInd);
    
    label interiorFacesOffset = splitAndUnsplitFacesInteriorPatchInd.size();
    label boundaryFacesOffset = interiorFacesOffset+splitAndUnsplitFacesBoundaryPatchInd.size();
    
    std::unordered_set<label> treatedPatches;
    label lastPatch = -1;
    DynamicList<label> boundaryPatchesStart;
    DynamicList<label> interfacePatchesStart;
    for(int i=0;i<facesToBoundaryPatchInd.size();i++)
    {
        if(i<interiorFacesOffset)
        {
            if(facesToBoundaryPatchInd[i]!=-1)
            {
                FatalErrorInFunction<<"Interior face must not be in a patch!"<< exit(FatalError);
            }
        }
        else// if(i<boundaryFacesOffset)
        {
            if(facesToBoundaryPatchInd[i]<lastPatch)
            {
                FatalErrorInFunction<<"Misorganized patches!"<< exit(FatalError);
            }
                
            while(facesToBoundaryPatchInd[i] != lastPatch)
            {
                lastPatch++;
                if(treatedPatches.find(lastPatch)!=treatedPatches.end())
                {
                    FatalErrorInFunction<<"Patch was already treated!"<< exit(FatalError);
                }
                if(Pstream::master())
                    std::cout<<"Patch "<<facesToBoundaryPatchInd[i]<<" start:"<<i<<"  last:"<<lastPatch<<std::endl;
                treatedPatches.insert(lastPatch);
                boundaryPatchesStart.append(i);
            }
        }
    }
    
    /*
    labelList boundaryPatchesStartSIZE(Pstream::nProcs(),0);
    boundaryPatchesStartSIZE[Pstream::myProcNo()] = boundaryPatchesStart.size();
    Pstream::gatherList(boundaryPatchesStartSIZE);
    if(Pstream::master())
        std::cout<<"--------boundaryPatchesStartSIZE:"<<Pstream::myProcNo()<<"----"<<boundaryPatchesStartSIZE[0]<<","<<boundaryPatchesStartSIZE[1]<<","<<boundaryPatchesStartSIZE[2]<<","<<boundaryPatchesStartSIZE[3]<<"----------"<<std::endl;
    
    labelList patchStartsSIZE(Pstream::nProcs(),0);
    patchStartsSIZE[Pstream::myProcNo()] = patchStarts.size();
    Pstream::gatherList(patchStartsSIZE);
    if(Pstream::master())
        std::cout<<"--------patchStartsSIZE:"<<Pstream::myProcNo()<<"----"<<patchStartsSIZE[0]<<","<<patchStartsSIZE[1]<<","<<patchStartsSIZE[2]<<","<<patchStartsSIZE[3]<<"----------"<<std::endl;
    */
    
//Barrier(false);
    
    if(boundaryPatchesStart.size()!=patchStarts.size())
    {
        Info<<"interiorFacesOffset:"<<interiorFacesOffset<<Foam::endl;
        Info<<"boundaryFacesOffset:"<<boundaryFacesOffset<<Foam::endl;
        Info<<"boundaryPatchesStart.size():"<<boundaryPatchesStart.size()<<Foam::endl;
        Info<<"patchStarts.size():"<<patchStarts.size()<<Foam::endl;
        FatalErrorInFunction<<"Patch size does not match!"<< exit(FatalError);
    }
    for(int i=0;i<patchStarts.size();i++)
    {
        if(treatedPatches.find(i)==treatedPatches.end())
            FatalErrorInFunction<<"Patch was not treated!"<< exit(FatalError);
    }
    patchStarts = boundaryPatchesStart;
    if(patchStarts[0] != interiorFacesOffset)
    {
        FatalErrorInFunction<<"Starting patch index does not match"<< exit(FatalError);
    }
    if(patchStarts[cutCellPatchIndex] != boundaryFacesOffset)
    {
        FatalErrorInFunction<<"Starting patch index does not match"<< exit(FatalError);
    }

    for(int i=0;i<patchStarts.size()-1;i++)
    {
        patchSizes[i] = patchStarts[i+1]-patchStarts[i];
        //Info<<"patchSizes["<<i<<"]:"<<patchSizes[i]<<Foam::endl;
    }
    patchSizes.last() = facesToBoundaryPatchInd.size()-patchStarts.last();        
    /*
    Info<<"patchSizes.last():"<<patchSizes.last()<<Foam::endl;
    Info<<"interiorFacesOffset:"<<interiorFacesOffset<<Foam::endl;
    Info<<"boundaryFacesOffset:"<<boundaryFacesOffset<<Foam::endl;
    Info<<"facesToBoundaryPatchInd.size():"<<facesToBoundaryPatchInd.size()<<Foam::endl;
    */
    
    /*
    deletedCellSIZE = labelList(Pstream::nProcs(),0);
    deletedCellSIZE[Pstream::myProcNo()] = deletedCell.size();
    Pstream::gatherList(deletedCellSIZE);
    if(Pstream::master())
        std::cout<<"--------deletedCellSIZE:"<<Pstream::myProcNo()<<"----"<<deletedCellSIZE[0]<<","<<deletedCellSIZE[1]<<","<<deletedCellSIZE[2]<<","<<deletedCellSIZE[3]<<"----------"<<std::endl;
    */
    

//Barrier(true);

    //std::function<faceList&(faceList&)> face_vertice_correction;
    auto face_vertice_correction = [&](faceList& faces)
    {
        for(face& oneFace : faces)
        {
            for(int i=0;i<oneFace.size();i++)
            {
                label reduced_index_vertice = reversePointMap[oneFace[i]];
                if(reduced_index_vertice<0)
                    FatalErrorInFunction<<"Invalid entry in reversePointMap!"<< exit(FatalError);
                oneFace[i] = reduced_index_vertice;
            }
        }
        return faces;
    };

    new_points = newMeshPoints_;
    
    new_faces.append(face_vertice_correction(splitAndUnsplitFacesInterior));
    new_faces.append(face_vertice_correction(splitAndUnsplitFacesBoundary));
    new_faces.append(face_vertice_correction(addedCutFaces));
    new_faces.append(face_vertice_correction(splitAndUnsplitFacesInteriorToBoundary));
    new_faces.append(face_vertice_correction(splitAndUnsplitFacesInterface));
    
    new_owner.append(splitAndUnsplitFacesInteriorOwner);
    new_owner.append(splitAndUnsplitFacesBoundaryOwner);
    new_owner.append(addedCutFacesOwner);
    new_owner.append(splitAndUnsplitFacesInteriorToBoundaryOwner);
    new_owner.append(splitAndUnsplitFacesInterfaceOwner);
    
    new_neighbour.append(splitAndUnsplitFacesInteriorNeighbor);
    new_neighbour.append(splitAndUnsplitFacesBoundaryNeighbor);
    new_neighbour.append(addedCutFacesNeighbor);
    new_neighbour.append(splitAndUnsplitFacesInteriorToBoundaryNeighbor);
    new_neighbour.append(splitAndUnsplitFacesInterfaceNeighbor);
    
    for(int i=0;i<owner.size();i++)
    {
        if(owner[i]<0)
        {
            Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
            Info<<"new_faces["<<i<<"]:"<<new_faces[i]<<endl;
            Info<<"new_owner[i]:"<<new_owner[i]<<endl;
            Info<<"new_neighbour[i]:"<<new_neighbour[i]<<endl;            
            FatalErrorInFunction<<"Owner fail stop"<< exit(FatalError); 
        }
        if(neighbour[i]<-1)
        {
            Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
            Info<<"new_faces["<<i<<"]:"<<new_faces[i]<<endl;
            Info<<"new_owner[i]:"<<new_owner[i]<<endl;
            Info<<"new_neighbour[i]:"<<new_neighbour[i]<<endl;       
            FatalErrorInFunction<<"Neighbour fail stop"<< exit(FatalError); 
        }
    }
    
    faceMap.setSize(0);
    faceMap.append(splitAndUnsplitFacesInteriorOldInd);
    faceMap.append(splitAndUnsplitFacesBoundaryOldInd);
    faceMap.append(labelList(addedCutFaces.size(),-1));
    faceMap.append(splitAndUnsplitFacesInteriorToBoundaryOldInd);
    faceMap.append(splitAndUnsplitFacesInterfaceOldInd);

    for(int newFacei=0;newFacei<faceMap.size();newFacei++)
    {
        if(faceMap[newFacei]!=-1)
        {
            reverseFaceMap[faceMap[newFacei]] = newFacei;
        }
    }

    DynamicList<DynamicList<nurbsReference>> meshPointNurbsReference_new;
    meshPointNurbsReference_new.setSize(new_points.size());
    scalarList pointDist_new;
    pointDist_new.setSize(new_points.size());
    for(int oldPointi=0;oldPointi<meshPointNurbsReference.size();oldPointi++)
    {
        label newPointi = reversePointMap[oldPointi];
        if(newPointi!=-1)
        {                
            meshPointNurbsReference_new[newPointi] = meshPointNurbsReference[oldPointi];
            pointDist_new[newPointi] = pointDist[oldPointi];
        }
    }
    meshPointNurbsReference = meshPointNurbsReference_new;
    pointDist = pointDist_new;
    
    //patchPointMap
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
        Info<<"  Face:"<<i<<" Owner:"<<addedCutFacesOwner[i]<<" Neighbor:"<<addedCutFacesNeighbor[i]<<" ";
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
        Info<<"  Face:"<<i<<" Owner:"<<splitAndUnsplitFacesInteriorOwner[i]<<" Neighbor:"<<splitAndUnsplitFacesInteriorNeighbor[i]<<" ";
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
        Info<<"  Face:"<<i<<" Owner:"<<splitAndUnsplitFacesBoundaryOwner[i]<<" Neighbor:"<<splitAndUnsplitFacesBoundaryNeighbor[i]<<" ";
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

Foam::vector crossProd(const Foam::vector& v1, const Foam::vector& v2)
{
    return Foam::vector(v1[1]*v2[2]-v1[2]*v2[1],
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
        {
            Info<<"meshCells["<<i<<"]:"<<meshCells[i]<<endl;
            Info<<"meshCells["<<i<<"].labels:"<<meshCells[i].labels(meshFaces)<<endl;
            Info<<"meshCells["<<i<<"].points:"<<meshCells[i].points(meshFaces,meshPoints)<<endl;
            FatalErrorInFunction<<"Cell with less than four faces!"<< exit(FatalError);
        }
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
    
    scalar maxCellSize = -1,avgCellSize=0,minCellSize=std::numeric_limits<double>::max();;
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
               crossProdsArea[k+1] >=partialThreeshold*maxFaceSize*(1.f/1e10) &&
               meshFaces[i].size()<=4
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
                Info<<"res:"<<res<<endl;
                
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
               crossProdsArea[k+1] >=partialThreeshold*maxFaceSize*(1.f/1e10) &&
               meshFaces[i].size()<=4
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

/*
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

    scalarList partialVolumeScale = scalarList(newCellVolume.size());
    //Info<<"new Cell Size: "<<newCellVolume.size()<<endl;

    label deletedCellsCount = 0;
    for(int i=0;i<deletedCell.size();i++)
    {
        if(deletedCell[i])
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
    
    
    Info<<endl<<"deletedCellsList: "<<deletedCell<<endl;
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
    
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    Info<< "took \t\t\t\t" << time_span.count() << " seconds."<<endl;
}
*/

void Foam::cutCellFvMesh::agglomerateSmallCells_MC33
(
    scalar partialThreeshold = 0.25
)
{
    std::unordered_map<label,DynamicList<label>> cell_to_faces;
    
    label maxCell=-1;
    for(int faceInd=0;faceInd<new_owner.size();faceInd++)
    {
        maxCell = std::max(maxCell,new_owner[faceInd]);
        cell_to_faces[new_owner[faceInd]].append(faceInd);
    }
    for(int faceInd=0;faceInd<new_neighbour.size();faceInd++)
    {
        maxCell = std::max(maxCell,new_neighbour[faceInd]);
        cell_to_faces[new_neighbour[faceInd]].append(faceInd);
    }
    
    cellList cells(maxCell+1);
    scalarList cellVolume(cells.size());
    for(int cellInd=0;cellInd<cells.size();cellInd++)
    {
        DynamicList<label>& cellFaces = cell_to_faces[cellInd];
        if(cellFaces.size()<4)
        {
            FatalErrorInFunction<<"Illformed cell!"<< exit(FatalError); 
        }
        cellList[cellInd] = cell(cellFaces);
        cellVolume[cellInd] = cellList[cellInd].mag(new_points,new_faces);
    }
    
    
    
    
    
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

    t1 = std::chrono::high_resolution_clock::now();
    scalarList partialVolumeScale = scalarList(newCellVolume.size());

    label deletedCellsCount = 0;
    for(int i=0;i<deletedCell.size();i++)
    {
        if(deletedCell[i])
            deletedCellsCount++;
    }
    if(newCellVolume.size()+deletedCellsCount != deletedCell.size())
    {
        Info<<"deletedCell.size():"<<deletedCell.size()<<endl;
        Info<<"newCellVolume.size():"<<newCellVolume.size()<<endl;
        Info<<"deletedCellsCount:"<<deletedCellsCount<<endl;
        Info<<"oldCellVolume.size():"<<oldCellVolume.size()<<endl;
        FatalErrorInFunction
        << "Must not happen!"
        << exit(FatalError); 
    }
    
    if(newCellVolume.size() != cellMap.size())
    {
        Info<<"newCellVolume.size():"<<newCellVolume.size()<<endl;
        Info<<"cellMap.size():"<<cellMap.size()<<endl;
        FatalErrorInFunction<< "Must not happen!"<< exit(FatalError);
    }
    
    for(int i=0;i<newCellVolume.size();i++)
    {
        if(cellMap[i] == -1)
            partialVolumeScale[i] = 1;
        else
            partialVolumeScale[i] = newCellVolume[i]/oldCellVolume[cellMap[i]];
    }
    
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
        mergeNecessary[i] = false;
        if((partialVolumeScale[i] < 1) && (partialVolumeScale[i] < partialThreeshold))
        {
            mergeNecessary[i] = true;
            
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
    
    for(int i=0;i<newCellVolume.size();i++)
    {
        //Info<<"i:"<<i<<endl;
        //Info<<"possibleMergeFaces[i]:"<<possibleMergeFaces[i]<<endl;
        for(int j=0;j<possibleMergeFaces[i].size();j++)
        {
            //Info<<"j:"<<j<<endl;
            
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
        //Info<<"End"<<endl;
    }
    
    // Sort possible merging cell by respect to face area biggest to smallest
    for(int i=0;i<possibleMergeFaceArea.size();i++)
    {
        //Info<<"Sort: "<<i<<endl;
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
    //label minCellInd = 0;
    scalar maxCellVol = newCells[0].mag(points,faces);
    //label maxCellInd = 0;
    scalar CellVolAvg = 0;

    scalar minCellVolOld = newCells[0].mag(points,faces);
    //label minCellIndOld = 0;
    scalar maxCellVolOld = newCells[0].mag(points,faces);
    //label maxCellIndOld = 0;
    scalar CellVolAvgOld = 0;
    
    scalar vol;
    for(int i=0;i<newCells.size();i++)
    {
        vol = newCells[i].mag(points,faces);
        CellVolAvg += vol;
        if(vol > maxCellVol)
        {
            maxCellVol = vol;
            //maxCellInd = i;
        }
        if(vol < minCellVol)
        {
            minCellVol = vol;
            //minCellInd = i;
        }
    }
    CellVolAvg /= newCells.size();
    
    for(int i=0;i<oldCellVolume.size();i++)
    {
        CellVolAvg += oldCellVolume[i];
        if(oldCellVolume[i] > maxCellVolOld)
        {
            maxCellVolOld = oldCellVolume[i];
            //maxCellIndOld = i;
        }
        if(oldCellVolume[i] < minCellVolOld)
        {
            minCellVolOld = oldCellVolume[i];
            //minCellIndOld = i;
        }
    }
    CellVolAvgOld /= oldCellVolume.size();
    
    for(int i=0;i<newCells.size();i++)
    {
        vol = newCells[i].mag(points,faces); 
        if((vol*factor) < minCellVolOld)
        {
            if(mergeNecessary[i] == false)
            {
                Info<<"i:"<<i<<endl;
                Info<<"newCells["<<i<<"]:"<<newCells[i]<<endl;
                Info<<"minCellVolOld:"<<minCellVolOld<<endl;
                Info<<"vol:"<<vol<<endl;
                Info<<"factor:"<<factor<<endl;
                Info<<"vol*factor:"<<vol*factor<<endl;
                FatalErrorInFunction<<"Not to merge but necessary"<<exit(FatalError);
            }
            if(possibleMergeCells[i].size() == 0)
            {
                FatalErrorInFunction
                << "No merge data available"
                << exit(FatalError);
            }            
        }
    }
    
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
        //label minCellInd = 0;
        scalar maxCellVol = cell[0].mag(point,face);
        //label maxCellInd = 0;
        scalar CellVolAvg = 0;
    
        scalar vol,neighbourVol;
        for(int i=0;i<cell.size();i++)
        {
            vol = cell[i].mag(point,face);
            CellVolAvg += vol;
            if(vol > maxCellVol)
            {
                maxCellVol = vol;
                //maxCellInd = i;
            }
            if(vol < minCellVol)
            {
                minCellVol = vol;
                //minCellInd = i;
            }
        }
        CellVolAvg /= cell.size();
        
        minCellVol = cell[0].mag(point,face);
        //minCellInd = 0;
        maxCellVol = cell[0].mag(point,face);
        //maxCellInd = 0;
        CellVolAvg = 0;
        //label mergeFace,mergeCell;
        for(int i=0;i<cell.size();i++)
        {
            if(mergeFaceOfCell[i][0] < 0)
                continue;

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
                neighbourVol = 0;
                for(int s=0;s<allCells.size();s++)
                {
                    neighbourVol += cell[allCells[s]].mag(points,faces);
                }

                CellVolAvg += neighbourVol+vol;
                if(neighbourVol+vol > maxCellVol)
                {
                    maxCellVol = neighbourVol+vol;
                    //maxCellInd = i;
                }
                if(neighbourVol+vol < minCellVol)
                {
                    minCellVol = neighbourVol+vol;
                    //minCellInd = i;
                }
            }
            else
            {
                CellVolAvg += vol;
                if(vol > maxCellVol)
                {
                    maxCellVol = vol;
                    //maxCellInd = i;
                }
                if(vol < minCellVol)
                {
                    minCellVol = vol;
                    //minCellInd = i;
                }
            }
        }
        CellVolAvg /= cell.size();
    }
//End:TestSection
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
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

/*
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
*/

/*
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
*/

/*
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
*/

/*
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
            Change = false;
        }
    


//Info<<"------"<<count<<"------"<<redIndToCell[count]<<"-------"<<"tryedCells["<<count<<"] = "<<tryedCells[count]<<"/"<<"possibleMergeCells_red["<<count<<"] = "<<possibleMergeCells_red[count].size()<<endl;

        if(mergeNecessary_red[count])

        {
            if(cellReserved.count(redIndToCell[count]) == 0)

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

                {
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
                    label bestTrackBackPoint = trackBackPoints[trackBackPoints.size()-1];

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
                            for(int l=0;l<blockedCells[k].size();l++)
                            {
                                cellReserved.erase(blockedCells[k][l]);
                            }
                        }
                        //Clean list of blockedCells
                        for(int k=backtrackingIndex+1;k<=count;k++)
                        {
                            tryedCells[k] = 0;
                            blockedCells[k].setSize(0);                                    
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
                    }
                    else
                    {
                        Info<<"Backtracking from cell: "<<redIndToCell[count]<<endl;
                        FatalErrorInFunction
                        << " Failed in Merging Selection! There is no found combination for all merging cells."<<endl
                        << exit(FatalError);
                    }
                }
                else
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
                    assignList[redIndToCell[count]] = mergeFace;
                    count++;
                }
            }
            else
            {
                DynamicList<label> temp;
                temp.append(-2);
                assignList[redIndToCell[count]] = temp;
                tryedCells[count] = 0;
                count++;
            }
        }
        else
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
    return assignList;
}
*/

/*
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
        {
            if(cellReserved.count(redIndToCell[count]) == 0 && !cellMergDone_red[count])
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
                {
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
Info<<"Sort: tBP:"<<trackBackPoints<<" tBPMI:"<<trackBackPointsMergeInd<<endl;

                    label bestTrackBackPoint = trackBackPoints[trackBackPoints.size()-1];
                    
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
                            for(int l=0;l<blockedCells[k].size();l++)
                            {
                                cellReserved.erase(blockedCells[k][l]);
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
                    }
                    else
                    {
                        Info<<"Backtracking from cell: "<<redIndToCell[count]<<endl;
                        FatalErrorInFunction
                        << " Failed in Merging Selection! There is no found combination for all merging cells."<<endl
                        << exit(FatalError);
                    }
                }
                else
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
                    count++;
                }
            }
            else
            {
                DynamicList<label> temp;
                temp.append(-2);
                assignList[redIndToCell[count]] = temp;
                tryedCells[count] = 0;
                count++;
            }
        }
        else
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
    return assignList;
}
*/

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
    Info<<"Correct0"<<endl;
    scalar maxFaceSize = -1,avgFaceSize=0,minFaceSize;
    Info<<"Correct01"<<endl;
    if(faces.size()>0) minFaceSize = faces[0].mag(points);
    Info<<"Correct02"<<endl;
    for(int i=0;i<faces.size();i++)
    {
        scalar thisFaceSize = faces[i].mag(points);
        avgFaceSize+=thisFaceSize;
        if(maxFaceSize<thisFaceSize)
            maxFaceSize=thisFaceSize;
        if(minFaceSize>thisFaceSize)
            minFaceSize=thisFaceSize;
    }
    Info<<"Correct03"<<endl;
    avgFaceSize /= faces.size();
    Info<<"owner.size()="<<owner.size()<<endl;
    Info<<"neighbour.size()="<<neighbour.size()<<endl;
    
    label numCellsOwner = 0;
    label numCellsNeighbour = 0;
    for(int i=0;i<owner.size();i++)
    {
        numCellsOwner = std::max(numCellsOwner,owner[i]);
    }
    for(int i=0;i<neighbour.size();i++)
    {
        numCellsNeighbour = std::max(numCellsNeighbour,neighbour[i]);
    }
    numCellsOwner++;
    numCellsNeighbour++;
    Info<<"numCellsOwner:"<<numCellsOwner<<endl;
    Info<<"numCellsNeighbour:"<<numCellsNeighbour<<endl;  

    Info<<"Correct1"<<endl;
    label numCells = std::max(numCellsOwner,numCellsNeighbour);
    List<DynamicList<label>> cellsFaces(numCells);

    Info<<"Correct2::"<<cellsFaces.size()<<endl;

    //Info<<"Correct2"<<endl;
    
    if(owner.size() != neighbour.size())
        FatalErrorInFunction<<"Can not happen"<< exit(FatalError); 
    
    for(int i=0;i<owner.size();i++)
    {
        if(owner[i]<0 || owner[i]>=numCells)
        {
            Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
            Info<<"face["<<i<<"]:"<<faces[i]<<endl;
            Info<<"owner[i]:"<<owner[i]<<endl;
            Info<<"neighbour[i]:"<<neighbour[i]<<endl;            
            FatalErrorInFunction<<"Owner fail stop"<< exit(FatalError); 
        }
        if(neighbour[i]<-1 || neighbour[i]>=numCells)
        {
            Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
            Info<<"face["<<i<<"]:"<<faces[i]<<endl;
            Info<<"owner[i]:"<<owner[i]<<endl;
            Info<<"neighbour[i]:"<<neighbour[i]<<endl;
            Info<<"numCells:"<<numCells<<endl;
            FatalErrorInFunction<<"Neighbour fail stop"<< exit(FatalError); 
        }       
    }
    for(int i=0;i<owner.size();i++)
    {
        cellsFaces[owner[i]].append(i);
        if(neighbour[i]!=-1)
            cellsFaces[neighbour[i]].append(i);
    }
    //FatalErrorInFunction<<"Temporary stop"<< exit(FatalError); 

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

List<bool> Foam::cutCellFvMesh::correctListFaceNormalDir
(
    const pointField& points,
    const faceList& faces,
    const labelList& owner,
    const labelList& neighbour
)
{
    Info<<"Correct0"<<endl;
    scalar maxFaceSize = -1,avgFaceSize=0,minFaceSize;
    Info<<"Correct01"<<endl;
    if(faces.size()>0) minFaceSize = faces[0].mag(points);
    Info<<"Correct02"<<endl;
    for(int i=0;i<faces.size();i++)
    {
        scalar thisFaceSize = faces[i].mag(points);
        avgFaceSize+=thisFaceSize;
        if(maxFaceSize<thisFaceSize)
            maxFaceSize=thisFaceSize;
        if(minFaceSize>thisFaceSize)
            minFaceSize=thisFaceSize;
    }
    Info<<"Correct03"<<endl;
    avgFaceSize /= faces.size();
    Info<<"owner.size()="<<owner.size()<<endl;
    Info<<"neighbour.size()="<<neighbour.size()<<endl;
    
    label numCellsOwner = 0;
    label numCellsNeighbour = 0;
    for(int i=0;i<owner.size();i++)
    {
        numCellsOwner = std::max(numCellsOwner,owner[i]);
    }
    for(int i=0;i<neighbour.size();i++)
    {
        numCellsNeighbour = std::max(numCellsNeighbour,neighbour[i]);
    }
    numCellsOwner++;
    numCellsNeighbour++;
    Info<<"numCellsOwner:"<<numCellsOwner<<endl;
    Info<<"numCellsNeighbour:"<<numCellsNeighbour<<endl;  

    Info<<"Correct1"<<endl;
    label numCells = std::max(numCellsOwner,numCellsNeighbour);
    List<DynamicList<label>> cellsFaces(numCells);

    Info<<"Correct2::"<<cellsFaces.size()<<endl;

    //Info<<"Correct2"<<endl;
    
    if(owner.size() != neighbour.size())
        FatalErrorInFunction<<"Can not happen"<< exit(FatalError); 

    /*
    for(int i=0;i<owner.size();i++)
        Info<<"i:"<<i<<" owner["<<i<<"]:"<<owner[i]<<" neighbour["<<i<<"]:"<<neighbour[i]<<" / "<<neighbour.size()<<" /-/ "<<cellsFaces.size()<<" numCells: "<<numCells<<endl;
    */
    Info<<"owner.size():"<<owner.size()<<endl;
    Info<<"cellsFaces.size():"<<cellsFaces.size()<<endl;
    
    for(int i=0;i<owner.size();i++)
    {
        if(owner[i]<-1 || owner[i]>=numCells)
        {
            Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
            Info<<"face["<<i<<"]:"<<faces[i]<<endl;
            Info<<"owner[i]:"<<owner[i]<<endl;
            Info<<"neighbour[i]:"<<neighbour[i]<<endl;            
            FatalErrorInFunction<<"Owner fail stop"<< exit(FatalError); 
        }
        if(neighbour[i]<-1 || neighbour[i]>=numCells)
        {
            Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
            Info<<"face["<<i<<"]:"<<faces[i]<<endl;
            Info<<"owner[i]:"<<owner[i]<<endl;
            Info<<"neighbour[i]:"<<neighbour[i]<<endl;
            Info<<"numCells:"<<numCells<<endl;
            FatalErrorInFunction<<"Neighbour fail stop"<< exit(FatalError); 
        }
        if(owner[i]<-0 && faces[i].size()!=0)
        {
            Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
            Info<<"face["<<i<<"]:"<<faces[i]<<endl;
            Info<<"owner[i]:"<<owner[i]<<endl;
            Info<<"neighbour[i]:"<<neighbour[i]<<endl;            
            FatalErrorInFunction<<"Owner fail stop"<< exit(FatalError); 
        }
    }
    Info<<"End"<<endl;
    for(int i=0;i<owner.size();i++)
    {
        Info<<"i:"<<i<<" owner[i]:"<<owner[i]<<" neighbour[i]:"<<neighbour[i]<<" / "<<neighbour.size()<<" /-/ "<<cellsFaces.size()<<" numCells: "<<numCells<<endl;
        if(owner[i]!=-1)
            cellsFaces[owner[i]].append(i);
        if(neighbour[i]!=-1)
            cellsFaces[neighbour[i]].append(i);
    }
    //FatalErrorInFunction<<"Temporary stop"<< exit(FatalError); 

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

    List<bool> flipFaces(faces.size(),false);
    
    for(int i=0;i<faces.size();i++)
    {
        if(faces[i].size()==0)
            continue;
        
        //Info<<"Test face "<<i<<endl;
        point centreFace = faces[i].centre(points);
        vector normalFace = faces[i].normal(points);
    
        vector faceCentreToOwnerCentre = cellList[owner[i]].centre(points,faces)-centreFace;
        if((faceCentreToOwnerCentre && normalFace)>=0)
        {
            flipFaces[i] = true;
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
        if(faces[i].size()==0)
            continue;
        
        Info<<"Test face "<<i<<endl;
        point centreFace = faces[i].centre(points);
        vector normalFace = faces[i].normal(points);
        scalar area = faces[i].mag(points);
    
        vector faceCentreToOwnerCentre = cellList[owner[i]].centre(points,faces)-centreFace;
        Info<<"faceCentreToOwnerCentre:"<<faceCentreToOwnerCentre<<endl;
        Info<<"normalFace:"<<normalFace<<endl;
        Info<<"faceCentreToOwnerCentre && normalFace:"<<(faceCentreToOwnerCentre&&normalFace)<<endl;
        if((faceCentreToOwnerCentre && normalFace)>=0 && !flipFaces[i])
        {
            Info<<"In 1"<<endl;
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
            
                if(neighbour[i]>=0)
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
        Info<<"End"<<endl;
    }
    Info<<"Comp End"<<endl;
    
    return flipFaces;
}
