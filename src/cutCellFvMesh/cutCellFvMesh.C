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

Foam::scalar det3x3
(
    const Foam::vector& c0,
    const Foam::vector& c1,
    const Foam::vector& c2
)
{
    return    c0[0]*c1[1]*c2[2]
           +  c1[0]*c2[1]*c0[2]
           +  c2[0]*c0[1]*c1[2]
           -  c0[2]*c1[1]*c2[0]
           -  c1[2]*c2[1]*c0[0]
           -  c2[2]*c0[1]*c1[0];
}

Foam::vector crossProd(const Foam::vector& v1, const Foam::vector& v2)
{
    return Foam::vector(v1[1]*v2[2]-v1[2]*v2[1],
                        v1[2]*v2[0]-v1[0]*v2[2],
                        v1[0]*v2[1]-v1[1]*v2[0]
                        );
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
        else if(lvlSet==0)
            pointsToSide[i] = 2;
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
        if(lvlSet > 0)
            pointsToSide_[i] = 1;
        else if(lvlSet==0)
            pointsToSide_[i] = 2;
        else if(lvlSet < 0)
            pointsToSide_[i] = -1;
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
    
    //Fill MC33Cube lists for every cut cell
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
        if(i==144059)
        {
            Info<<"----------------------------"<<Foam::endl;
            Info<<"oneCube.bitPattern:"<<oneCube.bitPattern<<Foam::endl;
            Info<<"oneCube.vertices:"<<oneCube.vertices<<Foam::endl;
            Info<<"oneCube.cubeCase:"<<oneCube.cubeCase<<Foam::endl;
            for(auto tuple : oneCube.cutTriangles)
            {
                Info<<"("<<std::get<0>(tuple)<<","<<std::get<1>(tuple)<<","<<std::get<2>(tuple)<<")"<<"  ";
            }
            Info<<Foam::endl;
            Info<<"oneCube.cutEdgeVerticeIndex:"<<oneCube.cutEdgeVerticeIndex<<Foam::endl;
            Info<<"----------------------------"<<Foam::endl;
        }
    }

    //Assign global edge Index to MC33Cube cell
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
    
    //Test for complete edge Ind assigment
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
    newMeshPointsInFunc.append(basisPoints);
    pointsToSide(basisPoints);
    
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
                        //Questionable if a valid fail criteria. Has count>2 any other consequences?
                        Info<<"neighborhood_count:"<<neighborhood_count<<Foam::endl;
                        FatalErrorInFunction<<"Triangle can only border maximaly two other triangles!"<< exit(FatalError);
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
    DynamicList<label> provisional_pointMap;
    removedOldPoints.setSize(nOldPoints,true);
    for(label oldInd=0; oldInd<nOldPoints; oldInd++)
    {
        if(pointsToSide_[oldInd]>0)
        {
            provisional_pointMap.append(oldInd);
            removedOldPoints[oldInd] = false;
        }
    }
    label addedPointsStartInd = provisional_pointMap.size();

    reversePointMap.setSize(nOldPoints,-1);
    for(label newInd=0; newInd<provisional_pointMap.size(); newInd++)
    {
        reversePointMap[provisional_pointMap[newInd]] = newInd;
    }

    pointToEgde_.setSize(basisPoints.size(),-1);
    edgeToPoint_.setSize(basisEdges.size(),-1);
    DynamicList<label> nearestPntOfAddedPnt;
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
                
                // Assure each point has an reference point. Does not match the mapPolyMesh
                label refIndex = -1;
                if(pointsToSide_[startLabel] > 0)
                    refIndex=startLabel;
                else if(pointsToSide_[endLabel] > 0)
                    refIndex=endLabel;
                else
                    FatalErrorInFunction<<"No refIndex found!"<< exit(FatalError);

                provisional_pointMap.append(-1);
                if(refIndex<0)
                    FatalErrorInFunction<<"No refIndex found!"<< exit(FatalError);
                nearestPntOfAddedPnt.append(refIndex);
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
                    mc33CutCellData[cellInd].cutVerticeIndexToEdge[cutPointInd] = mc33CubeLocalEdgeInd;
                }
            }
        }
    }
    
    for(int cellInd=0; cellInd<meshCells.size(); cellInd++)
    {
        if(mc33CutCellData[cellInd].cell!=-1)
        {
            MC33::MC33Cube& cube = mc33CutCellData[cellInd];
            std::unordered_set<label> cutEdges;
            for(unsigned int j=0;j<cube.cutTriangles.size();j++)
            {
                auto triangle = cube.cutTriangles[j];
                cutEdges.insert(std::get<0>(triangle));
                cutEdges.insert(std::get<1>(triangle));
                cutEdges.insert(std::get<2>(triangle));
            }
            for(int j=0;j<cube.cutEdgeVerticeIndex.size();j++)
            {
                if(cutEdges.find(j)!=cutEdges.end())
                {
                    if(cube.cutEdgeVerticeIndex[j]==-1)
                        FatalErrorInFunction<<"Missing Vertice assignment!"<< exit(FatalError);
                }
            }
            label countCutEdges = 0;
            vector centerPoint = Foam::zero();
            for(label vertice : cube.cutEdgeVerticeIndex)
            {
                if(vertice!=-1)
                {
                    countCutEdges++;
                    centerPoint += newMeshPointsInFunc[vertice];
                }
            }
            centerPoint /= countCutEdges;
            if(cube.cubeCase==MC33::Case::c73)
            {
                if(countCutEdges!=9)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
            }
            else if(cube.cubeCase==MC33::Case::c102_1262 || cube.cubeCase==MC33::Case::c102_1286)
            {
                if(countCutEdges!=8)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
            }
            else if(cube.cubeCase==MC33::Case::c122_1742 || cube.cubeCase==MC33::Case::c122_1646)
            {
                if(countCutEdges!=8)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
            }
            else if(cube.cubeCase==MC33::Case::c133_MostPos || cube.cubeCase==MC33::Case::c133_MostNeg)
            {
                if(countCutEdges!=9)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
            }
            else if(cube.cubeCase==MC33::Case::c134)
            {
                if(countCutEdges!=12)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
            }
            else
            {
                if(cutEdges.find(12)!=cutEdges.end() ||
                   cutEdges.find(13)!=cutEdges.end() ||
                   cutEdges.find(14)!=cutEdges.end() ||
                   cutEdges.find(15)!=cutEdges.end())
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);                    
                continue;
            }
            
            bool found;
            DynamicList<nurbsReference> reference;
            distToNurbs(centerPoint,found,reference);
            if(!found)
                FatalErrorInFunction<<"New Point have a dist to Nurbs!"<< exit(FatalError);
            
            // Assure each point has an reference point. Does not match the mapPolyMesh
            label refIndex = -1;
            for(label vertice : meshCells[cellInd])
            {
                if(pointsToSide_[vertice] > 0)
                {
                    refIndex=vertice;
                    break;
                }
            }
            if(refIndex==-1)
                FatalErrorInFunction<<"No refIndex found!"<< exit(FatalError);

            provisional_pointMap.append(-1);
            if(refIndex<0)
                FatalErrorInFunction<<"No refIndex found!"<< exit(FatalError);
            nearestPntOfAddedPnt.append(refIndex);
            label centerPointInd = newMeshPointsInFunc.size();
            newMeshPointsInFunc.append(centerPoint);
            meshPointNurbsReference.append(reference);            
            pointsToSide_.append(1);            
            pointToEgde_.append(-1);
            cube.centerPointInd = centerPointInd;
        }
    }
    
    for(int cellInd=0; cellInd<meshCells.size(); cellInd++)
    {
        if(mc33CutCellData[cellInd].cell!=-1)
        {          
            MC33::MC33Cube& cube = mc33CutCellData[cellInd];
            for(auto iter=cube.cutTriangles.begin(); iter!=cube.cutTriangles.end();)
            {
                label localFaceEdgeInd1 = std::get<0>(*iter);
                label cutEdge1VerticeInd = cube.cutEdgeVerticeIndex[localFaceEdgeInd1];
                
                label localFaceEdgeInd2 = std::get<1>(*iter);
                label cutEdge2VerticeInd = cube.cutEdgeVerticeIndex[localFaceEdgeInd2];
                
                label localFaceEdgeInd3 = std::get<2>(*iter);
                label cutEdge3VerticeInd = cube.cutEdgeVerticeIndex[localFaceEdgeInd3];
                
                if(cutEdge1VerticeInd==cutEdge2VerticeInd ||
                   cutEdge2VerticeInd==cutEdge3VerticeInd ||
                   cutEdge3VerticeInd==cutEdge1VerticeInd)
                    iter = cube.cutTriangles.erase(iter);
                else
                    iter++;
            }
        }
    }
    pointMap = provisional_pointMap;
    
    newMeshPoints_ = newMeshPointsInFunc;
    
    purePointMap = pointMap;
    label addedPntIndex=0;
    for(label pointInd=addedPointsStartInd; pointInd<pointMap.size(); pointInd++,addedPntIndex++)
    {
        if(pointMap[pointInd]!=-1)
            FatalErrorInFunction<<"Invalid!"<< exit(FatalError);
        if(!(addedPntIndex<nearestPntOfAddedPnt.size()))
            FatalErrorInFunction<<"Invalid!"<< exit(FatalError);
        pointMap[pointInd] = nearestPntOfAddedPnt[addedPntIndex];
    }
    
    for(label pnt : pointMap)
        if(pnt<0)
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
    
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
    Info<<"basisFaces.size():"<<basisFaces.size()<<Foam::endl;
    for(int faceInd=0;faceInd<basisFaces.size();faceInd++)
    {
        const face& thisFace = basisFaces[faceInd];
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
            const MC33::MC33Cube& thisCellCube = mc33CutCellData[cellInd];
            if(thisCellCube.cell!=-1)
            {
                if(cellInd != thisCellCube.cell)
                {
                    FatalErrorInFunction<<"Cells do not match!"<< exit(FatalError);
                }
                for(unsigned int k=0;k<thisCellCube.cutTriangles.size();k++)
                {
                    label localFaceEdgeInd1 = std::get<0>(thisCellCube.cutTriangles[k]);
                    label faceCutEdge1 = thisCellCube.edgeGlobalInd[localFaceEdgeInd1];
                    label cutEdge1VerticeInd = thisCellCube.cutEdgeVerticeIndex[localFaceEdgeInd1];
                    
                    label localFaceEdgeInd2 = std::get<1>(thisCellCube.cutTriangles[k]);
                    label faceCutEdge2 = thisCellCube.edgeGlobalInd[localFaceEdgeInd2];
                    label cutEdge2VerticeInd = thisCellCube.cutEdgeVerticeIndex[localFaceEdgeInd2];
                    
                    label localFaceEdgeInd3 = std::get<2>(thisCellCube.cutTriangles[k]);
                    label faceCutEdge3 = thisCellCube.edgeGlobalInd[localFaceEdgeInd3];
                    label cutEdge3VerticeInd = thisCellCube.cutEdgeVerticeIndex[localFaceEdgeInd3];
                    
                    if(cutEdge1VerticeInd==cutEdge2VerticeInd ||
                       cutEdge2VerticeInd==cutEdge3VerticeInd ||
                       cutEdge3VerticeInd==cutEdge1VerticeInd)
                    {
                        Info<<"cellInd:"<<cellInd<<Foam::endl;
                        Info<<"k:"<<k<<Foam::endl;
                        Info<<"cutEdge1VerticeInd:"<<cutEdge1VerticeInd<<Foam::endl;
                        Info<<"cutEdge2VerticeInd:"<<cutEdge2VerticeInd<<Foam::endl;
                        Info<<"cutEdge3VerticeInd:"<<cutEdge3VerticeInd<<Foam::endl;
                        FatalErrorInFunction<<"Triangle with duplicate points!"<< exit(FatalError);
                    }
                    
                    face triangleLocEdges(List<label>({localFaceEdgeInd1,localFaceEdgeInd2,localFaceEdgeInd3}));
                    face triangleEdges(List<label>({faceCutEdge1,faceCutEdge2,faceCutEdge3}));
                    face triangleVertices(List<label>({cutEdge1VerticeInd,cutEdge2VerticeInd,cutEdge3VerticeInd}));
                    
                    bool foundCutFaceEdges=false;
                    label which_Triangle_Vertice = -1;
                    label cutEdgeVertice0 = -1;
                    label cutEdgeVertice1 = -1;
                    for(label ip0=0; ip0<triangleLocEdges.size(); ip0++)
                    {
                        bool zeroPnt0InFace = thisFace.which(triangleVertices[ip0])!=-1;
                        bool cutEdgeInFacePnt0 = edgesOfFace.find(triangleEdges[ip0])!=edgesOfFace.end();
                        
                        label ip1 = triangleLocEdges.fcIndex(ip0);
                        bool zeroPnt1InFace = thisFace.which(triangleVertices[ip1])!=-1;
                        bool cutEdgeInFacePnt1 = edgesOfFace.find(triangleEdges[ip1])!=edgesOfFace.end();
                        
                        if((zeroPnt0InFace || cutEdgeInFacePnt0) && (zeroPnt1InFace || cutEdgeInFacePnt1))
                        {
                            if(foundCutFaceEdges)
                                FatalErrorInFunction<<"Multiple edges cut one face!"<< exit(FatalError);
                            foundCutFaceEdges = true;
                            which_Triangle_Vertice = ip0;
                            cutEdgeVertice0 = triangleVertices[ip0];
                            cutEdgeVertice1 = triangleVertices[ip1];
                        }
                    }
                    if(foundCutFaceEdges)
                    {
                        std::pair<std::pair<label,label>,std::pair<label,label>> keyValue =       
                                    {std::pair<label,label>(cutEdgeVertice0,cutEdgeVertice1),
                                     std::pair<label,label>(k,which_Triangle_Vertice)};
                        addedEdgesFromSide[j].insert(keyValue);
                        areThereCutsInFace[j] = true;
                    }
                    
                    if(faceInd==803686)
                    {
                        Info<<"localFaceEdgesInd:("<<localFaceEdgeInd1<<","<<localFaceEdgeInd2<<","<<localFaceEdgeInd3<<")  "<<foundCutFaceEdges<<Foam::endl;
                        Info<<"cutEdge1VerticeInd:"<<cutEdge1VerticeInd<<"   cutEdge2VerticeInd:"<<cutEdge2VerticeInd<<"   cutEdge3VerticeInd:"<<cutEdge3VerticeInd<<Foam::endl;
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
                label triangleNbr = cutEdgeIter0->second.first;
                label verticeNbr = cutEdgeIter0->second.second;
                
                auto& cutTriangle = thisCellCube.cutTriangles[triangleNbr];
                std::vector<label> cutTriangleVec = 
                            {std::get<0>(cutTriangle),std::get<1>(cutTriangle), std::get<2>(cutTriangle)};

                label mc33_loc_edg1 = cutTriangleVec[verticeNbr%3];
                label mc33_loc_edg2 = cutTriangleVec[(verticeNbr+1)%3];
                
                label mc33_edg_vert1 = thisCellCube.cutEdgeVerticeIndex[mc33_loc_edg1];
                label mc33_edg_vert2 = thisCellCube.cutEdgeVerticeIndex[mc33_loc_edg2];
                
                if(!((cutEdgeIter0->first.first==mc33_edg_vert1 || cutEdgeIter0->first.first==mc33_edg_vert2) &&
                    (cutEdgeIter0->first.second==mc33_edg_vert1 || cutEdgeIter0->first.second==mc33_edg_vert2)))
                    FatalErrorInFunction<<"Invalid!"<< exit(FatalError);
                
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
        
        if(faceInd==803686)
        {
            Info<<"cellsOfThisFace.size():"<<cellsOfThisFace.size()<<Foam::endl;
            Info<<"thisFace:"<<thisFace<<Foam::endl;
            Info<<"areThereCutsInFace[0]:"<<areThereCutsInFace[0]<<Foam::endl;
            
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
    
    newMeshEdges_.setCapacity(newMeshEdges_.size());
    edgesToSide_.setCapacity(edgesToSide_.size());
    edgeToFaces_.setCapacity(edgeToFaces_.size());
    faceToEdges_.setCapacity(faceToEdges_.size());
    //edgeToCells_.setCapacity(edgeToCells_.size());
    
    for(label edgInd=basisEdges.size(); edgInd<newMeshEdges_.size(); edgInd++)
    {
        const edge edg = newMeshEdges_[edgInd];
        for(const label* iter=edg.cbegin(); iter!=edg.cend(); iter++)
        {
            label vertice = *iter;
            if(vertice<nOldPoints)
            {
                if(reversePointMap[vertice]==-1)
                {
                    Info<<"vertice:"<<vertice<<" was deleted!"<<Foam::endl;
                    Info<<"pointsToSide_[vertice]:"<<pointsToSide_[vertice]<<Foam::endl;
                    FatalErrorInFunction<<"Error Stop!"<< exit(FatalError);
                }
            }
        }
    }

    Info<<"newMeshEdges_.size():"<<newMeshEdges_.size()<<Foam::endl;
    Info<<"nbrOfPrevEdges:"<<nbrOfPrevEdges<<Foam::endl;    
    //FatalErrorInFunction<<"Temp Stop!"<< exit(FatalError);
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

void Foam::cutCellFvMesh::newMeshFaces_MC33
(
)
{
    const cellList& meshCells = this->cells();
    const faceList& basisFaces = this->faces();
    
    nOldFaces = nbrOfPrevFaces = basisFaces.size();
    
    newMeshFaces_.append(basisFaces);

    facesToSide(newMeshFaces_);
    
    faceToCells_.setSize(basisFaces.size());

    cellToFaces_.setSize(meshCells.size());
    for(int cellInd=0;cellInd<meshCells.size();cellInd++)
    {
        const cell& thisCell = meshCells[cellInd];
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
                
                bool triangleInOneFace = false;
                for(label faceInd : thisCell)
                {
                    const face& oneFace = basisFaces[faceInd];
                    bool vertice0InFace = oneFace.which(vertice0)!=-1;
                    bool vertice1InFace = oneFace.which(vertice1)!=-1;
                    bool vertice2InFace = oneFace.which(vertice2)!=-1;
                    
                    if(vertice0InFace && vertice1InFace && vertice2InFace)
                    {
                        if(triangleInOneFace)
                            FatalErrorInFunction<<"Triangle completely in multiple faces!"<< exit(FatalError);
                        else
                            triangleInOneFace = true;
                    }
                }
                
                if(!triangleInOneFace)
                {
                    newMeshFaces_.append(addedFace);
        
                    facesToSide_.append(+2);
        
                    DynamicList<label> newFaceCell;;
                    newFaceCell.append(cellInd);
                    faceToCells_.append(newFaceCell);
        
                    cellToFaces_[cellInd].append(newMeshFaces_.size()-1);
                }
            }
            if(cellInd==144059)
            {
                Info<<"----------------------------"<<Foam::endl;
                Info<<"thisCellCube.bitPattern:"<<thisCellCube.bitPattern<<Foam::endl;
                Info<<"thisCellCube.vertices:"<<thisCellCube.vertices<<Foam::endl;
                Info<<"thisCellCube.cubeCase:"<<thisCellCube.cubeCase<<Foam::endl;
                for(auto tuple : thisCellCube.cutTriangles)
                {
                    Info<<"("<<std::get<0>(tuple)<<","<<std::get<1>(tuple)<<","<<std::get<2>(tuple)<<")"<<"  ";
                }
                Info<<Foam::endl;
                Info<<"thisCellCube.cutEdgeVerticeIndex:"<<thisCellCube.cutEdgeVerticeIndex<<Foam::endl;
                Info<<"----------------------------"<<Foam::endl;
            }
        }
    }
    
    newMeshFaces_.setCapacity(newMeshFaces_.size());
    faceToCells_.setCapacity(faceToCells_.size());
    cellToFaces_.setCapacity(cellToFaces_.size());
    
    for(label faceInd=basisFaces.size(); faceInd<newMeshFaces_.size(); faceInd++)
    {
        const face oneface = newMeshFaces_[faceInd];
        for(const label* iter=oneface.cbegin(); iter!=oneface.cend(); iter++)
        {
            label vertice = *iter;
            if(vertice<nOldPoints)
            {
                if(reversePointMap[vertice]==-1)
                {
                    Info<<"vertice:"<<vertice<<" was deleted!"<<Foam::endl;    
                    FatalErrorInFunction<<"Error Stop!"<< exit(FatalError);
                }
            }
        }
    }
    Info<<"newMeshFaces_.size():"<<newMeshFaces_.size()<<Foam::endl;
    Info<<"nOldFaces:"<<nOldFaces<<Foam::endl;    
    //FatalErrorInFunction<<"Temp Stop!"<< exit(FatalError);
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

    oldFacesToCutFaces_.setSize(meshFaces.size());
    Info<<"cutFaces_.size():"<<cutFaces_.size()<<Foam::endl;
    Info<<"oldFacesToCutFaces_.size():"<<oldFacesToCutFaces_.size()<<Foam::endl;
    
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
                edge addedEdge = newMeshEdges_[faceToEdges_[i][0]];
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
                        FatalErrorInFunction<<"Can not happen. Edges must be in order"<<exit(FatalError);
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
        
        if(i==803686)
        {
            Info<<"neighborCell:"<<neighborCell<<Foam::endl;
            Info<<"ownerCell:"<<ownerCell<<Foam::endl;
            Info<<"faceToEdges_["<<i<<"]:"<<faceToEdges_[i]<<Foam::endl;
            Info<<"oldFacesToCutFaces_["<<i<<"].size():"<<oldFacesToCutFaces_[i].size()<<Foam::endl;
        }
    }
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
    oldCellVolumesPtr = autoPtr<scalarField>(new scalarField(nOldCells,0));
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
                                else if(pointsToSide_[edgePoints[n]] > 0)
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
                                    if(pointsToSide_[otherPoint] > 0)
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
                       pointsToSide_[cutFace[k]] > 0 && 
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
 
    Info<<"newMeshPoints_.size():"<<newMeshPoints_.size()<<Foam::endl;
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
                
                for(int vertice : addedFace)
                {
                    if(vertice<0 || vertice>=newMeshPoints_.size())
                    {
                        Info<<i<<" "<<j<<" "<<"Vertice "<<vertice<<" does not fit"<<Foam::endl;
                        Info<<"newMeshPoints_.size():"<<newMeshPoints_.size()<<Foam::endl;
                    }
                }
                vector thisNormal = addedFace.normal(newMeshPoints_);
                //Info<<"This Normal: "<<thisNormal<<endl;
                //Info<<"addedFace:"<<addedFace<<Foam::endl;
                
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
                }
                else if(signFace1==2)
                {
                    FatalErrorInFunction<<"Splitted face can not be an all zero face" << exit(FatalError);
                    
                    splitAndUnsplitFacesInteriorToBoundary.append(face1);
                    splitAndUnsplitFacesInteriorToBoundaryPatchInd.append(cutCellPatchIndex);
                    splitAndUnsplitFacesInteriorToBoundaryOldInd.append(-1);
                    splitAndUnsplitFacesInteriorToBoundaryOriginInd.append(i);
                    
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
                    splitAndUnsplitFacesInterior.append(face1);
                    splitAndUnsplitFacesInteriorNeighbor.append(newNeighborCell);
                    splitAndUnsplitFacesInteriorOwner.append(newOwnerCell);
                    splitAndUnsplitFacesInteriorPatchInd.append(-1);
                    splitAndUnsplitFacesInteriorOldInd.append(-1);
                    splitAndUnsplitFacesInteriorOriginInd.append(i);
                }
            }
        }
        if(oldFacesToCutFaces_[i].size()==0)
        {
            //Info<<"GonetoElse"<<endl;
            // Interior uncut face on positive side is appended  without change
            if(facesToSide_[i] == 1)
            {
                splitAndUnsplitFacesInterior.append(meshFaces[i]);
                splitAndUnsplitFacesInteriorNeighbor.append(neighbour[i]);
                splitAndUnsplitFacesInteriorOwner.append(owner[i]);
                splitAndUnsplitFacesInteriorPatchInd.append(-1);
                splitAndUnsplitFacesInteriorOldInd.append(i);
                splitAndUnsplitFacesInteriorOriginInd.append(i);
            }
            else if(facesToSide_[i] == 0)
            {
                FatalErrorInFunction<<"A face with the side 0 was not cut! "<< exit(FatalError);
            }
            else if(facesToSide_[i] == 2)
            {              
                splitAndUnsplitFacesInteriorToBoundary.append(meshFaces[i]);
                splitAndUnsplitFacesInteriorToBoundaryPatchInd.append(cutCellPatchIndex);
                splitAndUnsplitFacesInteriorToBoundaryOldInd.append(i);
                splitAndUnsplitFacesInteriorToBoundaryOriginInd.append(i);
                
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
                    splitAndUnsplitFacesBoundaryOriginInd.append(i);
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
                    splitAndUnsplitFacesInterfaceOriginInd.append(i);
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
                    splitAndUnsplitFacesBoundary.append(meshFaces[i]);
                    splitAndUnsplitFacesBoundaryNeighbor.append(-1);
                    splitAndUnsplitFacesBoundaryOwner.append(owner[i]);
                    splitAndUnsplitFacesBoundaryPatchInd.append(oldFaceToPatchInd[i]);
                    splitAndUnsplitFacesBoundaryOldInd.append(i);
                    splitAndUnsplitFacesBoundaryOriginInd.append(i);
                }
                else if(oldFaceToPatchInd[i]>cutCellPatchIndex)
                {
                    splitAndUnsplitFacesInterface.append(meshFaces[i]);
                    splitAndUnsplitFacesInterfaceNeighbor.append(-1);
                    splitAndUnsplitFacesInterfaceOwner.append(owner[i]);
                    splitAndUnsplitFacesInterfacePatchInd.append(oldFaceToPatchInd[i]);
                    splitAndUnsplitFacesInterfaceOldInd.append(i);
                    splitAndUnsplitFacesInterfaceOriginInd.append(i);
                }
                else
                    FatalErrorInFunction<< "Old splitted boundary face can not be a cutCellPatch face" << exit(FatalError);
            }
            else if(facesToSide_[i] == 0)
            {
                Info<<"oldFacesToCutFaces_["<<i<<"].size():"<<oldFacesToCutFaces_[i].size()<<Foam::endl;
                Info<<"meshFaces["<<i<<"]:"<<meshFaces[i]<<Foam::endl;
                for(label vert : meshFaces[i])
                    Info<<"("<<pointsToSide_[vert]<<","<<pointDist[vert]<<")"<<Foam::endl;
                FatalErrorInFunction
                << "A zero side face can not appear if the face is not cut!"
                << exit(FatalError);
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
                        splitAndUnsplitFacesBoundary.append(meshFaces[i]);
                        splitAndUnsplitFacesBoundaryNeighbor.append(-1);
                        splitAndUnsplitFacesBoundaryOwner.append(owner[i]);
                        splitAndUnsplitFacesBoundaryPatchInd.append(oldFaceToPatchInd[i]);
                        splitAndUnsplitFacesBoundaryOldInd.append(i);
                        splitAndUnsplitFacesBoundaryOriginInd.append(i);
                    }
                    else if(oldFaceToPatchInd[i]>cutCellPatchIndex)
                    {
                        splitAndUnsplitFacesInterface.append(meshFaces[i]);
                        splitAndUnsplitFacesInterfaceNeighbor.append(-1);
                        splitAndUnsplitFacesInterfaceOwner.append(owner[i]);
                        splitAndUnsplitFacesInterfacePatchInd.append(oldFaceToPatchInd[i]);
                        splitAndUnsplitFacesInterfaceOldInd.append(i);
                        splitAndUnsplitFacesInterfaceOriginInd.append(i);
                    }
                    else
                        FatalErrorInFunction<< "Old splitted boundary face can not be a cutCellPatch face" << exit(FatalError);
                }
                else
                    FatalErrorInFunction<<"Cell side of owner or neighbor must be wrong" << exit(FatalError);
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
    cellMap.setSize(maxNumOfNewCellsCorr+1,-1);
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
    }
    patchSizes.last() = facesToBoundaryPatchInd.size()-patchStarts.last();

    //std::function<faceList&(faceList&)> face_vertice_correction;
    oldToNewPointInd = labelList(newMeshPointsInFunc.size(),-1);
    label maxOldPntNewIndex = 0;
    for(label oldPnti=0;oldPnti<reversePointMap.size();oldPnti++)
    {
        if(reversePointMap[oldPnti]!=-1)
        {
            oldToNewPointInd[oldPnti] = reversePointMap[oldPnti];
            maxOldPntNewIndex = std::max(reversePointMap[oldPnti],maxOldPntNewIndex);
        }
    }
    for(label newPnti=reversePointMap.size();newPnti<oldToNewPointInd.size();newPnti++)
    {
        maxOldPntNewIndex++;
        oldToNewPointInd[newPnti] = maxOldPntNewIndex;
    }

    new_points = pointField(maxOldPntNewIndex+1);
    for(label oldPnti=0;oldPnti<oldToNewPointInd.size();oldPnti++)
    {
        label newPnti = oldToNewPointInd[oldPnti];
        if(newPnti!=-1)
        {
            new_points[newPnti] = newMeshPoints_[oldPnti];
        }
    }
 
    auto face_vertice_correction = [&](faceList& faces)
    {
        for(face& oneFace : faces)
        {
            for(int i=0;i<oneFace.size();i++)
            {
                label reduced_index_vertice = oldToNewPointInd[oneFace[i]];
                if(reduced_index_vertice<0)
                {
                    Info<<"maxOldPntNewIndex:"<<maxOldPntNewIndex<<Foam::endl;
                    Info<<"newMeshPointsInFunc.size():"<<newMeshPointsInFunc.size()<<Foam::endl;
                    Info<<"this->points().size():"<<this->points().size()<<Foam::endl;
                    Info<<"nOldPoints:"<<nOldPoints<<Foam::endl;
                    Info<<"new_points.size():"<<new_points.size()<<Foam::endl;
                    Info<<"pointMap.size():"<<pointMap.size()<<Foam::endl;
                    Info<<"reversePointMap.size():"<<reversePointMap.size()<<Foam::endl;
                    Info<<"reduced_index_vertice:"<<reduced_index_vertice<<Foam::endl;
                    Info<<"pointsToSide_["<<oneFace[i]<<"]:"<<pointsToSide_[oneFace[i]]<<Foam::endl;
                    FatalErrorInFunction<<"Invalid entry in reversePointMap!"<< exit(FatalError);
                }
                oneFace[i] = reduced_index_vertice;
            }
        }
        return faces;
    };

    
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
        
    for(int i=0;i<new_owner.size();i++)
    {
        if(new_owner[i]<0)
        {
            Info<<"nbrOfPrevFaces:"<<nbrOfPrevFaces<<endl;
            Info<<"new_faces["<<i<<"]:"<<new_faces[i]<<endl;
            Info<<"new_owner[i]:"<<new_owner[i]<<endl;
            Info<<"new_neighbour[i]:"<<new_neighbour[i]<<endl;            
            FatalErrorInFunction<<"Owner fail stop"<< exit(FatalError); 
        }
        if(new_neighbour[i]<-1)
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
    
    DynamicList<label> faceToOriginInd;
    faceToOriginInd.append(splitAndUnsplitFacesInteriorOriginInd);
    faceToOriginInd.append(splitAndUnsplitFacesBoundaryOriginInd);
    faceToOriginInd.append(labelList(addedCutFaces.size(),-1));
    faceToOriginInd.append(splitAndUnsplitFacesInteriorToBoundaryOriginInd);
    faceToOriginInd.append(splitAndUnsplitFacesInterfaceOriginInd);
    if(faceToOriginInd.size()!=faceMap.size())
        FatalErrorInFunction<<"Error!"<< exit(FatalError);
    for(int newFacei=0; newFacei<faceToOriginInd.size(); newFacei++)
    {
        newFaceToOldFaceInds[newFacei] = faceToOriginInd[newFacei];
        oldFaceToNewFaceInds[faceToOriginInd[newFacei]].append(newFacei)
    }
    
    DynamicList<DynamicList<nurbsReference>> meshPointNurbsReference_new;
    meshPointNurbsReference_new.setSize(new_points.size());
    scalarList pointDist_new;
    pointDist_new.setSize(new_points.size());
    Info<<"oldToNewPointInd.size():"<<oldToNewPointInd.size()<<Foam::endl;
    Info<<"pointDist_new.size():"<<pointDist_new.size()<<Foam::endl;
    Info<<"pointDist.size():"<<pointDist.size()<<Foam::endl;
    Info<<"nOldPoints:"<<nOldPoints<<Foam::endl;
    Info<<"meshPointNurbsReference.size():"<<meshPointNurbsReference.size()<<Foam::endl;
    for(int oldPointi=0;oldPointi<meshPointNurbsReference.size();oldPointi++)
    {
        label newPointi = oldToNewPointInd[oldPointi];
        if(newPointi>=new_points.size())
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
        if(newPointi!=-1)
        {
            meshPointNurbsReference_new[newPointi] = meshPointNurbsReference[oldPointi];
            pointDist_new[newPointi] = pointDist[oldPointi];
        }
    }
    meshPointNurbsReference = meshPointNurbsReference_new;
    pointDist = pointDist_new;
    
    Info<<"End"<<Foam::endl;
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
    if(meshCells.size()>0) 
        minCellSize = meshCells[0].mag(meshPoints,meshFaces);
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
Barrier(true);

    
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

void Foam::cutCellFvMesh::checkPermutation()
{
    for(unsigned int i=0;i<permutationTable.size();i++)
    {
        for(unsigned int j=0;j<permutationTable.size();j++)
        {
            if(i!=j)
            {
                if(permutationTable[i]==permutationTable[j])
                {
                    FatalErrorInFunction<<"Illformed permutation Table!"<< exit(FatalError);
                }
            }
        }
    }
    
    List<std::unordered_set<label>> verticeToNeighbors(8);
    verticeToNeighbors[0].insert(1);
    verticeToNeighbors[0].insert(3);
    verticeToNeighbors[0].insert(4);
    
    verticeToNeighbors[1].insert(0);
    verticeToNeighbors[1].insert(2);
    verticeToNeighbors[1].insert(5);
    
    verticeToNeighbors[2].insert(1);
    verticeToNeighbors[2].insert(3);
    verticeToNeighbors[2].insert(6);
    
    verticeToNeighbors[3].insert(0);
    verticeToNeighbors[3].insert(2);
    verticeToNeighbors[3].insert(7);
    
    verticeToNeighbors[4].insert(0);
    verticeToNeighbors[4].insert(5);
    verticeToNeighbors[4].insert(7);
    
    verticeToNeighbors[5].insert(1);
    verticeToNeighbors[5].insert(4);
    verticeToNeighbors[5].insert(6);
    
    verticeToNeighbors[6].insert(2);
    verticeToNeighbors[6].insert(5);
    verticeToNeighbors[6].insert(7);
    
    verticeToNeighbors[7].insert(3);
    verticeToNeighbors[7].insert(4);
    verticeToNeighbors[7].insert(6);
    
    for(unsigned int i=0;i<permutationTable.size();i++)
    {
        label pnt0 = permutationTable[i][0];
        label pnt1 = permutationTable[i][1];
        label pnt2 = permutationTable[i][2];
        label pnt3 = permutationTable[i][3];
        label pnt4 = permutationTable[i][4];
        label pnt5 = permutationTable[i][5];
        label pnt6 = permutationTable[i][6];
        label pnt7 = permutationTable[i][7];
        
        // Edge 0,3,8
        if(verticeToNeighbors[pnt0].find(pnt1)!=verticeToNeighbors[pnt0].end() ||
           verticeToNeighbors[pnt0].find(pnt3)!=verticeToNeighbors[pnt0].end() ||
           verticeToNeighbors[pnt0].find(pnt4)!=verticeToNeighbors[pnt0].end() )
            FatalErrorInFunction<<"Illformed permutation Table!"<< exit(FatalError);
            
        // Edge 0,1,9
        if(verticeToNeighbors[pnt1].find(pnt0)!=verticeToNeighbors[pnt1].end() ||
           verticeToNeighbors[pnt1].find(pnt2)!=verticeToNeighbors[pnt1].end() ||
           verticeToNeighbors[pnt1].find(pnt5)!=verticeToNeighbors[pnt1].end() )
            FatalErrorInFunction<<"Illformed permutation Table!"<< exit(FatalError);
            
        // Edge 1,2,A
        if(verticeToNeighbors[pnt2].find(pnt1)!=verticeToNeighbors[pnt2].end() ||
           verticeToNeighbors[pnt2].find(pnt3)!=verticeToNeighbors[pnt2].end() ||
           verticeToNeighbors[pnt2].find(pnt6)!=verticeToNeighbors[pnt2].end() )
            FatalErrorInFunction<<"Illformed permutation Table!"<< exit(FatalError);
            
        // Edge 3,2,B
        if(verticeToNeighbors[pnt3].find(pnt0)!=verticeToNeighbors[pnt3].end() ||
           verticeToNeighbors[pnt3].find(pnt2)!=verticeToNeighbors[pnt3].end() ||
           verticeToNeighbors[pnt3].find(pnt7)!=verticeToNeighbors[pnt3].end() )
            FatalErrorInFunction<<"Illformed permutation Table!"<< exit(FatalError);
            
        // Edge 8,4,7
        if(verticeToNeighbors[pnt4].find(pnt0)!=verticeToNeighbors[pnt4].end() ||
           verticeToNeighbors[pnt4].find(pnt5)!=verticeToNeighbors[pnt4].end() ||
           verticeToNeighbors[pnt4].find(pnt7)!=verticeToNeighbors[pnt4].end() )
            FatalErrorInFunction<<"Illformed permutation Table!"<< exit(FatalError);
            
        // Edge 9,4,5
        if(verticeToNeighbors[pnt5].find(pnt1)!=verticeToNeighbors[pnt5].end() ||
           verticeToNeighbors[pnt5].find(pnt4)!=verticeToNeighbors[pnt5].end() ||
           verticeToNeighbors[pnt5].find(pnt6)!=verticeToNeighbors[pnt5].end() )
            FatalErrorInFunction<<"Illformed permutation Table!"<< exit(FatalError);
            
        // Edge A,5,6
        if(verticeToNeighbors[pnt6].find(pnt2)!=verticeToNeighbors[pnt6].end() ||
           verticeToNeighbors[pnt6].find(pnt5)!=verticeToNeighbors[pnt6].end() ||
           verticeToNeighbors[pnt6].find(pnt7)!=verticeToNeighbors[pnt6].end() )
            FatalErrorInFunction<<"Illformed permutation Table!"<< exit(FatalError);
            
        // Edge B,7,6
        if(verticeToNeighbors[pnt7].find(pnt3)!=verticeToNeighbors[pnt7].end() ||
           verticeToNeighbors[pnt7].find(pnt4)!=verticeToNeighbors[pnt7].end() ||
           verticeToNeighbors[pnt7].find(pnt6)!=verticeToNeighbors[pnt7].end() )
            FatalErrorInFunction<<"Illformed permutation Table!"<< exit(FatalError);       
    }
    
}

void Foam::cutCellFvMesh::triangulateCell
(
    const cell& thisCell
)
{
    DynamicList<label> remCell;
    DynamicList<face> faces;
    DynamicList<std::pair<label,bool>> faceMap;
    DynamicList<label> owner;
    DynamicList<label> neighbour;
    for(label locFaceInd=0; locFaceInd<thisCell.size(); locFaceInd++)
    {
        label faceInd = thisCell[locFaceInd];
        remCell.append(locFaceInd);
        faces.append(new_faces[faceInd]);
        faceMap.append({faceInd,true});
        owner.append(new_owner[faceInd]);
        neighbour.append(new_neighbour[faceInd]);
    }
    
    while(remCell.size()>0)
    {
        //Collect point to ordered face mapping
        std::unordered_map<label,DynamicList<label>> pointToFaces;
        for(label locFaceInd=0; locFaceInd<remCell.size(); locFaceInd++)
        {
            const face& thisFace = faces[remCell[locFaceInd]];
            for(const label pointInd : thisFace)
            {
                pointToFaces[pointInd].append(remCell[locFaceInd]);
            }
        }
        //Order faces in point to face mapping
        for(auto iter=pointToFaces.begin(); iter!=pointToFaces.end(); iter++)
        {
            label point = iter->first;
            DynamicList<label>& pointFaces = iter->second;
            label nbrFaces = pointFaces.size();
            std::list<label> pointFacesLis(pointFaces.begin(),pointFaces.end());
            pointFaces.clear();
            label faceFrontInd = pointFacesLis.front();
            pointFacesLis.pop_front();
            pointFaces.append(faceFrontInd);
            while(!pointFacesLis.empty())
            {
                const face& faceFront = faces[faceFrontInd];
                label locPntInd = faceFront.which(point);
                if(locPntInd==-1)
                    FatalErrorInFunction<<"Error"<< exit(FatalError);
                label nextPntInd = faceFront.nextLabel(locPntInd);
                label prevPntInd = faceFront.prevLabel(locPntInd);
                for(auto iter=pointFacesLis.cbegin(); ; )
                {
                    const face& faceCon = faces[*iter];
                    if(faceCon.which(point)==-1)
                        FatalErrorInFunction<<"Error"<< exit(FatalError);
                    bool nextPntInFaceCon = (faceCon.which(nextPntInd)!=-1);
                    bool prevPntInFaceCon = (faceCon.which(prevPntInd)!=-1);
                    
                    if((nextPntInFaceCon && !prevPntInFaceCon) ||
                       (!nextPntInFaceCon && prevPntInFaceCon))
                    {
                        pointFaces.append(*iter);
                        faceFrontInd = *iter;
                        pointFacesLis.erase(iter);
                        break;
                    }
                    else if(!nextPntInFaceCon && !prevPntInFaceCon)
                    {
                        iter++;
                        if(iter==pointFacesLis.end())
                            FatalErrorInFunction<<"Error"<< exit(FatalError);
                    }
                    else
                    {
                        Info<<"nextPntInFaceCon:"<<nextPntInFaceCon<<Foam::endl;
                        Info<<"prevPntInFaceCon:"<<prevPntInFaceCon<<Foam::endl;
                        FatalErrorInFunction<<"Error"<< exit(FatalError);
                    }
                }
            }
            const face& firstFace = faces[pointFaces.first()];
            const face& lastFace = faces[pointFaces.last()];
            label locPntIndfirstFace = firstFace.which(point);
            if(locPntIndfirstFace==-1)
                FatalErrorInFunction<<"Error"<< exit(FatalError);
            label locPntIndlastFace = lastFace.which(point);
            if(locPntIndlastFace==-1)
                FatalErrorInFunction<<"Error"<< exit(FatalError);
            label nextPntIndfirstFace = firstFace.nextLabel(locPntIndfirstFace);
            label prevPntIndfirstFace = firstFace.prevLabel(locPntIndfirstFace);
            label nextPntIndlastFace = lastFace.nextLabel(locPntIndlastFace);
            label prevPntIndlastFace = lastFace.prevLabel(locPntIndlastFace);
            if(!((nextPntIndfirstFace==prevPntIndlastFace) ||
                 (nextPntIndfirstFace==nextPntIndlastFace) ||
                 (prevPntIndfirstFace==prevPntIndlastFace) ||
                 (prevPntIndfirstFace==nextPntIndlastFace) ))
                FatalErrorInFunction<<"Error"<< exit(FatalError);
            if(nbrFaces!=pointFaces.size())
                FatalErrorInFunction<<"Error"<< exit(FatalError);

        }
        
        std::unordered_map<label,scalar> pointToMinAngle;
        break;
    }
}

std::unique_ptr<List<scalar>> Foam::cutCellFvMesh::faceToCellCenterRelation
(
    const cell& thisCell,
    const label cellInd
)
{
    vector thisCellCentre = thisCell.centre(new_points,new_faces);
    std::unique_ptr<List<scalar>>
    faceToCenterAngles(new List<scalar>());
    
    for(const label faceInd : thisCell)
    {
        const face& cellFace = new_faces[faceInd];
        vector normal = cellFace.normal(new_points);
        if(new_neighbour[faceInd]!=-1)
        {
            if(new_owner[faceInd]==cellInd)
            {}
            else if(new_neighbour[faceInd]==cellInd)
            {
                normal *= -1;
            }
            else
                FatalErrorInFunction<<"Error in face!"<< exit(FatalError);
        }
        vector faceCentre = cellFace.centre(new_points);
        vector centreToFaceCentre = faceCentre-thisCellCentre;
        
        scalar cosAngle = centreToFaceCentre & normal;
        cosAngle /= (norm2(centreToFaceCentre) * norm2(normal));
        faceToCenterAngles->append(cosAngle);
    }
    return faceToCenterAngles;
}

std::unique_ptr<List<DynamicList<std::pair<label,std::pair<label,label>>>>> Foam::cutCellFvMesh::compCellFaceEdgeGraph
(
    const cell& oneCell                                                                             
)
{
    std::unique_ptr<List<DynamicList<std::pair<label,std::pair<label,label>>>>> cellFaceEdgeGraph(new List<DynamicList<std::pair<label,std::pair<label,label>>>>(oneCell.size()));
    for(label locFaceInd=0;locFaceInd<oneCell.size();locFaceInd++)
    {
        const label oneFaceInd = oneCell[locFaceInd];
        const face& oneFace = new_faces[oneFaceInd];
        for(const label pntIni : oneFace)
        {
            label locVertInd = oneFace.which(pntIni);
            if(locVertInd==-1)
                FatalErrorInFunction<<"Error"<< exit(FatalError);
            label pntNext = oneFace.nextLabel(locVertInd);
            
            // Find the face connected to the two points pntIni,pntNext
            label j_Conn = -1;
            label smConnPnt = -1;
            label laConnPnt = -1;
            for(label j=0;j<oneCell.size();j++)
            {
                if(locFaceInd!=j)
                {
                    const label otherFaceInd = oneCell[j];
                    const face& otherFace = new_faces[otherFaceInd];
                    label otherLocVertInd = otherFace.which(pntIni);
                    if(otherLocVertInd!=-1)
                    {
                        label otherPntNext = otherFace.nextLabel(otherLocVertInd);
                        label otherPntPrev = otherFace.prevLabel(otherLocVertInd);
                        if(otherPntNext==pntNext || otherPntPrev==pntNext)
                        {
                            if(j_Conn!=-1)
                            {
                                Info<<"j_Conn:"<<j_Conn<<Foam::endl;
                                Info<<"j:"<<j<<Foam::endl;
                                FatalErrorInFunction<<"Error"<< exit(FatalError);
                            }
                            j_Conn = j;
                            if(pntIni<pntNext)
                            {
                                smConnPnt = pntIni;
                                laConnPnt = pntNext;
                            }
                            else if(pntIni>pntNext)
                            {
                                smConnPnt = pntNext;
                                laConnPnt = pntIni;
                            }
                            else
                                FatalErrorInFunction<<"Error"<< exit(FatalError);
                        }
                    }
                }
            }
            if(j_Conn==-1)
                FatalErrorInFunction<<"Error"<< exit(FatalError);
            (*cellFaceEdgeGraph)[locFaceInd].append({j_Conn,{smConnPnt,laConnPnt}});
        }
        if((*cellFaceEdgeGraph)[locFaceInd].size()!=oneFace.size())
            FatalErrorInFunction<<"Error"<< exit(FatalError);
    }
    return cellFaceEdgeGraph;
}

std::unique_ptr<List<List<DynamicList<label>>>> Foam::cutCellFvMesh::compCellFacePointGraph
(
    const cell& oneCell                                                                             
)
{
    std::unique_ptr<List<List<DynamicList<label>>>> cellFacePointGraph(new List<List<DynamicList<label>>>(oneCell.size()));
    for(label locFaceInd=0;locFaceInd<oneCell.size();locFaceInd++)
    {
        const label oneFaceInd = oneCell[locFaceInd];
        const face& oneFace = new_faces[oneFaceInd];
        (*cellFacePointGraph)[locFaceInd].setSize(oneFace.size());
        for(label locPntInd=0;locPntInd<oneFace.size();locPntInd++)
        {
            label pntIni = oneFace[locPntInd];
            label locVertInd = oneFace.which(pntIni);
            if(locVertInd==-1)
                FatalErrorInFunction<<"Error"<< exit(FatalError);
            
            for(label j=0;j<oneCell.size();j++)
            {
                if(locFaceInd!=j)
                {
                    const label otherFaceInd = oneCell[j];
                    const face& otherFace = new_faces[otherFaceInd];
                    label otherLocVertInd = otherFace.which(pntIni);
                    if(otherLocVertInd!=-1)
                    {
                        (*cellFacePointGraph)[locFaceInd][locPntInd].append(j);
                    }
                }
            }
            if((*cellFacePointGraph)[locFaceInd][locPntInd].size()<2)
                FatalErrorInFunction<<"Error"<< exit(FatalError);
        }
    }
    return cellFacePointGraph;
}

bool Foam::cutCellFvMesh::testCorrectFaceSizes
(
    const std::array<label,10> faceSizeNbrReq,
    const cell& thisCell,
    bool print
)
{
    if(faceSizeNbrReq[0]!=0 || faceSizeNbrReq[1]!=0 || faceSizeNbrReq[2]!=0)
        FatalErrorInFunction<<"Error!"<< exit(FatalError);
    
    std::array<label,10> faceSizeNbrGiven;
    faceSizeNbrGiven.fill(0);
    for(const label faceInd : thisCell)
    {
        const face& oneFace = new_faces[faceInd];
        if(static_cast<unsigned int>(oneFace.size())>=faceSizeNbrGiven.size())
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
        faceSizeNbrGiven[oneFace.size()]++;
    }
    if(print)
    {
        Info<<"{";
        for(label nbr : faceSizeNbrGiven)
            Info<<" "<<nbr;
        Info<<" }"<<Foam::endl;
    }
    return faceSizeNbrGiven==faceSizeNbrReq;
}

label Foam::cutCellFvMesh::memSecMC33VerticePntInd
(
    label mc33LocalPntInd_I,
    const MC33::MC33Cube& mc33cube,
    const std::unordered_set<label>& verticeOfCell
)
{
    if(mc33LocalPntInd_I<0 || mc33LocalPntInd_I>=8)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    label pntOldInd_I = mc33cube.vertices[mc33LocalPntInd_I];
    if(pntOldInd_I==-1 || pntOldInd_I>=oldToNewPointInd.size())
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    label pntInd_I = oldToNewPointInd[pntOldInd_I];
    if(pntInd_I==-1 || verticeOfCell.find(pntInd_I)==verticeOfCell.end() || pntInd_I>=new_points.size())
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    return pntInd_I;
}

label Foam::cutCellFvMesh::memSecMC33CutEdgePntInd
(
    label mc33LocalEdgeInd_I,
    const MC33::MC33Cube& mc33cube,
    const std::unordered_set<label>& verticeOfCell
)
{
    if(mc33LocalEdgeInd_I<0 || mc33LocalEdgeInd_I>=12)
        FatalErrorInFunction<<"Edge index out of range:"<<mc33LocalEdgeInd_I<< exit(FatalError);
    label cutPntEdgeOldInd_I = mc33cube.cutEdgeVerticeIndex[mc33LocalEdgeInd_I];
    if(cutPntEdgeOldInd_I==-1 || cutPntEdgeOldInd_I>=oldToNewPointInd.size())
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    label cutPntEdgeInd_I = oldToNewPointInd[cutPntEdgeOldInd_I];
    if(cutPntEdgeInd_I==-1 || verticeOfCell.find(cutPntEdgeInd_I)==verticeOfCell.end() || cutPntEdgeInd_I>=new_points.size())
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    return cutPntEdgeInd_I;
}

label Foam::cutCellFvMesh::getFaceIndFromPntList
(
    const cell& thisCell,
    std::vector<label> pointList
)
{
    label matchingFaceInd = -1;
    for(label locFaceInd=0; locFaceInd<thisCell.size(); locFaceInd++)
    {
        bool allInside = true;
        label faceInd = thisCell[locFaceInd];
        const face& thisFace = new_faces[faceInd];
        for(label pntInd : pointList)
        {
            if(thisFace.which(pntInd)==-1)
                allInside=false;
        }
        Info<<thisFace<<" "<<allInside<<Foam::endl;
        if(allInside)
        {
            if(matchingFaceInd!=-1)
                FatalErrorInFunction<<"Points match to multiple faces!"<< exit(FatalError);
            matchingFaceInd = thisCell[faceInd];
        }
    }
    if(matchingFaceInd==-1)
    {
        Info<<"pointList: "<<Foam::endl;
        for(label pnt : pointList)
        {
            Info<<" "<<pnt;
        }
        Info<<Foam::endl;
        
        Info<<"thisCellFaces:"<<Foam::endl;
        for(label faceInd : thisCell)
        {
            Info<<this->new_faces[faceInd]<<Foam::endl;
        }
        Info<<Foam::endl;
        FatalErrorInFunction<<"Points match to no face!"<< exit(FatalError);
    }
    return matchingFaceInd;
}

bool Foam::cutCellFvMesh::verticesAreIdentical
(
    const corrFaceVertice& vert1,
    const corrFaceVertice& vert2,
    const MC33::MC33Cube& mc33cube,
    corrFaceVertice& vertResult
)
{
    corrFaceVertice match;
    
    label globalPnt1;
    label globalPnt2;
    if(vert1.type==VType::cubePnt)
    {
        globalPnt1 = mc33cube.vertices[vert1.value];
        if(vert2.type==VType::cubePnt)
        {
            if(vert1.value==vert2.value)
            {
                vertResult.type = vert1.type;
                vertResult.value = vert1.value;
                return true;
            }
            else
                return false;
        }
        else if(vert2.type==VType::cutEdge)
        {
            globalPnt2 = mc33cube.cutEdgeVerticeIndex[vert2.value];
            if(globalPnt1==-1 || globalPnt2==-1)
                FatalErrorInFunction<<"Error!"<< exit(FatalError);
            if(globalPnt1==globalPnt2)
            {
                vertResult.type = vert1.type;
                vertResult.value = vert1.value;
                return true;
            }
            else
                return false;
        }
        else
            FatalErrorInFunction<<"Not made for this type of vertice"<< exit(FatalError);
    }
    else if(vert1.type==VType::cutEdge)
    {
        globalPnt1 = mc33cube.cutEdgeVerticeIndex[vert1.value];
        if(vert2.type==VType::cubePnt)
        {
            globalPnt2 = mc33cube.vertices[vert2.value];
            if(globalPnt1==-1 || globalPnt2==-1)
                FatalErrorInFunction<<"Error!"<< exit(FatalError);
            if(globalPnt1==globalPnt2)
            {
                vertResult.type = vert2.type;
                vertResult.value = vert2.value;
                return true;
            }
            else
                return false;
        }
        else if(vert2.type==VType::cutEdge)
        {
            globalPnt2 = mc33cube.cutEdgeVerticeIndex[vert2.value];
            if(globalPnt1==-1 || globalPnt2==-1)
                FatalErrorInFunction<<"Error!"<< exit(FatalError);
            if(globalPnt1==globalPnt2)
            {
                vertResult.type = vert2.type;
                vertResult.value = vert2.value;
                return true;
            }
            else
                return false;
        }
        else
            FatalErrorInFunction<<"Not made for this type of vertice"<< exit(FatalError);
    }
    else
        FatalErrorInFunction<<"Not made for this type of vertice"<< exit(FatalError);
    return false;
}

void Foam::cutCellFvMesh::getGlobalPoint
(
    const corrFaceVertice& vert,
    const MC33::MC33Cube& mc33cube,
    const pointField& points,
    label& verticeInd,
    vector& verticePoint
)
{   
    label globalPnt=-1;
    if(vert.type==VType::cubePnt)
        globalPnt = mc33cube.vertices[vert.value];
    else if(vert.type==VType::cutEdge)
        globalPnt = mc33cube.cutEdgeVerticeIndex[vert.value];
    else if(vert.type==VType::centerPnt)
        globalPnt = mc33cube.centerPointInd;
    else
        FatalErrorInFunction<<"Not made for this type of vertice"<< exit(FatalError);    
    if(globalPnt==-1)
        FatalErrorInFunction<<"Point does not exist in mc33Cube"<< exit(FatalError);
    
    globalPnt = oldToNewPointInd[globalPnt];
    if(globalPnt==-1)
        FatalErrorInFunction<<"Point does not exist as new point"<< exit(FatalError);
    
    verticeInd = globalPnt;
    verticePoint = points[globalPnt];
}

label Foam::cutCellFvMesh::getGlobalPoint
(
    const corrFaceVertice& vert,
    const MC33::MC33Cube& mc33cube,
    const std::unordered_map<std::uint8_t,label>& edgeToAddedPntInd,
    const pointField& points
)
{
    label globalPnt=-1;
    if(vert.type!=VType::edgePntAdded)
    {
        vector globalPntVec;
        getGlobalPoint(vert,mc33cube,points,globalPnt,globalPntVec);
    }
    else
    {
        auto iter = edgeToAddedPntInd.find(vert.value);
        if(iter==edgeToAddedPntInd.end())
            FatalErrorInFunction<<"Point does not exist"<< exit(FatalError);
        globalPnt = iter->second;
    }
    if(globalPnt==-1)
        FatalErrorInFunction<<"Invalid point label!"<< exit(FatalError);
    return globalPnt;    
}

face Foam::cutCellFvMesh::abstrFaceToGlobalFaceTransfer
(
    const corrFace& abstrFace,
    const MC33::MC33Cube& mc33cube,
    const std::unordered_map<std::uint8_t,label>& edgeToAddedPntInd,
    const faceList& faces,
    const pointField& points
)
{
    DynamicList<label> globalVertices;
    switch (abstrFace.type)
    {
        case FType::mc33Triangle:
        {
            if(abstrFace.faceData.size()!=3)
                FatalErrorInFunction<<"MC33 triangle must be sized 3!"<< exit(FatalError);
            if(abstrFace.origFace!=-2)
                FatalErrorInFunction<<"MC33 triangle have a original face of -2!"<< exit(FatalError);
            for(const corrFaceVertice& vert : abstrFace.faceData)
            {
                globalVertices.append(getGlobalPoint(vert,mc33cube,edgeToAddedPntInd,points));
            }
            break;
        }
        case FType::old:
        {
            if(abstrFace.origFace<0 || abstrFace.origFace>=nbrOfPrevFaces)
                FatalErrorInFunction<<"MC33 triangle have a valid original face number!"<< exit(FatalError);
            const DynamicList<label>& cutFaceInds = oldFacesToCutFaces_[abstrFace.origFace];
            if(cutFaceInds.size()<2)
                FatalErrorInFunction<<"Old face must have more cutFaces!"<< exit(FatalError);
            DynamicList<label> posCutFaceInds;
            for(label cutFaceInd : cutFaceInds)
            {
                label signFace = cutFacesToSide_[cutFaceInd];
                if(signFace>=0)
                    posCutFaceInds.append(signFace);
            }
            DynamicList<face> posCutFace;
            for(label posCutFaceInd : posCutFaceInds)
            {
                const face& thisFace = cutFaces_[posCutFaceInd];
                DynamicList<label> thisFaceNewPntsList;
                for(label oldVerticeInd : thisFace)
                {
                    label newVerticeInd = oldToNewPointInd[oldVerticeInd];
                    if(newVerticeInd==-1)
                        FatalErrorInFunction<<"Invalid new point indices!"<< exit(FatalError);
                    thisFaceNewPntsList.append(newVerticeInd);
                }
                posCutFace.append(face(thisFaceNewPntsList));
            }
            if(posCutFace.size()==0)
                FatalErrorInFunction<<"No remaining faces!"<< exit(FatalError);

            if(abstrFace.faceData.size()==0)
            {
                if(posCutFace.size()!=1)
                    FatalErrorInFunction<<"More than one possible face but no hint!"<< exit(FatalError);
                for(label vertices : posCutFace[0])
                     globalVertices.append(vertices);
            }
            else
            {
                List<std::unordered_set<label>> faceVertices(posCutFace.size());
                for(label faceInd=0; faceInd<posCutFace.size(); faceInd++)
                {
                    for(label vert : posCutFace[faceInd])
                    {
                        faceVertices[faceInd].insert(vert);
                    }
                }

                label matchInd = -1;
                for(label faceInd=0; faceInd<posCutFace.size(); faceInd++)
                {
                    bool allMatch = true;
                    for(const corrFaceVertice& vert : abstrFace.faceData)
                    {
                        auto iter = faceVertices[faceInd].find(getGlobalPoint(vert,mc33cube,edgeToAddedPntInd,points));
                        if(iter==faceVertices[faceInd].end())
                            allMatch=false;
                    }
                    if(allMatch)
                    {
                        if(matchInd!=-1)
                            FatalErrorInFunction<<"Hint matches with two faces!"<< exit(FatalError);
                        matchInd=faceInd;
                    }
                }
                if(matchInd==-1)
                    FatalErrorInFunction<<"Hint matches with no face!"<< exit(FatalError);
                for(label vertices : posCutFace[matchInd])
                     globalVertices.append(vertices);
            }
            break;
        }
        case FType::orig:
        {
            if(abstrFace.origFace<0 || abstrFace.origFace>=nbrOfPrevFaces)
                FatalErrorInFunction<<"MC33 triangle have a valid original face number!"<< exit(FatalError);
            const DynamicList<label>& cutFaceInds = oldFacesToCutFaces_[abstrFace.origFace];
            if(cutFaceInds.size()!=0)
                FatalErrorInFunction<<"Old face must not have cutFaces!"<< exit(FatalError);
            label origFaceInd = abstrFace.origFace;
            origFaceInd = reverseFaceMap[origFaceInd];
            if(origFaceInd<0 || origFaceInd>=faces.size())
                FatalErrorInFunction<<"New face index not valid!"<< exit(FatalError);
            globalVertices.append(faces[origFaceInd]);
            
            std::unordered_set<label> globalVerticesSet;
            for(label vert : globalVertices)
            {
                globalVerticesSet.insert(vert);
            }
            for(const corrFaceVertice& vert : abstrFace.faceData)
            {
                auto iter = globalVerticesSet.find(getGlobalPoint(vert,mc33cube,edgeToAddedPntInd,points));
                if(iter==globalVerticesSet.end())
                    FatalErrorInFunction<<"Hint does not match!"<< exit(FatalError);                
            }
            break;
        }
        case FType::splitOld:
        {
            if(abstrFace.origFace<0 || abstrFace.origFace>=nbrOfPrevFaces)
                FatalErrorInFunction<<"MC33 triangle have a valid original face number!"<< exit(FatalError);
            const DynamicList<label>& cutFaceInds = oldFacesToCutFaces_[abstrFace.origFace];
            if(cutFaceInds.size()<2)
                FatalErrorInFunction<<"Old face must have more cutFaces!"<< exit(FatalError);
            DynamicList<label> posCutFaceInds;
            for(label cutFaceInd : cutFaceInds)
            {
                label signFace = cutFacesToSide_[cutFaceInd];
                if(signFace>=0)
                    posCutFaceInds.append(signFace);
            }
            DynamicList<face> posCutFace;
            for(label posCutFaceInd : posCutFaceInds)
            {
                const face& thisFace = cutFaces_[posCutFaceInd];
                DynamicList<label> thisFaceNewPntsList;
                for(label oldVerticeInd : thisFace)
                {
                    label newVerticeInd = oldToNewPointInd[oldVerticeInd];
                    if(newVerticeInd==-1)
                        FatalErrorInFunction<<"Invalid new point indices!"<< exit(FatalError);
                    thisFaceNewPntsList.append(newVerticeInd);
                }
                posCutFace.append(face(thisFaceNewPntsList));
            }
            if(posCutFace.size()==0)
                FatalErrorInFunction<<"No remaining faces!"<< exit(FatalError);
            if(abstrFace.faceData.size()<3)
                FatalErrorInFunction<<"Split old face must be larger than 3!"<< exit(FatalError);
            for(const corrFaceVertice& vert : abstrFace.faceData)
            {
                globalVertices.append(getGlobalPoint(vert,mc33cube,edgeToAddedPntInd,points));
            }
            
            std::unordered_set<label> globalVerticesSet;
            for(label vert : globalVertices)
            {
                globalVerticesSet.insert(vert);
            }
            label allInOneFace = false;
            for(const face& oneFace : posCutFace)
            {
                bool allInThis = true;
                for(label vert : oneFace)
                {
                    auto iter = globalVerticesSet.find(vert);
                    if(iter==globalVerticesSet.end())
                        allInThis = false;
                }
                if(allInThis)
                {
                    if(allInOneFace)
                        FatalErrorInFunction<<"Double matching!"<< exit(FatalError);
                    else
                        allInOneFace = true;
                }
            }
            if(!allInOneFace)
                FatalErrorInFunction<<"No matching!"<< exit(FatalError);
            break;
        }
        case FType::splitOrig:
        {
            if(abstrFace.origFace<0 || abstrFace.origFace>=nbrOfPrevFaces)
                FatalErrorInFunction<<"MC33 triangle have a valid original face number!"<< exit(FatalError);
            const DynamicList<label>& cutFaceInds = oldFacesToCutFaces_[abstrFace.origFace];
            if(cutFaceInds.size()!=0)
                FatalErrorInFunction<<"Old face must not have cutFaces!"<< exit(FatalError);
            label origFaceInd = abstrFace.origFace;
            origFaceInd = reverseFaceMap[origFaceInd];
            if(origFaceInd<0 || origFaceInd>=faces.size())
                FatalErrorInFunction<<"New face index not valid!"<< exit(FatalError);
            if(abstrFace.faceData.size()<3)
                FatalErrorInFunction<<"Split old face must be larger than 3!"<< exit(FatalError);
            
            const face& origFace = faces[origFaceInd];
            for(const corrFaceVertice& vert : abstrFace.faceData)
            {
                globalVertices.append(getGlobalPoint(vert,mc33cube,edgeToAddedPntInd,points));
            }
            
            std::unordered_set<label> globalVerticesSet;
            for(label vert : globalVertices)
            {
                globalVerticesSet.insert(vert);
            }
            for(label vert : origFace)
            {
                auto iter = globalVerticesSet.find(vert);
                if(iter==globalVerticesSet.end())
                    FatalErrorInFunction<<"No matching!"<< exit(FatalError);
            }
            break;
        }
        case FType::added:
        {
            if(abstrFace.faceData.size()!=3)
                FatalErrorInFunction<<"Added face must be sized 3!"<< exit(FatalError);
            if(abstrFace.origFace!=-1)
                FatalErrorInFunction<<"Added face must have a original face of -1!"<< exit(FatalError);
            for(const corrFaceVertice& vert : abstrFace.faceData)
            {
                globalVertices.append(getGlobalPoint(vert,mc33cube,edgeToAddedPntInd,points));
            }
            break;
        }
    }

    return face(globalVertices);
}

Foam::cutCellFvMesh::FaceRelation Foam::cutCellFvMesh::facesOverlay
(
    const face& face1,
    const face& face2
)
{
    std::unordered_set<label> face1Set(face1.begin(),face1.end());
    std::unordered_set<label> face2Set(face2.begin(),face2.end());
    
    if(face1.size()==face2.size())
    {
        bool allMatch = true;
        label matching=0;
        for(label vert1 : face1)
        {
            auto iter =  face2Set.find(vert1);
            if(iter==face2Set.end())
                allMatch = false;
            else
                matching++;
        }
        if(allMatch)
            return FaceRelation::Ident;
        else
        {
            if(matching>2)
                return FaceRelation::mixedOverlay;
            else
                return FaceRelation::none;
        }
    }
    else if(face1.size()>face2.size())
    {
        bool allMatch = true;
        label matching=0;
        for(label vert2 : face2)
        {
            auto iter =  face1Set.find(vert2);
            if(iter==face1Set.end())
                allMatch = false;
            else
                matching++;
        }
        if(allMatch)
            return FaceRelation::twoInOne;
        else
        {
            if(matching>2)
                return FaceRelation::mixedOverlay;
            else
                return FaceRelation::none;
        }
    }
    else if(face2.size()>face1.size())
    {
        bool allMatch = true;
        label matching=0;
        for(label vert1 : face1)
        {
            auto iter =  face2Set.find(vert1);
            if(iter==face2Set.end())
                allMatch = false;
            else
                matching++;
        }
        if(allMatch)
            return FaceRelation::oneInTwo;
        else
        {
            if(matching>2)
                return FaceRelation::mixedOverlay;
            else
                return FaceRelation::none;
        }
    }
    return FaceRelation::none;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::generateCellSplitData
(
    const cell& thisCell,
    const label cellInd,
    corrConcaveCellData& corrData,
    std::function<bool(const corrConcaveCellData&, const List<face>&, label)> isCellValid
)
{
    label oldCellIndex = cellMap[cellInd];
    const MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
        
    auto newCellDataPtr = std::unique_ptr<CellSplitData>(new CellSplitData());
    CellSplitData& newCellData = *newCellDataPtr;
    newCellData.originalCell = cellInd;

    //Process possible added point
    bool addedPntsInConvexCorr = false;
    std::unordered_map<std::uint8_t,label> edgeToAddedPntInd;
    for(corrFace& oneFace : corrData.faces)
    {
        for(corrFaceVertice& oneVert : oneFace.faceData)
        {
            if(oneVert.type==VType::edgePntAdded)
            {
                addedPntsInConvexCorr=true;
                
                const addedPointPos& position = oneVert.position;
                const std::pair<VType,std::uint8_t> distStart = position.distStart;
                const std::pair<VType,std::uint8_t> distEnd = position.distEnd;
                const std::pair<VType,std::uint8_t> measureStart = position.measureStart;
                const std::pair<VType,std::uint8_t> measureEnd = position.measureEnd;
                corrFaceVertice distStartVert = {distStart.first,distStart.second};
                corrFaceVertice distEndVert = {distEnd.first,distEnd.second};
                corrFaceVertice measureStartVert = {measureStart.first,measureStart.second};
                corrFaceVertice measureEndVert = {measureEnd.first,measureEnd.second};
                corrFaceVertice match;
                if(verticesAreIdentical(distStartVert,distEndVert,mc33cube,match))
                {
                    if(match.type!=VType::cubePnt)
                        FatalErrorInFunction<<"Error"<< exit(FatalError);
                    if(edgeToAddedPntInd.find(oneVert.value)!=edgeToAddedPntInd.end())
                    {
                        edgeToAddedPntInd[oneVert.value] = measureStartVert.value;
                    }
                    oneVert = {position.measureStart.first,position.measureStart.second};
                }
                else
                {
                    label distStartInd;
                    vector distStartVec;
                    getGlobalPoint(oneVert,mc33cube,new_points,distStartInd,distStartVec);
                    label distEndInd;
                    vector distEndVec;
                    getGlobalPoint(oneVert,mc33cube,new_points,distEndInd,distEndVec);
                    label measureStartInd;
                    vector measureStartVec;
                    getGlobalPoint(oneVert,mc33cube,new_points,measureStartInd,measureStartVec);
                    label measureEndInd;
                    vector measureEndVec;
                    getGlobalPoint(oneVert,mc33cube,new_points,measureEndInd,measureEndVec);
                    
                    scalar distLen = norm2(distStartVec-distEndVec);
                    scalar measureLen = norm2(measureStartVec-measureEndVec);
                    
                    vector addedPoint = (distLen/measureLen) * (measureEndVec-measureStartVec) + measureStartVec;
                    
                    if(edgeToAddedPntInd.find(oneVert.value)!=edgeToAddedPntInd.end())
                    {
                        edgeToAddedPntInd[oneVert.value] = newCellData.addedPoints.size();
                        newCellData.addedPoints.append({addedPoint,edge(measureStartInd,measureEndInd)});
                    }
                }
            }
        }
    }
    
    //Reduce duplicate vertices in faces
    for(corrFace& oneFace : corrData.faces)
    {
        DynamicList<corrFaceVertice> faceData;
        label vertInd0 = 0;
        corrFaceVertice& vert0 = oneFace.faceData[vertInd0];
        faceData.append(vert0);
        do
        {
            label vertInd1 = (vertInd0+1)%oneFace.faceData.size();
            corrFaceVertice& vert1 = oneFace.faceData[vertInd1];
            corrFaceVertice match;
            if(!verticesAreIdentical(vert0,vert1,mc33cube,match))
            {
                faceData.append(vert1);
            }
            vert0 = vert1;
            vertInd0++;
        }
        while(vertInd0<oneFace.faceData.size());
        oneFace.faceData = faceData;
    }
    
    //Remove empty faces from cells
    for(List<label>& oneCell : corrData.cells)
    {
        DynamicList<label> reducedCell;
        for(label oneFaceInd : oneCell)
        {
            const corrFace& oneFace = corrData.faces[oneFaceInd];
            if(oneFace.faceData.size()>2)
                reducedCell.append(oneFaceInd);
        }
        oneCell = reducedCell;
    }
    
    //Create cell in newCellData
    DynamicList<face> nonConvexCellFaces;
    nonConvexCellFaces.setSize(corrData.faces.size())
    for(label faceInd=0; faceInd<corrData.faces.size(); faceInd++)
    {
        face& globFace = nonConvexCellFaces[faceInd];
        const corrFace& abstrFace = corrData.faces[faceInd];
        if(abstrFace.faceData.size()>2)
        {
            globFace = abstrFaceToGlobalFaceTransfer(abstrFace,mc33cube,edgeToAddedPntInd,new_faces,new_points);
        }
    }

    //Compute equal faces and faces as face subsets
    List<std::unordered_set<label>> equalFaces(nonConvexCellFaces.size());
    List<std::unordered_set<label>> faceSubsets(nonConvexCellFaces.size());
    List<std::unordered_set<label>> faceSupersets(nonConvexCellFaces.size());
    for(uint baseInd=0; baseInd<nonConvexCellFaces.size(); baseInd++)
    {
        face baseFace = nonConvexCellFaces[baseInd];
        std::unordered_set<label> baseFaceSet(baseFace.begin(),baseFace.end());
        for(uint compareInd=0; compareInd<nonConvexCellFaces.size(); compareInd++)
        {
            if(baseInd==compareInd)
                continue;

            face compFace = nonConvexCellFaces[compareInd];
            std::unordered_set<label> compFaceSet(compFace.begin(),compFace.end());
            bool allIn=true;
            for(label compVert : compFace)
            {
                if(baseFaceSet.find(compVert)==baseFaceSet.end())
                    allIn=false;
            }
            if(allIn)
            {
                if(baseFace.size()==compFace.size())
                {
                    equalFaces[baseInd].insert(compareInd);
                }
                else if(baseFace.size()>compFace.size())
                {
                    faceSubsets[baseInd].insert(compareInd);
                    faceSupersets[compareInd].insert(baseInd);
                }
                else
                    FatalErrorInFunction<<"Error"<< exit(FatalError);
            }
        }
    }
    for(uint faceInd=0; faceInd<nonConvexCellFaces.size(); faceInd++)
    {
        for(auto iter=equalFaces[faceInd].begin(); iter!=equalFaces[faceInd].end(); iter++)
        {
            label otherFace = *iter;
            if(equalFaces[otherFace].find(faceInd)==equalFaces[otherFace].end())
                FatalErrorInFunction<<"Inconsistent"<< exit(FatalError);
        }
        label nbrCategory=0;
        if(equalFaces[faceInd].size()>0)
            nbrCategory++;
        if(faceSubsets[faceInd].size()>0)
            nbrCategory++;
        if(faceSupersets[faceInd].size()>0)
            nbrCategory++;
        if(nbrCategory>1)
            FatalErrorInFunction<<"Face must be equal, a superset or a subset but not both!"<< exit(FatalError);
    }
    
    //Compute splitted faces for subsets
    List<bool> faceSubsetIsFilling(nonConvexCellFaces.size(),false);
    List<DynamicList<face>> faceSubsetsRemain(nonConvexCellFaces.size());
    for(uint faceInd=0; faceInd<faceSubsets.size(); faceInd++)
    {
        const face& baseFace = nonConvexCellFaces[faceInd];
        const std::unordered_set<label>& subFaces = faceSubsets[faceInd];
        DynamicList<face> remainingFaces;
        remainingFaces.append(baseFace);
        
        for(auto faceIndIter=subFaces.begin(); faceIndIter!=subFaces.end();)
        {
            const face& oneSubFace = nonConvexCellFaces[*faceIndIter];
            
            //Get remaining face to insert subface
            label remFaceInd = -1;
            for(label i=0; i<remainingFaces.size(); i++)
            {
                const face& oneRemainingFace = remainingFaces[i];
                //Test if all subface vertices are in one remaining faces
                bool allIn=true;
                for(label oneSubFaceInd : oneSubFace)
                    if(oneRemainingFace.which(oneSubFaceInd)==-1)
                        allIn=false;
                if(allIn)
                {
                    if(remFaceInd!=-1)
                        FatalErrorInFunction<<"Subface can not be in two remaining faces"<< exit(FatalError);
                    remFaceInd = i;
                }
            }
            if(remFaceInd==-1)
                FatalErrorInFunction<<"Subface must be in at least one remaining face"<< exit(FatalError);
            
            const face& remFace = remainingFaces[remFaceInd];
            if(oneSubFace.size()==remFace.size())
            {
                DynamicList<face> newRemainingFaces;
                for(label j=0; j<remainingFaces.size(); j++)
                    if(j!=remFaceInd)
                        newRemainingFaces.append(remainingFaces[j]);
                remainingFaces = newRemainingFaces;
                continue;
            }
            if(oneSubFace.size()>remFace.size())
                FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);

            //Search for subFace edges cutting the remFace
            DynamicList<label> oneSubFaceCuti0s;
            for(label i0=0; i0<oneSubFace.size(); i0++)
            {
                label i1 = oneSubFace.fcIndex(i0);
                label remainingFaceLocIndi0 = remFace.which(oneSubFace[i0]);
                if(!(remFace.nextLabel(remainingFaceLocIndi0)==oneSubFace[i1] ||
                     remFace.prevLabel(remainingFaceLocIndi0)==oneSubFace[i1]))
                {
                    oneSubFaceCuti0s.append(i0);
                }
            }
            if(oneSubFaceCuti0s.size()<1)
                FatalErrorInFunction<<"Subface must have at least one cutting edge"<< exit(FatalError);
            
            //Split remaining face
            DynamicList<face> newRemainingFaces;
            for(label j=0; j<oneSubFaceCuti0s.size(); j++)
            {
                label globIndSubFacei0 = oneSubFace[oneSubFaceCuti0s[j]];
                label globIndSubFacei1 = oneSubFace.nextLabel(oneSubFaceCuti0s[j]);
                
                label locIndRemFacei0 = remFace.which(globIndSubFacei0);                    
                label locIndRemFacei1 = remFace.which(globIndSubFacei1);
                if(locIndRemFacei0==-1 || locIndRemFacei1==-1)
                    FatalErrorInFunction<<"Cut edge not in remFace!"<< exit(FatalError);

                bool rotPValid = true;
                bool closedP = false;
                DynamicList<label> rotPInds;
                label locIndRemFacei0rotP = locIndRemFacei0;
                rotPInds.append(remFace[locIndRemFacei0rotP]);
                for(label k=0; k<remFace.size(); k++)
                {
                    locIndRemFacei0rotP = remFace.rcIndex(locIndRemFacei0rotP);
                    label globIndRemFace = remFace[locIndRemFacei0rotP];
                    rotPInds.append(globIndRemFace);
                    if(globIndSubFacei1==globIndRemFace)
                    {
                        closedP = true;
                        break;
                    }
                    else if(oneSubFace.which(globIndRemFace)!=-1)
                    {
                        rotPValid = false;
                    }
                }
                if(!closedP)
                    FatalErrorInFunction<<"Non closed cut face!"<< exit(FatalError);

                bool rotNValid = true;
                bool closedN = false;
                DynamicList<label> rotNInds;
                label locIndRemFacei0rotN = locIndRemFacei0;
                rotNInds.append(remFace[locIndRemFacei0rotN]);                
                for(label k=0; k<remFace.size(); k++)
                {
                    locIndRemFacei0rotN = remFace.fcIndex(locIndRemFacei0rotN);
                    label globIndRemFace = remFace[locIndRemFacei0rotN];
                    rotNInds.append(globIndRemFace);
                    if(globIndSubFacei1==globIndRemFace)
                    {
                        closedN = true;
                        break;
                    }
                    else if(oneSubFace.which(globIndRemFace)!=-1)
                    {
                        rotNValid = false;
                    }
                }
                if(!closedN)
                    FatalErrorInFunction<<"Non closed cut face!"<< exit(FatalError);
                
                if(rotPValid && !rotNValid)
                {
                    newRemainingFaces.append(face(rotNInds));
                }
                else if(!rotPValid && rotNValid)
                {
                    newRemainingFaces.append(face(rotPInds));
                }
                else
                    FatalErrorInFunction<<"Only one face can be valid"<< exit(FatalError);
            }
            
            //
            remainingFaces[remFaceInd] = newRemainingFaces[0];
            for(label j=1; j<newRemainingFaces.size(); j++)
            {
                remainingFaces.append(newRemainingFaces[j]);
            }
        }
        
        std::unordered_set<label> subFaceSet;
        for(auto faceIndIter=subFaces.begin(); faceIndIter!=subFaces.end();)
        {
            const face& oneSubFace = nonConvexCellFaces[*faceIndIter];
            for(label vert : oneSubFace)
                subFaceSet.insert(vert);
        }
        bool mustBeNonFilling = false;
        for(label vert : baseFace)
        {
            if(subFaceSet.find(vert)==subFaceSet.end())
                mustBeNonFilling = true;
        }
        
        if(remainingFaces.size()==0)
        {
            faceSubsetIsFilling[faceInd] = true;
            if(mustBeNonFilling)
                FatalErrorInFunction<<"Filling mismatch!"<< exit(FatalError);
        }
        else
        {
            if(!mustBeNonFilling)
                FatalErrorInFunction<<"Filling mismatch!"<< exit(FatalError);
            faceSubsetsRemain[faceInd].append(remainingFaces);
        }
    }

    // Redirect faces to their replacing faces
    List<bool> validFaces(nonConvexCellFaces.size(),true);
    List<std::pair<bool,DynamicList<label>>> faceReplacedByFaces(nonConvexCellFaces.size(),{false,DynamicList<label>()});
    for(uint faceInd=0; faceInd<validFaces.size(); faceInd++)
    {
        bool invalidBcEqual = false;
        if(equalFaces[faceInd].size()>0)
        {
            std::vector<label> equalFacesToThis(equalFaces[faceInd].begin(),equalFaces[faceInd].end());
            std::sort(equalFacesToThis.begin(),equalFacesToThis.end());
            if(equalFacesToThis[0]<faceInd)
            {
                validFaces[faceInd] = false;
                faceReplacedByFaces[faceInd].first = true;
                faceReplacedByFaces[faceInd].second.append(equalFacesToThis[0]);
            }
        }
        if(faceSubsets[faceInd].size()>0)
        {
            if(!validFaces[faceInd])
                FatalErrorInFunction<<"Face can not both be equal and have a subset!"<< exit(FatalError);
            validFaces[faceInd] = false;
            faceReplacedByFaces[faceInd].first = true;
            for(auto iter=faceSubsets[faceInd].begin(); iter!=faceSubsets[faceInd].end(); iter++)
            {
                if(faceSubsets[*iter].size()>0)
                    FatalErrorInFunction<<"Recursive face subset!"<< exit(FatalError);
                faceReplacedByFaces[faceInd].second.append(*iter);
            }
        }
    }
    
    // Reduce valid cells
    DynamicList<List<label>> trueCells;
    for(label cellInd=0; cellInd<corrData.cells.size(); cellInd++)
    {
        const List<label>& oneCellFaceInds = corrData.cells[cellInd];
        DynamicList<label> oneCellValidFaces;
        for(label faceInd : oneCellFaceInds)
            if(validFaces[faceInd])
                oneCellValidFaces.append(faceInd);
        
        //Compute neighborhood set for valid faces
        List<std::unordered_set<label>> neighborFaceOnEdge(oneCellValidFaces.size());
        List<bool> removed(oneCellValidFaces.size(),false);
        for(label iCell=0; iCell<oneCellValidFaces.size(); iCell++)
        // for all faces in cell
        {
            label oneFaceInd = oneCellValidFaces[iCell];
            const face& oneFace = nonConvexCellFaces[oneFaceInd];
            for(label i0Face=0; i0Face<oneFace.size(); i0Face++)
            // for all edges
            {
                label i1Face = oneFace.fcIndex(i0Face);
                label i0Vert = oneFace[i0Face];
                label i1Vert = oneFace[i1Face];
                bool edgeConnected=false;
                
                for(label iCellNei=0; iCellNei<oneCellValidFaces.size(); iCellNei++)
                {
                    if(iCell==iCellNei)
                        continue;
                    label oneFaceIndNei = oneCellValidFaces[iCellNei];
                    const face& oneFaceNei = nonConvexCellFaces[oneFaceIndNei];
                    
                    label neiFacei0 = oneFaceNei.which(i0Vert);
                    if(neiFacei0!=-1)
                    {
                        if(oneFaceNei.prevLabel(neiFacei0)==i1Vert ||
                           oneFaceNei.nextLabel(neiFacei0)==i1Vert)
                        {
                            if(neighborFaceOnEdge[iCell].find(iCellNei)!=neighborFaceOnEdge[iCell].end())
                                FatalErrorInFunction<<"Face borders same face twice!"<< exit(FatalError);
                            neighborFaceOnEdge[iCell].insert(iCellNei);
                        }
                    }
                }
            }
            if(neighborFaceOnEdge[iCell].size()<oneFace.size())
            {
                removed[iCell] = true;
            }
        }
        for(label iCell=0; iCell<oneCellValidFaces.size(); iCell++)
        // for all faces in cell
        {
            std::unordered_set<label>& neighborFaceOnEdgeSingle = neighborFaceOnEdge[iCell];
            for(label iCell1=0; iCell1<oneCellValidFaces.size(); iCell1++)
            {
                if(removed[iCell1])
                {
                    auto iter = neighborFaceOnEdgeSingle.find(iCell1);
                    neighborFaceOnEdgeSingle.erase(iter);
                }
            }
        }
        for(label iCell=0; iCell<oneCellValidFaces.size(); iCell++)
        // for all faces in cell
        {
            label oneFaceInd = oneCellValidFaces[iCell];
            const face& oneFace = nonConvexCellFaces[oneFaceInd];
            if(!removed[iCell])
            {
                if(neighborFaceOnEdge[iCell].size()<oneFace.size())
                {
                    FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
                }
            }
        }
        
        //List valid faces
        std::unordered_set<label> validFaceLocInd;
        for(label i=0;i<removed.size();i++)
            if(!removed[i])
                validFaceLocInd.insert(i);
        
        //Collect valid cells from valid faces
        List<std::unordered_set<label>> neighborFaceOnEdgeCopy = neighborFaceOnEdge;
        List<bool> treatedFace(neighborFaceOnEdgeCopy.size(),false);
        DynamicList<label> edgeFaceFront;
        DynamicList<label> trueCell;
        
        while(!validFaceLocInd.empty())
        // over all cells
        {
            if(edgeFaceFront.size()==0)
            {
                if(trueCell.size()>0)
                {
                    trueCells.append(trueCell);
                    trueCell.clear();
                }
                auto iter = validFaceLocInd.begin();
                edgeFaceFront.append(*iter);
                treatedFace[*iter] = true;
                validFaceLocInd.erase(iter);
            }
            
            // advance front
            DynamicList<label> new_edgeFaceFront;
            for(label locInd : edgeFaceFront)
            {
                trueCell.append(locInd);
                std::unordered_set<label>& neighborFaceOnEdgeSingle = neighborFaceOnEdgeCopy[locInd];
                for(auto iterNei=neighborFaceOnEdgeSingle.begin(); iterNei!=neighborFaceOnEdgeSingle.end(); iterNei++)
                {
                    if(!treatedFace[*iterNei])
                    {
                        new_edgeFaceFront.append(*iterNei);
                        treatedFace[*iterNei] = true;
                        auto iterValidFace = validFaceLocInd.find(*iterNei);
                        if(iterValidFace == validFaceLocInd.end())
                            FatalErrorInFunction<<"Must not happen!"<< exit(FatalError);
                        validFaceLocInd.erase(iterValidFace);
                    }
                }
            }
            edgeFaceFront = new_edgeFaceFront;
        }
    }
    
    //Create cells from valid faces alone
    List<DynamicList<label>> trueCellsTrueFaces(trueCells.size());
    for(label cellInd=0; cellInd<trueCells.size(); cellInd++)
    {
        const List<label>& cell = trueCells[cellInd];
        for(label faceInd : cell)
        {
            if(faceReplacedByFaces[faceInd].first)
            {
                trueCellsTrueFaces[cellInd].append(faceReplacedByFaces[faceInd].second);
            }
            else
            {
                trueCellsTrueFaces[cellInd].append(faceInd);
            }
        }
    }
    
    //Mark every face that is used for replacement of another face
    List<bool> faceUsedForReplacement(nonConvexCellFaces.size(),false);
    for(std::pair<bool,DynamicList<label>> faceReplacement : faceReplacedByFaces)
    {
        for(label replacementFace : faceReplacement.second)
        {
            if(faceUsedForReplacement[replacementFace])
                FatalErrorInFunction<<"Face used for replacement twice!"<< exit(FatalError);
            faceUsedForReplacement[replacementFace] = true;
        }
    }
    
    auto faceToVerticeCombNumber = [](const face& oneFace)
    {
        std::vector<label> verticeInds;
        for(label verticeInd : oneFace)
            verticeInds.push_back(verticeInd);
        std::sort(verticeInds.begin(),verticeInds.end());
        uint verticeCombNumber = 0;
        for(uint dec=0; dec<verticeInds.size()-1; dec++)
        {
            verticeCombNumber += verticeInds[dec];
            verticeCombNumber *= 10;
        }
        verticeCombNumber += verticeInds[verticeInds.size()-1];
        return verticeCombNumber;
    }
    
    std::unordered_map<uint,label> faceVerticesToInd;
    for(label faceInds : thisCell)
    {
        const face& oneFace = new_faces[faceInds];
        uint verticeCombNumber = faceToVerticeCombNumber(oneFace);        
        auto iter = faceVerticesToInd.find(verticeCombNumber);
        if(iter!=faceVerticesToInd.end())
            FatalErrorInFunction<<"Duplicate face!"<< exit(FatalError);
        faceVerticesToInd[verticeCombNumber] = faceInds;
    }
    
    // Compute new face reference index
    std::unordered_set<label> cellsFacesSet;
    cellsFacesSet.insert(thisCell.begin(),thisCell.end());
    List<label> newFaceInd(nonConvexCellFaces.size(),-1);
    for(label faceInd=0; faceInd<corrData.faces.size(); faceInd++)
    {
        const corrFace& oneFace = corrData.faces[faceInd];
        const face& globFace = nonConvexCellFaces[faceInd];
        label matchingNewFaceInd = -1;
        if(oneFace.type!=FType::mc33Triangle && oneFace.type!=FType::added)
        {
            label origFaceCellLocalInd = corrData.faces.origFace;
            label origFaceGlobalInd = mc33cube.origFaces[origFaceCellLocalInd];
            auto iterOldToNew = oldFaceToNewFaceInds.find(origFaceGlobalInd);
            if(iterOldToNew == oldFaceToNewFaceInds.end())
                FatalErrorInFunction<<"Old face not existing!"<< exit(FatalError);
            const DynamicList<label>& possNewFaceInds = iterOldToNew->second;
            for(label newFaceInd : possNewFaceInds)
            {
                const face& newFace = new_faces[newFaceInd];
                bool identicalFace = true;
                for(label vertice : globFace)
                {
                    if(newFace.which(vertice)==-1)
                    {
                        identicalFace = false;
                    }
                }
                if(identicalFace)
                {
                    if(matchingNewFaceInd!=-1)
                        FatalErrorInFunction<<"Identical to more than one face!"<< exit(FatalError);
                    matchingNewFaceInd = newFaceInd;
                }
            }
            if(matchingNewFaceInd==-1)
                FatalErrorInFunction<<"No identical face!"<< exit(FatalError);
            if(cellsFacesSet.find(matchingNewFaceInd)==cellsFacesSet.end())
                FatalErrorInFunction<<"Matching face not in cell!"<< exit(FatalError);
            newFaceInd[faceInd] = matchingNewFaceInd;
        }
        else if(oneFace.type==FType::mc33Triangle)
        {
            for(label newFaceInd : thisCell)
            {
                const face& newFace = new_faces[newFaceInd];
                bool identicalFace = true;
                for(label vertice : globFace)
                {
                    if(newFace.which(vertice)==-1)
                    {
                        identicalFace = false;
                    }
                }
                if(identicalFace)
                {
                    if(matchingNewFaceInd!=-1)
                        FatalErrorInFunction<<"Identical to more than one face!"<< exit(FatalError);
                    matchingNewFaceInd = newFaceInd;
                }
            }
            if(matchingNewFaceInd==-1)
                FatalErrorInFunction<<"No identical face!"<< exit(FatalError);
        }
    }
    
    List<std::pair<CSFaceType,label>> facesStoreInfo(nonConvexCellFaces.size());
    for(label faceInd=0; faceInd<corrData.faces.size(); faceInd++)
    {
        const corrFace& oneFace = corrData.faces[faceInd];
        const face& globFace = nonConvexCellFaces[faceInd];
        switch(oneFace.type)
        {
            case FType::old:
            case FType::orig:
            case FType::mc33Triangle:
                if(newFaceInd[faceInd]==-1)
                    FatalErrorInFunction<<"No new face identical!"<< exit(FatalError);
                facesStoreInfo[faceInd] = {CSFaceType::original,newFaceInd[faceInd]};                    
                break;
            case FType::added:
                facesStoreInfo[faceInd] = {CSFaceType::added,newCellData.addedFace.size()};
                newCellData.addedFace.append(oneFace);
                break;
            case FType::splitOld:
            case FType::splitOrig:
                if(newFaceInd[faceInd]==-1)
                    FatalErrorInFunction<<"No new face identical!"<< exit(FatalError);
                facesStoreInfo[faceInd] = {CSFaceType::splitted,newCellData.splittedFaces.size()};
                newCellData.splittedFaces.append({oneFace,newFaceInd[faceInd]});
                break;
            default:
                FatalErrorInFunction<<"Invalid face type!"<< exit(FatalError);
        }
    }
    
    enum FaceUsage {None,Partial,Complete};
    List<std::pair<FaceUsage,label>> faceUsedInSubCell(nonConvexCellFaces.size(),FaceUsage::None);
    for(label cellInd=0; cellInd<trueCellsTrueFaces.size(); cellInd++)
    {
        const List<label>& oneCell = trueCellsTrueFaces[cellInd];
        CellFaces oneCellData;
        
        DynamicList<label> originalFaceInds;
        DynamicList<label> splittedFaceInds;
        DynamicList<label> addedFaceInds;
        DynamicList<std::tuple<CSFaceType,label,CSFaceType,label>> replacingFaceInds;
        DynamicList<std::tuple<label,CSFaceType,label> replacingResidualFaceInds;
        
        for(label cellFaceInd : oneCell)
        {
            bool thisFaceValid = validFaces[cellFaceInd];
            const std::pair<bool,DynamicList<label>>& faceReplacement = faceReplacedByFaces[cellFaceInd];
            bool thisFaceIsReplaced = faceReplacement.first;
            bool thisFaceUsedForReplacement = faceUsedForReplacement[cellFaceInd];
            const DynamicList<face>& remainingFaces = faceSubsetsRemain[cellFaceInd];
            bool replacedFaceWithRemains = (remainingFaces.size()>0);
            
            if(thisFaceValid && faceReplacement.second.size()!=0)
                FatalErrorInFunction<<"Face valid but is replaced!"<< exit(FatalError);
            if(!thisFaceValid && faceReplacement.second.size()==0)
                FatalErrorInFunction<<"Face invalid but no replacement!"<< exit(FatalError);
            if(thisFaceUsedForReplacement && !thisFaceValid)
                FatalErrorInFunction<<"Face used for replacement but invalid!"<< exit(FatalError);
            if(!thisFaceIsReplaced && faceReplacement.second.size()>0)
                FatalErrorInFunction<<"Inconsistent face replacement data!"<< exit(FatalError);
            if(thisFaceIsReplaced && thisFaceUsedForReplacement)
                FatalErrorInFunction<<"Face replaced but also used for replacement"<< exit(FatalError);
            if(replacedFaceWithRemains && !thisFaceIsReplaced)
                FatalErrorInFunction<<"Face has remaining faces but not replaced"<< exit(FatalError);
            
            const corrFace& oneFace = corrData.faces[faceInd];
            const face& globFace = nonConvexCellFaces[faceInd];
            
            const std::pair<CSFaceType,label>& faceInfo = facesStoreInfo[cellFaceInd];
            label newFaceIndex = newFaceInd[cellFaceInd];
            
            if(thisFaceIsReplaced)
            {
                const DynamicList<label>& replacingFacesInd = faceReplacement.second;
                switch(oneFace.type)
                {
                    case FType::orig:
                    case FType::old:
                        bool replacedByMc33OrAdded = false;
                        for(label replacingFaceInd : replacingFacesInd)
                        {
                            const corrFace& replacingFace = corrData.faces[replacingFaceInd];
                            switch(replacingFace.type)
                            {
                                case FType::orig:
                                case FType::old:
                                case FType::splitOrig:
                                case FType::splitOld:
                                    FatalErrorInFunction<<"Invalid replacment of standard face by standard face!"<< exit(FatalError);
                                    break;
                                case FType::added:
                                case FType::mc33:
                                    replacedByMc33OrAdded = true;
                                    faceUsedInSubCell[replacingFaceInd] = Complete;
                                    break;
                                default:
                                    FatalErrorInFunction<<"Invalid face type!"<< exit(FatalError);
                            }
                        }
                        if(!replacedByMc33OrAdded);
                            FatalErrorInFunction<<"Original face can only be replaced by mc33 or added when part of a cell"<< exit(FatalError);
                        if(replacedFaceWithRemains && !replacedByMc33OrAdded)
                            FatalErrorInFunction<<"Invalid replacment of standard face by standard face!"<< exit(FatalError);
                        for(const face& remFace : remainingFaces)
                        {
                            splittedFaceInds.append(newCellData.splittedFaces.size());
                            newCellData.splittedFaces.append({remFace,faceInfo.second});
                        }
                        if(replacedFaceWithRemains)
                        {
                            faceUsedInSubCell[cellFaceInd].first = Partial;
                            faceUsedInSubCell[cellFaceInd].second++;
                        }
                        else
                            faceUsedInSubCell[cellFaceInd].first = Complete;
                        break;
                    case FType::splitOrig:
                    case FType::splitOld:
                        bool replacedByMc33OrAdded = false;
                        for(label replacingFaceInd : replacingFacesInd)
                        {
                            const corrFace& replacingFace = corrData.faces[replacingFaceInd];
                            switch(replacingFace.type)
                            {
                                case FType::orig:
                                case FType::old:
                                case FType::splitOrig:
                                case FType::splitOld:
                                    FatalErrorInFunction<<"Invalid replacment of standard face by standard face!"<< exit(FatalError);
                                case FType::added:
                                case FType::mc33:
                                    replacedByMc33OrAdded = true;
                                    faceUsedInSubCell[replacingFaceInd] = Complete;
                                default:
                                    FatalErrorInFunction<<"Invalid face type!"<< exit(FatalError);
                            }
                        }
                        if(!replacedByMc33OrAdded);
                            FatalErrorInFunction<<"Original face can only be replaced by mc33 or added when part of a cell"<< exit(FatalError);
                        if(replacedFaceWithRemains && !replacedByMc33OrAdded)
                            FatalErrorInFunction<<"Invalid replacment of standard face by standard face!"<< exit(FatalError);
                        for(const face& remFace : remainingFaces)
                        {
                            splittedFaceInds.append(newCellData.splittedFaces.size());
                            newCellData.splittedFaces.append({remFace,faceInfo.second});
                        }
                        if(replacedFaceWithRemains)
                        {
                            faceUsedInSubCell[cellFaceInd].first = Partial;
                            faceUsedInSubCell[cellFaceInd].second++;
                        }
                        else
                            faceUsedInSubCell[cellFaceInd].first = Complete;
                        break;
                    case FType::added:
                        for(label replacingFaceInd : replacingFacesInd)
                        {
                            const corrFace& replacingFace = corrData.faces[replacingFaceInd];
                            const std::pair<CSFaceType,label>& replacingFaceInfo = facesStoreInfo[replacingFaceInd];
                            switch(replacingFace.type)
                            {
                                case FType::orig:
                                case FType::old:
                                    replacingFaceInds.append({CSFaceType::added,faceInfo.second,
                                                              CSFaceType::original,replacingFaceInfo.second});
                                    faceUsedInCell[replacingFaceInd] = true;
                                    faceUsedInSubCell[replacingFaceInd] = Complete;
                                    break;
                                case FType::splitOrig:
                                case FType::splitOld:
                                    replacingFaceInds.append({CSFaceType::added,faceInfo.second,
                                                              CSFaceType::splitted,replacingFaceInfo.second});
                                    faceUsedInCell[replacingFaceInd] = true;
                                    faceUsedInSubCell[replacingFaceInd] = Complete;
                                    break;
                                case FType::added:
                                    FatalErrorInFunction<<"Invalid replacment of standard face by standard face!"<< exit(FatalError);
                                case FType::mc33:
                                    replacingFaceInds.append({CSFaceType::added,faceInfo.second,
                                                              CSFaceType::original,replacingFaceInfo.second});
                                    faceUsedInCell[replacingFaceInd] = true;
                                    faceUsedInSubCell[replacingFaceInd] = Complete;
                                    break;
                                default:
                                    FatalErrorInFunction<<"Invalid face type!"<< exit(FatalError);
                            }
                        }
                        for(const face& remFace : remainingFaces)
                        {
                            addedFaceInds.append(newCellData.addedFace.size());
                            newCellData.addedFaceCells[newCellData.addedFace.size()].append(cellInd);
                            newCellData.addedFace.append(remFace);
                        }
                        if(replacedFaceWithRemains)
                        {
                            faceUsedInSubCell[cellFaceInd].first = Partial;
                            faceUsedInSubCell[cellFaceInd].second++;
                        }
                        else
                            faceUsedInSubCell[cellFaceInd].first = Complete;
                        break;
                    case FType::mc33:
                        for(label replacingFaceInd : replacingFacesInd)
                        {
                            const corrFace& replacingFace = corrData.faces[replacingFaceInd];
                            switch(replacingFace.type)
                            {
                                case FType::orig:
                                case FType::old:
                                case FType::splitOrig:
                                case FType::splitOld:
                                    FatalErrorInFunction<<"Invalid replacment of standard face by standard face!"<< exit(FatalError);
                                case FType::added:
                                case FType::mc33:
                                    FatalErrorInFunction<<"Invalid replacment of standard face by added or mc33 face! No cell!"<< exit(FatalError);
                                default:
                                    FatalErrorInFunction<<"Invalid face type!"<< exit(FatalError);
                            }
                        }
                        if(remainingFaces.size>0)
                            FatalErrorInFunction<<"Invalid replacment of mc33 face!"<< exit(FatalError);
                        break;
                    default:
                        FatalErrorInFunction<<"Invalid face type!"<< exit(FatalError);
                }
            }
            else if(thisFaceUsedForReplacement)
            {
            }
            else
            {
                switch(faceInfo.first)
                {
                    case CSFaceType::original:                
                        originalFaceInds.append(faceInfo.second);
                        break;
                    case CSFaceType::added:
                        newCellData.addedFaceCells[faceInfo.second].append(cellInd);
                        addedFaceInds.append(faceInfo.second);
                        break;
                    case CSFaceType::splitted:
                        splittedFaceInds.append(faceInfo.second);
                        break;
                    default:
                        FatalErrorInFunction<<"Invalid face type!"<< exit(FatalError);
                }
                faceUsedInCell[cellFaceInd] = true;
            }
        }
        oneCellData.originalFaceInds = originalFaceInds;
        oneCellData.splittedFaceInds = splittedFaceInds;
        oneCellData.addedFaceInds = addedFaceInds;
        oneCellData.addedFaceCells = addedFaceCells;

        newCellData.cells.append(oneCellData);
    }
    
    
    for(label faceInd=0; faceInd<faceUsedInCell.size(); faceInd++)
    {
        bool usedForReplacement = faceUsedForReplacement[faceInd];
        std::pair<bool,DynamicList<label>>& beingReplaced = faceUsedForReplacement[faceInd];
        
        const DynamicList<face>& remainingFaces = faceSubsetsRemain[cellFaceInd];
        bool replacedFaceWithRemains = (remainingFaces.size()>0);

        if(faceUsedInSubCell[faceInd].first == None)
        {
            if(usedForReplacement && beingReplaced.first)
                FatalErrorInFunction<<"Circular replacement not possible!"<< exit(FatalError);
            if(beingReplaced.first)
            {
                const corrFace& oneFace = corrData.faces[faceInd];
                const face& globFace = nonConvexCellFaces[faceInd];
                
                const DynamicList<label>& replacingFacesInd = beingReplaced.second;
                switch(oneFace.type)
                {
                    case FType::orig:
                    case FType::old:
                        bool replacedByMc33OrAdded = false;
                        for(label replacingFaceInd : replacingFacesInd)
                        {
                            if(faceUsedInSubCell[replacingFaceInd].first != None)
                                FatalErrorInFunction<<"Not used replaced face is replaced by used face!"<< exit(FatalError);
                            const corrFace& replacingFace = corrData.faces[replacingFaceInd];
                            switch(replacingFace.type)
                            {
                                case FType::orig:
                                case FType::old:
                                case FType::splitOrig:
                                case FType::splitOld:
                                    FatalErrorInFunction<<"Invalid replacment of standard face by standard face!"<< exit(FatalError);
                                    break;
                                case FType::added:
                                    FatalErrorInFunction<<"Invalid replacment of orig face by added face!"<< exit(FatalError);
                                    break;
                                case FType::mc33:
                                    nonCellFaceChanges.append({CSFaceType::original,newFaceInd[faceInd],CSFaceType::original,newFaceInd[replacingFaceInd]});
                                    faceUsedInSubCell[replacingFaceInd] = Complete;
                                    break;
                                default:
                                    FatalErrorInFunction<<"Invalid face type!"<< exit(FatalError);
                            }
                        }
                        if(replacedFaceWithRemains)
                        {
                            FatalErrorInFunction<<"Off cell replaced face can not be replaced with remains!"<< exit(FatalError);
                        }
                        break;
                    case FType::splitOrig:
                    case FType::splitOld:
                        bool replacedByMc33OrAdded = false;
                        for(label replacingFaceInd : replacingFacesInd)
                        {
                            if(faceUsedInSubCell[replacingFaceInd].first != None)
                                FatalErrorInFunction<<"Not used replaced face is replaced by used face!"<< exit(FatalError);
                            const corrFace& replacingFace = corrData.faces[replacingFaceInd];
                            switch(replacingFace.type)
                            {
                                case FType::orig:
                                case FType::old:
                                case FType::splitOrig:
                                case FType::splitOld:
                                    FatalErrorInFunction<<"Invalid replacment of standard face by standard face!"<< exit(FatalError);
                                    break;
                                case FType::added:
                                    FatalErrorInFunction<<"Invalid replacment of orig face by added face!"<< exit(FatalError);
                                    break;
                                case FType::mc33:
                                    nonCellFaceChanges.append({CSFaceType::original,newFaceInd[faceInd],CSFaceType::original,newFaceInd[replacingFaceInd]});
                                    faceUsedInSubCell[replacingFaceInd] = Complete;
                                    break;
                                default:
                                    FatalErrorInFunction<<"Invalid face type!"<< exit(FatalError);
                            }
                        }
                        if(replacedFaceWithRemains)
                        {
                            FatalErrorInFunction<<"Off cell replaced face can not be replaced with remains!"<< exit(FatalError);
                        }
                        break;
                    case FType::added:
                        for(label replacingFaceInd : replacingFacesInd)
                        {
                            const corrFace& replacingFace = corrData.faces[replacingFaceInd];
                            const std::pair<CSFaceType,label>& replacingFaceInfo = facesStoreInfo[replacingFaceInd];
                            FatalErrorInFunction<<"Added face can not be noncell replaced!"<< exit(FatalError);
                            switch(replacingFace.type)
                            {
                                case FType::orig:
                                case FType::old:
                                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                                    break;
                                case FType::splitOrig:
                                case FType::splitOld:
                                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                                    break;
                                case FType::added:
                                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                                case FType::mc33:
                                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                                    break;
                                default:
                                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                            }
                        }
                        break;
                    case FType::mc33:
                        for(label replacingFaceInd : replacingFacesInd)
                        {
                            const corrFace& replacingFace = corrData.faces[replacingFaceInd];
                            switch(replacingFace.type)
                            {
                                case FType::orig:
                                case FType::old:
                                case FType::splitOrig:
                                case FType::splitOld:
                                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                                case FType::added:
                                case FType::mc33:
                                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                                default:
                                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                            }
                        }
                        if(remainingFaces.size>0)
                            FatalErrorInFunction<<"Invalid replacment of mc33 face!"<< exit(FatalError);
                        break;
                    default:
                        FatalErrorInFunction<<"Invalid face type!"<< exit(FatalError);
                }
            }
        }
    }
    
    return newCellDataPtr;
}

void Foam::cutCellFvMesh::testCellSplitData
(
    const cell& thisCell,
    const faceList& faces,
    const pointField& points,
    const CellSplitData& data
)
{
    std::unordered_set<label> facesSet(thisCell.begin(),thisCell.end());
    std::unordered_set<label> pointsSet;
    for(label faceInd : thisCell)
    {
        for(label pointInd : faces[faceInd])
        {
            pointsSet.insert(pointInd);
        }
    }
    
    if(data.addedPoints.size()>0)
    {
        if(data.addedPointsLimit!=nbrOfPrevPoints)
            FatalErrorInFunction<<"Mismatch in CellSplit Data structure!"<< exit(FatalError);
    }    
    std::unordered_set<label> pointSetAdded;
    for(label i=0;i<data.addedPoints.size();i++)
    {
        pointSetAdded.insert(i+data.addedPointsLimit);
    }
    
    for(const CellFaces& oneCell : data.cells)
    {
        for(label originalFace : oneCell.originalFaceInds)
        {
            if(facesSet.find(originalFace)!=facesSet.end())
                FatalErrorInFunction<<"Original face is not in cell!"<< exit(FatalError);
        }
        for(label splittedFaceInd : oneCell.splittedFaceInds)
        {
            const std::pair<face,label>& oneSplittedFace = data.splittedFaces[splittedFaceInd];
            if(facesSet.find(oneSplittedFace.second)!=facesSet.end())
                FatalErrorInFunction<<"Original face of splitted face is not in cell!"<< exit(FatalError);
            for(label pointInd : oneSplittedFace.first)
                if(pointsSet.find(pointInd)!=pointsSet.end() &&
                   pointSetAdded.find(pointInd)!=pointSetAdded.end())
                    FatalErrorInFunction<<"Split face point is not in cell!"<< exit(FatalError);
        }
        for(label addedFaceInd : oneCell.addedFaceInds)
        {
            const face& addedFace = data.addedFace[addedFaceInd];
            for(label pointInd : addedFace)
            {
                if(pointsSet.find(pointInd)!=pointsSet.end() &&
                   pointSetAdded.find(pointInd)!=pointSetAdded.end())
                    FatalErrorInFunction<<"Added point is not in cell!"<< exit(FatalError);
            }
        }        
    }
    
    for(const std::pair<vector,edge>& addedPoint : data.addedPoints)
    {
        label edgeStart = addedPoint.second.start();
        if(pointsSet.find(edgeStart)!=pointsSet.end())
            FatalErrorInFunction<<"Edge point is not in cell!"<< exit(FatalError);
        label edgeEnd = addedPoint.second.end();
        if(pointsSet.find(edgeEnd)!=pointsSet.end())
            FatalErrorInFunction<<"Edge point is not in cell!"<< exit(FatalError);
        vector addedPointVec = addedPoint.first;

        vector edgeVec = points[edgeEnd] - points[edgeStart];
        scalar edgeLen = norm2(edgeVec);
        
        vector partEdgeStart = addedPointVec - points[edgeStart];
        scalar partEdgeStartLen = norm2(partEdgeStart);
    
        vector partEdgeEnd = addedPointVec - points[edgeEnd];
        scalar partEdgeEndLen = norm2(partEdgeEnd);
        
        if(std::abs(edgeLen - (partEdgeStartLen+partEdgeEndLen))>=1e-8)
            FatalErrorInFunction<<"Added point not in edge!"<< exit(FatalError);
    }
}

bool Foam::cutCellFvMesh::redMarkSide(const MC33::MC33Cube& mc33cube)
{
    const std::array<unsigned short int,8>& permut = permutationTable[mc33cube.permutationTableIndex];
    const std::array<std::int8_t,8>& refSign = refSignVertice[mc33cube.cubeCase];
    std::array<std::int8_t,8> thisCellPattern;
    for(uint i=0;i<mc33cube.vertices.size();i++)
    {
        unsigned short int permInd = permut[i];
        label realVertIndPrevMesh = mc33cube.vertices[permInd];
        thisCellPattern[i] = (pointsToSide_[realVertIndPrevMesh]>=0)?+1:-1;
    }
    if(refSign == thisCellPattern)
        return true;
    for(std::int8_t& sign : thisCellPattern)
        sign = sign*-1;
    if(refSign == thisCellPattern)
        return false;
    FatalErrorInFunction<<"Non matching pattern!"<< exit(FatalError);
    return true;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC2
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C2"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    const MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c2)
        FatalErrorInFunction<<"Error"<< exit(FatalError);

    const markSide& c2convexCorr = convexCorrectionData[MC33::Case::c2];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
        corrData = c2convexCorr.redMarkPlus;
    else
        corrData = c2convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC32
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{   
    Info<<"Nonconvex correction C32"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c32)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c32convexCorr = convexCorrectionData[MC33::Case::c32];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
    {}
    else
        corrData = c32convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC412
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C412"<<Foam::endl;

    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c412)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c412convexCorr = convexCorrectionData[MC33::Case::c412];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
    {}
    else
        corrData = c412convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC5
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C5"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c5)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c5convexCorr = convexCorrectionData[MC33::Case::c5];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
        corrData = c5convexCorr.redMarkPlus;
    else
        corrData = c5convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC612
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C612"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c612)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c612convexCorr = convexCorrectionData[MC33::Case::c612];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
    {}
    else
        corrData = c612convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC62
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C62"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c62)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c62convexCorr = convexCorrectionData[MC33::Case::c62];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
    {}
    else
        corrData = c62convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC72
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C72"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c72_720 && 
       mc33cube.cubeCase!=MC33::Case::c72_760 && 
       mc33cube.cubeCase!=MC33::Case::c72_800)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c72convexCorr = convexCorrectionData[mc33cube.cubeCase];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
    {}
    else
        corrData = c72convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC73
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C73"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c73_840 && 
       mc33cube.cubeCase!=MC33::Case::c73_912 && 
       mc33cube.cubeCase!=MC33::Case::c73_984)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c73convexCorr = convexCorrectionData[mc33cube.cubeCase];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
        corrData = c73convexCorr.redMarkPlus;
    else
        corrData = c73convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC8
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C8"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c8)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c8convexCorr = convexCorrectionData[mc33cube.cubeCase];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
        corrData = c8convexCorr.redMarkPlus;
    else
        corrData = c8convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC1011
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C1011"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c1011_1190 && 
       mc33cube.cubeCase!=MC33::Case::c1011_1202)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    auto newCellDataPtr = std::unique_ptr<CellSplitData>(new CellSplitData());
    CellSplitData& newCellData = *newCellDataPtr;
    newCellData.originalCell = cellInd;

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC1211
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C1211"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c1211_1358 && 
       mc33cube.cubeCase!=MC33::Case::c1211_1406)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c1211convexCorr = convexCorrectionData[mc33cube.cubeCase];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
        corrData = c1211convexCorr.redMarkPlus;
    else
        corrData = c1211convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC122
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{   
    Info<<"Nonconvex correction C122"<<Foam::endl;
    
    FatalErrorInFunction<<"Not yet implemented!"<< exit(FatalError);

    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c122_1646 && 
       mc33cube.cubeCase!=MC33::Case::c122_1742)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c122convexCorr = convexCorrectionData[mc33cube.cubeCase];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
        corrData = c122convexCorr.redMarkPlus;
    else
        corrData = c122convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC132
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C132"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c132_MostNeg && 
       mc33cube.cubeCase!=MC33::Case::c132_MostPos)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c132convexCorr = convexCorrectionData[mc33cube.cubeCase];
    corrConcaveCellData corrData;
    if(mc33cube.cubeCase==MC33::Case::c132_MostNeg)
    {
        if(redMarkSide(mc33cube))
        {}
        else
            corrData = c132convexCorr.redMarkNotPlus;
    }
    else
    {
        if(redMarkSide(mc33cube))
            corrData = c132convexCorr.redMarkPlus;
        else
        {}
    }
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC133
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C133"<<Foam::endl;

    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c133_MostNeg && 
       mc33cube.cubeCase!=MC33::Case::c133_MostPos)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c133convexCorr = convexCorrectionData[mc33cube.cubeCase];
    corrConcaveCellData corrData;
    if(mc33cube.cubeCase==MC33::Case::c133_MostNeg)
    {
        if(redMarkSide(mc33cube))
            corrData = c133convexCorr.redMarkPlus;
        else
            corrData = c133convexCorr.redMarkNotPlus;
    }
    else
    {
        if(redMarkSide(mc33cube))
            corrData = c133convexCorr.redMarkPlus;
        else
            corrData = c133convexCorr.redMarkNotPlus;
    }
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC134
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C134"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c134)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c134convexCorr = convexCorrectionData[mc33cube.cubeCase];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
        corrData = c134convexCorr.redMarkPlus;
    else
        corrData = c134convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<Foam::cutCellFvMesh::CellSplitData> Foam::cutCellFvMesh::nonConvexCellSplitC14
(
    const cell& thisCell,
    const label cellInd,
    std::unique_ptr<List<scalar>> faceToCellCenterRelation
)
{
    Info<<"Nonconvex correction C14"<<Foam::endl;
    
    label oldCellIndex = cellMap[cellInd];
    MC33::MC33Cube& mc33cube = mc33CutCellData[oldCellIndex];
    if(mc33cube.cell!=oldCellIndex)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    if(mc33cube.cubeCase!=MC33::Case::c14)
        FatalErrorInFunction<<"Error"<< exit(FatalError);
    
    const markSide& c14convexCorr = convexCorrectionData[mc33cube.cubeCase];
    corrConcaveCellData corrData;
    if(redMarkSide(mc33cube))
        corrData = c14convexCorr.redMarkPlus;
    else
        corrData = c14convexCorr.redMarkNotPlus;
    
    std::unique_ptr<CellSplitData> newCellDataPtr = generateCellSplitData(thisCell,cellInd,corrData);

    return newCellDataPtr;
}

std::unique_ptr<cellList> Foam::cutCellFvMesh::createCells
(
    const faceList& faces,
    const labelList& owner,
    const labelList& neighbour 
)
{
    std::unordered_map<label,DynamicList<label>> cellsFacesMap;
    label cellNbr = -1;
    for(int faceInd=0; faceInd<faces.size(); faceInd++)
    {
        label thisOwner = owner[faceInd];
        label thisNeighbour = neighbour[faceInd];
        cellsFacesMap[thisOwner].append(faceInd);
        if(thisNeighbour!=-1)
            cellsFacesMap[thisNeighbour].append(faceInd);
        cellNbr = std::max(thisOwner,cellNbr);
        cellNbr = std::max(thisNeighbour,cellNbr);
    }
    cellNbr++;
    auto cells = std::unique_ptr<cellList>(new cellList(cellNbr));
    for(int cellInd=0; cellInd<cells->size(); cellInd++)
    {
        auto iter = cellsFacesMap.find(cellInd);
        if(iter==cellsFacesMap.end())
            FatalErrorInFunction<<"Missing cell index!"<< exit(FatalError);
        if(iter->second.size()<4)
            FatalErrorInFunction<<"Illformed cell has less than four faces!"<< exit(FatalError);
        
        (*cells)[cellInd] = cell(iter->second);
    }
    return cells;
}

void Foam::cutCellFvMesh::testCellClosure
(
    const faceList& faces,
    const cellList& cells
)
{
    for(const cell& oneCell : cells)
    {
        DynamicList<face> cellFaces;
        for(const label& faceInd : oneCell)
        {
            const face& oneFace = faces[faceInd];
            //Info<<"----------------"<<Foam::endl;
            //Info<<"oneFace:"<<oneFace<<Foam::endl;
            for(const label& vertice : faces[faceInd])
            {
                label locVertInd = oneFace.which(vertice);
                //Info<<" locVertInd:"<<locVertInd<<Foam::endl;
                if(locVertInd==-1)
                    FatalErrorInFunction<<"Error"<< exit(FatalError);
                label prevVertice = oneFace.prevLabel(locVertInd);
                //Info<<" prevVertice:"<<prevVertice<<Foam::endl;
                bool prevVerticeCon = false;
                label nextVertice = oneFace.nextLabel(locVertInd);
                //Info<<" nextVertice:"<<nextVertice<<Foam::endl;
                bool nextVerticeCon = false;
                DynamicList<label> cellConFaceInds;
                for(const label& othFaceInd : oneCell)
                {
                    if(othFaceInd == faceInd)
                        continue;
                    
                    const face& othFace = faces[othFaceInd];
                    if(othFace.which(vertice)!=-1)
                    {
                        cellConFaceInds.append(othFaceInd);
                    }
                }
                for(const label& conFaceInd : cellConFaceInds)
                {
                    const face& conFace = faces[conFaceInd];
                    //Info<<" conFace:"<<conFace<<Foam::endl;
                    label conLocVertInd = conFace.which(vertice);
                    if(conLocVertInd==-1)
                    {
                        Info<<" vertice:"<<vertice<<Foam::endl;
                        FatalErrorInFunction<<"Error"<< exit(FatalError);
                    }
                    label conPrevVertice = conFace.prevLabel(conLocVertInd);
                    label conNextVertice = conFace.nextLabel(conLocVertInd);
                    
                    if(prevVertice == conPrevVertice)
                    {
                        if(prevVerticeCon)
                            FatalErrorInFunction<<"Error"<< exit(FatalError);
                        prevVerticeCon = true;
                    }
                    else if(prevVertice == conNextVertice)
                    {
                        if(prevVerticeCon)
                            FatalErrorInFunction<<"Error"<< exit(FatalError);
                        prevVerticeCon = true;
                    }
                    
                    if(nextVertice == conPrevVertice)
                    {
                        if(nextVerticeCon)
                            FatalErrorInFunction<<"Error"<< exit(FatalError);
                        nextVerticeCon = true;
                    }
                    else if(nextVertice == conNextVertice)
                    {
                        if(nextVerticeCon)
                            FatalErrorInFunction<<"Error"<< exit(FatalError);
                        nextVerticeCon = true;
                    }
                }
                
                if(!prevVerticeCon || !nextVerticeCon)
                {
                    FatalErrorInFunction<<"Edges connected to vertice are not connected to other face!"<< exit(FatalError);
                }
            }
        }
    }
}

void Foam::cutCellFvMesh::correctFaceNormal
(
    faceList& faces,
    const labelList& owner,
    const labelList& neighbour,
    const cellList& cells
)
{
    for(int faceInd=0; faceInd<faces.size(); faceInd++)
    {
        label ownerCellInd = owner[faceInd];
        cell ownerCell = cells[ownerCellInd];
        //Info<<"ownerCell:"<<ownerCell<<Foam::endl;
        //Info<<"this->cells()[this->owner()[faceInd]]:"<<this->cells()[this->owner()[faceInd]]<<Foam::endl;
        //Info<<"this->cells()[this->owner()[faceInd]]:"<<this->cells()[this->owner()[faceInd]].points(this->faces(),this->points())<<Foam::endl;
        std::unordered_set<label> ownerCellFaceMap;
        for(const label& faceInd : ownerCell)
            ownerCellFaceMap.insert(faceInd);
        if(ownerCellFaceMap.find(faceInd)==ownerCellFaceMap.end())
            FatalErrorInFunction<<"Error"<< exit(FatalError);
        ownerCellFaceMap.erase(faceInd);
        
        vector centre = faces[faceInd].centre(new_points);
        vector normal = faces[faceInd].normal(new_points);        
        label posNormalHit = 0;
        DynamicList<vector> posNormalHitPoint;
        label negNormalHit = 0;
        DynamicList<vector> negNormalHitPoint;
        //Info<<"-----------------"<<Foam::endl;
        //Info<<"faceInd:"<<faceInd<<" -"<<faces[faceInd].points(new_points)<<Foam::endl;
        //Info<<"centre:"<<centre<<Foam::endl;
        //Info<<"normal:"<<normal<<Foam::endl;
        
        for(auto iter=ownerCellFaceMap.cbegin();
            iter!=ownerCellFaceMap.cend();
            iter++)
        {
            label othFaceInd = *iter;
            //Info<<" othFaceInd:"<<othFaceInd<<" "<<faces[othFaceInd].points(new_points)<<Foam::endl;
            const face& othFace = faces[othFaceInd];
            
            pointHit posNormRes= othFace.ray(centre,normal,new_points);
            if(posNormRes.hit() && posNormRes.distance()>0)
            {
                posNormalHit++;
                posNormalHitPoint.append(posNormRes.hitPoint());
                //Info<<"     Pos Hits "<<othFaceInd<<" using "<<normal<<" at dist:"<<posNormRes.distance()<<" and point:"<<posNormRes.hitPoint()<<Foam::endl;
            }
            
            pointHit negNormRes= othFace.ray(centre,-1*normal,new_points);
            if(negNormRes.hit() && negNormRes.distance()>0)
            {
                negNormalHit++;
                negNormalHitPoint.append(negNormRes.hitPoint());
                //Info<<"     Neg Hits "<<othFaceInd<<" using "<<-1*normal<<" at dist:"<<negNormRes.distance()<<" and point:"<<negNormRes.hitPoint()<<Foam::endl;
            }
        }
        
        //Info<<"posNormalHit:"<<posNormalHit<<Foam::endl;
        //Info<<"negNormalHit:"<<negNormalHit<<Foam::endl;
        
        //Info<<"posNormalHitPoint:"<<posNormalHitPoint<<Foam::endl;
        //Info<<"negNormalHitPoint:"<<negNormalHitPoint<<Foam::endl;
        
        List<bool> posEqualTreated(posNormalHitPoint.size(),false);
        for(label i=0;i<posNormalHitPoint.size();i++)
        {
            if(posEqualTreated[i])
                continue;
            
            posEqualTreated[i] = true;
            label equalsCount = 0;
            for(int j=i+1;j<posNormalHitPoint.size();j++)
            {
                if(posEqualTreated[j])
                    continue;
                
                scalar posHitPointDist = norm2(posNormalHitPoint[i]-posNormalHitPoint[j]);
                if(posHitPointDist<1e-10)
                {
                    if(posEqualTreated[j])
                        FatalErrorInFunction<<"Error"<< exit(FatalError);
                    equalsCount++;
                    posEqualTreated[j] = true;
                }
            }
            posNormalHit-=equalsCount;
        }
        List<bool> negEqualTreated(negNormalHitPoint.size(),false);
        for(label i=0;i<negNormalHitPoint.size();i++)
        {
            if(negEqualTreated[i])
                continue;
            
            negEqualTreated[i] = true;
            label equalsCount = 0;
            for(int j=i+1;j<negNormalHitPoint.size();j++)
            {
                if(negEqualTreated[j])
                    continue;
                
                scalar negHitPointDist = norm2(negNormalHitPoint[i]-negNormalHitPoint[j]);
                if(negHitPointDist<1e-10)
                {
                    if(negEqualTreated[j])
                    {
                        Info<<"-------------Error------------"<<Foam::endl;
                        Info<<"i:"<<i<<Foam::endl;
                        Info<<"j:"<<j<<Foam::endl;
                        Info<<"negNormalHitPoint:"<<negNormalHitPoint<<Foam::endl;
                        Info<<"negEqualTreated:"<<negEqualTreated<<Foam::endl;
                        FatalErrorInFunction<<"Error"<< exit(FatalError);
                    }
                    equalsCount++;
                    negEqualTreated[j] = true;
                }
            }
            negNormalHit-=equalsCount;
        }
        
        if(posNormalHit%2==0 && negNormalHit%2!=0)
        {
            //Normal direction correct
        }
        else if(posNormalHit%2!=0 && negNormalHit%2==0)
        {
            //new_faces[faceInd] = new_faces[faceInd].reverseFace();
        }
        else
        {
            Info<<"posNormalHit:"<<posNormalHit<<Foam::endl;
            Info<<"negNormalHit:"<<negNormalHit<<Foam::endl;
            Info<<(negNormalHitPoint[0]==negNormalHitPoint[1])<<Foam::endl;
            Info<<(negNormalHitPoint[0]-negNormalHitPoint[1])<<Foam::endl;
            FatalErrorInFunction<<"Error"<< exit(FatalError);
        }
        
    }
}

void Foam::cutCellFvMesh::splitFaceByEdge
(
    const face& oneFace,
    const edge& oneEdge,
    FixedList<face,2>& splitFaces
)
{
    if(oneFace.which(oneEdge.start())==-1 || oneFace.which(oneEdge.end())==-1)
        FatalErrorInFunction<<"Edge not in Face"<< exit(FatalError);   
    label startIndex = oneFace.which(oneEdge.start());
    if(oneFace.nextLabel(startIndex) == oneEdge.end() ||
       oneFace.prevLabel(startIndex) == oneEdge.end())
        FatalErrorInFunction<<"Edge does not split face"<< exit(FatalError);
    
    label index = startIndex;    
    DynamicList<label> face0;
    while(oneFace[index] != oneEdge.end())
    {
        face0.append(oneFace[index]);
        index = oneFace.fcIndex(index);
    }
    face0.append(oneFace[index]);
    
    index = startIndex;
    DynamicList<label> face1;
    while(oneFace[index] != oneEdge.end())
    {
        face1.append(oneFace[index]);
        index = oneFace.rcIndex(index);
    }
    face1.append(oneFace[index]);
    
    splitFaces[0] = face(face0);
    splitFaces[1] = face(face1);
}

bool Foam::cutCellFvMesh::edgeInFace
(
    const face& oneFace,
    const edge& oneEdge
)
{
    return (oneFace.which(oneEdge.start())!=-1 && oneFace.which(oneEdge.end())!=-1);
}

void Foam::cutCellFvMesh::facesEdgesIntersection
(
    const face& totalFace,
    const face& ownFace,
    const face& neiFace,
    List<List<bool>>& faceEdgeIntersections
)
{
    faceEdgeIntersections.setSize(ownFace.size(),List<bool>(neiFace.size(),false));
    for(label o=0; o<ownFace.size(); o++)
    {
        label oPnt0 = ownFace[o];
        label oPnt1 = ownFace[ownFace.fcIndex(o)];
        label totalFaceIndex = totalFace.which(oPnt0);
        for(label n=0; n<neiFace.size(); n++)
        {
            label nPnt0 = neiFace[n];
            label nPnt1 = neiFace[neiFace.fcIndex(n)];
            
            std::unordered_set<label> neiEdgeSet({nPnt0,nPnt1});
            bool openOwn = true;
            bool closedOwn = false;
            bool openNei = false;
            bool closedNei = false;
            
            label runTotalFaceIndex = totalFaceIndex;
            for(label k=0;k<totalFace.size();k++)
            {
                auto iter = neiEdgeSet.find(totalFace[runTotalFaceIndex]);
                if(iter!=neiEdgeSet.end())
                {
                    if(openNei)
                        closedNei=true;
                    else
                        openNei=true;
                    neiEdgeSet.erase(iter);
                }
                if(totalFace[runTotalFaceIndex]==oPnt1)
                    closedOwn=true;
                if(closedOwn)
                    break;
                runTotalFaceIndex = totalFace.fcIndex(runTotalFaceIndex);
            }
            if(openNei && !closedNei)
            {
                faceEdgeIntersections[o][n] = true;
            }
        }
    }
}

bool Foam::cutCellFvMesh::faceEdgesIntersect
(
    const List<List<bool>>& faceEdgeIntersections
)
{
    bool oneIntersection = false;
    for(const List<bool>& innerList : faceEdgeIntersections)
    {
        for(bool val : innerList)
        {
            oneIntersection = oneIntersection || val;
        }
    }
    return oneIntersection;
}

Foam::vector Foam::cutCellFvMesh::computeIntersectionPnt
(
    const edge edgeA,
    const edge edgeB
)
{
    vector AEdgeStartVec = new_points[edgeA.start()];
    vector AEdgeLineVec = edgeA.vec(new_points);

    vector BEdgeStartVec = new_points[edgeB.start()];
    vector BEdgeLineVec = edgeB.vec(new_points);

    vector n = crossProd(AEdgeLineVec,BEdgeLineVec);
    vector a = AEdgeLineVec;
    vector a0 = AEdgeStartVec;
    vector b = BEdgeLineVec;
    vector b0 = BEdgeStartVec;
    vector rhs = a0-b0;
    //Solve Linear Equation System [b,-a,-n](t,sx,sy)^T = a0-b0
    scalar detA = det3x3(b,-a,-n);
    if(std::abs(detA)<1e-10)
        FatalErrorInFunction<<"Edges are parrallel!"<< exit(FatalError);
    scalar t = det3x3(rhs,-a,-n)/detA;
    scalar sx = det3x3(b,rhs,-n)/detA;
    //scalar sy = det3x3(b,-a,rhs)/detA;
    if(t<0 || t>1 || sx<0 || sx>1)
        FatalErrorInFunction<<"Edges have no inner intersection!"<< exit(FatalError);

    vector AEdgeIntPnt = a0 + sx*a;
    vector BEdgeIntPnt = b0 + t*b;
    vector newPoint = (AEdgeIntPnt+BEdgeIntPnt)/2;
    return newPoint;
}

std::unique_ptr<std::vector<std::pair<label,label>>> Foam::cutCellFvMesh::matchingPoints
(
    const face& ownFace,
    const face& neiFace
)
{
    auto result = std::unique_ptr<std::vector<std::pair<label,label>>>(new std::vector<std::pair<label,label>>());
    for(label ownInd=0; ownInd<ownFace.size(); ownInd++)
    {
        label neiInd = neiFace.which(ownFace[ownInd]);
        if(neiInd!=-1)
        {
            result->push_back({ownInd,neiInd});
        }
    }
    return result;
}

std::array<std::pair<label,label>,2> Foam::cutCellFvMesh::matchingEdge
(
    const face& ownFace,
    const face& neiFace
)
{
    std::unique_ptr<std::vector<std::pair<label,label>>> matchPoints = matchingPoints(ownFace,neiFace);
    
    for(unsigned int i=0; i<matchPoints->size(); i++)
    {
        label ni = (i+1)%matchPoints->size();
        label neiFaceLoc_i = (*matchPoints)[i].second;
        label neiFaceLoc_ni = (*matchPoints)[ni].second;
        
        
        if(neiFaceLoc_ni == neiFace.fcIndex(neiFaceLoc_i) ||
           neiFaceLoc_ni == neiFace.rcIndex(neiFaceLoc_i))
        {
            return {(*matchPoints)[i],(*matchPoints)[ni]};
        }
    }
    FatalErrorInFunction<<"Error!"<< exit(FatalError);
    return {std::make_pair<>(0,0),std::make_pair<>(0,0)};
}

Foam::edge Foam::cutCellFvMesh::getEdge
(
    const face& theFace,
    label edgeInd
)
{
    label oPnt0 = theFace[edgeInd];
    label oPnt1 = theFace[theFace.fcIndex(edgeInd)];
    return edge(oPnt0,oPnt1);
}

List<label> Foam::cutCellFvMesh::faceIntersection
(
    const face& totalFace,
    const face& ownFace,
    const face& neiFace,
    edgToIntPntIndMap& edgesToAddPntInd,
    DynamicList<vector>& addPnt
)
{
    // Test for face edges intersection and check for correctness
    List<List<bool>> faceEdgeIntersections;
    facesEdgesIntersection(totalFace,ownFace,neiFace,faceEdgeIntersections);
    for(List<bool> oneEdgeInt : faceEdgeIntersections)
    {
        label numIntersec = 0;
        for(bool intersec : oneEdgeInt)
            if(intersec)
                numIntersec++;
        if(numIntersec>2)
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
    }
    for(label j=0;j<faceEdgeIntersections[0].size();j++)
    {
        label numIntersec = 0;
        for(label i=0;i<faceEdgeIntersections.size();i++)
        {
            if(faceEdgeIntersections[i][j])
                    numIntersec++;
        }
        if(numIntersec>2)
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
    }
    
    // Create new points when necessary
    List<List<label>> faceEdgeIntersectionsAddPntInd(ownFace.size(),List<label>(neiFace.size(),-1));
    label addedPntsIndex = this->new_points.size()+addPnt.size();
    for(label o=0; o<ownFace.size(); o++)
    {
        for(label n=0; n<neiFace.size(); n++)
        {
            if(faceEdgeIntersections[o][n])
            {
                label oPnt0 = ownFace[o];
                label oPnt1 = ownFace[ownFace.fcIndex(o)];
                edge oEdge;
                oEdge[0] = (oPnt0<oPnt1)?oPnt0:oPnt1;
                oEdge[1] = (oPnt0<oPnt1)?oPnt1:oPnt0;
               
                label nPnt0 = neiFace[n];
                label nPnt1 = neiFace[neiFace.fcIndex(n)];
                edge nEdge;             
                nEdge[0] = (nPnt0<nPnt1)?nPnt0:nPnt1;
                nEdge[1] = (nPnt0<nPnt1)?nPnt1:nPnt0;
                
                FixedList<label,4> key;
                if(oEdge[0]<nEdge[0])
                {
                    key[0] = oEdge[0]; key[1] = oEdge[1]; key[2] = nEdge[0]; key[3] = nEdge[1];
                }
                else
                {
                    key[0] = nEdge[0]; key[1] = nEdge[1]; key[2] = oEdge[0]; key[3] = oEdge[1];
                }
                
                auto iter = edgesToAddPntInd.find(key);
                if(iter==edgesToAddPntInd.end())
                {
                    edgesToAddPntInd.insert({key,addedPntsIndex});
                    addPnt.append(computeIntersectionPnt(oEdge,nEdge));
                    faceEdgeIntersectionsAddPntInd[o][n] = addedPntsIndex;
                    addedPntsIndex++;
                }
                else
                {
                    faceEdgeIntersectionsAddPntInd[o][n] = iter->second;
                }
            }
        }
    }

    bool doFacesEdgesIntersect = faceEdgesIntersect(faceEdgeIntersections);
    std::unique_ptr<std::vector<std::pair<label,label>>> matchingPntsOwnNei = matchingPoints(ownFace,neiFace);
    
    std::unordered_map<label,DynamicList<label>> pntNeighbors;
    for(label ownInd=0;ownInd<ownFace.size();ownInd++)
    {
        label currPnt = ownFace[ownInd];
        if(neiFace.which(currPnt)!=-1)
        {
            label prevPnt = ownFace.prevLabel(ownInd);
            if(neiFace.which(prevPnt)!=-1)
                pntNeighbors[currPnt].append(prevPnt);
            label nextPnt = ownFace.nextLabel(ownInd);
            if(neiFace.which(nextPnt)!=-1)
                pntNeighbors[currPnt].append(nextPnt);
        }
    }
    for(label o=0;o<faceEdgeIntersections.size();o++)
    {
        edge ownEdge = getEdge(ownFace,o);
        for(label n=0;n<faceEdgeIntersections[o].size();n++)
        {
            edge neiEdge = getEdge(neiFace,n);
            if(faceEdgeIntersections[o][n])
            {
                label intersecPnt = faceEdgeIntersectionsAddPntInd[o][n];
                label ownSecondInterse = -1;
                for(label n2=0;n2<faceEdgeIntersections[o].size();n2++)
                {
                    if(n2==n)
                        continue;
                    if(faceEdgeIntersections[o][n2])
                        ownSecondInterse = n2;
                }
                if(ownSecondInterse!=-1)
                {
                    label secOwnIntersecPnt = faceEdgeIntersectionsAddPntInd[o][ownSecondInterse];
                    pntNeighbors[intersecPnt].append(secOwnIntersecPnt);
                }
                else
                {
                    label otherOwnPnt=-1;
                    if(ownFace.which(ownEdge.start())!=-1 && neiFace.which(ownEdge.start())!=-1)
                        otherOwnPnt = ownEdge.start();
                    else if(ownFace.which(ownEdge.end())!=-1 && neiFace.which(ownEdge.end())!=-1)
                        otherOwnPnt = ownEdge.end();
                    else
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                    pntNeighbors[intersecPnt].append(otherOwnPnt);
                }
                
                label neiSecondIntersec = -1;
                for(label o2=0;o2<faceEdgeIntersections.size();o2++)
                {
                    if(o2==o)
                        continue;
                    if(faceEdgeIntersections[o2][n])
                        neiSecondIntersec = o2;
                }
                if(neiSecondIntersec!=-1)
                {
                    label secNeiIntersecPnt = faceEdgeIntersectionsAddPntInd[neiSecondIntersec][n];
                    pntNeighbors[intersecPnt].append(secNeiIntersecPnt);
                }
                else
                {
                    label otherNeiPnt=-1;
                    if(ownFace.which(neiEdge.start())!=-1 && neiFace.which(neiEdge.start())!=-1)
                        otherNeiPnt = neiEdge.start();
                    else if(ownFace.which(neiEdge.end())!=-1 && neiFace.which(neiEdge.end())!=-1)
                        otherNeiPnt = neiEdge.end();
                    else
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                    pntNeighbors[intersecPnt].append(otherNeiPnt);
                }
            }
        }
    }
    
    std::unordered_set<label> treatedPnts;
    DynamicList<label> facePnts;
    auto iter = pntNeighbors.begin();
    if(iter->second.size()!=2)
        FatalErrorInFunction<<"Error!"<< exit(FatalError);
    facePnts.append(iter->first);
    treatedPnts.insert(iter->first);
    facePnts.append(iter->second[0]);
    treatedPnts.insert(iter->second[0]);
    
    while(true)
    {
        label currPnt = facePnts[facePnts.size()-1];
        label prevPnt = facePnts[facePnts.size()-2];
        
        iter = pntNeighbors.find(currPnt);
        if(iter==pntNeighbors.end())
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
        label nextPnt = -1;
        if(iter->second[0]==currPnt)
            nextPnt = iter->second[1];
        else if(iter->second[1]==currPnt)
            nextPnt = iter->second[0];
        else
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
        
        if(nextPnt==facePnts[0])
        {
            break;
        }
        else
        {
            if(treatedPnts.find(nextPnt)!=treatedPnts.end())
                FatalErrorInFunction<<"Error!"<< exit(FatalError);
            if(facePnts.size()>20)
                FatalErrorInFunction<<"Face too large something is wrong here!"<< exit(FatalError);
            facePnts.append(nextPnt);
            treatedPnts.insert(nextPnt);
        }
    }
    
    for(label pnt : facePnts)
    {
        auto iterTr = treatedPnts.find(pnt);
        if(iterTr==treatedPnts.end())
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
        treatedPnts.erase(iterTr);
        
        auto iterMap = pntNeighbors.find(pnt);
        if(iterMap==pntNeighbors.end())
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
        pntNeighbors.erase(iterMap);
    }
    
    if(pntNeighbors.size()!=0 || treatedPnts.size()!=0)
        FatalErrorInFunction<<"Error!"<< exit(FatalError);
    
    /*
    if(doFacesEdgesIntersect)
    {
        DynamicList<label> facePnts;
        std::unordered_set<label> facePntsSet;
        List<List<bool>> visitedIntersec(ownFace.size(),List<bool>(neiFace.size(),false));
        label currIntersecEdgeOwn = -1;
        enum EdgeType {OwnEdgeInt, NeiEdgeInt, NoIntEdge};
        EdgeType currentEdge = EdgeType::OwnEdgeInt;
        bool endIntersec = false;
        for(label o=0;o<faceEdgeIntersections.size();o++)
        {
            label numIntersec = 0;
            for(bool intersec : faceEdgeIntersections[o])
                if(intersec)
                    numIntersec++;
            if(numIntersec>0)
            {
                if(numIntersec==1)
                    endIntersec=true;
                else
                    endIntersec=false;
                currIntersecEdgeOwn = o;
                break;
            }
        }
        if(currIntersecEdgeOwn==-1)
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
        
        label currIntersecEdgeNei = -1;
        for(label n=0;n<faceEdgeIntersections[currIntersecEdgeOwn].size();n++)
        {
            if(faceEdgeIntersections[currIntersecEdgeOwn][n])
            {
                currIntersecEdgeNei = n;
                visitedIntersec[currIntersecEdgeOwn][n] = true;
            }
        }
        if(currIntersecEdgeNei==-1)
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
        
        label pnt = faceEdgeIntersectionsAddPntInd[currIntersecEdgeOwn][currIntersecEdgeNei];
        facePnts.append(pnt);
        facePntsSet.insert(pnt);
        
        
        
        if(currentEdge != EdgeType::NoIntEdge)
        {
            if(endIntersec)
            {
                edge thisEdge;
                if(currentEdge==EdgeType::OwnEdgeInt)
                {
                    thisEdge = getEdge(ownFace,currIntersecEdgeOwn);
                }
                else
                {
                    thisEdge = getEdge(neiEdge,currIntersecEdgeNei);
                }
                if(ownFace.which(thisEdge.start())!=-1 &&
                   neiFace.which(thisEdge.start())!=-1)
                {
                    pnt = thisEdge.start();
                }
                else if(ownFace.which(thisEdge.end())!=-1 &&
                        neiFace.which(thisEdge.end())!=-1)
                {
                    pnt = thisEdge.end();
                }
                else
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                facePnts.append(pnt);
                facePntsSet.insert(pnt);
                
            }
            else
            {
            }
        }
        else
        {
            
        }
    }
    else
    {
        if(matchingPntsOwnNei->size()>=3)
        {
            face smallerFace = (ownFace.size()<neiFace.size()) ? ownFace : neiFace;
            face largerFace = (ownFace.size()<neiFace.size()) ? neiFace : ownFace;
            
            std::unordered_set<label> largerFaceSet(largerFace.begin(),largerFace.end());
            for(label pnt : smallerFace)
            {
                if(largerFaceSet.find(pnt)==largerFaceSet.end())
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
            }
            return smallerFace;
        }
        else
        { // Intersection is empty
        }
    }
    */
    return facePnts;
}

scalar Foam::cutCellFvMesh::computeCellSizeAndSizeLimit
(
    const List<cell>& cells,
    List<bool> tooSmallCells,
    List<scalar> cellSizes,
    const scalar limitPercentage,
    const scalar tooSmallFraction
)
{
    tooSmallCells.clear();
    tooSmallCells.resize(cells.size(),false);
    std::vector<scalar> cellSize(cells.size());
    for(label cellInd=0;cellInd<cells.size();cellInd++)
    {
        const cell& oneCell = cells[cellInd];
        cellSize[cellInd] = oneCell.mag(new_points,new_faces);
    }
    std::sort(cellSize.begin(),cellSize.end());
    
    scalar locMinCellSize = cellSize[cellSize.size()*limitPercentage];
    locMinCellSize *= tooSmallFraction;
    
    struct scalar_min
    {
        scalar operator()(const scalar a, const scalar b) const
        {
            return std::min<scalar>(a,b);
        }
    };
    scalar_min op;
    Pstream::gather(locMinCellSize,op);
    Pstream::scatter(locMinCellSize);
    
    cellSizes.setSize(cells.size());
    for(label cellInd=0;cellInd<cells.size();cellInd++)
    {
        const cell& oneCell = cells[cellInd];
        scalar cellSize = oneCell.mag(new_points,new_faces);
        cellSizes[cellInd] = cellSize;
        if(cellSize<locMinCellSize)
            tooSmallCells[cellInd]=true;
    }
    return locMinCellSize;
}

void Foam::cutCellFvMesh::findSmallCellMergeCandidate
(
    const List<cell>& cells,
    const List<bool> smallCell,
    List<DynamicList<std::pair<label,label>>>& smallCellMergeCand
)
{
    smallCellMergeCand.clear();
    smallCellMergeCand.resize(cells.size());
    for(label cellInd=0;cellInd<cells.size();cellInd++)
    {
        const cell& oneCell = cells[cellInd];
        scalar cellSize = oneCell.mag(new_points,new_faces);
        if(smallCell[cellInd])
        {
            for(label locFaceInd=0; locFaceInd<oneCell.size(); locFaceInd++)
            {
                label faceInd = oneCell[locFaceInd];
                if(faceInd>=new_neighbour.size() || faceInd>=new_owner.size())
                    FatalErrorInFunction<<"Error"<< exit(FatalError);

                label faceNeighborCell=-1;
                if(new_owner[faceInd]==cellInd)
                {
                    faceNeighborCell=new_neighbour[faceInd];
                }
                else if(new_neighbour[faceInd]==cellInd)
                {
                    faceNeighborCell=new_owner[faceInd];
                }
                else
                    FatalErrorInFunction<<"Error"<< exit(FatalError);
                
                if(faceNeighborCell!=-1)
                {
                    smallCellMergeCand[cellInd].append({locFaceInd,faceNeighborCell});
                }
            }
        }
    }
}

void Foam::cutCellFvMesh::testCellConvexity
(
    const List<cell>& cells,
    List<List<DynamicList<std::tuple<label,std::pair<label,label>,bool,scalar>>>>& cellFaceEdgeGraph,
    List<bool>& nonConvexCell
)
{
    cellFaceEdgeGraph.clear();
    cellFaceEdgeGraph.setSize(cells.size());
    nonConvexCell.clear();
    nonConvexCell.setSize(cells.size(),false);
    for(label cellInd=0;cellInd<cells.size();cellInd++)
    {
        const cell& oneCell = cells[cellInd];
        List<DynamicList<std::tuple<label,std::pair<label,label>,bool,scalar>>>& faceEdgeGraph = cellFaceEdgeGraph[cellInd];
        faceEdgeGraph.setSize(oneCell.size());
        for(label locFaceInd=0;locFaceInd<oneCell.size();locFaceInd++)
        {
            const label oneFaceInd = oneCell[locFaceInd];
            const face& oneFace = new_faces[oneFaceInd];
            for(const label pntIni : oneFace)
            {
                label locVertInd = oneFace.which(pntIni);
                if(locVertInd==-1)
                    FatalErrorInFunction<<"Error"<< exit(FatalError);
                label pntNext = oneFace.nextLabel(locVertInd);
                
                // Find the face connected to the two points pntIni,pntNext
                label j_Conn = -1;
                label smConnPnt = -1;
                label laConnPnt = -1;
                scalar cosAngle = 0;
                bool convex = false;
                for(label j=0;j<oneCell.size();j++)
                {
                    if(locFaceInd!=j)
                    {
                        const label otherFaceInd = oneCell[j];
                        const face& otherFace = new_faces[otherFaceInd];
                        label otherLocVertInd = otherFace.which(pntIni);
                        if(otherLocVertInd!=-1)
                        {
                            label otherPntNext = otherFace.nextLabel(otherLocVertInd);
                            label otherPntPrev = otherFace.prevLabel(otherLocVertInd);
                            if(otherPntNext==pntNext || otherPntPrev==pntNext)
                            {
                                if(j_Conn!=-1)
                                {
                                    Info<<"j_Conn:"<<j_Conn<<Foam::endl;
                                    Info<<"j:"<<j<<Foam::endl;
                                    FatalErrorInFunction<<"Error"<< exit(FatalError);
                                }
                                j_Conn = j;
                                if(pntIni<pntNext)
                                {
                                    smConnPnt = pntIni;
                                    laConnPnt = pntNext;
                                }
                                else if(pntIni>pntNext)
                                {
                                    smConnPnt = pntNext;
                                    laConnPnt = pntIni;
                                }
                                else
                                    FatalErrorInFunction<<"Error"<< exit(FatalError);
                                
                                label oneFaceOwnerCellInd = new_owner[oneFaceInd];
                                vector oneFaceNormal = oneFace.normal(new_points);
                                if(oneFaceOwnerCellInd!=cellInd)
                                    oneFaceNormal = -1*oneFaceNormal;
                                    
                                label otherFaceOwnerCellInd = new_owner[otherFaceInd];
                                vector otherFaceNormal = otherFace.normal(new_points);
                                if(otherFaceOwnerCellInd!=cellInd)
                                    otherFaceNormal = -1*otherFaceNormal;
                                
                                vector oneFaceCenter = oneFace.centre(new_points);
                                vector otherFaceCenter = otherFace.centre(new_points);
                                vector oneToOther = otherFaceCenter-oneFaceCenter;
                                vector otherToOne = oneFaceCenter-otherFaceCenter;
                                scalar oneToOtherToNormalAngle = oneToOther & otherFaceNormal;
                                scalar otherToOneToNormalAngle = otherToOne & oneFaceNormal;
                                
                                cosAngle = oneFaceNormal & otherFaceNormal;
                                if(oneToOtherToNormalAngle>0 && otherToOneToNormalAngle>0)
                                {
                                    convex = true;
                                }
                                else if(oneToOtherToNormalAngle<0 && otherToOneToNormalAngle<0)
                                {
                                    convex = false;
                                }
                                else
                                {
                                    scalar difference = std::abs(oneToOtherToNormalAngle-otherToOneToNormalAngle);
                                    if(difference<1e-10)
                                    {
                                        convex = true;
                                    }
                                    else
                                    {
                                        Info<<"cellInd:"<<cellInd<<Foam::endl;
                                        Info<<"oneFaceOwnerCellInd:"<<oneFaceOwnerCellInd<<Foam::endl;
                                        Info<<"otherFaceOwnerCellInd:"<<otherFaceOwnerCellInd<<Foam::endl;
                                        
                                        Info<<"locFaceInd:"<<locFaceInd<<" oneFace.points(new_points):"<<oneFace.points(new_points)<<Foam::endl;
                                        Info<<"j:"<<j<<" otherFace.points(new_points):"<<otherFace.points(new_points)<<Foam::endl;
                                        
                                        Info<<"oneFaceCenter:"<<oneFaceCenter<<Foam::endl;
                                        Info<<"otherFaceCenter:"<<otherFaceCenter<<Foam::endl;
                                        Info<<"oneToOther:"<<oneToOther<<Foam::endl;
                                        Info<<"otherToOne:"<<otherToOne<<Foam::endl;

                                        Info<<"otherFaceNormal:"<<otherFaceNormal<<Foam::endl;
                                        Info<<"oneFaceNormal:"<<oneFaceNormal<<Foam::endl;
                                        
                                        Info<<"otherFaceNormal:"<<otherFace.normal(new_points)<<Foam::endl;
                                        Info<<"oneFaceNormal:"<<oneFace.normal(new_points)<<Foam::endl;

                                        Info<<"oneToOtherToNormalAngle:"<<oneToOtherToNormalAngle<<Foam::endl;
                                        Info<<"otherToOneToNormalAngle:"<<otherToOneToNormalAngle<<Foam::endl;
                                        FatalErrorInFunction<<"Error"<< exit(FatalError);
                                    }
                                }
                            }
                        }
                    }
                }
                if(j_Conn==-1)
                    FatalErrorInFunction<<"Error"<< exit(FatalError);
                faceEdgeGraph[locFaceInd].append({j_Conn,{smConnPnt,laConnPnt},convex,cosAngle});
                nonConvexCell[cellInd] = true;
            }
            if(faceEdgeGraph[locFaceInd].size()!=oneFace.size())
                FatalErrorInFunction<<"Error"<< exit(FatalError);
        }
    }
}

void Foam::cutCellFvMesh::agglomerateSmallCells_MC33
(
    scalar partialThreeshold
)
{
    auto cellsPtr = createCells(new_faces,new_owner,new_neighbour);
    cellList& new_cells = *cellsPtr;
    
    Info<<"Test for closed cell!"<<Foam::endl;
    testCellClosure(new_faces,new_cells);
    
    Info<<"Test for correct normal direction!"<<Foam::endl;
    correctFaceNormal(new_faces,new_owner,new_neighbour,new_cells);
    
    //Testing for nonconvexivity in cell and correct them
    Info<<"Test for nonconvex cell!"<<Foam::endl;    
    List<List<DynamicList<std::tuple<label,std::pair<label,label>,bool,scalar>>>> cellFaceEdgeGraph;
    List<bool> nonConvexCell;
    testCellConvexity(new_cells,cellFaceEdgeGraph,nonConvexCell);
        
    Info<<"Test for cell sizes!"<<Foam::endl;
    Info<<"new_cells.size():"<<new_cells.size()<<Foam::endl;
    List<bool> smallCell;
    List<scalar> cellSizes;
    scalar locMinCellSize = computeCellSizeAndSizeLimit(new_cells,smallCell,cellSizes); 
        
    Info<<"Find merge cell candidates!"<<Foam::endl;
    List<DynamicList<std::pair<label,label>>> smallCellMergeCand;
    findSmallCellMergeCandidate(new_cells,smallCell,smallCellMergeCand);
    
    Info<<"Fix nonconvex cells"<<Foam::endl;
    List<std::unique_ptr<CellSplitData>> convexCorrectionData(new_cells.size());
    for(label cellInd=0;cellInd<new_cells.size();cellInd++)
    {
        Info<<"cellInd:"<<cellInd<<Foam::endl;

        label origCellIndex = cellMap[cellInd];
        if(origCellIndex==-1)
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
        
        Info<<"  origCellIndex:"<<origCellIndex<<Foam::endl;
        
        cell& thisCell = new_cells[cellInd];
        std::unique_ptr<List<scalar>> faceToCenterAngles = faceToCellCenterRelation(new_cells[cellInd],cellInd);
        
        Info<<"Got faceToCellCenterRelation"<<Foam::endl;

        MC33::MC33Cube& oneCube = mc33CutCellData[origCellIndex];
        
        Info<<"oneCube.cubeCase:"<<oneCube.cubeCase<<Foam::endl;
        
        if(oneCube.cell==-1)
            continue;
        
        Info<<"oneCube.cell:"<<oneCube.cell<<Foam::endl;
        
        if(oneCube.cell!=origCellIndex)
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
        
        Info<<"nonConvexCell.size():"<<nonConvexCell.size()<<Foam::endl;
        
        bool isNonConvex = nonConvexCell[cellInd];
        
        Info<<"Cases"<<Foam::endl;

        if(oneCube.cubeCase==MC33::c0)
            // Case 0 Triangle Number 0
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Can not be concave!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c1)
            // Case 1 Triangle Number 1
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Can not be concave!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c2 )
            // Case 2 Triangle Number 2
        {
            if(isNonConvex)
                convexCorrectionData[cellInd] = nonConvexCellSplitC2(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c31 )
            // Case 3.1 Triangle Number 2
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Can not be concave!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c32 )
            // Case 3.2 Triangle Number 4
        {
            convexCorrectionData[cellInd] = nonConvexCellSplitC32(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c411 )
            // Case 4.1.1 Triangle Number 2
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Can not be concave!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c412 )
            // Case 4.1.2 Triangle Number 6
        {
            if(isNonConvex)
                convexCorrectionData[cellInd] = nonConvexCellSplitC412(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c5 )
            // Case 5 Triangle Number 3
        {
            convexCorrectionData[cellInd] = nonConvexCellSplitC5(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c611 )
            // Case 6.1.1 Triangle Number 3
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Not yet implemented!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c612 )
            // Case 6.1.2 Triangle Number 7
        {
            if(!isNonConvex)
                FatalErrorInFunction<<"Must be concave!"<< exit(FatalError);
            convexCorrectionData[cellInd] = nonConvexCellSplitC612(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c62 )
            // Case 6.2 Triangle Number 5
        {
            if(!isNonConvex)
                FatalErrorInFunction<<"Must be concave!"<< exit(FatalError);
            convexCorrectionData[cellInd] = nonConvexCellSplitC62(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c71 )
            // Case 7.1 Triangle Number 3
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Can not be concave!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c72_720 || oneCube.cubeCase==MC33::c72_760 || oneCube.cubeCase==MC33::c72_800)
            // Case 7.2 Triangle Number 5
        {
            if(isNonConvex)
                convexCorrectionData[cellInd] = nonConvexCellSplitC72(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c73 )
            // Case 7.3 Triangle Number 9
        {
            if(!isNonConvex)
                FatalErrorInFunction<<"Must be concave!"<< exit(FatalError);
            convexCorrectionData[cellInd] = nonConvexCellSplitC73(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c741 )
            // Case 7.4.1 Triangle Number 5
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Can not be concave!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c742 )
            // Case 7.4.2 Triangle Number 9
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Can not be concave!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c8 )
            // Case 8 Triangle Number 2
        {            
            if(isNonConvex)
                convexCorrectionData[cellInd] = nonConvexCellSplitC8(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c9 )
            // Case 9 Triangle Number 4
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Can not be concave!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c1011_1190 || oneCube.cubeCase==MC33::c1011_1202)
            // Case 10.1.1 Triangle Number 4
        {
            if(isNonConvex)
                convexCorrectionData[cellInd] = nonConvexCellSplitC1011(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c1012_1214 || oneCube.cubeCase==MC33::c1012_1238)
            // Case 10.1.2 Triangle Number 8
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Can not be concave!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c102_1262 || oneCube.cubeCase==MC33::c102_1286)
            // Case 10.2 Triangle Number 8
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Not yet implemented!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c11 )
            // Case 11 Triangle Number 4
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Not yet implemented!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c14 )
            // Case 14 Triangle Number 4
        {
            if(!isNonConvex)
                FatalErrorInFunction<<"Must be concave!"<< exit(FatalError);
            convexCorrectionData[cellInd] = nonConvexCellSplitC14(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c1211_1358 || oneCube.cubeCase==MC33::c1211_1406)
            // Case 12.1.1 Triangle Number 4
        {
            convexCorrectionData[cellInd] = nonConvexCellSplitC1211(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c1212_1454 || oneCube.cubeCase==MC33::c1212_1550)
            // Case 12.1.2 Triangle Number 8
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Not yet implemented!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c122_1646 || oneCube.cubeCase==MC33::c122_1742)
            // Case 12.2 Triangle Number 8
        {
            if(!isNonConvex)
                FatalErrorInFunction<<"Must be concave!"<< exit(FatalError);
            convexCorrectionData[cellInd] = nonConvexCellSplitC122(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c131 )
            // Case 13.1 Triangle Number 4
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Not yet implemented!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c132_MostNeg || oneCube.cubeCase==MC33::c132_MostPos)
            // Case 13.2 Triangle Number 6
        {
            convexCorrectionData[cellInd] = nonConvexCellSplitC132(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c133_MostNeg || oneCube.cubeCase==MC33::c133_MostPos)
            // Case 13.3 Triangle Number 10
        {
            if(!isNonConvex)
                FatalErrorInFunction<<"Must be concave!"<< exit(FatalError);
            convexCorrectionData[cellInd] = nonConvexCellSplitC133(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c134 )
            // Case 13.4 Triangle Number 12
        {
            if(!isNonConvex)
                FatalErrorInFunction<<"Must be concave!"<< exit(FatalError);
            convexCorrectionData[cellInd] = nonConvexCellSplitC134(thisCell,cellInd,std::move(faceToCenterAngles));
        }
        else if(oneCube.cubeCase==MC33::c1352_2206 || oneCube.cubeCase==MC33::c1352_2246)
            // Case 13.5.2 Triangle Number 10
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Can not be concave!"<< exit(FatalError);
        }
        else if(oneCube.cubeCase==MC33::c1351)
            // Case 13.5.1 Triangle Number 6
        {
            if(isNonConvex)
                FatalErrorInFunction<<"Can not be concave!"<< exit(FatalError);
        }
        else
            FatalErrorInFunction<<"Error!"<< exit(FatalError);
    }
    
    FatalErrorInFunction<<"Temp Stop!"<< exit(FatalError);
    
    std::unordered_map<label,DynamicList<label>> pointFaces;
    for(label faceInd=0; faceInd<new_faces.size(); faceInd++)
    {
        for(label pntLocalInd=0; pntLocalInd<new_faces[faceInd].size(); pntLocalInd++)
        {
            pointFaces[new_faces[faceInd][pntLocalInd]].append(faceInd);
        }
    }
    
    //correct added points indices in CellSplitData
    //pointField new_points;
    DynamicList<vector> new_points_List;
    new_points_List.append(this->new_points);
    DynamicList<vector> added_points;
    std::unordered_map<label,DynamicList<std::pair<label,edge>>> faceToAddedPnts;    
    for(label cellInd=0;cellInd<new_cells.size();cellInd++)
    {
        if(convexCorrectionData[cellInd])
        {
            CellSplitData& cellSplit = *(convexCorrectionData[cellInd]);
            if(cellSplit.originalCell != cellMap[cellInd])
                FatalErrorInFunction<<"Error!"<< exit(FatalError);
            if(cellSplit.addedPointsLimit == -1)
            {
                if(cellSplit.addedPoints.size()!=0)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
            }
            else
            // cell has added points
            {
                if(cellSplit.addedPointsLimit!=new_points.size())
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                if(cellSplit.addedPoints.size()==0)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                std::unordered_set<label> thisCellFaces;
                thisCellFaces.insert(new_cells[cellInd].begin(),new_cells[cellInd].end());
                
                //Compute point index update
                std::unordered_map<label,label> pntIndTransfer;
                for(label i=0; i<cellSplit.addedPoints.size(); i++)
                {
                    pntIndTransfer[i+cellSplit.addedPointsLimit] = new_points.size()+added_points.size()+i;
                }
                
                //Correct added faces point indices
                for(face& addedFace : cellSplit.addedFace)
                {
                    for(label& pntInd : addedFace)
                    {
                        if(!(pntInd<cellSplit.addedPointsLimit))
                        {
                            auto map = pntIndTransfer.find(pntInd);
                            if(map!=pntIndTransfer.end())
                            {
                                pntInd = map->second;
                            }
                            else
                                FatalErrorInFunction<<"Error!"<< exit(FatalError);
                        }
                    }
                }
                std::unordered_map<label,DynamicList<std::pair<face,DynamicList<label>>>> origFaceAddPntToCutFaces;
                //Correct split faces point indices
                for(std::pair<face,label>& splitFaceData : cellSplit.splittedFaces)
                {
                    label origFaceInd = splitFaceData.second;
                    face& splitFace = splitFaceData.first;
                    DynamicList<label> addedPntInds;
                    for(label k=0; k<splitFace.size(); k++)
                    {
                        label& pntInd = splitFace[k];
                        if(!(pntInd<cellSplit.addedPointsLimit))
                        {
                            auto map = pntIndTransfer.find(pntInd);
                            if(map!=pntIndTransfer.end())
                            {
                                pntInd = map->second;
                                addedPntInds.append(map->first);
                            }
                            else
                                FatalErrorInFunction<<"Error!"<< exit(FatalError);
                        }
                    }
                    std::pair<face,DynamicList<label>> split = {splitFace,addedPntInds};
                    auto& dat = origFaceAddPntToCutFaces[origFaceInd];
                    dat.append(split);
                }
                
                // record added points for old faces
                if(origFaceAddPntToCutFaces.size() >= static_cast<long unsigned int>(new_cells[cellInd].size()))
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                for(auto iter=origFaceAddPntToCutFaces.cbegin(); iter!=origFaceAddPntToCutFaces.end(); iter++)
                {
                    label faceInd = iter->first;
                    face& thisFace = new_faces[faceInd];
                    
                    const DynamicList<std::pair<face,DynamicList<label>>>& cutFacesWithAddPnts = iter->second;
                    std::unordered_set<label> addedPntsForFace;
                    for(const std::pair<face,DynamicList<label>>& cutFace : cutFacesWithAddPnts)
                    {
                        addedPntsForFace.insert(cutFace.second.begin(),cutFace.second.end());
                    }
                    for(auto iter=addedPntsForFace.begin(); iter!=addedPntsForFace.end(); iter++)
                    {
                        label locAddedPntInd = *iter;
                        std::pair<vector,edge>& addPntDuo = cellSplit.addedPoints[locAddedPntInd];
                        edge& edgeOfAddedPnts = addPntDuo.second;
                        faceToAddedPnts[faceInd].append({pntIndTransfer[locAddedPntInd],edgeOfAddedPnts});
                        
                        label edgeStartPnt = edgeOfAddedPnts.start();
                        auto iterPntFaces = pointFaces.find(edgeStartPnt);
                        if(iterPntFaces==pointFaces.end())
                            FatalErrorInFunction<<"Error!"<< exit(FatalError);
                        DynamicList<label>& facesOfEdgeStartPnt = iterPntFaces->second;
                        if(facesOfEdgeStartPnt.size()<3)
                            FatalErrorInFunction<<"Error!"<< exit(FatalError);
                        std::unordered_set<label> facesOfEdgeStartPntSet;
                        facesOfEdgeStartPntSet.insert(facesOfEdgeStartPnt.begin(),facesOfEdgeStartPnt.end());
                        
                        label edgeEndPnt = edgeOfAddedPnts.end();
                        iterPntFaces = pointFaces.find(edgeEndPnt);
                        if(iterPntFaces==pointFaces.end())
                            FatalErrorInFunction<<"Error!"<< exit(FatalError);
                        DynamicList<label>& facesOfEdgeEndPnt = iterPntFaces->second;
                        
                        DynamicList<label> edgeAdjacentFaces;
                        for(label faceIndx : facesOfEdgeEndPnt)
                        {
                            auto iter = facesOfEdgeStartPntSet.find(faceIndx);
                            if(faceIndx!=faceInd && iter!=facesOfEdgeStartPntSet.end())
                            {
                                edgeAdjacentFaces.append(faceIndx);
                            }
                        }
                        if(edgeAdjacentFaces.size()<2)
                            FatalErrorInFunction<<"Error!"<< exit(FatalError);
                        
                        for(label faceIndx : edgeAdjacentFaces)
                        {
                            faceToAddedPnts[faceIndx].append({pntIndTransfer[locAddedPntInd],edgeOfAddedPnts});
                        }
                    }
                }
                
                for(std::pair<vector,edge>& addPnt : cellSplit.addedPoints)
                {
                    added_points.append(addPnt.first);
                }
            }
        }
    }
    new_points_List.append(added_points);    
    for(label faceInd=0;faceInd<new_faces.size();faceInd++)
    {
        auto iter = faceToAddedPnts.find(faceInd);
        if(iter!=faceToAddedPnts.end())
        {
            face& thisFace = new_faces[faceInd];
            DynamicList<std::pair<label,edge>>& cutPnts = iter->second;
            List<std::vector<label>> intersectionPnts(thisFace.size());
            for(std::pair<label,edge>& onePnt : cutPnts)
            {
                bool oneMatchFound = false;
                label thisPnt = onePnt.first;
                edge& thisEdge = onePnt.second;
                for(label i0=0;i0<thisFace.size();i0++)
                {
                    label i1 = thisFace.fcIndex(i0);
                    label p0 = thisFace[i0];
                    label p1 = thisFace[i1];
                    if(thisEdge.otherVertex(p0)==p1 || thisEdge.otherVertex(p1)==p0)
                    {
                        if(oneMatchFound)
                            FatalErrorInFunction<<"Error!"<< exit(FatalError);
                        oneMatchFound=true;
                        intersectionPnts[i0].push_back(thisPnt);
                    }
                }
                if(!oneMatchFound)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
            }
            for(label i=0;i<intersectionPnts.size();i++)
            {
                vector startPnt = new_points[thisFace[i]];
                std::vector<label>& intPnt = intersectionPnts[i];
                auto sortFunc = [&](label a, label b)
                {
                    scalar distA = norm2(startPnt-new_points[a]);
                    scalar distB = norm2(startPnt-new_points[b]);
                    return distA<distB;
                };
                std::sort(intPnt.begin(),intPnt.end(),sortFunc);
            }
            DynamicList<label> newFace;
            for(label i=0;i<thisFace.size();i++)
            {
                newFace.append(thisFace[i]);
                for(unsigned int j=0;j<intersectionPnts[i].size();j++)
                {
                    newFace.append(intersectionPnts[i][j]);
                }
            }
            thisFace = face(newFace);
        }
    }

    //Assign split and original faces to own and neighbor faces
    //std::unordered_set<label> removedFaces;
    //DynamicList<std::tuple<face,label,label>> addedFaces;
    //std::unordered_set<label> removedCell;
    DynamicList<label> cellIndStart({0});
    std::unordered_map<label,DynamicList<std::tuple<label,label,label>>> splittedFacesByNeighborCell;
    std::unordered_map<label,DynamicList<std::tuple<label,label,label>>> splittedFacesByOwnerCell;
    std::unordered_map<label,DynamicList<std::tuple<label,label,label>>> originalFacesByNeighborCell;
    std::unordered_map<label,DynamicList<std::tuple<label,label,label>>> originalFacesByOwnerCell;
    for(label cellInd=0; cellInd<new_cells.size(); cellInd++)
    {
        label cellMultiple = 1;
        if(convexCorrectionData[cellInd])
        {
            CellSplitData& cellSplit = *(convexCorrectionData[cellInd]);
            if(cellSplit.originalCell != cellMap[cellInd])
                FatalErrorInFunction<<"Error!"<< exit(FatalError);
            
            for(label locCell=0;locCell<cellSplit.cells.size();locCell++)
            {
                CellFaces& oneCell = cellSplit.cells[locCell];
                
                // handle original face indexes
                for(label locOFaceInd=0;locOFaceInd<oneCell.originalFaceInds.size();locOFaceInd++)
                {
                    label oFace = oneCell.originalFaceInds[locOFaceInd];
                    if(new_owner[oFace]==cellInd)
                    {
                        originalFacesByOwnerCell[oFace].append({cellInd,locCell,locOFaceInd});
                    }
                    else if(new_neighbour[oFace]==cellInd)
                    {
                        originalFacesByNeighborCell[oFace].append({cellInd,locCell,locOFaceInd});
                    }
                    else
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                }
                
                // handle splitted face indexes
                for(label locSpFaceInd=0;locSpFaceInd<oneCell.splittedFaceInds.size();locSpFaceInd++)
                {
                    label faceInds = oneCell.splittedFaceInds[locSpFaceInd];
                    std::pair<face,label>& splittedFace = cellSplit.splittedFaces[faceInds];
                    label origFace = splittedFace.second;
                    if(new_owner[origFace]==cellInd)
                    {
                        splittedFacesByOwnerCell[origFace].append({cellInd,locCell,locSpFaceInd});
                    }
                    else if(new_neighbour[origFace]==cellInd)
                    {
                        splittedFacesByNeighborCell[origFace].append({cellInd,locCell,locSpFaceInd});
                    }
                    else
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                }
                
                // handle replacing face indexes
                for(label locRepFaceInd=0;locRepFaceInd<oneCell.replacingFaceInds.size();locRepFaceInd++)
                {
                    label faceInds = oneCell.replacingFaceInds[locRepFaceInd];
                    const std::tuple<face,label,std::pair<bool,label>,bool>& replacingFace = cellSplit.replacingFaces[faceInds];

                    label replacingFaceInd = std::get<1>(replacingFace); // replacingFace
                    const std::pair<bool,label> replacedFaceData = std::get<2>(replacingFace);
                    bool replacedFaceIsNonAdded = replacedFaceInd.first;
                    label replacedFaceInd = replacedFaceInd.second;

                    if(new_owner[replacedFaceInd]==cellInd)
                    {
                        splittedFacesByOwnerCell[replacedFaceInd].append({cellInd,locCell,locSpFaceInd});
                    }
                    else if(new_neighbour[replacedFaceInd]==cellInd)
                    {
                        splittedFacesByNeighborCell[replacedFaceInd].append({cellInd,locCell,locSpFaceInd});
                    }
                    else
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                }
                


            }
            cellMultiple = cellSplit.cells.size();
        }
        cellIndStart.append(cellIndStart.last()+cellMultiple);
    }
    
    //
    List<DynamicList<face>> new_faces(this->new_faces.size());
    List<DynamicList<std::pair<label,label>>> new_owner(this->new_faces.size());
    List<DynamicList<std::pair<label,label>>> new_neighbour(this->new_faces.size());
    for(label faceInd=0; faceInd<new_faces.size(); faceInd++)
    {
        auto iterSpOwn = splittedFacesByOwnerCell.find(faceInd);
        auto iterSpNei = splittedFacesByNeighborCell.find(faceInd);
        auto iterOOwn = originalFacesByOwnerCell.find(faceInd);
        auto iterONei = originalFacesByNeighborCell.find(faceInd);
        DynamicList<std::tuple<label,label,label>>* singleSideCut;
        DynamicList<std::tuple<label,label,label>>* singleSideO;
        bool iterSpOwnB = iterSpOwn!=splittedFacesByOwnerCell.end();
        bool iterSpNeiB = iterSpNei!=splittedFacesByNeighborCell.end();
        bool iterOOwnB = iterOOwn!=originalFacesByOwnerCell.end();
        bool iterONeiB = iterONei!=originalFacesByNeighborCell.end();        

        std::pair<label,label> ownCell = {this->new_owner[faceInd],-1};
        std::pair<label,label> neiCell = {this->new_neighbour[faceInd],-1};
        if(!iterSpOwnB && !iterSpNeiB)
        {
            if(iterOOwnB)
            {
                singleSideO = &(iterOOwn->second);
                if(singleSideO->size()!=1)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                std::tuple<label,label,label>& oneSideO = (*singleSideO)[0];
                label cellInd = std::get<0>(oneSideO);
                if(cellInd!=ownCell.first)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                label locCell = std::get<1>(oneSideO);
                ownCell.second = locCell;
            }
            if(iterONeiB)
            {
                singleSideO = &(iterONei->second);
                if(singleSideO->size()!=1)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                std::tuple<label,label,label>& oneSideO = (*singleSideO)[0];
                label cellInd = std::get<0>(oneSideO);
                if(cellInd!=neiCell.first)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                label locCell = std::get<1>(oneSideO);
                neiCell.second = locCell;
            }
            new_faces[faceInd].append(this->new_faces[faceInd]);
            new_owner[faceInd].append(ownCell);
            new_neighbour[faceInd].append(neiCell);
        }
        else if(iterSpOwnB && iterSpNeiB)
        {
            if(iterOOwnB || iterONeiB)
                FatalErrorInFunction<<"Error!"<< exit(FatalError);
            
            edgToIntPntIndMap edgesToAddPntInd;
            
            const face& totalFace = this->new_faces[faceInd];
            
            if(iterSpOwn->first!=faceInd || iterSpNei->first!=faceInd)
                FatalErrorInFunction<<"Error!"<< exit(FatalError);            
            DynamicList<std::tuple<label,label,label>>& sideOwnSp = iterSpOwn->second;
            DynamicList<std::tuple<label,label,label>>& sideNeiSp = iterSpNei->second;

            for(std::tuple<label,label,label> spOwn : sideOwnSp)
            {
                label cellInd = std::get<0>(spOwn);
                if(cellInd!=this->new_owner[faceInd])
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                
                CellSplitData& cellSplitOwn = *(convexCorrectionData[cellInd]);
                label locCellOwn = std::get<1>(spOwn);
                CellFaces& oneCellOwn = cellSplitOwn.cells[locCellOwn];
                label locSpFaceIndOwn = std::get<2>(spOwn);
                label splitFaceIndOwn = oneCellOwn.splittedFaceInds[locSpFaceIndOwn];
                if(cellSplitOwn.splittedFaces[splitFaceIndOwn].second!=faceInd)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                
                const face& spOwnFace = cellSplitOwn.splittedFaces[splitFaceIndOwn].first;
                
                for(std::tuple<label,label,label> spNei : sideNeiSp)
                {
                    label cellInd = std::get<0>(spOwn);
                    if(cellInd!=this->new_neighbour[faceInd])
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                    
                    CellSplitData& cellSplitNei = *(convexCorrectionData[cellInd]);
                    label locCellNei = std::get<1>(spNei);
                    CellFaces& oneCellNei = cellSplitNei.cells[locCellNei];
                    label locSpFaceIndNei = std::get<2>(spNei);
                    label splitFaceIndNei = oneCellNei.splittedFaceInds[locSpFaceIndNei];
                    if(cellSplitNei.splittedFaces[splitFaceIndNei].second!=faceInd)
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                    
                    const face& spNeiFace = cellSplitNei.splittedFaces[splitFaceIndNei].first;
                    
                    List<label> intersecFacePnts = faceIntersection(totalFace,spOwnFace,spNeiFace,edgesToAddPntInd,new_points_List);
                    face intersecFace(intersecFacePnts);                    
                    
                    new_faces[faceInd].append(intersecFace);
                    std::pair<label,label> ownCell = {this->new_owner[faceInd],locCellOwn};
                    std::pair<label,label> neiCell = {this->new_neighbour[faceInd],locCellNei};
                    new_owner[faceInd].append(ownCell);
                    new_neighbour[faceInd].append(neiCell);
                }
            }
        }
        else
        {
            std::pair<label,label> othCell;
            if(iterSpOwnB)
            {
                othCell = {this->new_neighbour[faceInd],-1};
                if(iterOOwnB)
                {
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                }
                if(iterONeiB)
                {
                    singleSideO = &(iterONei->second);
                    if(singleSideO->size()!=1)
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                    std::tuple<label,label,label>& oneSideO = (*singleSideO)[0];
                    label cellInd = std::get<0>(oneSideO);
                    if(cellInd!=othCell.first)
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                    label locCell = std::get<1>(oneSideO);
                    othCell.second = locCell;
                }
                singleSideCut = &(iterSpOwn->second);
            }
            else
            {
                othCell = {this->new_owner[faceInd],-1};
                if(iterOOwnB)
                {
                    singleSideO = &(iterONei->second);
                    if(singleSideO->size()!=1)
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                    std::tuple<label,label,label>& oneSideO = (*singleSideO)[0];
                    label cellInd = std::get<0>(oneSideO);
                    if(cellInd!=othCell.first)
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                    label locCell = std::get<1>(oneSideO);
                    othCell.second = locCell;
                }
                if(iterONeiB)
                {
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                }
                singleSideCut = &(iterSpNei->second);
            }

            for(std::tuple<label,label,label>& oneSplitFace : *singleSideCut)
            {
                label cellInd = std::get<0>(oneSplitFace);
                label locCell = std::get<1>(oneSplitFace);
                label locSpFaceInd = std::get<2>(oneSplitFace);
                CellSplitData& cellSplit = *(convexCorrectionData[cellInd]);
                CellFaces& oneCell = cellSplit.cells[locCell];
                label splitFaceInd = oneCell.splittedFaceInds[locSpFaceInd];
                std::pair<face,label>& splitFace = cellSplit.splittedFaces[splitFaceInd];
                if(splitFace.second!=faceInd)
                    FatalErrorInFunction<<"Error!"<< exit(FatalError);
                
                new_faces[faceInd].append(splitFace.first);
                std::pair<label,label> ownCell = {this->new_owner[faceInd],-1};
                std::pair<label,label> neiCell = {this->new_neighbour[faceInd],-1};
                if(iterSpOwnB)
                {
                    ownCell.second = locCell;
                    neiCell = othCell;
                }
                else
                {
                    ownCell = othCell;
                    neiCell.second = locCell;
                }
                new_owner[faceInd].append(ownCell);
                new_neighbour[faceInd].append(neiCell);
            }
        }
    }
    
    pointField new_points;
    new_points.append(new_points_List);
    
    DynamicList<face> new_faces_flat;
    DynamicList<label> new_owner_flat;
    DynamicList<label> new_neighbour_flat;
    
    for(label i=0; i<new_faces.size(); i++)
    {
        for(label j=0; j<new_faces[i].size(); j++)
        {
            new_faces_flat.append(new_faces[i][j]);
            
            std::pair<label,label>& owner_dual = new_owner[i][j];
            label owner = cellIndStart[owner_dual.first];
            if(owner_dual.second!=-1)
            {
                owner+=owner_dual.second;
            }
            new_owner_flat.append(owner);
            
            std::pair<label,label>& neighbour_dual = new_neighbour[i][j];
            label neighbour = cellIndStart[neighbour_dual.first];
            if(neighbour_dual.second!=-1)
            {
                neighbour+=neighbour_dual.second;
            }
            new_neighbour_flat.append(neighbour);
        }
    }
    
    cellsPtr = createCells(new_faces_flat,new_owner_flat,new_neighbour_flat);
    new_cells = *cellsPtr;

    testCellClosure(new_faces_flat,new_cells);
    correctFaceNormal(new_faces_flat,new_owner_flat,new_neighbour_flat,new_cells);
    
    locMinCellSize =  computeCellSizeAndSizeLimit(new_cells,smallCell,cellSizes);    
    findSmallCellMergeCandidate(new_cells,smallCell,smallCellMergeCand);
    
    List<DynamicList<std::tuple<label,label,scalar,scalar>>> orderedMergeCand(new_cells.size());
    for(label cellInd=0; cellInd<new_cells.size(); cellInd++)
    {
        if(smallCell[cellInd])
        {
            scalar thisCellSize = cellSizes[cellInd];
            cell& thisCell = new_cells[cellInd];
            vector thisCellCentre = thisCell.centre(new_points,new_faces_flat);
            DynamicList<std::pair<label,label>>& mergeCandidates = smallCellMergeCand[cellInd];
            
            DynamicList<std::tuple<label,label,scalar,scalar>> mergeCandLargeEnough;            
            DynamicList<std::tuple<label,label,scalar,scalar>> mergeCandTooSmall;
            
            for(std::pair<label,label>& mergeCand : mergeCandidates)
            {
                label mergeCellInd = mergeCand.second;
                scalar mergeCellSize = cellSizes[mergeCellInd];
                cell& mergeCell = new_cells[mergeCellInd];
                vector mergeCellCentre = mergeCell.centre(new_points,new_faces_flat);
                
                scalar minAngle = std::numeric_limits<scalar>::max();
                scalar maxAngle = std::numeric_limits<scalar>::min();
                for(label locCellInd=0; locCellInd<thisCell.size(); locCellInd++)
                {
                    if(locCellInd==mergeCand.first)
                        continue;
                    
                    label locCellFaceInd = thisCell[locCellInd];
                    face& locCellFace = new_faces_flat[locCellFaceInd];
                    
                    vector locCellFaceCentre = locCellFace.centre(new_points);
                    vector mergeCellCentreToFaceCentre = locCellFaceCentre - mergeCellCentre;
                    
                    vector locCellFaceNormal = locCellFace.normal(new_points);
                    if(new_owner_flat[locCellFaceInd]==cellInd)
                    {
                        locCellFaceNormal *= -1;
                    }
                    else if(new_neighbour_flat[locCellFaceInd]==cellInd)
                    {
                    }
                    else
                        FatalErrorInFunction<<"Error!"<< exit(FatalError);
                    
                    scalar locCellFaceNormalLen = norm2(locCellFaceNormal);
                    scalar mergeCellCentreToFaceCentreLen = norm2(mergeCellCentreToFaceCentre);
                    scalar prod = locCellFaceNormal & mergeCellCentreToFaceCentre;
                    scalar angle = prod / (locCellFaceNormalLen*mergeCellCentreToFaceCentreLen);
                    
                    minAngle = std::min(minAngle,angle);
                    maxAngle = std::max(maxAngle,angle);
                }
                
                if(mergeCellSize+thisCellSize>=locMinCellSize)
                {
                    mergeCandLargeEnough.append({mergeCand.first,mergeCand.second,minAngle,maxAngle});
                }
                else
                {
                    mergeCandTooSmall.append({mergeCand.first,mergeCand.second,minAngle,maxAngle});
                }
            }

            auto compare = [](std::tuple<label,label,scalar,scalar> a, std::tuple<label,label,scalar,scalar> b)
            {
                return std::get<2>(a) > std::get<2>(b);
            };
            std::sort(mergeCandLargeEnough.begin(),mergeCandLargeEnough.end(),compare);
            std::sort(mergeCandTooSmall.begin(),mergeCandTooSmall.end(),compare);
            
            if(mergeCandLargeEnough.size()>0)
            {
                orderedMergeCand[cellInd] = mergeCandLargeEnough;
            }
            else
            {
                orderedMergeCand[cellInd] = mergeCandTooSmall;
            }
        }
    }
    
    // Decide for correct merge cell
    std::unordered_map<label,DynamicList<std::pair<label,label>>> cellToMergeCells;
    std::unordered_map<label,label> removedCellsToNewInd;
    std::unordered_set<label> removedFaces;
    for(label cellInd=0; cellInd<new_cells.size(); cellInd++)
    {
        if(smallCell[cellInd])
        {
            if(orderedMergeCand[cellInd].size()>0)
            {
                label locMergeFace = std::get<0>(orderedMergeCand[cellInd][0]);
                label mergeCell = std::get<1>(orderedMergeCand[cellInd][0]);
                label mergeFaceInd = new_cells[cellInd][locMergeFace];
                cellToMergeCells[mergeCell].append({mergeFaceInd,cellInd});
                removedCellsToNewInd[cellInd] = mergeCell;
                removedFaces.insert(mergeFaceInd);
            }
        }
    }
    
    List<label> cellIndUpdate(new_cells.size());
    label newCellInd = 0;
    for(label cellInd=0; cellInd<new_cells.size(); cellInd++)
    {
        auto iter = removedCellsToNewInd.find(cellInd);
        if(iter!=removedCellsToNewInd.end())
        {
            cellIndUpdate[cellInd] = iter->second;
        }
        else
        {
            cellIndUpdate[cellInd] = newCellInd;
            newCellInd++;
        }
    }
    
    if(newCellInd+1+removedFaces.size() != new_faces_flat.size())
        FatalErrorInFunction<<"Error!"<< exit(FatalError);
    
    this->new_faces.setSize(newCellInd+1);
    this->new_owner.setSize(newCellInd+1,-1);
    this->new_neighbour.setSize(newCellInd+1,-1);
    
    newCellInd=0;
    for(label faceInd=0;faceInd<new_faces_flat.size();faceInd++)
    {
        auto iter = removedFaces.find(faceInd);
        if(iter==removedFaces.end())
        {
            this->new_faces[newCellInd] = new_faces_flat[faceInd];
            this->new_owner[newCellInd] = new_owner_flat[faceInd];
            this->new_neighbour[newCellInd] = new_neighbour_flat[faceInd];
            this->new_owner[newCellInd] = cellIndUpdate[this->new_owner[newCellInd]];
            if(this->new_neighbour[newCellInd]!=-1)
            {
                this->new_neighbour[newCellInd] = cellIndUpdate[this->new_neighbour[newCellInd]];
            }
            newCellInd++;
        }
    }

    FatalErrorInFunction<<"Temp Stop!"<< exit(FatalError);
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

/*
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
        FatalErrorInFunction
        << "Cell size problem"
        << exit(FatalError);  
    }
}
*/

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
