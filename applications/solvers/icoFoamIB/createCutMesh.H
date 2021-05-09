Foam::Info
    << "Create cutCellMesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

Foam::Info<<"CUTCELLFVMESH"<<Foam::endl;    
    
int testdegree = 2;
scalarList knots(7);
scalarList weights(4);
List<vector> controlPoints(4);
List<std::shared_ptr<Nurbs>> items(0);
    
knots[0]=0; knots[1]=0; knots[2]=0; knots[3]=0.5;
knots[4]=1; knots[5]=1; knots[6]=1;    
weights[0] = 1;    weights[1] = 1;   weights[2] = 1;    weights[3] = 1; 
controlPoints[0]=vector(1,0,1); controlPoints[1]=vector(1,0,-0.1);
controlPoints[2]=vector(3,0,0.1); controlPoints[3]=vector(3,0,-1);
items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.2,4)));

knots[0]=0; knots[1]=0; knots[2]=0; knots[3]=0.5;
knots[4]=1; knots[5]=1; knots[6]=1;    
weights[0] = 1;    weights[1] = 1;   weights[2] = 1;    weights[3] = 1; 
controlPoints[0]=vector(2.4,0,1); controlPoints[1]=vector(2,0,0.33);
controlPoints[2]=vector(2,0,-0.33); controlPoints[3]=vector(1.6,0,-1);
items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.15,4)));

Info<<"Created Nurbs"<<endl;


Foam::cutCellFvMesh mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    ),
    items    
);

/*
Foam::fvMesh mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
*/
Foam::Info<<"CUTCELLFVMESH Test done"<<Foam::endl;

