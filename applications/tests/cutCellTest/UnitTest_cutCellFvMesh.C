#include "cutCellFvMesh.H"

void Foam::UnitTest_cutCellFvMesh(int argc, char *argv[],Time& runTime)
{
    Foam::Info<<"CUTCELLFVMESH"<<Foam::endl;    
    
    int testdegree = 2;
    scalarList knots;
    scalarList weights;
    List<vector> controlPoints;
    DynamicList<std::shared_ptr<Nurbs>> items(0);
    
    knots = scalarList(7);
    knots[0]=0; knots[1]=0; knots[2]=0; knots[3]=0.5;
    knots[4]=1; knots[5]=1; knots[6]=1;    
    weights = scalarList(4);
    weights[0] = 1;    weights[1] = 1;   weights[2] = 1;    weights[3] = 1; 
    controlPoints = List<vector>(4);
    
    scalar deltaX = 1;

// First row    
    controlPoints[0]=vector(1,-0.5,-0.5); controlPoints[1]=vector(1,-0.25,-0.375);
    controlPoints[2]=vector(1,-0.125,-0.25); controlPoints[3]=vector(1,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));
    
    controlPoints[0]=vector(1,-0.5,0.5); controlPoints[1]=vector(1,-0.25,0.375);
    controlPoints[2]=vector(1,-0.125,0.25); controlPoints[3]=vector(1,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));
    
    controlPoints[0]=vector(1,0.5,-0.5); controlPoints[1]=vector(1,0.25,-0.375);
    controlPoints[2]=vector(1,0.125,-0.25); controlPoints[3]=vector(1,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));
    
    controlPoints[0]=vector(1,0.5,0.5); controlPoints[1]=vector(1,0.25,0.375);
    controlPoints[2]=vector(1,0.125,0.25); controlPoints[3]=vector(1,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));

// Second row
    controlPoints[0]=vector(2,0,-1); controlPoints[1]=vector(2,0,-0.66);
    controlPoints[2]=vector(2,0,-0.33); controlPoints[3]=vector(2,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));
    
    controlPoints[0]=vector(2,-1,0); controlPoints[1]=vector(2,-0.66,0);
    controlPoints[2]=vector(2,-0.33,0); controlPoints[3]=vector(2,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));
    
    controlPoints[0]=vector(2,0,1); controlPoints[1]=vector(2,0,0.66);
    controlPoints[2]=vector(2,0,0.33); controlPoints[3]=vector(2,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));
    
    controlPoints[0]=vector(2,1,0); controlPoints[1]=vector(2,0.66,0);
    controlPoints[2]=vector(2,0.33,0); controlPoints[3]=vector(2,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));

// Third row
    controlPoints[0]=vector(4,-0.5,-0.5); controlPoints[1]=vector(3.75,-0.25,-0.375);
    controlPoints[2]=vector(3.5,-0.125,-0.25); controlPoints[3]=vector(3.25,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));
    
    controlPoints[0]=vector(4,-0.5,0.5); controlPoints[1]=vector(3.75,-0.25,0.375);
    controlPoints[2]=vector(3.5,-0.125,0.25); controlPoints[3]=vector(3.25,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));
    
    controlPoints[0]=vector(4,0.5,-0.5); controlPoints[1]=vector(3.75,0.25,-0.375);
    controlPoints[2]=vector(3.5,0.125,-0.25); controlPoints[3]=vector(3.25,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));
    
    controlPoints[0]=vector(4,0.5,0.5); controlPoints[1]=vector(3.75,0.25,0.375);
    controlPoints[2]=vector(3.5,0.125,0.25); controlPoints[3]=vector(3.25,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));
    
// Fourth part
    knots = scalarList(7);
    knots[0]=0; knots[1]=0; knots[2]=0; knots[3]=0.5;
    knots[4]=1; knots[5]=1; knots[6]=1;    
    weights = scalarList(4);
    weights[0] = 1;    weights[1] = 1;   weights[2] = 1;    weights[3] = 1; 
    controlPoints = List<vector>(4);    
    deltaX = 1;
    controlPoints[0]=vector(0,0,0); controlPoints[1]=vector(3.5,0,0);
    controlPoints[2]=vector(4,0,0); controlPoints[3]=vector(4.5,0,0);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.1,deltaX)));    
    Info<<"Created Nurbs"<<endl;

    /*
    std::unique_ptr<volScalarField> solidFraction;
     
    Foam::cutCellFvMesh nurbsMesh
    (
        Foam::IOobject
        (
            Foam::polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
            ),
        items,
        runTime,
        solidFraction
    );
    */
    
    Foam::cutCellFvMesh nurbsMesh
    (
        Foam::IOobject
        (
            Foam::polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
            ),
        items
     );
    
     Foam::Info<<"CUTCELLFVMESH Test done"<<Foam::endl;

}
