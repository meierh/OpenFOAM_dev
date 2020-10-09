#include "cutCellFvMesh.H"

void Foam::UnitTest_cutCellFvMesh(int argc, char *argv[],Time& runTime)
{
    Foam::Info<<"CUTCELLFVMESH"<<Foam::endl;    
    
    int testdegree = 2;
    scalarList knots;
    scalarList weights;
    List<vector> controlPoints;
    List<std::shared_ptr<Nurbs>> items(0);
    
    knots = scalarList(7);
    knots[0]=0; knots[1]=0; knots[2]=0; knots[3]=0.5;
    knots[4]=1; knots[5]=1; knots[6]=1;    
    weights = scalarList(4);
    weights[0] = 1;    weights[1] = 1;   weights[2] = 1;    weights[3] = 1; 
    controlPoints = List<vector>(4);
    controlPoints[0]=vector(1,0,0.5); controlPoints[1]=vector(1,0,0);
    controlPoints[2]=vector(3,0,0); controlPoints[3]=vector(3,0,-0.5);
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.2,4)));
    Info<<"Created Nurbs"<<endl;

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
    
    /*
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
     */
}
