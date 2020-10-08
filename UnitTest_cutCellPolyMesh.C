#include "cutCellPolyMesh.H"

void Foam::UnitTest_cutCellPolyMesh(int argc, char *argv[],Time& runTime)
{
    Foam::Info<<"CUTCELLPOLYMESH"<<Foam::endl;    
    
    int testdegree = 2;
    scalarList knots;
    scalarList weights;
    List<vector> controlPoints;
    List<std::shared_ptr<Nurbs>> items(0);
    
    knots = scalarList(6);
    knots[0]=0; knots[1]=0; knots[2]=0; knots[3]=1; knots[4]=1; knots[5]=1;    
    weights = scalarList(2);
    weights[0] = 1;    weights[1] = 1;    
    controlPoints = List<vector>(2);
    controlPoints[0]=vector(-1,-1,0); controlPoints[1]=vector(1,1,1);    
    items.append(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.3,4)));

    std::unique_ptr<volScalarField> solidFraction;
     
    /*
    Foam::cutCellPolyMesh nurbsMesh
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
    Foam::cutCellPolyMesh nurbsMesh
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
    
}
