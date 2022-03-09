#include "cutCellFvMesh.H"
#include <cmath>

Foam::List<Foam::List<Foam::vector>> simpleCurveMovementFunc(Time& runTime)
{
    scalar timeOutputValue = runTime.timeOutputValue();
    Foam::List<Foam::List<Foam::vector>> movingControlPoints(1,List<Foam::vector>(4));
    movingControlPoints[0][0] = Foam::vector(1,-0.6,0);
    movingControlPoints[0][1] = Foam::vector(1,0,0);
    movingControlPoints[0][2] = Foam::vector(3,0,0);
    movingControlPoints[0][3] = Foam::vector(3,0.6,0);
    
    Foam::List<Foam::vector> movementVector(4);
    movementVector[0] = Foam::vector(0,0,0);
    movementVector[1] = Foam::vector(0,0.1,0);
    movementVector[2] = Foam::vector(0,-0.1,0);
    movementVector[3] = Foam::vector(0,0,0);
    
    Foam::List<Foam::List<Foam::vector>> resultingCPs(1,List<Foam::vector>(4));
    resultingCPs[0] = movingControlPoints[0] +  movementVector * std::sin(std::atan(1)*timeOutputValue);

    Info<<"resultingCPs:"<<resultingCPs<<endl;
    return resultingCPs;
}

void Foam::UnitTest_cutCellFvMesh(int argc, char *argv[],Time& runTime)
{
    Foam::Info<<"CUTCELLFVMESH"<<Foam::endl;
    
    Foam::cutCellFvMesh nurbsMesh
    (
        Foam::IOobject
        (
            Foam::polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    nurbsMesh.write();

    
    /*
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        nurbsMesh.moveNurbsCurves(simpleCurveMovementFunc(runTime));
        nurbsMesh.moveTheMesh();
        Info<<"nurbsMesh.update()"<<endl;
        nurbsMesh.update();        
        Info<<"runTime.write()"<<endl;
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    */
    //nurbsMesh.write();

    Info<< "End\n" << endl;
    
    Foam::Info<<"CUTCELLFVMESH Test done"<<Foam::endl;
}
