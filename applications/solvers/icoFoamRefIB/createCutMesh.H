Foam::Info
    << "Create cutCellMesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

Foam::Info<<"CUTCELLFVMESH"<<Foam::endl;    

Foam::cutCellFvMesh mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    ),
    cutCellFvMesh::cutStatus::delNegMesh,
    0.1
);

Foam::Info<<"CUTCELLFVMESH Test done"<<Foam::endl;

