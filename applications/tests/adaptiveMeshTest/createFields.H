Foam::Info<< "Reading field refine\n" << Foam::endl;
Foam::volScalarField refine
(
    Foam::IOobject
    (
        "refine",
        runTime.name(),
        mesh,
        Foam::IOobject::MUST_READ,
        Foam::IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::volScalarField testField
(
    "test",refine
);
