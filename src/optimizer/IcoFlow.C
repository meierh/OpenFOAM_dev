#include "IcoFlow.H"

Foam::IcoFlow::IcoFlow
(
    int argc,
    char *argv[]
):
Application(argc,argv),
piso(mesh),
transportProperties(IOobject("transportProperties",runTime.constant(),mesh,IOobject::MUST_READ_IF_MODIFIED,IOobject::NO_WRITE)),
nu("nu",dimViscosity,transportProperties.lookup("nu")),
p(IOobject("p",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh),
U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh),
phi(IOobject("phi",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),fvc::flux(U))
{
    createPiso();
    createFields();
    initContinuityErrs();
}

void Foam::IcoFlow::createPiso()
{}

void Foam::IcoFlow::createFields()
{
    Info<< "Reading transportProperties\n" << endl;

    transportProperties = IOdictionary
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    nu = dimensionedScalar
    (
        "nu",
        dimViscosity,
        transportProperties.lookup("nu")
    );

    Info<< "Reading field p\n" << endl;
    p = volScalarField
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    U = volVectorField
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    phi = surfaceScalarField
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    );
    
    pRefCell = 0;
    pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());
}

void Foam::IcoFlow::initContinuityErrs()
{
    cumulativeContErr = 0;
}
