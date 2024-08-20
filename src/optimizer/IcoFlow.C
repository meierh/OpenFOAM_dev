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

void Foam::IcoFlow::run()
{
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        #include "CourantNo.H"

        // Momentum predictor
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
            + fvm::div(phi, U)
            - fvm::laplacian(nu, U)
        );
        
        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }
        
        if(velocityImmersedBoundary)
        {
            velocityImmersedBoundary.interpolateFluidVelocityToMarkers();
            velocityImmersedBoundary.computeCouplingForceOnMarkers();
            velocityImmersedBoundary.computeRodForceMoment();
            velocityImmersedBoundary.interpolateFluidForceField();
            Info<<"sum Force:"<<velocityImmersedBoundary.sumForces()<<Foam::endl;
        
            U = U + runTime.deltaT()*Uf;
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
                + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve();

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;
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
