#include "icoAdjointImmersedBoundary.H"

Foam::solvers::icoAdjointImmersedBoundary::icoAdjointImmersedBoundary
(
    fvMesh& mesh
):
icoImmersedBoundary(mesh),
adj_U_
(
    IOobject
    (
        "adj_U",runTime.name(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE
    ),
    mesh
)
{
    create_AdjointVelocityForcing();
    create_AdjointTemperature();
    create_AdjointTemperatureForcing(); 
}

void Foam::solvers::icoAdjointImmersedBoundary::create_AdjointVelocityForcing()
{
    IOobject adj_fU_IOobj
    (
        "adj_fU",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if(adj_fU_IOobj.filePath("",true).empty())
        adj_fU_ = std::make_unique<volVectorField>(adj_fU_IOobj,mesh);
}

void Foam::solvers::icoAdjointImmersedBoundary::create_AdjointTemperature()
{
    IOobject adj_T_IOobj
    (
        "adj_T",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if(adj_T_IOobj.filePath("",true).empty())
    {
        adj_T_ = std::make_unique<volScalarField>(adj_T_IOobj,mesh);
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::create_AdjointTemperatureForcing()
{
    IOobject adj_fT_IOobj
    (
        "adj_fT",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if(adj_fT_IOobj.filePath("",true).empty() && useTemperature)
        adj_fT_ = std::make_unique<volScalarField>(adj_fT_IOobj,mesh);
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_preSolve()
{
    volVectorField& U(U_);    
    phi = mesh.Sf() * U;
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_momentumPredictor()
{
    volVectorField& adj_U(adj_U_);
    surfaceScalarField& phi(phi_);

    tadj_UEqn =
    (
        fvm::ddt(adj_U)
        + fvm::div(-phi, adj_U) - fvc::grad(adj_Ua) & U
      + MRF.DDt(adj_U) // coriolis force
      + momentumTransport->divDevSigma(adj_U) // stress tensor
     ==
        fvModels().source(adj_U)
    );
    fvVectorMatrix& adj_UEqn = tadj_UEqn.ref();

    adj_UEqn.relax();

    fvConstraints().constrain(adj_UEqn);

    if (pimple.momentumPredictor())
    {
        if(useAdjointVelocityForcing)
        {
           volVectorField& adj_fU = *adj_fU_;
           if(useAdjointTemperature)
               adjUEqn_res = solve(UEqn == -fvc::grad(adj_p_) + T*fvc::grad(adj_T_) + adj_fU);
           else
               adjUEqn_res = solve(UEqn == -fvc::grad(adj_p_) + adj_fU);
           interaction_adj_fU->solve();
        }
        else
        {
            if(useAdjointTemperature)
                adjUEqn_res = solve(UEqn == -fvc::grad(adj_p_) + T*fvc::grad(adj_T_));
            else
                adjUEqn_res = solve(UEqn == -fvc::grad(adj_p_));
        }
        fvConstraints().constrain(adj_fU_);
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_thermophysicalPredictor()
{
    if(useAdjointTemperature)
    {
        volScalarField& adj_T = *adj_T_;
        fvScalarMatrix adj_TEqn(fvm::ddt(adj_T)+fvm::div(-phi,adj_T)+fvm::laplacian(alpha,adj_T));
        do
        {
            if(useTemperatureForcing)
            {
                volScalarField& adj_fT = *adj_fT_;
                adjTEqn_res = solve(adj_TEqn==adj_fT);
                
                interaction_ajd_fT->solve();
            }
            else
            {
                adjTEqn_res = solve(adj_TEqn); 
            }
        }
        while(adjTEqn_res.nIterations()>0);
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_pressureCorrector()
{
     while (pimple.correct())
     {
         adj_correctPressure();
     }
     tadj_UEqn.clear();
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_correctPressure()
{
    volScalarField& p(p_);
    volVectorField& U(U_);
    surfaceScalarField& phi(phi_);
    
    volScalarField& adj_p(adj_p_);
    volVectorField& adj_U(adj_U_);

    fvVectorMatrix& adj_UEqn = tadj_UEqn.ref();

    volScalarField rA_adjU(1.0/adj_UEqn.A());
    volVectorField adj_HbyA(constrainHbyA(rA_adjU*adj_UEqn.H(), adj_U, adj_p));
    surfaceScalarField adj_phiHbyA
    (
        "adj_phiHbyA",
        fvc::flux(adj_HbyA)
        + fvc::interpolate(rA_adjU) // *fvc::ddtCorr(U, phi, Uf)
    );

    MRF.makeRelative(adj_phiHbyA);

    if (p.needReference())
    {
        fvc::makeRelative(adj_phiHbyA, adj_U);
        adjustPhi(adj_phiHbyA, adj_U, adj_p);
        fvc::makeAbsolute(adj_phiHbyA, adj_U);
    }

    tmp<volScalarField> rA_adjtU(rA_adjU);

    if (pimple.consistent())
    {
        rA_adjtU = 1.0/max(1.0/rA_adjU - UEqn.H1(), 0.1/rA_adjU);
        adj_phiHbyA += fvc::interpolate(rA_adjtU() - rA_adjU)*fvc::snGrad(adj_p)*mesh.magSf();
        adj_HbyA -= (rA_adjU - rA_adjtU())*fvc::grad(adj_p);
    }

    if (pimple.nCorrPiso() <= 1)
    {
        tadj_UEqn.clear();
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(adj_p, adj_U, adj_phiHbyA, rA_adjtU(), MRF);

    // Evaluate any volume sources
    fvScalarMatrix p_rghEqnSource(fvModels().sourceProxy(adj_p));

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        std::unique_ptr<fvScalarMatrix> adj_pEqnPtr;
        if(useVelocityForcing)
        {
            volVectorField& adj_fU = *adj_fU_;
            if(useAdjointTemperature)
            {
                adj_pEqnPtr = std::make_unique<fvScalarMatrix>
                (
                    fvm::laplacian(rA_adjtU(), adj_p)
                    ==
                    fvc::div(adj_phiHbyA)
                    - p_rghEqnSource
                    + fvc::div(rA_adjtU()*adj_fU) + fvc::div(T*fvc::grad(adj_T_)) 
                );
            }
            else
            {
                adj_pEqnPtr = std::make_unique<fvScalarMatrix>
                (
                    fvm::laplacian(rA_adjtU(), adj_p)
                    ==
                    fvc::div(adj_phiHbyA)
                    - p_rghEqnSource
                    + fvc::div(rA_adjtU()*adj_fU) 
                );
            }
        }
        else
        {            
            if(useAdjointTemperature)
            {
                adj_pEqnPtr = std::make_unique<fvScalarMatrix>
                (
                    fvm::laplacian(rA_adjtU(), adj_p)
                    ==
                    fvc::div(adj_phiHbyA)
                    - p_rghEqnSource
                    + fvc::div(T*fvc::grad(adj_T_)) 
                );
            }
            else
            {
                adj_pEqnPtr = std::make_unique<fvScalarMatrix>
                (
                    fvm::laplacian(rA_adjtU(), adj_p)
                    ==
                    fvc::div(adj_phiHbyA)
                    - p_rghEqnSource
                );
            }
        }
        fvScalarMatrix& adj_pEqn = *adj_pEqnPtr;
        
        adj_pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        adj_pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = adj_phiHbyA - adj_pEqn.flux();
        }
    }

    continuityErrors();

    // Explicitly relax pressure for momentum corrector
    adj_p.relax();

    if(useVelocityForcing)
    {
        volVectorField& fU = *fU_;
        U = HbyA - rAtU*fvc::grad(p) + rAtU*fU;
    }
    else
    {
        U = HbyA - rAtU*fvc::grad(p);
    }
    U.correctBoundaryConditions();
    fvConstraints().constrain(U);

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi, MRF);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}
