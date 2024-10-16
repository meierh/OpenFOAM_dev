#include "icoAdjointImmersedBoundary.H"

Foam::solvers::icoAdjointImmersedBoundary::icoAdjointImmersedBoundary
(
    fvMesh& mesh,
    Time& time
):
icoImmersedBoundary(mesh),
time(time),
pimpleCtlr(incompressibleFluid::pimple),
adjPimpleCtlr(incompressibleFluid::pimple),
adj_U_
(
    IOobject
    (
        "adj_U",runTime.name(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE
    ),
    mesh
),
adj_p_
(
    IOobject
    (
        "adj_p",runTime.name(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE
    ),
    mesh
)
{
    create_AdjointVelocityForcing();
    create_AdjointTemperature();
    create_AdjointTemperatureForcing();
    
    // Set the initial time-step
    setDeltaT(time,*this);
    
    Info<<"--------------------------icoImmersedBoundary--------------------------"<<Foam::endl;
    Info<<"useAdjointVelocityForcing:"<<useAdjointVelocityForcing<<Foam::endl;
    Info<<"useAdjointTemperature:"<<useAdjointTemperature<<Foam::endl;
    Info<<"useAdjointTemperatureForcing:"<<useAdjointTemperatureForcing<<Foam::endl;
    Info<<"||||||||||||||||||||||||||icoImmersedBoundary||||||||||||||||||||||||||"<<Foam::endl;
}

void Foam::solvers::icoAdjointImmersedBoundary::create_AdjointVelocityForcing()
{
    if(useVelocityForcing)
    {
        Info<<"create_AdjointVelocityForcing"<<Foam::endl;
        useAdjointVelocityForcing = true;
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
    else
        useAdjointVelocityForcing = false;
}

void Foam::solvers::icoAdjointImmersedBoundary::create_AdjointTemperature()
{
    if(useTemperature)
    {
        Info<<"create_AdjointTemperature"<<Foam::endl;
        useAdjointTemperature = true;
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
    else
        useAdjointTemperature = false;
}

void Foam::solvers::icoAdjointImmersedBoundary::create_AdjointTemperatureForcing()
{
    if(useTemperatureForcing)
    {
        Info<<"create_AdjointTemperatureForcing"<<Foam::endl;
        useAdjointTemperatureForcing = true;        
        IOobject adj_fT_IOobj
        (
            "adj_fT",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );
        if(adj_fT_IOobj.filePath("",true).empty())
            adj_fT_ = std::make_unique<volScalarField>(adj_fT_IOobj,mesh);
    }
    else
        useAdjointTemperatureForcing = false;
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_preSolve()
{   
    if(useTemperature)
    {
        //volVectorField& U(U_);    
        //phi = mesh.Sf() * U;
        
        volScalarField& adj_T = *adj_T_;
        std::unique_ptr<fvScalarMatrix> adjTEqnPtr;
        if(steadyStateAdjoint)
            adjTEqnPtr = std::make_unique<fvScalarMatrix>(-fvm::div(phi,adj_T)+fvm::laplacian(alpha,adj_T));
        else
            adjTEqnPtr = std::make_unique<fvScalarMatrix>(fvm::ddt(adj_T)-fvm::div(phi,adj_T)+fvm::laplacian(alpha,adj_T));
        fvScalarMatrix& adjTEqn = *adjTEqnPtr;
        
        do
        {
            if(useTemperatureForcing)
            {
                volScalarField& adj_fT = *adj_fT_;
                adjTEqn_res = solve(adjTEqn==adj_fT);
                interaction_ajd_fT->solve();
            }
            else
            {
                adjTEqn_res = solve(adjTEqn); 
            }
        }
        while(adjTEqn_res.nIterations()>0);
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_momentumPredictor()
{
    volVectorField& U(U_);
    volVectorField& adj_U(adj_U_);
    surfaceScalarField& phi(phi_);

    volVectorField adjointTransposeConvection((fvc::grad(adj_U) & U));
    
    if(steadyStateAdjoint)
    {
        if(useTemperature && useAdjointTemperature)
        {
            volScalarField& adj_T = *adj_T_;
            volScalarField& T = *T_;
            if(useAdjointVelocityForcing)
            {
                volVectorField& adj_fU = *adj_fU_;
                tadj_UEqn = 
                (
                    fvm::div(-phi,adj_U) - adjointTransposeConvection + fvm::laplacian(nu,adj_U)
                    - T*fvc::grad(adj_T)
                    - adj_fU
                );
            }
            else
            {
                tadj_UEqn = 
                (
                    fvm::div(-phi,adj_U) - adjointTransposeConvection + fvm::laplacian(nu,adj_U)
                    - T*fvc::grad(adj_T)
                );
            }
        }
        else
        {
            if(useAdjointVelocityForcing)
            {
                volVectorField& adj_fU = *adj_fU_;
                tadj_UEqn = 
                (
                    fvm::div(-phi,adj_U) - adjointTransposeConvection + fvm::laplacian(nu,adj_U)
                    - adj_fU
                );
            }
            else
            {
                tadj_UEqn = 
                (
                    fvm::div(-phi,adj_U) - adjointTransposeConvection + fvm::laplacian(nu,adj_U)
                );
            }
        }
    }
    else
    {
        if(useTemperature && useAdjointTemperature)
        {
            volScalarField& adj_T = *adj_T_;
            volScalarField& T = *T_;
            if(useAdjointVelocityForcing)
            {
                volVectorField& adj_fU = *adj_fU_;
                tadj_UEqn = 
                (
                    fvm::ddt(adj_U)
                    + fvm::div(-phi,adj_U) - adjointTransposeConvection + fvm::laplacian(nu,adj_U)
                    - T*fvc::grad(adj_T)
                    - adj_fU
                );
            }
            else
            {
                tadj_UEqn = 
                (
                    fvm::ddt(adj_U)
                    + fvm::div(-phi,adj_U) - adjointTransposeConvection + fvm::laplacian(nu,adj_U)
                    - T*fvc::grad(adj_T)
                );
            }
        }
        else
        {
            if(useAdjointVelocityForcing)
            {
                volVectorField& adj_fU = *adj_fU_;
                tadj_UEqn = 
                (
                    fvm::ddt(adj_U)
                    + fvm::div(-phi,adj_U) - adjointTransposeConvection + fvm::laplacian(nu,adj_U)
                    - adj_fU
                );
            }
            else
            {
                tadj_UEqn = 
                (
                    fvm::ddt(adj_U)
                    + fvm::div(-phi,adj_U) - adjointTransposeConvection + fvm::laplacian(nu,adj_U)
                );
            }
        }
    }

    fvVectorMatrix& adj_UEqn = tadj_UEqn.ref();

    adj_UEqn.relax();
    fvConstraints().constrain(adj_UEqn);

    if (pimple.momentumPredictor())
    {
        if(useAdjointVelocityForcing)
        {
            adjUEqn_res = solve(adj_UEqn == -fvc::grad(adj_p_));
           interaction_adj_fU->solve();
        }
        else
        {
            adjUEqn_res = solve(adj_UEqn == -fvc::grad(adj_p_));
        }
        fvConstraints().constrain(adj_U_);
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_thermophysicalPredictor()
{
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
    //volScalarField& p(p_);
    //volVectorField& U(U_);
    surfaceScalarField& phi(phi_);
    
    volScalarField& adj_p(adj_p_);
    volVectorField& adj_U(adj_U_);

    fvVectorMatrix& adj_UEqn = tadj_UEqn.ref();

    volScalarField rA_adjU(1.0/adj_UEqn.A());
    volVectorField adj_HbyA(constrainHbyA(rA_adjU*adj_UEqn.H(),adj_U,adj_p));
    surfaceScalarField adj_phiHbyA
    (
        "adj_phiHbyA",
        fvc::flux(adj_HbyA)
        + fvc::interpolate(rA_adjU)*fvc::ddtCorr(adj_U,phi)
    );

    //MRF.makeRelative(adj_phiHbyA);
    /*
    if (p.needReference())
    {
        fvc::makeRelative(adj_phiHbyA, adj_U);
        adjustPhi(adj_phiHbyA, adj_U, adj_p);
        fvc::makeAbsolute(adj_phiHbyA, adj_U);
    }
    */

    tmp<volScalarField> rA_adjtU(rA_adjU);

    if (pimple.consistent())
    {
        rA_adjtU = 1.0/max(1.0/rA_adjU - adj_UEqn.H1(), 0.1/rA_adjU);
        adj_phiHbyA += fvc::interpolate(rA_adjtU() - rA_adjU)*fvc::snGrad(adj_p)*mesh.magSf();
        adj_HbyA -= (rA_adjU - rA_adjtU())*fvc::grad(adj_p);
    }
    if (pimple.nCorrPiso() <= 1)
    {
        tadj_UEqn.clear();
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(adj_p, adj_U, adj_phiHbyA, rA_adjtU()/*,MRF*/);

    // Evaluate any volume sources
    //fvScalarMatrix p_rghEqnSource(fvModels().sourceProxy(adj_p));

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix adj_pEqn
        (
            fvm::laplacian(rA_adjtU(), adj_p)
            ==
            fvc::div(adj_phiHbyA)
        );
        
        adj_pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        adj_pEqn.solve();

        /*
        if (pimple.finalNonOrthogonalIter())
        {
            phi = adj_phiHbyA - adj_pEqn.flux();
        }
        */
    }

    //continuityErrors();

    // Explicitly relax pressure for momentum corrector
    adj_p.relax();
    adj_U = adj_HbyA - rA_adjtU*fvc::grad(adj_p);
    adj_U.correctBoundaryConditions();
    fvConstraints().constrain(adj_U);

    // Correct Uf if the mesh is moving
    //fvc::correctUf(Uf, U, phi, MRF);

    // Make the fluxes relative to the mesh motion
    //fvc::makeRelative(phi, U);
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_postSolve()
{
    for(std::pair<Parameter,scalar>& singleParameter : gradient)
    {
        if(interaction_adj_fU)
            singleParameter.second += interaction_adj_fU->computeSensitivity(singleParameter.first);
        if(interaction_ajd_fT)
            singleParameter.second += interaction_ajd_fT->computeSensitivity(singleParameter.first);
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::solvePrimal()
{
    while (pimpleCtlr.run(time))
    {
        // Update PIMPLE outer-loop parameters if changed
        pimpleCtlr.read();
        
        Info<< " Presolve Time = " << runTime.userTimeName() << nl << endl;
        
        preSolve();

        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(time, *this);

        time++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        preMove();
        
        // PIMPLE corrector loop
        while (pimpleCtlr.loop())
        {
            moveMesh();
            motionCorrector();
            fvModels().correct();
            prePredictor();
            momentumPredictor();
            thermophysicalPredictor();
            pressureCorrector();
            postCorrector();
        }

        postSolve();

        runTime.write();

        Foam::Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << Foam::nl << endl;
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::solveAdjoint()
{
    while (adjPimpleCtlr.run(time))
    {
        // Update PIMPLE outer-loop parameters if changed
        adjPimpleCtlr.read();

        adj_preSolve();

        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(time, *this);

        time++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // PIMPLE corrector loop
        while (adjPimpleCtlr.loop())
        {
            adj_moveMesh();
            adj_motionCorrector();
            //adj_fvModels().correct();
            adj_prePredictor();
            adj_momentumPredictor();
            adj_thermophysicalPredictor();
            adj_pressureCorrector();
            adj_postCorrector();
        }

        adj_postSolve();

        runTime.write();

        Foam::Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << Foam::nl << endl;
    }
}
