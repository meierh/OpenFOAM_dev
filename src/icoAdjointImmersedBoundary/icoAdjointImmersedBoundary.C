#include "icoAdjointImmersedBoundary.H"

Foam::solvers::icoAdjointImmersedBoundary::icoAdjointImmersedBoundary
(
    fvMesh& mesh,
    Time& time
):
icoImmersedBoundary(mesh,time),
adjPimpleCtlr(incompressibleFluid::pimple,time),
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
    
    connectSolverPerformance(adjPimpleCtlr);
    
    steadyStateAdjoint = adjPimpleCtlr.isSteady();
    
    if(!steadyStateAdjoint)
        FatalErrorInFunction<<"Unsteady adjoint not implemented"<<exit(FatalError);
    
    checkDimensions();
    
    setAdjUBC();
    
    FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
    
    Info<<"--------------------------icoImmersedBoundary--------------------------"<<Foam::endl;
    Info<<"steadyStateAdjoint:"<<steadyStateAdjoint<<Foam::endl;
    Info<<"useAdjointVelocityForcing:"<<useAdjointVelocityForcing<<Foam::endl;
    Info<<"useAdjointTemperature:"<<useAdjointTemperature<<Foam::endl;
    Info<<"useAdjointTemperatureForcing:"<<useAdjointTemperatureForcing<<Foam::endl;
    Info<<"||||||||||||||||||||||||||icoImmersedBoundary||||||||||||||||||||||||||"<<Foam::endl;
}

void Foam::solvers::icoAdjointImmersedBoundary::create_AdjointVelocityForcing()
{
    if(false && useVelocityForcing)
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
        if(!adj_fU_IOobj.filePath("",true).empty())
            adj_fU_ = std::make_unique<volVectorField>(adj_fU_IOobj,mesh);
        else
            FatalErrorInFunction<<"Missing adj_fU file"<<exit(FatalError);
        
        ITstream rodMovementStream = structureDict->lookup("rodMovement");
        token rodMovementToken;
        rodMovementStream.read(rodMovementToken);
        if(!rodMovementToken.isWord())
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodType -- must be string"<<exit(FatalError);
        word rodMovementWord = rodMovementToken.wordToken();
        if(rodMovementWord == "StaticRod")
        {
            interaction_fU = std::make_unique<StaticVelocityPressureAction>(mesh,*structure,U_,*fU_,*structureDict,refinement_);
        }
        else if(rodMovementWord == "MovedRod")
        {
            interaction_fU = std::make_unique<ForcedMovementVelocityPressureAction>(mesh,*structure,U_,*fU_,*structureDict,refinement_);
        }
        else if(rodMovementWord == "FluidStructureRod")
        {
            FatalErrorInFunction<<"Not yet implemented"<<exit(FatalError);
        }
        else
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodMovement -- valid {StaticRod,MovedRod,FluidStructureRod}"<<exit(FatalError);
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
        if(!adj_T_IOobj.filePath("",true).empty())
            adj_T_ = std::make_unique<volScalarField>(adj_T_IOobj,mesh);
        else
            FatalErrorInFunction<<"Failed to create_AdjointTemperature"<<exit(FatalError);
    }
    else
        useAdjointTemperature = false;
}

void Foam::solvers::icoAdjointImmersedBoundary::create_AdjointTemperatureForcing()
{
    if(false && useTemperatureForcing)
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
        if(!adj_fT_IOobj.filePath("",true).empty())
            adj_fT_ = std::make_unique<volScalarField>(adj_fT_IOobj,mesh);
        else
            FatalErrorInFunction<<"Missing adj_fT file"<<exit(FatalError);
    }
    else
        useAdjointTemperatureForcing = false;
}

void Foam::solvers::icoAdjointImmersedBoundary::checkDimensions()
{
    //check adj_U
    if(adj_U_.dimensions() != dimensionSet(0,1,-1,0,0,0,0))
        FatalErrorInFunction<<"Wrong dimensions of adjU, given:"<<adj_U_.dimensions()<<" necessary:"<<dimensionSet(0,1,-1,0,0,0,0)<<exit(FatalError);
    
    //check adj_p
    if(adj_p_.dimensions() != dimensionSet(0,2,-2,0,0,0,0))
        FatalErrorInFunction<<"Wrong dimensions of adjP, given:"<<adj_p_.dimensions()<<" necessary:"<<dimensionSet(0,2,-2,0,0,0,0)<<exit(FatalError);

    //check adj_T_
    if(adj_T_)
    {
        if(adj_T_->dimensions() != dimensionSet(0,2,-2,-1,0,0,0))
            FatalErrorInFunction<<"Wrong dimensions of adjT, given:"<<adj_T_->dimensions()<<" necessary:"<<dimensionSet(0,2,-2,-1,0,0,0)<<exit(FatalError);
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::connectSolverPerformance
(
    pimpleAdjIBControl& pimpleCtlr
)
{
    icoImmersedBoundary::connectSolverPerformance(pimpleCtlr);
    pimpleCtlr.setAdjVelocityPerformance(&adjUEqn_res);
    pimpleCtlr.setAdjPressurePerformance(&adjPEqn_res);
    if(useAdjointTemperature)
        pimpleCtlr.setAdjTemperaturePerformance(&adjTEqn_res);
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_preSolve
(
    pimpleAdjIBControl& adjPimpleCtlr
)
{   
    Info<<"adj_preSolve  useAdjointTemperature:"<<useAdjointTemperature<<Foam::nl;
    if(useAdjointTemperature)
    {        
        volScalarField& adj_T = *adj_T_;
        std::unique_ptr<fvScalarMatrix> adjTEqnPtr;
        if(steadyStateAdjoint)
            adjTEqnPtr = std::make_unique<fvScalarMatrix>(-fvm::div(phi,adj_T)+fvm::laplacian(alpha,adj_T));
        else
            adjTEqnPtr = std::make_unique<fvScalarMatrix>(fvm::ddt(adj_T)-fvm::div(phi,adj_T)+fvm::laplacian(alpha,adj_T));
        fvScalarMatrix& adjTEqn = *adjTEqnPtr;

        while(adjPimpleCtlr.adjTemperatureLoop())
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
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_momentumPredictor
(
    pimpleAdjIBControl& adjPimpleCtlr
)
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

    if (adjPimpleCtlr.momentumPredictor())
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

void Foam::solvers::icoAdjointImmersedBoundary::adj_pressureCorrector
(
    pimpleAdjIBControl& adjPimpleCtlr
)
{
     while (adjPimpleCtlr.adjCorrect())
     {
         adj_correctPressure(adjPimpleCtlr);
     }
     tadj_UEqn.clear();
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_correctPressure
(
    pimpleAdjIBControl& adjPimpleCtlr
)
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

    if (adjPimpleCtlr.consistent())
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
    while (adjPimpleCtlr.correctNonOrthogonal())
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

void Foam::solvers::icoAdjointImmersedBoundary::adj_postSolve
(
    pimpleAdjIBControl& adjPimpleCtlr
)
{
    for(std::pair<Parameter,scalar>& singleParameter : gradient)
    {
        if(interaction_adj_fU)
            singleParameter.second += interaction_adj_fU->computeSensitivity(singleParameter.first);
        if(interaction_ajd_fT)
            singleParameter.second += interaction_ajd_fT->computeSensitivity(singleParameter.first);
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::oneAdjSteadyTimestep
(
    pimpleAdjIBControl& adjPimpleCtlr
)
{
    Info<<"--------------------------------------- Solve Adjoint ---------------------------------------"<<Foam::nl;
    adj_preSolve(adjPimpleCtlr);
    while (adjPimpleCtlr.adjMomentumLoop())
    {
        adj_moveMesh();
        adj_motionCorrector();
        //adj_fvModels().correct();
        adj_prePredictor();
        adj_momentumPredictor(adjPimpleCtlr);
        adj_thermophysicalPredictor();
        adj_pressureCorrector(adjPimpleCtlr);
        adj_postCorrector();
    }
    adj_postSolve(adjPimpleCtlr);
}

void Foam::solvers::icoAdjointImmersedBoundary::oneAdjSteadyTimestep()
{
    oneAdjSteadyTimestep(adjPimpleCtlr);
}

void Foam::solvers::icoAdjointImmersedBoundary::SolveSteadyAdjoint()
{
    
    while (adjPimpleCtlr.run(time))
    {
        
        Info<<"-------------------------------------- Presolve Time = "<<runTime.userTimeName()<<" --------------------------------------"<<nl;
        adjPimpleCtlr.read();
        preSolve(adjPimpleCtlr);
        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(time, *this);
        time++;
        Info<< "---------------------------------------- Time = "<<runTime.userTimeName()<<" -------------------------------------------"<<nl;        
        
        icoImmersedBoundary::oneTimestep(adjPimpleCtlr);
        oneAdjSteadyTimestep(adjPimpleCtlr);
        
        write_Analysis();
        runTime.write();
            
        Info<<"ExecutionTime = "<<runTime.elapsedCpuTime()<<" s"<<"  ClockTime = "<<runTime.elapsedClockTime()<<" s"<<nl<< nl;
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::setAdjUBC()
{
    const fvBoundaryMesh& boundaryMesh = mesh.boundary();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    GeometricBoundaryField<vector,fvPatchField,volMesh>& boundary = adj_U_.boundaryFieldRef();
    
    
    //Inlet
    label inletInd = boundaryMesh.findIndex("inlet");
    if(inletInd==-1)
        FatalErrorInFunction<<"Not found inlet patch"<<exit(FatalError);
    fvPatchField<vector>& inletBoundary = boundary[inletInd];
    const fvPatch& inletPatch = inletBoundary.patch();
    const labelList& inletPatchFaceCells = inletPatch.faceCells();
    vectorField inletAdjU(inletPatch.size(), vector(0,0,0));
    if(!dJgdp)
        FatalErrorInFunction<<"Objective function derivative dJgdp missing"<<exit(FatalError);
    for(label patchFaceI=0; patchFaceI<inletPatch.size(); patchFaceI++)
    {
        label faceInd = patchFaceI+inletPatch.start();
        vector faceNormal = faces[faceInd].normal(points);
        label cellInd = inletPatchFaceCells[patchFaceI];
        vector u = U_[cellInd];
        scalar p = p_[cellInd];
        scalar T = 0;
        if(T_)
            T = (*T_)[cellInd];
        inletAdjU[patchFaceI] = -1*faceNormal*dJgdp(faceNormal,u,p,T);
    }
    inletBoundary = inletAdjU;
    
    //Wall
    for(label patchI=0; patchI<boundary.types().size(); patchI++)
    {
        if(boundary.types()[patchI]=="noSlip")
        {
            fvPatchField<vector>& wallBoundary = boundary[inletInd];
            const fvPatch& wallPatch = outletBoundary.patch();
            const labelList& wallPatchFaceCells = wallPatch.faceCells();
            vectorField wallAdjU(wallPatch.size(), vector(0,0,0));
            if(!dJgdp)
                FatalErrorInFunction<<"Objective function derivative dJgdp missing"<<exit(FatalError);
            for(label patchFaceI=0; patchFaceI<wallPatch.size(); patchFaceI++)
            {
                label faceInd = patchFaceI+wallPatch.start();
                vector faceNormal = faces[faceInd].normal(points);
                label cellInd = wallPatchFaceCells[patchFaceI];
                vector u = U_[cellInd];
                scalar p = p_[cellInd];
                scalar T = 0;
                if(T_)
                    T = (*T_)[cellInd];
                wallAdjU[patchFaceI] = -1*faceNormal*dJgdp(faceNormal,u,p,T);
            }
            outletBoundary = outletAdjU;
        }
    }
    
    //Outlet
    label outletInd = boundaryMesh.findIndex("outlet");
    if(outletInd==-1)
        FatalErrorInFunction<<"Not found outlet patch"<<exit(FatalError);
    fvPatchField<vector>& outletBoundary = boundary[inletInd];
    const fvPatch& outletPatch = outletBoundary.patch();
    const labelList& outletPatchFaceCells = outletPatch.faceCells();
    vectorField outletAdjU(outletPatch.size(), vector(0,0,0));
    if(!dJgdp)
        FatalErrorInFunction<<"Objective function derivative dJgdp missing"<<exit(FatalError);
    for(label patchFaceI=0; patchFaceI<outletPatch.size(); patchFaceI++)
    {
        label faceInd = patchFaceI+outletPatch.start();
        vector faceNormal = faces[faceInd].normal(points);
        label cellInd = outletPatchFaceCells[patchFaceI];
        vector u = U_[cellInd];
        scalar p = p_[cellInd];
        scalar T = 0;
        if(T_)
            T = (*T_)[cellInd];
        outletAdjU[patchFaceI] = -1*faceNormal*dJgdp(faceNormal,u,p,T);
    }
    outletBoundary = outletAdjU;
}

void Foam::solvers::icoAdjointImmersedBoundary::setAdjPBC()
{
    const fvBoundaryMesh& boundaryMesh = mesh.boundary();
    GeometricBoundaryField<vector,fvPatchField,volMesh>& boundary = adj_U_.boundaryFieldRef();
}

void Foam::solvers::icoAdjointImmersedBoundary::setAdjTBC()
{
    const fvBoundaryMesh& boundaryMesh = mesh.boundary();
    GeometricBoundaryField<vector,fvPatchField,volMesh>& boundary = adj_U_.boundaryFieldRef();
}

void Foam::solvers::icoAdjointImmersedBoundary::set_DJgammaDp
(
    std::function<scalar(vector n, vector u, scalar p, scalar T)> dJgdp
)
{
    this->dJgdp=dJgdp;
}

void Foam::solvers::icoAdjointImmersedBoundary::set_DJgammaDun
(
    std::function<scalar(vector n, vector u, scalar p, scalar T)> dJgdun
)
{
    this->dJgdun=dJgdun;
}

void Foam::solvers::icoAdjointImmersedBoundary::set_DJgammaDut
(
    std::function<scalar(vector n, vector u, scalar p, scalar T)> dJgdut
)
{
    this->dJgdut=dJgdut;
}

void Foam::solvers::icoAdjointImmersedBoundary::set_DJgammaDT
(
    std::function<scalar(vector n, vector u, scalar p, scalar T)> dJgdT
)
{
    this->dJgdT=dJgdT;
}
