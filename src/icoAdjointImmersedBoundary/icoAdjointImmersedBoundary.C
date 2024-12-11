#include "icoAdjointImmersedBoundary.H"

Foam::solvers::icoAdjointImmersedBoundary::icoAdjointImmersedBoundary
(
    fvMesh& mesh,
    Time& time,
    std::vector<Parameter> parameters,
    objectiveFunction obj
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
/*
adj_phi_
(
    IOobject
    (
        "adj_phi",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(adj_U_) & mesh.Sf()
),
*/
adj_p_
(
    IOobject
    (
        "adj_p",runTime.name(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE
    ),
    mesh
)
{
    setupAdjoint();
    create_AdjointVelocityForcing();
    create_AdjointTemperature();
    create_AdjointTemperatureForcing();
    connectSolverPerformance(adjPimpleCtlr);
    steadyStateAdjoint = adjPimpleCtlr.isSteady();
            
    if(!steadyStateAdjoint)
        FatalErrorInFunction<<"Unsteady adjoint not implemented"<<exit(FatalError);
    
    checkDimensions();
            
    create_Objective(obj);
    setAdjUBC(obj.dJdp_Inlet,obj.dJdp_Wall,obj.dJdu_uOutlet);
    setAdjPBC(obj.dJdu_pOutlet);
    setAdjTBC(obj.dJdT_Outlet);
    J = obj.J;
    for(const Parameter& onePara : parameters)
        gradient.push_back({onePara,0});
    
    Info<<"--------------------------icoAdjointImmersedBoundary--------------------------"<<Foam::endl;
    Info<<"steadyStateAdjoint:"<<steadyStateAdjoint<<Foam::endl;
    Info<<"useAdjointVelocityForcing:"<<useAdjointVelocityForcing<<Foam::endl;
    Info<<"useAdjointTemperature:"<<useAdjointTemperature<<Foam::endl;
    Info<<"useAdjointTemperatureForcing:"<<useAdjointTemperatureForcing<<Foam::endl;
    Info<<"||||||||||||||||||||||||||icoAdjointImmersedBoundary||||||||||||||||||||||||||"<<Foam::endl;
}

Foam::solvers::icoAdjointImmersedBoundary::~icoAdjointImmersedBoundary()
{
}

void Foam::solvers::icoAdjointImmersedBoundary::setupAdjoint()
{
    IOobject fvSolutionIO("fvSolution","system",runTime,IOobject::MUST_READ,IOobject::NO_WRITE);
    if(!fvSolutionIO.filePath("",true).empty())
    {
        IOdictionary fvSolutionDict(fvSolutionIO);
        dictionary& adj_pimpleDict = fvSolutionDict.subDict("adj_PIMPLE");

        ITstream virtualMomentumTimestepStream = adj_pimpleDict.lookup("virtualMomentumTimestep");
        token virtualMomentumTimestepToken;
        virtualMomentumTimestepStream.read(virtualMomentumTimestepToken);
        if(!virtualMomentumTimestepToken.isScalar())
            FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/virtualMomentumTimestep -- must be scalar"<<exit(FatalError);
        virtualAdjMomentumTimestep = virtualMomentumTimestepToken.scalarToken();

        ITstream virtualTemperatureTimestepStream = adj_pimpleDict.lookup("virtualTemperatureTimestep");
        token virtualTemperatureTimestepToken;
        virtualTemperatureTimestepStream.read(virtualTemperatureTimestepToken);
        if(!virtualTemperatureTimestepToken.isScalar())
            FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/virtualMomentumTimestep -- must be scalar"<<exit(FatalError);
        virtualAdjTemperatureTimestep = virtualTemperatureTimestepToken.scalarToken();
    }
    else
        FatalErrorInFunction<<"Missing file in system/fvSolution"<<exit(FatalError);
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
        if(!adj_fU_IOobj.filePath("",true).empty())
            adj_fU_ = std::make_unique<volVectorField>(adj_fU_IOobj,mesh);
        else
        {
            Pout<<"adj_fU_IOobj.path():"<<adj_fU_IOobj.path(true)<<Foam::nl;
            FatalErrorInFunction<<"Missing adj_fU file"<<exit(FatalError);
        }
        
        ITstream rodMovementStream = structureDict->lookup("rodMovement");
        token rodMovementToken;
        rodMovementStream.read(rodMovementToken);
        if(!rodMovementToken.isWord())
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodType -- must be string"<<exit(FatalError);
        word rodMovementWord = rodMovementToken.wordToken();
        if(rodMovementWord == "StaticRod")
        {
            interaction_adj_fU = std::make_unique<SensitivityVelocityPressureForceInteraction>(mesh,*structure,*interaction_fU,adj_U_,*adj_fU_,*structureDict);
        }
        else if(rodMovementWord == "MovedRod")
        {
            FatalErrorInFunction<<"Not yet implemented"<<exit(FatalError);
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
        if(!adj_fT_IOobj.filePath("",true).empty())
            adj_fT_ = std::make_unique<volScalarField>(adj_fT_IOobj,mesh);
        else
            FatalErrorInFunction<<"Missing adj_fT file"<<exit(FatalError);
        
        ITstream rodHeatingStream = structureDict->lookup("rodHeating");
        token rodHeatingToken;
        rodHeatingStream.read(rodHeatingToken);
        if(!rodHeatingToken.isWord())
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodHeating -- must be word"<<exit(FatalError);
        word rodHeatingWord = rodHeatingToken.wordToken();
        if(rodHeatingWord == "FixedTemperature")
        {
            interaction_adj_fT = std::make_unique<SensitivityTemperatureInteraction>(mesh,*structure,*interaction_fT,*adj_T_,*adj_fT_,*structureDict);
        }
        else if(rodHeatingWord == "VariableTemperature")
        {
            FatalErrorInFunction<<"Not yet implemented"<<exit(FatalError);
        }
        else if(rodHeatingWord == "TemperatureInteractionRod")
        {
            FatalErrorInFunction<<"Not yet implemented"<<exit(FatalError);
        }
        else
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodMovement -- valid {FixedTemperature,VariableTemperature,TemperatureInteractionRod}"<<exit(FatalError);
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

void Foam::solvers::icoAdjointImmersedBoundary::create_Objective
(
    Foam::solvers::objectiveFunction& obj
)
{
    IOobject optimizationDictIO("optimizationDict","system",runTime,IOobject::MUST_READ,IOobject::NO_WRITE);
    if(!optimizationDictIO.filePath("",true).empty())
    {
        IOdictionary optimizationDict(optimizationDictIO);
        if(optimizationDict.found("objectiveFunction"))
        {
            ITstream objectiveFunctionStream = optimizationDict.lookup("objectiveFunction");
            token objectiveFunctionToken;
            objectiveFunctionStream.read(objectiveFunctionToken);
            if(!objectiveFunctionToken.isWord())
                FatalErrorInFunction<<"Invalid entry in system/optimizationDict/objectiveFunction -- must be word"<<exit(FatalError);
            word objectiveFunctionWord = objectiveFunctionToken.wordToken();
            if(objectiveFunctionWord=="totalPressureLoss")
                obj = createTotalPressureLoss();
            else if(objectiveFunctionWord=="totalPressureLossEnergy")
                FatalErrorInFunction<<"Not implemented"<<exit(FatalError);
            else if(objectiveFunctionWord=="heatGain")
                FatalErrorInFunction<<"Not implemented"<<exit(FatalError);
            else if(objectiveFunctionWord=="parameter")
            {
                if(obj.empty)
                    FatalErrorInFunction<<"Given empty objective function"<<exit(FatalError);
            }
            else
                FatalErrorInFunction<<"Invalid word in system/optimizationDict/objectiveFunction -- must be {totalPressureLoss,totalPressureLossEnergy,heatGain,parameter}"<<exit(FatalError);
        }
    }
    else
        FatalErrorInFunction<<"Missing file in system/optimizationDict"<<exit(FatalError);
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

void Foam::solvers::icoAdjointImmersedBoundary::initializeInteractions()
{
    if(useAdjointVelocityForcing)
    {
        interaction_adj_fU->solve(virtualAdjMomentumTimestep);
    }
    
    if(useAdjointTemperature)
    {        
        if(useAdjointTemperatureForcing)
        {
            interaction_adj_fT->solve(virtualAdjTemperatureTimestep);
        }
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_preSolve
(
    pimpleAdjIBControl& adjPimpleCtlr
)
{
    initializeInteractions();
    
    Info<<"J:"<<J(*this)<<Foam::endl;
    
    Info<<"adj_preSolve  useAdjointTemperature:"<<useAdjointTemperature<<Foam::nl;
    if(useAdjointTemperature)
    {        
        volScalarField& adj_T = *adj_T_;
        std::unique_ptr<fvScalarMatrix> adjTEqnPtr;

        auto generateTemperatureEqn = 
        [&]()
        {
            if(steadyStateAdjoint)
                adjTEqnPtr = std::make_unique<fvScalarMatrix>(fvm::div(-phi,adj_T)-fvm::laplacian(alpha,adj_T));
            else
                adjTEqnPtr = std::make_unique<fvScalarMatrix>(fvm::ddt(adj_T)+fvm::div(-phi,adj_T)-fvm::laplacian(alpha,adj_T));
            
            fvScalarMatrix& adjTEqn = *adjTEqnPtr;
            if(useAdjointTemperatureForcing)
            {
                volScalarField& adj_fT = *adj_fT_;
                adjTEqn += (-adj_fT);
            }
        };
        
        do
        {
            generateTemperatureEqn();
            fvScalarMatrix& adjTEqn = *adjTEqnPtr;            
            adjTEqn_res = solve(adjTEqn);
            
            if(useAdjointTemperatureForcing)
            {
                interaction_adj_fT->solve(virtualAdjTemperatureTimestep);

                adjTEqn_res = solve(adjTEqn);
            }
        }
        while(adjPimpleCtlr.adjTemperatureLoop());
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_momentumPredictor
(
    pimpleAdjIBControl& adjPimpleCtlr
)
{
    const volVectorField& U(U_);
    volVectorField& adj_U(adj_U_);
    const surfaceScalarField phi = linearInterpolate(U) & mesh.Sf();
    
    if(!steadyStateAdjoint)
        FatalErrorInFunction<<"SteadyStateAdjoint false"<<exit(FatalError);

    auto generateMomentumEqn = 
    [&]()
    {
        volVectorField adjTranspConv((fvc::grad(adj_U) & U));
        const labelList& inletCells = mesh.boundary()["inlet"].faceCells();
        for(label cellInd : inletCells)
            adjTranspConv[cellInd] = Foam::zero();

        if(useAdjointTemperature)
        {
            if(!T_)
                FatalErrorInFunction<<"T_ not a field"<<exit(FatalError);
            if(!adj_T_)
                FatalErrorInFunction<<"adj_T_ not a field"<<exit(FatalError);
            volScalarField& adj_T = *adj_T_;
            volScalarField& T = *T_;
            tadj_UEqn = fvm::div(-phi,adj_U) - adjTranspConv - fvm::laplacian(nu,adj_U) - T*fvc::grad(adj_T);
        }
        else
        {
            tadj_UEqn = fvm::div(-phi,adj_U) - adjTranspConv - fvm::laplacian(nu,adj_U);
        }
        
        if(useAdjointVelocityForcing)
        {
            volVectorField& adj_fU = *adj_fU_;
            fvVectorMatrix& adj_UEqn = tadj_UEqn.ref();
            adj_UEqn += (-adj_fU);
        }        
    };

    generateMomentumEqn();
    fvVectorMatrix& adj_UEqn = tadj_UEqn.ref();
    adj_UEqn.relax();
    fvConstraints().constrain(adj_UEqn);

    if (adjPimpleCtlr.adjMomentumPredictor())
    {
        adjUEqn_res = solve(adj_UEqn == -fvc::grad(adj_p_));
        if(useAdjointVelocityForcing)
        {
            interaction_adj_fU->solve(virtualAdjMomentumTimestep);
            //adjUEqn_res = solve(adj_UEqn == -fvc::grad(adj_p_));
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
    DynamicList<SolverPerformance<scalar>> pEqn_list;
    
    //volScalarField& p(p_);
    //volVectorField& U(U_);
    const surfaceScalarField& phi(phi_);
    
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
    while (adjPimpleCtlr.adjCorrectNonOrthogonal())
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

        SolverPerformance<scalar> pEqn = adj_pEqn.solve();
        pEqn_list.append(pEqn);
        
        if (pimple.finalNonOrthogonalIter())
        {
            //adj_phi_ = adj_phiHbyA - adj_pEqn.flux();
        }
        adj_p_.storePrevIter();
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
    if(pEqn_list.size()<1)
    {
        FatalErrorInFunction<<"Must be more than one iteration"<<exit(FatalError);
    }
    adjPEqn_res = pEqn_list[0];
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_postSolve
(
    pimpleAdjIBControl& adjPimpleCtlr
)
{
    Info<<"gradient.size():"<<gradient.size()<<Foam::nl;
    
    if(interaction_fU)
        interaction_fU->recomputeMarkerValues();
    if(interaction_fT)
        interaction_fT->recomputeMarkerValues();
    if(interaction_adj_fU)
        interaction_adj_fU->recomputeMarkerValues();
    if(interaction_adj_fT)
        interaction_adj_fT->recomputeMarkerValues();
    
    Info<<Foam::nl<<"computeSensitivity"<<Foam::nl;

    for(std::pair<Parameter,scalar>& singleParameter : gradient)
    {
        Info<<"----------------------------------------------"<<Foam::endl;
        singleParameter.second = 0;
        if(interaction_adj_fU)
            singleParameter.second += interaction_adj_fU->computeSensitivity(singleParameter.first);
        if(interaction_adj_fT)
            singleParameter.second += interaction_adj_fT->computeSensitivity(singleParameter.first);
        Info<<singleParameter.first.to_string()<<":"<<structure->getParameterValue(singleParameter.first)<<" -- "<<singleParameter.second<<Foam::nl;
        Info<<"Completed one parameter"<<Foam::endl;
        Info<<"||||||||||||||||||||||||||||||||||||||||||||||"<<Foam::endl;
    }
    Info<<Foam::nl<<"computeSensitivity done"<<Foam::nl;
}

void Foam::solvers::icoAdjointImmersedBoundary::oneAdjSteadyTimestep
(
    pimpleAdjIBControl& adjPimpleCtlr
)
{
    Info<<"--------------------------------------- Solve Adjoint ---------------------------------------"<<Foam::nl;
    Info<<"adj_U: "; printAvg(adj_U_); printMinMax(adj_U_);
    Info<<"adj_p: "; printAvg(adj_p_); printMinMax(adj_p_);
    adj_preSolve(adjPimpleCtlr);
    /*
    do
    {
        adj_moveMesh();
        adj_motionCorrector();
        //adj_fvModels().correct();
        adj_prePredictor();
        adj_momentumPredictor(adjPimpleCtlr);

        adj_thermophysicalPredictor();
        adj_pressureCorrector(adjPimpleCtlr);
        
        Info<<"adj_U: "; printAvg(adj_U_); printMinMax(adj_U_);
        Info<<"adj_p: "; printAvg(adj_p_); printMinMax(adj_p_);    
        adj_postCorrector();

        Info<<"ExecutionTime = "<<runTime.elapsedCpuTime()<<" s"<<"  ClockTime = "<<runTime.elapsedClockTime()<<" s"<<nl<< nl;
    }
    while (adjPimpleCtlr.adjMomentumLoop());
    */
    Info<<"--------------------------------------- Solved Adjoint ---------------------------------------"<<Foam::nl;
    Info<<"adj_U: "; printAvg(adj_U_); printMinMax(adj_U_);
    Info<<"adj_p: "; printAvg(adj_p_); printMinMax(adj_p_);
    runTime.write();
    adj_postSolve(adjPimpleCtlr);
}

void Foam::solvers::icoAdjointImmersedBoundary::oneAdjSteadyTimestepBase
(
    pimpleAdjIBControl& adjPimpleCtlr
)
{
    Info<<"--------------------------------------- Solve Adjoint ---------------------------------------"<<Foam::nl;
    adj_U_ = Foam::zero();
    //adj_phi_ = Foam::zero();
    adj_p_ = Foam::zero();
    
    label paRefCell = 0;
    scalar paRefValue = 0.0;
    //scalar cumulativeAdjointContErr = 0;
    
    Info<<"|||||||||||||||||||||||||||||||||||||||||||||"<<Foam::nl;
    Info<<"U: "; printAvg(U_); printMinMax(U_);
    Info<<"phi: "; printAvg(phi_); printMinMax(phi_);    
    Info<<"adj_U: "; printAvg(adj_U_); printMinMax(adj_U_);
    Info<<"adj_p: "; printAvg(adj_p_); printMinMax(adj_p_);
    Info<<"|||||||||||||||||||||||||||||||||||||||||||||"<<Foam::nl;
    
    const surfaceScalarField& phi(phi_);
    volScalarField& adj_p(adj_p_);
    volVectorField& adj_U(adj_U_);
    
    while (adjPimpleCtlr.adjMomentumLoop())
    {
        {
            volVectorField adjointTransposeConvection((fvc::grad(adj_U) & U));
            const labelList& inletCells = mesh.boundary()["inlet"].faceCells();
            forAll(inletCells,i)
                adjointTransposeConvection[inletCells[i]] = Foam::zero();

            tmp<fvVectorMatrix> tUaEqn
            (
                fvm::div(-phi, adj_U)
              - adjointTransposeConvection
              - fvm::laplacian(nu,adj_U)
              /*
              + turbulence->divDevSigma(Ua)
              + fvm::Sp(alpha, adj_U)
             ==
                fvModels.source
                */
            );
            fvVectorMatrix& UaEqn = tUaEqn.ref();

            UaEqn.relax();

            fvConstraints().constrain(UaEqn);

            solve(UaEqn == -fvc::grad(adj_p));

            fvConstraints().constrain(adj_U);

            volScalarField rAUa(1.0/UaEqn.A());
            volVectorField HbyAa("HbyAa", adj_U);
            HbyAa = rAUa*UaEqn.H();
            tUaEqn.clear();
            surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
            adjustPhi(phiHbyAa, adj_U, adj_p);

            // Non-orthogonal pressure corrector loop
            while (adjPimpleCtlr.correctNonOrthogonal())
            {
                fvScalarMatrix paEqn
                (
                    fvm::laplacian(rAUa, adj_p) == fvc::div(phiHbyAa)
                );

                paEqn.setReference(paRefCell, paRefValue);
                paEqn.solve();

                if (adjPimpleCtlr.finalNonOrthogonalIter())
                {
                    //adj_phi_ = phiHbyAa - paEqn.flux();
                }
            }
            
            /*
            {
                scalar sumLocalContErr = runTime.deltaTValue()*
                    mag(fvc::div(adj_phi_))().weightedAverage(mesh.V()).value();

                scalar globalContErr = runTime.deltaTValue()*
                    fvc::div(adj_phi_)().weightedAverage(mesh.V()).value();
                cumulativeAdjointContErr += globalContErr;

                Info<< "Adjoint continuity errors : sum local = " << sumLocalContErr
                    << ", global = " << globalContErr
                    << ", cumulative = " << cumulativeAdjointContErr
                    << endl;
            }
            */
    
            // Explicitly relax pressure for adjoint momentum corrector
            adj_p.relax();

            // Adjoint momentum corrector
            adj_U = HbyAa - rAUa*fvc::grad(adj_p);
            adj_U.correctBoundaryConditions();
            fvConstraints().constrain(adj_U);
        }
    Info<<"|||||||||||||||||||||||||||||||||||||||||||||"<<Foam::nl;
    Info<<"|| U: "; printAvg(U_); printMinMax(U_);
    Info<<"|| phi: "; printAvg(phi_); printMinMax(phi_);    
    Info<<"|| adj_U: "; printAvg(adj_U_); printMinMax(adj_U_);
    Info<<"|| adj_p: "; printAvg(adj_p_); printMinMax(adj_p_);
    Info<<"|||||||||||||||||||||||||||||||||||||||||||||"<<Foam::nl;
        
        //time += 1;
        //runTime.write();
        
    }
    adj_postSolve(adjPimpleCtlr);
}

void Foam::solvers::icoAdjointImmersedBoundary::oneAdjSteadyTimestep()
{
    oneAdjSteadyTimestep(adjPimpleCtlr);
}

void Foam::solvers::icoAdjointImmersedBoundary::SolvePrimal
(
    std::function<void(bool,const Time&,Time&)> writeData
)
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

        write_Analysis();
        writeData(false,runTime,time);
        Info<<"runTime:"<<runTime.toc()<<Foam::nl;
        Info<<"ExecutionTime = "<<runTime.elapsedCpuTime()<<" s"<<"  ClockTime = "<<runTime.elapsedClockTime()<<" s"<<nl<< nl;
    }
    writeData(true,runTime,time);
    Info<<"Primal final time = "<<runTime.elapsedCpuTime()<<" s"<<"  ClockTime = "<<runTime.elapsedClockTime()<<" s"<<nl<< nl;
}

void Foam::solvers::icoAdjointImmersedBoundary::SolveSteadyAdjoint()
{
    SolvePrimal();
    oneAdjSteadyTimestep(adjPimpleCtlr);
}

void Foam::solvers::icoAdjointImmersedBoundary::SolveFDGradient()
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

        write_Analysis();
        runTime.write();
        Info<<"runTime:"<<runTime.toc()<<Foam::nl;
        Info<<"ExecutionTime = "<<runTime.elapsedCpuTime()<<" s"<<"  ClockTime = "<<runTime.elapsedClockTime()<<" s"<<nl<< nl;
    }
    Info<<"Primal final time = "<<runTime.elapsedCpuTime()<<" s"<<"  ClockTime = "<<runTime.elapsedClockTime()<<" s"<<nl<< nl;

    oneAdjSteadyTimestep(adjPimpleCtlr);
    //oneAdjSteadyTimestepBase(adjPimpleCtlr);
    adj_U_.write();
    adj_p_.write();
    
    Info<<"SolveSteadyAdjoint done"<<nl<< nl;
}

void Foam::solvers::icoAdjointImmersedBoundary::setAdjUBC
(
    std::function<Field<scalar>(const icoAdjointVelocityInletBC&)> dJdp_Inlet,
    std::function<Field<scalar>(const icoAdjointVelocityWallBC&)> dJdp_Wall,
    std::function<Field<vector>(const icoAdjointVelocityOutletBC&)> dJdu_uOutlet
)
{
    DynamicList<label> inlet;
    DynamicList<label> wall;
    DynamicList<label> outlet;
    const fvBoundaryMesh& boundary = mesh.boundary();
    GeometricBoundaryField<vector,fvPatchField,volMesh>& boundaryField = adj_U_.boundaryFieldRef();
    for(label patchI=0; patchI<boundary.size(); patchI++)
    {
        fvPatch const& patch = boundary[patchI];
        if(patch.name()=="inlet")
        {
            inlet.append(patchI);
            continue;
        }
        else if(patch.name()=="outlet")
        {
            outlet.append(patchI);
            continue;
        }
        else if(patch.name().substr(0,12)!="procBoundary")
            wall.append(patchI);
    }
    for(label patchI : inlet)
    {
        fvPatchField<vector>& inletPatch = boundaryField[patchI];
        fvPatchField<vector>* inletPatchPtr = &inletPatch;
        icoAdjointVelocityInletBC* cast_inletPatchPtr = dynamic_cast<icoAdjointVelocityInletBC*>(inletPatchPtr);
        if(cast_inletPatchPtr==nullptr)
        {
            const fvPatch& patch = inletPatch.patch();
            Pout<<"patch.name:"<<patch.name()<<Foam::nl;
            FatalErrorInFunction<<"Inlet Patch is not icoAdjointVelocityInletBC"<<exit(FatalError);
        }
        cast_inletPatchPtr->set_dJdp_Inlet(dJdp_Inlet);
    }
    for(label patchI : wall)
    {
        fvPatchField<vector>& wallPatch = boundaryField[patchI];
        fvPatchField<vector>* wallPatchPtr = &wallPatch;
        icoAdjointVelocityWallBC* cast_wallPatchPtr = dynamic_cast<icoAdjointVelocityWallBC*>(wallPatchPtr);
        if(cast_wallPatchPtr==nullptr)
        {
            const fvPatch& patch = wallPatch.patch();
            Pout<<"patch.name:"<<patch.name()<<Foam::nl;
            FatalErrorInFunction<<"Inlet Patch is not icoAdjointVelocityWallBC"<<exit(FatalError);
        }
        cast_wallPatchPtr->set_dJdp_Wall(dJdp_Wall);
    }
    for(label patchI : outlet)
    {
        fvPatchField<vector>& outletPatch = boundaryField[patchI];
        fvPatchField<vector>* outletPatchPtr = &outletPatch;
        icoAdjointVelocityOutletBC* cast_outletPatchPtr = dynamic_cast<icoAdjointVelocityOutletBC*>(outletPatchPtr);
        if(cast_outletPatchPtr==nullptr)
        {
            const fvPatch& patch = outletPatch.patch();
            Pout<<"patch.name:"<<patch.name()<<Foam::nl;
            FatalErrorInFunction<<"Outlet Patch is not icoAdjointVelocityOutletBC"<<exit(FatalError);
        }
        cast_outletPatchPtr->set_dJdu_Outlet(dJdu_uOutlet);
        cast_outletPatchPtr->set_nu(nu);
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::setAdjPBC
(
    std::function<Field<vector>(const icoAdjointPressureOutletBC&)> dJdu_pOutlet
)
{
    DynamicList<label> outlet;
    const fvBoundaryMesh& boundary = mesh.boundary();
    GeometricBoundaryField<scalar,fvPatchField,volMesh>& boundaryField = adj_p_.boundaryFieldRef();
    for(label patchI=0; patchI<boundary.size(); patchI++)
    {
        fvPatch const& patch = boundary[patchI];
        if(patch.name()=="outlet")
        {
            outlet.append(patchI);
            continue;
        }
    }
    for(label patchI : outlet)
    {
        fvPatchField<scalar>& outletPatch = boundaryField[patchI];
        fvPatchField<scalar>* outletPatchPtr = &outletPatch;
        icoAdjointPressureOutletBC* cast_outletPatchPtr = dynamic_cast<icoAdjointPressureOutletBC*>(outletPatchPtr);
        if(cast_outletPatchPtr==nullptr)
            FatalErrorInFunction<<"Patch is not icoAdjointPressureOutletBC"<<exit(FatalError);
        cast_outletPatchPtr->set_dJdu_Outlet(dJdu_pOutlet);
        cast_outletPatchPtr->set_nu(nu);
        cast_outletPatchPtr->set_temperatureUsed(useAdjointTemperature & useTemperature);
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::setAdjTBC
(
    std::function<Field<scalar>(const icoAdjointTemperatureOutletBC&)> dJdT_Outlet
)
{
    if(adj_T_)
    {
        DynamicList<label> outlet;
        const fvBoundaryMesh& boundary = mesh.boundary();
        GeometricBoundaryField<scalar,fvPatchField,volMesh>& boundaryField = adj_T_->boundaryFieldRef();
        for(label patchI=0; patchI<boundary.size(); patchI++)
        {
            fvPatch const& patch = boundary[patchI];
            if(patch.name()=="outlet")
            {
                outlet.append(patchI);
                continue;
            }
        }
        for(label patchI : outlet)
        {
            fvPatchField<scalar>& outletPatch = boundaryField[patchI];
            fvPatchField<scalar>* outletPatchPtr = &outletPatch;
            icoAdjointTemperatureOutletBC* cast_outletPatchPtr = dynamic_cast<icoAdjointTemperatureOutletBC*>(outletPatchPtr);
            if(cast_outletPatchPtr==nullptr)
                FatalErrorInFunction<<"Patch is not icoAdjointTemperatureOutletBC"<<exit(FatalError);
            cast_outletPatchPtr->set_dJdT_Outlet(dJdT_Outlet);
            cast_outletPatchPtr->set_alpha(alpha);
        }
    }
}
