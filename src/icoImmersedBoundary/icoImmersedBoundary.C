#include "icoImmersedBoundary.H"

Foam::solvers::icoImmersedBoundary::icoImmersedBoundary
(
    fvMesh& mesh,
    Time& time
):
incompressibleFluid(mesh),
time(time),
pimpleCtlr(incompressibleFluid::pimple,time),
transportProperties
(
    IOobject
    (
        "transportProperties",runTime.constant(),mesh,IOobject::MUST_READ_IF_MODIFIED,IOobject::NO_WRITE
    )
),
nu("nu",dimensionSet(0,2,-1,0,0,0,0),transportProperties.lookup("nu")),
alpha("alpha",dimensionSet(0,2,-1,0,0,0,0),0)
{   
    IOobject structureIO("structureDict","constant",runTime,IOobject::MUST_READ,IOobject::NO_WRITE);
    Info<<"structureIO:"<<structureIO.name()<<Foam::endl;
    if(!structureIO.filePath("",true).empty())
    {
        structureDict = std::make_shared<IOdictionary>(structureIO);
        ITstream rodTypeStream = structureDict->lookup("rodType");
        token rodTypeToken;
        rodTypeStream.read(rodTypeToken);
        if(!rodTypeToken.isWord())
        {
            Info<<"rodTypeToken:"<<rodTypeToken<<Foam::endl;
            Info<<"rodTypeToken:"<<rodTypeToken.typeName()<<Foam::endl;
            FatalErrorInFunction<<"Invalid entry in constant/structureDict/rodType -- must be string"<<exit(FatalError);
        }
        word rodTypeWord = rodTypeToken.wordToken();
        if(rodTypeWord == "Line")
        {
            structure = std::make_unique<LineStructure>(mesh,structureDict);
            useStructure = true;
        }
        else if(rodTypeWord == "CrossSection")
        {
            structure = std::make_unique<CrossSectionStructure>(mesh,structureDict);
            useStructure = true;
        }
        else if(rodTypeWord == "None")
        {
            useStructure = false;
        }
        else
            FatalErrorInFunction<<"Invalid entry in constant/structureDict/rodType -- valid {Line,CrossSection,Empty}"<<exit(FatalError);
    }
    else
    {
        useStructure = false;
    }
    
    create_Refiner(mesh);
    create_VelocityForcing();
    create_Temperature();
    create_TemperatureForcing();
    
    setDeltaT(time,*this);
    
    connectSolverPerformance(pimpleCtlr);
    
    checkDimensions();
    
    Info<<"--------------------------icoImmersedBoundary--------------------------"<<Foam::endl;
    Info<<"useStructure:"<<useStructure<<Foam::endl;
    Info<<"useRefinement:"<<useRefinement<<Foam::endl;
    Info<<"useVelocityForcing:"<<useVelocityForcing<<Foam::endl;
    Info<<"useTemperature:"<<useTemperature<<Foam::endl;
    Info<<"useTemperatureForcing:"<<useTemperatureForcing<<Foam::endl;
    Info<<"||||||||||||||||||||||||||icoImmersedBoundary||||||||||||||||||||||||||"<<Foam::endl;
        
    create_Analysis();
}

void Foam::solvers::icoImmersedBoundary::create_VelocityForcing()
{
    if(useStructure)
    {
        useVelocityForcing = true;
        IOobject fU_IOobj
        (
            "fU",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );
        if(!fU_IOobj.filePath("",true).empty())
            fU_ = std::make_unique<volVectorField>(fU_IOobj,mesh);
        else
            fU_ = std::make_unique<volVectorField>("fU",U_);
        
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
        
        IOobject fvSolutionIO("fvSolution","system",runTime,IOobject::MUST_READ,IOobject::NO_WRITE);
        if(!fvSolutionIO.filePath("",true).empty())
        {            
            IOdictionary fvSolutionDict(fvSolutionIO);
            dictionary& pimpleDict = fvSolutionDict.subDict("PIMPLE");
            
            // Read delayedForcing
            if(pimpleDict.found("delayedForcing"))
            {
                ITstream delayedForcingStream = pimpleDict.lookup("delayedForcing");
                token delayedForcingToken;
                delayedForcingStream.read(delayedForcingToken);
                if(!delayedForcingToken.isWord())
                    FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/delayedForcing -- must be word"<<exit(FatalError);
                word delayedForcingTokenWord = delayedForcingToken.wordToken();
                if(delayedForcingTokenWord=="yes")
                {
                    delayedVelocityForcing = true;
                    delayedTemperatureForcing = true;
                }
                else if(delayedForcingTokenWord=="no")
                {
                    delayedVelocityForcing = false;
                    delayedTemperatureForcing = false;                    
                }
                else
                    FatalErrorInFunction<<"Invalid word in system/fvSolution/PIMPLE/delayedForcing -- must be {yes,no}"<<exit(FatalError);
                
                if(delayedForcingTokenWord=="yes")
                {
                    ITstream forcingStartStream = pimpleDict.lookup("forcingStart");
                    token forcingStartToken;
                    forcingStartStream.read(forcingStartToken);
                    if(!forcingStartToken.isScalar())
                        FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/forcingStart -- must be scalar"<<exit(FatalError);
                    velocityForcingStartTime = temperatureForcingStartTime = forcingStartToken.scalarToken();
                    
                    ITstream forcingRampUpStream = pimpleDict.lookup("forcingRampUp");
                    token forcingRampUpToken;
                    forcingRampUpStream.read(forcingRampUpToken);
                    if(!forcingRampUpToken.isScalar())
                        FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/forcingRampUp -- must be scalar"<<exit(FatalError);
                    velocityForcingRampUpTime = temperatureForcingRampUpTime = forcingRampUpToken.scalarToken();
                }
            }
        }
        else
            FatalErrorInFunction<<"Missing file in system/fvSolution"<<exit(FatalError);
    }
    else
        useVelocityForcing = false;
}

void Foam::solvers::icoImmersedBoundary::create_Temperature
(
    pimpleIBControl& pimpleCtlr
)
{
    IOobject T_IOobj
    (
        "T",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if(!T_IOobj.filePath("",true).empty())
    {
        T_ = std::make_unique<volScalarField>(T_IOobj,mesh);
        alpha.value() = dimensionedScalar(transportProperties.lookup("alpha")).value();
        pimpleCtlr.setTemperaturePerformance(&TEqn_res);
        useTemperature = true;
    }
    else
    {
        useTemperature = false;
    }
}

void Foam::solvers::icoImmersedBoundary::create_Temperature()
{
    create_Temperature(pimpleCtlr);
}

void Foam::solvers::icoImmersedBoundary::create_TemperatureForcing()
{   
    if(useStructure && useTemperature)
    {
        if(structureDict->found("rodHeating"))
        {
            if(!useTemperature)
                FatalErrorInFunction<<"Temperature forcing but no temperature fields!"<<exit(FatalError);
            useTemperatureForcing = true;
            
            IOobject fT_IOobj
            (
                "fT",
                runTime.name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );
            if(!fT_IOobj.filePath("",true).empty())
                fT_ = std::make_unique<volScalarField>(fT_IOobj,mesh);
            else
                fT_ = std::make_unique<volScalarField>("fT",p_);
            
            ITstream rodHeatingStream = structureDict->lookup("rodHeating");
            token rodHeatingToken;
            rodHeatingStream.read(rodHeatingToken);
            if(!rodHeatingToken.isWord())
                FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodHeating -- must be word"<<exit(FatalError);
            word rodHeatingWord = rodHeatingToken.wordToken();
            if(rodHeatingWord == "FixedTemperature")
            {
                interaction_fT = std::make_unique<FixedTemperatureAction>(mesh,*structure,*T_,*fT_,*structureDict);
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
        {
            Info<<"rodHeating not found"<<Foam::endl;
            useTemperatureForcing = false;
        }
    }
    else
        useTemperatureForcing = false;
}

void Foam::solvers::icoImmersedBoundary::create_Refiner
(
    fvMesh& mesh
)
{
    IOobject dynamicMeshDictIO("dynamicMeshDict","constant",runTime,IOobject::MUST_READ,IOobject::NO_WRITE);
    if(!dynamicMeshDictIO.filePath("",true).empty())
    {
        IOdictionary dynamicMeshDict(dynamicMeshDictIO);
        if(dynamicMeshDict.found("topoChanger"))
        {
            dictionary& topoChangerDict = dynamicMeshDict.subDict("topoChanger");
            ITstream topoChangerTypeStream = topoChangerDict.lookup("type");
            token topoChangerTypeToken;
            topoChangerTypeStream.read(topoChangerTypeToken);
            if(!topoChangerTypeToken.isWord())
            {
                Info<<"topoChangerTypeToken:"<<topoChangerTypeToken<<Foam::endl;
                FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChanger/type -- must be string"<<exit(FatalError);
            }
            word topoChangerTypeWord = topoChangerTypeToken.wordToken();
            if(topoChangerTypeWord!="refiner")
                FatalErrorInFunction<<"Invalid topoChanger/type -- valid {refiner}"<<exit(FatalError);
            
            word refineFieldName = "refineField";
            IOobject refineIO
            (
                refineFieldName,
                runTime.name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );
            if(!refineIO.filePath("",true).empty())
            {                
                refine_ = std::make_unique<volScalarField>(refineIO,mesh);
            }
            else
            {
                refine_ = std::make_unique<volScalarField>(refineFieldName,p_);
            }
            if(topoChangerDict.found("field")) topoChangerDict.set("upperRefineLevel",refineFieldName);
            else topoChangerDict.add("field",refineFieldName);
                        
            ITstream refinerTypeStream = topoChangerDict.lookup("refinerType");
            token refinerTypeToken;
            refinerTypeStream.read(refinerTypeToken);
            if(!refinerTypeToken.isWord())
                FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChangerDict/refinerType -- must be string"<<exit(FatalError);
            word refinerTypeWord = refinerTypeToken.wordToken();
            if(refinerTypeWord == "None")
            {
                useRefinement = false;
            }
            else if(refinerTypeWord == "MarkerOnly")
            {
                refinement_ = std::make_shared<MeshRefiner>(mesh,*structure,*refine_,dynamicMeshDict);
                useRefinement = true;
            }
            else if(refinerTypeWord == "MarkerAndFluid")
            {
                dictionary& fluidCriterionDict = topoChangerDict.subDict("fluidCriterion");
                ITstream typeStream = fluidCriterionDict.lookup("type");
                token typeToken;
                typeStream.read(typeToken);
                if(!typeToken.isWord())
                    FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChangerDict/fluidCriterion/type -- must be word"<<exit(FatalError);
                word type = typeToken.wordToken();
                
                if(type=="velocityMagnitude")
                {
                    refinement_ = std::make_shared<VelocityMagnitudeLaplacian>(mesh,*structure,*refine_,dynamicMeshDict,U_);
                }
                else if(type=="rmsVelocityCurvature")
                {
                    refinement_ = std::make_shared<RMSVelocityCurvature>(mesh,*structure,*refine_,dynamicMeshDict,U_);
                }
                else
                    FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChangerDict/fluidCriterion/type -- valid {velocityMagnitude}"<<exit(FatalError);                
                useRefinement = true;
            }
            else
                FatalErrorInFunction<<"Invalid entry in constant/dynamicMeshDict/topoChangerDict/refinerType -- valid {None,MarkerOnly,MarkerAndFluid}"<<exit(FatalError);
        }
        else
            useRefinement = false;
    }
    else
        useRefinement = false;    
}

void Foam::solvers::icoImmersedBoundary::checkDimensions()
{
    //check adj_U
    if(U_.dimensions() != dimensionSet(0,1,-1,0,0,0,0))
        FatalErrorInFunction<<"Wrong dimensions of U_, given:"<<U_.dimensions()<<" necessary:"<<dimensionSet(0,1,-1,0,0,0,0)<<exit(FatalError);
    
    //check adj_p
    if(p_.dimensions() != dimensionSet(0,2,-2,0,0,0,0))
        FatalErrorInFunction<<"Wrong dimensions of p, given:"<<p_.dimensions()<<" necessary:"<<dimensionSet(0,2,-2,0,0,0,0)<<exit(FatalError);

    //check adj_T_
    if(T_)
    {
        if(T_->dimensions() != dimensionSet(0,0,0,1,0,0,0))
            FatalErrorInFunction<<"Wrong dimensions of T, given:"<<T_->dimensions()<<" necessary:"<<dimensionSet(0,0,0,1,0,0,0)<<exit(FatalError);
    }
}

void Foam::solvers::icoImmersedBoundary::connectSolverPerformance
(
    pimpleIBControl& pimpleCtlr
)
{
    pimpleCtlr.setVelocityPerformance(&UEqn_res);
    pimpleCtlr.setPressurePerformance(&PEqn_res);
    if(useTemperature)
        pimpleCtlr.setTemperaturePerformance(&TEqn_res);
}

void Foam::solvers::icoImmersedBoundary::preSolve
(
    pimpleIBControl& pimpleCtlr
)
{
    incompressibleFluid::preSolve();
}

void Foam::solvers::icoImmersedBoundary::preMove
(
    pimpleIBControl& pimpleCtlr
)
{
    if(interaction_fU)
    {
        interaction_fU->preSolveMovement();
        interaction_fU->preSolveMarkerMeshAdaption();
    }
    if(delayedVelocityForcing)
    {
        if(velocityForcingStartTime+velocityForcingRampUpTime < runTime.value())
        {
            velocityForcingFactor = 1;
        }
        else if(velocityForcingStartTime < runTime.value())
        {
            scalar timeSinceForcingStart = runTime.value()-velocityForcingStartTime;
            velocityForcingFactor = timeSinceForcingStart/velocityForcingRampUpTime;
            Info<<"Ramp up velocityForcingFactor:"<<velocityForcingFactor<<Foam::nl;
        }
        else
            velocityForcingFactor = 0;
    }
    if(delayedTemperatureForcing)
    {
        if(temperatureForcingStartTime+temperatureForcingRampUpTime < runTime.value())
        {
            temperatureForcingFactor = 1;
        }
        else if(temperatureForcingStartTime < runTime.value())
        {
            scalar timeSinceForcingStart = runTime.value()-temperatureForcingStartTime;
            temperatureForcingFactor = timeSinceForcingStart/temperatureForcingRampUpTime;
            Info<<"Ramp up temperatureForcingFactor:"<<temperatureForcingFactor<<Foam::nl;
        }
        else
            temperatureForcingFactor = 0;
    }
}

void Foam::solvers::icoImmersedBoundary::moveMesh
(
    pimpleIBControl& pimpleCtlr
)
{
}

void Foam::solvers::icoImmersedBoundary::motionCorrector
(
    pimpleIBControl& pimpleCtlr
)
{
}

void Foam::solvers::icoImmersedBoundary::prePredictor
(
    pimpleIBControl& pimpleCtlr
)
{
    incompressibleFluid::prePredictor();
}

void Foam::solvers::icoImmersedBoundary::momentumPredictor
(
    pimpleIBControl& pimpleCtlr
)
{
    volVectorField& U(U_);
  
    if(useVelocityForcing)
    {
        volVectorField& fU = *fU_;
        tUEqn = (fvm::ddt(U)+fvm::div(phi,U)-fvm::laplacian(nu,U)-velocityForcingFactor*fU);
    }
    else
    {
        tUEqn = (fvm::ddt(U)+fvm::div(phi,U)-fvm::laplacian(nu,U));
    }

    fvVectorMatrix& UEqn = tUEqn.ref();
      
    UEqn.relax();
  
    fvConstraints().constrain(UEqn);
      
    if (pimpleCtlr.momentumPredictor())
    {
        UEqn_res = solve(UEqn == -fvc::grad(p));
        if(useVelocityForcing && runTime.timeIndex()>1)
        {
            interaction_fU->solve(pimpleCtlr.firstIter(),runTime.value(),pimpleCtlr.finalIter());
        }

        fvConstraints().constrain(U);
    }
}

void Foam::solvers::icoImmersedBoundary::pressureCorrector
(
    pimpleIBControl& pimpleCtlr
)
{
    while (pimpleCtlr.correct())
    {
        correctPressure(pimpleCtlr);
    }
    tUEqn.clear();
}

void Foam::solvers::icoImmersedBoundary::correctPressure
(
    pimpleIBControl& pimpleCtlr
)
{   
    volScalarField& p(p_);
    volVectorField& U(U_);
    surfaceScalarField& phi(phi_);

    fvVectorMatrix& UEqn = tUEqn.ref();

    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
        + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
    );

    //MRF.makeRelative(phiHbyA);
    /*
    if (p.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p);
        fvc::makeAbsolute(phiHbyA, U);
    }
    */

    tmp<volScalarField> rAtU(rAU);

    if (pimpleCtlr.consistent())
    {
        rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
        phiHbyA += fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    if (pimpleCtlr.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU()/*,MRF*/);

    // Evaluate any volume sources
    //fvScalarMatrix p_rghEqnSource(fvModels().sourceProxy(p));

    // Non-orthogonal pressure corrector loop
    while (pimpleCtlr.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(fvm::laplacian(rAtU(),p)==fvc::div(phiHbyA));
        
        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        pEqn.solve();

        if (pimpleCtlr.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    continuityErrors();

    // Explicitly relax pressure for momentum corrector
    p.relax();
    U = HbyA - rAtU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvConstraints().constrain(U);

    // Correct Uf if the mesh is moving
    //fvc::correctUf(Uf, U, phi, MRF);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
    
}

void Foam::solvers::icoImmersedBoundary::postCorrector
(
    pimpleIBControl& pimpleCtlr
)
{
    incompressibleFluid::postCorrector();
    if(interaction_fU)
    {
        interaction_fU->subTimestepStructureMovement();
        interaction_fU->subTimestepMarkerMeshAdaption();
    }
}

void Foam::solvers::icoImmersedBoundary::postSolve
(
    pimpleIBControl& pimpleCtlr
)
{    
    if(useTemperature)
    {
        volScalarField& T = *T_;
        fvScalarMatrix TEqn(fvm::ddt(T)+fvm::div(phi,T)-fvm::laplacian(alpha,T));
        
        do
        {
            if(useTemperatureForcing)
            {
                volScalarField& fT = *fT_;
                TEqn_res = solve(TEqn==temperatureForcingFactor*fT);
                interaction_fT->solve(pimpleCtlr.firstHeatIter(),runTime.value(),pimpleCtlr.finalHeatIter());
            }
            else
            {
                TEqn_res = solve(TEqn); 
            }
        }
        while(pimpleCtlr.temperatureLoop());
    }
}

void Foam::solvers::icoImmersedBoundary::store()
{
    scalar time = mesh.time().value();
    timeToIndex.insert({time,times.size()});
    times.push_back(time);
    
    //storage_Mesh[time] = mesh;
    
    U_.storeOldTimes();
    p_.storeOldTimes();
    if(useVelocityForcing)
        interaction_fU->store();
    
    if(useTemperature)
        T_->storeOldTimes();
    if(useTemperatureForcing)
        interaction_fT->store();
}

void Foam::solvers::icoImmersedBoundary::setToTime(scalar time)
{
    if(timeToIndex.find(time)==timeToIndex.end())
        FatalErrorInFunction<<"No data exists for "<<time<<exit(FatalError);
    label index = timeToIndex[time];
    
    U_.oldTimeRef(index);
    p_.oldTimeRef(index);
    if(useVelocityForcing)
        interaction_fU->setToTime(time);    

    if(useTemperature)
    {
        T_->oldTimeRef(index);
    }
    if(useTemperatureForcing)
        interaction_fT->setToTime(time);
}

void Foam::solvers::icoImmersedBoundary::oneTimestep
(
    pimpleIBControl& pimpleCtlr
)
{
    Info<<"--------------------------------------- Solve Primal ---------------------------------------"<<Foam::nl;
    preMove(pimpleCtlr);
    while (pimpleCtlr.momentumLoopOrig())
    {
        moveMesh(pimpleCtlr);
        motionCorrector(pimpleCtlr);
        //fvModels().correct();
        prePredictor(pimpleCtlr);
        momentumPredictor(pimpleCtlr);
        thermophysicalPredictor();
        pressureCorrector(pimpleCtlr);
        postCorrector(pimpleCtlr);
    }
    postSolve(pimpleCtlr);
    Info<<"OneTimestep done"<<Foam::nl;
}

void Foam::solvers::icoImmersedBoundary::oneTimestep()
{
    oneTimestep(pimpleCtlr);
}

void Foam::solvers::icoImmersedBoundary::Solve()
{   
    while (pimpleCtlr.run(time))
    {
        Info<< " Presolve Time = " << runTime.userTimeName() << nl << endl;
        pimpleCtlr.read();
        preSolve(pimpleCtlr);
        // Adjust the time-step according to the solver maxDeltaT
        adjustDeltaT(time, *this);
        time++;
        Info<< "Time = " << runTime.userTimeName() << nl << endl;
            
        oneTimestep();

        write_Analysis();
        runTime.write();

        Info<<"ExecutionTime = "<<runTime.elapsedCpuTime()<<" s"<<"  ClockTime = "<<runTime.elapsedClockTime()<<" s"<<nl<< nl;
    }
}

Foam::scalar Foam::solvers::icoImmersedBoundary::minDeltaT() const
{
    const IOdictionary& controlDict = runTime.controlDict();
    scalar minDeltaT = std::numeric_limits<scalar>::min();
    if(controlDict.found("minDeltaT"))
    {
        minDeltaT = controlDict.lookup<scalar>("minDeltaT");
    }
    return minDeltaT;
}

void Foam::solvers::icoImmersedBoundary::create_Analysis()
{
    IOobject cellSizes_IOobj
    (
        "cellSizes",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if(!cellSizes_IOobj.filePath("",true).empty())
        cellSizes = std::make_unique<volScalarField>(cellSizes_IOobj,mesh);
    
    IOobject cellMarkerCount_IOobj
    (
        "cellMarkerCount",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if(!cellMarkerCount_IOobj.filePath("",true).empty())
        cellMarkerCount = std::make_unique<volScalarField>(cellMarkerCount_IOobj,mesh);
    
    IOobject cellMarkerCharacSize_IOobj
    (
        "cellMarkerCharacSize",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if(!cellMarkerCharacSize_IOobj.filePath("",true).empty())
        cellMarkerCharacSize = std::make_unique<volScalarField>(cellMarkerCharacSize_IOobj,mesh);
    
}

void Foam::solvers::icoImmersedBoundary::write_Analysis()
{
    Info<<"write_Analysis"<<Foam::nl;
    if(structure)
    {
        if(cellSizes)
            structure->writeCellSizeField(*cellSizes);
        if(cellMarkerCount)
            structure->writeCellMarkerCountField(*cellMarkerCount);
        if(cellMarkerCharacSize)
            structure->writeCellMarkerCharacSizeField(*cellMarkerCharacSize);
    }
    Info<<"write_Analysis done"<<Foam::nl;
}
