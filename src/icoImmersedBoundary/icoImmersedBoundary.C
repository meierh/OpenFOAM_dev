#include "icoImmersedBoundary.H"

Foam::solvers::icoImmersedBoundary::icoImmersedBoundary
(
    fvMesh& mesh
):
incompressibleFluid(mesh),
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
    IOobject structureIO("structureDict","structure",runTime,IOobject::MUST_READ,IOobject::NO_WRITE);
    if(!structureIO.filePath("",true).empty())
    {
        structureDict = std::make_unique<IOdictionary>(structureIO);
        ITstream rodTypeStream = structureDict->lookup("rodType");
        token rodTypeToken;
        rodTypeStream.read(rodTypeToken);
        if(!rodTypeToken.isWord())
        {
            Info<<"rodTypeToken:"<<rodTypeToken<<Foam::endl;
            Info<<"rodTypeToken:"<<rodTypeToken.typeName()<<Foam::endl;
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodType -- must be string"<<exit(FatalError);
        }
        word rodTypeWord = rodTypeToken.wordToken();
        if(rodTypeWord == "Line")
        {
            structure = std::make_unique<LineStructure>(mesh,*structureDict);
            useStructure = true;
        }
        else if(rodTypeWord == "CrossSection")
        {
            structure = std::make_unique<CrossSectionStructure>(mesh,*structureDict);
            useStructure = true;
        }
        else if(rodTypeWord == "None")
        {
            useStructure = false;
        }
        else
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodType -- valid {Line,CrossSection,Empty}"<<exit(FatalError);
    }
    else
    {
        useStructure = false;
    }
    
    create_VelocityForcing();
    create_Temperature();
    create_TemperatureForcing();
    create_Refiner(mesh);
    
    Info<<"--------------------------icoImmersedBoundary--------------------------"<<Foam::endl;
    Info<<"useStructure:"<<useStructure<<Foam::endl;
    Info<<"useVelocityForcing:"<<useVelocityForcing<<Foam::endl;
    Info<<"useTemperature:"<<useTemperature<<Foam::endl;
    Info<<"useTemperatureForcing:"<<useTemperatureForcing<<Foam::endl;
    Info<<"useRefinement:"<<useRefinement<<Foam::endl;
    Info<<"||||||||||||||||||||||||||icoImmersedBoundary||||||||||||||||||||||||||"<<Foam::endl;
    FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
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
            interaction_fU = std::make_unique<StaticVelocityPressureAction>(mesh,*structure,U_,*fU_);
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
        useVelocityForcing = false;
}

void Foam::solvers::icoImmersedBoundary::create_Temperature()
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
        useTemperature = true;
    }
    else
    {
        useTemperature = false;
    }
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

void Foam::solvers::icoImmersedBoundary::create_Refiner(fvMesh& mesh)
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
            if(!topoChangerTypeToken.isString())
                FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodType -- must be string"<<exit(FatalError);
            word topoChangerTypeWord = topoChangerTypeToken.stringToken();
            if(topoChangerTypeWord!="refiner")
                FatalErrorInFunction<<"Invalid topoChanger/type -- valid {refiner}"<<exit(FatalError);
            
            refine_ = std::make_unique<volScalarField>
            (
                IOobject(topoChangerTypeWord,runTime.name(),mesh,Foam::IOobject::MUST_READ,Foam::IOobject::AUTO_WRITE),
                mesh
            );
            
            //Set refine/stay/unrefine field values
            if(topoChangerDict.found("upperRefineLevel")) topoChangerDict.set("upperRefineLevel",1.5);
            else topoChangerDict.add("upperRefineLevel",1.5);
            
            if(topoChangerDict.found("lowerRefineLevel")) topoChangerDict.set("lowerRefineLevel",0.5);
            else topoChangerDict.add("lowerRefineLevel",0.5);
            
            if(topoChangerDict.found("unrefineLevel")) topoChangerDict.set("unrefineLevel",-0.5);
            else topoChangerDict.add("unrefineLevel",-0.5);
            
            if(dynamicMeshDict.found("meshRefiner"))
            {
                dictionary meshRefinerDict = dynamicMeshDict.subDict("meshRefiner");
                ITstream refinerTypeStream = dynamicMeshDict.lookup("refinerType");
                token refinerTypeToken;
                refinerTypeStream.read(refinerTypeToken);
                if(!refinerTypeToken.isString())
                    FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodType -- must be string"<<exit(FatalError);
                word refinerTypeWord = refinerTypeToken.stringToken();
                if(refinerTypeWord == "None")
                {
                    useRefinement = false;
                }
                else if(refinerTypeWord == "MarkerOnly")
                {
                    refinement_ = std::make_unique<MeshRefiner>(mesh,*structure,*refine_,dynamicMeshDict);
                    useRefinement = true;
                }
                else
                    FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodType -- valid {None,MarkerOnly}"<<exit(FatalError);
            }
            else
                FatalErrorInFunction<<"Invalid dynamicMeshDict configuration: Missing meshRefiner entry"<<exit(FatalError);
        }
        else
            useRefinement = false;
    }
    else
        useRefinement = false;
}

void Foam::solvers::icoImmersedBoundary::preSolve()
{
    if(interaction_fU)
        interaction_fU->preSolveMovement();
}

void Foam::solvers::icoImmersedBoundary::momentumPredictor()
{
    volVectorField& U(U_);
  
    tUEqn =
    (
        fvm::ddt(U) + fvm::div(phi, U) - fvm::laplacian(nu,U)
      /*
      + MRF.DDt(U)
      + momentumTransport->divDevSigma(U)
     ==
        fvModels().source(U)
      */
    );
    fvVectorMatrix& UEqn = tUEqn.ref();
    
    Info<<"Created Ueqn"<<Foam::endl;
  
    UEqn.relax();
  
    fvConstraints().constrain(UEqn);
    
    Info<<"useVelocityForcing:"<<useVelocityForcing<<Foam::endl;
  
    if (pimple.momentumPredictor())
    {
        if(useVelocityForcing)
        {
            volVectorField& fU = *fU_;
            UEqn_res = solve(UEqn == -fvc::grad(p) + fU);
            
            interaction_fU->solve();
         }
         else
            UEqn_res = solve(UEqn == -fvc::grad(p));

        fvConstraints().constrain(U);
    }
    Info<<"done momentumPredictor"<<Foam::endl;
}

void Foam::solvers::icoImmersedBoundary::correctPressure()
{
    Info<<"correctPressure"<<Foam::endl;
    
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

    if (p.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p);
        fvc::makeAbsolute(phiHbyA, U);
    }

    tmp<volScalarField> rAtU(rAU);

    if (pimple.consistent())
    {
        rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
        phiHbyA += fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    if (pimple.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU(), MRF);

    // Evaluate any volume sources
    //fvScalarMatrix p_rghEqnSource(fvModels().sourceProxy(p));

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        std::unique_ptr<fvScalarMatrix> pEqnPtr;
        if(useVelocityForcing)
        {
            volVectorField& fU = *fU_;
            pEqnPtr = std::make_unique<fvScalarMatrix>
            (
                fvm::laplacian(rAtU(), p)
                ==
                fvc::div(phiHbyA)
                //- p_rghEqnSource
                + fvc::div(rAtU()*fU)
            );
        }
        else
        {
            pEqnPtr = std::make_unique<fvScalarMatrix>
            (
                fvm::laplacian(rAtU(), p)
                ==
                fvc::div(phiHbyA)
                //- p_rghEqnSource
            );
        }
        fvScalarMatrix& pEqn = *pEqnPtr;
        
        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    continuityErrors();

    // Explicitly relax pressure for momentum corrector
    p.relax();

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
    //fvc::correctUf(Uf, U, phi, MRF);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
    
}

void Foam::solvers::icoImmersedBoundary::postSolve()
{
    Info<<"postSolve"<<Foam::endl;
    Info<<"T_:"<<static_cast<bool>(T_)<<Foam::endl;
    Info<<"fT_:"<<static_cast<bool>(fT_)<<Foam::endl;
    
    if(useTemperature)
    {
        volScalarField& T = *T_;
        fvScalarMatrix TEqn(fvm::ddt(T)+fvm::div(phi,T)-fvm::laplacian(alpha,T));
        
        do
        {
            if(useTemperatureForcing)
            {
                volScalarField& fT = *fT_;
                TEqn_res = solve(TEqn==fT);
                
                interaction_fT->solve();
            }
            else
            {
                TEqn_res = solve(TEqn); 
            }
        }
        while(TEqn_res.nIterations()>0);
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
