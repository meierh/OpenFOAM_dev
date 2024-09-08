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
)
{
    create_VelocityForcing();
    create_Temperature();
    create_TemperatureForcing();
    create_Structure();
    create_Refiner(mesh);    
}

void Foam::solvers::icoImmersedBoundary::create_Structure()
{
    IOobject structureIO("structureDict","structure",runTime,IOobject::MUST_READ,IOobject::NO_WRITE);
    if(!structureIO.filePath("",true).empty())
    {
        stuctureDict = std::make_unique<IOdictionary>(structureIO);
        ITstream rodTypeStream = stuctureDict->lookup("rodType");
        token rodTypeToken;
        rodTypeStream.read(rodTypeToken);
        if(!rodTypeToken.isString())
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodType -- must be string"<<exit(FatalError);
        word rodTypeWord = rodTypeToken.stringToken();
        if(rodTypeWord == "Line")
            structure = std::make_unique<LineStructure>(mesh,*stuctureDict);
        else if(rodTypeWord == "CrossSection")
            structure = std::make_unique<CrossSectionStructure>(mesh,*stuctureDict);
        else
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodType -- valid {Line,CrossSection}"<<exit(FatalError);
    }
    else
        Info<<"No structure read"<<Foam::endl;
}

void Foam::solvers::icoImmersedBoundary::create_VelocityForcing()
{
    IOobject fU_IOobj
    (
        "fU",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if(fU_IOobj.filePath("",true).empty())
        fU_ = std::make_unique<volVectorField>(fU_IOobj,mesh);
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
    if(T_IOobj.filePath("",true).empty())
    {
        T_ = std::make_unique<volScalarField>(T_IOobj,mesh);
        alpha = dimensionedScalar(transportProperties.lookup("alpha"));
    }
}

void Foam::solvers::icoImmersedBoundary::create_TemperatureForcing()
{
    IOobject fT_IOobj
    (
        "fT",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if(fT_IOobj.filePath("",true).empty() && useTemperature)
        fT_ = std::make_unique<volScalarField>(fT_IOobj,mesh);
}

void Foam::solvers::icoImmersedBoundary::momentumPredictor()
{
    volVectorField& U(U_);
  
    tUEqn =
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + momentumTransport->divDevSigma(U)
     ==
        fvModels().source(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();
  
    UEqn.relax();
  
    fvConstraints().constrain(UEqn);
  
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
}

void Foam::solvers::icoImmersedBoundary::correctPressure()
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
        + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf)
    );

    MRF.makeRelative(phiHbyA);

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
    fvScalarMatrix p_rghEqnSource(fvModels().sourceProxy(p));

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
                - p_rghEqnSource
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
                - p_rghEqnSource
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
    fvc::correctUf(Uf, U, phi, MRF);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}

void Foam::solvers::icoImmersedBoundary::postSolve()
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
