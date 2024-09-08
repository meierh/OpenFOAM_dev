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
    
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_momentumPredictor()
{
    volVectorField& adj_U(adj_U_);

    tadj_UEqn =
    (
        fvm::ddt(adj_U) + fvm::div(phi, U)
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
        solve(UEqn == -fvc::grad(p));

        fvConstraints().constrain(U);
    }
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_thermophysicalPredictor()
{
     while (pimple.correct())
     {
         adj_correctPressure();
     }
     tadj_UEqn.clear();
}

void Foam::solvers::icoAdjointImmersedBoundary::adj_pressureCorrector()
{
     while (pimple.correct())
     {
         adj_correctPressure();
     }
     tadj_UEqn.clear();
}
