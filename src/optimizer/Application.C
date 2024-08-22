#include "Application.H"
#include <memory>

Foam::Application::Application
(
    int argc,
    char *argv[]
):
argc(argc),
argv(argv),
args(Foam::argList(argc, argv)),
runTime(Foam::Time::controlDictName, args),
mesh(std::make_unique<Foam::fvMesh>
        (
            Foam::IOobject(Foam::fvMesh::defaultRegion,
            runTime.name(),
            runTime,
            Foam::IOobject::MUST_READ)
        )
)
{
    setRootCaseLists();
    createTime();
    createMesh();
}

void Foam::Application::setRootCaseLists()
{
    //listOptions();
    setRootCase();
    //listOutput();
}
/*
void Foam::Application::listOptions()
{
    argList::addBoolOption
    (
        "listSwitches",
        "List switches declared in libraries but not set in etc/controlDict"
    );
    argList::addBoolOption
    (
        "listRegisteredSwitches",
        "List switches registered for run-time modification"
    );
    argList::addBoolOption
    (
        "listUnsetSwitches",
        "List switches declared in libraries but not set in etc/controlDict"
    );

    #ifdef fvPatchField_H
    argList::addBoolOption
    (
        "listScalarBCs",
        "List scalar field boundary conditions (fvPatchField<scalar>)"
    );
    argList::addBoolOption
    (
        "listVectorBCs",
        "List vector field boundary conditions (fvPatchField<vector>)"
    );
    #endif

    #ifdef functionObject_H
    argList::addBoolOption
    (
        "listFunctionObjects",
        "List functionObjects"
    );
    #endif

    #ifdef fvOption_H
    argList::addBoolOption
    (
        "listFvOptions",
        "List fvOptions"
    );
    #endif

    #if defined(turbulentTransportModel_H) || defined(turbulentFluidThermoModel_H)
    argList::addBoolOption
    (
        "listTurbulenceModels",
        "List turbulenceModels"
    );
    #endif
}
*/
void Foam::Application::setRootCase()
{
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }
}
/*
void Foam::Application::listOutput()
{
    _listOptions = false ;

    if
    (
        args.optionFound("listSwitches")
    )
    {
        debug::listSwitches(args.optionFound("includeUnsetSwitches"));
        _listOptions = true;
    }

    if
    (
        args.optionFound("listRegisteredSwitches")
    )
    {
        debug::listRegisteredSwitches(args.optionFound("includeUnsetSwitches"));
        _listOptions = true;
    }

    #ifdef fvPatchField_H
    if (args.optionFound("listScalarBCs"))
    {
        Info<< "scalarBCs"
            << fvPatchField<scalar>::dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        _listOptions = true;
    }

    if (args.optionFound("listVectorBCs"))
    {
        Info<< "vectorBCs"
            << fvPatchField<vector>::dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        _listOptions = true;
    }
    #endif

    #ifdef functionObject_H
    if (args.optionFound("listFunctionObjects"))
    {
        Info<< "functionObjects"
            << functionObject::dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        _listOptions = true;
    }
    #endif

    #ifdef fvOption_H
    if (args.optionFound("listFvOptions"))
    {
        Info<< "fvOptions"
            << fv::option::dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        _listOptions = true;
    }
    #endif

    #ifdef turbulentTransportModel_H
    if (args.optionFound("listTurbulenceModels"))
    {
        Info<< "Turbulence models"
            << incompressible::turbulenceModel::
            dictionaryConstructorTablePtr_->sortedToc()
            << endl;

        Info<< "RAS models"
            << incompressible::RASModel::
            dictionaryConstructorTablePtr_->sortedToc()
            << endl;

        Info<< "LES models"
            << incompressible::LESModel::
            dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        _listOptions = true;
    }
    #elif defined(turbulentFluidThermoModel_H)
    if (args.optionFound("listTurbulenceModels"))
    {
        Info<< "Turbulence models"
            << compressible::turbulenceModel::
            dictionaryConstructorTablePtr_->sortedToc()
            << endl;

        Info<< "RAS models"
            << compressible::RASModel::
            dictionaryConstructorTablePtr_->sortedToc()
            << endl;

        Info<< "LES models"
            << compressible::LESModel::
            dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        _listOptions = true;
    }
    #endif

    if (_listOptions)
    {
        FatalErrorInFunction<<"exit(0)"<<exit(FatalError);
    }
}
*/
void Foam::Application::createTime()
{
    Foam::Info<< "Create time\n" << Foam::endl;
}

void Foam::Application::createMesh()
{
    Foam::Info
        << "Create mesh for time = "
        << runTime.name() << Foam::nl << Foam::endl;

    mesh = std::make_unique<Foam::fvMesh>
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.name(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
}

void Foam::Application::createFields()
{
    
}
