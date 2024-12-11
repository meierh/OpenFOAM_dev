#include "SetupOpenFOAM.H"

 Foam::argList Foam::setRootCase(int argc, char *argv[])
 {
    Foam::argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }
    return args;
 }
 
std::unique_ptr<Foam::Time> Foam::createTime
(
    const Foam::argList& args
)
{
    Foam::Info<< "Create time\n" << Foam::endl;  
    auto timePtr = std::make_unique<Foam::Time>(Foam::Time::controlDictName, args);
    return timePtr;
}

std::unique_ptr<Foam::fvMesh> Foam::createMesh
(
    Foam::Time& runTime
)
{
    Foam::Info << "Create mesh for time = " << runTime.name() << Foam::endl;
    auto meshPtr = std::make_unique<Foam::fvMesh>(IOobject(fvMesh::defaultRegion,runTime.name(),runTime,IOobject::MUST_READ));
    return meshPtr;
}
 
