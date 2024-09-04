#include "icoImmersedBoundary.H"

Foam::solvers::icoImmersedBoundary::icoImmersedBoundary
(
    fvMesh& mesh
):
incompressibleFluid(mesh),
fU_(create_fU()),
T_(create_T()),
fT_(create_fT())
{
}

Foam::volVectorField Foam::solvers::icoImmersedBoundary::create_fU()
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
    {
        volVectorField fU(fU_IOobj,mesh);
        return fU;
    }
    else
    {
        volVectorField fU("fU",U_);
        return fU;
    }
}

Foam::volScalarField Foam::solvers::icoImmersedBoundary::create_T()
{
    IOobject fU_IOobj
    (
        "T",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if(fU_IOobj.filePath("",true).empty())
    {
        volScalarField T(fU_IOobj,mesh);
        return T;
    }
    else
    {
        volScalarField T("T",p_);
        return T;
    }
}

Foam::volScalarField Foam::solvers::icoImmersedBoundary::create_fT()
{
    IOobject fT_IOobj
    (
        "fT",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if(fT_IOobj.filePath("",true).empty() && computeTemperature)
    {
        volScalarField fT(fT_IOobj,mesh);
        return fT;
    }
    else
    {
        volScalarField fT("fT",T_);
        return fT;
    }
}

void Foam::solvers::icoImmersedBoundary::createStructure()
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
        word rodTypeWord = rodType.stringToken();
        if(rodTypeWord == "Line")
            structure = std::make_unique<LineStructure>(mesh,stuctureDict);
        else if(rodTypeWord == "CrossSection")
            structure = std::make_unique<CrossSectionStructure>(mesh,stuctureDict);
        else
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodType -- valid {Line,CrossSection}"<<exit(FatalError);
    }
    else
        Info<<"No structure read"<<Foam::endl;
}
