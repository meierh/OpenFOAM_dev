#include "icoEvalImmersedBoundary.H"

Foam::solvers::icoEvalImmersedBoundary::icoEvalImmersedBoundary
(
    argList& args,
    Parameter para,
    objectiveFunction obj
):
args(args),
para(para),
J(obj.J)
{
    std::unique_ptr<Time> timePtr = createTime(args);
    std::unique_ptr<fvMesh> meshPtr = createMesh(*timePtr);
    auto icoSolver = std::unique_ptr<icoAdjointImmersedBoundary>(new icoAdjointImmersedBoundary(*meshPtr,*timePtr,{para}));
    const std::unique_ptr<LineStructure>& structure = icoSolver->getStructure();
    if(!structure)
        FatalErrorInFunction<<"No structure set"<<exit(FatalError);
    scalarList values = structure->getParameterValue(para);
    if(values.size()!=1)
        FatalErrorInFunction<<"Parameter must be a single coefficient"<<exit(FatalError);
    parameterIniValue = values[0];
    
    if(Pstream::master())
    {
        recordFDFile = std::make_unique<std::ofstream>("JRecords");
        (*recordFDFile) << std::setprecision(20);
    }
        
    Info<<"--------------------------icoEvalImmersedBoundary--------------------------"<<Foam::endl;
    Info<<"parameter:"<<para.to_string()<<structure->getParameterValue(para)<<Foam::endl;
    Info<<"||||||||||||||||||||||||||icoEvalImmersedBoundary||||||||||||||||||||||||||"<<Foam::endl;
}

void Foam::solvers::icoEvalImmersedBoundary::Solve()
{
    std::unique_ptr<Time> timePtr = createTime(args);
    std::unique_ptr<fvMesh> meshPtr = createMesh(*timePtr);

    scalar epsilon = 0;
    IOobject structureIO("structureDict","constant",*timePtr,IOobject::MUST_READ,IOobject::NO_WRITE);
    if(!structureIO.filePath("",true).empty())
    {
        IOdictionary structureDict(structureIO);
        ITstream epsilonStream = structureDict.lookup("epsilon");
        token epsilonToken;
        epsilonStream.read(epsilonToken);
        if(!epsilonToken.isScalar())
        {
            Info<<"epsilonToken:"<<epsilonToken<<Foam::endl;
            Info<<"epsilonToken:"<<epsilonToken.typeName()<<Foam::endl;
            FatalErrorInFunction<<"Invalid entry in constant/structureDict/epsilon -- must be scalar"<<exit(FatalError);
        }
        epsilon = epsilonToken.scalarToken();
    }
    else
      FatalErrorInFunction<<"Missing structureDict file"<<exit(FatalError);

    auto icoSolver = std::unique_ptr<icoImmersedBoundary>(new icoImmersedBoundary(*meshPtr,*timePtr));

    std::unique_ptr<LineStructure>& structure = icoSolver->getStructure();
    structure->setParameterValue(para,{parameterIniValue+epsilon});
    structure->reInitializeMarkers(false,false);
    scalarList values = structure->getParameterValue(para);
    if(values.size()!=1)
        FatalErrorInFunction<<"Parameter must be a single coefficient"<<exit(FatalError);
    Info<<"Parameter set to:"<<values[0]<<Foam::nl;
    //icoSolver->SolvePrimalRepeated();
    icoSolver->Solve();
    scalar Jval = J(*icoSolver);

    if(Pstream::master())
    {
        Info<<"Write final"<<Foam::endl;
        std::vector<std::string> para = icoSolver->parametersToString();
        Info<<"Parameters"<<Foam::endl;
        (*recordFDFile)<<"Parameters"<<std::endl;
        for(auto& str : para)
        {
            Info<<str<<Foam::endl;
            (*recordFDFile)<<str<<std::endl;
        }
        Info<<Foam::endl;
        (*recordFDFile)<<std::endl;

        Info<<"val:"<<parameterIniValue<<"  epsilon:"<<epsilon<<"  J:"<<Jval<<Foam::endl;
        (*recordFDFile)<<"val:"<<parameterIniValue<<"  epsilon:"<<epsilon<<"  J:"<<Jval<<std::endl;
    }
}
