#include "icoFiniteDifferenceImmersedBoundary.H"

Foam::solvers::icoFiniteDifferenceImmersedBoundary::icoFiniteDifferenceImmersedBoundary
(
    argList& args,
    Parameter para,
    std::vector<scalar> percFD,
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
    
    std::sort(percFD.begin(),percFD.end());
    for(uint i=0; i<percFD.size(); i++)
        epsilons.push_back(parameterIniValue*(percFD[i]/100));
    
    if(Pstream::master())
    {
        recordFDFile = std::make_unique<std::ofstream>("fdRecords");
        (*recordFDFile) << std::setprecision(20);
    }
        
    Info<<"--------------------------icoFiniteDifferenceImmersedBoundary--------------------------"<<Foam::endl;
    Info<<"parameter:"<<para.to_string()<<structure->getParameterValue(para)<<Foam::endl;
    Info<<"||||||||||||||||||||||||||icoFiniteDifferenceImmersedBoundary||||||||||||||||||||||||||"<<Foam::endl;
}

void Foam::solvers::icoFiniteDifferenceImmersedBoundary::Solve()
{
    scalar fdGradient;
    std::vector<scalar> J_eps;
    label resultCount = 1;
    
    for(scalar eps : epsilons)
    {      
        if(Pstream::master())
            Info<<"Start eps:"<<eps<<Foam::nl;

        for(label sign : {1,-1})
        {
            std::unique_ptr<Time> timePtr = createTime(args);
            std::unique_ptr<fvMesh> meshPtr = createMesh(*timePtr);
            auto icoSolver = std::unique_ptr<icoAdjointImmersedBoundary>(new icoAdjointImmersedBoundary(*meshPtr,*timePtr,{para}));
            std::unique_ptr<LineStructure>& structure = icoSolver->getStructure();
            structure->setParameterValue(para,{parameterIniValue+(static_cast<scalar>(sign)*eps)});
            icoSolver->SolvePrimal([&](bool completed, const Time& runTime, Time& time)
            {
                if(completed)
                {
                    time.setTime(static_cast<scalar>(resultCount),resultCount);
                    runTime.write();
                }
            });
            scalarList values = structure->getParameterValue(para);
            if(values.size()!=1)
                FatalErrorInFunction<<"Parameter must be a single coefficient"<<exit(FatalError);
            if(Pstream::master())
                Info<<"Done "<<sign<<" para:"<<values[0]<<Foam::nl;
            J_eps.push_back(J(*icoSolver));
            resultCount++;
        }
        fdGradient = (J_eps[0]-J_eps[1])/(2*eps);
        
        if(Pstream::master())
        {
            Info<<"val:"<<parameterIniValue<<"  eps:"<<eps<<"  +eps J:"<<J_eps[0]<<"  -eps J:"<<J_eps[1]<<"  fdgrad:"<<fdGradient<<Foam::endl;
            (*recordFDFile)<<"val:"<<parameterIniValue<<"  eps:"<<eps<<"  +eps J:"<<J_eps[0]<<"  -eps J:"<<J_eps[1]<<"  fdgrad:"<<fdGradient<<std::endl;

        }
    }
}
