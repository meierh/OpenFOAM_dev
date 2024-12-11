#include "icoFiniteDifferenceImmersedBoundary.H"

Foam::solvers::icoFiniteDifferenceImmersedBoundary::icoFiniteDifferenceImmersedBoundary
(
    int argc,
    char *argv[],
    Parameter para,
    std::vector<scalar> percFD,
    objectiveFunction obj
):
para(para),
J(obj.J)
{
    std::unique_ptr<Time> timePtr = createTime(setRootCase(argc,argv));
    std::unique_ptr<fvMesh> meshPtr = createMesh(*timePtr);
    icoSolver = std::unique_ptr<icoAdjointImmersedBoundary>(new icoAdjointImmersedBoundary(*meshPtr,*timePtr,{para}));
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
    
    recordFDFile = std::make_unique<std::ofstream>("fdRecords");
    
    icoSolver->checkIOObjects();
    
    
    
    
    Info<<"--------------------------icoFiniteDifferenceImmersedBoundary--------------------------"<<Foam::endl;
    Info<<"parameter:"<<para.to_string()<<structure->getParameterValue(para)<<Foam::endl;
    Info<<"||||||||||||||||||||||||||icoFiniteDifferenceImmersedBoundary||||||||||||||||||||||||||"<<Foam::endl;

    icoSolver.release();
    meshPtr.release();
    timePtr.release();
}

void Foam::solvers::icoFiniteDifferenceImmersedBoundary::Solve(int argc,char *argv[])
{
    scalar fdGradient;
    std::vector<scalar> J_eps;
    label resultCount = 1;
    
    for(scalar eps : epsilons)
    {       
        Info<<"Start eps:"<<eps<<Foam::nl;
        
        for(label sign : {1,-1})
        {
            // plus epsilon
            std::unique_ptr<Time> timePtr = createTime(setRootCase(argc,argv));
            std::unique_ptr<fvMesh> meshPtr = createMesh(*timePtr);
            auto icoSolver = std::unique_ptr<icoAdjointImmersedBoundary>(new icoAdjointImmersedBoundary(*meshPtr,*timePtr,{para}));
            {
                std::unique_ptr<LineStructure>& structure_plus = icoSolver->getStructure();
                structure_plus->setParameterValue(para,{parameterIniValue+(static_cast<scalar>(sign)*eps)});
            }
            icoSolver->SolvePrimal([&](bool completed, const Time& runTime, Time& time)
            {
                if(completed)
                {
                    time.setTime(static_cast<scalar>(resultCount),resultCount);
                    runTime.write();
                }
            });
            J_eps.push_back(J(*icoSolver));
            resultCount++;
            Info<<"Done "<<sign<<" eps:"<<eps<<Foam::nl;
        }
        fdGradient = (J_eps[0]-J_eps[1])/(2*eps);
        
        (*recordFDFile)<<"val:"<<parameterIniValue<<"  eps:"<<eps<<"  +eps J:"<<J_eps[0]<<"  -eps J:"<<J_eps[1]<<"  fdgrad:"<<fdGradient<<std::endl;
    }
}
