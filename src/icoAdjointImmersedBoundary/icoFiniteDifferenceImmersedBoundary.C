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
    
    if(Pstream::master())
        recordFDFile = std::make_unique<std::ofstream>("fdRecords");

    Info<<"--------------------------icoFiniteDifferenceImmersedBoundary--------------------------"<<Foam::endl;
    Info<<"parameter:"<<para.to_string()<<structure->getParameterValue(para)<<Foam::endl;
    Info<<"||||||||||||||||||||||||||icoFiniteDifferenceImmersedBoundary||||||||||||||||||||||||||"<<Foam::endl;

    icoSolver.release();
    meshPtr.release();
    timePtr.release();
}

void Foam::solvers::icoFiniteDifferenceImmersedBoundary::Solve()
{
    scalar fdGradient;
    std::vector<scalar> J_eps;
    scalar writeTime = 1;
    
    Info<<"epsilons:";
    for(scalar eps : epsilons)
        Info<<eps<<" ";

    std::unique_ptr<Time> timePtr = createTime(args);
    std::unique_ptr<fvMesh> meshPtr = createMesh(*timePtr);
    auto icoSolver = std::unique_ptr<icoAdjointImmersedBoundary>(new icoAdjointImmersedBoundary(*meshPtr,*timePtr,{para}));

    for(scalar eps : epsilons)
    {
        Info<<"Start eps:"<<eps<<Foam::nl;
        for(label sign : {1,-1})
        {
            std::unique_ptr<LineStructure>& structure = icoSolver->getStructure();
            structure->setParameterValue(para,{parameterIniValue+(static_cast<scalar>(sign)*eps)});
            structure->reInitializeMarkers(false,false);
            scalarList values = structure->getParameterValue(para);
            if(values.size()!=1)
                FatalErrorInFunction<<"Parameter must be a single coefficient"<<exit(FatalError);
            Info<<"Parameter set to:"<<values[0]<<Foam::nl;
            icoSolver->SolvePrimalRepeated([&](bool completed, const Time& runTime, Time& time)
            {
                /*
                if(false && completed)
                {
                    scalar convTime = time.userTimeValue();
                    writeTime = std::max(convTime,writeTime);
                    label intWriteTime = writeTime;
                    intWriteTime+=1;
                    writeTime = intWriteTime;
                    scalar deltaTMissing = writeTime-convTime;
                    time+=(deltaTMissing);
                    Pout<<"Print into timestep set to:"<<writeTime<<" being:"<<time.userTimeValue()<<" from "<<convTime<<Foam::endl;
                    runTime.write();
                }
                */
            });
            J_eps.push_back(J(*icoSolver));
            Info<<"Done single instance "<<sign<<" eps:"<<eps<<Foam::nl;
        }
        fdGradient = (J_eps[0]-J_eps[1])/(2*eps);
        if(Pstream::master())
            (*recordFDFile)<<"val:"<<parameterIniValue<<"  eps:"<<eps<<"  +eps J:"<<J_eps[0]<<"  -eps J:"<<J_eps[1]<<"  fdgrad:"<<fdGradient<<std::endl;
    }
}
