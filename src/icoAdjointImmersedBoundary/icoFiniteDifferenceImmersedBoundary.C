#include "icoFiniteDifferenceImmersedBoundary.H"

Foam::solvers::icoFiniteDifferenceImmersedBoundary::icoFiniteDifferenceImmersedBoundary
(
    fvMesh& mesh,
    Time& time,
    Parameter para,
    std::vector<scalar> percFD,
    objectiveFunction obj
):
mesh(mesh),
time(time),
J(obj.J)
{    
    icoSolver = std::unique_ptr<icoAdjointImmersedBoundary>(new icoAdjointImmersedBoundary(mesh,time,{}));
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
    
    Info<<"--------------------------icoFiniteDifferenceImmersedBoundary--------------------------"<<Foam::endl;
    Info<<"parameter:"<<para.to_string()<<structure->getParameterValue(para)<<Foam::endl;
    Info<<"||||||||||||||||||||||||||icoFiniteDifferenceImmersedBoundary||||||||||||||||||||||||||"<<Foam::endl;    
}

void Foam::solvers::icoFiniteDifferenceImmersedBoundary::Solve()
{
    scalar J_plus,J_minus,fdGradient;
    
    for(scalar eps : epsilons)
    {       
        Info<<"Start eps:"<<eps<<Foam::nl;
        
        // plus epsilon
        icoSolver = std::unique_ptr<icoAdjointImmersedBoundary>(new icoAdjointImmersedBoundary(mesh,time,{}));
        {
            std::unique_ptr<LineStructure>& structure_plus = icoSolver->getStructure();
            structure_plus->setParameterValue(para,{parameterIniValue+eps});
        }
        icoSolver->SolvePrimal();
        J_plus = J(*icoSolver);
        
        Info<<"Done eps:"<<eps<<" J_plus"<<Foam::nl;
        
        // minus epsilon
        icoSolver = std::unique_ptr<icoAdjointImmersedBoundary>(new icoAdjointImmersedBoundary(mesh,time,{}));
        {
            std::unique_ptr<LineStructure>& structure_minus = icoSolver->getStructure();
            structure_minus->setParameterValue(para,{parameterIniValue-eps});
        }
        icoSolver->SolvePrimal();
        J_minus = J(*icoSolver);
        
        Info<<"Done eps:"<<eps<<" J_minus"<<Foam::nl;
        
        fdGradient = (J_plus-J_minus)/(2*eps);
        
        (*recordFDFile)<<"val:"<<parameterIniValue<<"  eps:"<<eps<<"  +eps J:"<<J_plus<<"  -eps J:"<<J_minus<<"  fdgrad:"<<fdGradient<<std::endl;
    }
}
