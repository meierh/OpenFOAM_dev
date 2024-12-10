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
    
    Info<<"--------------------------icoFiniteDifferenceImmersedBoundary--------------------------"<<Foam::endl;
    Info<<"parameter:"<<para.to_string()<<Foam::endl;
    Info<<"||||||||||||||||||||||||||icoFiniteDifferenceImmersedBoundary||||||||||||||||||||||||||"<<Foam::endl;    
}

void Foam::solvers::icoFiniteDifferenceImmersedBoundary::Solve()
{
    for(scalar eps : epsilons)
    {
        // plus epsilon
        icoSolver = std::unique_ptr<icoAdjointImmersedBoundary>(new icoAdjointImmersedBoundary(mesh,time,{}));
        std::unique_ptr<LineStructure>& structure = icoSolver->getStructure();
        structure->setParameterValue(para,{parameterIniValue+eps});
        icoSolver->SolvePrimal();
        scalar J_plus = J(*icoSolver);
        
        // minus epsilon
        icoSolver = std::unique_ptr<icoAdjointImmersedBoundary>(new icoAdjointImmersedBoundary(mesh,time,{}));
        std::unique_ptr<LineStructure>& structure = icoSolver->getStructure();
        structure->setParameterValue(para,{parameterIniValue-eps});
        icoSolver->SolvePrimal();
        scalar J_minus = J(*icoSolver);
        
        scalar fdGradient = (J_plus-J_minus)/(2*eps);
    }
}
