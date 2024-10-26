#include "pimpleIBControl.H"

Foam::solvers::pimpleIBControl::pimpleIBControl
(
    pimpleNoLoopControl& pimple,
    bool delayedInitialConvergence
):
pimpleSingleRegionControl(pimple),
delayedInitialConvergence(delayedInitialConvergence)
{
    timeIteration = 0;
}

bool Foam::solvers::pimpleIBControl::run
(
    Time& time
)
{
    timeIteration++;
    innerLoopIteration=0;
    return pimpleSingleRegionControl::run(time);
}

bool Foam::solvers::pimpleIBControl::loop()
{
    innerLoopIteration++;
    bool convergence = evalConvergence();
    bool pimpleSingleRegionControlLoop = pimpleSingleRegionControl::loop();
    bool combinedLoop = convergence | pimpleSingleRegionControlLoop;
    Info<<"||||||||||||||||||||||||||||||||||||||||||||||||||Loop  time:"<<timeIteration<<" (conv:"<<convergence<<" , pSRCL:"<<pimpleSingleRegionControlLoop<<") -> "<<combinedLoop<<Foam::endl;
    return combinedLoop;
}

bool Foam::solvers::pimpleIBControl::evalConvergence()
{
    return false;
    innerLoopIteration++;
    
    if(delayedInitialConvergence)
    {
        if(innerLoopIteration>timeIteration)
        {
            Info<<"-----------------------------------delayedInitialConvergence------------------------"<<Foam::endl;
            return false;
        }
    }
    
    bool allEqnsConverged = true;
    for(const SolverPerformance<vector>* perf : vecEqns)
        allEqnsConverged &= perf->converged();
    for(const SolverPerformance<scalar>* perf : scaEqns)
        allEqnsConverged &= perf->converged();
    //Info<<"allEqnsConverged:"<<allEqnsConverged<<Foam::endl;
    if(!allEqnsConverged)
        return true;
    
    bool allEqnsFinal = true;
    for(const SolverPerformance<vector>* perf : vecEqns)
    {
        vector nIterations = perf->nIterations();
        if(nIterations[0]>1)
            allEqnsFinal = false;
        if(nIterations[1]>1)
            allEqnsFinal = false;
        if(nIterations[2]>1)
            allEqnsFinal = false;
    }
    for(const SolverPerformance<scalar>* perf : scaEqns)
    {
        label nIterations = perf->nIterations();
        if(nIterations>1)
            allEqnsFinal = false;
    }
    //Info<<"allEqnsFinal:"<<allEqnsFinal<<Foam::endl;
    if(!allEqnsFinal)
        return true;   
    
    return false;
}

void Foam::solvers::pimpleIBControl::addEqnPerformance
(
    SolverPerformance<vector>* eqn
)
{
    vecEqns.append(eqn);
}

void Foam::solvers::pimpleIBControl::addEqnPerformance
(
    SolverPerformance<scalar>* eqn
)
{
    scaEqns.append(eqn);
}
