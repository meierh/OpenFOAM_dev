#include "nonStopPimpleAdjIBControl.H"

Foam::solvers::nonStopPimpleAdjIBControl::nonStopPimpleAdjIBControl
(
    pimpleNoLoopControl& pimple,
    Time& runTime
):
pimpleAdjIBControl(pimple,runTime)
{
}

bool Foam::solvers::nonStopPimpleAdjIBControl::run(Time& time)
{
    if(nonStop)
        return false;

    bool converged = pimple_.converged();
    if(converged)
    {
        time.writeNow();
        nonStop = false;
    }
    else
    {
        pimple_.storePrevIterFields();
    }
    return time.run();
}
