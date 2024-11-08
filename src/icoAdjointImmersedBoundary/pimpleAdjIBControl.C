#include "pimpleAdjIBControl.H"

Foam::solvers::pimpleAdjIBControl::pimpleAdjIBControl
(
    pimpleNoLoopControl& pimple,
    Time& runTime
):
pimpleIBControl(pimple,runTime)
{
    readFromDict();
}

void Foam::solvers::pimpleAdjIBControl::readFromDict()
{
    pimpleIBControl::readFromDict();
    
    IOobject fvSolutionIO("fvSolution","system",runTime,IOobject::MUST_READ,IOobject::NO_WRITE);
    if(!fvSolutionIO.filePath("",true).empty())
    {
        IOdictionary fvSolutionDict(fvSolutionIO);
        
        dictionary& solverDict = fvSolutionDict.subDict("solvers");
        
        // Read velocityTolerance
        dictionary& adj_uDict = solverDict.subDict("adj_U");
        ITstream adj_uDictToleranceStream = adj_uDict.lookup("tolerance");
        token adj_uDictToleranceToken;
        adj_uDictToleranceStream.read(adj_uDictToleranceToken);
        if(!adj_uDictToleranceToken.isScalar())
            FatalErrorInFunction<<"Invalid entry in system/fvSolution/solvers/adj_U/tolerance -- must be scalar"<<exit(FatalError);
        adjVelocityTolerance = adj_uDictToleranceToken.scalarToken();
        
        // Read pressureTolerance
        dictionary& adj_pDict = solverDict.subDict("adj_p");
        ITstream adj_pDictToleranceStream = adj_pDict.lookup("tolerance");
        token adj_pDictToleranceToken;
        adj_pDictToleranceStream.read(adj_pDictToleranceToken);
        if(!adj_pDictToleranceToken.isScalar())
            FatalErrorInFunction<<"Invalid entry in system/fvSolution/solvers/adj_p/tolerance -- must be scalar"<<exit(FatalError);
        adjPressureTolerance = adj_pDictToleranceToken.scalarToken();
        
        // Read temperatureTolerance
        if(solverDict.found("adj_T"))
        {
            dictionary& adj_tDict = solverDict.subDict("adj_T");
            ITstream adj_tDictToleranceStream = adj_tDict.lookup("tolerance");
            token adj_tDictToleranceToken;
            adj_tDictToleranceStream.read(adj_tDictToleranceToken);
            if(!adj_tDictToleranceToken.isScalar())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/solvers/adj_T/tolerance -- must be scalar"<<exit(FatalError);
            adjTemperatureTolerance = adj_tDictToleranceToken.scalarToken();
            adjTemperatureToleranceSet = true;
        }
        
        dictionary& adj_pimpleDict = fvSolutionDict.subDict("adj_PIMPLE");
        
        // Read nOuterCorrectors
        if(adj_pimpleDict.found("nOuterCorrectors"))
        {
            ITstream adjnOuterCorrectorsStream = adj_pimpleDict.lookup("nOuterCorrectors");
            token adjnOuterCorrectorsToken;
            adjnOuterCorrectorsStream.read(adjnOuterCorrectorsToken);
            if(!adjnOuterCorrectorsToken.isLabel())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/adjnOuterCorrectors -- must be label"<<exit(FatalError);
            adjMinOuterIterations = adjnOuterCorrectorsToken.labelToken();
        }
        
        // Read minMomentumIterations
        if(adj_pimpleDict.found("minMomentumIterations"))
        {
            ITstream adjMinMomentumIterationsStream = adj_pimpleDict.lookup("minMomentumIterations");
            token adjMinMomentumIterationsToken;
            adjMinMomentumIterationsStream.read(adjMinMomentumIterationsToken);
            if(!adjMinMomentumIterationsToken.isLabel())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/adjMinMomentumIterations -- must be label"<<exit(FatalError);
            adjMinMomentumIterations = adjMinMomentumIterationsToken.labelToken();
            if(adjMinMomentumIterations<1)
                FatalErrorInFunction<<"adjMinMomentumIterations must be >= 1"<<exit(FatalError);
        }
        
        // Read momentumIterTolFrac
        if(adj_pimpleDict.found("momentumIterTolFrac"))
        {
            ITstream adjMomentumIterTolFracStream = adj_pimpleDict.lookup("momentumIterTolFrac");
            token adjMomentumIterTolFracToken;
            adjMomentumIterTolFracStream.read(adjMomentumIterTolFracToken);
            if(!adjMomentumIterTolFracToken.isLabel())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/adjMomentumIterTolFrac -- must be label"<<exit(FatalError);
            adjMomentumIterTolFrac = adjMomentumIterTolFracToken.labelToken();
        }
        
        // Read momentumIterTolFrac
        if(adj_pimpleDict.found("minTemperatureIterations"))
        {
            ITstream adjMinTemperatureIterationsStream = adj_pimpleDict.lookup("minTemperatureIterations");
            token adjMinTemperatureIterationsToken;
            adjMinTemperatureIterationsStream.read(adjMinTemperatureIterationsToken);
            if(!adjMinTemperatureIterationsToken.isLabel())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/adjMinTemperatureIterations -- must be label"<<exit(FatalError);
            adjMinTemperatureIterations = adjMinTemperatureIterationsToken.labelToken();
            if(adjMinTemperatureIterations<1)
                FatalErrorInFunction<<"adjMinTemperatureIterations must be >= 1"<<exit(FatalError);
        }
        
        if(adj_pimpleDict.found("temperatureIterTolFrac"))
        {
            ITstream adjTemperatureIterTolFracStream = adj_pimpleDict.lookup("temperatureIterTolFrac");
            token adjTemperatureIterTolFracToken;
            adjTemperatureIterTolFracStream.read(adjTemperatureIterTolFracToken);
            if(!adjTemperatureIterTolFracToken.isLabel())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/adjTemperatureIterTolFrac -- must be label"<<exit(FatalError);
            adjTemperatureIterTolFrac = adjTemperatureIterTolFracToken.labelToken();
        }
        
        if(adj_pimpleDict.found("delayedInitialIterations"))
        {
            ITstream adjDelayedInitialIterationsStream = adj_pimpleDict.lookup("delayedInitialIterations");
            token adjDelayedInitialIterationsToken;
            adjDelayedInitialIterationsStream.read(adjDelayedInitialIterationsToken);
            if(!adjDelayedInitialIterationsToken.isWord())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/delayedInitialIterations -- must be label"<<exit(FatalError);
            word adjDelayedInitialIterationsWord = adjDelayedInitialIterationsToken.wordToken();
            if(adjDelayedInitialIterationsWord=="yes")
                adjDelayedInitialConvergence = true;
            else if(adjDelayedInitialIterationsWord=="no")
                adjDelayedInitialConvergence = false;
            else
                FatalErrorInFunction<<"Invalid word in system/fvSolution/adj_PIMPLE/adjDelayedInitialIterations -- must be {yes,no}"<<exit(FatalError);
        }
    }
    else
        FatalErrorInFunction<<"Missing file in system/fvSolution"<<exit(FatalError);
}

bool Foam::solvers::pimpleAdjIBControl::run(Time& time)
{
    return pimpleIBControl::run(time);
}

bool Foam::solvers::pimpleAdjIBControl::adjMomentumLoop()
{    
    bool contLoop;
    if(adjDelayedInitialConvergence && adjMomentumInnerIteration>timeIteration)
    {
        contLoop = false;
    }
    else if(adjMomentumInnerIteration<adjMinMomentumIterations)
    {
        contLoop = true;
    }
    else
    {
        if(adjVelocityEqns==nullptr)
            FatalErrorInFunction<<"adjVelocityEqns not set!"<<exit(FatalError);
        bool adjUEqnConverged = adjVelocityEqns->converged();
        vector adjUInitialRes = adjVelocityEqns->initialResidual();
        bool adjUIterationConverged = true;
        for(label d=0; d<3; d++)
           adjUIterationConverged &= (adjUInitialRes[d] < adjVelocityTolerance*adjMomentumIterTolFrac);
        
        if(adjPressureEqns==nullptr)
            FatalErrorInFunction<<"adjPressureEqns not set!"<<exit(FatalError);
        bool adjPEqnConverged = adjPressureEqns->converged();
        scalar adjPInitialRes = adjPressureEqns->initialResidual();
        bool adjPIterationConverged = (adjPInitialRes < adjPressureTolerance*adjMomentumIterTolFrac);
        
        contLoop = adjUEqnConverged && adjUIterationConverged && adjPEqnConverged && adjPIterationConverged;
    }
    
    adjMomentumInnerIteration++;
    return contLoop;
}

bool Foam::solvers::pimpleAdjIBControl::adjMomentumPredictor()
{
    return pimple.momentumPredictor();
}

bool Foam::solvers::pimpleAdjIBControl::adjConsistent()
{
    return pimple.consistent();
}

Foam::label Foam::solvers::pimpleAdjIBControl::adjNCorrPiso()
{
    return pimple.nCorrPiso();
}

bool Foam::solvers::pimpleAdjIBControl::adjCorrect()
{
    return pimple.correct();
}

bool Foam::solvers::pimpleAdjIBControl::adjCorrectNonOrthogonal()
{
    return pimple.correctNonOrthogonal();
}

bool Foam::solvers::pimpleAdjIBControl::adjFinalNonOrthogonalIter()
{
    return pimple.finalNonOrthogonalIter();
}

bool Foam::solvers::pimpleAdjIBControl::adjTemperatureLoop()
{
    if(!adjTemperatureUsed)
        FatalErrorInFunction<<"adjTemperature not used but temperature loop called"<<exit(FatalError);
    
    if(!adjTemperatureToleranceSet)
        FatalErrorInFunction<<"adjTemperatureTolerance not set"<<exit(FatalError);
    
    bool contLoop;
    if(adjDelayedInitialConvergence && adjTemperatureInnerIteration>timeIteration)
    {
        contLoop = false;
    }
    if(adjTemperatureInnerIteration<adjMinTemperatureIterations)
    {
        contLoop = true;
    }
    else
    {        
        if(adjTemperatureEqns==nullptr)
            FatalErrorInFunction<<"adjTemperatureEqns not set!"<<exit(FatalError);
        bool adjTEqnConverged = adjTemperatureEqns->converged();
        scalar adjTInitialRes = adjTemperatureEqns->initialResidual();
        bool adjTIterationConverged = (adjTInitialRes < adjTemperatureTolerance*adjTemperatureIterTolFrac);
        
        contLoop = adjTEqnConverged && adjTIterationConverged;
    }
    
    adjTemperatureInnerIteration++;
    return contLoop;
}

void Foam::solvers::pimpleAdjIBControl::setAdjVelocityPerformance
(
    SolverPerformance<vector>* adjUEqn
)
{
    adjVelocityEqns = adjUEqn;
}

void Foam::solvers::pimpleAdjIBControl::setAdjPressurePerformance
(
    SolverPerformance<scalar>* adjPEqn
)
{
    adjPressureEqns = adjPEqn;
}

void Foam::solvers::pimpleAdjIBControl::setAdjTemperaturePerformance
(
    SolverPerformance<scalar>* adjTEqn
)
{
    adjTemperatureEqns = adjTEqn;
    adjTemperatureUsed = true;
}
