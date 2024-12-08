#include "pimpleAdjIBControl.H"

Foam::solvers::pimpleAdjIBControl::pimpleAdjIBControl
(
    pimpleNoLoopControl& pimple,
    Time& runTime
):
pimpleIBControl(pimple,runTime)
{
    readFromFvSolution();
}

void Foam::solvers::pimpleAdjIBControl::readFromFvSolution()
{
    pimpleIBControl::readFromFvSolution();
    
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

        // Read adjointMomentum
        if(adj_pimpleDict.found("adjointMomentum"))
        {
            dictionary& adjointHeatDict = adj_pimpleDict.subDict("adjointMomentum");
            ITstream toleranceStream = adjointHeatDict.lookup("tolerance");
            token toleranceToken;
            toleranceStream.read(toleranceToken);
            if(!toleranceToken.isScalar())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/adjointMomentum/tolerance -- must be scalar"<<exit(FatalError);
            adjMomentumTolerance = toleranceToken.scalarToken();
            if(adjMomentumTolerance < adjVelocityTolerance)
                FatalErrorInFunction<<"adjMomentumTolerance lower than adjVelocityTolerance: ("<<adjMomentumTolerance<<"<"<<adjVelocityTolerance<<")"<<exit(FatalError);
            if(adjMomentumTolerance < adjPressureTolerance)
                FatalErrorInFunction<<"adjMomentumTolerance lower than adjPressureTolerance: ("<<adjMomentumTolerance<<"<"<<adjPressureTolerance<<")"<<exit(FatalError);
        }
        else
            adjMomentumTolerance = std::max<scalar>(adjPressureTolerance,adjVelocityTolerance);
        
        if(adj_pimpleDict.found("nCorrectors"))
        {
            ITstream nCorrectorsStream = adj_pimpleDict.lookup("nCorrectors");
            token nCorrectorsToken;
            nCorrectorsStream.read(nCorrectorsToken);
            if(!nCorrectorsToken.isLabel())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/nCorrectors -- must be label"<<exit(FatalError);
            adjnCorrectors = nCorrectorsToken.labelToken();
            adjnCorrectors_set = true;
            if(adjnCorrectors<0)
                FatalErrorInFunction<<"adjnCorrectors lower than zero"<<exit(FatalError);
        }
        else
            adjnCorrectors_set = false;

        if(adj_pimpleDict.found("nNonOrthogonalCorrectors"))
        {
            ITstream nNonOrthogonalCorrectorsStream = adj_pimpleDict.lookup("nNonOrthogonalCorrectors");
            token nNonOrthogonalCorrectorsToken;
            nNonOrthogonalCorrectorsStream.read(nNonOrthogonalCorrectorsToken);
            if(!nNonOrthogonalCorrectorsToken.isLabel())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/adj_PIMPLE/nNonOrthogonalCorrectors -- must be label"<<exit(FatalError);
            adjnNonOrthogonalCorrectorsToken = nNonOrthogonalCorrectorsToken.labelToken();
            adjnNonOrthogonalCorrectorsToken_set = true;
            if(adjnCorrectors<0)
                FatalErrorInFunction<<"adjnCorrectors lower than zero"<<exit(FatalError);
        }
        else
            adjnCorrectors_set = false;
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
    if(adjVelocityEqns==nullptr)
        FatalErrorInFunction<<"adjVelocityEqns not set"<<exit(FatalError);
    if(adjPressureEqns==nullptr)
        FatalErrorInFunction<<"adjPressureEqns not set"<<exit(FatalError);
    
    bool contLoop = false;
    vector adj_u_iniTol = adjVelocityEqns->initialResidual();
    scalar adj_p_iniTol = adjPressureEqns->initialResidual();
    
    for(label dim=0; dim<3; dim++)
    {
        if(adj_u_iniTol[dim]>=adjMomentumTolerance)
            contLoop = true;
    }
    if(adj_p_iniTol>=adjMomentumTolerance)
        contLoop = true;

    /*
    Info<<"pimpleAdjIBControl::adjMomentumLoop"<<Foam::nl;
    
    bool contLoop;
    Info<<"||--adjMomentumInnerIteration:"<<adjMomentumInnerIteration<<" timeIteration:"<<timeIteration<<" adjMinMomentumIterations:"<<adjMinMomentumIterations<<" adjDelayedInitialConvergence:"<<adjDelayedInitialConvergence<<Foam::nl;
    if(adjMomentumInnerIteration==0)
    {
        contLoop = true;
        Info<<"||--First iteration"<<Foam::nl;
    }
    else if(adjDelayedInitialConvergence && adjMomentumInnerIteration>timeIteration)
    {
        contLoop = false;
        Info<<"||--Delayed stop"<<Foam::nl;
    }
    else if(adjMomentumInnerIteration<adjMinMomentumIterations)
    {
        contLoop = true;
        Info<<"||--Lower than minimum"<<Foam::nl;
    }
    else if(steady)
    {
        contLoop = false;
        Info<<"||--Steady"<<Foam::nl;        
    }
    else
    {
        if(adjVelocityEqns==nullptr)
            FatalErrorInFunction<<"adjVelocityEqns not set!"<<exit(FatalError);
        bool adjUEqnConverged = adjVelocityEqns->converged();
        Vector<label> adjUIterations = adjVelocityEqns->nIterations();
        bool adjUIterationConverged = (adjUIterations[0]<1) && (adjUIterations[1]<1) && (adjUIterations[2]<1);
        
        if(adjPressureEqns==nullptr)
            FatalErrorInFunction<<"adjPressureEqns not set!"<<exit(FatalError);
        bool adjPEqnConverged = adjPressureEqns->converged();
        label adjPIterations = adjPressureEqns->nIterations();
        bool adjPIterationConverged = (adjPIterations<1);
        
        contLoop = !(adjUEqnConverged && adjUIterationConverged && adjPEqnConverged && adjPIterationConverged);
        
        if(contLoop)
            Info<<"||--Not Converged  adjUIterations:"<<adjUIterations<<" adjPIterations:"<<adjPIterations<<Foam::nl;
        else
            Info<<"||--Converged  adjUIterations:"<<adjUIterations<<" adjPIterations:"<<adjPIterations<<Foam::nl;
    }
    
    adjMomentumInnerIteration++;
    */
    
    momentumIteration++;
    
    if(contLoop)
        Info<<"momentum iteration:"<<momentumIteration<<" adj_U iniRes "<<adj_u_iniTol<<" adj_p iniRes:"<<adj_p_iniTol<<" limit:"<<adjMomentumTolerance<<" not converged"<<Foam::nl<<Foam::nl;
    else
        Info<<"momentum iteration:"<<momentumIteration<<" adj_U iniRes "<<adj_u_iniTol<<" adj_p iniRes:"<<adj_p_iniTol<<" limit:"<<adjMomentumTolerance<<" converged"<<Foam::nl<<Foam::nl;
    
    adjnCorrectors_count = 0;
    
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
    adjnNonOrthogonalCorrectorsToken_count = 0;
    if(adjnCorrectors_set)
    {
        adjnCorrectors_count++;
        if(adjnCorrectors_count<=adjnCorrectors)
            return true;
        else 
            return false;
    }
    else
        return pimple.correct();
}

bool Foam::solvers::pimpleAdjIBControl::adjCorrectNonOrthogonal()
{
    if(adjnNonOrthogonalCorrectorsToken_set)
    {
        adjnNonOrthogonalCorrectorsToken_count++;
        if(adjnNonOrthogonalCorrectorsToken_count<=adjnNonOrthogonalCorrectorsToken)
            return true;
        else 
            return false;
    }
    else
        return pimple.correct();
    
    return pimple.correctNonOrthogonal();
}

bool Foam::solvers::pimpleAdjIBControl::adjFinalNonOrthogonalIter()
{
    return pimple.finalNonOrthogonalIter();
}

bool Foam::solvers::pimpleAdjIBControl::adjTemperatureLoop()
{
    return false;
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
}
