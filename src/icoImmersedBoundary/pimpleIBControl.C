#include "pimpleIBControl.H"

Foam::solvers::pimpleIBControl::pimpleIBControl
(
    pimpleNoLoopControl& pimple,
    Time& runTime
):
pimpleSingleRegionControl(pimple),
pimple(pimple),
runTime(runTime),
timeIteration(0)
{
    readFromFvSolution();
}

void Foam::solvers::pimpleIBControl::readFromFvSolution()
{
    IOobject fvSolutionIO("fvSolution","system",runTime,IOobject::MUST_READ,IOobject::NO_WRITE);
    if(!fvSolutionIO.filePath("",true).empty())
    {
        IOdictionary fvSolutionDict(fvSolutionIO);
        
        dictionary& solverDict = fvSolutionDict.subDict("solvers");
        
        // Read velocityTolerance
        dictionary& uDict = solverDict.subDict("U");
        ITstream uDictToleranceStream = uDict.lookup("tolerance");
        token uDictToleranceToken;
        uDictToleranceStream.read(uDictToleranceToken);
        if(!uDictToleranceToken.isScalar())
            FatalErrorInFunction<<"Invalid entry in system/fvSolution/solvers/U/tolerance -- must be scalar"<<exit(FatalError);
        velocityTolerance = uDictToleranceToken.scalarToken();
        
        // Read pressureTolerance
        dictionary& pDict = solverDict.subDict("p");
        ITstream pDictToleranceStream = pDict.lookup("tolerance");
        token pDictToleranceToken;
        pDictToleranceStream.read(pDictToleranceToken);
        if(!pDictToleranceToken.isScalar())
            FatalErrorInFunction<<"Invalid entry in system/fvSolution/solvers/p/tolerance -- must be scalar"<<exit(FatalError);
        pressureTolerance = pDictToleranceToken.scalarToken();
        
        // Read temperatureTolerance
        /*
        if(solverDict.found("T"))
        {
            dictionary& tDict = solverDict.subDict("T");
            ITstream tDictToleranceStream = tDict.lookup("tolerance");
            token tDictToleranceToken;
            tDictToleranceStream.read(tDictToleranceToken);
            if(!tDictToleranceToken.isScalar())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/solvers/T/tolerance -- must be scalar"<<exit(FatalError);
            temperatureTolerance = tDictToleranceToken.scalarToken();
            temperatureUsed = true;
            temperatureToleranceSet = true;
        }
        */
        
        dictionary& pimpleDict = fvSolutionDict.subDict("PIMPLE");
        
        // Read nOuterCorrectors
        if(pimpleDict.found("nOuterCorrectors"))
        {
            ITstream nOuterCorrectorsStream = pimpleDict.lookup("nOuterCorrectors");
            token nOuterCorrectorsToken;
            nOuterCorrectorsStream.read(nOuterCorrectorsToken);
            if(!nOuterCorrectorsToken.isLabel())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/nOuterCorrectors -- must be label"<<exit(FatalError);
            minOuterIterations = nOuterCorrectorsToken.labelToken();
            nOuterTemperatureCorrector = minOuterIterations;
            nOuterTemperatureCorrectorSet = true;
        }
        
        // Read minMomentumIterations
        if(pimpleDict.found("minMomentumIterations"))
        {
            ITstream minMomentumIterationsStream = pimpleDict.lookup("minMomentumIterations");
            token minMomentumIterationsToken;
            minMomentumIterationsStream.read(minMomentumIterationsToken);
            if(!minMomentumIterationsToken.isLabel())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/minMomentumIterations -- must be label"<<exit(FatalError);
            minMomentumIterations = minMomentumIterationsToken.labelToken();
            if(minMomentumIterations<1)
                FatalErrorInFunction<<"minMomentumIterations must be >= 1"<<exit(FatalError);
        }
        
        // Read momentumIterTolFrac
        if(pimpleDict.found("momentumIterTolFrac"))
        {
            ITstream momentumIterTolFracStream = pimpleDict.lookup("momentumIterTolFrac");
            token momentumIterTolFracToken;
            momentumIterTolFracStream.read(momentumIterTolFracToken);
            if(!momentumIterTolFracToken.isScalar())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/momentumIterTolFrac -- must be scalar"<<exit(FatalError);
            momentumIterTolFrac = momentumIterTolFracToken.scalarToken();
        }
        
        // Read momentumIterTolFrac
        /*
        if(pimpleDict.found("minTemperatureIterations"))
        {
            ITstream minTemperatureIterationsStream = pimpleDict.lookup("minTemperatureIterations");
            token minTemperatureIterationsToken;
            minTemperatureIterationsStream.read(minTemperatureIterationsToken);
            if(!minTemperatureIterationsToken.isLabel())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/minTemperatureIterations -- must be label"<<exit(FatalError);
            minTemperatureIterations = minTemperatureIterationsToken.labelToken();
            if(minTemperatureIterations<1)
                FatalErrorInFunction<<"minTemperatureIterations must be >= 1"<<exit(FatalError);
        }
        */
        
        /*
        if(pimpleDict.found("temperatureIterTolFrac"))
        {
            ITstream temperatureIterTolFracStream = pimpleDict.lookup("temperatureIterTolFrac");
            token temperatureIterTolFracToken;
            temperatureIterTolFracStream.read(temperatureIterTolFracToken);
            if(!temperatureIterTolFracToken.isScalar())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/temperatureIterTolFrac -- must be scalar"<<exit(FatalError);
            temperatureIterTolFrac = temperatureIterTolFracToken.scalarToken();
        }
        */
        
        if(pimpleDict.found("advancedInitialConvergence"))
        {
            ITstream advancedInitialConvergenceStream = pimpleDict.lookup("advancedInitialConvergence");
            token advancedInitialConvergenceToken;
            advancedInitialConvergenceStream.read(advancedInitialConvergenceToken);
            if(!advancedInitialConvergenceToken.isWord())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/advancedInitialConvergence -- must be word"<<exit(FatalError);
            word advancedInitialConvergenceWord = advancedInitialConvergenceToken.wordToken();
            if(advancedInitialConvergenceWord=="yes")
                advancedInitialConvergence = true;
            else if(advancedInitialConvergenceWord=="no")
                advancedInitialConvergence = false;
            else
                FatalErrorInFunction<<"Invalid word in system/fvSolution/PIMPLE/advancedInitialConvergence -- must be {yes,no}"<<exit(FatalError);
        }
        
        ITstream solutionTypeStream = pimpleDict.lookup("solutionType");
        token solutionTypeToken;
        solutionTypeStream.read(solutionTypeToken);
        if(!solutionTypeToken.isWord())
            FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/solutionType -- must be label"<<exit(FatalError);
        word solutionTypeWord = solutionTypeToken.wordToken();
        if(solutionTypeWord=="steady")
            steady = true;
        else if(solutionTypeWord=="unsteady")
            steady = false;
        else
            FatalErrorInFunction<<"Invalid word in system/fvSolution/PIMPLE/solutionType -- must be {steady,unsteady}"<<exit(FatalError);
        
        if(pimpleDict.found("convergenceDelay"))
        {
            ITstream convergenceDelayStream = pimpleDict.lookup("convergenceDelay");
            token convergenceDelayToken;
            convergenceDelayStream.read(convergenceDelayToken);
            if(!convergenceDelayToken.isScalar())
                FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/convergenceDelay -- must be scalar"<<exit(FatalError);
            convergenceDelay = convergenceDelayToken.scalarToken();            
            convergenceDelaySet = true;
        }
        else
            convergenceDelaySet = false;
        
        if(pimpleDict.found("outerCorrectorResidualControl"))
        {
            dictionary& outerResCtrl = pimpleDict.subDict("outerCorrectorResidualControl");
            if(outerResCtrl.found("T"))
            {
                dictionary& TDict = outerResCtrl.subDict("T");
                ITstream TTimestepToleranceStream = TDict.lookup("tolerance");
                token TTimestepToleranceToken;
                TTimestepToleranceStream.read(TTimestepToleranceToken);
                if(!TTimestepToleranceToken.isScalar())
                    FatalErrorInFunction<<"Invalid entry in system/fvSolution/PIMPLE/outerCorrectorResidualControl/T/tolerance -- must be scalar"<<exit(FatalError);
                TTimestepTolerance = TTimestepToleranceToken.scalarToken();
                TTimestepToleranceSet = true;
            }
        }
    }
    else
        FatalErrorInFunction<<"Missing file in system/fvSolution"<<exit(FatalError);    
}

bool Foam::solvers::pimpleIBControl::run(Time& time)
{
    momentumInnerIteration = 0;
    iOuterTemperatureCorrector = 0;
    return pimpleSingleRegionControl::run(time);
}

bool Foam::solvers::pimpleIBControl::runNonStop(Time& time)
{
    momentumInnerIteration = 0;
    iOuterTemperatureCorrector = 0;
    if(resetActive)
    {
        pimple.storePrevIterFields();
        resetActive = false;
        return true;
    }
    else
    {
        timeIteration++;
        bool converged = pimple.converged();
        if(converged)
        {
            time.writeNow();
            return false;
        }
        else
        {
            pimple.storePrevIterFields();
            return time.loop();
        }
    }
}

bool Foam::solvers::pimpleIBControl::momentumLoopOrig()
{
    return pimpleSingleRegionControl::loop();
}

bool Foam::solvers::pimpleIBControl::momentumLoop()
{    
    bool contLoop;
    Info<<"||--momentumInnerIteration:"<<momentumInnerIteration<<" timeIteration:"<<timeIteration<<" minMomentumIterations:"<<minMomentumIterations<<" advancedInitialConvergence:"<<advancedInitialConvergence<<Foam::endl;
    if(momentumInnerIteration==0)
    {
        contLoop = true;
        Info<<"||--First iteration"<<Foam::nl;
    }
    else if(advancedInitialConvergence && momentumInnerIteration>timeIteration)
    {
        contLoop = false;
        Info<<"||--Advanced stop"<<Foam::nl;
    }
    else if(momentumInnerIteration<minMomentumIterations)
    {
        contLoop = true;
        Info<<"||--Lower than minimum"<<Foam::nl;
    }
    else if(steady)
    {
        contLoop = false;
        Info<<"||--Steady"<<Foam::nl;        
    }
    else if(convergenceDelaySet && convergenceDelay < runTime.value())
    {
        contLoop = false;
        Info<<"||--Convergence delay"<<Foam::nl;   
    }
    else
    {
        if(velocityEqns==nullptr)
            FatalErrorInFunction<<"velocityEqns not set!"<<exit(FatalError);
        bool uEqnConverged = velocityEqns->converged();
        Vector<label> uIterations = velocityEqns->nIterations();
        bool uIterationConverged = (uIterations[0]<1) && (uIterations[1]<1) && (uIterations[2]<1);
        
        if(pressureEqns==nullptr)
            FatalErrorInFunction<<"pressureEqns not set!"<<exit(FatalError);
        bool pEqnConverged = pressureEqns->converged();
        label pIterations = pressureEqns->nIterations();
        bool pIterationConverged = (pIterations<1);
        
        contLoop = !(uEqnConverged && uIterationConverged && pEqnConverged && pIterationConverged);
        
        if(contLoop)
            Info<<"||--Not Converged  uIterations:"<<uIterations<<" pIterations:"<<pIterations<<Foam::nl;
        else
            Info<<"||--Converged  uIterations:"<<uIterations<<" pIterations:"<<pIterations<<Foam::nl;
    }
    
    momentumInnerIteration++;   
    return contLoop;
}

bool Foam::solvers::pimpleIBControl::momentumPredictor()
{
    return pimple.momentumPredictor();
}

bool Foam::solvers::pimpleIBControl::consistent()
{
    return pimple.consistent();
}

Foam::label Foam::solvers::pimpleIBControl::nCorrPiso()
{
    return pimple.nCorrPiso();
}

bool Foam::solvers::pimpleIBControl::correct()
{
    bool pimpleCorrect = pimple.correct();
    return pimpleCorrect;
}

bool Foam::solvers::pimpleIBControl::correctNonOrthogonal()
{
    bool pimpleCorrectNonOrthogonal = pimple.correctNonOrthogonal();
    return pimpleCorrectNonOrthogonal;
}

bool Foam::solvers::pimpleIBControl::finalNonOrthogonalIter()
{
    return pimple.finalNonOrthogonalIter();
}

bool Foam::solvers::pimpleIBControl::temperatureLoop()
{
    if(!temperatureUsed)
        FatalErrorInFunction<<"Temperature not used but temperature loop called"<<exit(FatalError);
    if(!nOuterTemperatureCorrectorSet)
        FatalErrorInFunction<<"nOuterTemperatureCorrector not set"<<exit(FatalError);
    
    if(iOuterTemperatureCorrector<nOuterTemperatureCorrector)
    {
        iOuterTemperatureCorrector++;
        if(TTimestepToleranceSet && iOuterTemperatureCorrector>1)
        {
            return true;
            if(temperatureEqns==nullptr)
                FatalErrorInFunction<<"temperatureEqns not set!"<<exit(FatalError);
            //bool tEqnConverged = temperatureEqns->converged();
            scalar tInitialRes = temperatureEqns->initialResidual();
            bool tIterationConverged = (tInitialRes<TTimestepTolerance);
            if(tIterationConverged)
                return false;
            else
                return true;
        }
        else
            return true;
    }
    else
    {
        return false;
    }
}

/*
bool Foam::solvers::pimpleIBControl::temperatureLoop()
{
    Info<<"pimpleIBControl::temperatureLoop()"<<Foam::nl;
    
    if(!temperatureUsed)
        FatalErrorInFunction<<"Temperature not used but temperature loop called"<<exit(FatalError);
    
    if(!temperatureToleranceSet)
        FatalErrorInFunction<<"TemperatureTolerance not set"<<exit(FatalError);
    
    bool contLoop = false;
    if(temperatureInnerIteration==0)
    {
        contLoop = true;
        Info<<"||--First iteration"<<Foam::nl;
    }
    else if(advancedInitialConvergence && temperatureInnerIteration>timeIteration)
    {
        contLoop = false;
        Info<<"||--Advanced stop"<<Foam::nl;
    }
    else if(temperatureInnerIteration<minTemperatureIterations)
    {
        contLoop = true;
        Info<<"||--Lower than minimum"<<Foam::nl;
    }
    else if(steady)
    {
        contLoop = false;
        Info<<"||--Steady"<<Foam::nl;        
    }
    else if(convergenceDelaySet && convergenceDelay < runTime.value())
    {
        contLoop = false;
        Info<<"||--Convergence delay"<<Foam::nl;   
    }
    else
    {        
        if(temperatureEqns==nullptr)
            FatalErrorInFunction<<"temperatureEqns not set!"<<exit(FatalError);
        bool tEqnConverged = temperatureEqns->converged();
        scalar tInitialRes = temperatureEqns->initialResidual();
        bool tIterationConverged = (tInitialRes < temperatureTolerance*temperatureIterTolFrac);
        
        contLoop = tEqnConverged && tIterationConverged;
    }
    
    temperatureInnerIteration++;
    return contLoop;
}
*/

void Foam::solvers::pimpleIBControl::setVelocityPerformance(SolverPerformance<vector>* uEqn)
{
    velocityEqns = uEqn;
}

void Foam::solvers::pimpleIBControl::setPressurePerformance(SolverPerformance<scalar>* pEqn)
{
    pressureEqns = pEqn;
}

void Foam::solvers::pimpleIBControl::setTemperaturePerformance(SolverPerformance<scalar>* tEqn)
{
    temperatureEqns = tEqn;
    temperatureUsed = true;
}

void Foam::solvers::pimpleIBControl::reset()
{
    timeIteration=0;
    resetActive=true;
    corr_=0;
    converged_=false;
}
