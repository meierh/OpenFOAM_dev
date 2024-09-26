#include "Optimizer.H"

bool Foam::Optimizer::singleton = false;
nlopt::vfunc Foam::Optimizer::core_objFunc = nullptr;
bool Foam::Optimizer::masterDone = false;


bool Foam::OptimizerStop::activated = false;
bool Foam::OptimizerStop::singleton = false;

Foam::OptimizerStop::OptimizerStop
():
std::runtime_error("Master process done with optimizer")
{
    if(OptimizerStop::singleton)
        FatalErrorInFunction<<"Can not have two activations of OptimizerStop"<<exit(FatalError);
    OptimizerStop::singleton=true;
    OptimizerStop::activated=true;
}

bool Foam::OptimizerStop::active()
{
    return activated;
}

void Foam::OptimizerStop::reset()
{
    activated=false;
}

Foam::Optimizer::Optimizer
(
    nlopt::vfunc obj,
    label n,
    void* f_data
):
n(n),
opt(nlopt::LD_MMA, n)
{
    if(singleton)
        FatalErrorInFunction<<"Can not have two instances of optimizer"<<exit(FatalError);
    Optimizer::singleton = true;

    Optimizer::core_objFunc = obj;
    opt.set_min_objective(&(Optimizer::nlopt_objFunc), f_data);
    opt.set_xtol_rel(1e-4);
}

Foam::scalar Foam::Optimizer::run()
{
    Info<<"Run"<<Foam::endl;
    double minf;
    if(static_cast<label>(x_initial.size())!=n)
        FatalErrorInFunction<<"Dimension mismatch for initial variables"<<exit(FatalError);
    
    bool validEnd = false;
    nlopt::result result;
    while(!validEnd)
    {
        bool slaveProcStopped = false;
        try
        {
            result = opt.optimize(x_initial, minf);
        }
        catch(std::exception &e)
        {
            //Slave process stopped due to master done
            if(OptimizerStop::active())
            {
                if(Pstream::master())
                    FatalErrorInFunction<<"Master can not be stopped"<<exit(FatalError);
                slaveProcStopped = true;
                OptimizerStop::reset();
                break;
            }
            
            //Other exception
            FatalErrorInFunction<<"Optimizer internal fail"<<exit(FatalError);
        }
        
        if(Pstream::master())
        {
            Optimizer::masterDone = true;
        }
        else
        {
            if(!slaveProcStopped)
                continue;
        }
        
        Pstream::scatter(Optimizer::masterDone);
        if(!Optimizer::masterDone)
            FatalErrorInFunction<<"Master not done can not appear here"<<exit(FatalError);
        if(Optimizer::masterDone)
            break;
    }
    
    return minf;
}

void Foam::Optimizer::setBounds
(
    const std::vector<scalar>& lower,
    const std::vector<scalar>& upper
)
{
    if(lower.size()!=0)
    {
        if(static_cast<label>(lower.size())!=n)
            FatalErrorInFunction<<"Dimension mismatch for lower bounds"<<exit(FatalError);
        opt.set_lower_bounds(lower);
    }
    if(upper.size()!=0)
    {
        if(static_cast<label>(upper.size())!=n)
            FatalErrorInFunction<<"Dimension mismatch for upper bounds"<<exit(FatalError);
        opt.set_upper_bounds(upper);
    }
}

void Foam::Optimizer::setInitial
(
    const std::vector<scalar>& initial
)
{
    if(static_cast<label>(initial.size())!=n)
        FatalErrorInFunction<<"Dimension mismatch for initial values"<<exit(FatalError);
    x_initial = initial;
}

Foam::scalar Foam::Optimizer::nlopt_objFunc
(
    const std::vector<double> &x,
    std::vector<double> &grad,
    void *my_func_data
)
{
    if(Pstream::master())
        Optimizer::masterDone = false;
    else
        Optimizer::masterDone = true;
    Pstream::scatter(Optimizer::masterDone);
    if(!Pstream::master())
    {
        if(Optimizer::masterDone)
        {
            throw OptimizerStop();
        }
    }
    
    return Optimizer::core_objFunc(x,grad,my_func_data);
}
