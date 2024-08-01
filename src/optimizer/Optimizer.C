#include "Optimizer.H"

scalar Foam::Optimizer::run()
{
    double minf;
    if(x_initial.size()!=n)
        FatalErrorInFunction<<"Dimension mismatch for initial variables"<<exit(FatalError);
    try
    {
        nlopt::result result = opt.optimize(x, minf);
    }
    catch(std::exception &e)
    {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    return minf;
}

void Foam::Optimizer::setBounds
(
    std::vector<scalar> lower,
    std::vector<scalar> upper
)
{
    if(lower.size()==0)
    {
        if(lower.size()!=n)
            FatalErrorInFunction<<"Dimension mismatch for lower bounds"<<exit(FatalError);
        opt.set_lower_bounds(lower);
    }
    if(upper.size()==0)
    {
        if(upper.size()!=n)
            FatalErrorInFunction<<"Dimension mismatch for upper bounds"<<exit(FatalError);
        opt.set_lower_bounds(upper);
    }
}

