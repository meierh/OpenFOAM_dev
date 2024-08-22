#ifndef OPTIMIZER_H
#define OPTIMIZER_H
#include <nlopt.hpp>
#include "fvMesh.H"
#include <stdexcept>

namespace Foam
{
    
class OptimizerStop : public std::runtime_error
{
    public:
        explicit OptimizerStop();
        static bool active();
        static void reset();
        
    private:
        static bool activated;        
        static bool singleton;
};
    
class Optimizer
{
    public:
        explicit Optimizer(nlopt::vfunc obj,label n);
        virtual ~Optimizer()=default;
        Optimizer(const Optimizer&) = delete;
        
        scalar run();
        
        void setBounds(const std::vector<scalar>& lower, const std::vector<scalar>& upper);
        void setInitial(const std::vector<scalar>& initial);

    protected:
        label n;
        nlopt::opt opt;
        std::vector<scalar> x_initial;
        static scalar nlopt_objFunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);
        static nlopt::vfunc core_objFunc;
        
    private:
        static bool singleton;
        static bool masterDone;
};

}
#endif
