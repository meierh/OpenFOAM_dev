#ifndef OPTIMIZER_H
#define OPTIMIZER_H
#include <nlopt.hpp>

namespace Foam
{
    
class Optimizer
{
    public:
        explicit Optimizer
        (
            std::function<void(void)> 
        );
        ~Optimizer();
        Optimizer(const Optimizer&) = delete;
        
        scalar run();

    protected:
        label n;
        std::vector<scalar> x_initial;
        nlopt::opt opt;
        
        virtual scalar objectiveFunction
        (
            const std::vector<scalar>& x,
            std::vector<double>& grad,
            void *data
        );
        
        void setBounds(std::vector<scalar> lower, std::vector<scalar> upper);
        
        
        
        virtual
        
        
};

}
#endif
