/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/
#include "LineStructure.H"
#include <nlopt.hpp>
#include "Optimizer.H"
#include <memory>
#include "fvCFD.H"
#include "dynamicRefineFvMesh.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

double objFunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{    
    if (!grad.empty()) {
        grad[0] = 2*x[0];
        grad[1] = 2*x[1];
    }
    return x[0]*x[0]+x[1]*x[1];
}

class Quadratic : public Optimizer
{
public:
    Quadratic():Optimizer(objFunc,2){}
};

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"    
    pisoControl piso(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"

    nlopt::opt opt(nlopt::LD_MMA, 2);
    std::vector<double> lb(2);
    lb[0] = -HUGE_VAL; lb[1] = 0;
    opt.set_lower_bounds(lb);
    opt.set_min_objective(objFunc, NULL);
    opt.set_xtol_rel(1e-4);
    std::vector<double> x(2);
    x[0] = 1.234; x[1] = 5.678;
    double minf;
    try
    {
        nlopt::result result = opt.optimize(x, minf);
        std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
            << minf << std::endl;
    }
    catch(std::exception &e)
    {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    
    x[0] = 1.234; x[1] = 5.678;
    Quadratic optCl;
    optCl.setBounds(lb,{});
    optCl.setInitial(x);
    optCl.run();
    
    std::vector<bool> cont = {true,true,true,false,false,false,false};
    if(!Pstream::master())
        cont = {true,false,true,false,true,false,true};
    
    DynamicList<word> message;
    
    bool validEnd = false;
    bool masterDone;
    label funcInvoc = 0;
    message.append("Start");
    std::size_t i = 0;
    while(!validEnd)
    {
        bool slaveProcStopped = false;
        try
        {
            for(; i<cont.size(); i++)
            {
                // Pre function
                if(Pstream::master())
                    masterDone = false;
                else
                    masterDone = true;
                Pstream::scatter(masterDone);
                if(!Pstream::master())
                {
                    if(masterDone)
                    {
                        message.append("\t\t\tOptimizerStop invocation");
                        throw OptimizerStop();
                    }
                }
                
                message.append("\t\tFunction invocation: "+std::to_string(i)+"-"+std::to_string(cont[i]));
                funcInvoc++;
            
                //Post function
                if(!cont[i])
                {
                    message.append("\t\t\tRegular stop");
                    i++;
                    break;
                }
            }
        }
        catch(OptimizerStop& stop)
        {
            if(Pstream::master())
                FatalErrorInFunction<<"Master can not be stopped"<<exit(FatalError);
            slaveProcStopped = true;
            message.append("\tOptimizerStop catched");
            break;
        }
        catch(std::exception &e)
        {
            FatalErrorInFunction<<"Optimizer internal fail"<<exit(FatalError);
        }
        
        if(Pstream::master())
        {
            message.append("\tMaster done");
            masterDone = true;
        }
        else
        {
            if(!slaveProcStopped)
                continue;
        }
        
        Pstream::scatter(masterDone);
        if(!masterDone)
            FatalErrorInFunction<<"Master not done can not appear here"<<exit(FatalError);
        if(masterDone)
            break;
    }
    message.append("Completed:"+std::to_string(funcInvoc));
    
    if(Pstream::master())
        Pout<<"message:"<<message<<Foam::endl;
    
    Barrier(false);
    
    if(!Pstream::master())
        Pout<<"message:"<<message<<Foam::endl;
    
    
    return 0;
}


// ************************************************************************* //
