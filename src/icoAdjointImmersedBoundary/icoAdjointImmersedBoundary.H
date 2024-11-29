/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

Class
    Foam::solvers::icoImmersedBoundary

Description
    Solver module for steady or transient turbulent flow of incompressible
    isothermal fluids with optional mesh motion and change.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, constraining or limiting
    the solution.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    icoImmersedBoundary.C

See also
    Foam::solvers::fluidSolver
    Foam::solvers::isothermalFluid

\*---------------------------------------------------------------------------*/

#ifndef icoAdjointImmersedBoundary_H
#define icoAdjointImmersedBoundary_H

#include "icoImmersedBoundary.H"
#include "SensitivityVelocityPressureForceMarkerInteraction.H"
#include "SensitivityTemperatureInteraction.H"

#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "fvcMeshPhi.H"
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvmLaplacian.H"

#include "icoAdjointVelocityInletBC.H"
#include "icoAdjointVelocityWallBC.H"
#include "icoAdjointVelocityOutletBC.H"
#include "icoAdjointPressureOutletBC.H"
#include "icoAdjointTemperatureOutletBC.H"

#include "pimpleAdjIBControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class incompressibleFluid Declaration
\*---------------------------------------------------------------------------*/

class icoAdjointImmersedBoundary : public icoImmersedBoundary
{
    public:
    // Objective function structure
        typedef struct
        {
            std::function<Field<scalar>(const icoAdjointVelocityInletBC&)> dJdp_Inlet;
            std::function<Field<scalar>(const icoAdjointVelocityWallBC&)> dJdp_Wall;
            std::function<Field<vector>(const icoAdjointVelocityOutletBC&)> dJdu_uOutlet;
            std::function<Field<vector>(const icoAdjointPressureOutletBC&)> dJdu_pOutlet;
            std::function<Field<scalar>(const icoAdjointTemperatureOutletBC&)> dJdT_Outlet;
            std::function<scalar(const icoAdjointImmersedBoundary&)> J;
        } objectiveFunction;
        
    // Constructors
        //- Construct from region mesh
        icoAdjointImmersedBoundary
        (
            fvMesh& mesh,
            Time& time,
            objectiveFunction obj
        );
        //- Disallow default bitwise copy construction
        icoAdjointImmersedBoundary(const icoImmersedBoundary&) = delete;
        //- Destructor
        virtual ~icoAdjointImmersedBoundary()=default;
        //- Disallow default bitwise assignment
        void operator=(const icoAdjointImmersedBoundary&) = delete;
    
    // High level functions
        void SolveSteadyAdjoint();
    
    protected:
        pimpleAdjIBControl adjPimpleCtlr;
        bool steadyStateAdjoint;
        void setupAdjoint();
        
        // Adjoint Veclocity
        volVectorField adj_U_;
        void setAdjUBC
        (
            std::function<Field<scalar>(const icoAdjointVelocityInletBC&)> dJdp_Inlet,
            std::function<Field<scalar>(const icoAdjointVelocityWallBC&)> dJdp_Wall,
            std::function<Field<vector>(const icoAdjointVelocityOutletBC&)> dJdu_uOutlet
        );
        surfaceScalarField adj_phi_;
        
        // Adjoint Pressure
        volScalarField adj_p_;
        void setAdjPBC
        (
            std::function<Field<vector>(const icoAdjointPressureOutletBC&)> dJdu_pOutlet
        );
            
        //- Adjoint Velocity forcing
        bool useAdjointVelocityForcing;
        std::unique_ptr<volVectorField> adj_fU_;
        std::unique_ptr<SensitivityVelocityPressureForceInteraction> interaction_adj_fU;
        void create_AdjointVelocityForcing();
        
        //- Adjoint Temperature
        bool useAdjointTemperature;
        std::unique_ptr<volScalarField> adj_T_;
        void create_AdjointTemperature();
        void setAdjTBC
        (
            std::function<Field<scalar>(const icoAdjointTemperatureOutletBC&)> dJdT_Outlet
        );
        
        //- Adjoint Temperature forcing
        bool useAdjointTemperatureForcing;
        std::unique_ptr<volScalarField> adj_fT_;
        std::unique_ptr<SensitivityTemperatureInteraction> interaction_ajd_fT;
        void create_AdjointTemperatureForcing();
        
        //Set objective function boundary part derivatives
        std::function<scalar(vector n, vector u, scalar p, scalar T)> dJgdp;
        std::function<scalar(vector n, vector u, scalar p, scalar T)> dJgdun;
        std::function<scalar(vector n, vector u, scalar p, scalar T)> dJgdut;
        std::function<scalar(vector n, vector u, scalar p, scalar T)> dJgdT;
        
        // Check adjoint field dimension
        void checkDimensions();

        //Solution process
        bool delayedSolution = false;
        scalar solutionStartTime = 0;

        // Adjoint U equation
        tmp<fvVectorMatrix> tadj_UEqn;
        
        //- Solver Performances
        SolverPerformance<vector> adjUEqn_res;
        SolverPerformance<scalar> adjPEqn_res;
        SolverPerformance<scalar> adjTEqn_res;
        void connectSolverPerformance(pimpleAdjIBControl& pimpleCtlr);
        
        //Gradient Computation
        std::vector<std::pair<Parameter,scalar>> gradient;
        
        //Objective function
        std::function<scalar(const icoAdjointImmersedBoundary&)> J;

    // Member Functions

        //- Load stored sets of 
        void adj_preSolve(pimpleAdjIBControl& adjPimpleCtlr);

        //- Called at the start of the PIMPLE loop to move the mesh
        void adj_moveMesh(){}

        //- Corrections that follow mesh motion
        void adj_motionCorrector(){}

        //- Called at the start of the PIMPLE loop
        void adj_prePredictor(){};
        
        //- Construct and optionally solve the momentum equation
        void adj_momentumPredictor(pimpleAdjIBControl& adjPimpleCtlr);

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        void adj_thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        void adj_pressureCorrector(pimpleAdjIBControl& adjPimpleCtlr);

        //- Correct the momentum and thermophysical transport modelling
        void adj_postCorrector(){};
        
        void adj_correctPressure(pimpleAdjIBControl& adjPimpleCtlr);
        
        //- Computes gradient after one timestep
        void adj_postSolve(pimpleAdjIBControl& adjPimpleCtlr);
    
    // High level functions
        void oneAdjSteadyTimestep(pimpleAdjIBControl& adjPimpleCtlr);
        void oneAdjSteadyTimestep();        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
template<class Type,template<class> class PatchField,class Mesh>
void printBCs(GeometricField<Type,PatchField,Mesh> const& field)
{
    GeometricBoundaryField<Type,PatchField,Mesh> const& boundary = field.boundaryField();
    for(label i=0; i<boundary.size(); i++)
    {
        fvPatchField<Type> const& boundaryField = boundary[i];
        Info<<"-------------------------"<<Foam::nl;
        Info<<"| type: "<<boundaryField.type()<<Foam::nl;
        Info<<"| size: "<<boundaryField.size()<<Foam::nl;
        Info<<"| value: "<<dynamic_cast<List<Type> const&>(boundaryField)<<Foam::nl;
        Type val = Foam::zero();
        for(Type const& v : boundaryField)
            val += v;
        val /= boundaryField.size();
        Info<<"| average value: "<<val<<Foam::nl;
        Info<<"-------------------------"<<Foam::nl;
    }
}
template<class Type,template<class> class PatchField,class Mesh>
void printAvg(GeometricField<Type,PatchField,Mesh> const& field)
{
    Type val = Foam::zero();
    for(label i=0; i<field.size(); i++)
        val += field[i];
    val /= field.size();
    Info<<"Average: "<<val;
}

template<class Type,template<class> class PatchField,class Mesh>
void printMinMax(GeometricField<Type,PatchField,Mesh> const& field);

template<>
void printMinMax<scalar,fvsPatchField,surfaceMesh>(surfaceScalarField const& field)
{
    if(field.size()<1)
        FatalErrorInFunction<<"Error"<<exit(FatalError);
    scalar minVal = field[0];
    scalar maxVal = field[0];
    for(label i=0; i<field.size(); i++)
    {
        minVal = std::min(minVal,field[i]);
        maxVal = std::max(maxVal,field[i]);
    }
    Info<<"  Min: "<<minVal<<" Max: "<<maxVal<<Foam::nl;
}

template<>
void printMinMax<scalar,fvPatchField,volMesh>(volScalarField const& field)
{
    if(field.size()<1)
        FatalErrorInFunction<<"Error"<<exit(FatalError);
    scalar minVal = field[0];
    scalar maxVal = field[0];
    for(label i=0; i<field.size(); i++)
    {
        minVal = std::min(minVal,field[i]);
        maxVal = std::max(maxVal,field[i]);
    }
    Info<<"  Min: "<<minVal<<" Max: "<<maxVal<<Foam::nl;
}

template<>
void printMinMax<vector,fvPatchField,volMesh>(volVectorField const& field)
{
    if(field.size()<1)
        FatalErrorInFunction<<"Error"<<exit(FatalError);
    vector minVal = field[0];
    vector maxVal = field[0];
    for(label i=0; i<field.size(); i++)
    {
        for(label d=0; d<3; d++)
        {
            minVal[d] = std::min(minVal[d],field[i][d]);
            maxVal[d] = std::max(maxVal[d],field[i][d]);
        }
    }
    Info<<"  Min: "<<minVal<<" Max: "<<maxVal<<Foam::nl;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
