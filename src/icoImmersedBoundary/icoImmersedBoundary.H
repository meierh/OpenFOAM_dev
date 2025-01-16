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

#ifndef icoImmersedBoundary_H
#define icoImmersedBoundary_H

#include "incompressibleFluid.H"
#include "CrossSectionStructure.H"
#include "StaticVelocityPressureAction.H"
#include "ForcedMovementVelocityPressureAction.H"
#include "FixedTemperatureAction.H"
#include "MeshRefiner.H"
#include "VelocityMagnitudeLaplacian.H"
#include "RMSVelocityCurvature.H"
#include "VelocityJump.H"

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

#include "pimpleIBControl.H"
#include "setDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class incompressibleFluid Declaration
\*---------------------------------------------------------------------------*/

class icoImmersedBoundary : public incompressibleFluid
{
    public:
        //- Construct from region mesh
        icoImmersedBoundary(fvMesh& mesh, Time& time);
        //- Disallow default bitwise copy construction
        icoImmersedBoundary(const icoImmersedBoundary&) = delete;
        //- Destructor
        virtual ~icoImmersedBoundary() = default;
        //- Disallow default bitwise assignment
        void operator=(const icoImmersedBoundary&) = delete;
        
        void Solve();
        
        void checkIOObjects(){Info<<"checkIOObjects"<<mesh.toc()<<Foam::nl;}
        
        scalar minDeltaT() const;
        
        const std::unique_ptr<LineStructure>& getStructure() const {return structure;}
        std::unique_ptr<LineStructure>& getStructure() {return structure;}


    protected:
        Time& time;
    private:
        pimpleIBControl pimpleCtlr;
    protected:
        IOdictionary transportProperties;
        dimensionedScalar nu;
        
        enum ImmersionMethod {addForce, substrValue};   

        //- RodStructure
        bool useStructure;
        std::shared_ptr<IOdictionary> structureDict;
        std::unique_ptr<LineStructure> structure;
        
        //- Velocity forcing
        bool useVelocityForcing;
        std::unique_ptr<volVectorField> fU_;
        std::unique_ptr<VelocityPressureForceInteraction> interaction_fU;
        scalar velocityForcingFactor=1;
        bool delayedVelocityForcing=false;
        scalar velocityForcingStartTime=0;
        scalar velocityForcingRampUpTime=0;
        void create_VelocityForcing();
        
        //- Temperature
        bool useTemperature;
        std::unique_ptr<volScalarField> T_;
        dimensionedScalar alpha;
        void create_Temperature(pimpleIBControl& pimpleCtlr);
        void create_Temperature();
        
        //- Temperature forcing
        bool useTemperatureForcing;
        std::unique_ptr<volScalarField> fT_;
        std::unique_ptr<TemperatureInteraction> interaction_fT;
        scalar temperatureForcingFactor=1;
        bool delayedTemperatureForcing=false;
        scalar temperatureForcingStartTime=0;
        scalar temperatureForcingRampUpTime=0;
        void create_TemperatureForcing();
        
        //- Refinement
        bool useRefinement;
        std::shared_ptr<MeshRefiner> refinement_;
        std::unique_ptr<volScalarField> refine_;
        void create_Refiner(fvMesh& mesh);
        
        // Check field dimension
        void checkDimensions();
        
        //- Solver Performances
        SolverPerformance<vector> UEqn_res;
        SolverPerformance<scalar> PEqn_res;
        SolverPerformance<scalar> TEqn_res;
        void connectSolverPerformance(pimpleIBControl& pimpleCtlr);
        
        //Storage
        std::map<scalar,fvMesh> storage_Mesh;
        std::map<scalar,volVectorField> storage_U;
        std::map<scalar,volScalarField> storage_p;
        std::map<scalar,volScalarField> storage_T;
        std::vector<scalar> times;
        std::unordered_map<scalar,label> timeToIndex;
        void store();
        void setToTime(scalar time);
    
    // Member Functions
        
        //- Called at the start of the time-step, before the PIMPLE loop
        void preSolve(pimpleIBControl& pimpleCtlr);
        void preSolve()override{FatalErrorInFunction<<"Invalid call preSolve"<<exit(FatalError);};
        
        void preMove(pimpleIBControl& pimpleCtlr);
        //void preMove()override{};
        
        //- Called at the start of the PIMPLE loop to move the mesh
        void moveMesh(pimpleIBControl& pimpleCtlr);
        void moveMesh()override{FatalErrorInFunction<<"Invalid call moveMesh"<<exit(FatalError);};

        //- Corrections that follow mesh motion
        void motionCorrector(pimpleIBControl& pimpleCtlr);
        void motionCorrector()override{FatalErrorInFunction<<"Invalid call motionCorrector"<<exit(FatalError);};

        //- Called at the start of the PIMPLE loop
        void prePredictor(pimpleIBControl& pimpleCtlr);
        void prePredictor()override{FatalErrorInFunction<<"Invalid call prePredictor"<<exit(FatalError);};

        //- Construct and optionally solve the momentum equation
        void momentumPredictor(pimpleIBControl& pimpleCtlr);
        void momentumPredictor()override{FatalErrorInFunction<<"Invalid call momentumPredictor"<<exit(FatalError);};

    /*
        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();
    */
        //- Construct and solve the pressure equation in the PISO loop
        void pressureCorrector(pimpleIBControl& pimpleCtlr);
        void pressureCorrector()override{FatalErrorInFunction<<"Invalid call pressureCorrector"<<exit(FatalError);};
           
        void correctPressure(pimpleIBControl& pimpleCtlr);
        void correctPressure()override{FatalErrorInFunction<<"Invalid call correctPressure"<<exit(FatalError);};
      
        //- Correct the momentum and thermophysical transport modelling
        void postCorrector(pimpleIBControl& pimpleCtlr);
        void postCorrector()override{FatalErrorInFunction<<"Invalid call postCorrector"<<exit(FatalError);};
        
        //- Solve temperature equation after one timestep pimple loop
        void postSolve(pimpleIBControl& pimpleCtlr);
        void postSolve()override{FatalErrorInFunction<<"Invalid call postSolve"<<exit(FatalError);};
        
    // High level functions
        void oneTimestep(pimpleIBControl& pimpleCtlr);
        void oneTimestep();
        
        void create_Analysis();
        void write_Analysis();
            
    private:
        std::unique_ptr<volScalarField> cellSizes;
        std::unique_ptr<volScalarField> cellMarkerCount;
        std::unique_ptr<volScalarField> cellMarkerCharacSize;        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
