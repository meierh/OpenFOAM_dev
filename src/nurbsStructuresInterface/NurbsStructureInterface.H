#ifndef NURBSSTRUCTURE_H
#define NURBSSTRUCTURE_H

#include <stdio.h>
#include <iostream>
//#include <direct.h>
#include <random>
#include <chrono>
#include <memory>

#include "activeRodMesh.h"
//#include "rodOpt.h"
//#include "rodTools.h"
//#include "rodLatticePBC.h"
#include "rodImport.h"
#include "rodCScircle.h"
//#include "rodCSrectangle.h"
#include "timing.h"
#include "rodMesh.h"
//#include "bvh.h"

//#include <dirent.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "cutCellFvMesh.H"

namespace Foam
{
class NurbsStructureInterface
{
    public:
        explicit NurbsStructureInterface
        (
            Time& runTime,
            const dimensionedScalar& alpha,
            const volScalarField& T,
            const volScalarField& p,
            const volVectorField& U,
            cutCellFvMesh& mesh,
            const dimensionedScalar nu,
            const word& IBpatchName = "cutCell"
        );
        ~NurbsStructureInterface();
        NurbsStructureInterface(const NurbsStructureInterface&) = delete;
        
        void solveOneStep();
        
    private:
        Time& runTime;
        fileName runDirectory;
        fileName caseName;
        fileName xmlPath;
        word name;
        word getXMLPath();
        int loadRodsFromXML();
        void createNurbsStructure();
        void createNurbsBoundary();
        void setSolverOptions();
        void createDeformationCurve();
        std::unique_ptr<ActiveRodMesh::rodMesh> myMesh;
        std::vector< ActiveRodMesh::rodCrossSection* > Geo;
        std::vector< ActiveRodMesh::rodCosserat* > Rods;
        std::vector< gsNurbsBasis<double>* > BasisRef;
        std::vector<gsNurbs<double>> rodsList;
        const int nR; // Number of Rods
        std::string folder;
        
        //******************* Parameters *******************//
            double    latScale;// = 1;			// Scaling factor for XML geometry
            const int	latDir = 2;				// Lattice direction (e.g. for gravity load)
            gsVector<double,3> latSize;
        
            // * Parameters for IGA discretization
            const bool	use_mixed = 1;			// Use mixed method?
            const int	p_sim = 6;				// Refine degree of splines to p
            const int	el_sim = 16;			// Refine knot vector of splines to el elements

            // * Parameters for cross-sections
            const double	geoR0 = 0.8;			// Radius [mm]
            const double	geoE0 = 60;				// Young's modulus [MPa]
            const double	geoNu = 0.45;			// Poisson's ratio [.]
            const double	geoRho = 1100e-9;		// Density [kg/mm設
            const double	Temp0 = 25.;			// Ambient temperature [蚓]
            
            // * Parameters for loads and boundary conditions
            const int	bcMode = 0;				// 0:Clamp in latDir
            const bool	bcRot = 1;				// Fix rotations on clamped boundaries
            const bool	applyGravity = 1;		// Apply gravity load
            // Apply non-zero DBC along latDir
            const double	loadStrain = 0.1;		// ... equivalent to an effective strain of ...
            const int	loadDir = +1;			// ... in +1:tension or -1:compression
            const int	loadSteps = 10;			// ... in ... load steps
                
            // * Parameters for simulation control
            const int	omp_nthreads = 16;		// Use multithreading with OpenMP 
            const int	plot_vtk_geo = 1;		// Write VTK file with initial geometry
            const bool	plot_vtk_solid = 0;		// Write VTK files with volume geometry
            const int	plot_vtk = 15;			// Write VTK each ... load steps
            const int	plot_n = 11;			// No. of points along axis in VTK files
            const int	plot_c = 9;				// No. of points around cross-section in VTK files

            // * Parameters for contact
            const bool	isContact = 0;            
            ActiveRodMesh::ContactOptions cntOpt;
            std::vector<int> pp_rid;
            
            //Solver options
            ActiveRodMesh::SolveOptions solveOpt;


        //******************* Parameters *******************//
        const dimensionedScalar& alpha;
        const volScalarField& T;
        const volScalarField& p;
        const volVectorField& U;
        cutCellFvMesh& mesh;
        label IBPatchID;
        const word& IBpatchName;
        const dimensionedScalar nu;
        const DynamicList<DynamicList<cutCellFvMesh::nurbsReference>>& meshPointNurbsReference;
        const std::shared_ptr<std::vector<Nurbs1D>> Curves;
        std::unordered_multimap<label,label> boundaryPntToFaces;
        List<std::multimap<scalar,label>> nurbsParameterToPnt;
        
        std::vector<std::vector<scalar>> nurbsToLimits;
        std::vector<std::vector<std::vector<label>>> nurbsToLimitsFaces;
        std::vector<std::vector<std::vector<scalar>>> nurbsToLimitsFacesWeights;
        
        void computeIBHeatFlux();
        void computeIBTemperature();
        void assignBoundaryFacesToNurbsCurves();
        template<typename Tensor_Type>
        std::unique_ptr<std::vector<std::vector<Tensor_Type>>> computeDistributedLoad
        (
            const Field<Tensor_Type> immersedBoundaryField
        );
        
        void assignForceOnCurve();        
        std::vector<std::unique_ptr<gismo::gsNurbs<double>>> forceCurveStorage;
        void moveNurbs();        
};
}
#endif
