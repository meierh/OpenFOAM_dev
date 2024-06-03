#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <stdio.h>
#include <iostream>
//#include <direct.h>
#include <random>
#include <chrono>
#include <memory>

#include "activeRodMesh.h"
//#include "rodOpt.h"
#include "rodTools.h"
//#include "rodLatticePBC.h"
#include "rodImport.h"
#include "rodCScircle.h"
//#include "rodCSrectangle.h"
#include "timing.h"
#include "rodMesh.h"
//#include "bvh.h"
#include <unordered_set>
#include <unordered_map>

#include <vector>
#if defined( _MSC_VER )
        #if !defined( _CRT_SECURE_NO_WARNINGS )
                #define _CRT_SECURE_NO_WARNINGS         // This test file is not intended to be secure.
        #endif
#endif

#include "tinyxml2.h"
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <memory>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

#include <dirent.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include "fvCFD.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "dynamicRefineFvMesh.H"
#include <array>

namespace Foam
{
class Structure
{
    public:
        explicit Structure
        (
            dynamicRefineFvMesh& mesh,
            const Time& runTime
        );
        ~Structure();
        Structure(const Structure&) = delete;
        
        label getNumberRods();
        label getMaxDegree(const ActiveRodMesh::rodCosserat* oneRod);
        
        void solveOneStep();
        
        static gsNurbs<scalar> createNurbs
        (
            std::vector<scalar> knots,
            uint degree,
            std::vector<scalar> weights,
            std::vector<scalar> coefficients
        )
        {
            gsKnotVector<scalar> cKnots(knots,degree);
            Info<<"Generated gsKnotVector"<<Foam::endl;
            gsMatrix<scalar> cWeight(weights.size(),1);
            for(label i=0; i<weights.size(); i++)
                cWeight.at(i) = weights[i];
            Info<<"Generated weights"<<Foam::endl;
            gsMatrix<scalar> cCoeff(coefficients.size(),1);
            for(label i=0; i<coefficients.size(); i++)
                cCoeff.at(i) = coefficients[i];
            Info<<"Generated coefficients"<<Foam::endl;
            return gsNurbs<scalar>(cKnots,cWeight,cCoeff);
        }
        
        static gsNurbs<scalar> createNurbs
        (
            std::vector<scalar> knots,
            uint degree,
            std::vector<scalar> weights,
            std::vector<vector> coefficients
        )
        {
            gsKnotVector<scalar> cKnots(knots,degree);
            gsMatrix<scalar> cWeight(weights.size(),1);
            for(label i=0; i<weights.size(); i++)
                cWeight.at(i) = weights[i];
            gsMatrix<scalar> cCoeff(coefficients.size(),3);
            for(int i=0; i<coefficients.size(); i++)
            {
                for(label d=0;d<3;d++)
                {
                    cCoeff(i,d) = coefficients[i][d];
                }
            }            
            return gsNurbs<scalar>(cKnots,cWeight,cCoeff);
        }
        static vector rodDerivEval(const ActiveRodMesh::rodCosserat* rod, scalar parameter);
        static void rodEval(const ActiveRodMesh::rodCosserat* rod, scalar parameter, vector& d1, vector& d2, vector& d3, vector& r);
        static void rodEval(const ActiveRodMesh::rodCosserat* rod, scalar parameter, vector& r);
        
    protected:
        const Time& runTime;
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
            const double	geoR0 = 200;			// Radius [mm]
            const double	geoE0 = 0.006;			// Young's modulus [MPa]
            const double	geoNu = 0.45;			// Poisson's ratio [.]
            const double	geoRho = 1100e-9;		// Density [kg/mm設
            const double	Temp0 = 25.;			// Ambient temperature [蚓]
            
            // * Parameters for loads and boundary conditions
            const int	bcMode = 0;				// 0:Clamp in latDir
            const bool	bcRot = 1;				// Fix rotations on clamped boundaries
            const bool	applyGravity = 0;		// Apply gravity load
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
        /*
        const dimensionedScalar& alpha;
        const volScalarField& T;
        const volScalarField& p;
        const volVectorField& U;
        label IBPatchID;
        const dimensionedScalar nu;
        */
        dynamicRefineFvMesh& mesh;

        
        /*
        void computeIBHeatFlux();
        void computeIBTemperature();
        void assignBoundaryFacesToNurbsCurves();
        template<typename Tensor_Type>
        std::unique_ptr<List<List<Tensor_Type>>> computeDistributedLoad
        (
            const Field<Tensor_Type> immersedBoundaryField
        );
        
        void assignForceOnCurve();        
        std::vector<std::unique_ptr<gismo::gsNurbs<double>>> forceCurveStorage;
        void moveNurbs();
        */
        
        
        
        void updateRodCoordinateSystem();
        
        /*********************************Getter: Collect halo cells***********************************/
        const std::unordered_map<label,DynamicList<label>>& getSelfHaloCellToNeighboursMap()
        {
            return selfHaloCellToNeighboursMap;
        }
        const DynamicList<label>& getSelfHaloList_Sorted() {return selfHaloList_Sorted;}
        const List<DynamicList<label>>& getGlobalSelfHaloList_Sorted() {return globalSelfHaloList_Sorted;}
        const std::unordered_map<label,label>& getSelfHaloList_Sorted_IndexMap() {return selfHaloList_Sorted_IndexMap;}
        /*************************************Collect halo cells***************************************/
        
        
        std::pair<scalar,scalar> minMaxSpan (const cell& thisCell);
        
        void cellDistances
        (
            const cell& thisCell,
            std::vector<scalar>&
        );
        
        scalar supportDomainMinSize(const DynamicList<label>& supportDomainCells);
        
    private:
        /*************************************Collect halo cells***************************************/
        // cell -> [neighbourProcesses]
        std::unordered_map<label,DynamicList<label>> selfHaloCellToNeighboursMap;
        // DynamicList< cell >
        DynamicList<label> selfHaloList_Sorted;
        List<DynamicList<label>> globalSelfHaloList_Sorted;
        // cell -> index
        std::unordered_map<label,label> selfHaloList_Sorted_IndexMap;
        
        /*************************************Collect halo cells***************************************/
        
        typedef struct
        {
            label neighbourProcess;
            label neighbourProcessPatchIndex;
            label neighbourProcessFaceLocalIndex;
            DynamicList<DynamicList<label>> cells;
            DynamicList<DynamicList<vector>> cellsCentre;
            DynamicList<DynamicList<scalar>> cellsMag;
        } singleFaceHalo;
        
        // faceInd -> halo
        std::unordered_map<label,singleFaceHalo> haloData;
                

        
        void collectMeshHalos(label iterations=1);
};
}
#endif
