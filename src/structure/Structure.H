#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <stdio.h>
#include <iostream>
//#include <direct.h>
#include <random>
#include <chrono>
#include <memory>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "activeRodMesh.h"
#include "rodTools.h"
#include "rodImport.h"
#include "rodCScircle.h"
#include "timing.h"
#include "rodMesh.h"
#pragma GCC diagnostic pop 

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
struct stdPairHash
{
    label operator()(const std::pair<label, label> &p) const
    {
        return std::hash<label>{}(p.first)^std::hash<label>{}(p.second);
    }
};
    
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
        );
        
        static gsNurbs<scalar> createNurbs
        (
            std::vector<scalar> knots,
            uint degree,
            std::vector<scalar> weights,
            std::vector<vector> coefficients
        );

        static vector rodDerivEval(const ActiveRodMesh::rodCosserat* rod, scalar parameter);
        static void rodEval(const ActiveRodMesh::rodCosserat* rod, scalar parameter, vector& d1, vector& d2, vector& d3, vector& r);
        static void rodEval(const ActiveRodMesh::rodCosserat* rod, scalar parameter, vector& r);
        const scalar initialMeshSpacing;
        
        const std::unordered_set<label>& getNeighbourProcesses(){return neighbourhoodProcesses;}
        
        const List<List<Pair<label>>>& getMeshGraph(){return meshGraph;}
        
    protected:
        const Time& runTime;
        const fileName runDirectory;
        const fileName caseName;
        const fileName xmlPath;
        const word name;
        word getXMLPath();
        word getName();
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
        dynamicRefineFvMesh& mesh;

        void updateRodCoordinateSystem();

        std::pair<scalar,scalar> minMaxSpan (const cell& thisCell);
        
        void cellDistances
        (
            const cell& thisCell,
            std::vector<scalar>&
        );
        
        scalar supportDomainMinSize(const DynamicList<label>& supportDomainCells);

    protected:
        /*****************************Getter: Collect halo cells*******************************/
        typedef struct
        {
            label index;
            vector centre;
            scalar volume;
        } CellDescription;
        
    public:
        const DynamicList<CellDescription>& getHaloCellList(label process) const
        {
            return globalHaloCellList_Sorted[process];
        }
        const std::unordered_map<label,label>& getHaloCellToIndexMap(label process) const
        {
            return globalHaloCellToIndexMap[process];
        }
        const std::unordered_map<label,List<DynamicList<std::pair<label,label>>>>& getPatchFaceToCellMap() const
        {
            return patchFaceToCellMap;
        }
        const List<List<Pair<label>>>& getHaloMeshGraph(label process) const
        {
            return globalHaloMeshGraph[process];
        }
        
        template<typename T>
        std::unique_ptr<List<List<T>>> broadcastHaloFields
        (
            const GeometricField<T,fvPatchField,volMesh>& fieldData
        );
        /*********************************Collect halo cells***********************************/
    private:
        /*********************************Collect halo cells***********************************/
        // [proc] -> processes' halo cells 
        List<DynamicList<CellDescription>> globalHaloCellList_Sorted;
        // [proc] -> cell -> index
        List<std::unordered_map<label,label>> globalHaloCellToIndexMap;
        //faceInd -> [{neighbor proc,cell}]
        std::unordered_map<label,List<DynamicList<std::pair<label,label>>>> patchFaceToCellMap;
        // [proc] -> processes' halo cells -> edges -> {proc,cell}
        List<List<List<Pair<label>>>> globalHaloMeshGraph;
        /*********************************Collect halo cells***********************************/
        
        std::unordered_set<label> neighbourhoodProcesses;
        
        List<List<Pair<label>>> meshGraph;
        
        void generateMeshGraph();
        
        void collectMeshHaloData(label iterations=1);
    public:
        scalar initialSpacingFromMesh(const dynamicRefineFvMesh& mesh, label cellInd=0);
};
void Barrier(bool stop);
}
#endif
