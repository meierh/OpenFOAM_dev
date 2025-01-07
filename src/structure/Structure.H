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
#include <chrono>

#include <sys/stat.h>
#include <sys/types.h>

#include <dirent.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include "fvMesh.H"
#include "meshSearch.H"
#include "processorPolyPatch.H"
#include "GeometricField.H"
#include "volFieldsFwd.H"
#include "volMesh.H"
#include "fvPatchField.H"
#include "fvFieldSource.H"
#include <array>
#include "BoundingBox.H"
#include "LinearSystem.H"
#include "OFstream.H"

namespace Foam
{

template<typename A,typename B>
struct stdPairHash
{
    label operator()(const std::pair<A, B> &p) const
    {
        return std::hash<A>{}(p.first)^std::hash<B>{}(p.second);
    }
};

template<typename A>
struct foamPairHash
{
    label operator()(const Pair<A> &p) const
    {
        return std::hash<A>{}(p.first())^std::hash<A>{}(p.second());
    }
};

class Rotation;
class Quaternion
{
    public:
        Quaternion(){};
        Quaternion(scalar x,scalar y,scalar z,scalar w);
        Quaternion(const gsMatrix<scalar>& gsQuaternion);
        
        Quaternion operator*(Quaternion const& q) const;
        Quaternion operator/(Quaternion const& q) const;
        Quaternion operator-(Quaternion const& q) const;
        Quaternion invert() const;
        scalar len() const;
        scalar distanceNorm2(Quaternion const& q) const;
        void normalize();

        inline scalar qx() const {return x;}
        inline scalar qy() const {return y;}
        inline scalar qz() const {return z;}
        inline scalar qw() const {return w;}
        
        scalar operator[](uint index) const;
        scalar& operator[](uint index);
        
    private:
        scalar x;
        scalar y;
        scalar z;
        scalar w;
        
    friend Ostream& operator<<(Ostream& os, Quaternion const& m);
    friend Rotation;
};
Ostream& operator<<(Ostream& os, Quaternion const& q);

class Rotation
{
    public:
        Rotation(){};
        Rotation(vector d1, vector d2, vector d3);
        Rotation(const Quaternion& q);
        Rotation operator-(Rotation const& R) const;
        Rotation operator+(Rotation const& R) const;
        bool operator!=(Rotation const& R) const;
        Rotation operator/(scalar alpha) const;
        scalar distanceNorm2(Rotation const& R) const;
        scalar norm2() const;
        inline const vector& get_d1() const {return T[0];}
        inline const vector& get_d2() const {return T[1];}
        inline const vector& get_d3() const {return T[2];}
        
        static Rotation compute_dRdX(const Quaternion& dqdX, const Quaternion& q);
        static FixedList<FixedList<vector,4>,3> compute_dRdq(const Quaternion& q);
        static Rotation compute_d2RdX(const Quaternion& d2qdX, const Quaternion& dqdX, const Quaternion& q);
        static FixedList<FixedList<vector,4>,3> compute_d2Rdq(const Quaternion& q);
        
    private:
        FixedList<vector,3> T;
        
    friend Ostream& operator<<(Ostream& os, Rotation const& m);
};
Ostream& operator<<(Ostream& os, Rotation const& m);


class Structure
{
    public:
        explicit Structure
        (
            const fvMesh& mesh,
            const Time& runTime
        );
        explicit Structure
        (
            const fvMesh& mesh,
            const std::shared_ptr<IOdictionary> structureDict,
            const Time& runTime
        );
        virtual ~Structure();
        void cleanupActiveRodMesh();
        Structure(const Structure&) = delete;
        
        label getNumberRods() const;
        scalar domainStart (label rodNumber) const;
        scalar domainEnd (label rodNumber) const;
        label getMaxDegree(const ActiveRodMesh::rodCosserat* oneRod) const;
        
        void solveOneStep();
        
        static scalar matrixDistance(const FixedList<vector,3>& m1, const FixedList<vector,3>& m2);
        static scalar vectorDistance(const vector& v1, const vector& v2);
        
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
        
        static gsNurbs<scalar> createNurbs
        (
            scalar domainStart,
            scalar domainEnd,
            uint degree,
            std::vector<scalar> coefficients
        );
        
        static gsNurbs<scalar> createNurbs
        (
            scalar domainStart,
            scalar domainEnd,
            uint degree,
            std::vector<vector> coefficients
        );
        
        static std::vector<scalar> computeUniformKnots
        (
            uint degree,
            uint nbrCoefficients,
            scalar domainStart,
            scalar domainEnd
        );
        
        static void nurbsMinusX
        (
            const gsMatrix<scalar>& x,
            const gsMatrix<scalar>& s,
            const gsNurbs<scalar>& nurbs,
            gsMatrix<scalar>& N_x
        );
                
        static void dNurbsdCoeff
        (
            const gsMatrix<scalar>& s,
            const gsNurbs<scalar>& nurbs,
            List<gsMatrix<scalar>>& dNdC
        );
        
        static void intNurbsMinusX
        (
            const gsMatrix<scalar>& x,
            const gsMatrix<scalar>& s,
            const gsNurbs<scalar>& nurbs,
            vector& intNmX
        );
        
        static void gradIntNurbsMinX
        (
            const gsMatrix<scalar>& x,
            const gsMatrix<scalar>& s,
            const gsNurbs<scalar>& nurbs,
            gsMatrix<scalar>& N_x,
            List<gsMatrix<scalar>>& dNdC,
            gsMatrix<scalar>& gradIntN_x
        );
        
        static void minGradStep
        (
            gsNurbs<scalar>& nurbs,
            const gsMatrix<scalar>& gradient,
            vector stepsize
        );
        
        static void linearLeastSquare
        (
            const List<vector>& points,
            const List<scalar>& parameters,
            const gsNurbs<scalar>& nurbs,
            gsMatrix<scalar>& fittedCoeffs
        );
        
        static void fitNurbsCoeffsToPoints
        (
            const List<vector>& points,
            const gsNurbs<scalar>& nurbs,
            gsMatrix<scalar>& fittedCoeffs,
            scalar epsilon = 1e-4
        );
        
        static void fitNurbsCoeffsToPoints
        (
            const List<vector>& points,
            const List<scalar>& parameters,
            const gsNurbs<scalar>& nurbs,
            gsMatrix<scalar>& fittedCoeffs,
            scalar epsilon = 1e-4
        );
        
    // --------------- Rod setting method ---------------
        void getCurveCoeffs(List<List<vector>>& coeffs) const;
        void setCurveCoeffs(const List<List<vector>>& coeffs);
        scalar getCurveCoeff(label rodNumber, label derivCoeffNumber, label dimension) const;
        void setCurveCoeff(label rodNumber, label derivCoeffNumber, label dimension, scalar value);
    // --------------- ------------------ ---------------
        
    // --------------- Rod movement methods ---------------
    private:
        using timedDef = std::vector<std::unique_ptr<std::pair<gsNurbs<scalar>,scalar>>>;
        using timedRot = std::vector<std::unique_ptr<std::pair<gsNurbs<scalar>,scalar>>>;
        std::pair<timedDef,timedRot> prevState;
        std::pair<timedDef,timedRot> currentState;
    protected:
        const std::pair<gsNurbs<scalar>,scalar>* readPrevRodDeformation(label rodNumber) const;
        const std::pair<gsNurbs<scalar>,scalar>* readPrevRodRotation(label rodNumber) const;
    public:
        void pushBackDeformationState();
        virtual void setDeformation(const List<List<vector>>& deformationCoeffs);
        void setDeformation(label rodNumber, const List<vector>& deformationCoeffs);
        virtual vector evaluateRodVelocity(label rodNumber, scalar parameter, scalar angle, scalar radiusFrac) = 0;
        const gsNurbs<scalar>& getDeformation(label rodNumber) const;
    // ----------------- ------------------------- -----------------
        
    // ----------------- --------------- -----------------
        
    // ----------------- Storage methods -----------------
    private:
        std::map<scalar,std::vector<gsNurbs<scalar>>> storage;
    public:
        virtual void store();
        virtual void setToTime(scalar time);
    // ----------------- --------------- -----------------
        
    // ----------------- Rod evaluation methods -----------------
        static vector rodEval(const ActiveRodMesh::rodCosserat* rod, scalar parameter);
        static vector rodDerivEval(const ActiveRodMesh::rodCosserat* rod, scalar parameter);
        static vector rodDeriv2Eval(const ActiveRodMesh::rodCosserat* rod, scalar parameter);
        
        void rodEval(label rodNumber, scalar parameter, vector& d1, vector& d2, vector& d3, vector& r);
        static void rodEval(const ActiveRodMesh::rodCosserat* rod, scalar parameter, vector& d1, vector& d2, vector& d3, vector& r);
        void rodDerivEval(label rodNumber, scalar parameter, vector& d1, vector& d2, vector& d3, vector& r);
        static void rodDerivEval(const ActiveRodMesh::rodCosserat* rod, scalar parameter, vector& dd1dp, vector& dd2dp, vector& dd3dp, vector& drdp);
        void rodDeriv2Eval(label rodNumber, scalar parameter, vector& d1, vector& d2, vector& d3, vector& r);
        static void rodDeriv2Eval(const ActiveRodMesh::rodCosserat* rod, scalar parameter, vector& d2d1dp, vector& d2d2dp, vector& d2d3dp, vector& d2rdp);
        
        static void rodEval(const ActiveRodMesh::rodCosserat* rod, scalar parameter, vector& r);
        static void rodEval(const gsNurbs<scalar>& curve, const gsNurbs<scalar>& def, scalar parameter, vector& r);
        static void rodEval(const gsNurbs<scalar>& curve, const gsNurbs<scalar>& def, const gsNurbs<scalar>& rot, scalar parameter, vector& d1, vector& d2, vector& d3, vector& r);
        static void rodDerivEval(const gsNurbs<scalar>& curve, const gsNurbs<scalar>& def, const gsNurbs<scalar>& rot, scalar parameter, vector& dd1dp, vector& dd2dp, vector& dd3dp, vector& drdp);
        static void rodDeriv2Eval(const gsNurbs<scalar>& curve, const gsNurbs<scalar>& def, const gsNurbs<scalar>& rot, scalar parameter, vector& d2d1dp, vector& d2d2dp, vector& d2d3dp, vector& d2rdp);
        
        Quaternion m_Rot_Eval(label rodNumber, scalar parameter);
        static Quaternion m_Rot_Eval(const ActiveRodMesh::rodCosserat* rod, scalar parameter);
        const gsNurbs<scalar>& getQuaternions(label rodNumber) const;
        List<Rotation> getRotationCurveCoefs(label rodNumber) const;

        static gsMatrix<scalar> evalNurbs(const gsNurbs<scalar>& nurbs,scalar parameter);
        static gsMatrix<scalar> evalNurbsDeriv(const gsNurbs<scalar>& nurbs,scalar parameter);
        static gsMatrix<scalar> evalNurbsDeriv2(const gsNurbs<scalar>& nurbs,scalar parameter);
    // ----------------- ---------------------- -----------------
        
    // ----------------- Rod coefficient derivation methods -----------------
    protected:
        bool constructedCoeffDerivedData = false;
        void constructCoeffDerivedData();
    private:
        //[rodNumber] -> [coeffNumber] -> [dims] -> nurbs
        std::vector<std::vector<std::vector<gsNurbs<scalar>>>> coeffDerivedCenterline;
        //[rodNumber] -> nurbs
        std::vector<gsNurbs<scalar>> initialRotation;
        //[rodNumber] -> [coeffNumber] -> [dims] -> nurbs
        std::vector<std::vector<std::vector<gsNurbs<scalar>>>> coeffDerivedQuaternions;
    public:
        label numberCurveCoeffs(label rodNumber) const;
        label numberDeformationCoeffs(label rodNumber) const;
        vector get_drdC(label rodNumber, label derivCoeffNumber, label derivDimension, scalar parameter) const;
        Quaternion get_iniQuaternions(label rodNumber, scalar parameter) const;
        Quaternion get_Quaternions(label rodNumber, scalar parameter) const;
        Quaternion get_coeffDerivQuaternions
        (
            label rodNumber,
            label derivCoeffNumber,
            label derivDimension,
            scalar parameter
        ) const;
        const gsNurbs<scalar>& getCoeffDerivedQuaternions
        (
            label rodNumber,
            label derivCoeffNumber,
            label derivDimension
        ) const;
        void rodEvalDerivCoeff
        (
            label rodNumber,
            label derivCoeffNumber,
            label derivDimension,
            scalar parameter,
            vector& d1dC,
            vector& d2dC,
            vector& d3dC,
            vector& rdC
        );
        Quaternion m_Rot_Eval_Deriv(label rodNumber, label derivCoeffNumber, label derivDimension, scalar parameter);
    // ----------------- ---------------------------------- -----------------
        
        void printCurves()
        {
            std::cout<<"--------------------------------------------------------------------------------------"<<std::endl;
            std::cout<<"--------------------------------------------------------------------------------------"<<std::endl;
            std::cout<<"--------------------------------------------------------------------------------------"<<std::endl;
            std::cout<<"--------------------------------------------------------------------------------------"<<std::endl;
            std::cout<<"Rods[0]->m_Curve:"<<Rods[0]->m_Curve<<std::endl;
            std::cout<<"Rods[0]->m_Rot:"<<Rods[0]->m_Rot<<std::endl;
            std::cout<<"--------------------------------------------------------------------------------------"<<std::endl;
            std::cout<<"--------------------------------------------------------------------------------------"<<std::endl;
            std::cout<<"--------------------------------------------------------------------------------------"<<std::endl;
            std::cout<<"--------------------------------------------------------------------------------------"<<std::endl;
        }
        
        const scalar initialMeshSpacing;
        
    private:
        List<std::unordered_map<scalar,std::array<vector,4>>> rodEvalBuffer;
        List<std::unordered_map<scalar,std::array<vector,4>>> rodDerivEvalBuffer;
        List<std::unordered_map<scalar,std::array<vector,4>>> rodDeriv2EvalBuffer;
        
    protected:

    // -------------------------------- ActiveRodMesh Setup --------------------------------
        const Time& runTime;
        const fileName runDirectory;
        const fileName caseName;
        const fileName xmlPath;
        const word name;
        word getXMLPath();
        word xmlFromDict(const IOdictionary& stuctureDict);
        word getName();
        int loadRodsFromXML();
        void createRodScaling();
        void createNurbsStructure();
        void createNurbsBoundary();
        void setSolverOptions();
        //void createDeformationCurve();
        void updateRodCoordinateSystem();
        void checkActiveRodMesh();
        std::unique_ptr<ActiveRodMesh::rodMesh> myMesh;
        std::vector< ActiveRodMesh::rodCrossSection* > Geo;
        std::vector< ActiveRodMesh::rodCosserat* > Rods;
        std::vector< gsNurbsBasis<double>* > BasisRef;
        std::vector<gsNurbs<double>> rodsList;
        const int nR; // Number of Rods
        std::string folder;
        void setupActiveRodMesh();
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
    // -------------------------------- ------------------ --------------------------------

    protected:
        const fvMesh& mesh;
        const std::shared_ptr<IOdictionary> structureDict;
        BoundingBox meshBoundingBox;
        BoundingBox computeMeshBoundingBox();
        
        std::pair<scalar,scalar> minMaxSpan (const cell& thisCell);
        void cellDistances
        (
            const cell& thisCell,
            std::vector<scalar>&
        );
        scalar supportDomainMinSize(const DynamicList<label>& supportDomainCells);
        
        /*********************************Parallel topology***********************************/
    private:
        //Neighbour processes of this process (MPI local)
        std::unordered_set<label> neighbourProcesses;
        //Map of process pairs to comm: {procA,procB} -> comm (MPI global)
        std::unordered_map<Pair<label>,label,foamPairHash<label>> proc_ProcToComm;
        //List of comms ordered (MPI global)
        DynamicList<label> orderedComms;
        //List of existing neighborhood data (MPI local)
        List<bool> myProcExistingData;
        // Generates : (neighbourProcesses,proc_ProcToComm,orderedComms)
        void createParallelTopology();
        bool exchangeSpecific = true;
    public:
        const std::unordered_set<label>& getNeighbourProcesses() const {return neighbourProcesses;}
        inline bool neighbourDataExists(label proc) const {return myProcExistingData[proc];}
        template<typename T> void exchangeBetweenAll(List<T>& globalData) const;
        template<typename T> static void exchangeBetweenTwo(List<T>& globalData, const label comm);
        /*********************************Parallel topology***********************************/
        
        /****************************Global(neighbours) mesh data*****************************/
    private:
        //faceInd -> [{proc,cell}]
        std::unordered_map<label,Pair<label>> patchFaceToNeighCellMap;
        // [proc] -> cell -> neighbours -> {proc,cellInd}
        List<List<List<Pair<label>>>> globalMeshGraph;
        // localCell -> neighbours -> {proc,cellInd}
        List<List<Pair<label>>>* localMeshGraph;        
        // [proc] -> cell -> {centre}
        List<List<vector>> globalMeshCentres;
        // [proc] -> cell -> {magnitude}
        List<List<scalar>> globalMeshMagn;
        // Generates : (patchFaceToNeighCellMap)
        void computePatchFaceToCellMap();
        // Generates:(globalMeshGraph,localMeshGraph,globalMeshCentres,globalMeshMagn), Calls:(computePatchFaceToCellMap)
        void generateMeshGraph();
    public:
        const List<List<List<Pair<label>>>>& getMeshGraph() const {return globalMeshGraph;}
        const List<List<Pair<label>>>& getMeshGraph(label proc) const;
        inline vector getCellCentre(label proc, label cell) const {return globalMeshCentres[proc][cell];}
        inline scalar getCellVolume(label proc, label cell) const {return globalMeshMagn[proc][cell];}
        /****************************Global(neighbours) mesh data*****************************/
        /*
    private:
        std::unordered_map<Pair<label>,std::pair<vector,scalar>,foamPairHash<label>> cellDataBuffer;
    public:
        std::unordered_map<Pair<label>,std::pair<vector,scalar>,foamPairHash<label>>& getCellDataBuffer()
        */
        
        /*************************+++++++****Halo cells************++++++++*******************/
    private:
        const label haloWidth = 4;
    public:
        typedef struct{ label index; vector centre; scalar volume; } CellDescription;
    private:
        // [proc] -> processes' halo cells 
        List<DynamicList<CellDescription>> globalHaloCellList_Sorted;
        // [proc] -> cell -> index
        List<std::unordered_map<label,label>> globalHaloCellToIndexMap;
        // Generates: (globalHaloCellList_Sorted,globalHaloCellToIndexMap), Calls:(generateMeshGraph)
        void generateMeshHaloData(label iterations=1);
        void selfCheckHalos();
    public:
        const DynamicList<CellDescription>& getGlobalHaloCellList_Sorted(label proc) const;
        const std::unordered_map<label,label>& getHaloCellToIndexMap(label proc) const;
        
        //Exchange fields on halo cells
        template<typename T>
        std::unique_ptr<List<List<T>>> broadcastHaloFields(const Field<T>& fieldData);
        /*************************+++++++****Halo cells************++++++++*******************/
        
        /********************************Mesh communication***********+++++*******************/
    public:
        // Calls:(computeHaloData)
        void computeMeshSetup();
        /********************************Mesh communication***********+++++*******************/
    
    public:
        static scalar spacingFromMesh(const fvMesh& mesh, label cellInd=0);
        static void neighbourCells(const fvMesh& mesh, label cellInd, DynamicList<label>& neighbours);
        
        label findCell(vector position, label cellHint=-1) const;
        
        void writeCellSizeField(volScalarField& field) const;
        
        void quaternionCheck();
        void rotationCheck();
        void transformationCheck();
        virtual void parameterGradientCheck();
        virtual void selfCheck();
};
void Barrier(bool stop);



#include "Structure.T.H"

}
#endif
