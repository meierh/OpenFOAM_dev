#ifndef CUTCELLFVMESH_H
#define CUTCELLFVMESH_H

//#include "fvCFD.H"
#include "dynamicRefineFvMesh.H"
#include <functional>
#include <memory>
#include "BsTree.H"
#include "KdTree.H"
#include "Nurbs1D.H"
#include "NurbsReader.H"
#include <unordered_set>
#include <unordered_map>
#include <list>
#include <algorithm>
#include <chrono>
#include <queue>
#include <cassert>
//#include <sys/types.h>
//#include <dirent.h>
//#include <string.h>
//#include <stdio.h>

#include "motionSolver.H"
#include "points0MotionSolver.H"
#include "displacementMotionSolver.H"
#include "valuePointPatchField.H"
#include "fixedValuePointPatchField.H"
#include "processorPointPatchField.H"
#include "motionDiffusivity.H"
#include "displacementLaplacianFvMotionSolver.H"
#include "patchWave.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"

typedef double GRD_data_type;
typedef double MC33_real;

namespace Foam
{
//class motionSolver;
//class points0MotionSolver;
//class displacementMotionSolver;
    
class cutCellFvMesh   :   public dynamicRefineFvMesh //fvMesh
{
    //Public attributes
    public:
        enum cutStatus 
        {
            delNegMesh = 0,
            internalCut
        };
        
        struct nurbsReference
        {
            label nurbsInd=-1;
            scalar nurbsPara=-1;
        };
        
        class MC33{
            public:
                explicit MC33(cutCellFvMesh& mesh);
                
                typedef struct{
                    label cell=-1;
                    List<label> vertices = List<label>(8,-1);
                    List<edge> edges = List<edge>(12,edge(-1,-1));
                    List<label> edgeGlobalInd = labelList(12,-1);
                    List<label> cutEdgeVerticeIndex = labelList(12,-1);
                    //List<label> faces = List<label>(6);
                    std::vector<std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>> cutTriangles;
                    unsigned int bitPattern;
                } MC33Cube;
                
                MC33Cube computeCutCell(int cellInd);
                std::string getBitPattern(const unsigned short int x);
                
            private:
                cutCellFvMesh& mesh;
                const unsigned short int table[2310];
                
                MC33_real v[8];
                //Procedures
                int face_tests(int *, int) const;
                int face_test1(int) const;
                int interior_test(int, int) const;
          
                void setBit(unsigned short int & bitField, std::uint8_t bitIndex);
                bool  getBit(const unsigned short int  bitField, std::uint8_t bitIndex);
                const unsigned short int* getTriangleCase(unsigned int verticePattern);
                std::vector<std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>> collectTriangles
                (
                    const unsigned short int* triangleCase
                );
                MC33Cube generateMC33Cube(int cellInd);
                unsigned int computeSignBitPattern(const MC33Cube& cell);
                bool cellIsCut(int cellInd);
        };

    //Private Attributes
    private:
        // cellMinDimToStructureDim
        scalar cellDimToStructureDimLimit;
        
        // nurbs intersection radius
        scalar intersectionRadius = 0;
        
        //Motion Solver pointer
        autoPtr<motionSolver> motionPtr_;
        
        // Type of cut immersed boundary algorithm
        cutStatus ibAlgorithm;
        
        // Threshold of cell size difference
        scalar partialThreeshold = 0.25;
        
        // Level-Set Function for Boundary Computation
        std::function<scalar(const vector)> levelSet;
        
        // distance to Surface for all Points
        scalarList pointDist;
        
        // point vector to dist to reduce complexity
        std::unordered_map<label,label> newToOldPointIndMap;
        scalarList pointDistMap;
        
        
        //old Cell Volume 
        scalarList oldCellVolume;
        
        //new Cell Volume
        scalarList newCellVolume;
                
        //Points to Side (+1,0,-1)
        DynamicList<label> pointsToSide_;
        
        //Edge to Side (+1,0,-1)
        DynamicList<label> edgesToSide_;
        
        //face to Side (+1,0,-1)
        DynamicList<label> facesToSide_;
        
        //cell to Side (+1,0,-1)
        labelList cellsToSide_;
        
        //Set of all points on cut cell patch
        std::unordered_set<label> cutPatchPoint;
        
        //cutfaces to Side(+1,0,-1)
        DynamicList<label> cutFacesToSide_;
        
      //Points and data added by newMeshPoints function
      //Begin
            MC33 marchingCubesAlgorithm;
        
            //Cut data for every cell
            DynamicList<MC33::MC33Cube> mc33CutCellData;

            //Nbr of previous existing points
            label nbrOfPrevPoints;
            // new PointField; previous points + added Points
            DynamicList<point> newMeshPointsInFunc;
            pointField newMeshPoints_;
        
            //Reference nurbs curves for zero points
            DynamicList<DynamicList<nurbsReference>> meshPointNurbsReference;
            DynamicList<point> meshPoint;
            
            //label of edge for every added Point and -1 for intersection point
            //labelList pointToEgde_;
            DynamicList<label> pointToEgde_;
            //label of added Point for every edge or -1 if edge has no added point
            //labelList edgeToPoint_;
            DynamicList<label> edgeToPoint_;
            
            //label of face for every added Point or zero side point and null for intersection non-cut point
            //labelListList pointToFaces_;
            //DynamicList<DynamicList<label>> pointToFaces_;
            //label of added points for every face or null if no point added
            //labelListList faceToPoints_;
            //DynamicList<DynamicList<label>> faceToPoints_;
            
            //label of cells for every added point and null for every intersection point
            //labelListList pointToCells_;
            //DynamicList<DynamicList<label>> pointToCells_;
            //label of added point for every cell or null if no added Points
            //labelListList cellToPoints_;
            //DynamicList<DynamicList<label>> cellToPoints_;
      //End
            
      //Edges and data added by addMeshEdges function
            //Begin
            //Nbr of previous existing edges
            label nbrOfPrevEdges;
            // new PointField; previous points + added Points
            DynamicList<edge> newMeshEdges_;
        
            //label of face for every added edge or zero side edge and null every other edge
            DynamicList<DynamicList<label>> edgeToFaces_;
            //labelListList edgeToFaces_;
            //label of added points for every face or null if no point added
            DynamicList<DynamicList<label>> faceToEdges_;
            //labelListList faceToEdges_;
        
            //label of cells for every added point and null for every intersection point
            //DynamicList<DynamicList<label>> edgeToCells_;
            //labelListList edgeToCells_;
            //label of added edge or zero side edge for every cell or zero size
            //DynamicList<DynamicList<label>> cellToEdges_;
            //labelListList cellToEdges_;
            
            //true if face is problematic due to zero points -> 4,3 or 2 zero points
            //List<bool> problematicFace;
            //number of total zero points in problematic face
            //List<label> problematicFacePoints;
            //number of total new zero points in problematic face
            //List<label> problematicFaceNewPoints;
      //End
            
      //Faces and data added by addMeshFaces function
            //Begin
            //Nbr of previous existing edges
            label nbrOfPrevFaces;
            // new PointField; previous points + added Points
            DynamicList<face> newMeshFaces_;
            //faceList newMeshFaces_;
        
            //label of cells for every added face or zero side edge and null for every intersection point
            DynamicList<DynamicList<label>> faceToCells_;
            //labelListList faceToCells_;
            //label of added point for every cell or null if no added Points
            //DynamicList<DynamicList<label>> cellToFaces_;
            DynamicList<DynamicList<label>> cellToFaces_;
      //End
            
      //Split faces and data added by cutOldFaces
            //Begin
            // list of all cut faces replacing old faces
            DynamicList<face>    cutFaces_;
            // list of all faces containing a list of cut Faces or empty
            DynamicList<DynamicList<label>>   oldFacesToCutFaces_;
      //End
            
      //Data to construct new mesh
            //Begin
            // starting face for each patch
            labelList patchStarts;
            // nbr of faces for each patch
            labelList   patchSizes;
            
            //insert cell split faces
            // added cell split faces
            DynamicList<face>  addedCutFaces;
            // added split faces owner
            DynamicList<label> addedCutFacesOwner;
            // added split faces neighbor
            DynamicList<label> addedCutFacesNeighbor;
            // added split faces patch Index
            DynamicList<label> addedCutFacesPatchInd;
            
            //insert interior splitted faces
            // added interior split faces
            DynamicList<face>  splitAndUnsplitFacesInterior;
            // added interior split faces owner
            DynamicList<label> splitAndUnsplitFacesInteriorOwner;
            // added interior split faces neighbor
            DynamicList<label> splitAndUnsplitFacesInteriorNeighbor;
            // added interior split faces patch Index
            DynamicList<label> splitAndUnsplitFacesInteriorPatchInd;
            // added interior split faces old Index
            DynamicList<label> splitAndUnsplitFacesInteriorOldInd;
            // faceInd to sub list index
            //std::unordered_map<label,label> sAUFI_NewToOldMap;
            //std::unordered_multimap<label,label> sAUFI_OldToNewMap;
            
            //insert interior to boundary splitted faces
            // added interior to boundary split faces
            DynamicList<face>  splitAndUnsplitFacesInteriorToBoundary;
            // added interior to boundary split faces owner
            DynamicList<label> splitAndUnsplitFacesInteriorToBoundaryOwner;
            // added interior to boundary split faces neighbor
            DynamicList<label> splitAndUnsplitFacesInteriorToBoundaryNeighbor;
            // added interior to boundary split faces patch Index
            DynamicList<label> splitAndUnsplitFacesInteriorToBoundaryPatchInd;
            // added interior to boundary split faces old Index
            DynamicList<label> splitAndUnsplitFacesInteriorToBoundaryOldInd;
            // faceInd to sub list index
            //std::unordered_map<label,label> sAUFITB_NewToOldMap;
            //std::unordered_multimap<label,label> sAUFITB_OldToNewMap;
            
            //insert boundary splitted faces
            // added boundary split faces
            DynamicList<face>  splitAndUnsplitFacesBoundary;
            // added boundary split faces owner
            DynamicList<label> splitAndUnsplitFacesBoundaryOwner;
            // added boundary split faces neighbor
            DynamicList<label> splitAndUnsplitFacesBoundaryNeighbor;
            // added boundary split faces patch Index
            DynamicList<label> splitAndUnsplitFacesBoundaryPatchInd;
            // added boundary split faces old index
            DynamicList<label> splitAndUnsplitFacesBoundaryOldInd;
            // faceInd to sub list index
            //std::unordered_map<label,label> sAUFB_NewToOldMap;
            //std::unordered_multimap<label,label> sAUFB_OldToNewMap;
            
            //insert interface splitted faces
            // added interface split faces
            DynamicList<face>  splitAndUnsplitFacesInterface;
            // added interface split faces owner
            DynamicList<label> splitAndUnsplitFacesInterfaceOwner;
            // added interface split faces neighbor
            DynamicList<label> splitAndUnsplitFacesInterfaceNeighbor;
            // added interface split faces patch Index
            DynamicList<label> splitAndUnsplitFacesInterfacePatchInd;
            // added interface split faces old Ind
            DynamicList<label>  splitAndUnsplitFacesInterfaceOldInd;
            // faceInd to sub list index
            //std::unordered_map<label,label> sAUFIf_NewToOldMap;
            //std::unordered_multimap<label,label> sAUFIf_OldToNewMap;
            
            //list of indexes of minus side cell of splitted cell or empty if cell is not split
            List<DynamicList<label>> oldSplittedCellToNewMinusCell;
            //list of indexes of plus side cell of splitted cell or empty if cell is not split
            List<DynamicList<label>> oldSplittedCellToNewPlusCell;
            //list of all cell, true for the ones that have been deleted
            List<bool> deletedCell;
            //list of labels to map old cell index to new plus side cut cell index; deleted Cells are set to -1
            DynamicList<DynamicList<label>> mapOldCellsToNewCells;
            //list of labels to map new cell index to old cells
            //labelList mapNewCellsToOldCells;
            
            //mapPolyMesh data
            label nOldPoints;
            labelList pointMap; 
            labelList reversePointMap; 
            List<objectMap> pointsFromPoints; // not written
            pointField preMotionPoints;
            labelListList pointZoneMap; // not written
            
            label nOldFaces;
            labelList faceMap;
            labelList reverseFaceMap;
            List<objectMap> facesFromPoints; // not written
            List<objectMap> facesFromEdges; // not written
            List<objectMap> facesFromFaces; // not written
            labelListList faceZonePointMap; // not written
            labelListList faceZoneFaceMap; // not written
            labelHashSet flipFaceFlux; // not written
            
            label nOldCells;
            labelList cellMap;
            labelList reverseCellMap;
            List<objectMap> cellsFromPoints; // not written
            List<objectMap> cellsFromEdges; // not written
            List<objectMap> cellsFromFaces; // not written
            List<objectMap> cellsFromCells; // not written
            autoPtr<scalarField> oldCellVolumesPtr;
            labelListList cellZoneMap; // not written
            
            labelListList patchPointMap;
            labelList oldPatchStarts;
            labelList oldPatchNMeshPoints;
            
            //mapPolyMesh meshCuttingMap;
        //End
            
        //Data for solidFraction Computation
            //index of faces for all cut cell faces
            labelListList cutCellsMinusAndPlus;
            // uncut cell to minus side cut cell face collection
            labelList oldCellToMinusCutCell;
            //uncut cell to plus side cut cell face collection
            labelList oldCellToPlusCutCell;
        //End
            
            //New mesh data
            pointField new_points;
            faceList new_faces;
            labelList new_owner;
            labelList new_neighbour;
            
        //List<std::shared_ptr<Nurbs>> Curves;
        std::shared_ptr<std::vector<Nurbs1D>> Curves;
        //std::shared_ptr<std::vector<Nurbs1D>> Deformations;

        
        std::unique_ptr<KdTree> MainTree;
        List<std::unique_ptr<BsTree>> NurbsTrees;
        
        //std::unordered_map<label,scalar> pointToNurbsDist;
        
        
         
    //Public Memberfunctions    
    public:
        //Constructor
        explicit cutCellFvMesh
        (
            const IOobject& io,
            std::function<scalar(const vector)> levelSet,
            cutStatus state = delNegMesh
        );
        
        explicit cutCellFvMesh
        (
            const IOobject& io,
            std::shared_ptr<std::vector<Nurbs1D>> Curves,
            cutStatus state = delNegMesh
        );
        
        explicit cutCellFvMesh
        (
            const IOobject& io,
            std::shared_ptr<std::vector<Nurbs1D>> Curves,
            Time& runTime,
            std::unique_ptr<volScalarField>& solidFraction
        );
        
        explicit cutCellFvMesh
        (
            const IOobject& io,
            cutStatus state = delNegMesh,
            scalar cellDimToStructureDimLimit = 0.51
        );
        
        explicit cutCellFvMesh
        (
            const IOobject& io,
            scalar cellDimToStructureDimLimit = 0.51
        );
        
        void cutTheImmersedBoundary();
        void cutTheImmersedBoundary_MC33();
        void refineTheImmersedBoundary();
        void computeSolidFraction_MC33(std::unique_ptr<volScalarField>& solidFraction);

        
        // Print Functions for debugg purposes
        void printAddedPoints();
        void printAddedEdges();
        void printAddedFaces();
        void printCutFaces();
        void printNewMeshData();
        void printMesh();
        void selfTestMesh();
        
        //- Return the motionSolver
        const motionSolver& motion() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
        
        // Set Defomration Curve
        void setInitialDeformationCurve
        (
            std::vector<scalarList>& nurbs_to_knots,
            std::vector<List<vector>>& nurbs_to_controlPoints,
            std::vector<scalarList>& nurbs_to_weights,
            std::vector<label>& nurbs_to_degree
        );
        
        //Compute immersed boundary mesh points movement vector
        void moveTheMesh();
        
        // move the Nurbs
        void moveNurbsCurves(std::vector<List<List<vector>>>& movedDeformationControlPoints);

        //- Write the underlying polyMesh and other data
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool write = true
        ) const;
        
        const DynamicList<DynamicList<nurbsReference>>& getMeshPointNurbsReference() const;
        
        const std::shared_ptr<std::vector<Nurbs1D>> getCurves() const;

    //Private Memberfunctions
    private:
        //Create pointsToSide Data
        void pointsToSide();
        void pointsToSide(const pointField& points);
        void findCutCellPnts();
        
        //Create edgesToSide Data
        void edgesToSide();
        void edgesToSide(const edgeList& edges);
        
        //Create facesToSide Data
        void facesToSide();
        void facesToSide(const faceList& faces);
        
        //Create cellsToSide Data
        void cellsToSide();
        void cellsToSide(const cellList& cells);
        
        //Create Marching Cube data
        void executeMarchingCubes();
        
        // Add points for cutting the cells and create dependency data for creating cut cells
        // Adds to points_.
        void newMeshPoints();
        void newMeshPoints_MC33();
        
        scalar distToNurbs(point pnt,bool& foundFlag,DynamicList<nurbsReference>& reference);
        scalar distToNurbs(point pnt,bool& foundFlag);

        label sideToNurbs(point pnt,bool& foundFlag);
        
        //std::unique_ptr<mapPolyMesh> CutCellPolyMeshMap;
        
        // Add edges for cutting the cells and create dependency data for creating cut cells
        // Adds to edges_;
        void newMeshEdges_MC33();
        
        
        bool equalEdges
        (
            std::unordered_set<label> setA,
            std::unordered_set<label> setB
        );
        /*
        void findCycles
        (
            int i,
            label startPoint,
            label nextPoint,
            label prevPoint,
            DynamicList<label> cyclePath,
            DynamicList<label> cycleEdgePath,
            std::unordered_set<label> coveredPoints,
            std::unordered_set<label> usedEdges,
            std::unordered_map<label,std::unordered_set<label>>& pointGraphData,
            DynamicList<std::pair<label,label>>& pointEdgeComb,
            DynamicList<DynamicList<label>>& closedCyclePoints,
            DynamicList<std::unordered_set<label>>& closedCycleEdges,
            DynamicList<DynamicList<label>>& closedCycleEdgesList
        );
        */
        vector vectorToNurbs
        (
            point pnt,
            bool& foundFlag
        );
        List<vector> vectorsToNurbsOfEdge
        (
            point startPoint,
            point endPoint,
            label nbrOfVectors = 10
        );
        
        List<scalar> distToNursOfEdge
        (
            point startPoint,
            point endPoint,
            label nbrOfPoints = 10
        );

        scalar nearestNurbsIndexPara
        (
            point pnt,
            bool& foundFlag,
            label& nurbsInd,
            scalar& nurbsPara
        );
        
        /*
        void computeClosedFaceFront
        (
            label centerPointInd,
            label faceInd,
            labelList& problematicFacePoints,
            const labelList& cellsInd,
            List<DynamicList<face>>& facesInCellsIn,
            List<DynamicList<std::unordered_set<label>>>& facesMapInCellsIn,
            DynamicList<bool>& closedFaceFrontOutFaceInFace,
            DynamicList<DynamicList<face>>& closedFaceFrontOut,
            DynamicList<DynamicList<std::unordered_set<label>>>& closedFaceFrontMapOut
        );
        */
        
        /*
        void groupEdgesToFaces
        (
            List<DynamicList<DynamicList<label>>>& cellNonconnectedEdges,
            List<DynamicList<DynamicList<DynamicList<label>>>>& cellNonConnectedMultiPoints,
            List<DynamicList<DynamicList<DynamicList<label>>>>& cellNonConnectedMultiEdges,
            List<DynamicList<DynamicList<face>>>& cellNonConnectedMultiFaces,
            List<DynamicList<DynamicList<std::unordered_set<label>>>>& cellNonConnectedMultiPointMap,
            List<DynamicList<DynamicList<std::unordered_set<label>>>>& cellNonConnectedMultiEdgeMap
        );
        */
        
        // Add faces for cutting the cells and create dependency data for creating cut cells
        // Adds to faces_;
        void newMeshFaces_MC33();
        
        // Splits all old faces being cut into too
        void cutOldFaces_MC33();

        //rebuild mesh deleting negative side
        void createNewMeshData_MC33();
        
        // project distance from Nurbs Surface on Points
        void projectNurbsSurface(bool reset=false);
        
        // project level set function on Points
        void projectLevelSet();
        
        // correct the point order in faces to correct the normal dir 
        void correctFaceNormalDir(const pointField& points,faceList& faces,const labelList& owner,const labelList& neighbour);
        List<bool> correctListFaceNormalDir(const pointField& points,const faceList& faces,const labelList& owner,const labelList& neighbour);

        //cell agglomeration
        /*
        void agglomerateSmallCells_cutNeg
        (
            scalarList& newCellVolume,
            scalarList& oldCellVolume,
            scalar partialThreeshold = 0.25
        );
        */
        
        void agglomerateSmallCells_MC33
        (
            scalar partialThreeshold = 0.25
        );
        
        void agglomerateSmallCells_cutNeg_plus
        (
            scalarList& newCellVolume,
            scalarList& oldCellVolume,
            scalar partialThreeshold = 0.25
        );
        
        /*
        labelList searchDown
        (
            scalarListList& possibleMergeFaceArea,
            labelListList& possibleMergeFaces,
            labelListList& possibleMergeCells,
            boolList& oneMergeFaceSufficient,
            boolList& mergeNecessary,
            label count,
            std::unordered_set<label> cellReserved
        );
        */
        
        /*
        labelList searchDown_rec
        (
            DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
            DynamicList<DynamicList<label>>& possibleMergeFaces,
            DynamicList<DynamicList<label>>& possibleMergeCells,
            DynamicList<bool>& oneMergeFaceSufficient,
            DynamicList<bool>& mergeNecessary,
            label count,
            DynamicList<DynamicList<label>>& blockedCells,
            std::unordered_set<label>& cellReserved
        );
        */
        
        /*
        List<DynamicList<label>> searchDown_iter
        (
            DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
            DynamicList<bool>& oneMergeFaceSufficient,
            DynamicList<bool>& mergeNecessary
        );
        */
        
        /*
        List<DynamicList<label>> searchDown_iter_preBlock
        (
            const labelList& owner,
            const labelList& neighbour,
            DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
            DynamicList<bool>& oneMergeFaceSufficient,
            DynamicList<bool>& mergeNecessary
        );
        */
        
        List<DynamicList<label>> assignMergeFaces
        (
            const labelList& owner,
            const labelList& neighbour,
            DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
            DynamicList<bool>& oneMergeFaceSufficient,
            DynamicList<bool>& mergeNecessary,
            DynamicList<DynamicList<scalar>>& possibleMergeCellsPartialSize,
            scalar partialThreeshold
        );
        
        //test data for new Mesh creation
        void testNewMeshData
        (
            const faceList& newFaces,
            const labelList& newFaceOwner,
            const labelList& newFaceNeighbor,
            const labelList& patchStarts,
            const labelList& patchSizes
        );
        
        void testForCellSize
        (
            DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
            DynamicList<bool>& oneMergeFaceSufficient,
            DynamicList<bool>& mergeNecessary,
            labelList& mergeFaceOfCell,
            scalar partialThreeshold
        );
        
        void checkForHexCellsInCutArea();
                
        scalar minNurbsRadius();
        
        void Barrier(bool stop=false);
};
void UnitTest_cutCellFvMesh(int argc, char *argv[], Time& runtime);
}

#endif
