#ifndef CUTCELLFVMESH_H
#define CUTCELLFVMESH_H

#include "fvCFD.H"
#include <functional>
#include <memory>
#include "BsTree.H"
#include "KdTree.H"
#include "Nurbs.H"
#include <unordered_set>
#include <unordered_map>
#include <list>
#include <algorithm>
#include <chrono>
#include <queue>

namespace Foam
{    
class cutCellFvMesh   :   public fvMesh
{
    //Public attributes
    public:
        enum cutStatus 
        {
            delNegMesh = 0,
            internalCut
        };


    //Private Attributes
    private:
        // Threshold of cell size difference
        scalar partialThreeshold = 0.25;
        
        // Level-Set Function for Boundary Computation
        std::function<scalar(const vector)> levelSet;
        
        // distance to Surface for all Points
        scalarList pointDist;
        
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
        
        //cutfaces to Side(+1,0,-1)
        DynamicList<label> cutFacesToSide_;
        
      //Points and data added by newMeshPoints function
      //Begin
            //Nbr of previous existing points
            label nbrOfPrevPoints;
            // new PointField; previous points + added Points
            DynamicList<point> newMeshPointsInFunc;
            pointField newMeshPoints_;
        
            //label of edge for every added Point and -1 for intersection point
            //labelList pointToEgde_;
            DynamicList<label> pointToEgde_;
            //label of added Point for every edge or -1 if edge has no added point
            //labelList edgeToPoint_;
            DynamicList<label> edgeToPoint_;
            
            //label of face for every added Point or zero side point and null for intersection non-cut point
            //labelListList pointToFaces_;
            DynamicList<DynamicList<label>> pointToFaces_;
            //label of added points for every face or null if no point added
            //labelListList faceToPoints_;
            DynamicList<DynamicList<label>> faceToPoints_;
            
            //label of cells for every added point and null for every intersection point
            //labelListList pointToCells_;
            DynamicList<DynamicList<label>> pointToCells_;
            //label of added point for every cell or null if no added Points
            //labelListList cellToPoints_;
            DynamicList<DynamicList<label>> cellToPoints_;
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
            DynamicList<DynamicList<label>> edgeToCells_;
            //labelListList edgeToCells_;
            //label of added edge or zero side edge for every cell or zero size
            DynamicList<DynamicList<label>> cellToEdges_;
            //labelListList cellToEdges_;
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
            DynamicList<face> addedCutFaces;
            // added split faces owner
            DynamicList<label> addedCutFaceOwner;
            // added split faces neighbor
            DynamicList<label> addedCutFaceNeighbor;
            
            //insert interior splitted faces
            // added interior split faces
            DynamicList<face> splitAndUnsplitFacesInterior;
            // added interior split faces owner
            DynamicList<label> splitAndUnsplitFaceInteriorOwner;
            // added interior split faces neighbor
            DynamicList<label> splitAndUnsplitFaceInteriorNeighbor;
            
            //insert interior splitted faces
            // added interior split faces
            DynamicList<face> splitAndUnsplitFacesInteriorToBoundary;
            // added interior split faces owner
            DynamicList<label> splitAndUnsplitFaceInteriorToBoundaryOwner;
            // added interior split faces neighbor
            DynamicList<label> splitAndUnsplitFaceInteriorToBoundaryNeighbor;
            
            //insert boundary splitted faces
            // added boundary split faces
            DynamicList<face> splitAndUnsplitFacesBoundary;
            // added boundary split faces owner
            DynamicList<label> splitAndUnsplitFaceBoundaryOwner;
            // added boundary split faces neighbor
            DynamicList<label> splitAndUnsplitFaceBoundaryNeighbor;
            
            //list of indexes of minus side cell of splitted cell or empty if cell is not split
            List<DynamicList<label>> oldSplittedCellToNewMinusCell;
            //list of indexes of plus side cell of splitted cell or empty if cell is not split
            List<DynamicList<label>> oldSplittedCellToNewPlusCell;
            //list of all cells that have been deleted
            labelList deletedCellsList;
            //list of labels to map old cell index to new plus side cut cell index; deleted Cells are set to -1
            DynamicList<DynamicList<label>> mapOldCellsToNewCells;
            //list of labels to map new cell index to old cells
            labelList mapNewCellsToOldCells;
      //End
            
        //Data for solidFraction Computation
            //index of faces for all cut cell faces
            labelListList cutCellsMinusAndPlus;
            // uncut cell to minus side cut cell face collection
            labelList oldCellToMinusCutCell;
            //uncut cell to plus side cut cell face collection
            labelList oldCellToPlusCutCell;
        
        List<std::shared_ptr<Nurbs>> Curves;

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
            List<std::shared_ptr<Nurbs>> Curves,
            cutStatus state = delNegMesh
        );
        
        explicit cutCellFvMesh
        (
            const IOobject& io,
            List<std::shared_ptr<Nurbs>> Curves,
            Time& runTime,
            std::unique_ptr<volScalarField>& solidFraction
        );
        
        // Print Functions for debugg purposes
        void printAddedPoints();
        void printAddedEdges();
        void printAddedFaces();
        void printCutFaces();
        void printNewMeshData();
        void printMesh();
        void selfTestMesh();
    

    //Private Memberfunctions
    private:
        //Create pointsToSide Data
        void pointsToSide();
        void pointsToSide(const pointField& points);
        
        //Create edgesToSide Data
        void edgesToSide();
        void edgesToSide(const edgeList& edges);
        
        //Create facesToSide Data
        void facesToSide();
        void facesToSide(const faceList& faces);
        
        //Create cellsToSide Data
        void cellsToSide();
        void cellsToSide(const cellList& cells);
        
        // Add points for cutting the cells and create dependency data for creating cut cells
        // Adds to points_.
        void newMeshPoints();
        
        scalar distToNurbs(point pnt,bool& foundFlag);
        label sideToNurbs(point pnt,bool& foundFlag);
        
        // Add edges for cutting the cells and create dependency data for creating cut cells
        // Adds to edges_;
        void newMeshEdges();
        bool equalEdges
        (
            std::unordered_set<label> setA,
            std::unordered_set<label> setB
        );
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
        
        // Add faces for cutting the cells and create dependency data for creating cut cells
        // Adds to faces_;
        void newMeshFaces();
        void newMeshFaces_plus();

        
        // Splits all old faces being cut into too
        void cutOldFaces();
        void cutOldFaces_plus();

        
        //rebuild mesh by using cut data
        void createNewMeshData();
        
        //rebuild mesh deleting negative side
        void createNewMeshData_cutNeg();
        void createNewMeshData_cutNeg_plus();
        
        // project distance from Nurbs Surface on Points
        void projectNurbsSurface();
        
        // project level set function on Points
        void projectLevelSet();
        
        // correct the point order in faces to correct the normal dir 
        void correctFaceNormalDir(const pointField& points,faceList& faces,const labelList& owner,const labelList& neighbour);
        
        //cell agglomeration
        void agglomerateSmallCells_cutNeg
        (
            scalarList& newCellVolume,
            scalarList& oldCellVolume,
            scalar partialThreeshold = 0.25
        );
        
        void agglomerateSmallCells_cutNeg_plus
        (
            scalarList& newCellVolume,
            scalarList& oldCellVolume,
            scalar partialThreeshold = 0.25
        );
        
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
        
        List<DynamicList<label>> searchDown_iter
        (
            DynamicList<DynamicList<scalar>>& possibleMergeFaceArea,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeFaces,
            DynamicList<DynamicList<DynamicList<label>>>& possibleMergeCells,
            DynamicList<bool>& oneMergeFaceSufficient,
            DynamicList<bool>& mergeNecessary
        );
        
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
};
void UnitTest_cutCellFvMesh(int argc, char *argv[], Time& runtime);
}

#endif
