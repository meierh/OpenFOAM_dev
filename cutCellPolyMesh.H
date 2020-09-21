#ifndef CUTCELLPOLYMESH_H
#define CUTCELLPOLYMESH_H

#include "fvCFD.H"
#include <functional>
#include <memory>
#include "BsTree.H"
#include "KdTree.H"
#include "Nurbs.H"
#include <unordered_map>

namespace Foam
{    
class cutCellPolyMesh   :   public polyMesh
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
        // Level-Set Function for Boundary Computation
        std::function<scalar(const vector)> levelSet;
        
        // distance to Surface for all Points
        scalarList pointDist;
        
        //old Cell Volume 
        scalarList oldCellVolume;
                
        //Points to Side (+1,0,-1)
        labelList pointsToSide_;
        
        //Edge to Side (+1,0,-1)
        labelList edgesToSide_;
        
        //face to Side (+1,0,-1)
        labelList facesToSide_;
        
        //cell to Side (+1,0,-1)
        labelList cellsToSide_;
        
        //cutfaces to Side(+1,0,-1)
        labelList   cutFacesToSide_;
        
      //Points and data added by newMeshPoints function
      //Begin
            //Nbr of previous existing points
            label nbrOfPrevPoints;
            // new PointField; previous points + added Points
            pointField newMeshPoints_;
        
            //label of edge for every added Point and -1 for intersection point
            labelList pointToEgde_;
            //label of added Point for every edge or -1 if edge has no added point
            labelList edgeToPoint_;
        
            //label of face for every added Point or zero side point and null for intersection non-cut point
            labelListList pointToFaces_;
            //label of added points for every face or null if no point added
            labelListList faceToPoints_;
        
            //label of cells for every added point and null for every intersection point
            labelListList pointToCells_;
            //label of added point for every cell or null if no added Points
            labelListList cellToPoints_;
      //End
            
      //Edges and data added by addMeshEdges function
            //Begin
            //Nbr of previous existing edges
            label nbrOfPrevEdges;
            // new PointField; previous points + added Points
            edgeList newMeshEdges_;
        
            //label of face for every added edge or zero side edge and null every other edge
            labelListList edgeToFaces_;
            //label of added points for every face or null if no point added
            labelListList faceToEdges_;
        
            //label of cells for every added point and null for every intersection point
            labelListList edgeToCells_;
            //label of added edge or zero side edge for every cell or zero size
            labelListList cellToEdges_;
      //End
            
      //Faces and data added by addMeshFaces function
            //Begin
            //Nbr of previous existing edges
            label nbrOfPrevFaces;
            // new PointField; previous points + added Points
            faceList newMeshFaces_;
        
            //label of cells for every added face or zero side edge and null for every intersection point
            labelListList faceToCells_;
            //label of added point for every cell or null if no added Points
            labelListList cellToFaces_;
      //End
            
      //Split faces and data added by cutOldFaces
            //Begin
            // list of all cut faces replacing old faces
            faceList    cutFaces_;
            // list of all faces containing a list of cut Faces or empty
            labelListList   oldFacesToCutFaces_;
      //End
            
      //Data to construct new mesh
            //Begin
            // starting face for each patch
            labelList patchStarts;
            // nbr of faces for each patch
            labelList   patchSizes;
            
            //insert cell split faces
            // added cell split faces
            faceList addedCutFaces;
            // added split faces owner
            labelList addedCutFaceOwner;
            // added split faces neighbor
            labelList addedCutFaceNeighbor;
            
            //insert interior splitted faces
            // added interior split faces
            faceList splitAndUnsplitFacesInterior;
            // added interior split faces owner
            labelList splitAndUnsplitFaceInteriorOwner;
            // added interior split faces neighbor
            labelList splitAndUnsplitFaceInteriorNeighbor;
            
            //insert boundary splitted faces
            // added boundary split faces
            faceList splitAndUnsplitFacesBoundary;
            // added boundary split faces owner
            labelList splitAndUnsplitFaceBoundaryOwner;
            // added boundary split faces neighbor
            labelList splitAndUnsplitFaceBoundaryNeighbor;
            
            //index of minus side cell of splitted cell or -1 if cell is not split
            labelList oldSplittedCellToNewMinusCell;
            //index of plus side cell of splitted cell or -1 if cell is not split
            labelList oldSplittedCellToNewPlusCell;
      //End
        
        List<std::shared_ptr<Nurbs>> Curves;

        std::unique_ptr<KdTree> MainTree;
        List<std::unique_ptr<BsTree>> NurbsTrees;
        
        std::unordered_map<label,scalar> pointToNurbsDist;
        
    //Public Memberfunctions    
    public:
        //Constructor
        explicit cutCellPolyMesh
        (
            const IOobject& io,
            std::function<scalar(const vector)> levelSet,
            cutStatus state = delNegMesh
        );
        
        explicit cutCellPolyMesh
        (
            const IOobject& io,
            List<std::shared_ptr<Nurbs>> Curves,
            cutStatus state = delNegMesh
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
        
        // Add edges for cutting the cells and create dependency data for creating cut cells
        // Adds to edges_;
        void newMeshEdges();
        
        // Add faces for cutting the cells and create dependency data for creating cut cells
        // Adds to faces_;
        void newMeshFaces();
        
        // Splits all old faces being cut into too
        void cutOldFaces();
        
        //rebuild mesh by using cut data
        void createNewMeshData();
        
        //rebuild mesh deleting negative side
        void createNewMeshData_cutNeg();
        
        // project distance from Nurbs Surface on Points
        void projectNurbsSurface();
        
        // project level set function on Points
        void projectLevelSet();

};
void UnitTest_cutCellPolyMesh(int argc, char *argv[], Time& runtime);
}

#endif
