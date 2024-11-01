#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H

#include <memory>
#include <unordered_set>
#include "fvMesh.H"
#include "activeRodMesh.h"

namespace Foam
{
    
class BoundingBox
{
    public:
        BoundingBox(){}
        BoundingBox(vector smaller, vector larger);
        BoundingBox(const fvMesh& mesh, const cell& cell);
        BoundingBox operator+(const BoundingBox& rhs) const;
        //Compute the largest necessary bounding box around both boxes -- There must be no intersection
        BoundingBox boundingBoxUnion(const BoundingBox& rhs) const;
        //Computes the  Bounding box around both boxes -- There must be no intersection
        BoundingBox boundingBoxIntersection(const BoundingBox& rhs) const;
        scalar distance(const BoundingBox& rhs) const;
        bool boundingBoxOverlap(const BoundingBox& rhs) const;
        void enlarge(scalar size);
        bool inside(vector point, scalar eps=0) const;
        scalar innerSize();
        bool isEmpty(){return empty;}
        FixedList<vector,8> allVertices();
        void print() const {Pout<<smaller<<" -> "<<larger<<Foam::endl;}
        Pair<vector> getBB() const {return {smaller,larger};}
        vector getSmaller() const {return smaller;}
        vector getLarger() const {return larger;}
        
        static BoundingBox boundsOfCoefficients(const gsMatrix<scalar>& coefs);
        static BoundingBox boundsOfNurbs(const gsNurbs<scalar>& curve);
        static BoundingBox boundsOfNurbs(gsNurbs<scalar> curve, scalar start, scalar end);
        
    private:
        vector smaller;
        vector larger;
        bool empty = true;
};
    
template<class Key,class Value>
class BinaryTree
{
    public:
        BinaryTree& operator=(const BinaryTree& rhs);
        const Value* find(const Key& key);
        class Node
        {
            public:
                Node(){}
                Node(const Node& other);
                Key key;
                Value value;
                std::unique_ptr<Node> leftChild;
                std::unique_ptr<Node> rightChild;
        };
        
    protected:
        std::unique_ptr<Node> root;
};

class BoundingBoxTree : public BinaryTree<scalar,BoundingBox>
{
    public:
        BoundingBoxTree& operator=(const BoundingBoxTree& rhs);
        void findPointParameters(std::vector<scalar>& parameters, vector point) const;
        std::unique_ptr<Node>& getRoot(){return root;}
        const BoundingBox& rootBox(){return root->value;}
        void printRoot(){root->value.print();}
};
#include "BoundingBox.T.H"
}

#endif
