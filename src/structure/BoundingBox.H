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
        BoundingBox(vector smaller, vector larger):smaller(smaller),larger(larger){}
        BoundingBox operator+(const BoundingBox& rhs) const;
        //Union of BoundingBoxes
        BoundingBox U(const BoundingBox& rhs) const;
        void enlarge(scalar size);
        bool inside(vector point, scalar eps=0) const;
        scalar innerSize();
        FixedList<vector,8> allVertices();
        
        static BoundingBox boundsOfCoefficients(const gsMatrix<scalar>& coefs);
        static BoundingBox boundsOfNurbs(const gsNurbs<scalar>& curve);
        static BoundingBox boundsOfNurbs(gsNurbs<scalar> curve, scalar start, scalar end);
        
    private:
        vector smaller;
        vector larger;
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
};

}

#endif
