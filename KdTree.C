#include "KdTree.H"
#include <math.h> 

Foam::KdTree::KdTree
(
    std::unique_ptr<List<Nurbs*>> Items
)
{
    this->Items = std::unique_ptr<List<Nurbs*>>(std::move(Items));
    
    listMinMaxBoxes = List<BoundingBox>(this->Items->size());
    
    for(int i=0;i<this->Items->size();i++)
    {
        listMinMaxBoxes[i] = (*(this->Items))[i]->computeBoundingBox();
    }
    Info<<listMinMaxBoxes[0].Min<<" "<<listMinMaxBoxes[0].Max<<endl;
    
    labelList nurbsCurves(this->Items->size());
    for(int i=0;i<nurbsCurves.size();i++)
    {
        nurbsCurves[i] = i;
    }
    
    _nil =  new Node();
    root = newNode(nurbsCurves,_nil);
    if(nurbsCurves.size() >= 1)
    {
        constructTree(root,nurbsCurves,0);
    }
}


Foam::KdTree::Node *Foam::KdTree::newNode
(
    Node* parent
)
{

    Node* newNodeItem = new Node();
    
    newNodeItem->left = _nil;
    newNodeItem->right = _nil;
    newNodeItem->parent = parent;
    
    Foam::BoundingBox MinMaxBox = listMinMaxBoxes[nurbsCurves[0]];
    for(int i=0; i<nurbsCurves.size();i++)
    {
        BoundingBox temp = listMinMaxBoxes[nurbsCurves[0]];
        for(int d=0;d<3;d++)
        {
            if(MinMaxBox.Min[d] > temp.Min[d])
                MinMaxBox.Min[d] = temp.Min[d];
            if(MinMaxBox.Max[d] < temp.Max[d])
                MinMaxBox.Max[d] = temp.Max[d];
        }
    }
    newNodeItem->MinMaxBox = MinMaxBox;
    
    return newNodeItem;
}

void Foam::KdTree::constructTree
(
    Node* thisNode,
    labelList nurbsCurves,
    label treeHeight,
    labelList firstLevel,
    labelList secondLevel,
)
{
    if(nurbsCurves.size() <= 1)
    {
        FatalErrorInFunction
        << " Node with zero Nurbs inside forbidden!"<<endl
        << abort(FatalError);
    }
    
    //Create the Bounding Box for the Subtree
    Foam::BoundingBox MinMaxBox = listMinMaxBoxes[thisNode->nurbsCurves[0]];
    for(int i=0; i<nurbsCurves.size();i++)
    {
        BoundingBox temp = listMinMaxBoxes[nurbsCurves[i]];
        for(int d=0;d<3;d++)
        {
            if(MinMaxBox.Min[d] > temp.Min[d])
                MinMaxBox.Min[d] = temp.Min[d];
            if(MinMaxBox.Max[d] < temp.Max[d])
                MinMaxBox.Max[d] = temp.Max[d];
        }
    }
    for(int i=0; i<firstLevel.size();i++)
    {
        BoundingBox temp = listMinMaxBoxes[firstLevel[i]];
        for(int d=0;d<3;d++)
        {
            if(MinMaxBox.Min[d] > temp.Min[d])
                MinMaxBox.Min[d] = temp.Min[d];
            if(MinMaxBox.Max[d] < temp.Max[d])
                MinMaxBox.Max[d] = temp.Max[d];
        }
    }
    for(int i=0; i<secondLevel.size();i++)
    {
        BoundingBox temp = listMinMaxBoxes[secondLevel[i]];
        for(int d=0;d<3;d++)
        {
            if(MinMaxBox.Min[d] > temp.Min[d])
                MinMaxBox.Min[d] = temp.Min[d];
            if(MinMaxBox.Max[d] < temp.Max[d])
                MinMaxBox.Max[d] = temp.Max[d];
        }
    }
    this->MinMaxBox = MinMaxBox;
    // Bounding Box for Subtree created
    
    label divideDim = treeHeight % 3;
    scalar divideBound = (thisNode->MinMaxBox.Min[divideDim]+thisNode->MinMaxBox.Max[divideDim])/2;
    
    labelList leftSide(0);
    labelList rightSide(0);
    labelList bothSides(0);
    
    labelList nextFirstLevel(0);
    labelList nextSecondLevel(0);
    labelList nextThirdLevel(0);
    
    for(int i=0;i<nurbsCurves.size();i++)
    {
        BoundingBox temp = listMinMaxBoxes[nurbsCurves[i]];
        if(temp.Min[divideDim] < divideBound && temp.Max[divideDim] < divideBound)
        {
            leftSide.append(nurbsCurves[i]);
        }
        else if(temp.Min[divideDim] > divideBound && temp.Max[divideDim] > divideBound)
        {
            rightSide.append(nurbsCurves[i]);
        }
        else
        {
            nextFirstLevel.append(nurbsCurves[i]);
        }
    }
    for(int i=0;i<firstLevel.size();i++)
    {
        BoundingBox temp = listMinMaxBoxes[firstLevel[i]];
        if(temp.Min[divideDim] < divideBound && temp.Max[divideDim] < divideBound)
        {
            leftSide.append(firstLevel[i]);
        }
        else if(temp.Min[divideDim] > divideBound && temp.Max[divideDim] > divideBound)
        {
            rightSide.append(firstLevel[i]);
        }
        else
        {
            nextSecondLevel.append(firstLevel[i]);
        }
    }
    for(int i=0;i<secondLevel.size();i++)
    {
        BoundingBox temp = listMinMaxBoxes[secondLevel[i]];
        if(temp.Min[divideDim] < divideBound && temp.Max[divideDim] < divideBound)
        {
            leftSide.append(secondLevel[i]);
        }
        else if(temp.Min[divideDim] > divideBound && temp.Max[divideDim] > divideBound)
        {
            rightSide.append(secondLevel[i]);
        }
        else
        {
            thisNode->nurbsCurves.append(secondLevel[i]);
        }
    }
    
    thisNode->nurbsCurves = bothSides;
    // Limit BoundaryBox to bothSides 

    if
    (
        leftSide.size()+rightSide.size()+nextfirstLevel.size()
        +nextsecondLevel.size()+thisNode->nurbsCurves.size()
        >=maxCurvesPerNode
    )
    {
        if(leftSide.size() >= 1)
        {
            thisNode->left = newNode(thisNode);
            constructTree(thisNode->left,leftSide,treeHeight++);
        }
        if(rightSide.size() >= 1)
        {
            thisNode->right = newNode(thisNode);
            constructTree(thisNode->right,rightSide,treeHeight++);
        }
    }
    else
    {
        thisNode->nurbsCurves.append(leftSide);
        thisNode->nurbsCurves.append(rightSide);
        thisNode->nurbsCurves.append(nextFirstLevel);
        thisNode->nurbsCurves.append(nextSecondLevel);
    }
}
