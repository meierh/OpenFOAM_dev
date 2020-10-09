#include "KdTree.H"
#include <math.h>
#include <set>

Foam::KdTree::KdTree
(
    List<std::shared_ptr<Nurbs>> Items,
    label maxCurvesPerNode
):
Items(Items),
maxCurvesPerNode(maxCurvesPerNode)
{
    listMinMaxBoxes = List<BoundingBox>(this->Items.size());
    
    for(int i=0;i<this->Items.size();i++)
    {
        listMinMaxBoxes[i] = this->Items[i]->computeBoundingBox();
    }
    //Info<<listMinMaxBoxes[0].Min<<" "<<listMinMaxBoxes[0].Max<<endl;
    
    labelList nurbsCurves(this->Items.size());
    for(int i=0;i<nurbsCurves.size();i++)
    {
        nurbsCurves[i] = i;
    }
    
    _nil =  new Node();
    root = newNode(_nil);
    if(nurbsCurves.size() >= 1)
    {
        constructTree(root,nurbsCurves,0,labelList(0),labelList(0),labelList(0));
    }
}

Foam::KdTree::~KdTree()
{
    recursiveNodeDeleter(root);
    delete _nil;
}

void Foam::KdTree::recursiveNodeDeleter(Node* thisNode)
{
    if(thisNode->left != _nil)
        recursiveNodeDeleter(thisNode->left);
    if(thisNode->right != _nil)
        recursiveNodeDeleter(thisNode->right);
    delete thisNode;     
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

    return newNodeItem;
}

void Foam::KdTree::constructTree
(
    Node* thisNode,
    labelList nurbsCurves,
    label treeHeight,
    labelList firstLevel,
    labelList secondLevel,
    labelList thirdLevel
)
{
    /*
    Info<<"constructTree on height: "<<treeHeight<<" with "
        <<nurbsCurves.size()+firstLevel.size()+secondLevel.size()+thirdLevel.size()
        <<" Curves and split in "<<treeHeight % 3<<endl;
    */
    if(nurbsCurves.size()+firstLevel.size()+secondLevel.size()+thirdLevel.size() < 1)
    {
        FatalErrorInFunction
        << " Node with zero Nurbs inside forbidden!"<<endl
        << abort(FatalError);
    }
    
    BoundingBox MinMaxBox;
    //Info<<"MinMaxBox on height "<<treeHeight<<" is "<<MinMaxBox.Min<<"->"<<MinMaxBox.Max<<endl;
    //Create the Bounding Box for the Subtree
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
    for(int i=0; i<thirdLevel.size();i++)
    {
        BoundingBox temp = listMinMaxBoxes[thirdLevel[i]];
        for(int d=0;d<3;d++)
        {
            if(MinMaxBox.Min[d] > temp.Min[d])
                MinMaxBox.Min[d] = temp.Min[d];
            if(MinMaxBox.Max[d] < temp.Max[d])
                MinMaxBox.Max[d] = temp.Max[d];
        }
    }
    thisNode->MinMaxBox = MinMaxBox;
    //Info<<"MinMaxBox on height "<<treeHeight<<" is "<<MinMaxBox.Min<<"->"<<MinMaxBox.Max<<endl;
    // Bounding Box for Subtree created
    
    label divideDim = treeHeight % 3;
    scalar divideBound = (thisNode->MinMaxBox.Min[divideDim]+thisNode->MinMaxBox.Max[divideDim])/2;
    thisNode->divideDim = divideDim;
    thisNode->divideBound = divideBound;
    
    labelList leftSide(0);
    labelList rightSide(0);
    
    labelList nextFirstLevel(0);
    
    labelList nextSecondLevelLeft(0);
    labelList nextSecondLevelRight(0);

    labelList nextThirdLevelLeft(0);
    labelList nextThirdLevelRight(0);
    
    for(int i=0;i<nurbsCurves.size();i++)
    {
        BoundingBox temp = listMinMaxBoxes[nurbsCurves[i]];
        if(temp.Min[divideDim] < divideBound && temp.Max[divideDim] < divideBound)
        {
            //Info<<nurbsCurves[i]<<" from main to left"<<endl;
            leftSide.append(nurbsCurves[i]);
        }
        else if(temp.Min[divideDim] > divideBound && temp.Max[divideDim] > divideBound)
        {
            //Info<<nurbsCurves[i]<<" from main to right"<<endl;
            rightSide.append(nurbsCurves[i]);
        }
        else
        {
            //Info<<nurbsCurves[i]<<" from main to firstLevel"<<endl;
            nextFirstLevel.append(nurbsCurves[i]);
        }
    }
    for(int i=0;i<firstLevel.size();i++)
    {
        BoundingBox temp = listMinMaxBoxes[firstLevel[i]];
        if(temp.Min[divideDim] < divideBound && temp.Max[divideDim] < divideBound)
        {
            //Info<<firstLevel[i]<<" from firstLevel to secondLevelLeft"<<endl;
            nextSecondLevelLeft.append(firstLevel[i]);
        }
        else if(temp.Min[divideDim] > divideBound && temp.Max[divideDim] > divideBound)
        {
            //Info<<firstLevel[i]<<" from firstLevel to secondLevelRight"<<endl;
            nextSecondLevelRight.append(firstLevel[i]);
        }
        else
        {
            //Info<<firstLevel[i]<<" from firstLevel to secondLevelLeft/Right"<<endl;
            nextSecondLevelLeft.append(firstLevel[i]);
            nextSecondLevelRight.append(firstLevel[i]);
        }
    }
    for(int i=0;i<secondLevel.size();i++)
    {
        BoundingBox temp = listMinMaxBoxes[secondLevel[i]];
        if(temp.Min[divideDim] < divideBound && temp.Max[divideDim] < divideBound)
        {
            //Info<<secondLevel[i]<<" from secondLevel to thirdLevelLeft"<<endl;
            nextThirdLevelLeft.append(secondLevel[i]);
        }
        else if(temp.Min[divideDim] > divideBound && temp.Max[divideDim] > divideBound)
        {
            //Info<<secondLevel[i]<<" from secondLevel to thirdLevelRight"<<endl;
            nextThirdLevelRight.append(secondLevel[i]);
        }
        else
        {
            //Info<<secondLevel[i]<<" from secondLevel to thirdLevelLeft/Right"<<endl;
            nextThirdLevelLeft.append(secondLevel[i]);
            nextThirdLevelRight.append(secondLevel[i]);
        }
    }
    for(int i=0;i<thirdLevel.size();i++)
    {
        //Info<<thirdLevel[i]<<" from thirdLevel to nurbsCurves"<<endl;
        thisNode->nurbsCurves.append(thirdLevel[i]);
    }
    
    if
    (
        leftSide.size()+rightSide.size()+nextFirstLevel.size()+
        nextSecondLevelRight.size()+nextSecondLevelLeft.size()+
        nextThirdLevelRight.size()+nextThirdLevelLeft.size()+
        thisNode->nurbsCurves.size()
        >maxCurvesPerNode
    )
    {
        if(leftSide.size()+nextFirstLevel.size()+nextSecondLevelLeft.size()+
            nextThirdLevelLeft.size() >= 1)
        {
            thisNode->left = newNode(thisNode);
            constructTree(thisNode->left,leftSide,treeHeight+1,
                          nextFirstLevel,nextSecondLevelLeft,nextThirdLevelLeft);
        }
        else
        {
            thisNode->left = _nil;
        }
        if(rightSide.size()+nextFirstLevel.size()+nextSecondLevelRight.size()+
            nextThirdLevelRight.size() >= 1)
        {            
            thisNode->right = newNode(thisNode);
            constructTree(thisNode->right,rightSide,treeHeight+1,
                          nextFirstLevel,nextSecondLevelRight,nextThirdLevelRight);
        }
        else
        {
            thisNode->right = _nil;
        }
    }
    else
    {
        thisNode->nurbsCurves.append(leftSide);
        thisNode->nurbsCurves.append(rightSide);
        thisNode->nurbsCurves.append(nextFirstLevel);
        thisNode->nurbsCurves.append(nextSecondLevelLeft);
        thisNode->nurbsCurves.append(nextSecondLevelRight);
        thisNode->nurbsCurves.append(nextThirdLevelLeft);
        thisNode->nurbsCurves.append(nextThirdLevelRight);
    }
}

std::unique_ptr<labelList> Foam::KdTree::nearNurbsCurves
(
    vector point
) const
{
    std::unique_ptr<labelList> foundNurbsCurves = std::unique_ptr<labelList>(new labelList(0));
    traverseKdTree(root,point,foundNurbsCurves);
    return foundNurbsCurves;
}

void Foam::KdTree::traverseKdTree
(
    Node* currentNode,
    vector& point,
    std::unique_ptr<labelList>& foundNurbsCurves
) const
{
    if(currentNode == _nil || !currentNode->MinMaxBox.isInside(point))
        return;
    
    for(int i=0;i<currentNode->nurbsCurves.size();i++)
        if(listMinMaxBoxes[currentNode->nurbsCurves[i]].isInside(point))
            foundNurbsCurves->append(currentNode->nurbsCurves[i]);
    
    if(point[currentNode->divideDim] < currentNode->divideBound)
    {
        traverseKdTree(currentNode->left,point,foundNurbsCurves);
    }
    else
    {
        traverseKdTree(currentNode->right,point,foundNurbsCurves);
    }
}

labelList Foam::KdTree::printPath(vector point)
{
    Info<<"Path for: "<<point<<endl;
    labelList foundNurbsCurves(0);
    printTraversedKdTree(root,point,foundNurbsCurves);
    return foundNurbsCurves;
}

void Foam::KdTree::printTraversedKdTree
(
    Node* currentNode,
    vector point,
    labelList& foundNurbsCurves
) const
{
    if(currentNode == _nil)
    {
        Info<<"End of path because of nil"<<endl<<endl;
        return;
    }
        
    if(!currentNode->MinMaxBox.isInside(point))
    {
        Info<<"End of path because of MinMaxBox:"<<currentNode->MinMaxBox.Min<<"->"<<currentNode->MinMaxBox.Max<<endl<<endl;
        return;
    }
    
    Info<<"Append Curves of Node (MinMaxBox:"<<currentNode->MinMaxBox.Min<<"->"<<currentNode->MinMaxBox.Max<<"): ";
    for(int i=0;i<currentNode->nurbsCurves.size();i++)
        if(listMinMaxBoxes[currentNode->nurbsCurves[i]].isInside(point))
        {
            Info<<currentNode->nurbsCurves[i]<<" ";
            foundNurbsCurves.append(currentNode->nurbsCurves[i]);
        }
    
    Info<<"continuing dim_"<<currentNode->divideDim<<":"<<currentNode->divideBound;
    if(point[currentNode->divideDim] < currentNode->divideBound)
    {
        Info<<" left"<<endl;
        printTraversedKdTree(currentNode->left,point,foundNurbsCurves);
    }
    else
    {
        Info<<" right"<<endl;
        printTraversedKdTree(currentNode->right,point,foundNurbsCurves);
    }    
}
