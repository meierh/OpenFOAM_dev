#include "cutCellPolyMesh.H"

Foam::cutCellPolyMesh::cutCellPolyMesh
(
    const IOobject& io,
    List<std::shared_ptr<Nurbs>> Curves,
    cutStatus state
):
polyMesh(io),
Curves(std::move(Curves)),
MainTree(new KdTree(this->Curves)),
NurbsTrees(List<std::unique_ptr<BsTree>>(this->Curves.size()))
{
    for(int i=0;i<this->Curves.size();i++)
    {
        NurbsTrees[i] = std::move(std::unique_ptr<BsTree>(new BsTree(this->Curves[i])));
    }
    
    const pointField& points = this->points();
    labelList firstOrderNearNurbs;
    for(int i=0;i<points.size();i++)
    {
        std::unique_ptr<labelList> firstOrderNearNurbs = MainTree->nearNurbsCurves(points[i]);
        if(firstOrderNearNurbs->size() == 0)
            continue;
        
        minDistanceToPoint(points[i];
        //Reoganize the 
    }
    
}
