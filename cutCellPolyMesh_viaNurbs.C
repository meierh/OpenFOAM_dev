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
    Info<<"Prepared all Data"<<endl;
    
    projectNurbsSurface();
}
 
void Foam::cutCellPolyMesh::projectNurbsSurface()
{
    const pointField& points = this->points();
    pointDist = scalarList(points.size());
    
    for(int i=0;i<points.size();i++)
    {
        Info<<"Worked on Point: "<<points[i]<<endl;
        
        std::unique_ptr<labelList> firstOrderNearNurbs = MainTree->nearNurbsCurves(points[i]);
        if(firstOrderNearNurbs->size() == 0)
            continue;
        Info<<"Got list size:"<<firstOrderNearNurbs->size()<<endl;
        
        scalarList distToNurbsSurface(0);
        scalarList paraToNurbsSurface(0);
        for(int k=0;k<firstOrderNearNurbs->size();k++)
        {
            label thisNurbs = (*firstOrderNearNurbs)[k];
            Info<<"Index of nurbs:"<<thisNurbs<<endl;
            scalar thisNodePara = NurbsTrees[thisNurbs]->closestParaOnNurbsToPoint(points[i]);
            if(thisNodePara < this->Curves[thisNurbs]->min_U())
            {
                pointDist[i] = 1;
                continue;
            }
            paraToNurbsSurface.append(thisNodePara);
            distToNurbsSurface.append(this->Curves[thisNurbs]->distanceToNurbsSurface(thisNodePara,points[i]));
        }
        
        scalar minDistToNurbsSurface = std::numeric_limits<scalar>::max();
        for(int k=0;k<distToNurbsSurface.size();k++)
        {
            if(distToNurbsSurface[k] < minDistToNurbsSurface)
                minDistToNurbsSurface = distToNurbsSurface[k];
        }
        pointDist[i] = minDistToNurbsSurface;
        
        Info<<"Finished working on Point: "<<points[i]<<" "<<pointDist[i]<<endl;
    }
    
    Info<<"Final point writing"<<endl;
    for(int i=0;i<points.size();i++)
    {
        Info<<points[i]<<" "<<pointDist[i]<<endl;
    }
}
