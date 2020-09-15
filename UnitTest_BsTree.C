#include "Nurbs.H"
#include "BsTree.H"

void Foam::UnitTest_BsTree()
{
    Foam::Info<<"BsTree"<<Foam::endl;    
    std::unique_ptr<scalarList> knots(new scalarList(6));
    (*knots)[0] = 0;    (*knots)[1] = 0;    (*knots)[2] = 0;    (*knots)[3] = 1;
    (*knots)[4] = 1;    (*knots)[5] = 1;    
    //Info<<"Knoten"<<endl;    
    int testdegree = 2;    
    std::unique_ptr<scalarList> weights(new scalarList(3));
    (*weights)[0] = 1;    (*weights)[1] = 1;    (*weights)[2] = 2;
    //Info<<"Gewichte"<<endl;    
    std::unique_ptr<List<vector>> controlPoints(new List<vector>(3));
    (*controlPoints)[0] = vector(1,0,0);    
    (*controlPoints)[1] = vector(1,1,0);    
    (*controlPoints)[2] = vector(0,1,0);
    Info<<"Kontrollpunkte"<<endl;    
    std::unique_ptr<Nurbs> QuarterCircle(new Nurbs(
        std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.2,0.2));
    
    BsTree Tree(std::move(QuarterCircle),10);
    
    vector testPointOutside(0.1,0.1,0);
    
    scalarList resOut = Tree.nearestPoint(testPointOutside);
    Info<<resOut.size()<<endl;
}
