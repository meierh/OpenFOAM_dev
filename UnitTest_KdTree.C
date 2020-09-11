#include "KdTree.H"

void Foam::UnitTest_KdTree()
{
    Foam::Info<<"KDTREE"<<Foam::endl;
    
    int testdegree = 2;
    std::unique_ptr<scalarList> knots;
    std::unique_ptr<scalarList> weights;
    std::unique_ptr<List<vector>> controlPoints;
    std::unique_ptr<List<Nurbs*>> items(new List<Nurbs*>());

    //0) X[0,1]Y0Z0
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(0,0,0); (*controlPoints)[1]=vector(1,0,0);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));

    //1) X1Y[0,1]Z0
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(1,0,0); (*controlPoints)[1]=vector(1,1,0);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));
    
    //2) X[1,0]Y1Z0
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(1,1,0); (*controlPoints)[1]=vector(0,1,0);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));

    //3) X0Y[1,0]Z0
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(0,1,0); (*controlPoints)[1]=vector(0,0,0);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));

    //4) X0Y0Z[0,1]
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(0,0,0); (*controlPoints)[1]=vector(0,0,1);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));

    //5) X1Y0Z[0,1]
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(1,0,0); (*controlPoints)[1]=vector(1,0,1);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));

    //6) X1Y1Z[0,1]
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(1,1,0); (*controlPoints)[1]=vector(1,1,1);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));

    //7) X0Y1Z[0,1]
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(0,1,0); (*controlPoints)[1]=vector(0,1,1);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));

    //8) X[0,1]Y1Z1
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(0,0,1); (*controlPoints)[1]=vector(1,0,1);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));

    //9) X1Y[0,1]Z1
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(1,0,1); (*controlPoints)[1]=vector(1,1,1);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));
    
    //10) X[1,0]Y1Z1
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(1,1,1); (*controlPoints)[1]=vector(0,1,1);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));

    //11) X0Y[1,0]Z1
    knots = std::unique_ptr<scalarList>(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    weights = std::unique_ptr<scalarList>(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;    
    controlPoints = std::unique_ptr<List<vector>>(new List<vector>(2));
    (*controlPoints)[0]=vector(0,1,1); (*controlPoints)[1]=vector(0,0,1);    
    items->append(new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.1,0.1));

    KdTree testTree(std::move(items),1);
}
