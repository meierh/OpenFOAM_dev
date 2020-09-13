#include "KdTree.H"
#include <random>
#include <time.h>

bool listContainsUnique(List<int> list, int number)
{
    bool contains = false;
    for(int i=0;i<list.size();i++)
    {
        if(contains == true && list[i] == number)
        {
            contains = false;
            break;
        }
        if(contains == false && list[i] == number)
            contains = true;
    }
    return contains;
}

double randFrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

Nurbs* generateRandomNurbs(scalar minCoord, scalar maxCoord, scalar diameter, scalar delta_X, int degree)
{
    std::unique_ptr<scalarList> knots(new scalarList(6));
    (*knots)[0]=0; (*knots)[1]=0; (*knots)[2]=0; (*knots)[3]=1; (*knots)[4]=1; (*knots)[5]=1;    
    
    std::unique_ptr<scalarList> weights(new scalarList(2));
    (*weights)[0] = 1;    (*weights)[1] = 1;
    
    std::unique_ptr<List<vector>> controlPoints(new List<vector>(2));
    time_t t;
    srand(static_cast<unsigned>(time(&t)));
    vector start;
    for(int d=0;d<3;d++)
        start[d] = randFrom(minCoord,maxCoord);
    vector end;
    for(int d=0;d<3;d++)
        end[d] = randFrom(minCoord,maxCoord);
    (*controlPoints)[0]=start; (*controlPoints)[1]=end;
    
    return new Nurbs(std::move(knots),std::move(controlPoints),std::move(weights),degree,diameter,delta_X);
}

labelList getNearNurbsBruteForce(const vector point, const List<Foam::BoundingBox>& listMinMaxBoxes)
{
    labelList nearNurbs(0);
    for(int i=0;i<listMinMaxBoxes.size();i++)
    {
        if(listMinMaxBoxes[i].isInside(point))
            nearNurbs.append(i);
    }
    return nearNurbs;
}


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

    //8) X[0,1]Y0Z1
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
    int correctNodes = 0;
    
    KdTree::Node* _nil = testTree._nil;
    KdTree::Node* testTreeRoot = testTree.root;
    if(testTreeRoot->nurbsCurves.size()==0) correctNodes++;
    else    Info<<"Root contains: "<<testTreeRoot->nurbsCurves.size()<<" != 0 curves"<<endl;
    
    KdTree::Node* lvl1_l = testTree.root->left;
    KdTree::Node* lvl1_r = testTree.root->right;
    if(lvl1_l->nurbsCurves.size()==0)   correctNodes++;
    else    Info<<"lvl1_l contains: "<<lvl1_l->nurbsCurves.size()<<" != 0 curves"<<endl;
    if(lvl1_r->nurbsCurves.size()==0)   correctNodes++;
    else    Info<<"lvl1_r contains: "<<lvl1_r->nurbsCurves.size()<<" != 0 curves"<<endl;
    
    KdTree::Node* lvl2_ll = lvl1_l->left;
    KdTree::Node* lvl2_lr = lvl1_l->right;
    KdTree::Node* lvl2_rl = lvl1_r->left;
    KdTree::Node* lvl2_rr = lvl1_r->right;
    if(lvl2_ll->nurbsCurves.size()==0)   correctNodes++;
    else    Info<<"lvl2_ll contains: "<<lvl2_ll->nurbsCurves.size()<<" != 0 curves"<<endl;
    if(lvl2_lr->nurbsCurves.size()==0)   correctNodes++;
    else    Info<<"lvl2_lr contains: "<<lvl2_lr->nurbsCurves.size()<<" != 0 curves"<<endl;
    if(lvl2_rl->nurbsCurves.size()==0)   correctNodes++;
    else    Info<<"lvl2_rl contains: "<<lvl2_rl->nurbsCurves.size()<<" != 0 curves"<<endl;
    if(lvl2_rr->nurbsCurves.size()==0)   correctNodes++;
    else    Info<<"lvl2_rr contains: "<<lvl2_rr->nurbsCurves.size()<<" != 0 curves"<<endl;
    
    KdTree::Node* lvl3_lll = lvl2_ll->left;
    KdTree::Node* lvl3_llr = lvl2_ll->right;
    KdTree::Node* lvl3_lrl = lvl2_lr->left;
    KdTree::Node* lvl3_lrr = lvl2_lr->right;
    KdTree::Node* lvl3_rll = lvl2_rl->left;
    KdTree::Node* lvl3_rlr = lvl2_rl->right;
    KdTree::Node* lvl3_rrl = lvl2_rr->left;
    KdTree::Node* lvl3_rrr = lvl2_rr->right;
    if(lvl3_lll->nurbsCurves.size()==1 && listContainsUnique(lvl3_lll->nurbsCurves,0))   correctNodes++;
    else    Info<<"lvl3_lll contains: "<<lvl3_lll->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl3_llr->nurbsCurves.size()==1 && listContainsUnique(lvl3_llr->nurbsCurves,8))   correctNodes++;
    else    Info<<"lvl3_llr contains: "<<lvl3_llr->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl3_lrl->nurbsCurves.size()==1 && listContainsUnique(lvl3_lrl->nurbsCurves,2))   correctNodes++;
    else    Info<<"lvl3_lrl contains: "<<lvl3_lrl->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl3_lrr->nurbsCurves.size()==1 && listContainsUnique(lvl3_lrr->nurbsCurves,10))   correctNodes++;
    else    Info<<"lvl3_lrr contains: "<<lvl3_lrr->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl3_rll->nurbsCurves.size()==1 && listContainsUnique(lvl3_rll->nurbsCurves,0))   correctNodes++;
    else    Info<<"lvl3_rll contains: "<<lvl3_rll->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl3_rlr->nurbsCurves.size()==1 && listContainsUnique(lvl3_rlr->nurbsCurves,8))   correctNodes++;
    else    Info<<"lvl3_rlr contains: "<<lvl3_rlr->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl3_rrl->nurbsCurves.size()==1 && listContainsUnique(lvl3_rrl->nurbsCurves,2))   correctNodes++;
    else    Info<<"lvl3_rrl contains: "<<lvl3_rrl->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl3_rrr->nurbsCurves.size()==1 && listContainsUnique(lvl3_rrr->nurbsCurves,10))   correctNodes++;
    else    Info<<"lvl3_rrr contains: "<<lvl3_rrr->nurbsCurves.size()<<" != 1 curves"<<endl;
    
    KdTree::Node* lvl4_llll = lvl3_lll->left;
    KdTree::Node* lvl4_lllr = lvl3_lll->right;
    KdTree::Node* lvl4_llrl = lvl3_llr->left;
    KdTree::Node* lvl4_llrr = lvl3_llr->right;
    KdTree::Node* lvl4_lrll = lvl3_lrl->left;
    KdTree::Node* lvl4_lrlr = lvl3_lrl->right;
    KdTree::Node* lvl4_lrrl = lvl3_lrr->left;
    KdTree::Node* lvl4_lrrr = lvl3_lrr->right;
    KdTree::Node* lvl4_rlll = lvl3_rll->left;
    KdTree::Node* lvl4_rllr = lvl3_rll->right;
    KdTree::Node* lvl4_rlrl = lvl3_rlr->left;
    KdTree::Node* lvl4_rlrr = lvl3_rlr->right;
    KdTree::Node* lvl4_rrll = lvl3_rrl->left;
    KdTree::Node* lvl4_rrlr = lvl3_rrl->right;
    KdTree::Node* lvl4_rrrl = lvl3_rrr->left;
    KdTree::Node* lvl4_rrrr = lvl3_rrr->right;
    if(lvl4_llll->nurbsCurves.size()==1 && listContainsUnique(lvl4_llll->nurbsCurves,3))   correctNodes++;
    else    Info<<"lvl4_llll contains: "<<lvl4_llll->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl4_lllr == _nil)   correctNodes++;
    else    Info<<"lvl4_lllr is not _nil"<<endl;
    if(lvl4_llrl->nurbsCurves.size()==1 && listContainsUnique(lvl4_llrl->nurbsCurves,11))   correctNodes++;
    else    Info<<"lvl4_llrl contains: "<<lvl4_llrl->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl4_llrr == _nil)   correctNodes++;
    else    Info<<"lvl4_llrr is not _nil"<<endl;
    if(lvl4_lrll->nurbsCurves.size()==1 && listContainsUnique(lvl4_lrll->nurbsCurves,3))   correctNodes++;
    else    Info<<"lvl4_lrll contains: "<<lvl4_lrll->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl4_lrlr == _nil)   correctNodes++;
    else    Info<<"lvl4_lrlr  is not _nil"<<endl;
    if(lvl4_lrrl->nurbsCurves.size()==1 && listContainsUnique(lvl4_lrrl->nurbsCurves,11))   correctNodes++;
    else    Info<<"lvl4_lrrl contains: "<<lvl4_lrrl->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl4_lrrr == _nil)   correctNodes++;
    else    Info<<"lvl4_lrrr is not _nil"<<endl;
    if(lvl4_rlll == _nil)   correctNodes++;
    else    Info<<"lvl4_rlll contains: "<<lvl3_lll->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl4_rllr->nurbsCurves.size()==1 && listContainsUnique(lvl4_rllr->nurbsCurves,1))   correctNodes++;
    else    Info<<"lvl4_rllr contains: "<<lvl4_rllr->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl4_rlrl == _nil)   correctNodes++;
    else    Info<<"lvl4_rlrl is not _nil"<<endl;
    if(lvl4_rlrr->nurbsCurves.size()==1 && listContainsUnique(lvl4_rlrr->nurbsCurves,9))   correctNodes++;
    else    Info<<"lvl4_rlrr contains: "<<lvl4_rlrr->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl4_rrll == _nil)   correctNodes++;
    else    Info<<"lvl4_rrll is not _nil"<<endl;
    if(lvl4_rrlr->nurbsCurves.size()==1 && listContainsUnique(lvl4_rrlr->nurbsCurves,1))   correctNodes++;
    else    Info<<"lvl4_rrlr contains: "<<lvl4_rrlr->nurbsCurves.size()<<" != 1 curves"<<endl;
    if(lvl4_rrrl == _nil)   correctNodes++;
    else    Info<<"lvl4_rrrl is not _nil"<<endl;
    if(lvl4_rrrr->nurbsCurves.size()==1 && listContainsUnique(lvl4_rrrr->nurbsCurves,9))   correctNodes++;
    else    Info<<"lvl4_rrrr contains: "<<lvl4_rrrr->nurbsCurves.size()<<" != 1 curves"<<endl;
    
    //No testing of lowest level    
    Info<<"UnitTest KdTree Structure Done:"<<correctNodes<<"/31 Nodes correct"<<endl;
    
    
    label NUM_RUNS = 100;
    scalar minCoord = -1;
    scalar maxCoord = 1;
    scalar diameter = 0.1;
    scalar delta_X = 0.1;
    int degree = 2;
    label MAX_NUM_NURBS = 100;
    label NUM_TESTPOINTS_PER_RUN = 1000;
    
    for(int run=0;run<NUM_RUNS;run++)
    {
        std::unique_ptr<List<Nurbs*>> thisRunNurbsCurves(new List<Nurbs*>());
        for(int numNurbs=0;numNurbs<MAX_NUM_NURBS;numNurbs++)
        {
            thisRunNurbsCurves->append(generateRandomNurbs(minCoord, maxCoord, diameter, delta_X, degree));
        }
        KdTree Tree(std::move(thisRunNurbsCurves));
        
        for(int numP=0;numP<NUM_TESTPOINTS_PER_RUN;numP++)
        {
            vector point;
            for(int d=0;d<3;d++)
                point[d] = randFrom(minCoord,maxCoord);
        
            labelList BruteForceRes = getNearNurbsBruteForce(point,Tree.listMinMaxBoxes);            
        }
    }
}
