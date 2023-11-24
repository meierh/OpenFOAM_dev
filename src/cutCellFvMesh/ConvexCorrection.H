enum VType {cubePnt=0,cutEdge=1,centerPnt=2,edgePntAdded=3};
class addedPointPos
{
public:
    std::pair<VType,std::uint8_t> distStart;
    std::pair<VType,std::uint8_t> distEnd;
    std::pair<VType,std::uint8_t> measureStart;
    std::pair<VType,std::uint8_t> measureEnd;
};
class corrFaceVertice
{
public:
    VType type;
    std::uint8_t value;
    addedPointPos position;
    corrFaceVertice(VType type,std::uint8_t value,addedPointPos position):type(type),value(value),position(position){}
    corrFaceVertice(VType type,std::uint8_t value):type(type),value(value){}
    corrFaceVertice(){}
};

enum FType {mc33Triangle=-1,added=0,splitOld=1,splitOrig=2,old=3,orig=4};
//added: new Face
//splitOld: old Face splitted
//splitOrig: orig Face that is the same as in the starting mesh splitted
//old: face stays the same
//orig: face stays the same and is the same as in the starting mesh
class corrFace
{
public:
    FType type;
    List<corrFaceVertice> faceData;
    label origFace;
    corrFace(FType type,List<corrFaceVertice> faceData,label origFace):
            type(type),faceData(faceData),origFace(origFace) {}
    corrFace();
};

class corrConcaveCellData
{
public:
    List<List<label>> cells;
    List<corrFace> faces;
};

class markSide
{
public:
    corrConcaveCellData redMarkPlus;
    corrConcaveCellData redMarkNotPlus;
};

List<markSide> convexCorrectionData = List<markSide>(32);

void createConvexCorrectionData()
{
// Case c0
//---

// Case c1
//---

// Case c2
{
    corrFace mc33TriUp =       {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,5},{VType::cutEdge,9}},-2};
    corrFace mc33TriLow =      {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,8},{VType::cutEdge,9}},-2};
    // RedMarkPlus
        corrFace addedRedMark={FType::added,{{VType::cutEdge,7},{VType::cutEdge,9},{VType::cubePnt,1},{VType::edgePntAdded,3}},-1};
        //Upper Cell {3,4,5,6,7}
        corrFace split3Up =   {FType::splitOld,{{VType::cutEdge,7},{VType::cubePnt,7},{VType::cubePnt,3},{VType::edgePntAdded,3}},3};
        corrFace split4Up  =  {FType::splitOrig,{{VType::cubePnt,1},{VType::edgePntAdded,3},{VType::cubePnt,3},{VType::cubePnt,2}},4};
        corrFace old2 =       {FType::orig,{},2};
        corrFace old5 =       {FType::old,{},5};
        corrFace old1 =       {FType::old,{},1};
        //Lower Cell {8,9,10}
        corrFace split3Low  = {FType::splitOld,{{VType::cutEdge,7},{VType::cubePnt,7},{VType::cubePnt,3},{VType::edgePntAdded,3}},3};
        corrFace split4Low  = {FType::splitOrig,{{VType::cutEdge,7},{VType::cubePnt,7},{VType::cubePnt,3},{VType::edgePntAdded,3}},4};
        corrFace old0 =       {FType::old,{},0};
        corrConcaveCellData redMarkSide = {{{0,2,3,4,5,6,7},{1,2,8,9,10}},{mc33TriUp,mc33TriLow,addedRedMark,split3Up,split4Up,old2,old5,old1,split3Low,split4Low,old0}};
        
    // RedMarkNotPlus
        corrFace addedNonRedMark = {FType::added,{{VType::cutEdge,9},{VType::cubePnt,5},{VType::cutEdge,7}},-1};
        //Upper Cell {3,4}
        corrFace split5Up =   {FType::splitOld,{{VType::cubePnt,5},{VType::cutEdge,7},{VType::cutEdge,5}},5};
        old1 =                {FType::old,{},1};
        //Lower Cell {5,6,7}
        corrFace split5Low  = {FType::splitOld,{{VType::cubePnt,5},{VType::cubePnt,4},{VType::cutEdge,7}},5};
        old0 =                {FType::old,{},0};
        corrFace old3 =       {FType::old,{},3};
        corrConcaveCellData nonRedMarkSide = {{{0,2,3,4},{1,2,5,6,7}},{mc33TriUp,mc33TriLow,addedNonRedMark,split5Up,old1,split5Low,old0,old3}};
        
    convexCorrectionData[MC33::Case::c2] = {redMarkSide,nonRedMarkSide};
}
// Case c31
//---

// Case c32
{
    corrFace mc33TriLeft =    {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,6},{VType::cutEdge,8}},-2};
    corrFace mc33TriLeftCent ={FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,6},{VType::cutEdge,8}},-2};
    corrFace mc33TriRightCent={FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,5},{VType::cutEdge,8}},-2};
    corrFace mc33TriRight =   {FType::mc33Triangle,{{VType::cutEdge,4},{VType::cutEdge,5},{VType::cutEdge,8}},-2};
    // RedMarkPlus
        //---
    // RedMarkNotPlus
        corrFace addedNonRedMark = {FType::added,{{VType::cutEdge,10},{VType::cubePnt,2},{VType::cubePnt,0},{VType::cutEdge,8}},-1};
        //Upper Cell {5,6,7,8}
        corrFace old2 =       {FType::old,{},2};
        corrFace old3 =       {FType::old,{},3};
        corrFace old5Up =     {FType::old,{{VType::cubePnt,7}},5};
        corrFace split4Up =   {FType::splitOrig,{{VType::cubePnt,2},{VType::cubePnt,3},{VType::cubePnt,0}},4};
        //Lower Cell {9,10,11,12}
        corrFace old0 =       {FType::old,{},0};
        corrFace old1 =       {FType::old,{},1};
        corrFace old5Low =    {FType::old,{{VType::cubePnt,5}},5};
        corrFace split4Low =  {FType::splitOrig,{{VType::cubePnt,0},{VType::cubePnt,1},{VType::cubePnt,2}},4};
        corrConcaveCellData nonRedMarkSide = {{{0,1,4,5,6,7,8},{2,3,4,9,10,11,12}},
                                             {mc33TriLeft,mc33TriLeftCent,mc33TriRightCent,mc33TriRight,addedNonRedMark,
                                              old2,old3,old5Up,split4Up,old0,old1,old5Low,split4Low}
                                            };
        
    convexCorrectionData[MC33::Case::c32] = {{{{}},{}},nonRedMarkSide};
}
// Case c411
//---

// Case c412
{
    corrFace mc33TriA27 = {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,2},{VType::cutEdge,7}},-2};
    corrFace mc33Tri74A = {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,4},{VType::cutEdge,10}},-2};
    corrFace mc33TriA14 = {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,1},{VType::cutEdge,4}},-2};//2
    corrFace mc33Tri841 = {FType::mc33Triangle,{{VType::cutEdge,8},{VType::cutEdge,4},{VType::cutEdge,1}},-2};
    corrFace mc33Tri821 = {FType::mc33Triangle,{{VType::cutEdge,8},{VType::cutEdge,2},{VType::cutEdge,1}},-2};//4
    corrFace mc33Tri782 = {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,8},{VType::cutEdge,2}},-2};
    // RedMarkPlus
        //---
    // RedMarkNotPlus
        corrFace addedAtEdgeA   ={FType::added,{{VType::cutEdge,10},{VType::cubePnt,6},{VType::cutEdge,7}},-1};//6
        corrFace addedAtEdge2   ={FType::added,{{VType::cutEdge,2} ,{VType::cubePnt,3},{VType::cutEdge,7}},-1};
        corrFace addedAtEdge8   ={FType::added,{{VType::cutEdge,8} ,{VType::cubePnt,0},{VType::cutEdge,2}},-1};//10
        corrFace addedAtEdge4   ={FType::added,{{VType::cutEdge,4} ,{VType::cubePnt,5},{VType::cutEdge,10}},-1};
        corrFace addedAtEdge1to4={FType::added,{{VType::cutEdge,1} ,{VType::cubePnt,1},{VType::cutEdge,4}},-1};//12
        corrFace addedAtEdge1to8={FType::added,{{VType::cutEdge,1} ,{VType::cubePnt,1},{VType::cutEdge,8}},-1};

/*
        //Cell 0 {14,15}
        corrFace split1Part1= {FType::,{},1};
        corrFace split5Part1= {FType::,{},5};
        //Cell 1 {16,17,18}
        corrFace split5Part2= {FType::,{},5};
        corrFace split3Part1= {FType::,{},3};
        corrFace old2 = {FType::,{},2};
        //Cell 2 {19,20}
        corrFace split3Part2= {FType::,{},3};
        corrFace split4Part1= {FType::,{},4};
        //Cell 3 {21,22}
        corrFace split1Part2= {FType::,{},1};
        corrFace split0Part1= {FType::,{},0};
        //Cell 4 {23,24}
        corrFace split0Part3= {FType::,{},0};
        corrFace split4Part2= {FType::,{},4};
        //Cell 5 {25}
        corrFace split0Part2= {FType::,{},0};
        
        corrConcaveCellData nonRedMarkSide = {{{1, 6,11, },
                                               {0, },
                                               {5, },
                                               {2, },
                                               {4, },
                                               {3, }
                                                   
        },
                                             {}
                                            };
        
    convexCorrectionData[MC33::Case::c32] = {{{{}},{}},nonRedMarkSide};
*/
}

// Case c5
//---

// Case c611
//---

// Case c612
//---

// Case c71
//---

// Case c72
//---

// Case c73
//---

// Case c741
//---

// Case c742
//---

// Case c8
//---

// Case c9
//---

// Case c1011
//---

// Case c1012
//---

// Case c102
//---

// Case c11
//---

// Case c1211
//---

// Case c1212
//---

// Case c122
//---

// Case c131
//---

// Case c132
//---

// Case c133
//---

// Case c134
//---

// Case c1351
//---

// Case c1352
//---

// Case c14
//---
};
