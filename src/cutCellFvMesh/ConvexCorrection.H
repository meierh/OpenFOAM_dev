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
    std::int8_t value;
    addedPointPos position;
    corrFaceVertice(VType type,std::int8_t value,addedPointPos position):type(type),value(value),position(position){}
    corrFaceVertice(VType type,std::int8_t value):type(type),value(value){}
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
        corrFace addedRedMark = {FType::added,{{VType::cutEdge,9},{VType::cubePnt,5},{VType::cutEdge,7}},-1};
        //Upper Cell {3,4}
        corrFace split5Up =   {FType::splitOld,{{VType::cubePnt,5},{VType::cutEdge,7},{VType::cutEdge,5}},5};
        corrFace old1 =                {FType::old,{},1};
        //Lower Cell {5,6,7}
        corrFace split5Low  = {FType::splitOld,{{VType::cubePnt,5},{VType::cubePnt,4},{VType::cutEdge,7}},5};
        corrFace old0 =                {FType::old,{},0};
        corrFace old3 =       {FType::old,{},3};
        corrConcaveCellData redMarkSide = {{{0,2,3,4},{1,2,5,6,7}},{mc33TriUp,mc33TriLow,addedRedMark,split5Up,old1,split5Low,old0,old3}};
        
    // RedMarkNotPlus
        corrFace addedNonRedMark={FType::added,{{VType::cutEdge,7},{VType::cutEdge,9},{VType::cubePnt,1},{VType::edgePntAdded,3}},-1};
        //Upper Cell {3,4,5,6,7}
        corrFace split3Up =   {FType::splitOld,{{VType::cutEdge,7},{VType::cubePnt,7},{VType::cubePnt,3},{VType::edgePntAdded,3}},3};
        corrFace split4Up  =  {FType::splitOrig,{{VType::cubePnt,1},{VType::edgePntAdded,3},{VType::cubePnt,3},{VType::cubePnt,2}},4};
        corrFace old2 =       {FType::orig,{},2};
        corrFace old5 =       {FType::old,{},5};
        old1 =       {FType::old,{},1};
        //Lower Cell {8,9,10}
        corrFace split3Low  = {FType::splitOld,{{VType::cutEdge,7},{VType::cubePnt,7},{VType::cubePnt,3},{VType::edgePntAdded,3}},3};
        corrFace split4Low  = {FType::splitOrig,{{VType::cutEdge,7},{VType::cubePnt,7},{VType::cubePnt,3},{VType::edgePntAdded,3}},4};
        old0 =       {FType::old,{},0};
        corrConcaveCellData nonRedMarkSide = {{{0,2,3,4,5,6,7},{1,2,8,9,10}},{mc33TriUp,mc33TriLow,addedNonRedMark,split3Up,split4Up,old2,old5,old1,split3Low,split4Low,old0}};
        
    convexCorrectionData[MC33::Case::c2] = {redMarkSide,nonRedMarkSide};
}
// Case c31
//---

// Case c32
{
    corrFace mc33TriLeft =    {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,6},{VType::cutEdge,8}},-2};
    corrFace mc33TriLeftCent ={FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,6},{VType::cutEdge,8}},-2};
    corrFace mc33TriRightCent={FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,4},{VType::cutEdge,8}},-2};
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

// Case c412 Further work necessary
{
    corrFace mc33TriA27 = {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,2},{VType::cutEdge,7}},-2};
    corrFace mc33Tri74A = {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,4},{VType::cutEdge,10}},-2};
    corrFace mc33TriA14 = {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,1},{VType::cutEdge,4}},-2};
    corrFace mc33Tri841 = {FType::mc33Triangle,{{VType::cutEdge,8},{VType::cutEdge,4},{VType::cutEdge,1}},-2};
    corrFace mc33Tri821 = {FType::mc33Triangle,{{VType::cutEdge,8},{VType::cutEdge,2},{VType::cutEdge,1}},-2};//4
    corrFace mc33Tri782 = {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,8},{VType::cutEdge,2}},-2};
    // RedMarkPlus
        //---
    // RedMarkNotPlus
        corrFace addedAtEdgesA7  = {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,7} ,{VType::cubePnt,6}},-1};
        corrFace addedAtEdges27  = {FType::mc33Triangle,{{VType::cutEdge,2} ,{VType::cutEdge,7} ,{VType::cubePnt,3}},-1};
        corrFace addedAtEdges28  = {FType::mc33Triangle,{{VType::cutEdge,2}, {VType::cutEdge,8} ,{VType::cubePnt,0}},-1};//8
        corrFace addedAtEdges81  = {FType::mc33Triangle,{{VType::cutEdge,8} ,{VType::cutEdge,1} ,{VType::cubePnt,1}},-1};
        corrFace addedAtEdges41  = {FType::mc33Triangle,{{VType::cutEdge,4} ,{VType::cutEdge,1} ,{VType::cubePnt,1}},-1};
        corrFace addedAtEdgesA4  = {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,4} ,{VType::cubePnt,5}},-1};

        //Cell Top {12,13,14}
        corrFace old2 =       {FType::old,{},2}
        corrFace split3Top  = {FType::splitOld,{{VType::cubePnt,3},{VType::cubePnt,7},{VType::cutEdge,7}},3};
        corrFace split5Top  = {FType::splitOld,{{VType::cubePnt,6},{VType::cubePnt,7},{VType::cutEdge,7}},5};
        //Cell Back 1 Down {15,16}
        corrFace split3Down = {FType::splitOld,{{VType::cubePnt,3},{VType::cutEdge,7},{VType::cutEdge,8},{VType::cubePnt,0}},3};
        corrFace split4Left = {FType::splitOld,{{VType::cubePnt,0},{VType::cubePnt,3},{VType::cutEdge,2}},4};
        //Cell Front 1 Down {17,18}
        corrFace split5Down = {FType::splitOld,{{VType::cubePnt,4},{VType::cubePnt,5},{VType::cubePnt,6},{VType::cutEdge,7}},5};
        corrFace split1Front= {FType::splitOld,{{VType::cutEdge,10},{VType::cubePnt,6},{VType::cubePnt,5}},1};
        //Cell Back 2 Down {19,20}
        corrFace split4Right = {FType::splitOld,{{VType::cubePnt,0},{VType::cubePnt,1},{VType::cutEdge,1},{VType::cutEdge,2}},4};
        corrFace split0Back= {FType::splitOld,{{VType::cutEdge,8},{VType::cubePnt,1},{VType::cubePnt,0}},0};
        //Cell Front 2 Down {21,22}
        corrFace split1Back = {FType::splitOld,{{VType::cubePnt,5},{VType::cubePnt,1},{VType::cutEdge,1},{VType::cutEdge,10}},1};
        corrFace split0Front= {FType::splitOld,{{VType::cutEdge,4},{VType::cubePnt,5},{VType::cubePnt,1}},0};
        //Cell 5 {23}
        corrFace split0Middle={FType::splitOld,{{VType::cutEdge,8},{VType::cutEdge,4},{VType::cubePnt,1}},0};
        
        corrConcaveCellData nonRedMarkSide = {{{0, 6,7,  12,13,14},
                                               {5, 7,8,  15,16},
                                               {1, 6,11, 17,18},
                                               {4, 8,9,  19,20},
                                               {2,10,11, 21,22},
                                               {3, 9,10, 23}},
                                             {mc33TriA27,mc33Tri74A,mc33TriA14,mc33Tri841,mc33Tri821,mc33Tri782,
                                              addedAtEdgesA7,addedAtEdges27,addedAtEdges28,addedAtEdges81,addedAtEdges41,addedAtEdgesA4,
                                              old2,split3Top,split5Top,split3Down,split4Left,split5Down,split1Front,split4Right,split0Back,split1Back,split0Front,split0Middle
                                            }
                                            };
        
    convexCorrectionData[MC33::Case::c32] = {{{{}},{}},nonRedMarkSide};
*/
}
// Case c5
{
    corrFace mc33TriE513 =  {FType::mc33Triangle,{{VType::cutEdge,5},{VType::cutEdge,1},{VType::cutEdge,3}},-2};
    corrFace mc33TriE538 =  {FType::mc33Triangle,{{VType::cutEdge,5},{VType::cutEdge,3},{VType::cutEdge,8}},-2};
    corrFace mc33TriE845 =  {FType::mc33Triangle,{{VType::cutEdge,8},{VType::cutEdge,4},{VType::cutEdge,5}},-2};
    // RedMarkPlus
        corrFace addedRedMark = {FType::added,{{VType::cutEdge,5},{VType::cubePnt,5},{VType::cubePnt,0},{VType::cutEdge,3}},-1};
        //Right Cell {4,5,6}
        corrFace old1 =      {FType::old,{},1};
        corrFace old4 =      {FType::old,{},4};
        corrFace split0Right={FType::splitOld,{{VType::cubePnt,5},{VType::cubePnt,1},{VType::cubePnt,0}},0};
        //Left Cell {7,8,9}
        corrFace old5 =      {FType::old,{},5};
        corrFace old3 =      {FType::old,{},3};
        corrFace split0Left ={FType::splitOld,{{VType::cubePnt,0},{VType::cubePnt,5},{VType::cutEdge,4},{VType::cutEdge,8}},0};
        corrConcaveCellData redMarkSide =   {{{0,3,4,5,6},{1,2,3,7,8,9}},
                                             {mc33TriE513,mc33TriE538,mc33TriE845,addedRedMark,
                                              old1,old4,split0Right,old5,old3,split0Left}
                                            };        
    // RedMarkNotPlus        
        corrFace addedNonRedMark = {FType::added,{{VType::cutEdge,5},{VType::cubePnt,6},{VType::cubePnt,3},{VType::cutEdge,3}},-1};
        //Right Cell {4,5,6}
        old1 =      {FType::old,{},1};
        old4 =      {FType::old,{},4};
        corrFace split2Right={FType::splitOld,{{VType::cubePnt,6},{VType::cubePnt,2},{VType::cubePnt,3}},2};
        //Left Cell {7,8,9,10}
        old3 =       {FType::old,{},3};
        old5 =       {FType::old,{},5};
        corrFace old0 =       {FType::old,{},0};
        corrFace split2Left={FType::splitOld,{{VType::cubePnt,6},{VType::cubePnt,3},{VType::cubePnt,7}},2};
        corrConcaveCellData nonRedMarkSide = {{{0,3,4,5,6},{1,2,3,7,8,9,10}},
                                              {mc33TriE513,mc33TriE538,mc33TriE845,addedRedMark,
                                               old1,old4,split2Right,old3,old5,old0,split2Left}
                                             };
        
    convexCorrectionData[MC33::Case::c5] = {redMarkSide,nonRedMarkSide};
}
// Case c611
{
    corrFace mc33TriUp =       {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,5},{VType::cutEdge,9}},-2};
    corrFace mc33TriLow =      {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,8},{VType::cutEdge,9}},-2};
    corrFace mc33TriOther=     {FType::mc33Triangle,{{VType::cutEdge,1},{VType::cutEdge,2},{VType::cutEdge,10}},-2};
    // RedMarkPlus
        corrFace addedRedMark = {FType::added,{{VType::cutEdge,9},{VType::cubePnt,5},{VType::cutEdge,7}},-1};
        //Upper Cell {3,4}
        corrFace split5Up =   {FType::splitOld,{{VType::cubePnt,5},{VType::cutEdge,7},{VType::cutEdge,5}},5};
        corrFace old1 =                {FType::old,{},1};
        //Lower Cell {5,6,7}
        corrFace split5Low  = {FType::splitOld,{{VType::cubePnt,5},{VType::cubePnt,4},{VType::cutEdge,7}},5};
        corrFace old0 =                {FType::old,{},0};
        corrFace old3 =       {FType::old,{},3};
        corrConcaveCellData redMarkSide = {{{0,2,3,4},{1,2,5,6,7}},{mc33TriUp,mc33TriLow,addedRedMark,split5Up,old1,split5Low,old0,old3}};
        
    // RedMarkNotPlus
        corrFace addedNonRedMark={FType::added,{{VType::cutEdge,7},{VType::cutEdge,9},{VType::cubePnt,1},{VType::edgePntAdded,3}},-1};
        //Upper Cell {3,4,5,6,7}
        corrFace split3Up =   {FType::splitOld,{{VType::cutEdge,7},{VType::cubePnt,7},{VType::cubePnt,3},{VType::edgePntAdded,3}},3};
        corrFace split4Up  =  {FType::splitOld,{{VType::cubePnt,1},{VType::edgePntAdded,3},{VType::cubePnt,3},{VType::cutEdge,2},{VType::cutEdge,1}},4};
        corrFace old2 =       {FType::old,{},2};
        corrFace old5 =       {FType::old,{},5};
        old1 =       {FType::old,{},1};
        //Lower Cell {8,9,10}
        corrFace split3Low  = {FType::splitOld,{{VType::cutEdge,7},{VType::cubePnt,7},{VType::cubePnt,3},{VType::edgePntAdded,3}},3};
        corrFace split4Low  = {FType::splitOld,{{VType::cutEdge,7},{VType::cubePnt,7},{VType::cubePnt,3},{VType::edgePntAdded,3}},4};
        old0 =       {FType::old,{},0};
        corrConcaveCellData nonRedMarkSide = {{{0,2,3,4,5,6,7,11},{1,2,8,9,10}},{mc33TriUp,mc33TriLow,addedNonRedMark,split3Up,split4Up,old2,old5,old1,split3Low,split4Low,old0,mc33TriOther}};
        
    convexCorrectionData[MC33::Case::c611] = {redMarkSide,nonRedMarkSide};
}
// Case c612
{
    corrFace mc33Tri57A =   {FType::mc33Triangle,{{VType::cutEdge,5},{VType::cutEdge,7},{VType::cutEdge,10}},-2};
    corrFace mc33Tri7A2 =   {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,10},{VType::cutEdge,2}},-2};
    corrFace mc33Tri728 =   {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,2},{VType::cutEdge,8}},-2};
    corrFace mc33Tri281 =   {FType::mc33Triangle,{{VType::cutEdge,2},{VType::cutEdge,8},{VType::cutEdge,1}},-2};
    corrFace mc33Tri819 =   {FType::mc33Triangle,{{VType::cutEdge,8},{VType::cutEdge,1},{VType::cutEdge,9}},-2}; //4
    
    // RedMarkPlus
        // --         
    // RedMarkNotPlus
        corrFace addedNonRedMarkA77 ={FType::added,{{VType::cutEdge,10},{VType::cubePnt,7},{VType::cutEdge,7}},-1};
        corrFace addedNonRedMark237 ={FType::added,{{VType::cutEdge,2},{VType::cubePnt,3},{VType::cutEdge,7}},-1};
        corrFace addedNonRedMark238 ={FType::added,{{VType::cutEdge,2},{VType::cubePnt,3},{VType::cutEdge,8}},-1};
        corrFace addedNonRedMark108 ={FType::added,{{VType::cutEdge,1},{VType::cubePnt,0},{VType::cutEdge,8}},-1}; //8

        //Split Cell 1 {9,10,11}
        corrFace old5 =         {FType::old,{},5};
        corrFace old1Front =    {FType::old,{{VType::cubePnt,6}},1};
        corrFace split2Front=   {FType::splitOld,{{VType::cubePnt,6},{VType::cutEdge,10},{VType::cubePnt,7}},2};
        //Split Cell 2 {12,13}
        corrFace split2Back=    {FType::splitOld,{{VType::cubePnt,7},{VType::cutEdge,10},{VType::cutEdge,2},{VType::cubePnt,3}},2};
        corrFace split3High=    {FType::splitOld,{{VType::cubePnt,7},{VType::cubePnt,3},{VType::cutEdge,7}},3};
        //Split Cell 3 {14}
        corrFace split3Middle = {FType::splitOld,{{VType::cutEdge,7},{VType::cubePnt,3},{VType::cutEdge,8}},3};
        //Split Cell 4 {15,16}
        corrFace split4High =   {FType::splitOld,{{VType::cutEdge,2},{VType::cutEdge,1},{VType::cubePnt,0},{VType::cubePnt,3}},4};
        corrFace split3Low =    {FType::splitOld,{{VType::cutEdge,8},{VType::cubePnt,3},{VType::cubePnt,0}},3};
        //Split Cell 5 {17,18,19}
        corrFace old0 =         {FType::old,{},0};
        corrFace old1Back =     {FType::old,{{VType::cubePnt,1}},1};
        corrFace split4Low =    {FType::splitOld,{{VType::cutEdge,1},{VType::cubePnt,1},{VType::cubePnt,0}},4};
        
        corrConcaveCellData nonRedMarkSide = {{{0,5,9,10,11},{1,5,6,12,13},{2,6,7,14},{3,7,8,15,16},{4,8,17,18,19}},
                            {mc33Tri57A,mc33Tri7A2,mc33Tri728,mc33Tri281,mc33Tri819,
                             addedNonRedMarkA77,addedNonRedMark237,addedNonRedMark238,addedNonRedMark108,
                             old5,old1Front,split2Front,split2Back,split3High,split3Middle,split4High,split3Low,old0,old1Back,split4Low
                            }};
        
    convexCorrectionData[MC33::Case::c612] = {{{{}},{}},nonRedMarkSide};
}
// Case c62
{
    corrFace mc33Tri57A =   {FType::mc33Triangle,{{VType::cutEdge,5},{VType::cutEdge,7},{VType::cutEdge,10}},-2};
    corrFace mc33Tri7A2 =   {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,10},{VType::cutEdge,2}},-2};
    corrFace mc33Tri728 =   {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,2},{VType::cutEdge,8}},-2};
    corrFace mc33Tri281 =   {FType::mc33Triangle,{{VType::cutEdge,2},{VType::cutEdge,8},{VType::cutEdge,1}},-2};
    corrFace mc33Tri819 =   {FType::mc33Triangle,{{VType::cutEdge,8},{VType::cutEdge,1},{VType::cutEdge,9}},-2}; //4
    
    // RedMarkPlus
        // --         
    // RedMarkNotPlus
        corrFace addedNonRedMarkA77 ={FType::added,{{VType::cutEdge,10},{VType::cubePnt,7},{VType::cutEdge,7}},-1};
        corrFace addedNonRedMark237 ={FType::added,{{VType::cutEdge,2},{VType::cubePnt,3},{VType::cutEdge,7}},-1};
        corrFace addedNonRedMark238 ={FType::added,{{VType::cutEdge,2},{VType::cubePnt,3},{VType::cutEdge,8}},-1};
        corrFace addedNonRedMark108 ={FType::added,{{VType::cutEdge,1},{VType::cubePnt,0},{VType::cutEdge,8}},-1}; //8

        //Split Cell 1 {9,10,11}
        corrFace old5 =         {FType::old,{},5};
        corrFace old1Front =    {FType::old,{{VType::cubePnt,6}},1};
        corrFace split2Front=   {FType::splitOld,{{VType::cubePnt,6},{VType::cutEdge,10},{VType::cubePnt,7}},2};
        //Split Cell 2 {12,13}
        corrFace split2Back=    {FType::splitOld,{{VType::cubePnt,7},{VType::cutEdge,10},{VType::cutEdge,2},{VType::cubePnt,3}},2};
        corrFace split3High=    {FType::splitOld,{{VType::cubePnt,7},{VType::cubePnt,3},{VType::cutEdge,7}},3};
        //Split Cell 3 {14}
        corrFace split3Middle = {FType::splitOld,{{VType::cutEdge,7},{VType::cubePnt,3},{VType::cutEdge,8}},3};
        //Split Cell 4 {15,16}
        corrFace split4High =   {FType::splitOld,{{VType::cutEdge,2},{VType::cutEdge,1},{VType::cubePnt,0},{VType::cubePnt,3}},4};
        corrFace split3Low =    {FType::splitOld,{{VType::cutEdge,8},{VType::cubePnt,3},{VType::cubePnt,0}},3};
        //Split Cell 5 {17,18,19}
        corrFace old0 =         {FType::old,{},0};
        corrFace old1Back =     {FType::old,{{VType::cubePnt,1}},1};
        corrFace split4Low =    {FType::splitOld,{{VType::cutEdge,1},{VType::cubePnt,1},{VType::cubePnt,0}},4};
        
        corrConcaveCellData nonRedMarkSide = {{{0,5,9,10,11},{1,5,6,12,13},{2,6,7,14},{3,7,8,15,16},{4,8,17,18,19}},
                            {mc33Tri57A,mc33Tri7A2,mc33Tri728,mc33Tri281,mc33Tri819,
                             addedNonRedMarkA77,addedNonRedMark237,addedNonRedMark238,addedNonRedMark108,
                             old5,old1Front,split2Front,split2Back,split3High,split3Middle,split4High,split3Low,old0,old1Back,split4Low
                            }};
        
    convexCorrectionData[MC33::Case::c62] = {{{{}},{}},nonRedMarkSide};
}
// Case c71
//---

// Case c72_720
{
    corrFace mc33TriLeftLeft = {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,4},{VType::cutEdge,9}},-2};
    corrFace mc33TriLeft =     {FType::mc33Triangle,{{VType::cutEdge,9},{VType::cutEdge,7},{VType::cutEdge,11}},-2};
    corrFace mc33TriRight =    {FType::mc33Triangle,{{VType::cutEdge,9},{VType::cutEdge,11},{VType::cutEdge,5}},-2};
    corrFace mc33TriRightRight={FType::mc33Triangle,{{VType::cutEdge,11},{VType::cutEdge,5},{VType::cutEdge,6}},-2};    
    corrFace mc33TriOther=     {FType::mc33Triangle,{{VType::cutEdge,1},{VType::cutEdge,2},{VType::cutEdge,10}},-2}; //4
    // RedMarkPlus
        //---
    // RedMarkNotPlus
        corrFace addedNonRedMark={FType::added,{{VType::cutEdge,11},{VType::cubePnt,3},{VType::cubePnt,1},{VType::cutEdge,9}},-1};
        //Top Cell {6,7,8}
        corrFace old1 =       {FType::old,{},1};
        corrFace old2 =       {FType::old,{},2};
        corrFace split4Top= {FType::splitOld,{{VType::cutEdge,2},{VType::cubePnt,3},{VType::cubePnt,1},{VType::cutEdge,1}},4};
        //Bottom Cell {9,10,11}
        corrFace old3 =       {FType::old,{},3};
        corrFace old0 =       {FType::old,{},0};
        corrFace split4Down = {FType::splitOld,{{VType::cubePnt,3},{VType::cubePnt,0},{VType::cubePnt,1}},4};
        
        corrConcaveCellData nonRedMarkSide = {{{2,3,4,5,6,7,8},
                                               {0,1,5,9,10,11}},{mc33TriLeftLeft,mc33TriLeft,mc33TriRight,mc33TriRightRight,mc33TriOther,
                                                addedNonRedMark,
                                                old1,old2,split4Top,old3,old0,split4Down}};
        
    convexCorrectionData[MC33::Case::c72_720] = {{{{}},{}},nonRedMarkSide};
}
// Case c72_760
{
    corrFace mc33TriBackBack = {FType::mc33Triangle,{{VType::cutEdge,9},{VType::cutEdge,1},{VType::cutEdge,2}},-2};
    corrFace mc33TriBack =     {FType::mc33Triangle,{{VType::cutEdge,4},{VType::cutEdge,9},{VType::cutEdge,2}},-2};
    corrFace mc33TriFront =    {FType::mc33Triangle,{{VType::cutEdge,2},{VType::cutEdge,4},{VType::cutEdge,10}},-2};
    corrFace mc33TriFrontFront={FType::mc33Triangle,{{VType::cutEdge,4},{VType::cutEdge,5},{VType::cutEdge,10}},-2};    
    corrFace mc33TriOther=     {FType::mc33Triangle,{{VType::cutEdge,6},{VType::cutEdge,7},{VType::cutEdge,11}},-2}; //4
    // RedMarkPlus
        //---
    // RedMarkNotPlus
        corrFace addedNonRedMark={FType::added,{{VType::cutEdge,2},{VType::cubePnt,3},{VType::cubePnt,4},{VType::cutEdge,4}},-1};
        //Front Cell {6,7,8}
        corrFace old2 =       {FType::old,{},2};
        corrFace old5 =       {FType::old,{},5};
        corrFace split3Front= {FType::splitOld,{{VType::cutEdge,11},{VType::cubePnt,3},{VType::cubePnt,4},{VType::cutEdge,7}},3};
        //Back Cell {9,10,11}
        corrFace old0 =       {FType::old,{},0};
        corrFace old4 =       {FType::old,{},4};
        corrFace split3Back = {FType::splitOld,{{VType::cubePnt,0},{VType::cubePnt,4},{VType::cubePnt,3}},3};
        
        corrConcaveCellData nonRedMarkSide = {{{2,3,4,5,6,7,8},
                                               {0,1,5,9,10,11}},
                                               {mc33TriBackBack,mc33TriBack,mc33TriFront,mc33TriFrontFront,mc33TriOther,
                                                addedNonRedMark,
                                                old2,old5,split3Front,old0,old4,split3Back}};
        
    convexCorrectionData[MC33::Case::c72_760] = {{{{}},{}},nonRedMarkSide};
}
// Case c72_800
{
    corrFace mc33TriBackBack = {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,11},{VType::cutEdge,2}},-2};
    corrFace mc33TriBack =     {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,1},{VType::cutEdge,2}},-2};
    corrFace mc33TriFront =    {FType::mc33Triangle,{{VType::cutEdge,7},{VType::cutEdge,6},{VType::cutEdge,1}},-2};
    corrFace mc33TriFrontFront={FType::mc33Triangle,{{VType::cutEdge,1},{VType::cutEdge,6},{VType::cutEdge,10}},-2};    
    corrFace mc33TriOther=     {FType::mc33Triangle,{{VType::cutEdge,4},{VType::cutEdge,5},{VType::cutEdge,9}},-2}; //4
    // RedMarkPlus
        //---
    // RedMarkNotPlus
        corrFace addedNonRedMark={FType::added,{{VType::cutEdge,7},{VType::cubePnt,4},{VType::cubePnt,1},{VType::cutEdge,1}},-1};
        //Front Cell {6,7,8}
        corrFace old1 =       {FType::old,{},1};
        corrFace old5 =       {FType::old,{},5};
        corrFace split0Front= {FType::splitOld,{{VType::cutEdge,9},{VType::cubePnt,1},{VType::cubePnt,4},{VType::cutEdge,4}},0};
        //Back Cell {9,10,11}
        corrFace old3 =       {FType::old,{},3};
        corrFace old4 =       {FType::old,{},4};
        corrFace split0Back = {FType::splitOld,{{VType::cubePnt,4},{VType::cubePnt,0},{VType::cubePnt,1}},0};
        
        corrConcaveCellData nonRedMarkSide = {{{2,3,4,5,6,7,8},
                                               {0,1,5,9,10,11}},{mc33TriBackBack,mc33TriBack,mc33TriFront,mc33TriFrontFront,mc33TriOther,
                                                addedNonRedMark,
                                                old1,old5,split0Front,old3,old4,split0Back}};
        
    convexCorrectionData[MC33::Case::c72_800] = {{{{}},{}},nonRedMarkSide};
}
// Case c73_840
{
    corrFace mc33TriB7 =   {FType::mc33Triangle,{{VType::cutEdge,11},{VType::cutEdge,7}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri74 =   {FType::mc33Triangle,{{VType::cutEdge,7}, {VType::cutEdge,4}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri49 =   {FType::mc33Triangle,{{VType::cutEdge,4}, {VType::cutEdge,9}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri91 =   {FType::mc33Triangle,{{VType::cutEdge,9}, {VType::cutEdge,1}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri12 =   {FType::mc33Triangle,{{VType::cutEdge,1}, {VType::cutEdge,2}, {VType::centerPnt,-1}},-2}; //4
    corrFace mc33Tri2A =   {FType::mc33Triangle,{{VType::cutEdge,2}, {VType::cutEdge,10},{VType::centerPnt,-1}},-2}; 
    corrFace mc33TriA5 =   {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,5}, {VType::centerPnt,-1}},-2}; 
    corrFace mc33Tri56 =   {FType::mc33Triangle,{{VType::cutEdge,5}, {VType::cutEdge,6}, {VType::centerPnt,-1}},-2}; 
    corrFace mc33Tri6B =   {FType::mc33Triangle,{{VType::cutEdge,6}, {VType::cutEdge,11},{VType::centerPnt,-1}},-2}; //8
    
    // RedMarkPlus
        corrFace addedRedMark54 ={FType::added,{{VType::cutEdge,5},{VType::cutEdge,4},{VType::centerPnt,-1}},-1};
        corrFace addedRedMark59 ={FType::added,{{VType::cutEdge,5},{VType::cutEdge,9},{VType::centerPnt,-1}},-1};
        //Cell Left {11,12,13}
        corrFace old3 =        {FType::old,{},3};
        corrFace old2Left =    {FType::old,{{VType::cubePnt,7}},2};
        corrFace split5Left=   {FType::splitOld,{{VType::cutEdge,5},{VType::cutEdge,6},{VType::cubePnt,7},{VType::cutEdge,7},{VType::cutEdge,4}},5};
         //Cell Right {14,15,16}
        corrFace old0 =        {FType::old,{},0};
        corrFace split5Right=  {FType::splitOld,{{VType::cutEdge,5},{VType::cubePnt,5},{VType::cutEdge,4}},5};
        corrFace split1Front=  {FType::splitOld,{{VType::cutEdge,5},{VType::cubePnt,5},{VType::cutEdge,9}},1};
        //Cell Right {17,18,19}
        corrFace old4 =        {FType::old,{},4};
        corrFace old2Right=    {FType::old,{{VType::cubePnt,2}},2};
        corrFace split1Back=   {FType::splitOld,{{VType::cutEdge,5},{VType::cutEdge,10},{VType::cubePnt,2},{VType::cutEdge,1},{VType::cutEdge,9}},1};
        
        corrConcaveCellData redMarkSide = {{{7,8,0,1,9,11,12,13},
                                            {2,9,10,14,15,16},
                                            {3,4,5,6,10,17,18,19}},
                            {mc33TriB7,mc33Tri74,mc33Tri49,mc33Tri91,mc33Tri12,mc33Tri2A,mc33TriA5,mc33Tri56,mc33Tri6B
                             addedRedMark54,addedRedMark59,
                             old3,old2Left,split5Left,old0,split5Right,split1Front,old4,old2Right,split1Back}
                            };
                
    // RedMarkNotPlus
        corrFace addedNonRedMarkB2 ={FType::added,{{VType::cutEdge,11},{VType::cutEdge,2},{VType::centerPnt,-1}},-1};
        //Cell Top {10,11,12}
        corrFace old5Top =     {FType::old,{{VType::cubePnt,6}},5};
        corrFace old1Top =     {FType::old,{{VType::cubePnt,6}},1};
        corrFace split2Front=  {FType::splitOld,{{VType::cutEdge,2},{VType::cutEdge,10},{VType::cubePnt,6},{VType::cutEdge,6},{VType::cutEdge,11}},2};
        //Cell Back {13,14,15,16,17,18}
        old3 =                 {FType::old,{},3};
        old4 =                 {FType::old,{},4};
        old0 =                 {FType::old,{},0};        
        corrFace old5Left =    {FType::old,{{VType::cubePnt,4}},5};
        corrFace old1Back =    {FType::old,{{VType::cubePnt,3}},1};
        corrFace split2Back =  {FType::splitOld,{{VType::cutEdge,2},{VType::cubePnt,3},{VType::cutEdge,11}},2};

        
        corrConcaveCellData nonRedMarkSide = {{{5,6,7,8,9,10,11,12},
                                               {0,1,2,3,4,9,13,14,15,16,17,18}},
                            {mc33TriB7,mc33Tri74,mc33Tri49,mc33Tri91,mc33Tri12,mc33Tri2A,mc33TriA5,mc33Tri56,mc33Tri6B
                             addedNonRedMarkB2,
                             old5Top,old1Top,split2Front,old3,old4,old0,old5Left,old1Back,split2Back
                            }};
        
    convexCorrectionData[MC33::Case::c73_840] = {redMarkSide,nonRedMarkSide};
}
//Case c73_912
{
    corrFace mc33Tri65 =   {FType::mc33Triangle,{{VType::cutEdge,6}, {VType::cutEdge,5}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri59 =   {FType::mc33Triangle,{{VType::cutEdge,5}, {VType::cutEdge,9}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri94 =   {FType::mc33Triangle,{{VType::cutEdge,9}, {VType::cutEdge,4}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri47 =   {FType::mc33Triangle,{{VType::cutEdge,4}, {VType::cutEdge,7}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri7B =   {FType::mc33Triangle,{{VType::cutEdge,7}, {VType::cutEdge,11},{VType::centerPnt,-1}},-2}; //4
    corrFace mc33TriB2 =   {FType::mc33Triangle,{{VType::cutEdge,11},{VType::cutEdge,2}, {VType::centerPnt,-1}},-2}; 
    corrFace mc33Tri21 =   {FType::mc33Triangle,{{VType::cutEdge,2}, {VType::cutEdge,1}, {VType::centerPnt,-1}},-2}; 
    corrFace mc33Tri1A =   {FType::mc33Triangle,{{VType::cutEdge,1}, {VType::cutEdge,10},{VType::centerPnt,-1}},-2}; 
    corrFace mc33TriA6 =   {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,6}, {VType::centerPnt,-1}},-2}; //8
    
    // RedMarkPlus
        corrFace addedRedMark67 ={FType::added,{{VType::cutEdge,6},{VType::cutEdge,7},{VType::centerPnt,-1}},-1};
        corrFace addedRedMark6B ={FType::added,{{VType::cutEdge,6},{VType::cutEdge,11},{VType::centerPnt,-1}},-1};
        //Cell Bottom {11,12,13}
        corrFace old0 =        {FType::old,{},0};
        corrFace old1Front=    {FType::old,{{VType::cubePnt,5}},1};
        corrFace split5Right=  {FType::splitOld,{{VType::cutEdge,6},{VType::cutEdge,5},{VType::cubePnt,5},{VType::cutEdge,4},{VType::cutEdge,7}},5};
         //Cell Middle {14,15,16}
        corrFace old3 =        {FType::old,{},3};
        corrFace split5Left =  {FType::splitOld,{{VType::cutEdge,6},{VType::cubePnt,7},{VType::cutEdge,7}},5};
        corrFace split2Left =  {FType::splitOld,{{VType::cutEdge,6},{VType::cubePnt,7},{VType::cutEdge,11}},2};
        //Cell Right {17,18,19}
        corrFace old4 =        {FType::old,{},4};
        corrFace old1Back =    {FType::old,{{VType::cubePnt,2}},1};
        corrFace split2Right=  {FType::splitOld,{{VType::cutEdge,6},{VType::cutEdge,10},{VType::cubePnt,2},{VType::cutEdge,2},{VType::cutEdge,11}},2};
        
        corrConcaveCellData redMarkSide = {{{0,1,2,3,9,11,12,13},
                                            {4,9,10,14,15,16},
                                            {5,6,7,8,10,17,18,19}},
                            {mc33Tri65,mc33Tri59,mc33Tri94,mc33Tri47,mc33Tri7B,mc33TriB2,mc33Tri21,mc33Tri1A,mc33TriA6
                             addedRedMark67,addedRedMark6B,
                             old0,old1Front,split5Right,old3,split5Left,split2Left,old4,old1Back,split2Right}
                            };
                
    // RedMarkNotPlus
        corrFace addedNonRedMark19 ={FType::added,{{VType::cutEdge,1},{VType::cutEdge,9},{VType::centerPnt,-1}},-1};
        //Cell Right {10,11,12}
        corrFace old2Front =   {FType::old,{{VType::cubePnt,6}},2};
        corrFace old5Right =   {FType::old,{{VType::cubePnt,6}},5};
        corrFace split1Front=  {FType::splitOld,{{VType::cutEdge,1},{VType::cutEdge,10},{VType::cubePnt,6},{VType::cutEdge,5},{VType::cutEdge,9}},1};
        //Cell Back {13,14,15,16,17,18}
        old3 =                 {FType::old,{},3};
        old4 =                 {FType::old,{},4};
        old0 =                 {FType::old,{},0};        
        corrFace old5Left =    {FType::old,{{VType::cubePnt,4}},5};
        corrFace old2Back =    {FType::old,{{VType::cubePnt,3}},2};
        corrFace split1Back =  {FType::splitOld,{{VType::cutEdge,9},{VType::cubePnt,1},{VType::cutEdge,1}},1};

        
        corrConcaveCellData nonRedMarkSide = {{{0,1,7,8,9,10,11,12},
                                               {2,3,4,5,6,9,13,14,15,16,17,18}},
                            {mc33Tri65,mc33Tri59,mc33Tri94,mc33Tri47,mc33Tri7B,mc33TriB2,mc33Tri21,mc33Tri1A,mc33TriA6
                             addedNonRedMark19,
                             old2Front,old5Right,split1Front,old3,old4,old0,old5Left,old2Back,split1Back
                            }};
        
    convexCorrectionData[MC33::Case::c73_840] = {redMarkSide,nonRedMarkSide};
}
// Case c73_984
{
    corrFace mc33TriA6 =   {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,6}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri67 =   {FType::mc33Triangle,{{VType::cutEdge,6}, {VType::cutEdge,7}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri7B =   {FType::mc33Triangle,{{VType::cutEdge,7}, {VType::cutEdge,11},{VType::centerPnt,-1}},-2};
    corrFace mc33TriB2 =   {FType::mc33Triangle,{{VType::cutEdge,11},{VType::cutEdge,2}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri21 =   {FType::mc33Triangle,{{VType::cutEdge,2}, {VType::cutEdge,1}, {VType::centerPnt,-1}},-2}; //4
    corrFace mc33Tri19 =   {FType::mc33Triangle,{{VType::cutEdge,1}, {VType::cutEdge,9}, {VType::centerPnt,-1}},-2}; 
    corrFace mc33Tri94 =   {FType::mc33Triangle,{{VType::cutEdge,9}, {VType::cutEdge,4}, {VType::centerPnt,-1}},-2}; 
    corrFace mc33Tri45 =   {FType::mc33Triangle,{{VType::cutEdge,4}, {VType::cutEdge,5}, {VType::centerPnt,-1}},-2}; 
    corrFace mc33Tri5A =   {FType::mc33Triangle,{{VType::cutEdge,5}, {VType::cutEdge,10},{VType::centerPnt,-1}},-2}; //8
    
    // RedMarkPlus
        corrFace addedRedMarkA2 ={FType::added,{{VType::cutEdge,10},{VType::cutEdge,2},{VType::centerPnt,-1}},-1};
        corrFace addedRedMarkA1 ={FType::added,{{VType::cutEdge,10},{VType::cutEdge,1},{VType::centerPnt,-1}},-1};
        //Cell Top {11,12,13}
        corrFace old3 =        {FType::old,{},3};
        corrFace old5Top =     {FType::old,{{VType::cubePnt,7}},5};
        corrFace split2Left=   {FType::splitOld,{{VType::cutEdge,2},{VType::cutEdge,11},{VType::cubePnt,7},{VType::cutEdge,6},{VType::cutEdge,10}},2};
         //Cell Middle {14,15,16}
        corrFace old4 =        {FType::old,{},4};
        corrFace split1Back =  {FType::splitOld,{{VType::cutEdge,10},{VType::cubePnt,2},{VType::cutEdge,1}},1};
        corrFace split2Right = {FType::splitOld,{{VType::cutEdge,10},{VType::cubePnt,2},{VType::cutEdge,2}},2};
        //Cell Right {17,18,19}
        corrFace old0       =  {FType::old,{},0};
        corrFace split1Front=  {FType::old,{{VType::cutEdge,10},{VType::cutEdge,5},{VType::cubePnt,5},{VType::cutEdge,9},{VType::cutEdge,1}},1};
        corrFace old5Bottom =  {FType::old,{{VType::cubePnt,5}},5};
        
        corrConcaveCellData redMarkSide = {{{0,1,2,3,9,11,12,13},
                                            {4,9,10,14,15,16},
                                            {5,6,7,8,10,17,18,19}},
                            {mc33TriA6,mc33Tri67,mc33Tri7B,mc33TriB2,mc33Tri21,mc33Tri19,mc33Tri94,mc33Tri45,mc33Tri5A,
                             addedRedMarkA2,addedRedMarkA1,
                             old3,old5Top,split2Left,old4,split1Back,split2Right,old0,split1Front,old5Bottom
                            }};
                
    // RedMarkNotPlus
        corrFace addedNonRedMark74 ={FType::added,{{VType::cutEdge,7},{VType::cutEdge,4},{VType::centerPnt,-1}},-1};
        //Cell Front {10,11,12}
        corrFace old2Front =   {FType::old,{{VType::cubePnt,6}},2};
        corrFace old1Top =     {FType::old,{{VType::cubePnt,6}},1};
        corrFace split5Right=  {FType::splitOld,{{VType::cutEdge,4},{VType::cutEdge,5},{VType::cubePnt,6},{VType::cutEdge,6},{VType::cutEdge,7}},5};
        //Cell Back {13,14,15,16,17,18}
        old3       =           {FType::old,{},3};
        old4 =                 {FType::old,{},4};
        old0 =                 {FType::old,{},0};        
        corrFace old2Back =    {FType::old,{{VType::cubePnt,3}},2};
        corrFace old1Bottom=   {FType::old,{{VType::cubePnt,1}},1};
        corrFace split5Left =  {FType::splitOld,{{VType::cutEdge,4},{VType::cubePnt,4},{VType::cutEdge,7}},5};
        
        corrConcaveCellData nonRedMarkSide = {{{0,1,7,8,9,10,11,12},
                                               {2,3,4,5,6,9,13,14,15,16,17,18}},
                            {mc33TriA6,mc33Tri57,mc33Tri7B,mc33TriB2,mc33Tri21,mc33Tri19,mc33Tri94,mc33Tri45,mc33Tri5A,
                             addedNonRedMark74,
                             old2Front,old1Top,old5,old3,old4,old0,old2Back,old1Bottom
                            }};
        
    convexCorrectionData[MC33::Case::c73_984] = {redMarkSide,nonRedMarkSide};
}
// Case c741
//---

// Case c742
//---

// Case c8
{
    corrFace mc33TriFront =   {FType::mc33Triangle,{{VType::cutEdge,5},{VType::cutEdge,7},{VType::cutEdge,3}},-2};
    corrFace mc33TriBack =    {FType::mc33Triangle,{{VType::cutEdge,5},{VType::cutEdge,1},{VType::cutEdge,5}},-2};
    // RedMarkPlus
        corrFace addedRedMark={FType::added,{{VType::cutEdge,5},{VType::cubePnt,5},{VType::cubePnt,0},{VType::cutEdge,3}},-1};
        //Front Cell {3,4,5}
        corrFace old3 =       {FType::old,{},3};
        corrFace old5 =       {FType::old,{},5};
        corrFace split0Front= {FType::splitOld,{{VType::cubePnt,5},{VType::cubePnt,4},{VType::cubePnt,0}},0};
        //Back Cell {6,7,8}
        corrFace old1 =       {FType::old,{},1};
        corrFace old4 =       {FType::old,{},4};
        corrFace split0Back = {FType::splitOld,{{VType::cubePnt,0},{VType::cubePnt,1},{VType::cubePnt,5}},0};
        
        corrConcaveCellData redMarkSide = {{{0,2,3,4,5},
                                            {1,2,6,7,8}},
                                    {mc33TriFront,mc33TriBack,
                                     addedRedMark,
                                     old3,old5,split0Front,old1,old4,split0Back}};

    // RedMarkNotPlus
        corrFace addedNonRedMark={FType::added,{{VType::cutEdge,5},{VType::cubePnt,5},{VType::cubePnt,3},{VType::cutEdge,3}},-1};
        //Front Cell {3,4,5}
        old3 =       {FType::old,{},3};
        old5 =       {FType::old,{},5};
        corrFace split2Front= {FType::splitOld,{{VType::cubePnt,6},{VType::cubePnt,7},{VType::cubePnt,3}},2};
        //Back Cell {6,7,8}
        old1 =       {FType::old,{},1};
        old4 =       {FType::old,{},4};
        corrFace split2Back = {FType::splitOld,{{VType::cubePnt,3},{VType::cubePnt,2},{VType::cubePnt,6}},2};
        
        corrConcaveCellData nonRedMarkSide = {{{0,2,3,4,5},
                                               {1,2,6,7,8}},
                                    {mc33TriFront,mc33TriBack,
                                     addedRedMark,
                                     old3,old5,split2Front,old1,old4,split2Back}};

    convexCorrectionData[MC33::Case::c8] = {redMarkSide,nonRedMarkSide};
}
// Case c9
//---

// Case c1011
{
    /*
    corrFace mc33TriA20  =   {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,2},{VType::cutEdge,0}},-2};
    corrFace mc33TriA09  =   {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,0},{VType::cutEdge,9}},-2};
    corrFace mc33TriB68  =   {FType::mc33Triangle,{{VType::cutEdge,11},{VType::cutEdge,6},{VType::cutEdge,8}},-2};
    corrFace mc33Tri648  =   {FType::mc33Triangle,{{VType::cutEdge,6},{VType::cutEdge,4},{VType::cutEdge,8}},-2};
    // RedMarkPlus
        corrFace addedRedMark={FType::added,{{VType::cutEdge,7},{VType::cubePnt,4},{VType::cubePnt,1},{VType::cutEdge,1}},-1};
        //Front Cell {3,4,5}
        corrFace old1 =       {FType::old,{},1};
        corrFace old5 =       {FType::old,{},5};
        corrFace split0Front= {FType::splitOld,{{VType::cubePnt,4},{VType::cubePnt,5},{VType::cubePnt,1}},0};
        //Back Cell {6,7,8}
        corrFace old3 =       {FType::old,{},3};
        corrFace old4 =       {FType::old,{},4};
        corrFace split0Back = {FType::splitOld,{{VType::cubePnt,1},{VType::cubePnt,0},{VType::cubePnt,4}},0};
        
        corrConcaveCellData redMarkSide = {{{0,2,3,4,5},{1,2,6,7,8}},{mc33TriFront,mc33TriBack,addedRedMark,old1,old5,split0Front,old3,old4,split0Back}};
         
    // RedMarkNotPlus
        corrFace addedNonRedMark={FType::added,{{VType::cutEdge,7},{VType::cubePnt,7},{VType::cubePnt,2},{VType::cutEdge,1}},-1};
        //Front Cell {3,4,5}
        old1 =       {FType::old,{},5};
        old5 =       {FType::old,{},1};
        corrFace split2Front= {FType::splitOld,{{VType::cubePnt,2},{VType::cubePnt,6},{VType::cubePnt,7}},2};
        //Back Cell {6,7,8}
        old3 =       {FType::old,{},3};
        old4 =       {FType::old,{},4};
        corrFace split2Back = {FType::splitOld,{{VType::cubePnt,2},{VType::cubePnt,3},{VType::cubePnt,7}},2};
        
        corrConcaveCellData nonRedMarkSide = {{{0,2,3,4,5},{1,2,6,7,8}},{mc33TriFront,mc33TriBack,addedNonRedMark,old1,old5,split2Front,old3,old4,split2Back}};
        
    convexCorrectionData[MC33::Case::c8] = {redMarkSide,nonRedMarkSide};
    */
}
// Case c1012
//---
// Case c102_1262
{
    corrFace mc33Tri64  =   {FType::mc33Triangle,{{VType::cutEdge,6}, {VType::cutEdge,4}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri48  =   {FType::mc33Triangle,{{VType::cutEdge,4}, {VType::cutEdge,8}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri8B  =   {FType::mc33Triangle,{{VType::cutEdge,8}, {VType::cutEdge,11},{VType::centerPnt,-1}},-2};
    corrFace mc33TriB2  =   {FType::mc33Triangle,{{VType::cutEdge,11},{VType::cutEdge,2}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri20  =   {FType::mc33Triangle,{{VType::cutEdge,2}, {VType::cutEdge,0}, {VType::centerPnt,-1}},-2};//4
    corrFace mc33Tri09  =   {FType::mc33Triangle,{{VType::cutEdge,0}, {VType::cutEdge,9}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri9A  =   {FType::mc33Triangle,{{VType::cutEdge,9}, {VType::cutEdge,10},{VType::centerPnt,-1}},-2};
    corrFace mc33TriA6  =   {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,6}, {VType::centerPnt,-1}},-2};

    // RedMarkPlus
        corrFace addedRedMark49={FType::added,{{VType::cutEdge,4},{VType::cutEdge,9},{VType::centerPnt,-1}},-1};//8
        corrFace addedRedMark08={FType::added,{{VType::cutEdge,0},{VType::cutEdge,8},{VType::centerPnt,-1}},-1};
        //Front Cell {10,11,12,13}
        corrFace old2Front=   {FType::old,{{VType::cubePnt,6}},2};
        corrFace old5 =       {FType::old,{},5};
        corrFace old1 =       {FType::old,{},1};
        corrFace split0Front= {FType::splitOld,{{VType::cutEdge,4},{VType::cubePnt,5},{VType::cutEdge,9}},0};
        //Middle Cell {14}
        corrFace split0Middle={FType::splitOld,{{VType::cutEdge,4},{VType::cutEdge,9},{VType::cutEdge,0},{VType::cutEdge,8}},0};
        //Right Cell {15,16,17,18}
        corrFace old2Back  =  {FType::old,{{VType::cubePnt,3}},2};
        corrFace old3 =       {FType::old,{},3};
        corrFace old4 =       {FType::old,{},4};
        corrFace split0Back = {FType::splitOld,{{VType::cutEdge,8},{VType::cubePnt,0},{VType::cutEdge,0}},0};
        
        corrConcaveCellData redMarkSide = {{{0,6,7,8,10,11,12,13},
                                            {1,5,8,9,14},
                                            {2,3,4,9,15,16,17,18}}, 
                            {mc33Tri64,mc33Tri48,mc33Tri8B,mc33TriB2,mc33Tri20,mc33Tri09,mc33Tri9A,mc33TriA6,
                             addedRedMark49,addedRedMark08,
                             old2Front,old5,old1,split0Front,split0Middle,old2Back,old3,old4,split0Back
                    }};
         
    // RedMarkNotPlus
        corrFace addedRedMarkB6={FType::added,{{VType::cutEdge,11},{VType::cutEdge,6},{VType::centerPnt,-1}},-1};//8
        corrFace addedRedMark2A={FType::added,{{VType::cutEdge,2},{VType::cutEdge,10},{VType::centerPnt,-1}},-1};
        //Left Cell {10,11,12,13}
        corrFace old0Left =   {FType::old,{{VType::cubePnt,4}},0};
        old3 =                {FType::old,{},3};
        old5 =                {FType::old,{},5};
        corrFace split2Left = {FType::splitOld,{{VType::cubePnt,7},{VType::cutEdge,11},{VType::cutEdge,6}},2};
        //Middle Cell {14}
        corrFace split2Middle={FType::splitOld,{{VType::cutEdge,6},{VType::cutEdge,10},{VType::cutEdge,2},{VType::cutEdge,11}},2};
        //Back Cell {15,16,17,18}
        corrFace old0Right=   {FType::old,{{VType::cubePnt,1}},2};
        old1 =                {FType::old,{},1};
        old4 =                {FType::old,{},4};
        corrFace split2Right= {FType::splitOld,{{VType::cubePnt,2},{VType::cutEdge,10},{VType::cutEdge,2}},2};
        
        corrConcaveCellData nonRedMarkSide = {{{0,1,2,8,10,11,12,13},
                                               {3,7,8,9,14},
                                               {4,5,6,9,15,16,17,18}}, 
                            {mc33Tri64,mc33Tri48,mc33Tri8B,mc33TriB2,mc33Tri20,mc33Tri09,mc33Tri9A,mc33TriA6,
                            addedRedMarkB6,addedRedMark2A,
                            old0Left,old3,old5,split2Left,split2Middle,old0Right,old1,old4,split2Right
                            }};
        
    convexCorrectionData[MC33::Case::c102_1262] = {redMarkSide,nonRedMarkSide};
}
// Case c102_1286
{
    corrFace mc33Tri49  =   {FType::mc33Triangle,{{VType::cutEdge,4}, {VType::cutEdge,9}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri9A  =   {FType::mc33Triangle,{{VType::cutEdge,9}, {VType::cutEdge,10},{VType::centerPnt,-1}},-2};
    corrFace mc33TriA2  =   {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,2}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri20  =   {FType::mc33Triangle,{{VType::cutEdge,2}, {VType::cutEdge,0}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri08  =   {FType::mc33Triangle,{{VType::cutEdge,0}, {VType::cutEdge,8}, {VType::centerPnt,-1}},-2};//4
    corrFace mc33Tri8B  =   {FType::mc33Triangle,{{VType::cutEdge,8}, {VType::cutEdge,11},{VType::centerPnt,-1}},-2};
    corrFace mc33TriB6  =   {FType::mc33Triangle,{{VType::cutEdge,11},{VType::cutEdge,6}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri64  =   {FType::mc33Triangle,{{VType::cutEdge,6}, {VType::cutEdge,4}, {VType::centerPnt,-1}},-2};

    // RedMarkPlus
        corrFace addedRedMark6A={FType::added,{{VType::cutEdge,6},{VType::cutEdge,10},{VType::centerPnt,-1}},-1};//8
        corrFace addedRedMarkB2={FType::added,{{VType::cutEdge,11},{VType::cutEdge,2},{VType::centerPnt,-1}},-1};
        //Front Cell {10,11,12,13}
        corrFace old0Front=   {FType::old,{{VType::cubePnt,5}},0};
        corrFace old5 =       {FType::old,{},5};
        corrFace old1 =       {FType::old,{},1};
        corrFace split2Front= {FType::splitOld,{{VType::cubePnt,6},{VType::cutEdge,10},{VType::cutEdge,6}},2};
        //Middle Cell {14}
        corrFace split2Middle={FType::splitOld,{{VType::cutEdge,6},{VType::cutEdge,10},{VType::cutEdge,2},{VType::cutEdge,11}},2};
        //Back Cell {15,16,17,18}
        corrFace old0Back =   {FType::old,{{VType::cubePnt,0}},2};
        corrFace old4 =       {FType::old,{},4};
        corrFace old3 =       {FType::old,{},3};
        corrFace split2Back = {FType::splitOld,{{VType::cutEdge,3},{VType::cutEdge,11},{VType::cutEdge,2}},2};
        corrConcaveCellData redMarkSide = {{{0,1,7,8,10,11,12,13},
                                            {2,6,8,9,14},
                                            {3,4,5,9,15,16,17,18}}, 
                    {mc33Tri49,mc33Tri9A,mc33TriA2,mc33Tri20,mc33Tri08,mc33Tri8B,mc33TriB6,mc33Tri64,
                     addedRedMark6A,addedRedMarkB2,
                     old0Front,old5,old1,split2Front,split2Middle,old0Back,old4,old3,split2Back
                    }};
        
    // RedMarkNotPlus
        corrFace addedRedMark48={FType::added,{{VType::cutEdge,4},{VType::cutEdge,8},{VType::centerPnt,-1}},-1};//8
        corrFace addedRedMark09={FType::added,{{VType::cutEdge,0},{VType::cutEdge,9},{VType::centerPnt,-1}},-1};
        //Left Cell {10,11,12,13}
        corrFace old2Left =   {FType::old,{{VType::cubePnt,7}},2};
        corrFace old5 =       {FType::old,{},5};
        corrFace old3 =       {FType::old,{},3};
        corrFace split0Left = {FType::splitOld,{{VType::cubePnt,4},{VType::cutEdge,4},{VType::cutEdge,8}},0};
        //Middle Cell {14}
        corrFace split0Middle={FType::splitOld,{{VType::cutEdge,4},{VType::cutEdge,9},{VType::cutEdge,0},{VType::cutEdge,8}},0};
        //Right Cell {15,16,17,18}
        corrFace old2Right =  {FType::old,{{VType::cubePnt,2}},2};
        corrFace old1 =       {FType::old,{},1};
        corrFace old4 =       {FType::old,{},4};
        corrFace split0Right= {FType::splitOld,{{VType::cubePnt,1},{VType::cutEdge,0},{VType::cutEdge,9}},0};
        
        corrConcaveCellData nonRedMarkSide = {{{5,6,7,8,10,11,12,13},
                                               {0,4,8,9,14},
                                               {1,2,3,9,15,16,17,18}}, 
                    {mc33Tri49,mc33Tri9A,mc33TriA2,mc33Tri20,mc33Tri08,mc33Tri8B,mc33TriB6,mc33Tri64,
                     addedRedMark48,addedRedMark09,
                     old2Left,old5,old3,split0Left,split0Middle,old2Right,old1,old4,split0Right
                    }};
                    
    convexCorrectionData[MC33::Case::c102_1286] = {redMarkSide,nonRedMarkSide};
}
// Case c11
//---

// Case c1211
{
    corrFace mc33TriE513 =  {FType::mc33Triangle,{{VType::cutEdge,5},{VType::cutEdge,1},{VType::cutEdge,3}},-2};
    corrFace mc33TriE538 =  {FType::mc33Triangle,{{VType::cutEdge,5},{VType::cutEdge,3},{VType::cutEdge,8}},-2};
    corrFace mc33TriE845 =  {FType::mc33Triangle,{{VType::cutEdge,8},{VType::cutEdge,4},{VType::cutEdge,5}},-2};
    corrFace mc33TriE67B =  {FType::mc33Triangle,{{VType::cutEdge,6},{VType::cutEdge,7},{VType::cutEdge,11}},-2};
    // RedMarkPlus
        corrFace addedRedMark = {FType::added,{{VType::cutEdge,5},{VType::cubePnt,5},{VType::cubePnt,0},{VType::cutEdge,3}},-1};
        //Right Cell {4,5,6}
        corrFace old1 =      {FType::old,{},1};
        corrFace old4 =      {FType::old,{},4};
        corrFace split0Right={FType::splitOld,{{VType::cubePnt,5},{VType::cubePnt,1},{VType::cubePnt,0}},0};
        //Left Cell {7,8,9}
        corrFace old5 =      {FType::old,{},5};
        corrFace old3 =      {FType::old,{},3};
        corrFace split0Left ={FType::splitOld,{{VType::cubePnt,0},{VType::cubePnt,5},{VType::cutEdge,4},{VType::cutEdge,8}},0};
        corrConcaveCellData redMarkSide =   {{{0,3,4,5,6},{1,2,3,7,8,9}},
                                             {mc33TriE513,mc33TriE538,mc33TriE845,addedRedMark,
                                              old1,old4,split0Right,old5,old3,split0Left}
                                            };        
    // RedMarkNotPlus        
        corrFace addedNonRedMark = {FType::added,{{VType::cutEdge,5},{VType::cubePnt,6},{VType::cubePnt,3},{VType::cutEdge,3}},-1};
        //Right Cell {4,5,6}
        old1 =      {FType::old,{},1};
        old4 =      {FType::old,{},4};
        corrFace split2Right={FType::splitOld,{{VType::cubePnt,6},{VType::cubePnt,2},{VType::cubePnt,3}},2};
        //Left Cell {7,8,9,10}
        old3 =       {FType::old,{},3};
        old5 =       {FType::old,{},5};
        corrFace old0 =     {FType::old,{},0};
        corrFace split2Left={FType::splitOld,{{VType::cubePnt,6},{VType::cubePnt,3},{VType::cutEdge,11},{VType::cutEdge,6}},2};
        corrConcaveCellData nonRedMarkSide = {{{0,3,4,5,6},{1,2,3,7,8,9,10,11}},
                                              {mc33TriE513,mc33TriE538,mc33TriE845,addedRedMark,
                                               old1,old4,split2Right,old3,old5,old0,split2Left,mc33TriE67B}
                                             };
        
    convexCorrectionData[MC33::Case::c1211] = {redMarkSide,nonRedMarkSide};
}
// Case c1212
{
    corrFace mc33TriE156 =  {FType::mc33Triangle,{{VType::cutEdge,1},{VType::cutEdge,5},{VType::cutEdge,6}},-2};
    corrFace mc33TriE16B =  {FType::mc33Triangle,{{VType::cutEdge,1},{VType::cutEdge,6},{VType::cutEdge,11}},-2};
    corrFace mc33TriE13B =  {FType::mc33Triangle,{{VType::cutEdge,1},{VType::cutEdge,3},{VType::cutEdge,11}},-2};
    // RedMarkPlus
        // --
    // RedMarkNotPlus        
        corrFace addedNonRedMark621 = {FType::added,{{VType::cutEdge,6},{VType::cubePnt,2},{VType::cutEdge,1}},-1};//3
        corrFace addedNonRedMarkB21 = {FType::added,{{VType::cutEdge,11},{VType::cubePnt,2},{VType::cutEdge,1}},-1};
        //Right Cell {5,6,7}
        corrFace old1 =      {FType::old,{},1};
        corrFace old5Up =    {FType::old,{{VType::cubePnt,6}},5};
        corrFace split2Front={FType::splitOld,{{VType::cubePnt,6},{VType::cubePnt,2},{VType::cutEdge,6}},2};
        //Middle Cell {8}
        corrFace split2Middle={FType::splitOld,{{VType::cutEdge,6},{VType::cutEdge,11},{VType::cubePnt,2}},2};
        //Left Cell {9,10,11}
        corrFace old4 =      {FType::old,{},4};
        corrFace old3Up =    {FType::old,{{VType::cubePnt,3}},3};
        corrFace split2Back= {FType::splitOld,{{VType::cubePnt,2},{VType::cutEdge,11},{VType::cubePnt,3}},2};
        
        corrConcaveCellData nonRedMarkSide = {{{0,3,5,6,7},{1,3,4,8},{2,4,9,10,11}},
                                              {mc33TriE156,mc33TriE16B,mc33TriE13B,
                                               addedNonRedMark621,addedNonRedMarkB21,
                                               old1,old5Up,split2Front,split2Middle,old4,old3Up,split2Back}
                                             };
        
    convexCorrectionData[MC33::Case::c1212] = {{{{}},{}},nonRedMarkSide};
}
// Case c122
{
    corrFace mc33Tri65  =   {FType::mc33Triangle,{{VType::cutEdge,4}, {VType::cutEdge,9}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri51  =   {FType::mc33Triangle,{{VType::cutEdge,9}, {VType::cutEdge,10},{VType::centerPnt,-1}},-2};
    corrFace mc33Tri13  =   {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,2}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri38  =   {FType::mc33Triangle,{{VType::cutEdge,2}, {VType::cutEdge,0}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri84  =   {FType::mc33Triangle,{{VType::cutEdge,0}, {VType::cutEdge,8}, {VType::centerPnt,-1}},-2};//4
    corrFace mc33Tri47  =   {FType::mc33Triangle,{{VType::cutEdge,8}, {VType::cutEdge,11},{VType::centerPnt,-1}},-2};
    corrFace mc33Tri7B  =   {FType::mc33Triangle,{{VType::cutEdge,11},{VType::cutEdge,6}, {VType::centerPnt,-1}},-2};
    corrFace mc33TriB6  =   {FType::mc33Triangle,{{VType::cutEdge,6}, {VType::cutEdge,4}, {VType::centerPnt,-1}},-2};

    // RedMarkPlus
        corrFace addedRedMark45={FType::added,{{VType::cutEdge,4},{VType::cutEdge,5},{VType::centerPnt,-1}},-1};//8
        //Front Cell {9,10,11}
        corrFace old2    =    {FType::old,{},2};
        corrFace old3Front=   {FType::old,{{VType::cubePnt,7}},3};
        corrFace split5Left = {FType::splitOld,{{VType::cubePnt,7},{VType::cutEdge,7},{VType::cutEdge,4},{VType::cutEdge,5},{VType::cutEdge,6}},5};
        //Right Cell {15,16,17,18,19}
        corrFace old0 =       {FType::old,{},0};
        corrFace old1 =       {FType::old,{},1};
        corrFace old3Back =   {FType::old,{{VType::cubePnt,0}},3};
        corrFace old4 =       {FType::old,{},4};
        corrFace split5Right= {FType::splitOld,{{VType::cubePnt,5},{VType::cutEdge,5},{VType::cutEdge,4}},5};
        
        corrConcaveCellData redMarkSide = {{{5,6,7,0,8,9,10,11},{1,2,3,4,8,15,16,17,18,19}}, 
                    {mc33Tri65,mc33Tri51,mc33Tri13,mc33Tri38,mc33Tri84,mc33Tri47,mc33Tri7B,mc33TriB6,
                     addedRedMark45,
                     old2,old3Front,split5Left,old0,old1,old3Back,old4,split5Right
                    }};
         
    // RedMarkNotPlus
        corrFace addedNonRedMark3B={FType::added,{{VType::cutEdge,3},{VType::cutEdge,11},{VType::centerPnt,-1}},-1};//8
        //Front Cell {9,10,11,12,13}
        old1 =                {FType::old,{},1};
        old2 =                {FType::old,{},2};
        old4 =                {FType::old,{},4};
        corrFace old5Top  =   {FType::splitOld,{{VType::cubePnt,6}},5};
        corrFace split3Top =  {FType::splitOld,{{VType::cubePnt,3},{VType::cutEdge,3},{VType::cutEdge,11}},3};
        //Back Cell {15,16,17}
        old0 =                {FType::old,{},0};
        corrFace old5Bottom = {FType::old,{{VType::cubePnt,4}},5};
        corrFace split3Bottom={FType::splitOld,{{VType::cubePnt,4},{VType::cutEdge,7},{VType::cutEdge,11},{VType::cutEdge,3},{VType::cutEdge,8}},3};
        
        corrConcaveCellData nonRedMarkSide = {{{7,0,1,2,8,9,10,11,12,13},{3,4,5,6,8,15,16,17}}, 
                    {mc33Tri65,mc33Tri51,mc33Tri13,mc33Tri38,mc33Tri84,mc33Tri47,mc33Tri7B,mc33TriB6,
                     addedNonRedMark3B,
                     old1,old2,old4,old5Top,split3Top,old0,old5Bottom,split3Bottom
                    }};
        
    convexCorrectionData[MC33::Case::c122] = {redMarkSide,nonRedMarkSide};
}
// Case c131
//---

// Case c132
{
    corrFace mc33TriLeftLeft = {FType::mc33Triangle,{{VType::cutEdge,6},{VType::cutEdge,5},{VType::cutEdge,11}},-2};
    corrFace mc33TriLeft     = {FType::mc33Triangle,{{VType::cutEdge,5},{VType::cutEdge,3},{VType::cutEdge,11}},-2};
    corrFace mc33TriRight =    {FType::mc33Triangle,{{VType::cutEdge,5},{VType::cutEdge,3},{VType::cutEdge,10}},-2};
    corrFace mc33TriRightRight={FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,2},{VType::cutEdge,3}},-2};
    
    corrFace mc33TriOtherLeft= {FType::mc33Triangle,{{VType::cutEdge,4},{VType::cutEdge,7},{VType::cutEdge,8}},-2};
    corrFace mc33TriOtherRight={FType::mc33Triangle,{{VType::cutEdge,0},{VType::cutEdge,1},{VType::cutEdge,9}},-2};
    // RedMarkPlus
        //---
    // RedMarkNotPlus
        corrFace addedNonRedMark={FType::added,{{VType::cutEdge,5},{VType::cubePnt,5},{VType::cubePnt,0},{VType::cutEdge,3}},-1};
        //Right Cell {7,8,9,10}
        corrFace old1 =       {FType::old,{},1};
        corrFace old4 =       {FType::old,{},4};
        corrFace split0Right= {FType::splitOld,{{VType::cutEdge,9},{VType::cutEdge,0},{VType::cubePnt,0},{VType::cubePnt,5}},0};
        corrFace old2Right =  {FType::old,{{VType::cubePnt,2}},2};
        //Left Cell {11,12,13,14}
        corrFace old3 =       {FType::old,{},3};
        corrFace old5 =       {FType::old,{},5};
        corrFace split0Left = {FType::splitOld,{{VType::cubePnt,5},{VType::cubePnt,0},{VType::cutEdge,8},{VType::cutEdge,4}},0};
        corrFace old2Left  =  {FType::old,{{VType::cubePnt,7}},2};
        
        corrConcaveCellData nonRedMarkSide = {{{0,1,4,6,7,8,9,10},{2,3,5,11,12,13,14}},{mc33TriLeftLeft,mc33TriLeft,mc33TriRight,mc33TriRightRight,mc33TriOtherLeft,mc33TriOtherRight,
         addedNonRedMark,old1,old4,split0Right,old2Right,old3,old5,split0Left,old2Left}};
        
    convexCorrectionData[MC33::Case::c132] = {{{{}},{}},nonRedMarkSide};
}
// Case c133
{
    corrFace mc33Tri65  =   {FType::mc33Triangle,{{VType::cutEdge,6}, {VType::cutEdge,5}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri5A  =   {FType::mc33Triangle,{{VType::cutEdge,5}, {VType::cutEdge,10},{VType::centerPnt,-1}},-2};
    corrFace mc33TriA2  =   {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,2}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri23  =   {FType::mc33Triangle,{{VType::cutEdge,2}, {VType::cutEdge,3}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri38  =   {FType::mc33Triangle,{{VType::cutEdge,3}, {VType::cutEdge,8}, {VType::centerPnt,-1}},-2};//4
    corrFace mc33Tri84  =   {FType::mc33Triangle,{{VType::cutEdge,8}, {VType::cutEdge,4}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri47  =   {FType::mc33Triangle,{{VType::cutEdge,4}, {VType::cutEdge,7}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri7B  =   {FType::mc33Triangle,{{VType::cutEdge,7}, {VType::cutEdge,11},{VType::centerPnt,-1}},-2};
    corrFace mc33TriB6  =   {FType::mc33Triangle,{{VType::cutEdge,11},{VType::cutEdge,6}, {VType::centerPnt,-1}},-2};//8
    
    corrFace mc33Tri019  =  {FType::mc33Triangle,{{VType::cutEdge,0}, {VType::cutEdge,1}, {VType::cutEdge,9}},-2};

    // RedMarkPlus
        corrFace addedNonRedMarkB2={FType::added,{{VType::cutEdge,11},{VType::cutEdge,2},{VType::centerPnt,-1}},-1};//10
        corrFace addedNonRedMarkB3={FType::added,{{VType::cutEdge,11},{VType::cutEdge,5},{VType::centerPnt,-1}},-1};
        //Left Center Cell {12,13,14}
        corrFace old5Left  =  {FType::old,{{VType::cubePnt,4}},5};
        corrFace old0Left =   {FType::old,{{VType::cubePnt,4}},0};
        corrFace split3Front= {FType::splitOld,{{VType::cubePnt,4},{VType::cutEdge,7},{VType::cutEdge,11},{VType::cutEdge,3},{VType::cutEdge,8}},3};
        //Middle Upper Cell {15,16,17}
        corrFace old4Up =     {FType::old,{{VType::cubePnt,3}},4};        
        corrFace split2Back = {FType::splitOld,{{VType::cubePnt,3},{VType::cutEdge,11},{VType::cutEdge,2}},2};
        corrFace split3Back = {FType::splitOld,{{VType::cubePnt,3},{VType::cutEdge,3},{VType::cutEdge,11}},3};
        //Right Upper Cell {18,19,20}
        corrFace old5Right  = {FType::old,{{VType::cubePnt,6}},5};
        corrFace old1Up   =   {FType::old,{{VType::cubePnt,6}},1};
        corrFace split2Front= {FType::splitOld,{{VType::cutEdge,11},{VType::cutEdge,2},{VType::cutEdge,10},{VType::cubePnt,6},{VType::cutEdge,6}},2};
        
        corrConcaveCellData redMarkSide = {{{4,5,6,7,11,12,13,14},{3,10,11,15,16,17},                 
                                            {8,0,1,2,10,18,19,20}}, 
                    {mc33Tri65,mc33Tri5A,mc33TriA2,mc33Tri23,mc33Tri38,mc33Tri84,mc33Tri47,mc33Tri7B,mc33TriB6,mc33Tri019,
                    addedNonRedMarkB2,addedNonRedMarkB3,
                    old5Left,old0Left,split3Front,old4Up,split2Back,split3Back,old5Right,old1Up,split2Front
                    }};
         
    // RedMarkNotPlus
        corrFace addedNonRedMark55={FType::added,{{VType::cutEdge,5},{VType::cubePnt,5},{VType::centerPnt,-1}},-1};//10
        corrFace addedNonRedMark45={FType::added,{{VType::cutEdge,4},{VType::cubePnt,5},{VType::centerPnt,-1}},-1};
        //Front Cell {12,13,14}
        corrFace old3Front =  {FType::old,{{VType::cubePnt,7}},3};
        corrFace old2Front =  {FType::old,{{VType::cubePnt,7}},2};
        corrFace old5 =       {FType::old,{},5};
        //Back Cell {15,16,17}
        corrFace old3Back  =  {FType::old,{{VType::cubePnt,0}},3};
        corrFace old2Back  =  {FType::old,{{VType::cubePnt,2}},2};
        corrFace old4 =       {FType::old,{},4};
        corrFace old1 =       {FType::old,{},1};
        corrFace old0 =       {FType::old,{},0};
        
        corrConcaveCellData nonRedMarkSide = {{{6,7,8,0,10,11,12,13,14},{1,2,3,4,5,10,11,15,16,17}},
                    {mc33Tri65,mc33Tri5A,mc33TriA2,mc33Tri23,mc33Tri38,mc33Tri84,mc33Tri47,mc33Tri7B,mc33TriB6,mc33Tri019,
                    addedNonRedMarkB2,addedNonRedMarkB3,
                    old3Front,old2Front,old5,old3Back,old2Back,old4,old1,old0
                    }};
        
    convexCorrectionData[MC33::Case::c133] = {redMarkSide,nonRedMarkSide};
}
// Case c134
{
    corrFace mc33Tri65  =   {FType::mc33Triangle,{{VType::cutEdge,6}, {VType::cutEdge,5}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri5A  =   {FType::mc33Triangle,{{VType::cutEdge,5}, {VType::cutEdge,10},{VType::centerPnt,-1}},-2};
    corrFace mc33TriA2  =   {FType::mc33Triangle,{{VType::cutEdge,10},{VType::cutEdge,2}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri21  =   {FType::mc33Triangle,{{VType::cutEdge,2}, {VType::cutEdge,1}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri19  =   {FType::mc33Triangle,{{VType::cutEdge,1}, {VType::cutEdge,9}, {VType::centerPnt,-1}},-2};//4
    corrFace mc33Tri90  =   {FType::mc33Triangle,{{VType::cutEdge,9}, {VType::cutEdge,0}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri03  =   {FType::mc33Triangle,{{VType::cutEdge,0}, {VType::cutEdge,3}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri38  =   {FType::mc33Triangle,{{VType::cutEdge,3}, {VType::cutEdge,8}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri84  =   {FType::mc33Triangle,{{VType::cutEdge,8}, {VType::cutEdge,4}, {VType::centerPnt,-1}},-2};//8
    corrFace mc33Tri47  =   {FType::mc33Triangle,{{VType::cutEdge,4}, {VType::cutEdge,7}, {VType::centerPnt,-1}},-2};
    corrFace mc33Tri7B  =   {FType::mc33Triangle,{{VType::cutEdge,7}, {VType::cutEdge,11},{VType::centerPnt,-1}},-2};
    corrFace mc33TriB6  =   {FType::mc33Triangle,{{VType::cutEdge,11},{VType::cutEdge,6}, {VType::centerPnt,-1}},-2};
    
    // RedMarkPlus
        corrFace addedNonRedMark23={FType::added,{{VType::cutEdge,5},{VType::cubePnt,5},{VType::centerPnt,-1}},-1};//12
        corrFace addedNonRedMark33={FType::added,{{VType::cutEdge,4},{VType::cubePnt,5},{VType::centerPnt,-1}},-1};
        corrFace addedNonRedMarkB3={FType::added,{{VType::cutEdge,9},{VType::cubePnt,5},{VType::centerPnt,-1}},-1};
        //Left Center Cell {15,16,17}
        corrFace old5Left  =  {FType::old,{{VType::cubePnt,4}},5};
        corrFace old0Left =   {FType::old,{{VType::cubePnt,4}},0};
        corrFace old3 =       {FType::old,{},5};
        //Middle Lower Cell {18,19,20}
        corrFace old2 =       {FType::old,{},2};
        corrFace old5Right =  {FType::old,{{VType::cubePnt,6}},5};
        corrFace old1Front =  {FType::old,{{VType::cubePnt,6}},1};
        //Back Center Cell {21,22,23}
        corrFace old4 =       {FType::old,{},4};
        corrFace old0Right  = {FType::old,{{VType::cubePnt,1}},0};
        corrFace old1Back  =  {FType::old,{{VType::cubePnt,1}},1};
        
        corrConcaveCellData redMarkSide = {{{7,8,9,10,14,13,15,16,17},{11,0,1,2,12,14,18,19,20},                 
                                            {3,4,5,6,12,13,21,22,23}}, 
                    {mc33Tri65,mc33Tri5A,mc33TriA2,mc33Tri21,mc33Tri19,mc33Tri90,mc33Tri03,mc33Tri38,mc33Tri84,mc33Tri47,mc33Tri7B,mc33TriB6,
                    addedNonRedMark23,addedNonRedMark33,addedNonRedMarkB3,
                    old5Left,old0Left,old3,old2,old5Right,old1Front,old4,old0Right,old1Back
                    }};
         
    // RedMarkNotPlus
        corrFace addedNonRedMark55={FType::added,{{VType::cutEdge,5},{VType::cubePnt,5},{VType::centerPnt,-1}},-1};//12
        corrFace addedNonRedMark45={FType::added,{{VType::cutEdge,4},{VType::cubePnt,5},{VType::centerPnt,-1}},-1};
        corrFace addedNonRedMark95={FType::added,{{VType::cutEdge,9},{VType::cubePnt,5},{VType::centerPnt,-1}},-1};
        //Front Center Cell {15,16,17}
        corrFace old3Front =  {FType::old,{{VType::cubePnt,7}},3};
        corrFace old2Front =  {FType::old,{{VType::cubePnt,7}},2};
        corrFace old5 =       {FType::old,{},5};
        //Middle Lower Cell {18,19,20}
        corrFace old0 =       {FType::old,{},0};
        corrFace old3Back  =  {FType::old,{{VType::cubePnt,0}},3};
        corrFace old4Left  =  {FType::old,{{VType::cubePnt,0}},4};
        //Right Center Cell {21,22,23}
        corrFace old1 =       {FType::old,{},1};
        corrFace old2Back  =  {FType::old,{{VType::cubePnt,2}},2};
        corrFace old4Right =  {FType::old,{{VType::cubePnt,2}},4};
        
        corrConcaveCellData nonRedMarkSide = {{{9,10,11,0,12,13,15,16,17},{5,6,7,8,13,14,18,19,20},                 
                                               {1,2,3,4,12,14,21,22,23}}, 
                    {mc33Tri65,mc33Tri5A,mc33TriA2,mc33Tri21,mc33Tri19,mc33Tri90,mc33Tri03,mc33Tri38,mc33Tri84,mc33Tri47,mc33Tri7B,mc33TriB6,
                    addedNonRedMark55,addedNonRedMark45,addedNonRedMark95,
                    old3Front,old2Front,old5,old0,old3Back,old4Left,old1,old2Back,old4Right
                    }};
        
    convexCorrectionData[MC33::Case::c134] = {redMarkSide,nonRedMarkSide};
}
// Case c1351
//---

// Case c1352
//---

// Case c14
{
    corrFace mc33Tri12B = {FType::mc33Triangle,{{VType::cutEdge,1},{VType::cutEdge,2},{VType::cutEdge,11}},-2};
    corrFace mc33Tri14B = {FType::mc33Triangle,{{VType::cutEdge,1},{VType::cutEdge,4},{VType::cutEdge,11}},-2};
    corrFace mc33Tri145 = {FType::mc33Triangle,{{VType::cutEdge,1},{VType::cutEdge,4},{VType::cutEdge,5}},-2};
    corrFace mc33Tri48B = {FType::mc33Triangle,{{VType::cutEdge,4},{VType::cutEdge,8},{VType::cutEdge,11}},-2};

    // RedMarkPlus
        corrFace addedRedMark={FType::added,{{VType::cutEdge,1},{VType::cubePnt,1},{VType::cutEdge,4}},-1};
        //Front Cell {5,6,7}
        corrFace old1 =       {FType::old,{},1};
        corrFace old5 =       {FType::old,{},5};
        corrFace split0Front ={FType::splitOld,{{VType::cutEdge,4},{VType::cubePnt,1},{VType::cubePnt,5}},0};
        //Back Cell {8,9,10,11}
        corrFace old2 =       {FType::old,{},2};
        corrFace old3 =       {FType::old,{},3};
        corrFace old4 =       {FType::old,{},4};
        corrFace split0Back = {FType::splitOld,{{VType::cutEdge,4},{VType::cubePnt,1},{VType::cubePnt,0},{VType::cutEdge,8}},0};
        
        corrConcaveCellData redMarkSide = {{{2,4,5,6,7},{0,1,3,4,8,9,10,11}},{mc33Tri12B,mc33Tri14B,mc33Tri145,mc33Tri48B,addedRedMark,
         old1,old5,split0Front,old2,old3,old4,split0Back}};
        
    // RedMarkNotPlus
        corrFace addedNonRedMark12B={FType::added,{{VType::cutEdge,1},{VType::cubePnt,2},{VType::cutEdge,11},},-1};
        corrFace addedNonRedMark47B={FType::added,{{VType::cutEdge,4},{VType::cubePnt,7},{VType::cutEdge,11},},-1};
        //Right Cell {6,7}
        old4 =                  {FType::old,{},4};
        corrFace split2Back =   {FType::splitOld,{{VType::cutEdge,2},{VType::cutEdge,11},{VType::cubePnt,2}},2};
        //Center Cell {8,9}
        corrFace split2Front =  {FType::splitOld,{{VType::cutEdge,11},{VType::cubePnt,7},{VType::cubePnt,6},{VType::cubePnt,2}},2};
        corrFace split5Right =  {FType::splitOld,{{VType::cutEdge,4},{VType::cutEdge,5},{VType::cubePnt,6},{VType::cubePnt,7}},5};
        //Left Cell {10,11}
        old3 =                  {FType::old,{},3};
        corrFace old0 =         {FType::old,{},0};
        corrFace split5Left =   {FType::splitOld,{{VType::cutEdge,4},{VType::cubePnt,7},{VType::cubePnt,4}},5};
        
        corrConcaveCellData nonRedMarkSide = {{{0,4,6,7},{1,2,4,5,8,9},{3,5,10,11}},
                                    {mc33Tri12B,mc33Tri14B,mc33Tri145,mc33Tri48B,
                                     addedNonRedMark12B,addedNonRedMark47B,
                                     old4,split2Back,split2Front,split5Right,old3,old0,split5Left}};    
            
    convexCorrectionData[MC33::Case::c14] = {redMarkSide,nonRedMarkSide};
}
};
