// Permutation of a Cube given in MC33_LookUpTable.h
{
  {0,1,2,3,4,5,6,7},//  Top Face 2: Rot E
  {1,5,6,2,0,4,7,3},//              Rot Rz*E
  {5,4,7,6,1,0,3,2},//              Rot Rz²*E
  {4,0,3,7,5,1,2,6},//              Rot Rz³*E

  {4,5,1,0,7,6,2,3},//  Top Face 4: Rot Ry
  {5,6,2,1,4,7,3,0},//              Rot Rz*Ry
  {6,7,3,2,5,4,0,1},//              Rot Rz²*Ry
  {7,4,0,3,6,5,1,2},//              Rot Rz³*Ry
  
  {7,6,5,4,3,2,1,0},//  Top Face 0: Rot Ry²
  {6,2,1,5,7,3,0,4},//              Rot Rz*Ry²
  {2,3,0,1,6,7,4,5},//              Rot Rz²*Ry²
  {3,7,4,0,2,6,5,1},//              Rot Rz³*Ry²

  {3,2,6,7,0,1,5,4},//  Top Face 5: Rot Ry³
  {2,1,5,6,3,0,4,7},//              Rot Rz*Ry³
  {1,0,4,5,2,3,7,6},//              Rot Rz²*Ry³
  {0,3,7,4,1,2,6,5},//              Rot Rz³*Ry³
  
  {3,0,1,2,7,4,5,6},//  Top Face 1: Rot Rx
  {0,4,5,1,3,7,6,2},//              Rot Rz*Rx
  {4,7,6,5,0,3,2,1},//              Rot Rz²*Rx
  {7,3,2,6,4,0,1,5},//              Rot Rz³*Rx
  
  {1,2,3,0,5,6,7,4},//  Top Face 3: Rot -Rx
  {2,6,7,3,1,5,4,0},//              Rot Rz*-Rx
  {6,5,4,7,2,1,0,3},//              Rot Rz²*-Rx
  {5,1,0,4,6,2,3,7} //              Rot Rz³*-Rx
}
