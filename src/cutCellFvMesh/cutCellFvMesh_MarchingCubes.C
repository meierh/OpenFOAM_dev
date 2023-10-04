#include "cutCellFvMesh.H"

inline float invSqrt(float f) {
	return 1.0/sqrt(f);
}

inline unsigned int permuteBitPattern(std::array<unsigned short int,8>& permutation, unsigned int bitPattern)
{
	unsigned int resPattern = 0;
	for(int vertI=0;vertI<permutation.size();vertI++)
	{
		bool bitSet = getBit(bitPattern,7-permutation[vertI]);
		if(bitSet)
			setBit(resPattern,7-vertI);
	}
	return resPattern;
}

void Foam::cutCellFvMesh::MC33::setBit(unsigned short int& bitField, std::uint8_t bitIndex)
{
    if(bitIndex>15)
        FatalErrorInFunction<<"Bit Index must be inside [0,15]"<< exit(FatalError);

    unsigned short int mask = 1;
    mask = mask << bitIndex;
    
    bitField = bitField | mask;
}

bool Foam::cutCellFvMesh::MC33::getBit(const unsigned short int bitField, std::uint8_t bitIndex)
{
    if(bitIndex>15)
        FatalErrorInFunction<<"Bit Index must be inside [0,15]"<< exit(FatalError);

    unsigned short int mask = 1;
    mask = mask << bitIndex;
    
    return bitField & mask;
}

std::string Foam::cutCellFvMesh::MC33::getBitPattern(const unsigned short int x)
{
	std::string pattern;
	for(int i=15;i>=0;i--)
	{
		pattern.append(std::to_string(getBit(x,i)));
		if(i%4==0)
			pattern.append(" ");
	}
	return pattern;
}

Foam::cutCellFvMesh::MC33::MC33
(
    cutCellFvMesh& mesh
):
mesh(mesh),
table(mesh.caseTable),
permutations(mesh.permutationTable)
{
}

Foam::cutCellFvMesh::MC33::MC33Cube Foam::cutCellFvMesh::MC33::generateMC33Cube
(
	int cellInd
)
{
	//Info<<Foam::endl<<cellInd<<Foam::endl;
	const cellList& meshCells = mesh.cells();
	const faceList& meshFaces = mesh.faces();
    const cellShapeList& meshShapes = mesh.cellShapes();
	
	/*
	Info<<"meshCells points:"<<meshCells[cellInd].labels(meshFaces)<<Foam::endl<<"(";
	for(label pnt : meshCells[cellInd].labels(meshFaces))
		Info<<mesh.pointDist[pnt]<<" ";
	Info<<")"<<Foam::endl;
	*/
	
	cell thisCell = meshCells[cellInd];
	edgeList thisEdges = thisCell.edges(meshFaces);
	cellShape thisCellShape = meshShapes[cellInd];
    cellModel thisCellModel = thisCellShape.model();
	
	//Info<<"meshFaces points:"<<meshFaces[thisCell[0]]<<Foam::endl;
	
	bool isHex = true;
	isHex = isHex && (thisCell.nFaces()==6);
	labelList vertices = thisCell.labels(meshFaces);
	isHex = isHex && (vertices.size()==8);
	edgeList edges = thisCell.edges(meshFaces);
	isHex = isHex && (edges.size()==12);
	isHex = isHex && (thisCellModel.name()=="hex");
	isHex = isHex && (thisCell.size()==6);
	
	if(!isHex)
	{
		Info<<thisCellModel.name()<<Foam::endl;
		FatalErrorInFunction<<"Can not happen!"<< exit(FatalError);
	}

	
	MC33Cube oneCube;
	
	//Info<<"Vertices:"<<oneCube.vertices<<Foam::endl;
	
	oneCube.cell = cellInd;
	oneCube.vertices[0] = meshFaces[thisCell[0]][0];
	oneCube.vertices[1] = meshFaces[thisCell[0]][1];
	oneCube.vertices[2] = meshFaces[thisCell[0]][2];
	oneCube.vertices[3] = meshFaces[thisCell[0]][3];
	oneCube.edges[0] = edge(oneCube.vertices[0],oneCube.vertices[1]);
	oneCube.edges[1] = edge(oneCube.vertices[1],oneCube.vertices[2]);
	oneCube.edges[2] = edge(oneCube.vertices[2],oneCube.vertices[3]);
	oneCube.edges[3] = edge(oneCube.vertices[3],oneCube.vertices[0]);
	//oneCube.faces[4] = thisCell[0];
	
	std::unordered_set<label> startingFaceVert = {oneCube.vertices[0],oneCube.vertices[1],
												  oneCube.vertices[2],oneCube.vertices[3]};
	std::unordered_set<label> oppositeFaceVert;
	for(label vert:vertices)
		if(startingFaceVert.find(vert)!=startingFaceVert.end())
			oppositeFaceVert.insert(vert);
	
	//Info<<"Vertices:"<<oneCube.vertices<<Foam::endl;
	
	//Info<<"Edges:"<<thisEdges<<Foam::endl;	
	for(label vertI=0;vertI<4;vertI++)
	{
		label seenFaceVertice = oneCube.vertices[vertI];
		label oppositeVertice = -1;
		for(edge edg : thisEdges)
		{
			label edgOtherVert = edg.otherVertex(seenFaceVertice);
			if(edgOtherVert!=-1 && startingFaceVert.find(edgOtherVert)==startingFaceVert.end())
			{
				if(oppositeVertice!=-1)
				{
					
					Info<<Foam::endl;
					Info<<"vertI:"<<vertI<<Foam::endl;
					Info<<"oppositeVertice:"<<oppositeVertice<<Foam::endl;
					Info<<"edgOtherVert:"<<edgOtherVert<<Foam::endl;
					Info<<"Vertices:"<<oneCube.vertices<<Foam::endl;
					
					FatalErrorInFunction<<"Double assignment!"<< exit(FatalError);
				}
				else
					oppositeVertice = edgOtherVert;
			}
		}
		oneCube.vertices[4+vertI] = oppositeVertice;
		oneCube.edges[8+vertI] = edge(oneCube.vertices[vertI],oneCube.vertices[4+vertI]);
		//Info<<vertI<<" Vertices:"<<oneCube.vertices<<Foam::endl;
	}
	oneCube.edges[4] = edge(oneCube.vertices[4],oneCube.vertices[5]);
	oneCube.edges[5] = edge(oneCube.vertices[5],oneCube.vertices[6]);
	oneCube.edges[6] = edge(oneCube.vertices[6],oneCube.vertices[7]);
	oneCube.edges[7] = edge(oneCube.vertices[7],oneCube.vertices[4]);
	
	std::vector<std::unordered_set<label>> faceVertices(6);
	faceVertices[4] = startingFaceVert;
	faceVertices[5] = oppositeFaceVert;
	faceVertices[0] = {oneCube.vertices[0],oneCube.vertices[1],oneCube.vertices[4],oneCube.vertices[5]};
	faceVertices[1] = {oneCube.vertices[1],oneCube.vertices[2],oneCube.vertices[5],oneCube.vertices[6]};
	faceVertices[2] = {oneCube.vertices[2],oneCube.vertices[3],oneCube.vertices[6],oneCube.vertices[7]};
	faceVertices[3] = {oneCube.vertices[0],oneCube.vertices[3],oneCube.vertices[4],oneCube.vertices[7]};

	return oneCube;
}

unsigned int Foam::cutCellFvMesh::MC33::computeSignBitPattern
(
	const MC33Cube& cell
)
{
	unsigned int bitPattern = 0;
	for(int vertI=0;vertI<cell.vertices.size();vertI++)
	{
		label pntInd = cell.vertices[vertI];
		scalar iso = mesh.pointDist[pntInd];
		if(iso>=0)
		{
			unsigned int mask = 1;
			mask = mask << (7-vertI);
			bitPattern = bitPattern | mask;
		}
		v[vertI] = iso;
	}
	return bitPattern;
}

Foam::cutCellFvMesh::MC33::MC33Cube Foam::cutCellFvMesh::MC33::computeCutCell(int cellInd)
{
	MC33Cube mc33Cube = generateMC33Cube(cellInd);
	unsigned int bitPattern = computeSignBitPattern(mc33Cube);
	if(bitPattern!=0 && bitPattern!=255) 
	{
		mc33Cube.bitPattern = bitPattern;
		const unsigned short int* triangleCase = getTriangleCase(bitPattern);
		label arrayIndex = triangleCase-table;
		mc33Cube.cutTriangles = collectTriangles(triangleCase);
		
		std::vector<std::uint8_t> validPermutations;
		for(label permutInd=0;permutInd<permutations.size();permutInd++)
		{
			if(bitPattern==permuteBitPattern(permutations[permutInd],bitPattern))
				validPermutations.push_back(permutInd);
		}
	
		{
			if(arrayIndex<0)
				FatalErrorInFunction<<"ArrayIndex must not be smaller than zero!"<< exit(FatalError);
			else if(arrayIndex<128)
				FatalErrorInFunction<<"ArrayIndex not a case!"<< exit(FatalError);
			else if(arrayIndex<136)
				// Case 1 Triangle Number 1
				mc33Cube.cubeCase = c1;
			else if(arrayIndex<160)
				// Case 2 Triangle Number 2
				mc33Cube.cubeCase = c2;
			else if(arrayIndex<184)
				// Case 3.1 Triangle Number 2
				mc33Cube.cubeCase = c31;
			else if(arrayIndex<232)
				// Case 3.2 Triangle Number 4
				mc33Cube.cubeCase = c32;
			else if(arrayIndex<240)
				// Case 4.1.1 Triangle Number 2
				mc33Cube.cubeCase = c411;
			else if(arrayIndex<264)
				// Case 4.1.2 Triangle Number 6
				mc33Cube.cubeCase = c412;
			else if(arrayIndex<336)
				// Case 5 Triangle Number 3
				mc33Cube.cubeCase = c5;
			else if(arrayIndex<408)
				// Case 6.1.1 Triangle Number 3
				mc33Cube.cubeCase = c611;
			else if(arrayIndex<576)
				// Case 6.1.2 Triangle Number 7
				mc33Cube.cubeCase = c612;
			else if(arrayIndex<696)
				// Case 6.2 Triangle Number 5
				mc33Cube.cubeCase = c62;
			else if(arrayIndex<720)
				// Case 7.1 Triangle Number 3
				mc33Cube.cubeCase = c71;
			else if(arrayIndex<840)
				// Case 7.2 Triangle Number 5
				mc33Cube.cubeCase = c72;
			else if(arrayIndex<1056)
				// Case 7.3 Triangle Number 9
				mc33Cube.cubeCase = c73;
			else if(arrayIndex<1096)
				// Case 7.4.1 Triangle Number 5
				mc33Cube.cubeCase = c741;
			else if(arrayIndex<1168)
				// Case 7.4.2 Triangle Number 9
				mc33Cube.cubeCase = c742;
			else if(arrayIndex<1174)
				// Case 8 Triangle Number 2
				mc33Cube.cubeCase = c8;
			else if(arrayIndex<1190)
				// Case 9 Triangle Number 4
				mc33Cube.cubeCase = c9;
			else if(arrayIndex<1214)
				// Case 10.1.1 Triangle Number 4
				mc33Cube.cubeCase = c1011;
			else if(arrayIndex<1262)
				// Case 10.1.2 Triangle Number 8
				mc33Cube.cubeCase = c1012;
			else if(arrayIndex<1310)
				// Case 10.2 Triangle Number 8
				mc33Cube.cubeCase = c102;
			else if(arrayIndex<1334)
				// Case 11 Triangle Number 4
				mc33Cube.cubeCase = c11;
			else if(arrayIndex<1358)
				// Case 14 Triangle Number 4
				mc33Cube.cubeCase = c14;
			else if(arrayIndex<1454)
				// Case 12.1.1 Triangle Number 4
				mc33Cube.cubeCase = c1211;
			else if(arrayIndex<1646)
				// Case 12.1.2 Triangle Number 8
				mc33Cube.cubeCase = c1212;
			else if(arrayIndex<1838)
				// Case 12.2 Triangle Number 8
				mc33Cube.cubeCase = c122;
			else if(arrayIndex<1846)
				// Case 13.1 Triangle Number 4
				mc33Cube.cubeCase = c131;
			else if(arrayIndex<1918)
				// Case 13.2 Triangle Number 6
				mc33Cube.cubeCase = c132;
			else if(arrayIndex<2158)
				// Case 13.3 Triangle Number 10
				mc33Cube.cubeCase = c133;
			else if(arrayIndex<2206)
				// Case 13.4 Triangle Number 12
				mc33Cube.cubeCase = c134;
			else if(arrayIndex<2286)
				// Case 13.5.2 Triangle Number 10
				mc33Cube.cubeCase = c1352;
			else if(arrayIndex<2310)
				// Case 13.5.1 Triangle Number 6
				mc33Cube.cubeCase = c1351;
			else
				FatalErrorInFunction<<"ArrayIndex must not be larger than 2309!"<< exit(FatalError);
			
			if(mc33Cube.cubeCase==c1){
				// Case 1 Triangle Number 1
				if(mc33Cube.cutTriangles.size()!=1)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c2 ){
				// Case 2 Triangle Number 2
				if(mc33Cube.cutTriangles.size()!=2)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c31 ){
				// Case 3.1 Triangle Number 2
				if(mc33Cube.cutTriangles.size()!=2)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c32 ){
				// Case 3.2 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c411 ){
				// Case 4.1.1 Triangle Number 2
				if(mc33Cube.cutTriangles.size()!=2)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c412 ){
				// Case 4.1.2 Triangle Number 6
				if(mc33Cube.cutTriangles.size()!=6)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c5 ){
				// Case 5 Triangle Number 3
				if(mc33Cube.cutTriangles.size()!=3)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c611 ){
				// Case 6.1.1 Triangle Number 3
				if(mc33Cube.cutTriangles.size()!=3)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c612 ){
				// Case 6.1.2 Triangle Number 7
				if(mc33Cube.cutTriangles.size()!=7)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c62 ){
				// Case 6.2 Triangle Number 5
				if(mc33Cube.cutTriangles.size()!=5)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c71 ){
				// Case 7.1 Triangle Number 3
				if(mc33Cube.cutTriangles.size()!=3)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c72 ){
				// Case 7.2 Triangle Number 5
				if(mc33Cube.cutTriangles.size()!=5)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c73 ){
				// Case 7.3 Triangle Number 9
				if(mc33Cube.cutTriangles.size()!=9)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c741 ){
				// Case 7.4.1 Triangle Number 5
				if(mc33Cube.cutTriangles.size()!=5)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c742 ){
				// Case 7.4.2 Triangle Number 9
				if(mc33Cube.cutTriangles.size()!=9)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c8 ){
				// Case 8 Triangle Number 2
				if(mc33Cube.cutTriangles.size()!=2)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c9 ){
				// Case 9 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c1011 ){
				// Case 10.1.1 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c1012 ){
				// Case 10.1.2 Triangle Number 8
				if(mc33Cube.cutTriangles.size()!=8)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c102 ){
				// Case 10.2 Triangle Number 8
				if(mc33Cube.cutTriangles.size()!=8)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c11 ){
				// Case 11 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c14 ){
				// Case 14 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c1211 ){
				// Case 12.1.1 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c1212 ){
				// Case 12.1.2 Triangle Number 8
				if(mc33Cube.cutTriangles.size()!=8)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c122 ){
				// Case 12.2 Triangle Number 8
				if(mc33Cube.cutTriangles.size()!=8)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c131 ){
				// Case 13.1 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c132 ){
				// Case 13.2 Triangle Number 6
				if(mc33Cube.cutTriangles.size()!=6)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c133 ){
				// Case 13.3 Triangle Number 10
				if(mc33Cube.cutTriangles.size()!=10)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c134 ){
				// Case 13.4 Triangle Number 12
				if(mc33Cube.cutTriangles.size()!=12)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c1352){
				// Case 13.5.2 Triangle Number 10
				if(mc33Cube.cutTriangles.size()!=10)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else if(mc33Cube.cubeCase==c1351){
				// Case 13.5.1 Triangle Number 6
				if(mc33Cube.cutTriangles.size()!=6)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);}
			else
				FatalErrorInFunction<<"Error!"<< exit(FatalError);
		}
		
		/*
		const cellList& meshCells = mesh.cells();
		const faceList& meshFaces = mesh.faces();
		const labelList& cellPoints = meshCells[cellInd].labels(meshFaces);
		Info<<"cellInd:"<<cellInd<<Foam::endl<<"(";
		for(label vertice: cellPoints)
			Info<<" | "<<mesh.pointDist[vertice];
		Info<<")"<<Foam::endl;
		Info<<"bitPattern:"<<getBitPattern(bitPattern)<<Foam::endl;
		Info<<"Triangle number:"<<mc33Cube.cutTriangles.size()<<Foam::endl;
		*/
		
		return mc33Cube;
	}
	else
	{
		mc33Cube.cell=-1;
		return mc33Cube;		
	}
	
}

const unsigned short int* Foam::cutCellFvMesh::MC33::getTriangleCase(unsigned int verticePattern)
{
	union { // memory saving
		int f[6];//for the face tests
		double r[6];//for intercept and normal coordinates
	};
	const unsigned short int *pcase = table;

	unsigned int k, c, m, n;
	//double t;
	if (verticePattern&0x80) {
		c = pcase[verticePattern^0xFF];
		m = (c&0x800) == 0;
		n = !m;
	} else {
		c = pcase[verticePattern];
		n = (c&0x800) == 0;
		m = !n;
	}
	k = c&0x7FF;
	switch (c>>12) { //find the MC33 case
		case 0: // cases 1, 2, 5, 8, 9, 11 and 14
			pcase += k;
			break;
		case 1: // case 3
			pcase += ((m? verticePattern: verticePattern^0xFF)&face_test1(k>>2)? 183 + (k<<1): 159 + k);
			break;
		case 2: // case 4
			pcase += (interior_test(k,0)? 239 + 6*k: 231 + (k<<1));
			break;
		case 3: // case 6
			if ((m? verticePattern: verticePattern^0xFF)&face_test1(k%6))
				pcase += 575 + 5*k; //6.2
			else
				pcase += (interior_test(k/6,0)? 407 + 7*k: 335 + 3*k); //6.1
			break;
		case 4: // case 7
			switch (face_tests(f,(m? verticePattern: verticePattern^0xFF))) {
			case -3:
				pcase += 695 + 3*k; //7.1
				break;
			case -1: //7.2
				pcase += (f[4] + f[5] < 0? (f[0] + f[2] < 0? 759: 799): 719) + 5*k;
				break;
			case 1: //7.3
				pcase += (f[4] + f[5] < 0? 983: (f[0] + f[2] < 0? 839: 911)) + 9*k;
				break;
			default: //7.4
				pcase += (interior_test(k>>1,0)? 1095 + 9*k: 1055 + 5*k);
			}
			break;
		case 5: // case 10
			switch (face_tests(f,(m? verticePattern: verticePattern^0xFF))) {
			case -2:
				if (k&2? interior_test(0,0): interior_test(0,0)||interior_test(k? 1: 3,0))
					pcase += 1213 + (k<<3); //10.1.2
				else
					pcase += 1189 + (k<<2); //10.1.1
				break;
			case 0: //10.2
				pcase += (f[2 + k] < 0? 1261: 1285) + (k<<3);
				break;
			default:
				if (k&2? interior_test(1,0): interior_test(2,0)||interior_test(k? 3: 1,0))
					pcase += 1237 + (k<<3); //10.1.2
				else
					pcase += 1201 + (k<<2); //10.1.1
			}
			break;
		case 6: // case 12
			switch (face_tests(f,(m? verticePattern: verticePattern^0xFF))) {
			case -2: //12.1
				pcase += (interior_test((0xDA010C>>(k<<1))&3,0)? 1453 + (k<<3): 1357 + (k<<2));
				break;
			case 0: //12.2
				pcase += (f[k>>1] < 0? 1645: 1741) + (k<<3);
				break;
			default: //12.1
				pcase += (interior_test((0xA7B7E5>>(k<<1))&3,0)? 1549 + (k<<3): 1405 + (k<<2));
			}
			break;
		default: // case 13
			switch (abs(face_tests(f, 165))) {
			case 0:
				k = ((f[1] < 0)<<1)|(f[5] < 0);
				if (f[0]*f[1] == f[5]) //13.4
					pcase += 2157 + 12*k;
				else {
					c = interior_test(k, 1); // 13.5.1 if c == 0 else 13.5.2
					pcase += 2285 + (c? 10*k - 40*c: 6*k);
				}
				break;
			case 2: //13.3
				pcase += 1917 + 10*((f[0] < 0? f[2] > 0: 12 + (f[2] < 0)) + (f[1] < 0? f[3] < 0: 6 + (f[3] > 0)));
				if (f[4] > 0)
					pcase += 30;
				break;
			case 4: //13.2
				k = 21 + 11*f[0] + 4*f[1] + 3*f[2] + 2*f[3] + f[4];
				if (k >> 4)
					k -= (k&32? 20: 10);
				pcase += 1845 + 3*k;
				break;
			default: //13.1
				pcase += 1839 + 2*f[0];
			}
	}
	
	return ++pcase;
}

std::vector<std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>> Foam::cutCellFvMesh::MC33::collectTriangles
(
	const unsigned short int* triangleCase
)
{
	std::vector<std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>> triangles;
	do
	//while first hex in case is 1
	{
		std::uint8_t vertice1 = ((*triangleCase)&0x0F00)>>8;
		std::uint8_t vertice2 = ((*triangleCase)&0x00F0)>>4;
		std::uint8_t vertice3 = ((*triangleCase)&0x000F);
		triangles.push_back(std::tie(vertice1,vertice2,vertice3));
		if(triangles.size()>12)
			FatalErrorInFunction<<"Error. No cases with more than 12 triangles"<< exit(FatalError);  
	}
	while(*(triangleCase++)&0x1000);
	
	return triangles;
}

/******************************************************************
Vertices:           Faces:
    3 __________2        ___________
   /|          /|      /|          /|
  / |         / |     / |   2     / |
7/__________6/  |    /  |     4  /  |
|   |       |   |   |¯¯¯¯¯¯¯¯¯¯¯| 1 |     z
|   0_______|___1   | 3 |_______|___|     |
|  /        |  /    |  /  5     |  /      |____y
| /         | /     | /     0   | /      /
4/__________5/      |/__________|/      x


This function returns a vector with all six test face results (face[6]). Each
result value is 1 if the positive face vertices are joined, -1 if the negative
vertices are joined, and 0 (unchanged) if the test should not be applied. The
return value of this function is the the sum of all six results.*/
int Foam::cutCellFvMesh::MC33::face_tests(int *face, int i) const {
	if (i&0x80) { //vertex 0
		face[0] = ((i&0xCC) == 0x84? (v[0]*v[5] < v[1]*v[4]? -1: 1): 0);//0x84 = 10000100, vertices 0 and 5
		face[3] = ((i&0x99) == 0x81? (v[0]*v[7] < v[3]*v[4]? -1: 1): 0);//0x81 = 10000001, vertices 0 and 7
		face[4] = ((i&0xF0) == 0xA0? (v[0]*v[2] < v[1]*v[3]? -1: 1): 0);//0xA0 = 10100000, vertices 0 and 2
	} else {
		face[0] = ((i&0xCC) == 0x48? (v[0]*v[5] < v[1]*v[4]? 1: -1): 0);//0x48 = 01001000, vertices 1 and 4
		face[3] = ((i&0x99) == 0x18? (v[0]*v[7] < v[3]*v[4]? 1: -1): 0);//0x18 = 00011000, vertices 3 and 4
		face[4] = ((i&0xF0) == 0x50? (v[0]*v[2] < v[1]*v[3]? 1: -1): 0);//0x50 = 01010000, vertices 1 and 3
	}
	if (i&0x02) { //vertex 6
		face[1] = ((i&0x66) == 0x42? (v[1]*v[6] < v[2]*v[5]? -1: 1): 0);//0x42 = 01000010, vertices 1 and 6
		face[2] = ((i&0x33) == 0x12? (v[3]*v[6] < v[2]*v[7]? -1: 1): 0);//0x12 = 00010010, vertices 3 and 6
		face[5] = ((i&0x0F) == 0x0A? (v[4]*v[6] < v[5]*v[7]? -1: 1): 0);//0x0A = 00001010, vertices 4 and 6
	} else {
		face[1] = ((i&0x66) == 0x24? (v[1]*v[6] < v[2]*v[5]? 1: -1): 0);//0x24 = 00100100, vertices 2 and 5
		face[2] = ((i&0x33) == 0x21? (v[3]*v[6] < v[2]*v[7]? 1: -1): 0);//0x21 = 00100001, vertices 2 and 7
		face[5] = ((i&0x0F) == 0x05? (v[4]*v[6] < v[5]*v[7]? 1: -1): 0);//0x05 = 00000101, vertices 5 and 7
	}
	return face[0] + face[1] + face[2] + face[3] + face[4] + face[5];
}

/* Faster function for the face test, the test is applied to only one face
(int face). This function is only used for the cases 3 and 6 of MC33*/
int Foam::cutCellFvMesh::MC33::face_test1(int face) const {
	switch (face) {
	case 0:
		return (v[0]*v[5] < v[1]*v[4]? 0x48: 0x84);
	case 1:
		return (v[1]*v[6] < v[2]*v[5]? 0x24: 0x42);
	case 2:
		return (v[3]*v[6] < v[2]*v[7]? 0x21: 0x12);
	case 3:
		return (v[0]*v[7] < v[3]*v[4]? 0x18: 0x81);
	case 4:
		return (v[0]*v[2] < v[1]*v[3]? 0x50: 0xA0);
	default:
		return (v[4]*v[6] < v[5]*v[7]? 0x05: 0x0A);
	}
}

// Silence dereferencing type-punned pointer warning in GCC
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif
// an ugly signbit:
#if MC33_double_precision
inline unsigned int signbf(double x) {
	return ((*(reinterpret_cast<unsigned long long int*>(&x))>>32)&0x80000000);
}
#else
inline unsigned int signbf(float x) {
	return (*(reinterpret_cast<unsigned int*>(&x))&0x80000000);
}
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

/******************************************************************
Interior test function. If the test is positive, the function returns a value
different from 0. The integer i must be 0 to test if the vertices 0 and 6 are
joined. 1 for vertices 1 and 7, 2 for vertices 2 and 4, and 3 for 3 and 5.
For case 13, the integer flag13 must be 1, and the function returns 2 if one
of the vertices 0, 1, 2 or 3 is joined to the center point of the cube (case
13.5.2), returns 1 if one of the vertices 4, 5, 6 or 7 is joined to the
center point of the cube (case 13.5.2 too), and it returns 0 if the vertices
are not joined (case 13.5.1)*/
int Foam::cutCellFvMesh::MC33::interior_test(int i, int flag13) const {
	//Signs of cube vertices were changed to use signbit function in calc_isosurface
	//A0 = -v[0], B0 = -v[1], C0 = -v[2], D0 = -v[3]
	//A1 = -v[4], B1 = -v[5], C1 = -v[6], D1 = -v[7]
	//But the function still works
	double At = v[4] - v[0], Bt = v[5] - v[1],
				Ct = v[6] - v[2], Dt = v[7] - v[3];
	double t = At*Ct - Bt*Dt;//the "a" value.
	if (signbf(t)) {
		if (i&0x01) return 0;
	} else {
		if (!(i&0x01) || t == 0) return 0;
	}
	t = 0.5f*(v[3]*Bt + v[1]*Dt - v[2]*At - v[0]*Ct)/t;//t = -b/2a

	if (t > 0 && t < 1) {
		At = v[0] + At*t;
		Bt = v[1] + Bt*t;
		Ct = v[2] + Ct*t;
		Dt = v[3] + Dt*t;
		Ct *= At;
		Dt *= Bt;
		if (i&0x01) {
			if (Ct < Dt && signbf(Dt) == 0)
				return (signbf(Bt) == signbf(v[i])) + flag13;
		} else {
			if (Ct > Dt && signbf(Ct) == 0)
				return (signbf(At) == signbf(v[i])) + flag13;
		}
	}
	return 0;
}

#undef FF
#ifdef _MSC_VER
#pragma warning( pop )
#endif

