#include "cutCellFvMesh.H"

inline float invSqrt(float f) {
	return 1.0/sqrt(f);
}

bool Foam::cutCellFvMesh::MC33::detectInFaceTriangles
(
	std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>& triangle
)
{
	bool allInOneFace = false;
	for(auto iter=faceToEdges.cbegin(); iter!=faceToEdges.cend(); iter++)
	{
		std::uint8_t edge0 = std::get<0>(triangle);
		std::uint8_t edge1 = std::get<1>(triangle);
		std::uint8_t edge2 = std::get<2>(triangle);
		
		const std::unordered_set<std::uint8_t>& faceEdges = iter->second;
		
		bool vertice0In = faceEdges.find(edge0)!=faceEdges.end();
		bool vertice1In = faceEdges.find(edge1)!=faceEdges.end();
		bool vertice2In = faceEdges.find(edge2)!=faceEdges.end();
		
		if(vertice0In && vertice1In && vertice2In)
		{
			if(allInOneFace)
			{
				Info<<"edge0:"<<edge0<<Foam::endl;
				Info<<"edge1:"<<edge1<<Foam::endl;
				Info<<"edge2:"<<edge2<<Foam::endl;
				FatalErrorInFunction<<"Error!"<< exit(FatalError);
			}
			else
				allInOneFace = true;
		}
	}
	return allInOneFace;
}

unsigned int Foam::cutCellFvMesh::MC33::permuteBitPattern
(
	const std::array<unsigned short int,8>& permutation,
	unsigned int bitPattern
)
{
	unsigned short int resPattern = 0;
	for(unsigned int vertI=0;vertI<permutation.size();vertI++)
	{
		bool bitSet = getBit(bitPattern,7-permutation[vertI]);
		if(bitSet)
			setBit(resPattern,7-vertI);
	}
	return resPattern;
}

std::array<std::int8_t,6> permuteFacePattern
(
	const std::array<unsigned short int,6>& permutation,
	const std::array<std::int8_t,6>& facePattern
)
{
	std::array<std::int8_t,6> permutedFacePattern;
	for(uint faceInd=0; faceInd<facePattern.size(); faceInd++)
	{
		permutedFacePattern[faceInd] = facePattern[permutation[faceInd]];
	}
	return permutedFacePattern;
}

void Foam::cutCellFvMesh::MC33::setBit
(
	unsigned short int& bitField,
	std::uint8_t bitIndex
)
{
    if(bitIndex>15)
        FatalErrorInFunction<<"Bit Index must be inside [0,15]"<< exit(FatalError);

    unsigned short int mask = 1;
    mask = mask << bitIndex;
    
    bitField = bitField | mask;
}

bool Foam::cutCellFvMesh::MC33::getBit
(
	const unsigned short int bitField,
	std::uint8_t bitIndex
)
{
    if(bitIndex>15)
        FatalErrorInFunction<<"Bit Index must be inside [0,15]"<< exit(FatalError);

    unsigned short int mask = 1;
    mask = mask << bitIndex;
    
    return bitField & mask;
}

std::string Foam::cutCellFvMesh::MC33::getBitPattern
(
	const unsigned short int x
)
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

void Foam::cutCellFvMesh::MC33::computeFacePermutations()
{
	for(unsigned int permutationInd=0; permutationInd<facePermutations.size(); permutationInd++)
	{
		const std::array<unsigned short int,8>& onePointPermutation = permutations[permutationInd];
		std::array<unsigned short int,6>& oneFacePermutation = facePermutations[permutationInd];
		
		for(unsigned int faceInd=0; faceInd<oneFacePermutation.size(); faceInd++)
		{
			std::array<uint,4> facePnts = faceToPnts[faceInd];
			for(uint& pnt : facePnts)
			{
				pnt = onePointPermutation[pnt];
			}
			std::sort(facePnts.begin(),facePnts.end());
			
			uint num = 0;
			for(uint decInd=0; decInd<facePnts.size(); decInd++)
			{
				num *= 10;
				num += facePnts[decInd];
			}
			auto iter = pntsToFace.find(num);
			if(iter==pntsToFace.end())
			{
				for(auto iter=pntsToEdge.begin(); iter!=pntsToEdge.end(); iter++)
					Info<<iter->first<<"/"<<iter->second<<Foam::endl;
				Info<<"---------"<<Foam::endl;
				for(auto iter=pntsToFace.begin(); iter!=pntsToFace.end(); iter++)
					Info<<iter->first<<"/"<<iter->second<<Foam::endl;
				FatalErrorInFunction<<"Error num:"<<num<< exit(FatalError);
			}
			oneFacePermutation[faceInd] = iter->second;
		}
		std::bitset<6> refFaces;
		for(short unsigned int permFace : oneFacePermutation)
			refFaces[permFace] = true;
		if(!refFaces.all())
		{
			Info<<"[";
			for(label perm : oneFacePermutation)
				Info<<perm<<",";
			Info<<"]"<<Foam::endl;
			std::cout<<refFaces<<std::endl;
			FatalErrorInFunction<<"Mismatch face permutations"<< exit(FatalError);
		}
	}
}

Foam::cutCellFvMesh::MC33::MC33
(
    cutCellFvMesh& mesh
):
mesh(mesh),
table(mesh.caseTable),
permutations(mesh.permutationTable),
edgeToPnts(mesh.edgeToPnts),
pntsToEdge(mesh.pntsToEdge),
faceToPnts(mesh.faceToPnts),
pntsToFace(mesh.pntsToFace),
faceToEdges(mesh.faceToEdges)
{
	computeFacePermutations();
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
		if(startingFaceVert.find(vert)==startingFaceVert.end())
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
	
	for(uint mc33faceInd=0; mc33faceInd<faceVertices.size(); mc33faceInd++)
	{
		label matchMc33FaceInd = -1;
		for(label locFaceInd=0; locFaceInd<thisCell.size(); locFaceInd++)
		{
			const face& thisFace = meshFaces[thisCell[locFaceInd]];
			bool allIn = true;
			for(label vertice : thisFace)
			{
				auto iter=faceVertices[mc33faceInd].find(vertice);
				if(iter==faceVertices[mc33faceInd].end())
					allIn = false;
			}
			if(allIn)
			{
				if(matchMc33FaceInd!=-1)
					FatalErrorInFunction<<"Duplicate face assignment!"<< exit(FatalError);
				matchMc33FaceInd = locFaceInd;
			}
		}
		oneCube.origFaces[mc33faceInd] = thisCell[matchMc33FaceInd];
	}
	std::unordered_set<label> globFaces;
	for(label oFace : oneCube.origFaces)
	{
		auto iter = globFaces.find(oFace);
		if(iter!=globFaces.end())
		{
			Info<<"startingFaceVert:";
			for(auto iter=startingFaceVert.begin(); iter!=startingFaceVert.end(); iter++)
				Info<<*iter<<" ";
			Info<<Foam::endl;
			Info<<"oppositeFaceVert:";
			for(auto iter=oppositeFaceVert.begin(); iter!=oppositeFaceVert.end(); iter++)
				Info<<*iter<<" ";
			Info<<Foam::endl;
			Info<<"vertices:"<<Foam::endl;
			for(auto vert : vertices)
			{
				Info<<vert<<" ";
			}
			Info<<Foam::endl;
			Info<<"faceVertices:"<<Foam::endl;
			for(auto set : faceVertices)
			{
				for(auto iter=set.begin(); iter!=set.end(); iter++)
					Info<<*iter<<" ";
				Info<<Foam::endl;
			}
			Info<<"oneCube.vertices:"<<oneCube.vertices<<Foam::endl;
			Info<<"oneCube.origFaces:"<<oneCube.origFaces<<Foam::endl;
			for(auto faceInd : oneCube.origFaces)
			{
				Info<<faceInd<<": ";
				for(label vert : meshFaces[faceInd])
					Info<<vert<<" ";
				Info<<Foam::endl;
			}
			Info<<"thisCell:"<<thisCell<<Foam::endl;
			for(auto iter=globFaces.begin(); iter!=globFaces.end(); iter++)
				Info<<*iter<<" ";
			Info<<Foam::endl;
			Info<<" :"<<*iter<<Foam::endl;
			FatalErrorInFunction<<"Duplicate face assignment!"<< exit(FatalError);
		}
		globFaces.insert(oFace);
	}
	//FatalErrorInFunction<<"Temp Stop!"<< exit(FatalError);
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

std::array<std::int8_t,6> Foam::cutCellFvMesh::MC33::computeFacePattern()
{
	auto crossedFace = [](std::array<MC33_real,4> facePnts)
	{
		std::bitset<4> signs;
		//label count=0;
		for(uint i=0; i<signs.size(); i++)
		{
			signs[i] = (facePnts[i]>=0)?true:false;
		}
		if(signs[0]==signs[2] && signs[1]==signs[3])
		{
			if(signs[0]!=signs[1] && signs[2]!=signs[3])
				return true;
			else
				return false;
		}
		else
			return false;
	};
	std::bitset<8> signs;
	for(uint i=0; i<signs.size(); i++)
	{
		signs[i] = (v[i]>=0)?true:false;
	}
	
	std::array<std::int8_t,6> pattern;

	std::array<std::array<MC33_real,4>,6> faceV;
	faceV[0] = {v[4],v[5],v[1],v[0]};
	faceV[1] = {v[5],v[6],v[2],v[1]};
	faceV[2] = {v[3],v[2],v[6],v[7]};
	faceV[3] = {v[0],v[3],v[7],v[4]};
	faceV[4] = {v[0],v[1],v[2],v[3]};
	faceV[5] = {v[5],v[4],v[7],v[6]};

	if (signs[0]){ //i&0x80) { //vertex 0
		pattern[0] = (crossedFace(faceV[0])? (v[0]*v[5] < v[1]*v[4]? -1: 1): 0);//vertices 0 and 5
		pattern[3] = (crossedFace(faceV[3])? (v[0]*v[7] < v[3]*v[4]? -1: 1): 0);//vertices 0 and 7
		pattern[4] = (crossedFace(faceV[4])? (v[0]*v[2] < v[1]*v[3]? -1: 1): 0);//vertices 0 and 2
	} else {
		pattern[0] = (crossedFace(faceV[0])? (v[0]*v[5] < v[1]*v[4]? 1: -1): 0);//vertices 1 and 4
		pattern[3] = (crossedFace(faceV[3])? (v[0]*v[7] < v[3]*v[4]? 1: -1): 0);//vertices 3 and 4
		pattern[4] = (crossedFace(faceV[4])? (v[0]*v[2] < v[1]*v[3]? 1: -1): 0);//vertices 1 and 3
	}
	if (signs[6]){//i&0x02) { //vertex 6
		pattern[1] = (crossedFace(faceV[1])? (v[1]*v[6] < v[2]*v[5]? -1: 1): 0);//vertices 1 and 6
		pattern[2] = (crossedFace(faceV[2])? (v[3]*v[6] < v[2]*v[7]? -1: 1): 0);//vertices 3 and 6
		pattern[5] = (crossedFace(faceV[5])? (v[4]*v[6] < v[5]*v[7]? -1: 1): 0);//vertices 4 and 6
	} else {
		pattern[1] = (crossedFace(faceV[1])? (v[1]*v[6] < v[2]*v[5]? 1: -1): 0);//vertices 2 and 5
		pattern[2] = (crossedFace(faceV[2])? (v[3]*v[6] < v[2]*v[7]? 1: -1): 0);//vertices 2 and 7
		pattern[5] = (crossedFace(faceV[5])? (v[4]*v[6] < v[5]*v[7]? 1: -1): 0);//vertices 5 and 7
	}
	return pattern;
}

bool equalTriangles
(
	const std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>& mc33TriA,
	const std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>& mc33TriB
)
{
	std::unordered_set<std::uint8_t> mc33TriASet = {std::get<0>(mc33TriA),std::get<1>(mc33TriA),std::get<2>(mc33TriA)};
	bool equal = mc33TriASet.find(std::get<0>(mc33TriB))!=mc33TriASet.end() &&
				 mc33TriASet.find(std::get<1>(mc33TriB))!=mc33TriASet.end() &&
				 mc33TriASet.find(std::get<2>(mc33TriB))!=mc33TriASet.end();
	return equal;
};

using Triangle = std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>;
struct TriangleHasher
{
	std::size_t operator()(const Triangle& tri) const
	{
		return 	std::hash<std::uint8_t>()(std::get<0>(tri)) ^ 
				std::hash<std::uint8_t>()(std::get<1>(tri)) ^
				std::hash<std::uint8_t>()(std::get<2>(tri));
	}
};
struct TriangleEqual
{
	std::size_t operator()(const Triangle& triA, const Triangle& triB) const
	{
		return 	equalTriangles(triA,triB);
	}
};

void Foam::cutCellFvMesh::getCubeCaseMC33Triangles
(
	MC33::Case cubeCase,
    bool redMarkPlus,
    DynamicList<std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>>& mc33Triangles
)
{
	const markSide& thisCaseCorrData = convexCorrectionData[cubeCase];

	std::unordered_set<Triangle,TriangleHasher,TriangleEqual> mc33TrianglesSet;
	for(const corrFace& oneCorrFace : thisCaseCorrData.redMarkPlus.faces)
	{
		//Info<<"Iterate red :"<<oneCorrFace.to_string();
		if(oneCorrFace.type==FType::mc33Triangle)
		{
			Triangle oneMc33Triangle;
			if(oneCorrFace.faceData.size()!=3)
				FatalErrorInFunction<<"MC33 must be sized 3!"<< exit(FatalError);
			std::get<0>(oneMc33Triangle) = static_cast<std::uint8_t>(oneCorrFace.faceData[0].value);
			std::get<1>(oneMc33Triangle) = static_cast<std::uint8_t>(oneCorrFace.faceData[1].value);
			std::get<2>(oneMc33Triangle) = static_cast<std::uint8_t>(oneCorrFace.faceData[2].value);
			if(mc33TrianglesSet.find(oneMc33Triangle)==mc33TrianglesSet.end())
			{
				mc33Triangles.append(oneMc33Triangle);
				mc33TrianglesSet.insert(oneMc33Triangle);
				//Info<<" insert";
			}
		}
		//Info<<Foam::endl;
	}
	for(const corrFace& oneCorrFace : thisCaseCorrData.redMarkNotPlus.faces)
	{
		//Info<<"Iterate non red :"<<oneCorrFace.to_string();
		if(oneCorrFace.type==FType::mc33Triangle)
		{
			Triangle oneMc33Triangle;
			if(oneCorrFace.faceData.size()!=3)
				FatalErrorInFunction<<"MC33 must be sized 3!"<< exit(FatalError);
			std::get<0>(oneMc33Triangle) = static_cast<std::uint8_t>(oneCorrFace.faceData[0].value);
			std::get<1>(oneMc33Triangle) = static_cast<std::uint8_t>(oneCorrFace.faceData[1].value);
			std::get<2>(oneMc33Triangle) = static_cast<std::uint8_t>(oneCorrFace.faceData[2].value);
			if(mc33TrianglesSet.find(oneMc33Triangle)==mc33TrianglesSet.end())
			{
				mc33Triangles.append(oneMc33Triangle);
				mc33TrianglesSet.insert(oneMc33Triangle);
				//Info<<" insert";
			}
		}
		//Info<<Foam::endl;
	}
}

Foam::cutCellFvMesh::MC33::MC33Cube Foam::cutCellFvMesh::MC33::computeCutCell
(
	int cellInd
)
{
	MC33Cube mc33Cube = generateMC33Cube(cellInd);
	unsigned int bitPattern = computeSignBitPattern(mc33Cube);
	auto compEdgePermFromPointPerm = [&](const std::array<unsigned short int,8>& pointPermutation)
	{
		std::array<std::uint8_t,12> edgePermutation;
		for(unsigned int edgeInd=0; edgeInd<edgePermutation.size(); edgeInd++)
		{
			std::array<uint,2> edgePnts = edgeToPnts[edgeInd];
			for(uint& pnt : edgePnts)
			{
				pnt = pointPermutation[pnt];
			}
			std::sort(edgePnts.begin(),edgePnts.end());
			
			uint num = 0;
			for(uint decInd=0; decInd<edgePnts.size(); decInd++)
			{
				num *= 10;
				num += edgePnts[decInd];
			}
			auto iter = pntsToEdge.find(num);
			if(iter==pntsToEdge.end())
			{				
				FatalErrorInFunction<<"Error num:"<<num<< exit(FatalError);
			}
			edgePermutation[edgeInd] = iter->second;
		}
		std::bitset<12> refEdges;
		for(std::uint8_t permEdg : edgePermutation)
			refEdges[permEdg] = true;
		if(!refEdges.all())
			FatalErrorInFunction<<"Mismatch edge permutations"<< exit(FatalError);
		return edgePermutation;
	};
	if(bitPattern!=0 && bitPattern!=255) 
	{
		mc33Cube.bitPattern = bitPattern;
		mc33Cube.facePattern = computeFacePattern();
		auto dominantSign = [](std::array<std::int8_t,6> facePattern)
		{
			std::uint8_t posCount=0;
			std::uint8_t negCount=0;
			for(std::int8_t ind : facePattern)
			{
				if(ind>0)
					posCount++;
				if(ind<0)
					negCount++;
			}
			if(negCount>posCount)
				return -1;
			else if(posCount>negCount)
				return +1;
			else
				return 0;
		};
		const unsigned short int* triangleCase = getTriangleCase(bitPattern);
		label arrayIndex = triangleCase-table;
		mc33Cube.cutTriangles = collectTriangles(triangleCase);

		unsigned int referenceBitPattern=0;
		std::array<std::int8_t,6> referenceFacePattern;
		std::vector<std::uint8_t> validPermutations;
		{
		// Get starting array index in MC33 list
			if(arrayIndex<0)
				FatalErrorInFunction<<"ArrayIndex must not be smaller than zero!"<< exit(FatalError);
			else if(arrayIndex<128)
				FatalErrorInFunction<<"ArrayIndex not a case!"<< exit(FatalError);
			else if(arrayIndex<136)
				// Case 1 Triangle Number 1
				mc33Cube.cubeCase = MC33::Case::c1;
			else if(arrayIndex<160)
				// Case 2 Triangle Number 2
				mc33Cube.cubeCase = MC33::Case::c2;
			else if(arrayIndex<184)
				// Case 3.1 Triangle Number 2
				mc33Cube.cubeCase = MC33::Case::c31;
			else if(arrayIndex<232)
				// Case 3.2 Triangle Number 4
				mc33Cube.cubeCase = MC33::Case::c32;
			else if(arrayIndex<240)
				// Case 4.1.1 Triangle Number 2
				mc33Cube.cubeCase = MC33::Case::c411;
			else if(arrayIndex<264)
				// Case 4.1.2 Triangle Number 6
				mc33Cube.cubeCase = MC33::Case::c412;
			else if(arrayIndex<336)
				// Case 5 Triangle Number 3
				mc33Cube.cubeCase = MC33::Case::c5;
			else if(arrayIndex<408)
				// Case 6.1.1 Triangle Number 3
				mc33Cube.cubeCase = MC33::Case::c611;
			else if(arrayIndex<576)
				// Case 6.1.2 Triangle Number 7
				mc33Cube.cubeCase = MC33::Case::c612;
			else if(arrayIndex<696)
				// Case 6.2 Triangle Number 5
				mc33Cube.cubeCase = MC33::Case::c62;
			else if(arrayIndex<720)
				// Case 7.1 Triangle Number 3
				mc33Cube.cubeCase = MC33::Case::c71;
			else if(arrayIndex<760)
				// Case 7.2 Triangle Number 5
				mc33Cube.cubeCase = MC33::Case::c72_720;
			else if(arrayIndex<800)
				// Case 7.2 Triangle Number 5
				mc33Cube.cubeCase = MC33::Case::c72_760;
			else if(arrayIndex<840)
				// Case 7.2 Triangle Number 5
				mc33Cube.cubeCase = MC33::Case::c72_800;
			else if(arrayIndex<912)
				// Case 7.3 Triangle Number 9
				mc33Cube.cubeCase = MC33::Case::c73_840;
			else if(arrayIndex<984)
				// Case 7.3 Triangle Number 9
				mc33Cube.cubeCase = MC33::Case::c73_912;
			else if(arrayIndex<1056)
				// Case 7.3 Triangle Number 9
				mc33Cube.cubeCase = MC33::Case::c73_984;
			else if(arrayIndex<1096)
				// Case 7.4.1 Triangle Number 5
				mc33Cube.cubeCase = MC33::Case::c741;
			else if(arrayIndex<1168)
				// Case 7.4.2 Triangle Number 9
				mc33Cube.cubeCase = MC33::Case::c742;
			else if(arrayIndex<1174)
				// Case 8 Triangle Number 2
				mc33Cube.cubeCase = MC33::Case::c8;
			else if(arrayIndex<1190)
				// Case 9 Triangle Number 4
				mc33Cube.cubeCase = MC33::Case::c9;
			else if(arrayIndex<1202)
				// Case 10.1.1 Triangle Number 4
				mc33Cube.cubeCase = MC33::Case::c1011_1190;
			else if(arrayIndex<1214)
				// Case 10.1.1 Triangle Number 4
				mc33Cube.cubeCase = MC33::Case::c1011_1202;
			else if(arrayIndex<1238)
				// Case 10.1.2 Triangle Number 8
				mc33Cube.cubeCase = MC33::Case::c1012_1214;
			else if(arrayIndex<1262)
				// Case 10.1.2 Triangle Number 8
				mc33Cube.cubeCase = MC33::Case::c1012_1238;
			else if(arrayIndex<1286)
				// Case 10.2 Triangle Number 8
				mc33Cube.cubeCase = MC33::Case::c102_1262;
			else if(arrayIndex<1310)
				// Case 10.2 Triangle Number 8
				mc33Cube.cubeCase = MC33::Case::c102_1286;
			else if(arrayIndex<1334)
				// Case 11 Triangle Number 4
				mc33Cube.cubeCase = MC33::Case::c11;
			else if(arrayIndex<1358)
				// Case 14 Triangle Number 4
				mc33Cube.cubeCase = MC33::Case::c14;
			else if(arrayIndex<1406)
				// Case 12.1.1 Triangle Number 4
				mc33Cube.cubeCase = MC33::Case::c1211_1358;
			else if(arrayIndex<1454)
				// Case 12.1.1 Triangle Number 4
				mc33Cube.cubeCase = MC33::Case::c1211_1406;
			else if(arrayIndex<1550)
				// Case 12.1.2 Triangle Number 8
				mc33Cube.cubeCase = MC33::Case::c1212_1454;
			else if(arrayIndex<1646)
				// Case 12.1.2 Triangle Number 8
				mc33Cube.cubeCase = MC33::Case::c1212_1550;
			else if(arrayIndex<1742)
				// Case 12.2 Triangle Number 8
				mc33Cube.cubeCase = MC33::Case::c122_1646;
			else if(arrayIndex<1838)
				// Case 12.2 Triangle Number 8
				mc33Cube.cubeCase = MC33::Case::c122_1742;
			else if(arrayIndex<1846)
				// Case 13.1 Triangle Number 4
				mc33Cube.cubeCase = MC33::Case::c131;
			else if(arrayIndex<1918)
			{
				// Case 13.2 Triangle Number 6
				//mc33Cube.cubeCase = c132;
				std::int8_t domSign = dominantSign(mc33Cube.facePattern);
				if(domSign==+1)
					mc33Cube.cubeCase = MC33::Case::c132_MostPos;
				else if(domSign==-1)
					mc33Cube.cubeCase = MC33::Case::c132_MostNeg;
				else
					FatalErrorInFunction<<"False Sign!"<< exit(FatalError);
			}
			else if(arrayIndex<2158)
			{
				// Case 13.3 Triangle Number 10
				//mc33Cube.cubeCase = c133;
				std::int8_t domSign = dominantSign(mc33Cube.facePattern);
				if(domSign==+1)
					mc33Cube.cubeCase = MC33::Case::c133_MostPos;
				else if(domSign==-1)
					mc33Cube.cubeCase = MC33::Case::c133_MostNeg;
				else
					FatalErrorInFunction<<"False Sign!"<< exit(FatalError);
			}
			else if(arrayIndex<2206)
				// Case 13.4 Triangle Number 12
				mc33Cube.cubeCase = MC33::Case::c134;
			else if(arrayIndex<2246)
				// Case 13.5.2 Triangle Number 10
				mc33Cube.cubeCase = MC33::Case::c1352_2206;
			else if(arrayIndex<2286)
				// Case 13.5.2 Triangle Number 10
				mc33Cube.cubeCase = MC33::Case::c1352_2246;
			else if(arrayIndex<2310)
				// Case 13.5.1 Triangle Number 6
				mc33Cube.cubeCase = MC33::Case::c1351;
			else
				FatalErrorInFunction<<"ArrayIndex must not be larger than 2309!"<< exit(FatalError);

		// Getting referenceBitPattern for every MC33 case
			if(mc33Cube.cubeCase==MC33::Case::c1){
				referenceBitPattern = 0x0008;
				mc33Cube.permMeth = permutationMethod::Point;
			}
			else if(mc33Cube.cubeCase==MC33::Case::c2){
				referenceBitPattern = 0x000C;
				mc33Cube.permMeth = permutationMethod::Point;
			}
			else if(mc33Cube.cubeCase==MC33::Case::c31 || mc33Cube.cubeCase==MC33::Case::c32){
				referenceBitPattern = 0x000A;
				mc33Cube.permMeth = permutationMethod::Point;
			}
			else if(mc33Cube.cubeCase==MC33::Case::c411 || mc33Cube.cubeCase==MC33::Case::c412){
				referenceBitPattern = 0x0028;
				mc33Cube.permMeth = permutationMethod::Point;
			}
			else if(mc33Cube.cubeCase==MC33::Case::c5){
				referenceBitPattern = 0x00C4;
				mc33Cube.permMeth = permutationMethod::Point;
			}
			else if(mc33Cube.cubeCase==MC33::Case::c611 || mc33Cube.cubeCase==MC33::Case::c612 || mc33Cube.cubeCase==MC33::Case::c62){
				referenceBitPattern = 0x002C;
				mc33Cube.permMeth = permutationMethod::Point;
			}
			else if(mc33Cube.cubeCase==MC33::Case::c71 || mc33Cube.cubeCase==MC33::Case::c72_720 || mc33Cube.cubeCase==MC33::Case::c72_760 || mc33Cube.cubeCase==MC33::Case::c72_800){
				referenceBitPattern = 0x0025;
				mc33Cube.permMeth = permutationMethod::Both;
				if(mc33Cube.cubeCase==MC33::Case::c72_720)
					referenceFacePattern = {0,-1,-1,0,0,+1};
				else if(mc33Cube.cubeCase==MC33::Case::c72_760)
					referenceFacePattern = {0,+1,-1,0,0,-1};
				else if(mc33Cube.cubeCase==MC33::Case::c72_800)
					referenceFacePattern = {0,-1,+1,0,0,-1};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c73_840 || mc33Cube.cubeCase==MC33::Case::c73_912 || mc33Cube.cubeCase==MC33::Case::c73_984){
				referenceBitPattern = 0x0025;
				mc33Cube.permMeth = permutationMethod::Both;
				if(mc33Cube.cubeCase==MC33::Case::c73_840)
					referenceFacePattern = {0,+1,-1,0,0,+1};
				else if(mc33Cube.cubeCase==MC33::Case::c73_912)
					referenceFacePattern = {0,-1,+1,0,0,+1};
				else if(mc33Cube.cubeCase==MC33::Case::c73_984)
					referenceFacePattern = {0,+1,+1,0,0,-1};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c741 || mc33Cube.cubeCase==MC33::Case::c742){
				referenceBitPattern = 0x0025;
				mc33Cube.permMeth = permutationMethod::Point;
			}
			else if(mc33Cube.cubeCase==MC33::Case::c8){
				referenceBitPattern = 0x00CC;
				mc33Cube.permMeth = permutationMethod::Point;
			}
			else if(mc33Cube.cubeCase==MC33::Case::c9){
				referenceBitPattern = 0x00D8;
				mc33Cube.permMeth = permutationMethod::Point;
			}
			else if(mc33Cube.cubeCase==MC33::Case::c1011_1190 || mc33Cube.cubeCase==MC33::Case::c1011_1202){
				referenceBitPattern = 0x0096;
				mc33Cube.permMeth = permutationMethod::Both;
				if(mc33Cube.cubeCase==MC33::Case::c1011_1190)
					referenceFacePattern = {-1,0,-1,0,0,0};
				else if(mc33Cube.cubeCase==MC33::Case::c1011_1202)
					referenceFacePattern = {+1,0,+1,0,0,0};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c1012_1214 || mc33Cube.cubeCase==MC33::Case::c1012_1238){
				referenceBitPattern = 0x0096;
				mc33Cube.permMeth = permutationMethod::Both;
				if(mc33Cube.cubeCase==MC33::Case::c1012_1214)
					referenceFacePattern = {-1,0,-1,0,0,0};
				else if(mc33Cube.cubeCase==MC33::Case::c1012_1238)
					referenceFacePattern = {+1,0,+1,0,0,0};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c102_1262 || mc33Cube.cubeCase==MC33::Case::c102_1286){
				referenceBitPattern = 0x0096;
				mc33Cube.permMeth = permutationMethod::Both;
				if(mc33Cube.cubeCase==MC33::Case::c102_1262)
					referenceFacePattern = {+1,0,-1,0,0,0};
				else if(mc33Cube.cubeCase==MC33::Case::c102_1286)
					referenceFacePattern = {-1,0,+1,0,0,0};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c11){
				referenceBitPattern = 0x00E8;
				mc33Cube.permMeth = permutationMethod::Point;
			}
			else if(mc33Cube.cubeCase==MC33::Case::c1211_1358 || mc33Cube.cubeCase==MC33::Case::c1211_1406){
				referenceBitPattern = 0x00A5;
				mc33Cube.permMeth = permutationMethod::Both;
				if(mc33Cube.cubeCase==MC33::Case::c1211_1358)
					referenceFacePattern = {0,0,0,-1,0,-1};
				else if(mc33Cube.cubeCase==MC33::Case::c1211_1406)
					referenceFacePattern = {0,0,0,+1,0,+1};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c1212_1454 || mc33Cube.cubeCase==MC33::Case::c1212_1550){
				referenceBitPattern = 0x00A5;
				mc33Cube.permMeth = permutationMethod::Both;
				if(mc33Cube.cubeCase==MC33::Case::c1212_1454)
					referenceFacePattern = {0,0,0,-1,0,-1};
				else if(mc33Cube.cubeCase==MC33::Case::c1212_1550)
					referenceFacePattern = {0,0,0,+1,0,+1};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c122_1646 || mc33Cube.cubeCase==MC33::Case::c122_1742){
				referenceBitPattern = 0x00A5;
				mc33Cube.permMeth = permutationMethod::Both;
				if(mc33Cube.cubeCase==MC33::Case::c1212_1454)
					referenceFacePattern = {0,0,0,+1,0,-1};
				else if(mc33Cube.cubeCase==MC33::Case::c1212_1550)
					referenceFacePattern = {0,0,0,-1,0,+1};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c131){
				referenceBitPattern = 0x00A5;
				mc33Cube.permMeth = permutationMethod::Point; // False but irrelavant
			}
			else if(mc33Cube.cubeCase==MC33::Case::c132_MostNeg){
				referenceBitPattern = 0x00A5;
				mc33Cube.permMeth = permutationMethod::Face;
				referenceFacePattern = {-1,-1,+1,-1,-1,-1};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c132_MostPos){
				referenceBitPattern = 0x00A5;
				mc33Cube.permMeth = permutationMethod::Face;
				referenceFacePattern = {+1,+1,-1,+1,+1,+1};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c133_MostPos){
				referenceBitPattern = 0x00A5;
				mc33Cube.permMeth = permutationMethod::Face;
				referenceFacePattern = {+1,+1,-1,-1,+1,+1};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c133_MostNeg){
				referenceBitPattern = 0x00A5;
				mc33Cube.permMeth = permutationMethod::Face;
				referenceFacePattern = {-1,-1,+1,+1,-1,-1};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c134){
				referenceBitPattern = 0x00A5;
				mc33Cube.permMeth = permutationMethod::Face;
				referenceFacePattern = {-1,-1,+1,+1,-1,+1};
			}
			else if(mc33Cube.cubeCase==MC33::Case::c1351){
				referenceBitPattern = 0x00A5;
				mc33Cube.permMeth = permutationMethod::Point; // False but irrelevant
			}
			else if(mc33Cube.cubeCase==MC33::Case::c1352_2206  || mc33Cube.cubeCase==MC33::Case::c1352_2246){
				referenceBitPattern = 0x00A5;
				mc33Cube.permMeth = permutationMethod::Point; // False but irrelevant
			}
			else if(mc33Cube.cubeCase==MC33::Case::c14){
				referenceBitPattern = 0x00D4;
				mc33Cube.permMeth = permutationMethod::Point;
			}
			else
				FatalErrorInFunction<<"Error!"<< exit(FatalError);
			
		//Compute permutations for MC33 case and assign the correct one
			auto findPointPermutation = [&]()
			{
				std::vector<std::uint8_t> validPermutationsRef;
				std::vector<std::uint8_t> validPermutationsInvert;
				for(unsigned int permutInd=0;permutInd<permutations.size();permutInd++)
				{
					if(referenceBitPattern==permuteBitPattern(permutations[permutInd],bitPattern))
					{
						validPermutationsRef.push_back(permutInd);
					}
				}
				referenceBitPattern = referenceBitPattern^0x00FF;
				for(unsigned int permutInd=0;permutInd<permutations.size();permutInd++)
				{
					if(referenceBitPattern==permuteBitPattern(permutations[permutInd],bitPattern))
					{
						validPermutationsInvert.push_back(permutInd);
					}
				}
				if(validPermutationsRef.size()>0 && validPermutationsInvert.size()==0)
				{
					mc33Cube.redMarkIsPlusSide = true;
					validPermutations = validPermutationsRef;
				}
				else if(validPermutationsRef.size()==0 && validPermutationsInvert.size()>0)
				{
					mc33Cube.redMarkIsPlusSide = false;
					validPermutations = validPermutationsInvert;
				}
				else if(validPermutationsRef.size()>0 && validPermutationsInvert.size()>0)
				{
					if(mc33Cube.cubeCase==c8 ||
					mc33Cube.cubeCase==c9 ||
					mc33Cube.cubeCase==c1011_1190 || mc33Cube.cubeCase==c1011_1202 || 
					mc33Cube.cubeCase==c1012_1214 || mc33Cube.cubeCase==c1012_1238 ||
					mc33Cube.cubeCase==c102_1262 || mc33Cube.cubeCase==c102_1286 ||
					mc33Cube.cubeCase==c11 ||
					mc33Cube.cubeCase==c1211_1358 || mc33Cube.cubeCase==c1211_1406 ||
					mc33Cube.cubeCase==c1212_1454 || mc33Cube.cubeCase==c1212_1550|| 
					mc33Cube.cubeCase==c122_1646 || mc33Cube.cubeCase==c122_1742 ||
					mc33Cube.cubeCase==c131 || mc33Cube.cubeCase==c132_MostNeg || 
					mc33Cube.cubeCase==c132_MostPos || mc33Cube.cubeCase==c133_MostNeg ||
					mc33Cube.cubeCase==c133_MostPos ||
					mc33Cube.cubeCase==c134 || mc33Cube.cubeCase==c1351 ||
					mc33Cube.cubeCase==c1352_2206 || mc33Cube.cubeCase==c1352_2246 ||
					mc33Cube.cubeCase==c14)
					{
						mc33Cube.redMarkIsPlusSide = true;
						validPermutations = validPermutationsRef;
					}
				}
				else
				{
					Info<<"mc33Cube.cubeCase:"<<mc33Cube.cubeCase<<Foam::endl;
					Info<<"referenceBitPattern:"<<referenceBitPattern<<Foam::endl;
					Info<<"validPermutationsRef.size():"<<validPermutationsRef.size()<<Foam::endl;
					Info<<"validPermutationsInvert.size():"<<validPermutationsInvert.size()<<Foam::endl;
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				}
				
				using Triangle = std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>;

				DynamicList<Triangle> mc33TrianglesInCorrection;
				mesh.getCubeCaseMC33Triangles
				(
					mc33Cube.cubeCase,
					mc33Cube.redMarkIsPlusSide,
					mc33TrianglesInCorrection
				);
				const std::vector<Triangle>& cutTriangles = mc33Cube.cutTriangles;
				std::unordered_set<Triangle,TriangleHasher,TriangleEqual> cutTrianglesSet(cutTriangles.cbegin(),cutTriangles.cend());
				
				/*
				Info<<"cutTriangles:"<<Foam::endl;
				for(Triangle tri : cutTriangles)
				{
					Info<<"   ("<<std::get<0>(tri)<<","<<std::get<1>(tri)<<","<<std::get<2>(tri)<<")"<<Foam::endl;
				}
				Info<<"mc33TrianglesInCorrection:"<<Foam::endl;
				for(Triangle tri : mc33TrianglesInCorrection)
				{
					Info<<"   ("<<std::get<0>(tri)<<","<<std::get<1>(tri)<<","<<std::get<2>(tri)<<")"<<Foam::endl;
				}
				*/
				
				std::vector<std::uint8_t> validPermutationsMatchTri;
				for(label onePermutation : validPermutations)
				{
					//Info<<"onePermutation:"<<onePermutation<<" (";
					std::array<unsigned short int,8> pointPerm = permutations[onePermutation];
					std::array<std::uint8_t,12> edgePerm = compEdgePermFromPointPerm(pointPerm);
					
					/*
					for(std::uint8_t pedg : edgePerm)
					{
						Info<<" "<<pedg<<" ";
					}
					Info<<")"<<Foam::endl;
					*/
					
					List<Triangle> mc33TrianglesInCorrPerm(mc33TrianglesInCorrection.size());
					for(label triInd=0; triInd<mc33TrianglesInCorrection.size(); triInd++)
					{
						Triangle mc33Tri = mc33TrianglesInCorrection[triInd];
						std::get<0>(mc33Tri) = edgePerm[std::get<0>(mc33Tri)];
						std::get<1>(mc33Tri) = edgePerm[std::get<1>(mc33Tri)];
						std::get<2>(mc33Tri) = edgePerm[std::get<2>(mc33Tri)];
						mc33TrianglesInCorrPerm[triInd] = mc33Tri;
					}
					
					bool allMatch = true;
					for(const Triangle& corrMc33 : mc33TrianglesInCorrPerm)
					{
						if(cutTrianglesSet.find(corrMc33)==cutTrianglesSet.end())
							allMatch = false;;
					}
					if(allMatch)
					{
						validPermutationsMatchTri.push_back(onePermutation);
					}
				}
				if(validPermutationsMatchTri.size()<1)
					FatalErrorInFunction<<"No matching permutation!"<< exit(FatalError);
				validPermutations = validPermutationsMatchTri;
				return validPermutations[0];
			};
			auto findFacePermutation = [&]()
			{
				std::vector<std::uint8_t> validPermutations;
				for(unsigned int permutInd=0;permutInd<permutations.size();permutInd++)
				{
					if(referenceFacePattern==permuteFacePattern(facePermutations[permutInd],mc33Cube.facePattern))
					{
						validPermutations.push_back(permutInd);
					}
				}

				if(validPermutations.size()==0)
				{
					Info<<"mc33Cube.cubeCase:"<<mc33Cube.cubeCase<<Foam::endl;
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				}
				return validPermutations[0];
			};
			
			std::uint8_t permutationTableIndex = 0;
			if(mc33Cube.permMeth==permutationMethod::Point)
			{
				permutationTableIndex = findPointPermutation();
			}
			else if(mc33Cube.permMeth==permutationMethod::Both)
			{
				permutationTableIndex = findPointPermutation();
			}
			else if(mc33Cube.permMeth == permutationMethod::Face)
			{
				permutationTableIndex = findFacePermutation();
				unsigned int referenceBitPatternInv = referenceBitPattern^0x00FF;
				if(referenceBitPattern==permuteBitPattern(permutations[permutationTableIndex],bitPattern))
				{
					mc33Cube.redMarkIsPlusSide = true;
				}
				else if(referenceBitPatternInv==permuteBitPattern(permutations[permutationTableIndex],bitPattern))
				{
					mc33Cube.redMarkIsPlusSide = false;
				}
				else
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
			}
			mc33Cube.permutationTableIndex = permutationTableIndex;
			
			/*
			if(mc33Cube.cubeCase==c1){
				// Case 1 Triangle Number 1
				if(mc33Cube.cutTriangles.size()!=1)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=3)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c2 ){
				// Case 2 Triangle Number 2
				if(mc33Cube.cutTriangles.size()!=2)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=2)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c31 ){
				// Case 3.1 Triangle Number 2
				if(mc33Cube.cutTriangles.size()!=2)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=2)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c32 ){
				// Case 3.2 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=2)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c411 ){
				// Case 4.1.1 Triangle Number 2
				if(mc33Cube.cutTriangles.size()!=2)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=3)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c412 ){
				// Case 4.1.2 Triangle Number 6
				if(mc33Cube.cutTriangles.size()!=6)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=3)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c5 ){
				// Case 5 Triangle Number 3
				if(mc33Cube.cutTriangles.size()!=3)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c611 ){
				// Case 6.1.1 Triangle Number 3
				if(mc33Cube.cutTriangles.size()!=3)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c612 ){
				// Case 6.1.2 Triangle Number 7
				if(mc33Cube.cutTriangles.size()!=7)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c62 ){
				// Case 6.2 Triangle Number 5
				if(mc33Cube.cutTriangles.size()!=5)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c71 ){
				// Case 7.1 Triangle Number 3
				if(mc33Cube.cutTriangles.size()!=3)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c72_720 ||  mc33Cube.cubeCase==c72_760 || mc33Cube.cubeCase==c72_800){
				// Case 7.2 Triangle Number 5
				if(mc33Cube.cutTriangles.size()!=5)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c73_840 || mc33Cube.cubeCase==c73_912 || mc33Cube.cubeCase==c73_984){
				// Case 7.3 Triangle Number 9
				if(mc33Cube.cutTriangles.size()!=9)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c741 ){
				// Case 7.4.1 Triangle Number 5
				if(mc33Cube.cutTriangles.size()!=5)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c742 ){
				// Case 7.4.2 Triangle Number 9
				if(mc33Cube.cutTriangles.size()!=9)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c8 ){
				// Case 8 Triangle Number 2
				if(mc33Cube.cutTriangles.size()!=2)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=4)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c9){
				// Case 9 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=3)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c1011_1190 || mc33Cube.cubeCase==c1011_1202){
				// Case 10.1.1 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=4)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c1012_1214 || mc33Cube.cubeCase==c1012_1238){
				// Case 10.1.2 Triangle Number 8
				if(mc33Cube.cutTriangles.size()!=8)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=4)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c102_1262 || mc33Cube.cubeCase==c102_1286){
				// Case 10.2 Triangle Number 8
				if(mc33Cube.cutTriangles.size()!=8)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=4)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c11 ){
				// Case 11 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c14 ){
				// Case 14 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c1211_1358 || mc33Cube.cubeCase==c1211_1406){
				// Case 12.1.1 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c1212_1454 || mc33Cube.cubeCase==c1212_1550){
				// Case 12.1.2 Triangle Number 8
				if(mc33Cube.cutTriangles.size()!=8)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c122_1646 || mc33Cube.cubeCase==c122_1742){
				// Case 12.2 Triangle Number 8
				if(mc33Cube.cutTriangles.size()!=8)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=1)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c131 ){
				// Case 13.1 Triangle Number 4
				if(mc33Cube.cutTriangles.size()!=4)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=4)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c132 ){
				// Case 13.2 Triangle Number 6
				if(mc33Cube.cutTriangles.size()!=6)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=4)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c133 ){
				// Case 13.3 Triangle Number 10
				if(mc33Cube.cutTriangles.size()!=10)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=4)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c134 ){
				// Case 13.4 Triangle Number 12
				if(mc33Cube.cutTriangles.size()!=12)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=4)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c1352_2206 || mc33Cube.cubeCase==c1352_2246){
				// Case 13.5.2 Triangle Number 10
				if(mc33Cube.cutTriangles.size()!=10)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=4)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else if(mc33Cube.cubeCase==c1351){
				// Case 13.5.1 Triangle Number 6
				if(mc33Cube.cutTriangles.size()!=6)
					FatalErrorInFunction<<"Error!"<< exit(FatalError);
				if(validPermutations.size()!=4)
					FatalErrorInFunction<<"Error! validPermutations.size():"<<validPermutations.size()<< exit(FatalError);
			}
			else
				FatalErrorInFunction<<"Error!"<< exit(FatalError);
			*/
		}

		/*
		if(mc33Cube.cubeCase == MC33::Case::c8)
		{
			Info<<"referenceBitPattern:"<<referenceBitPattern<<Foam::endl;
			Info<<"bitPattern:"<<bitPattern<<Foam::endl;
			Info<<"validPermutations: (";
			for(auto val : validPermutations)
				Info<<" "<<val<<" ";
			Info<<")"<<Foam::endl;
			FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);
		}
		*/
		
		mc33Cube.pointPermutation = permutations[mc33Cube.permutationTableIndex];
		std::bitset<8> refPoints;
		for(short unsigned int permPnt : mc33Cube.pointPermutation)
			refPoints[permPnt] = true;
		if(!refPoints.all())
			FatalErrorInFunction<<"Mismatch point permutations"<< exit(FatalError);

		for(unsigned int edgeInd=0; edgeInd<mc33Cube.edgePermutation.size(); edgeInd++)
		{
			std::array<uint,2> edgePnts = edgeToPnts[edgeInd];
			for(uint& pnt : edgePnts)
			{
				pnt = mc33Cube.pointPermutation[pnt];
			}
			std::sort(edgePnts.begin(),edgePnts.end());
			
			uint num = 0;
			for(uint decInd=0; decInd<edgePnts.size(); decInd++)
			{
				num *= 10;
				num += edgePnts[decInd];
			}
			auto iter = pntsToEdge.find(num);
			if(iter==pntsToEdge.end())
			{				
				FatalErrorInFunction<<"Error num:"<<num<< exit(FatalError);
			}
			mc33Cube.edgePermutation[edgeInd] = iter->second;
		}
		mc33Cube.edgePermutation = compEdgePermFromPointPerm(mc33Cube.pointPermutation);
		
		for(unsigned int faceInd=0; faceInd<mc33Cube.facePermutation.size(); faceInd++)
		{
			std::array<uint,4> facePnts = faceToPnts[faceInd];
			for(uint& pnt : facePnts)
			{
				pnt = mc33Cube.pointPermutation[pnt];
			}
			std::sort(facePnts.begin(),facePnts.end());
			
			uint num = 0;
			for(uint decInd=0; decInd<facePnts.size(); decInd++)
			{
				num *= 10;
				num += facePnts[decInd];
			}
			auto iter = pntsToFace.find(num);
			if(iter==pntsToFace.end())
			{
				for(auto iter=pntsToEdge.begin(); iter!=pntsToEdge.end(); iter++)
					Info<<iter->first<<"/"<<iter->second<<Foam::endl;
				Info<<"---------"<<Foam::endl;
				for(auto iter=pntsToFace.begin(); iter!=pntsToFace.end(); iter++)
					Info<<iter->first<<"/"<<iter->second<<Foam::endl;
				FatalErrorInFunction<<"Error num:"<<num<< exit(FatalError);
			}
			mc33Cube.facePermutation[faceInd] = iter->second;
		}
		std::bitset<6> refFaces;
		for(short unsigned int permFace : mc33Cube.facePermutation)
			refFaces[permFace] = true;
		if(!refFaces.all())
		{
			Info<<"[";
			for(label perm : mc33Cube.facePermutation)
				Info<<perm<<",";
			Info<<"]"<<Foam::endl;
			std::cout<<refFaces<<std::endl;
			FatalErrorInFunction<<"Mismatch face permutations"<< exit(FatalError);
		}
		
		std::vector<std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>>& triangles = mc33Cube.cutTriangles;
		for(auto iter=triangles.begin(); iter!=triangles.end(); )
		{
			std::tuple<std::uint8_t,std::uint8_t,std::uint8_t>& triangle = *iter;
			if(detectInFaceTriangles(triangle))
			{
				
				if(mc33Cube.cubeCase!=c612 && mc33Cube.cubeCase!=c742)
				{
					Info<<"("<<std::get<0>(triangle)<<","<<std::get<1>(triangle)<<","<<std::get<2>(triangle)<<")"<<Foam::endl;
					FatalErrorInFunction<<"Triangles in face in case where not possible: "<<mc33Cube.cubeCase<< exit(FatalError);
				}
				iter = triangles.erase(iter);
			}
			else
				iter++;
		}
		
		return mc33Cube;
	}
	else
	{
		mc33Cube.cell=-1;
		return mc33Cube;
	}
	
}

void Foam::cutCellFvMesh::MC33::mc33Cube_print
(
	const MC33Cube& cube,
	const labelList& oldToNewPointInd
)
{
	Info<<"Cell:"<<cube.cell<<Foam::endl;
	Info<<"vertices:"<<cube.vertices<<Foam::endl;
	Info<<"new vertices: 8(";
	for(label vert : cube.vertices)
		Info<<" "<<oldToNewPointInd[vert]<<" ";
	Info<<")"<<Foam::endl;
	List<label> permVertices(cube.vertices.size(),-1);
	Info<<"new perm vertices: 8(";
	for(label i=0; i<permVertices.size(); i++)
	{
		label perm = cube.pointPermutation[i];
		Info<<" "<<oldToNewPointInd[cube.vertices[perm]]<<" ";
	}
	Info<<")"<<Foam::endl;
	Info<<"bitPattern:"<<cube.bitPattern<<Foam::endl;
	//Info<<"edges:"<<cube.edges<<Foam::endl;
	//Info<<"edgeGlobalInd:"<<cube.edgeGlobalInd<<Foam::endl;
	Info<<"cutEdgeVerticeIndex:"<<cube.cutEdgeVerticeIndex<<Foam::endl;
	Info<<"centerPointInd:"<<cube.centerPointInd<<Foam::endl;
	Info<<"origFaces:"<<cube.origFaces<<Foam::endl;
	Info<<"cubeCase:"<<cube.cubeCase<<Foam::endl;
	Info<<"permutationTableIndex:"<<cube.permutationTableIndex<<Foam::endl;
	Info<<"redMarkIsPlusSide:"<<cube.redMarkIsPlusSide<<Foam::endl;
	Info<<"pointPermutation: (";
	for(int pointPerm : cube.pointPermutation)
		Info<<" "<<pointPerm;
	Info<<" )"<<Foam::endl;
	Info<<"edgePermutation: (";
	for(int edgePerm : cube.edgePermutation)
		Info<<" "<<edgePerm;
	Info<<" )"<<Foam::endl;
	Info<<"facePermutation: (";
	for(int facePerm : cube.facePermutation)
		Info<<" "<<facePerm;
	Info<<" )"<<Foam::endl;
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

