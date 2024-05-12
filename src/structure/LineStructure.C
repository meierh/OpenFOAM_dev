#include "LineStructure.H"

Foam::LagrangianMarker::LagrangianMarker
(
    const dynamicRefineFvMesh& mesh,
    const label rodNumber
    const ActiveRodMesh::rodCosserat* baseRod,
    const scalar markerParameter
):
mesh(mesh),
rodNumber(rodNumber),
baseRod(baseRod),
markerParameter(markerParameter),
markerPosition(LineStructure::evaluateRodPos(baseRod,markerParameter)),
markerCell(mesh.findCell(marker.markerPosition))
{
    computeSupport();
    minMaxSupportWidth();
    dilationFactors();
}

std::string Foam::LagrangianMarker::to_string()
{
    return std::to_string(markerParameter);
}

vector Foam::LagrangianMarker::getMarkerVelocity()
{
    return vector(0,0,0);
}

scalar Foam::LagrangianMarker::getMarkerTemperature()
{
    return 1;
}

void Foam::LineStructure::computeSupport
(
    label iterations
)
{
    const cellList& cellList = mesh.cells();
    const faceList& facesList = mesh.faces();
    const labelList& owners = mesh.owner();
    const labelList& neighbours = mesh.neighbour();
    std::unordered_set<label> supportCells;
    this->supportCells.resize(0);
    if(markerCell!=-1)
    {
        if(markerCell<0 || markerCell>=cellList.size())
            FatalErrorInFunction<<"Invalid cell index"<< exit(FatalError);
        supportCells.insert(markerCell);
        
        for(label iter=0; iter<iterations; iter++)
        {
            DynamicList<label> newCells;
            for(auto cellIter=supportCells.begin(); cellIter!=supportCells.end(); cellIter++)
            {
                label cellInd = *cellIter;
                const cell& thisCell = cellList[cellInd];
                for(label faceInd : thisCell)
                {
                    supportCells.insert(owners[faceInd]);
                    if(faceInd<neighbours.size())
                        supportCells.insert(neighbours[faceInd]);                    
                }
            }
            supportCells.insert(newCells.begin(),newCells.end());
        }
    }
    for(auto iterCells=supportCells.begin(); iterCells!=supportCells.end(); iterCells++)
    {
        this->supportCells.append(*iterCells);
    }
}

void Foam::LagrangianMarker::minMaxSupportWidth()
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    scalar min = std::numeric_limits<scalar>::min();
    scalar max = std::numeric_limits<scalar>::max();
    vector minSpan(max,max,max);
    vector maxSpan(min,min,min);

    for(label i=0; i<supportCells.size(); i++)
    {
        vector cellCentreA = cells[supportCells[i]].centre(faces,points);
        for(label j=0; j<supportCells.size(); j++)
        {
            if(i!=j)
            {
                vector cellCentreB = cells[supportCells[j]].centre(faces,points);
                vector conn = cellCentreA-cellCentreB;
                for(label d=0; d<3; d++)
                {
                    conn[d] = std::abs(conn[d]);
                    minSpan[d] = std::min(minSpan[d],conn[d]);
                    maxSpan[d] = std::max(maxSpan[d],conn[d]);
                }
            }
        }
    }
    h_plus = maxSpan;
    h_minus = minSpan;
}

void Foam::LagrangianMarker::dilationFactors()
{
    vector minSpan = h_minus;
    vector maxSpan = h_plus;
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    dilation = 5.0/6.0 * maxSpan + 1.0/6.0 * minSpan + eps
}

Foam::FieldMarkerStructureInteraction::FieldMarkerStructureInteraction
(
    dynamicRefineFvMesh& mesh,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
mesh(mesh),
h(std::cbrt(mesh.cells()[0].mag(mesh.points(),mesh.faces()))),
modusFieldToMarker(modusFieldToMarker),
modusMarkerToField(modusMarkerToField)
{
}

Foam::scalar Foam::FieldMarkerStructureInteraction::phiFunction(Foam::scalar r)
{
    //Info<<"phiFunction -- r:"<<r<<Foam::endl;
    
    scalar abs_r = std::abs(r);
    if(abs_r < 0.5)
    {
        scalar result = (1.0/3.0)*(1+std::sqrt(-3*(abs_r*abs_r)+1));
        //Info<<"abs_r:"<<abs_r<<" -> "<<result<<Foam::endl;
        return result;
    }
    else if(abs_r <= 1.5)
    {
        scalar result = (1.0/6.0)*(5.0-3.0*abs_r-std::sqrt(-3.0*((1-abs_r)*(1-abs_r))+1));
        //Info<<"abs_r:"<<abs_r<<" -> "<<result<<Foam::endl;
        return result;        
    }
    else
    {
        //Info<<"abs_r:"<<abs_r<<" -> "<<0<<Foam::endl;
        return 0;
    }
}

Foam::scalar Foam::FieldMarkerStructureInteraction::deltaDirac
(
    Foam::vector X,
    Foam::vector x,
    Foam::scalar h
)
{
    return deltaDirac(X,x,vector(h,h,h));
}

Foam::scalar Foam::FieldMarkerStructureInteraction::deltaDirac
(
    Foam::vector X,
    Foam::vector x,
    Foam::vector h
)
{   
    vector sigma_d;
    for(label dim=0; dim<3; dim++)
    {
        scalar X_i_x_i = X[dim]-x[dim];
        scalar r = X_i_x_i / h[dim];
        sigma_d[dim] = phiFunction(r);
    }
    scalar deltaDir = sigma_d[0]*sigma_d[1]*sigma_d[2];
    deltaDir /= (h[0]*h[1]*h[2]);
    //Info<<h<<"  X:"<<X<<"  x:"<<x<<" -> "<<deltaDir<<Foam::endl;
    return deltaDir;
}

scalar Foam::FieldMarkerStructureInteraction::correctedDeltaDirac
(
    vector X,
    vector x,
    scalar h,
    std::array<scalar,10> b
)
{
    vector conn = x-X;
    scalar correctionFactor =   b[0] + 
                                conn[0]*b[1] + conn[1]*b[2] + conn[2]*b[3] +
                                conn[0]*conn[1]*b[4] + conn[1]*conn[2]*b[5] + conn[2]*conn[0]*b[6] +
                                conn[0]*conn[0]*b[7] + conn[1]*conn[1]*b[8] + conn[2]*conn[2]*b[9];
    return correctionFactor*deltaDirac(X,x,h);
}

scalar Foam::FieldMarkerStructureInteraction::computeMoment
(
    const LagrangianMarker& marker,
    vector indices
)
{
    X = marker.markerPosition;
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    scalar moment = 0;
    for(label i=0; i<marker.supportCells.size(); i++)
    {
        label cellInd = supportCells[i];
        vector x = cells[cellInd].centre(faces,points);
        vector conn = x-X;
        vector coeff(1,1,1);
        for(label dim=0; dim<3; dim++)
        {
            for(label index=0; index<indices[dim]; index++)
            {
                coeff[dim]*=coeff[dim];
            }
        }
        scalar volume = cells[*suppIter].mag(faces,points);
        scalar dirac = deltaDirac(X,x,marker.getDilation());
        moment += coeff[0]*coeff[1]*coeff[2]*dirac*volume;
    }
    return moment
};

std::unique_ptr<gismo::gsMatrix<scalar>> Foam::FieldMarkerStructureInteraction::computeMomentMatrix
(
    const LagrangianMarker& marker
)
{
    auto moments3DPtr = std::unique_ptr<gismo::gsMatrix<scalar>>(new gismo::gsMatrix<scalar>(10,10);
    gismo::gsMatrix<scalar>& moments3D = *moments3DPtr;
    
    //Diagonal
    moments3D[0,0] = computeMoment(marker,vector(0,0,0));
    moments3D[1,1] = computeMoment(marker,vector(2,0,0));
    moments3D[2,2] = computeMoment(marker,vector(0,2,0));
    moments3D[3,3] = computeMoment(marker,vector(0,0,2));
    moments3D[4,4] = computeMoment(marker,vector(2,2,0));
    moments3D[5,5] = computeMoment(marker,vector(0,2,2));
    moments3D[6,6] = computeMoment(marker,vector(2,0,2));
    moments3D[7,7] = computeMoment(marker,vector(4,0,0));
    moments3D[8,8] = computeMoment(marker,vector(0,4,0));
    moments3D[9,9] = computeMoment(marker,vector(0,0,4));
    
    //Row / Column 0
    moments3D[0,1] = moments3D[1,0] = computeMoment(marker,vector(1,0,0));
    moments3D[0,2] = moments3D[2,0] = computeMoment(marker,vector(0,1,0));
    moments3D[0,3] = moments3D[3,0] = computeMoment(marker,vector(0,0,1));
    moments3D[0,4] = moments3D[4,0] = computeMoment(marker,vector(1,1,0));
    moments3D[0,5] = moments3D[5,0] = computeMoment(marker,vector(0,1,1));
    moments3D[0,6] = moments3D[6,0] = computeMoment(marker,vector(1,0,1));
    moments3D[0,7] = moments3D[7,0] = computeMoment(marker,vector(2,0,0));
    moments3D[0,8] = moments3D[8,0] = computeMoment(marker,vector(0,2,0));
    moments3D[0,9] = moments3D[9,0] = computeMoment(marker,vector(0,0,2));

    //Row / Column 1
    moments3D[1,2] = moments3D[2,1] = computeMoment(marker,vector(1,1,0));
    moments3D[1,3] = moments3D[3,1] = computeMoment(marker,vector(1,0,1));
    moments3D[1,4] = moments3D[4,1] = computeMoment(marker,vector(2,1,0));
    moments3D[1,5] = moments3D[5,1] = computeMoment(marker,vector(1,1,1));
    moments3D[1,6] = moments3D[6,1] = computeMoment(marker,vector(2,0,1));
    moments3D[1,7] = moments3D[7,1] = computeMoment(marker,vector(3,0,0));
    moments3D[1,8] = moments3D[8,1] = computeMoment(marker,vector(1,2,0));
    moments3D[1,9] = moments3D[9,1] = computeMoment(marker,vector(1,0,2));
    
    //Row / Column 2
    moments3D[2,3] = moments3D[3,2] = computeMoment(marker,vector(0,1,1));
    moments3D[2,4] = moments3D[4,2] = computeMoment(marker,vector(1,2,0));
    moments3D[2,5] = moments3D[5,2] = computeMoment(marker,vector(0,2,1));
    moments3D[2,6] = moments3D[6,2] = computeMoment(marker,vector(1,1,1));
    moments3D[2,7] = moments3D[7,2] = computeMoment(marker,vector(2,1,0));
    moments3D[2,8] = moments3D[8,2] = computeMoment(marker,vector(0,3,0));
    moments3D[2,9] = moments3D[9,2] = computeMoment(marker,vector(0,1,2));
    
    //Row / Column 3
    moments3D[3,4] = moments3D[4,3] = computeMoment(marker,vector(1,1,1));
    moments3D[3,5] = moments3D[5,3] = computeMoment(marker,vector(0,1,2));
    moments3D[3,6] = moments3D[6,3] = computeMoment(marker,vector(1,0,2));
    moments3D[3,7] = moments3D[7,3] = computeMoment(marker,vector(2,0,1));
    moments3D[3,8] = moments3D[8,3] = computeMoment(marker,vector(0,2,1));
    moments3D[3,9] = moments3D[9,3] = computeMoment(marker,vector(0,0,3));

    //Row / Column 4
    moments3D[3,4] = moments3D[4,3] = computeMoment(marker,vector(1,1,1));
    moments3D[3,5] = moments3D[5,3] = computeMoment(marker,vector(0,1,2));
    moments3D[3,6] = moments3D[6,3] = computeMoment(marker,vector(1,0,2));
    moments3D[3,7] = moments3D[7,3] = computeMoment(marker,vector(2,0,1));
    moments3D[3,8] = moments3D[8,3] = computeMoment(marker,vector(0,2,1));
    moments3D[3,9] = moments3D[9,3] = computeMoment(marker,vector(0,0,3));

    //Row / Column 5
    moments3D[4,5] = moments3D[5,4] = computeMoment(marker,vector(1,2,1));
    moments3D[4,6] = moments3D[6,4] = computeMoment(marker,vector(2,1,1));
    moments3D[4,7] = moments3D[7,4] = computeMoment(marker,vector(3,1,0));
    moments3D[4,8] = moments3D[8,4] = computeMoment(marker,vector(1,3,0));
    moments3D[4,9] = moments3D[9,4] = computeMoment(marker,vector(1,1,2));
    
    //Row / Column 6
    moments3D[5,6] = moments3D[6,5] = computeMoment(marker,vector(1,1,2));
    moments3D[5,7] = moments3D[7,5] = computeMoment(marker,vector(2,1,1));
    moments3D[5,8] = moments3D[8,5] = computeMoment(marker,vector(0,3,1));
    moments3D[5,9] = moments3D[9,5] = computeMoment(marker,vector(0,1,3));
    
    //Row / Column 7
    moments3D[6,7] = moments3D[7,6] = computeMoment(marker,vector(3,0,1));
    moments3D[6,8] = moments3D[8,6] = computeMoment(marker,vector(1,2,1));
    moments3D[6,9] = moments3D[9,6] = computeMoment(marker,vector(1,0,3));
    
    //Row / Column 8
    moments3D[7,8] = moments3D[8,7] = computeMoment(marker,vector(2,2,0));
    moments3D[7,9] = moments3D[9,7] = computeMoment(marker,vector(2,0,2));
    
    //Row / Column 9
    moments3D[8,9] = moments3D[9,8] = computeMoment(marker,vector(0,2,2));
    
    return moments3DPtr;
}

std::unique_ptr<std::array<scalar,10>> Foam::FieldMarkerStructureInteraction::rescalingDiagonal
(
    const LagrangianMarker& marker
)
{
    auto diagPtr = std::unique_ptr<std::array<scalar,10>>(new std::array<scalar,10>());
    std::array<scalar,10>& diag = *diagPtr;
    diag[0] = 1;
    diag[1] = 1/dilation[0];
    diag[2] = 1/dilation[1];
    diag[3] = 1/dilation[2];
    diag[4] = 1/(dilation[0]*dilation[1]);
    diag[5] = 1/(dilation[1]*dilation[2]);
    diag[6] = 1/(dilation[0]*dilation[2]);
    diag[7] = 1/(dilation[0]*dilation[0]);
    diag[8] = 1/(dilation[1]*dilation[1]);
    diag[9] = 1/(dilation[2]*dilation[2]);
    return diagPtr;
}

std::unique_ptr<std::array<scalar,10>> Foam::FieldMarkerStructureInteraction::computeCorrectionWeights
(
    const LagrangianMarker& marker
)
{
    auto bPtr = std::unique_ptr<std::array<scalar,10>>(new std::array<scalar,10>());
    std::array<scalar,10>& b = *bPtr;
    
    std::unique_ptr<gismo::gsMatrix<scalar>> M = computeMomentMatrix(marker);
    std::unique_ptr<std::array<scalar,10>> H = rescalingDiagonal(marker);
    gismo::gsMatrix<scalar> e(10,1),c(10,1);
    for(label i=0; i<10; i++)
    {
        (*M)(i,i) = (*M)(i,i) * (*H)[i];
        e(i,0) = 0;
        c(i,0) = 0;
    }
    e(0,0) = 1

    gismo::gsGMRes<scalar> HM(*M);
    HM.solve(e,c);
    
    for(label i=0; i<10; i++)
        b[i] = c(i,0) * (*H)[i];

    return bPtr;
}


void Foam::FieldMarkerStructureInteraction::scatterNurbs
(
    std::pair<gsNurbs<scalar>,label> in,
    std::pair<gsNurbs<scalar>,label>& out
)
{
    List<List<scalar>> nurbsData(7);
    // degree
    // gsKnotVector
    // weights
    // cols // nbrControlPoints
    // rows // dimensionControlPoints
    // coefficients
    // number
    if(Pstream::master())
    {
        gismo::gsKnotVector<scalar> knots = in.first.knots();
        label degree = knots.degree();
        nurbsData[0] = List<scalar>(1);
        nurbsData[0][0] = degree;
        nurbsData[1] = List<scalar>(knots.size());
        for(label i=0;i<knots.size();i++)
            nurbsData[1][i] = knots[i];
        
        gismo::gsMatrix<scalar> weights = in.first.weights();
        if(weights.rows()!=1)
            FatalErrorInFunction<<"Wrong weight dimensions"<< exit(FatalError);
        nurbsData[2] = List<scalar>(weights.cols());
        for(label i=0;i<weights.cols();i++)
            nurbsData[2][i] = weights(0,i);

        gismo::gsMatrix<scalar> coefs = in.first.coefs();
        nurbsData[3] = List<scalar>(1);
        nurbsData[3][0] = coefs.cols();
        nurbsData[4] = List<scalar>(1);
        nurbsData[4][0] = coefs.rows();
        nurbsData[5] = List<scalar>(coefs.cols()*coefs.rows());
        label index = 0;
        for(int n=0;n<coefs.cols();n++)
        {
            for(int d=0;d<coefs.rows();d++)
            {
                nurbsData[5][index] = coefs(d,n);
                index++;
            }
        }
        
        nurbsData[6] = List<scalar>(1);
        nurbsData[6][0] = in.second;
    }
    Pstream::scatter(nurbsData);
    if(Pstream::master())
    {
        out = in;
    }
    else
    {
        std::vector<double> knotContainer(nurbsData[1].size());
        for(label i=0;i<nurbsData[1].size();i++)
            knotContainer[i] = nurbsData[1][i];
        gismo::gsKnotVector<scalar> knots(knotContainer,nurbsData[0][0]);

        gismo::gsMatrix<scalar> weights(nurbsData[2].size(),1);
        for(label i=0;i<nurbsData[2].size();i++)
            weights(i,0) = nurbsData[2][i];

        gismo::gsMatrix<scalar> coefs(nurbsData[4][0],nurbsData[3][0]);
        nurbsData[5] = List<scalar>(coefs.cols()*coefs.rows());
        label index = 0;
        for(int n=0;n<coefs.cols();n++)
        {
            for(int d=0;d<coefs.rows();d++)
            {
                coefs(d,n) = nurbsData[5][index];;
                index++;
            }
        }
       
        out.first = gismo::gsNurbs<double>(knots,weights,coefs);
        
        out.second = nurbsData[6][0];
    }
}

Foam::LineStructure::LineStructure
(
    dynamicRefineFvMesh& mesh,
    const List<scalar> crossSecArea
):
Structure(mesh,mesh.time()/*,alpha,T,p,U,nu*/),
crossSecArea(crossSecArea)
{
    
}

void Foam::LineStructure::transferMarkers(FieldMarkerStructureInteraction& connector)
{
    if(!myMesh)
        FatalErrorInFunction<<"Rod Mesh not set!"<<exit(FatalError);
    if(myMesh->m_nR!=rodMarkers.size())
    {
        rodMarkers.resize(myMesh->m_nR);
    }
    connector.markers.resize(0);
    
    if(crossSecArea.size()!=rodMarkers.size())
        FatalErrorInFunction<<"Mismatch in size of crossSecArea and rodMarkers"<<exit(FatalError);
    
    for(uint rodIndex=0; rodIndex<rodMarkers.size(); rodIndex++)
    {
        std::unique_ptr<std::vector<LagrangianMarker>>& oneRodMarkers = rodMarkers[rodIndex];
        if(!oneRodMarkers)
        {
            myMesh->m_Rods[rodIndex]->m_Rot.setCoefs(myMesh->m_Rods[rodIndex]->m_Rot_coefs.transpose());
            oneRodMarkers = constructMarkerSet(rodIndex,myMesh->m_Rods[rodIndex],crossSecArea[rodIndex]);
        }
        for(LagrangianMarker& marker : *oneRodMarkers)
        {
            connector.markers.push_back(&marker);
        }
    }
}

std::unique_ptr<std::vector<scalar>> Foam::LineStructure::createSpacedPointsOnRod
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar spacing
)
{
    auto pointsPtr = std::unique_ptr<std::vector<scalar>>(new std::vector<scalar>());
    std::vector<scalar>& pointsVec = *pointsPtr;
    
    std::list<scalar> points;
    points.push_back(oneRod->m_Curve.domainStart());
    points.push_back(oneRod->m_Curve.domainEnd());
    
    bool refined=true;
    while(refined)
    {
        //Info<<"Refine "<<refinementCount<<" size:"<<markers.size()<<Foam::endl;
        refined = false;
        bool cond = true;
        auto pntsIter0 = points.begin();
        auto pntsIter1 = ++(points.begin());
        for( ; pntsIter1!=points.end() ; )
        {
            scalar dist = distance(oneRod,*pntsIter0,*pntsIter1);
            if(dist>spacing)
            {
                scalar middlePar = 0.5*(*pntsIter0 + *pntsIter1);
                auto inserted = markers.insert(pntsIter1,middlePar);
                refined=true;
            }
            markersIter0 = markersIter1;
            markersIter1++;
        }
    }
    
    for(scalar para : points)
        pointsVec.push_back(para);
    
    return pointsPtr;
}

std::pair<double,double> Foam::LineStructure::minMaxSpan
(
    const cell& thisCell
)
{
    std::vector<double> spans;
    cellDistances(thisCell,spans);
    std::sort(spans.begin(),spans.end());
    if(spans.size()<24)
        FatalErrorInFunction<<"Cell can not have less than 24 spans"<< exit(FatalError);
    return {spans.front(),spans.back()};
}

void Foam::LineStructure::cellDistances
(
    const cell& thisCell,
    std::vector<double>& spans
)
{
    const pointField& cellPoints = thisCell.points(mesh.faces(),mesh.points());
    for(label pntI=0; pntI<cellPoints.size(); pntI++)
    {
        for(label pntJ=0; pntJ<cellPoints.size(); pntJ++)
        {
            if(pntI!=pntJ)
            {
                vector distanceVec = cellPoints[pntI]-cellPoints[pntJ];
                scalar distanceVal = Foam::mag(distanceVec);
                spans.push_back(distanceVal);
            }
        }
    }    
}

scalar Foam::LineStructure::supportDomainMinSize
(
    label cellInd
)
{
    //Info<<"     supportDomainMinSize:"<<cellInd<<Foam::endl;
    const cellList& cells = mesh.cells();
    
    double minSuppSize = std::numeric_limits<double>::max(); 
    std::unordered_set<label> neighbourhoodCells;
    getSupportDomain(cellInd,neighbourhoodCells);
    //Info<<"     neighbourhoodCells.size():"<<neighbourhoodCells.size()<<Foam::endl;
    for(auto cellIter=neighbourhoodCells.begin(); cellIter!=neighbourhoodCells.end(); cellIter++)
    {
        //Info<<"     cell:"<<*cellIter<<Foam::endl;
        const cell& suppCell = cells[*cellIter];
        double minSpan = minMaxSpan(suppCell).first;
        minSuppSize = std::min<double>(minSuppSize,minSpan);
    }
    //Info<<"     supportDomainMinSize done"<<Foam::endl;
    return minSuppSize;
}

scalar Foam::LineStructure::supportDomainMinSize
(
    const std::unordered_set<label>& supportDomainCells
)
{    
    const cellList& cells = mesh.cells();    
    double minSuppSize = std::numeric_limits<double>::max(); 
    for(auto cellIter=supportDomainCells.begin(); cellIter!=supportDomainCells.end(); cellIter++)
    {
        const cell& suppCell = cells[*cellIter];
        double minSpan = minMaxSpan(suppCell).first;
        minSuppSize = std::min<double>(minSuppSize,minSpan);
    }
    return minSuppSize;
}

std::unique_ptr<std::vector<LagrangianMarker>> Foam::LineStructure::constructMarkerSet
(
    label rodNumber,
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar crossSecArea
)
{
    Info<<"constructMarkerSet"<<Foam::endl;
    
    std::list<LagrangianMarker> markers;
    
    // Insert start marker
    scalar startPar = oneRod->m_Curve.domainStart();
    markers.push_back(createLagrangianMarker(startPar,rodNumber,oneRod));
    Info<<"Start marker created: "<<startPar<<Foam::endl;
    Info<<"markers.back():"<<markers.back().markerParameter<<" - "<<markers.back().markerPosition<<Foam::endl;
    
    // Insert end marker
    scalar endPar = oneRod->m_Curve.domainEnd();
    markers.push_back(createLagrangianMarker(endPar,rodNumber,oneRod));
    Info<<"End marker created: "<<endPar<<Foam::endl;
    Info<<"markers.back():"<<markers.back().markerParameter<<" - "<<markers.back().markerPosition<<Foam::endl;
    
    bool refined=true;
    label refinementCount = 0;
    while(refined)
    {
        //Info<<"Refine "<<refinementCount<<" size:"<<markers.size()<<Foam::endl;
        refined = false;
        bool cond = true;
        auto markersIter0 = markers.begin();
        auto markersIter1 = ++(markers.begin());
        for( ; markersIter1!=markers.end() ; )
        {
            bool subdivide = doSubdivision(*markersIter0, *markersIter1);
            if(refinementCount<minRefinement)
                subdivide = true;
            if(subdivide)
            {
                scalar middlePar = 0.5*(markersIter0->markerParameter+markersIter1->markerParameter);
                auto inserted = markers.insert(markersIter1,createLagrangianMarker(middlePar,rodNumber,oneRod));
                refined=true;
                //Info<<"   Subdivision:"<<markersIter0->markerParameter<<"--"<<markersIter1->markerParameter<<" / "<<middlePar<<Foam::endl;
            }
            else
            {
                //Info<<"   Non Subdivision:"<<markersIter0->markerParameter<<"--"<<markersIter1->markerParameter<<Foam::endl;
            }
            markersIter0 = markersIter1;
            markersIter1++;
        }
        refinementCount++;
    }
        
    std::list<LagrangianMarker>::iterator iterPrev = markers.end();
    std::list<LagrangianMarker>::iterator iterNext;
    for(auto iter=markers.begin(); iter!=markers.end(); iter++)
    {
        iterNext = iter;
        iterNext++;
        
        scalar spanStart;
        if(iterPrev!=markers.end())
        {
            spanStart = iterPrev->markerParameter;
            spanStart = spanStart + iter->markerParameter;
            spanStart /= 2;
        }
        else
            spanStart = iter->markerParameter;

        scalar spanEnd;
        if(iterNext!=markers.end())
        {            
            spanEnd = iterNext->markerParameter;
            spanEnd = spanEnd + iter->markerParameter;
            spanEnd /= 2;
        }
        else
            spanEnd = iter->markerParameter;
        
        if(iterPrev==markers.end() && iterNext==markers.end())
            FatalErrorInFunction<<"Marker with no predecessor and no succesor"<<exit(FatalError);
        
        scalar span = spanEnd-spanStart;
        iter->markerVolume = span*crossSecArea;
        
        iterPrev = iter;
    }

    std::unique_ptr<std::vector<LagrangianMarker>> markersPtr(new std::vector<LagrangianMarker>());
    for(auto iter=markers.begin(); iter!=markers.end(); iter++)
    {
        Info<<iter->markerParameter<<" -- "<<iter->markerPosition<<" "<<iter->markerCell<<Foam::endl;
        if(iter->markerCell!=-1)
            markersPtr->push_back(*iter);
    }
    
    Info<<"Created "<<markersPtr->size()<<" markers!"<<Foam::endl;
    return markersPtr;
}

Foam::vector Foam::LineStructure::evaluateRodPos
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter
)
{
    vector r;
    rodEval(oneRod,parameter,r);
    return r;
}

/*
Foam::vector Foam::LineStructure::evaluateRodDeriv
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter
)
{
    gsMatrix<double,1,1> parMat(1,1);
    parMat.at(0) = parameter;
    
    gsMatrix<double> basePntVec;
    gsNurbs<double> base = oneRod->m_Curve;
    base.deriv_into(parMat,basePntVec);

    gsMatrix<double> defPntVec;
    gsNurbs<double> def = oneRod->m_Def;    
    def.deriv_into(parMat,defPntVec);
    
    gsMatrix<double> resultMat = basePntVec+defPntVec;
    return vector(resultMat(0,0),resultMat(0,1),resultMat(0,2));    
}
*/

Foam::LagrangianMarker Foam::LineStructure::createLagrangianMarker
(
    scalar markerParameter,
    label rodNumber,
    const ActiveRodMesh::rodCosserat* oneRod
)
{
    //Info<<" createLagrangianMarker"<<Foam::endl;
    LagrangianMarker marker(markerParameter,evaluateRodPos(oneRod,markerParameter),rodNumber,oneRod);
    label cellOfMarker = mesh.findCell(marker.markerPosition);
    if(cellOfMarker!=-1)
        getSupportDomain(cellOfMarker,marker.supportCells);
    marker.markerCell = mesh.findCell(marker.markerPosition);
    /*
    for(auto iter=marker.supportCells.begin(); iter!=marker.supportCells.end(); iter++)
        if(mesh.pointInCell(marker.markerPosition,*iter))
        {
            Info<<*iter<<":"<<marker.markerPosition<<Foam::endl;
            if(marker.markerCell!=-1)
                FatalErrorInFunction<<"Double assignment marker cell"<<exit(FatalError);
            marker.markerCell = *iter;
        }
    */
    return marker;
}

template<typename T>
static T integrateRodwise
(
    const ActiveRodMesh::rodCosserat* oneRod,
    std::function<T(scalar)> function
)
{
    return integrateRodwise<T>(oneRod,oneRod->knots().first(),oneRod->knots().last(),function);
}

template<typename T>
static T integrateRodwise
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parA,
    scalar parB,
    std::function<T(scalar)> function
)
{
    if(parA>parB)
        FatalErrorInFunction<<"Integration from higher to lower value"<<exit(FatalError);
    
    std::set<scalar> knotSet;
    for(scalar knot : oneRod->knots())
        knotSet.insert(knot);
    
    label degree = oneRod->knots().degree();
    auto start = knotSet.insert(parA);
    auto end = knotSet.insert(parB);
    
    label numberOfAbcissa = (degree+1)*(degree+1);
    
    T totalValue = Foam::zero();
    
    std::set<scalar>::iterator currNode = start.first;
    std::set<scalar>::iterator nextNode = ++(start.first);
    for(; currNode!=knotSet.end() && nextNode!=knotSet.end() && currNode!=end.first ; )
    {
        scalar startPar = *currNode;
        scalar endPar = *nextNode;
        scalar dist = endPar-startPar;
        scalar step = dist/numberOfAbcissa;
        scalar initialStep = startPar+step/2;
        T summedValue = Foam::zero();
        for(label i=0;i<numberOfAbcissa;i++)
        {
            summedValue += function(initialStep)*step;
            initialStep+=step;
        }
        totalValue+=summedValue;
    }
    return totalValue;
}


scalar Foam::LineStructure::distance
(
    const ActiveRodMesh::rodCosserat* oneRod,   
    scalar parA,
    scalar parB
)
{
    std::function<scalar(scalar)> curveLen = 
    [rodPtr = oneRod](scalar parameter)
    {
        vector derivRod = rodDerivEval(rodPtr,parameter);
        return std::sqrt(derivRod&derivRod);
    };    
    return integrateRodwise<scalar>(oneRod,parA,parB,curveLen);
}

scalar Foam::LineStructure::distance
(
    const LagrangianMarker& A,
    const LagrangianMarker& B
)
{
    if(A.baseRod!=B.baseRod)
        FatalErrorInFunction<<"Distance can not be computed between points on different rods"<<exit(FatalError);
    
    return distance(A.baseRod,A.markerParameter,B.markerParameter);
}

Foam::scalar Foam::LineStructure::evaluateRodArcLen
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parA,
    scalar parB
)
{
    vector parAVec = evaluateRodPos(oneRod,parA);
    vector parBVec = evaluateRodPos(oneRod,parB);
    vector connec = parAVec-parBVec;
    return Foam::mag(connec);
}

bool Foam::LineStructure::doSubdivision
(
    const LagrangianMarker& smallerSide,
    const LagrangianMarker& largerSide
)
{
    bool supportDomainOverlap = false;
    for(auto iterSm=smallerSide.supportCells.begin();
        iterSm!=smallerSide.supportCells.end(); iterSm++)
    {
        for(auto iterLa=largerSide.supportCells.begin();
            iterLa!=largerSide.supportCells.end(); iterLa++)
        {
            
            if(*iterSm == *iterLa)
                supportDomainOverlap = true;
        }
    }
    bool smallerSideHasSupp = smallerSide.supportCells.size()!=0;
    bool largerSideHasSupp = largerSide.supportCells.size()!=0;
    bool bothSidesHaveSupp = largerSideHasSupp&&smallerSideHasSupp;
    
    scalar distOnNurbs = distance(smallerSide,largerSide);
    scalar distDirect = Foam::mag(smallerSide.markerPosition-largerSide.markerPosition);
    
    scalar smallerSideDomainMinSize = std::numeric_limits<scalar>::max();
    if(smallerSideHasSupp)
        smallerSideDomainMinSize = supportDomainMinSize(smallerSide.supportCells);
    
    scalar largerSideDomainMinSize = std::numeric_limits<scalar>::max();
    if(largerSideHasSupp)
        largerSideDomainMinSize = supportDomainMinSize(largerSide.supportCells);
    
    scalar overallDomainMinSide = std::min(smallerSideDomainMinSize,largerSideDomainMinSize);
    
    //Info<<" doSubdivision:"<<distOnNurbs<<"|"<<distDirect<<"/"<<overallDomainMinSide<<smallerSide.markerPosition<<"->"<<largerSide.markerPosition<<" -- "<<distance(smallerSide,largerSide)<<Foam::endl;
    
    if(bothSidesHaveSupp)
    {
        if(distOnNurbs<overallDomainMinSide)
        {
            if(!supportDomainOverlap)
            {
                Info<<Foam::endl;
                std::cout<<"smallerSide.markerParameter:"<<smallerSide.markerParameter<<std::endl;
                std::cout<<"largerSide.markerParameter:"<<largerSide.markerParameter<<std::endl;
                Info<<"distOnNurbs:"<<distOnNurbs<<Foam::endl;
                Info<<"distDirect:"<<distDirect<<Foam::endl;
                Info<<"smallerSideDomainMinSize:"<<smallerSideDomainMinSize<<Foam::endl;
                Info<<"largerSideDomainMinSize:"<<largerSideDomainMinSize<<Foam::endl;
                Info<<"overallDomainMinSide:"<<overallDomainMinSide<<Foam::endl;
                Info<<"smallerSideHasSupp:"<<smallerSideHasSupp<<Foam::endl;
                Info<<"largerSideHasSupp:"<<largerSideHasSupp<<Foam::endl;
                FatalErrorInFunction<<"Dist matches but no overlap"<<exit(FatalError);
            }
            return false;
        }
        else
            return true;
    }
    else if(smallerSideHasSupp || largerSideHasSupp)
    {
        if(distOnNurbs<overallDomainMinSide)
            return false;
        else
            return true;
    }
    else
    {
        if(distOnNurbs<overallDomainMinSide)
            return false;
        else
            return true;        
    }
}
