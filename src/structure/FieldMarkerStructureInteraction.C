#include "LineStructure.H"

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
    Info<<"Completed FieldMarkerStructureInteraction setup"<<Foam::endl;
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

void Foam::FieldMarkerStructureInteraction::computeMarkerEpsilon
(
    const dynamicRefineFvMesh& mesh,
    const DynamicList<LagrangianMarker*>& markers
)
{
    std::unique_ptr<gismo::gsMatrix<scalar>> markerAreaMatrix = computeMarkerEpsilonMatrix(mesh,markers);
    
    gismo::gsMatrix<scalar> ones(markers.size(),1);
    gismo::gsMatrix<scalar> epsilon(markers.size(),1);
    for(label I=0; I<markers.size(); I++)
    {
        ones(I,0) = 1;
    }

    gismo::gsGMRes HM(*markerAreaMatrix);
    HM.solve(ones,epsilon);

    for(label I=0; I<markers.size(); I++)
    {
        markers[I]->setMarkerVolume(markers[I]->getMarkerVolume()*epsilon(I,0));
    }
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

std::unique_ptr<gismo::gsMatrix<scalar>> Foam::FieldMarkerStructureInteraction::computeMarkerEpsilonMatrix
(
    const dynamicRefineFvMesh& mesh,
    const DynamicList<LagrangianMarker*>& markers
)
{
    auto result = std::unique_ptr<gismo::gsMatrix<scalar>>
    (
        new gismo::gsMatrix<scalar>(markers.size(),markers.size())
    );
    gismo::gsMatrix<scalar>& resultM = *result;
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    for(label I=0; I<markers.size(); I++)
    {
        const LagrangianMarker& markerI = *(markers[I]);
        vector XI = markerI.getMarkerPosition();
        vector dilationI = markerI.getDilation();
        scalar dilationIMax = std::max<scalar>(dilationI[0],dilationI[1]);
        dilationIMax = std::max<scalar>(dilationIMax,dilationI[2]);
        const DynamicList<std::tuple<bool,label,label>>& supportI = markerI.getSupportCells();
        std::unordered_set<label> supportSetI;
        for(auto tupl : supportI)
            supportSetI.insert(std::get<0>(tupl));
            
        
        for(label K=0; K<markers.size(); K++)
        {
            const LagrangianMarker& markerK = *(markers[K]);
            vector XK = markerK.getMarkerPosition();
            vector dilationK = markerK.getDilation();
            scalar dilationKMax = std::max<scalar>(dilationK[0],dilationK[1]);
            dilationKMax = std::max<scalar>(dilationKMax,dilationK[2]);
            scalar markerKVol = markerK.getMarkerVolume();
            const DynamicList<std::tuple<bool,label,label>>& supportK = markerK.getSupportCells();
            
            vector distVec = XI-XK;
            scalar distMag = std::sqrt(distVec&distVec);
            
            scalar matrixEntry = 0;
            if(distMag < 10*dilationIMax || distMag < 10*dilationKMax)
            {
                for(auto cellKT : supportK)
                {
                    label cellK = std::get<2>(cellKT);
                    if(supportSetI.find(cellK)!=supportSetI.end())
                    {
                        const cell& overlapCell = cells[cellK];
                        scalar overlapCellVol = overlapCell.mag(points,faces);
                        vector overlapCellCentre = overlapCell.centre(points,faces);
                        
                        using FMSI=FieldMarkerStructureInteraction;
                        scalar weightI = FMSI::deltaDirac(XI,overlapCellCentre,dilationI);
                        scalar weightK = FMSI::deltaDirac(XK,overlapCellCentre,dilationK);

                        matrixEntry += weightK*weightI*overlapCellVol;
                    }
                }
                matrixEntry *= markerKVol;
            }
            resultM(I,K) = matrixEntry;
        }
    }
    return result;
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
    const dynamicRefineFvMesh& mesh,
    const LagrangianMarker& marker,
    vector indices
)
{
    vector X = marker.getMarkerPosition();
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    const DynamicList<std::tuple<bool,label,label>>& supportCells = marker.getSupportCells();
    
    scalar moment = 0;
    for(label i=0; i<supportCells.size(); i++)
    {
        label cellInd = std::get<2>(supportCells[i]);
        vector x = cells[cellInd].centre(points,faces);
        vector conn = x-X;
        vector coeff(1,1,1);
        for(label dim=0; dim<3; dim++)
        {
            for(label index=0; index<indices[dim]; index++)
            {
                coeff[dim]*=coeff[dim];
            }
        }
        scalar volume = cells[cellInd].mag(points,faces);
        scalar dirac = deltaDirac(X,x,marker.getDilation());
        moment += coeff[0]*coeff[1]*coeff[2]*dirac*volume;
    }
    return moment;
};

std::unique_ptr<gismo::gsMatrix<scalar>> Foam::FieldMarkerStructureInteraction::computeMomentMatrix
(
    const dynamicRefineFvMesh& mesh,
    const LagrangianMarker& marker
)
{
    auto moments3DPtr = std::unique_ptr<gismo::gsMatrix<scalar>>(new gismo::gsMatrix<scalar>(10,10));
    gismo::gsMatrix<scalar>& moments3D = *moments3DPtr;
    
    //Diagonal
    moments3D(0,0) = computeMoment(mesh,marker,vector(0,0,0));
    moments3D(1,1) = computeMoment(mesh,marker,vector(2,0,0));
    moments3D(2,2) = computeMoment(mesh,marker,vector(0,2,0));
    moments3D(3,3) = computeMoment(mesh,marker,vector(0,0,2));
    moments3D(4,4) = computeMoment(mesh,marker,vector(2,2,0));
    moments3D(5,5) = computeMoment(mesh,marker,vector(0,2,2));
    moments3D(6,6) = computeMoment(mesh,marker,vector(2,0,2));
    moments3D(7,7) = computeMoment(mesh,marker,vector(4,0,0));
    moments3D(8,8) = computeMoment(mesh,marker,vector(0,4,0));
    moments3D(9,9) = computeMoment(mesh,marker,vector(0,0,4));
    
    //Row / Column 0
    moments3D(0,1) = moments3D(1,0) = computeMoment(mesh,marker,vector(1,0,0));
    moments3D(0,2) = moments3D(2,0) = computeMoment(mesh,marker,vector(0,1,0));
    moments3D(0,3) = moments3D(3,0) = computeMoment(mesh,marker,vector(0,0,1));
    moments3D(0,4) = moments3D(4,0) = computeMoment(mesh,marker,vector(1,1,0));
    moments3D(0,5) = moments3D(5,0) = computeMoment(mesh,marker,vector(0,1,1));
    moments3D(0,6) = moments3D(6,0) = computeMoment(mesh,marker,vector(1,0,1));
    moments3D(0,7) = moments3D(7,0) = computeMoment(mesh,marker,vector(2,0,0));
    moments3D(0,8) = moments3D(8,0) = computeMoment(mesh,marker,vector(0,2,0));
    moments3D(0,9) = moments3D(9,0) = computeMoment(mesh,marker,vector(0,0,2));

    //Row / Column 1
    moments3D(1,2) = moments3D(2,1) = computeMoment(mesh,marker,vector(1,1,0));
    moments3D(1,3) = moments3D(3,1) = computeMoment(mesh,marker,vector(1,0,1));
    moments3D(1,4) = moments3D(4,1) = computeMoment(mesh,marker,vector(2,1,0));
    moments3D(1,5) = moments3D(5,1) = computeMoment(mesh,marker,vector(1,1,1));
    moments3D(1,6) = moments3D(6,1) = computeMoment(mesh,marker,vector(2,0,1));
    moments3D(1,7) = moments3D(7,1) = computeMoment(mesh,marker,vector(3,0,0));
    moments3D(1,8) = moments3D(8,1) = computeMoment(mesh,marker,vector(1,2,0));
    moments3D(1,9) = moments3D(9,1) = computeMoment(mesh,marker,vector(1,0,2));
    
    //Row / Column 2
    moments3D(2,3) = moments3D(3,2) = computeMoment(mesh,marker,vector(0,1,1));
    moments3D(2,4) = moments3D(4,2) = computeMoment(mesh,marker,vector(1,2,0));
    moments3D(2,5) = moments3D(5,2) = computeMoment(mesh,marker,vector(0,2,1));
    moments3D(2,6) = moments3D(6,2) = computeMoment(mesh,marker,vector(1,1,1));
    moments3D(2,7) = moments3D(7,2) = computeMoment(mesh,marker,vector(2,1,0));
    moments3D(2,8) = moments3D(8,2) = computeMoment(mesh,marker,vector(0,3,0));
    moments3D(2,9) = moments3D(9,2) = computeMoment(mesh,marker,vector(0,1,2));
    
    //Row / Column 3
    moments3D(3,4) = moments3D(4,3) = computeMoment(mesh,marker,vector(1,1,1));
    moments3D(3,5) = moments3D(5,3) = computeMoment(mesh,marker,vector(0,1,2));
    moments3D(3,6) = moments3D(6,3) = computeMoment(mesh,marker,vector(1,0,2));
    moments3D(3,7) = moments3D(7,3) = computeMoment(mesh,marker,vector(2,0,1));
    moments3D(3,8) = moments3D(8,3) = computeMoment(mesh,marker,vector(0,2,1));
    moments3D(3,9) = moments3D(9,3) = computeMoment(mesh,marker,vector(0,0,3));

    //Row / Column 4
    moments3D(3,4) = moments3D(4,3) = computeMoment(mesh,marker,vector(1,1,1));
    moments3D(3,5) = moments3D(5,3) = computeMoment(mesh,marker,vector(0,1,2));
    moments3D(3,6) = moments3D(6,3) = computeMoment(mesh,marker,vector(1,0,2));
    moments3D(3,7) = moments3D(7,3) = computeMoment(mesh,marker,vector(2,0,1));
    moments3D(3,8) = moments3D(8,3) = computeMoment(mesh,marker,vector(0,2,1));
    moments3D(3,9) = moments3D(9,3) = computeMoment(mesh,marker,vector(0,0,3));

    //Row / Column 5
    moments3D(4,5) = moments3D(5,4) = computeMoment(mesh,marker,vector(1,2,1));
    moments3D(4,6) = moments3D(6,4) = computeMoment(mesh,marker,vector(2,1,1));
    moments3D(4,7) = moments3D(7,4) = computeMoment(mesh,marker,vector(3,1,0));
    moments3D(4,8) = moments3D(8,4) = computeMoment(mesh,marker,vector(1,3,0));
    moments3D(4,9) = moments3D(9,4) = computeMoment(mesh,marker,vector(1,1,2));
    
    //Row / Column 6
    moments3D(5,6) = moments3D(6,5) = computeMoment(mesh,marker,vector(1,1,2));
    moments3D(5,7) = moments3D(7,5) = computeMoment(mesh,marker,vector(2,1,1));
    moments3D(5,8) = moments3D(8,5) = computeMoment(mesh,marker,vector(0,3,1));
    moments3D(5,9) = moments3D(9,5) = computeMoment(mesh,marker,vector(0,1,3));
    
    //Row / Column 7
    moments3D(6,7) = moments3D(7,6) = computeMoment(mesh,marker,vector(3,0,1));
    moments3D(6,8) = moments3D(8,6) = computeMoment(mesh,marker,vector(1,2,1));
    moments3D(6,9) = moments3D(9,6) = computeMoment(mesh,marker,vector(1,0,3));
    
    //Row / Column 8
    moments3D(7,8) = moments3D(8,7) = computeMoment(mesh,marker,vector(2,2,0));
    moments3D(7,9) = moments3D(9,7) = computeMoment(mesh,marker,vector(2,0,2));
    
    //Row / Column 9
    moments3D(8,9) = moments3D(9,8) = computeMoment(mesh,marker,vector(0,2,2));
    
    return moments3DPtr;
}

std::unique_ptr<std::array<scalar,10>> Foam::FieldMarkerStructureInteraction::rescalingDiagonal
(
    const LagrangianMarker& marker
)
{
    auto diagPtr = std::unique_ptr<std::array<scalar,10>>(new std::array<scalar,10>());
    std::array<scalar,10>& diag = *diagPtr;
    vector dilation = marker.getDilation();
    
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
    const dynamicRefineFvMesh& mesh,
    const LagrangianMarker& marker
)
{
    auto bPtr = std::unique_ptr<std::array<scalar,10>>(new std::array<scalar,10>());
    std::array<scalar,10>& b = *bPtr;
    
    std::unique_ptr<gismo::gsMatrix<scalar>> M = computeMomentMatrix(mesh,marker);
    std::unique_ptr<std::array<scalar,10>> H = rescalingDiagonal(marker);
    gismo::gsMatrix<scalar> e(10,1),c(10,1);
    for(label i=0; i<10; i++)
    {
        (*M)(i,i) = (*M)(i,i) * (*H)[i];
        e(i,0) = 0;
        c(i,0) = 0;
    }
    e(0,0) = 1;

    gismo::gsGMRes HM(*M);
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
