#include "NurbsStructureInterface.H"

Foam::NurbsStructureInterface::NurbsStructureInterface
(
    Time& runTime,
    const dimensionedScalar& alpha,
    const volScalarField& T,
    const volScalarField& p,
    const volVectorField& U,
    cutCellFvMesh& mesh,
    const dimensionedScalar nu,
    const word& IBpatchName
):
runTime(runTime),
runDirectory(runTime.rootPath()),
caseName(runTime.caseName()),
xmlPath(getXMLPath()),
nR(loadRodsFromXML()),
alpha(alpha),
T(T),
p(p),
U(U),
mesh(mesh),
IBpatchName(IBpatchName),
nu(nu),
meshPointNurbsReference(mesh.getMeshPointNurbsReference()),
Curves(mesh.getCurves())
{    
    cntOpt.ptsType = 2;
    cntOpt.ptsN = 2;
    cntOpt.csFac = 2.0;
    cntOpt.preg = 0.02 * geoR0; 
    cntOpt.kc = 0.2 * geoE0*geoR0*geoR0; 
    cntOpt.initOut = 0;
    
    folder = runDirectory+"/"+caseName;

    createNurbsStructure();
    createNurbsBoundary();
    setSolverOptions();
    
    auto cutCellBound = this->mesh.Sf().boundaryField();
    const fvBoundaryMesh& bound = mesh.boundary();
    IBPatchID = bound.findPatchID(IBpatchName);
    
    if(Curves)
        Info<<"Exists"<<Foam::endl;
    else
        Info<<"Does not exist"<<Foam::endl;
    
    createDeformationCurve();    
    assignBoundaryFacesToNurbsCurves();
    Info<<"Completed Structure setup"<<Foam::endl;
}

void Foam::NurbsStructureInterface::assignBoundaryFacesToNurbsCurves()
{
    const pointField& points = mesh.points();
    const fvBoundaryMesh& boundary = mesh.boundary();
    const fvPatch& nurbsBoundary = boundary[IBPatchID];
    const faceList& faces = mesh.faces();
    label i=0;
    label faceInd=nurbsBoundary.start();
    //Info<<"Collect faces"<<Foam::endl;
    
    nurbsToLimits.resize(Curves->size());
    for(label nurbsInd = 0; nurbsInd<Curves->size(); nurbsInd++)
    {
        const scalar minU = (*Curves)[nurbsInd].min_U();
        const scalar maxU = (*Curves)[nurbsInd].max_U();
        gismo::gsMatrix<double> collocationPts = (myMesh->m_Rods[i])->m_collPts;
        int rows = collocationPts.rows();
        if(rows!=1)
            FatalErrorInFunction<<"Row direction of Span of limits is too small"<<exit(FatalError);
        int cols = collocationPts.cols();
        nurbsToLimits[nurbsInd].resize(cols+1);
        if(nurbsToLimits[nurbsInd].size()<2)
            FatalErrorInFunction<<"Span of limits is too small"<<exit(FatalError);
        nurbsToLimits[nurbsInd].front() = minU;
        nurbsToLimits[nurbsInd].back() = maxU;
        for(int k=1;k<nurbsToLimits[nurbsInd].size()-1;k++)
        {
            double lower = collocationPts.at(k-1);
            double upper = collocationPts.at(k);
            double middle = 0.5*(lower+upper);
            nurbsToLimits[nurbsInd][k] = middle;
        }
    }

    for(label i=0;i<nurbsBoundary.size();i++,faceInd++)
    {
        face thisFace = faces[faceInd];
        for(label j=0;j<thisFace.size();j++)
        {
            if(meshPointNurbsReference[thisFace[j]].size()==0)
            {
                Info<<"pnt:"<<thisFace[j]<<Foam::endl;
                FatalErrorInFunction<<"Not having a Nurbs reference"<<exit(FatalError);
            }
            boundaryPntToFaces.insert(std::pair<label,label>(thisFace[j],faceInd));
        }
    }

    nurbsParameterToPnt = List<std::multimap<scalar,label>>(Curves->size());
    for(auto iterPnt=boundaryPntToFaces.begin(); iterPnt!=boundaryPntToFaces.end(); ++iterPnt)
    {
        label pntLabel = iterPnt->first;
        const DynamicList<cutCellFvMesh::nurbsReference>&  refs = meshPointNurbsReference[pntLabel];
        const cutCellFvMesh::nurbsReference& ref = refs[0];
        nurbsParameterToPnt[ref.nurbsInd].insert(std::pair<scalar,label>(ref.nurbsPara,pntLabel));
    }

    nurbsToLimitsFaces.resize(Curves->size());
    nurbsToLimitsFacesWeights.resize(Curves->size());
    for(int nurbsInd=0; nurbsInd<Curves->size(); nurbsInd++)
    {
        std::vector<scalar>& limits = nurbsToLimits[nurbsInd];
        std::vector<std::unordered_map<label,scalar>> facesInLimitsWithWeight;
        
        //Info<<"Create limits list"<<Foam::endl;
        facesInLimitsWithWeight.resize(limits.size()-1);
        auto iterPara = nurbsParameterToPnt[nurbsInd].lower_bound(limits[0]);
        auto iterPara2 = std::max_element(nurbsParameterToPnt[nurbsInd].begin(),nurbsParameterToPnt[nurbsInd].end());
        /*
        Info<<"limits:"<<limits.back()<<Foam::endl;
        Info<<"lower:"<<iterPara->first<<"  "<<iterPara->second<<Foam::endl;
        Info<<"upper:"<<iterPara2->first<<"  "<<iterPara2->second<<Foam::endl;
        */
        for(int i=0;i<limits.size()-1;i++)
        {
            scalar lower = limits[i];
            scalar upper = limits[i+1];
            
            for(;iterPara!=nurbsParameterToPnt[nurbsInd].end() && iterPara->first<upper;iterPara++)
            {
                label pntInd = iterPara->second;
                for(auto iterPnt=boundaryPntToFaces.find(pntInd);
                    iterPnt!=boundaryPntToFaces.end() && iterPnt->first==pntInd;
                    iterPnt++)
                {
                    label faceInd = iterPnt->second;
                    face thisFace = faces[faceInd];
                    facesInLimitsWithWeight[i][faceInd] += 1.0/(scalar)thisFace.size();
                }
            }
        }
        
        //Info<<"Insert in limits"<<Foam::endl;
        nurbsToLimits[nurbsInd] = limits;
        nurbsToLimitsFaces[nurbsInd].resize(limits.size()-1);
        nurbsToLimitsFacesWeights[nurbsInd].resize(limits.size()-1);
        for(int i=0;i<limits.size()-1;i++)
        {
            /*
            Info<<"limits:"<<i<<"/"<<limits.size();
            scalar lower = limits[i];
            scalar upper = limits[i+1];
            Info<<" ("<<lower<<","<<upper<<") "<<facesInLimitsWithWeight[i].size()<<" [";
            */
            for(auto iter=facesInLimitsWithWeight[i].begin();
                iter!=facesInLimitsWithWeight[i].end();
                iter++)
            {
                nurbsToLimitsFaces[nurbsInd][i].push_back(iter->first);
                nurbsToLimitsFacesWeights[nurbsInd][i].push_back(iter->second);
                
                //Info<<"-["<<nurbsToLimitsFaces[nurbsInd][i].back()<<";"<<nurbsToLimitsFacesWeights[nurbsInd][i].back()<<"]-";
            }
            //Info<<"]"<<Foam::endl;

        }
    }
    //Info<<"End"<<Foam::endl;
}

template<typename Tensor_Type>
std::unique_ptr<std::vector<std::vector<Tensor_Type>>> Foam::NurbsStructureInterface::computeDistributedLoad
(
    const Field<Tensor_Type> immersedBoundaryField
)
{
    const fvBoundaryMesh& boundary = mesh.boundary();
    const fvPatch& nurbsBoundary = boundary[IBPatchID];
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    const label nurbsBoundaryStart = nurbsBoundary.start();
        
    auto result = std::unique_ptr<std::vector<std::vector<Tensor_Type>>>
    (
        new std::vector<std::vector<Tensor_Type>>(Curves->size())
    );
    for(int nurbsInd=0; nurbsInd<Curves->size(); nurbsInd++)
    {
        (*result)[nurbsInd].resize(nurbsToLimits[nurbsInd].size()-1);
        for(int i=0;i<(*result)[nurbsInd].size();i++)
        {
            std::vector<label>& facesInLimits = nurbsToLimitsFaces[nurbsInd][i];
            std::vector<scalar>& weightsInLimits = nurbsToLimitsFacesWeights[nurbsInd][i];
            
            Foam::zero null;
            Tensor_Type thisSpanValue = Tensor_Type(null);
            //Info<<thisSpanValue<<"  ";
            for(int j=0;j<facesInLimits.size();j++)
            {
                face thisFace = faces[facesInLimits[j]];
                Tensor_Type thisFaceField = immersedBoundaryField[facesInLimits[j]-nurbsBoundaryStart];
                thisSpanValue += weightsInLimits[j]*thisFaceField*thisFace.mag(points);                
            }
            //Info<<nurbsBoundaryStart<<"  "<<facesInLimits.size()<<" "<<weightsInLimits.size()<<thisSpanValue<<Foam::endl;

            Tensor_Type valuePerSpan = thisSpanValue/(nurbsToLimits[nurbsInd][i+1]-nurbsToLimits[nurbsInd][i]);
            (*result)[nurbsInd][i] = valuePerSpan;
        }
    }
    return std::move(result);
}

void Foam::NurbsStructureInterface::assignForceOnCurve()
{
    tmp<GeometricField<Tensor<double>,fvPatchField,volMesh>> gU = fvc::grad(U);
    GeometricField<Tensor<double>,fvPatchField,volMesh>& gradU = gU.ref();
    GeometricField<SymmTensor<double>,fvPatchField,volMesh> totalStress = symm(-p*tensor::one + nu*(gradU + gradU.T()));
    
    const vectorField& Sfp = mesh.Sf().boundaryField()[IBPatchID];
    const scalarField& magSfp = mesh.magSf().boundaryField()[IBPatchID];
    const symmTensorField& totalStressIB = totalStress.boundaryField()[IBPatchID];
    
    Field<vector> ibWallForces = (-Sfp/magSfp) & totalStressIB;
    
    vector avgVec;
    for(vector currVec : ibWallForces)
        avgVec += currVec;
    avgVec /= ibWallForces.size();
    Info<<"ibWallForces sum:"<<avgVec<<Foam::endl;
    
    auto distrLoad = computeDistributedLoad<vector>(ibWallForces);
    
    for(int nurbsInd=0;nurbsInd<distrLoad->size();nurbsInd++)
    {
        const std::vector<vector>& oneCurveDistrLoad = (*distrLoad)[nurbsInd];      
        
        std::vector<double> knotContainer(nurbsToLimits[nurbsInd].size()*2);
        int degree=1;
        double delta = nurbsToLimits[nurbsInd][1]-nurbsToLimits[nurbsInd][0];
        knotContainer[0] = knotContainer[1] = nurbsToLimits[nurbsInd][0] - delta/200;
        for(int k=1;k<nurbsToLimits[nurbsInd].size()-1;k++)
        {
            double delta = nurbsToLimits[nurbsInd][k+1]-nurbsToLimits[nurbsInd][k-1];
            double epsilon = delta/200;
            knotContainer[2*k] = nurbsToLimits[nurbsInd][k]-epsilon;
            knotContainer[2*k+1] = nurbsToLimits[nurbsInd][k]+epsilon;
        }
        knotContainer[knotContainer.size()-2] = knotContainer[knotContainer.size()-1] = nurbsToLimits[nurbsInd].back() + delta/200;
        
        /*
        Info<<"delta:"<<delta<<Foam::endl;
        Info<<"delta/200:"<<delta/200<<Foam::endl;
        Info<<"nurbsToLimits[nurbsInd].back():"<<nurbsToLimits[nurbsInd].back()<<Foam::endl;
        Info<<"nurbsToLimits[nurbsInd].back()+delta/200:"<<nurbsToLimits[nurbsInd].back()+delta/2000<<Foam::endl;
        */
        
        gismo::gsKnotVector<double> knotVector(knotContainer,degree);
        
        gismo::gsMatrix<double> w(nurbsToLimits[nurbsInd].size()*2-2,1);
        for(int i=0;i<nurbsToLimits[nurbsInd].size()*2-2;i++)
            w(i,0) = 1;
        
        gismo::gsMatrix<double> Pcoeff(nurbsToLimits[nurbsInd].size()*2-2,3);
        int k=0;
        for(int i=0;i<nurbsToLimits[nurbsInd].size()*2-2;i+=2,k++)
        {
            for(label d=0;d<3;d++)
            {
                Pcoeff(i,d) = oneCurveDistrLoad[k][d];
                Pcoeff(i+1,d) = oneCurveDistrLoad[k][d];
            }
        }

        //Info<<"RodStart:"<<Rods[nurbsInd]->m_Curve.domainStart()<<Foam::endl;
        //Info<<"RodEnd:"<<Rods[nurbsInd]->m_Curve.domainEnd()<<Foam::endl;
        
        gismo::gsNurbs<double>* forceNurbs = new gismo::gsNurbs<double>(knotVector,w,Pcoeff);
        //Info<<"domainStart:"<<forceNurbs->domainStart()<<Foam::endl;
        //Info<<"domainEnd:"<<forceNurbs->domainEnd()<<Foam::endl;
        gismo::gsKnotVector<double>& knotVectorRead = forceNurbs->knots();
        
        /*
        for(int j=0;j<knotVectorRead.size();j++)
        {
            Info<<knotVectorRead[j]<<"  ";
        }        
        Info<<Foam::endl;
        for(double start=forceNurbs->domainStart(); start<forceNurbs->domainEnd(); start+=0.1)
        {
            gismo::gsMatrix<double> u(1,1);
            u(0,0) = start;
            gismo::gsMatrix<double> f = forceNurbs->eval(u);
            Info<<"f.rows:"<<f.rows()<<Foam::endl;
            Info<<"f.cols:"<<f.cols()<<Foam::endl;
            Info<<u(0,0)<<": ("<<f(0,0)<<","<<f(1,0)<<","<<f(2,0)<<")"<<Foam::endl;
        }
        */
        forceCurveStorage[nurbsInd].reset(forceNurbs);

        (myMesh->m_Rods[nurbsInd])->set_force_lG(forceCurveStorage[nurbsInd].get());
    }
}

void Foam::NurbsStructureInterface::computeIBHeatFlux()
{
    tmp<GeometricField<double,fvsPatchField,surfaceMesh>> gradT = fvc::snGrad(T,IBpatchName);
    GeometricField<double,fvsPatchField,surfaceMesh>& gTF = gradT.ref();
    const GeometricField<double,fvsPatchField,surfaceMesh>::Boundary& gradTBoundary = gTF.boundaryField();
    const fvsPatchField<double> ibgradT = gradTBoundary[IBPatchID];
    
    
    
    gismo::gsNurbs<double> heatFlux;
    
    
    
    Info<<"gradT size:"<<gradTBoundary.size()<<endl; 
    Info<<"ibgradT size:"<<ibgradT.size()<<endl;
}

Foam::word Foam::NurbsStructureInterface::getXMLPath()
{
    Info<<"lateScale:"<<latScale<<endl;
    Info<<"latDir:"<<latDir<<endl;

    fileName caseDirectory = runDirectory+"/"+caseName;
    fileName constantDirectory = caseDirectory;

    DIR  *dir = NULL;
    const char *pathConstantDirectory = constantDirectory.c_str();
    dir = opendir(pathConstantDirectory);
    if(dir==NULL)
        FatalIOError<<"Error reading Nurbs file. Could not open the directory!"<<exit(FatalIOError);
    DynamicList<word> directoryFiles;
    struct dirent *dp;
    while((dp=readdir(dir))!=NULL)
    {
        if(dp->d_name==NULL)
            FatalIOError<<"Reading existing files name as null pointer!"<<exit(FatalIOError);
        directoryFiles.append(word(dp->d_name));
    }
    DynamicList<word> xmlFiles;
    for(int i=0;i<directoryFiles.size();i++)
    {
        std::size_t fileEndingStart = directoryFiles[i].rfind(".");
        if(fileEndingStart==std::string::npos)
            continue;
        word fileEnding = directoryFiles[i].substr(fileEndingStart,directoryFiles[i].size()-fileEndingStart);
        if(fileEnding.compare(".xml")==0)
            xmlFiles.append(directoryFiles[i]);
    }
    if(xmlFiles.size()==0)
        FatalIOError<<"No Nurbs file found!"<<exit(FatalIOError);
    if(xmlFiles.size()>1)
        Info<<"Multiple Nurbs files found. First one will be used!"<<endl;
    word fullPath = constantDirectory+"/"+xmlFiles[0];    
    name = xmlFiles[0].erase(xmlFiles[0].length()-4,std::string::npos);
    return fullPath;
}

int Foam::NurbsStructureInterface::loadRodsFromXML()
{
    printf("loadRodsFromXML\n");
    std::string rodsXMLFilePath = xmlPath;
    bool importSuccess = ActiveRodMesh::import_xmlCrv(rodsList, rodsXMLFilePath, 3, 1, 0);
    if(!importSuccess)
    {
        FatalIOError<<"Importing of Nurbs into rodMesh failed"<<exit(FatalIOError);
    }
	const int  nR = rodsList.size();
    Info<<"rodsList.size():"<<rodsList.size()<<endl;
    Info<<"nR:"<<nR<<endl;
    for(int i=0;i<nR;i++)
    {
        Info<<"---Curve:"<<i<<endl;
        gismo::gsKnotVector<double> knots = rodsList[i].knots();
        Info<<"Knots [";
        for(double knot: knots)
            Info<<knot<<", ";
        Info<<"]"<<endl;
        gismo::gsMatrix<double> coefs = rodsList[i].coefs();
        Info<<"Coefs ["<<endl;
        for(int j=0;j<coefs.rows();j++)
        {
            Info<<"("<<coefs(j,0)<<","<<coefs(j,1)<<","<<coefs(j,2)<<")"<<endl;
        }
        Info<<"]"<<endl;

    }

    double	x0 = 1e6, x1 = -1e6, y0 = 1e6, y1 = -1e6, z0 = 1e6, z1 = -1e6;
    for(int i=0; i<nR; i++)
    {
        for (int j=0; j<rodsList[i].coefs().rows(); j++)
        {
            for (int k=0; k<3; k++)
            {
                rodsList[i].coefs()(j, k) = 1e-6 * round(1e6 * rodsList[i].coefs()(j, k));
            }
        }
        x0 = std::fmin(x0, rodsList[i].coefs().topRows(0).coeff(0, 0));
        x0 = std::fmin(x0, rodsList[i].coefs().bottomRows(1).coeff(0, 0));
        x1 = std::fmax(x1, rodsList[i].coefs().topRows(0).coeff(0, 0));
        x1 = std::fmax(x1, rodsList[i].coefs().bottomRows(1).coeff(0, 0));
        y0 = std::fmin(y0, rodsList[i].coefs().topRows(0).coeff(0, 1));
        y0 = std::fmin(y0, rodsList[i].coefs().bottomRows(1).coeff(0, 1));
        y1 = std::fmax(y1, rodsList[i].coefs().topRows(0).coeff(0, 1));
        y1 = std::fmax(y1, rodsList[i].coefs().bottomRows(1).coeff(0, 1));
        z0 = std::fmin(z0, rodsList[i].coefs().topRows(0).coeff(0, 2));
        z0 = std::fmin(z0, rodsList[i].coefs().bottomRows(1).coeff(0, 2));
        z1 = std::fmax(z1, rodsList[i].coefs().topRows(0).coeff(0, 2));
        z1 = std::fmax(z1, rodsList[i].coefs().bottomRows(1).coeff(0, 2));
    }
    Info<<"x0:"<<x0<<endl;
    Info<<"x1:"<<x1<<endl;
    Info<<"y0:"<<y0<<endl;
    Info<<"y1:"<<y1<<endl;
    Info<<"z0:"<<z0<<endl;
    Info<<"z1:"<<z1<<endl;
    Info<<"lateScale:"<<latScale<<endl;
    Info<<"latDir:"<<latDir<<endl;
    latSize << x1 - x0, y1 - y0, z1 - z0;
    
    latScale=1;
    gsVector<double,3> dX;
    dX << -x0, -y0, -z0;
    for (int i = 0; i < nR; i++)
    {
        rodsList[i].knots().transform(0., 1.);	// re-scale knot vector
        rodsList[i].translate(dX);				// translate to 0
        rodsList[i].scale(latScale);
    }
    Info<<"lateScale:"<<latScale<<endl;
    latSize *= latScale;
    printf("Rods:  %i\n", nR);
    printf("Dimensions: %4.1fx%4.1fx%4.1f mm\n", latSize[0], latSize[1], latSize[2]);
    
    //FatalIOError<<"Temp Stop"<<exit(FatalIOError);
    
    return nR;
}

void Foam::NurbsStructureInterface::createNurbsStructure()
{
    printf("Create Nurbs Structure ... \n");

        // * Base NURBS curve - straight line as B-Spline
    const int	el = 1;
    const int	p = 1;
    const int	n = p + el;
    gsKnotVector<double> knots(0, 1, el - 1, p + 1, 1, p);
    gsMatrix<double> weights(n, 1);
    weights.setOnes();

    // * CS input vectors
    gsVector<double> csE(1);
    csE << geoE0;
    gsVector<double> csR(1);
    csR << geoR0;

    // * CS basis for optimization
    const int ep = 1;
    const int ee = 1;
    const int ekn_n = ep + ee;
    const int ekn_mi = ee - 1;
    gsKnotVector<double> ekn(0., 1., ekn_mi, ep + 1, 1, ep);
    gsMatrix<double> ewts(ekn_n, 1);
    ewts.setOnes();

    // * Basic classes
    gsVector<double, 3> gVec;
    gVec.setZero();
    gVec(latDir) = -9.8;
    gsConstantFunction<double>* loadG;
    loadG = new gsConstantFunction<double>(gVec, 1);
    gsConstantFunction<double>* twist = new gsConstantFunction<double>(0., 1);
    Geo = std::vector< ActiveRodMesh::rodCrossSection* >(nR);
    Rods = std::vector< ActiveRodMesh::rodCosserat* >(nR);
    BasisRef = std::vector< gsNurbsBasis<double>* >(nR);
    
    // * Make rods
    printf("Make rods ... \n");
    #ifdef _OPENMP
    omp_set_num_threads(omp_nthreads);
    #endif	

    #pragma omp parallel for  
    for (int i = 0; i < nR; i++)
    {
        // CS Geometry
        Geo[i] = new ActiveRodMesh::rodCScircle(csE, csR, geoNu, geoRho, Temp0, 0, plot_c);

        // Basis refine
        //rodsList[i]->degreeElevate(p_sim - rodsList[i]->degree());
        //rodsList[i]->uniformRefine(el_sim - rodsList[i]->numElements());
        //Rods[i] = new ActiveRodMesh::rodCosserat(rodsList[i], *Geo[i], twist, 2, 0);
        BasisRef[i] = rodsList[i].basis().clone();
        BasisRef[i]->degreeElevate(p_sim - BasisRef[i]->degree());
        BasisRef[i]->uniformRefine(el_sim - BasisRef[i]->numElements());

        // Make rod
        if (use_mixed)
            Rods[i] = new ActiveRodMesh::rodCosseratMixed(rodsList[i], *BasisRef[i], *Geo[i], twist, 2, 0);
        else
            Rods[i] = new ActiveRodMesh::rodCosserat(rodsList[i], *BasisRef[i], *Geo[i], twist, 2, 0);

        // Cross-section basis for optimization
        Rods[i]->resetEbasis(ekn, ewts);

        // Gravity 
        if (applyGravity)
            Rods[i]->set_force_lG(loadG);
    }
    forceCurveStorage.resize(nR);
}

void Foam::NurbsStructureInterface::createNurbsBoundary()
{
    printf("Boundary conditions ... \n");

    std::vector<bool> uBC0(3, false);
    std::vector<bool> uBCxz(3, false);
    std::vector<bool> uBCyz(3, false);
    std::vector<bool> uBCx(3, true);
    std::vector<bool> uBCy(3, true);
    std::vector<bool> uBCz(3, true);
    std::vector<bool> uBCdir(3, true);
    uBCxz[1] = true;
    uBCyz[0] = true;
    uBCx[0] = false;
    uBCy[1] = false;
    uBCz[2] = false;
    uBCdir[latDir] = false;
    bool rBC0 = false;
    const double bceps = 1e-5;
    const int dir1 = (latDir + 1) % 3;
    const int dir2 = (latDir + 2) % 3;
    gsVector<double, 3> uBCxV;
    uBCxV.setZero();

    pp_rid = std::vector<int>(0);
    std::vector<bool>	pp_eid(0);
    std::vector<double>		pp_wts(0);
    
    for (int i = 0; i < nR; i++)
    {
        Rods[i]->setBCu(0, uBC0); 
        Rods[i]->setBCr(0, rBC0); 
        Rods[i]->setBCu(1, uBC0); 
        Rods[i]->setBCr(1, rBC0); 
    }
    
/*
    if (bcMode == 0)
    {
        for (int i = 0; i < nR; i++)
        {
            if (Rods[i]->m_eval_Crv.leftCols(1)(latDir, 0) < bceps * latSize(latDir))
            {
                Rods[i]->setBCu(0, uBC0); 
                if (bcRot)
                    Rods[i]->setBCr(0, false);
            }
            if (Rods[i]->m_eval_Crv.rightCols(1)(latDir, 0) < bceps * latSize(latDir))
            {
                Rods[i]->setBCu(1, uBC0); 
                if (bcRot)
                    Rods[i]->setBCr(1, false);
            }
        }

        uBCxV(latDir) = loadDir * latSize(latDir) * loadStrain;
        for (int i = 0; i < nR; i++)
        {
            if (Rods[i]->m_eval_Crv.leftCols(1)(latDir, 0) > (1. - bceps) * latSize(latDir))
            {
                Rods[i]->setBCu(0, uBC0, uBCxV); //uBCx
                if (bcRot)
                    Rods[i]->setBCr(0, false);
                pp_rid.push_back(i);
                pp_eid.push_back(0);
                pp_wts.push_back(1.);
            }

            if (Rods[i]->m_eval_Crv.rightCols(1)(latDir, 0) > (1. - bceps) * latSize(latDir))
            {
                Rods[i]->setBCu(1, uBC0, uBCxV); //uBCx
                if (bcRot)
                    Rods[i]->setBCr(1, false);
                pp_rid.push_back(i);
                pp_eid.push_back(1);
                pp_wts.push_back(1.);
            }
        }
    }
    */
}

void Foam::NurbsStructureInterface::setSolverOptions()
{
    printf("Rod mesh ... \n");

    ActiveRodMesh::RodMeshOptions meshOpt;
    meshOpt.name = name;
    myMesh = std::unique_ptr<ActiveRodMesh::rodMesh>(new ActiveRodMesh::rodMesh(Rods, meshOpt));
    myMesh->setTemp(Temp0, Temp0);

    // * Output mesh
    printf("Set ParaView Options ... \n");
    ParaViewOptions vtkOpt;
    vtkOpt.name = name;
    vtkOpt.folder = folder;
    vtkOpt.nPoints = plot_n;
    vtkOpt.nodeThicken = 0;
    vtkOpt.solidGeo = plot_vtk_solid;

    if (plot_vtk_geo)
    {
        ParaViewOptions vtkOptGeo = vtkOpt;
        vtkOptGeo.name = name + "_geo";
        vtkOptGeo.folder = folder + "\\geo";
        vtkOptGeo.unDeformed = 1;
        vtkOptGeo.solidGeo = (plot_vtk > 1);
        /*
        if (_mkdir(vtkOptGeo.folder.c_str()) != 0)
            std::cout << "Problem making folder!\n";
        */
        //myMesh->writeParaView(vtkOptGeo);
        //myMesh.writeXml(vtkOptGeo.name, vtkOptGeo.folder);
    }
    //std::cout << "Total length: " << myMesh.m_initLength << "\n";

    // * Enable contact
    if (isContact)
    {
        printf("Initialize contact ... \n");
        myMesh->contact_init(cntOpt);
    }
// **** End **** //

// **** Setup simulation parameters **** //
// * Parameters
    printf("Setup solve ... \n");

    int solveOK;
    uint64 tg1, tg2;
    myMesh->setNewtonParams(1e-5, 17, 7., 3);
    myMesh->setLinearSolver(1);
    myMesh->setOMPnThreads(omp_nthreads);

    
    solveOpt.stats  = 1;	// Write to command line
    solveOpt.vtkOut = 1;	// Write VTK every X load steps
    solveOpt.mmOut  = 0;	// Write matrices & vectors in 1:first / 2:every Newton step 
    solveOpt.tempCtrl  = 0;	// Temperature controlled simulation (for 4DP)
    solveOpt.vtk_n = plot_n;	// vtk eval. points
    solveOpt.folder = folder;
    solveOpt.name    = "mesh";	// Name for vtk files
    solveOpt.time_method = 1;	// 0: static, 1: Backward Euler, 2: Crank-Nicolson
    solveOpt.time_dt    = 0.01;	// Time step size
    solveOpt.useLoadFun  = 0;	// Use loadFun also for static case

    std::string evalFile = folder + "\\eval";
    myMesh->setSolveEvalOut(std::vector<int>(0), std::vector<double>(0), evalFile);
    
    // * Load steps
    std::vector<double> load_steps;
    int nls = loadSteps;
    std::string loadStr;
    double t1 = 1.;
    double dt = t1 / (1. * nls);

    for (int i=1; i<=nls; i++)
        load_steps.push_back(dt * i);

    // * Pre-part output files
    const int pp_n = pp_rid.size();
    gsVector<double, 3> pp_sum_f;
    gsVector<double, 3> pp_cf, pp_cu;
    std::ofstream outfile1, outfile2, outfile3;
    outfile1.open(folder + "\\output_u.dat");
    outfile2.open(folder + "\\output_f.dat");
    outfile3.open(folder + "\\output_eval.dat");
    for (int i=0; i<pp_n; i++)
    {
        //outfile1 << pp_eid[i] << "\t" << pp_eid[i] << "\t" << pp_eid[i] << "\t";
        outfile1 << "x" << pp_rid[i] << "\ty" << pp_rid[i] << "\tz" << pp_rid[i] << "\t";
        outfile2 << "x" << pp_rid[i] << "\ty" << pp_rid[i] << "\tz" << pp_rid[i] << "\t";
    }
    outfile1 << "\n";
    outfile2 << "\n";
    outfile3 << "t\tlf\tfx\tfy\tfz\tenEl\tenU\tenQ\tenVi\tenHa\tenTot\tdiss\n";
}

Foam::NurbsStructureInterface::~NurbsStructureInterface()
{
    for (int i = 0; i < nR; i++)
	{
		delete Geo[i];
		delete Rods[i];
	}
	Geo.clear();
	Rods.clear();
}

void Foam::NurbsStructureInterface::solveOneStep()
{
    Info<<"Solve Nurbs structure"<<Foam::endl;
    assignForceOnCurve();
    myMesh->solve(1.0,solveOpt);
    moveNurbs();
}

void Foam::NurbsStructureInterface::moveNurbs()
{
    Info<<"Mesh assign Deformation curve"<<Foam::endl;
    label nbrNurbs = myMesh->m_Rods.size();
    std::vector<List<List<vector>>> nurbs_to_new_controlPoints;
    
    for(int i=0;i<nbrNurbs;i++)
    {
        gsMatrix<double> new_CPs = myMesh->m_Rods[i]->m_Def_coefs;
        label nbrControlPoints = new_CPs.cols();
        nurbs_to_new_controlPoints.push_back(List<List<vector>>(1,List<vector>(nbrControlPoints)));
        for(int n=0;n<nbrControlPoints;n++)
        {
            for(int d=0;d<3;d++)
            {
                nurbs_to_new_controlPoints.back()[0][n][d] = new_CPs(d,n);
            }
        }
    }
    mesh.moveNurbsCurves(nurbs_to_new_controlPoints);
    mesh.moveTheMesh();
    mesh.update();
}

void Foam::NurbsStructureInterface::createDeformationCurve()
{
    label nbrNurbs = myMesh->m_Rods.size();
    label nbrNurbs2 = myMesh->m_nR;
    Info<<"nbrNurbs:"<<nbrNurbs<<Foam::endl;
    Info<<"nbrNurbs2:"<<nbrNurbs2<<Foam::endl;
    //FatalIOError<<"Temp Stop"<<exit(FatalIOError);

    std::vector<scalarList> nurbs_to_knots;
    std::vector<List<vector>> nurbs_to_controlPoints;
    std::vector<scalarList> nurbs_to_weights;
    std::vector<label> nurbs_to_degree;

    for(int i=0;i<nbrNurbs;i++)
    {
        gsMatrix<double> P = myMesh->m_Rods[i]->m_Def_coefs;
        gsNurbs<double> curve = myMesh->m_Rods[i]->m_Def;
        gismo::gsKnotVector<double> knots = curve.knots();
        label nbrKnots = knots.size();
        label degree = knots.degree();
        label nbrControlPoints = P.cols();
        if(P.rows()!=3)
            FatalErrorInFunction<<"Wrong Deformation curve dimension"<< exit(FatalError); 

        nurbs_to_controlPoints.push_back(List<vector>(nbrControlPoints));
        nurbs_to_weights.push_back(scalarList(nbrControlPoints));
        for(int n=0;n<nbrControlPoints;n++)
        {
            for(int d=0;d<3;d++)
            {
                nurbs_to_controlPoints.back()[n][d] = P(d,n);
            }
            nurbs_to_weights.back()[n] = 1;
        }
        nurbs_to_degree.push_back(degree);
        nurbs_to_knots.push_back(scalarList(nbrKnots));
        for(int n=0;n<nbrKnots;n++)
        {
            nurbs_to_knots.back()[n] = knots[n];
        }
    }
    mesh.setInitialDeformationCurve(nurbs_to_knots,nurbs_to_controlPoints,nurbs_to_weights,nurbs_to_degree);
}
