#include "Structure.H"

Foam::Structure::Structure
(
    dynamicRefineFvMesh& mesh,
    const Time& runTime
):
runTime(runTime),
runDirectory(runTime.rootPath()),
caseName(runTime.caseName()),
xmlPath(getXMLPath()),
name(getName()),
nR(loadRodsFromXML()),
mesh(mesh),
initialMeshSpacing(initialSpacingFromMesh(mesh))
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
       
    collectMeshHaloData(2);
}

Foam::Structure::~Structure()
{
    for (int i = 0; i < nR; i++)
	{
		delete Geo[i];
		delete Rods[i];
	}
	Geo.clear();
	Rods.clear();
}

Foam::word Foam::Structure::getXMLPath()
{
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
    return fullPath;
}

Foam::word Foam::Structure::getName()
{
    word name = xmlPath;
    label lastSlash = name.rfind("/");
    label lastDot = name.rfind(".xml");
    return name.substr(lastSlash+1,lastDot-lastSlash-1);
}

int Foam::Structure::loadRodsFromXML()
{
    Info<<"loadRodsFromXML"<<Foam::endl;
    //printf("loadRodsFromXML\n");
    std::string rodsXMLFilePath = xmlPath;
    bool importSuccess = ActiveRodMesh::import_xmlCrv(rodsList, rodsXMLFilePath, 3, 1, 0);
    if(!importSuccess)
    {
        FatalIOError<<"Importing of Nurbs into rodMesh failed"<<exit(FatalIOError);
    }
	const int  nR = rodsList.size();
    Info<<"rodsList.size():"<<rodsList.size()<<endl;

    //Info<<"nR:"<<nR<<endl;
    for(int i=0;i<nR;i++)
    {
        Info<<"---Curve:"<<i;
        gismo::gsKnotVector<double> knots = rodsList[i].knots();
        Info<<"  Knots [";
        for(double knot: knots)
            Info<<knot<<", ";
        Info<<"]  ";
        gismo::gsMatrix<double> coefs = rodsList[i].coefs();
        Info<<"  Coefs [";
        for(int j=0;j<coefs.rows();j++)
        {
            Info<<"("<<coefs(j,0)<<","<<coefs(j,1)<<","<<coefs(j,2)<<")";
        }
        Info<<"]";

    }
    Info<<endl;

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
    Info<<"Bounding Box (x:["<<x0<<"-"<<x1<<"], y:["<<y0<<"-"<<y1<<"], z:["<<z0<<"-"<<z1<<"])"<<endl;
    //Info<<"lateScale:"<<latScale<<endl;
    //Info<<"latDir:"<<latDir<<endl;
    latSize << x1 - x0, y1 - y0, z1 - z0;
    
    latScale=1;
    gsVector<double,3> dX;
    dX << -x0, -y0, -z0;
    std::cout<<"dX:"<<dX<<std::endl;
    for (int i = 0; i < nR; i++)
    {
        rodsList[i].knots().transform(0., 1.);	// re-scale knot vector
        //rodsList[i].translate(dX);				// translate to 0
        rodsList[i].scale(latScale);
    }
    //Info<<"lateScale:"<<latScale<<endl;
    latSize *= latScale;
    //printf("Rods:  %i\n", nR);
    //printf("Dimensions: %4.1fx%4.1fx%4.1f mm\n", latSize[0], latSize[1], latSize[2]);
    
    //FatalIOError<<"Temp Stop"<<exit(FatalIOError);
    
    return nR;
}

void Foam::Structure::createNurbsStructure()
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
        
        std::cout<<"rodsList[i].coeffs:"<<rodsList[i].coefs()<<std::endl;
        std::cout<<"Rods[i]->m_Curve.coeffs:"<<Rods[i]->m_Curve.coefs()<<std::endl;
        
        Rods[i]->m_Rot.setCoefs(Rods[i]->m_init_Rot.transpose());
        std::cout<<"Rods[i]->m_Curve.coeffs:"<<Rods[i]->m_Curve.coefs()<<std::endl;

        // Cross-section basis for optimization
        Rods[i]->resetEbasis(ekn, ewts);

        // Gravity 
        if (applyGravity)
            Rods[i]->set_force_lG(loadG);
    }
    
    for(ActiveRodMesh::rodCosserat* rodPtr : Rods)
        std::cout<<rodPtr->m_Curve.coefs()<<std::endl;
    
    printf("Rod mesh ... \n");

    ActiveRodMesh::RodMeshOptions meshOpt;
    meshOpt.name = name;
    myMesh = std::unique_ptr<ActiveRodMesh::rodMesh>(new ActiveRodMesh::rodMesh(Rods, meshOpt));
    myMesh->setTemp(Temp0, Temp0);
}

void Foam::Structure::createNurbsBoundary()
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

void Foam::Structure::setSolverOptions()
{
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

void Foam::Structure::updateRodCoordinateSystem()
{
    for(uint rodIndex=0; rodIndex<myMesh->m_Rods.size(); rodIndex++)
    {
        myMesh->m_Rods[rodIndex]->m_Rot.setCoefs(myMesh->m_Rods[rodIndex]->m_Rot_coefs.transpose());
    }
}

label Foam::Structure::getNumberRods()
{
    return myMesh->m_Rods.size();
}

label Foam::Structure::getMaxDegree
(
    const ActiveRodMesh::rodCosserat* oneRod
)
{
    const gsNurbs<double>& base = oneRod->m_Curve;
    label baseDegreeX = base.degree();
    //label baseDegreeY = base.degree(1);
    //label baseDegreeZ = base.degree(2);
    const gsNurbs<double>& def = oneRod->m_Def;
    label defDegreeX = def.degree();
    //label defDegreeY = def.degree(1);
    //label defDegreeZ = def.degree(2);
    /*
    label maxDegree = std::max<label>(baseDegreeX,baseDegreeY);
    maxDegree = std::max(maxDegree,baseDegreeY);
    maxDegree = std::max(maxDegree,baseDegreeZ);
    maxDegree = std::max(maxDegree,defDegreeX);
    maxDegree = std::max(maxDegree,defDegreeY);
    maxDegree = std::max(maxDegree,defDegreeZ);
    */
    return std::max<label>(baseDegreeX,defDegreeX);
}

vector Foam::Structure::rodDerivEval
(
    const ActiveRodMesh::rodCosserat* rod,
    scalar parameter
)
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    gsMatrix<scalar> baseDeriv;
    const gsNurbs<scalar>& curve = rod->m_Curve;
    curve.deriv_into(parMat,baseDeriv);

    gsMatrix<scalar> defDeriv;
    const gsNurbs<scalar>& def = rod->m_Def;
    def.deriv_into(parMat,defDeriv);
    
    gsMatrix<scalar> pnt = baseDeriv+defDeriv;    
    return vector(pnt(0,0),pnt(1,0),pnt(2,0));    
}

void Foam::Structure::rodEval
(
    const ActiveRodMesh::rodCosserat* rod,
    scalar parameter,
    vector& d1,
    vector& d2,
    vector& d3,
    vector& r
)
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
        
    gsMatrix<scalar> basePnt;
    const gsNurbs<scalar>& curve = rod->m_Curve;
    curve.eval_into(parMat,basePnt);

    gsMatrix<scalar> defPnt;
    const gsNurbs<scalar>& def = rod->m_Def;
    def.eval_into(parMat,defPnt);
    
    gsMatrix<scalar> pnt = basePnt+defPnt;    
    r = vector(pnt(0,0),pnt(1,0),pnt(2,0));

    gsMatrix<scalar> rotQuat;
    const gsNurbs<scalar>& quat = rod->m_Rot;
    quat.eval_into(parMat,rotQuat);
    gsMatrix<scalar,3,3> R;
    ActiveRodMesh::quat_R(rotQuat,R);
    
    d1 = vector(R(0,0),R(1,0),R(2,0));
    d2 = vector(R(0,1),R(1,1),R(2,1));
    d3 = vector(R(0,2),R(1,2),R(2,2));
}

void Foam::Structure::rodEval
(
    const ActiveRodMesh::rodCosserat* rod,
    scalar parameter,
    vector& r
)
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
        
    gsMatrix<scalar> basePnt;
    const gsNurbs<scalar>& curve = rod->m_Curve;
    curve.eval_into(parMat,basePnt);

    gsMatrix<scalar> defPnt;
    const gsNurbs<scalar>& def = rod->m_Def;
    def.eval_into(parMat,defPnt);
    
    gsMatrix<scalar> pnt = basePnt+defPnt;    
    r = vector(pnt(0,0),pnt(1,0),pnt(2,0));
}

std::pair<double,double> Foam::Structure::minMaxSpan
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

void Foam::Structure::cellDistances
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

scalar Foam::Structure::supportDomainMinSize
(
    const DynamicList<label>& supportDomainCells
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

void Foam::Structure::collectMeshHaloData
(
    label iterations
)
{
    if(iterations<1)
        FatalErrorInFunction<<"Must be at least one iterations"<<exit(FatalError);
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const pointField& points = mesh.points();
    const polyBoundaryMesh& boundaries = mesh.boundaryMesh();
    
    // DynamicList<std::tuple<neighborProcess,start of faces,size of border>>
    DynamicList<std::tuple<label,label,label>> neighborProcesses;
    
    // Collect halo data for own process along the processor boundary patches
    label nbrPatchFaces = 0;
    List<DynamicList<label>> procPatchOwner(Pstream::nProcs());
    List<DynamicList<label>> procPatchNeighbour(Pstream::nProcs());
    List<DynamicList<label>> procPatchIndex(Pstream::nProcs());
    List<DynamicList<label>> procPatchFaceLocalIndex(Pstream::nProcs());
    List<DynamicList<List<DynamicList<label>>>> procCellInds(Pstream::nProcs());
    List<DynamicList<List<DynamicList<vector>>>> procCellCentres(Pstream::nProcs());
    List<DynamicList<List<DynamicList<scalar>>>> procCellMags(Pstream::nProcs());
    
    for(label patchIndex=0; patchIndex<boundaries.size(); patchIndex++)
    {
        const polyPatch& patch = boundaries[patchIndex];
        if(isA<processorPolyPatch>(patch))
        {
            const processorPolyPatch* pPP = dynamic_cast<const processorPolyPatch*>(&patch);
            label neighborProcess = pPP->neighbProcNo();
            label patchStartFace = pPP->start();
            label patchSizeFaces = pPP->size();

            neighborProcesses.append({neighborProcess,patchStartFace,patchSizeFaces});

            for(label locFaceInd=0; locFaceInd<patchSizeFaces; locFaceInd++)
            {
                label faceInd = locFaceInd+patchStartFace;
                if(faceInd<neighbour.size())
                    FatalErrorInFunction<<"Boundary face has a neighbour"<<exit(FatalError);
                
                procPatchOwner[Pstream::myProcNo()].append(Pstream::myProcNo());
                procPatchNeighbour[Pstream::myProcNo()].append(neighborProcess);
                procPatchIndex[Pstream::myProcNo()].append(patchIndex);
                procPatchFaceLocalIndex[Pstream::myProcNo()].append(locFaceInd);
                
                procCellInds[Pstream::myProcNo()].append(List<DynamicList<label>>());                procCellCentres[Pstream::myProcNo()].append(List<DynamicList<vector>>());
                procCellMags[Pstream::myProcNo()].append(List<DynamicList<scalar>>());
                
                List<DynamicList<label>>& thisFaceCellInds = procCellInds[Pstream::myProcNo()].last();
                thisFaceCellInds.resize(iterations);
                List<DynamicList<vector>>& thiFaceCellCentres = procCellCentres[Pstream::myProcNo()].last();
                thiFaceCellCentres.resize(iterations);
                List<DynamicList<scalar>>& thiFaceCellMags = procCellMags[Pstream::myProcNo()].last();
                thiFaceCellMags.resize(iterations);

                std::unordered_set<label> borderCells;
                label cellInd = owner[faceInd];
                
                auto addCellData = [&](label cellInd,label k)
                {
                    thisFaceCellInds[k].append(cellInd);
                    borderCells.insert(cellInd);
                    thiFaceCellCentres[k].append(cells[cellInd].centre(points,faces));
                    thiFaceCellMags[k].append(cells[cellInd].mag(points,faces));
                };
                
                addCellData(cellInd,0);
                for(label k=1; k<iterations; k++)
                {
                    for(label frontCellInd : thisFaceCellInds[k-1])
                    {
                        const cell& frontCell = cells[frontCellInd];
                        for(label faceInd : frontCell)
                        {
                            if(owner[faceInd]!=frontCellInd)
                            {
                                if(borderCells.find(owner[faceInd])==borderCells.end())
                                {
                                    addCellData(owner[faceInd],k);
                                }
                            }
                            if(faceInd<neighbour.size())
                            {
                                if(neighbour[faceInd]!=frontCellInd)
                                {
                                    if(borderCells.find(neighbour[faceInd])==borderCells.end())
                                    {
                                        addCellData(neighbour[faceInd],k);
                                    }
                                }
                            }
                        }
                    }
                }
                
                nbrPatchFaces++;
            }
        }
    }
        
    //Exchange data
    Pstream::gatherList(procPatchOwner);
    Pstream::gatherList(procPatchNeighbour);
    Pstream::gatherList(procPatchIndex);
    Pstream::gatherList(procPatchFaceLocalIndex);
    Pstream::gatherList(procCellInds);
    Pstream::gatherList(procCellCentres);
    Pstream::gatherList(procCellMags);
    
    Pstream::scatterList(procPatchOwner);
    Pstream::scatterList(procPatchNeighbour);
    Pstream::scatterList(procPatchIndex);
    Pstream::scatterList(procPatchFaceLocalIndex);
    Pstream::scatterList(procCellInds);
    Pstream::scatterList(procCellCentres);
    Pstream::scatterList(procCellMags);
    
    /*
    if(Pstream::myProcNo()==0)
    {
        Pout<<"procPatchOwner"<<Foam::endl;
        Pout<<procPatchOwner.size()<<Foam::endl;
        for(label proc=0; proc<procPatchOwner.size(); proc++)
            for(label cellInd=0; cellInd<procPatchOwner[proc].size(); cellInd++)
                Pout<<proc<<":"<<cellInd<<":"<<procPatchOwner[proc][cellInd]<<" | "<<procPatchNeighbour[proc][cellInd]<<" | "<<procPatchIndex[proc][cellInd]<<" | "<<procPatchFaceLocalIndex[proc][cellInd]<<" | "<<procCellInds[proc][cellInd]<<procCellCentres[proc][cellInd]<<Foam::endl;
    }
    */
    
    globalHaloCellList_Sorted.resize(Pstream::nProcs());
    globalHaloCellToIndexMap.resize(Pstream::nProcs());
    for(label process=0; process<Pstream::nProcs(); process++)
    {
        const DynamicList<label>& patchOwner = procPatchOwner[process];
        const DynamicList<label>& patchNeighbour = procPatchNeighbour[process];
        const DynamicList<label>& patchIndex = procPatchIndex[process];
        const DynamicList<label>& patchFaceLocalIndex = procPatchFaceLocalIndex[process];
        const DynamicList<List<DynamicList<label>>>& patchCellInds = procCellInds[process];
        const DynamicList<List<DynamicList<vector>>>& patchCellCentres = procCellCentres[process];
        const DynamicList<List<DynamicList<scalar>>>& patchCellMags = procCellMags[process];
        
        std::unordered_map<label,std::tuple<label,label,label>> uniqueHaloCells;
        for(label patchFaceInd=0; patchFaceInd<patchOwner.size(); patchFaceInd++)
        {
            label owner = patchOwner[patchFaceInd];
            label neighbour = patchNeighbour[patchFaceInd];
            label index = patchIndex[patchFaceInd];
            label faceLocalIndex = patchFaceLocalIndex[patchFaceInd];
            const List<DynamicList<label>>& cellInds = patchCellInds[patchFaceInd];
            //const List<DynamicList<vector>>& cellCentres = patchCellCentres[patchFaceInd];
            //const List<DynamicList<scalar>>& cellMags = patchCellMags[patchFaceInd];
            
            if(owner!=process)
                FatalErrorInFunction<<"Invalid data"<<exit(FatalError);
            
            for(label iteration=0; iteration<cellInds.size(); iteration++)
            {
                const DynamicList<label>& iterCellInds = cellInds[iteration];
                //const DynamicList<vector>& iterCellCentres = cellCentres[iteration];
                //const DynamicList<scalar>& iterCellMags = cellMags[iteration];
                
                for(label i=0; i<iterCellInds.size(); i++)
                {
                    label cellInd = iterCellInds[i];
                    if(uniqueHaloCells.find(cellInd)==uniqueHaloCells.end())
                    {
                        uniqueHaloCells[cellInd] = {patchFaceInd,iteration,i};
                    }
                }
            }
        }
        
        DynamicList<CellDescription>& haloCellList_Sorted = globalHaloCellList_Sorted[process];
        for(auto iter=uniqueHaloCells.begin(); iter!=uniqueHaloCells.end(); iter++)
        {
            label cellInd = iter->first;
            label patchFaceInd = std::get<0>(iter->second);
            label iterNbr = std::get<1>(iter->second);
            label i = std::get<2>(iter->second);
            if(cellInd!=patchCellInds[patchFaceInd][iterNbr][i])
                FatalErrorInFunction<<"Invalid data"<<exit(FatalError);
            CellDescription oneCellData =
            {
                patchCellInds[patchFaceInd][iterNbr][i],
                patchCellCentres[patchFaceInd][iterNbr][i],
                patchCellMags[patchFaceInd][iterNbr][i]
            };
            haloCellList_Sorted.append(oneCellData);
        }
        std::sort(haloCellList_Sorted.begin(),haloCellList_Sorted.end(),
                  [](const CellDescription& a, const CellDescription& b)
                    {
                        return a.index < b.index;
                    }
                );
        for(label index=0; index<haloCellList_Sorted.size(); index++)
        {
            globalHaloCellToIndexMap[process][haloCellList_Sorted[index].index] = index;
        }
    }
    
    // [proc] -> neighbour -> patchInd -> locFaceInd -> index
    List<std::unordered_map<label,std::unordered_map<label,std::map<label,label>>>> globalProcessToPatchIndToLocInd(Pstream::nProcs());
    for(label process=0; process<Pstream::nProcs(); process++)
    {
        const DynamicList<label>& patchOwner = procPatchOwner[process];
        const DynamicList<label>& patchNeighbour = procPatchNeighbour[process];
        const DynamicList<label>& patchIndex = procPatchIndex[process];
        const DynamicList<label>& patchFaceLocalIndex = procPatchFaceLocalIndex[process];
        const DynamicList<List<DynamicList<label>>>& patchCellInds = procCellInds[process];
        
        for(label patchFaceInd=0; patchFaceInd<patchOwner.size(); patchFaceInd++)
        {
            label owner = patchOwner[patchFaceInd];
            label neighbour = patchNeighbour[patchFaceInd];
            label index = patchIndex[patchFaceInd];
            label faceLocalIndex = patchFaceLocalIndex[patchFaceInd];
            const List<DynamicList<label>>& cellInds = patchCellInds[patchFaceInd];
            
            if(owner!=process)
                FatalErrorInFunction<<"Data mismatch"<<exit(FatalError);
            
            globalProcessToPatchIndToLocInd[process][neighbour][index][faceLocalIndex] = patchFaceInd;
        }
    }    
    
    for(label patchIndex=0; patchIndex<boundaries.size(); patchIndex++)
    {
        const polyPatch& patch = boundaries[patchIndex];
        if(isA<processorPolyPatch>(patch))
        {
            const processorPolyPatch* pPP = dynamic_cast<const processorPolyPatch*>(&patch);
            label neighborProcess = pPP->neighbProcNo();
            label patchStartFace = pPP->start();
            label patchSizeFaces = pPP->size();
            
            auto iter = globalProcessToPatchIndToLocInd[neighborProcess].find(Pstream::myProcNo());
            if(iter==globalProcessToPatchIndToLocInd[neighborProcess].end())
                FatalErrorInFunction<<"Data mismatch"<<exit(FatalError);
                
            // patchInd -> locFaceInd -> index
            const std::unordered_map<label,std::map<label,label>>& patchIndToLocInd = iter->second;
            if(patchIndToLocInd.size()!=1)
                FatalErrorInFunction<<"Neighbour must have only one patch type"<<exit(FatalError);
            const std::map<label,label>& locFaceIndToListInd = patchIndToLocInd.cbegin()->second;
            
            const DynamicList<label>& patchOwner = procPatchOwner[neighborProcess];
            const DynamicList<label>& patchNeighbour = procPatchNeighbour[neighborProcess];
            const DynamicList<label>& patchIndex = procPatchIndex[neighborProcess];
            const DynamicList<label>& patchFaceLocalIndex = procPatchFaceLocalIndex[neighborProcess];
            const DynamicList<List<DynamicList<label>>>& patchCellInds = procCellInds[neighborProcess];

            for(label locFacePatchInd=0; locFacePatchInd<patchSizeFaces; locFacePatchInd++)
            {
                auto iter = locFaceIndToListInd.find(locFacePatchInd);
                if(iter==locFaceIndToListInd.end())
                    FatalErrorInFunction<<"Missing local face indice in patch"<<exit(FatalError);
                label patchFaceInd = iter->second;             

                label owner = patchOwner[patchFaceInd];
                label neighbour = patchNeighbour[patchFaceInd];
                label index = patchIndex[patchFaceInd];
                label faceLocalIndex = patchFaceLocalIndex[patchFaceInd];
                const List<DynamicList<label>>& cellInds = patchCellInds[patchFaceInd];
                
                if(owner!=neighborProcess)
                    FatalErrorInFunction<<"Data mismatch"<<exit(FatalError);
                if(neighbour!=Pstream::myProcNo())
                    FatalErrorInFunction<<"Data mismatch"<<exit(FatalError);
                if(faceLocalIndex!=locFacePatchInd)
                    FatalErrorInFunction<<"Data mismatch"<<exit(FatalError);

                label faceInd = patchStartFace+locFacePatchInd;
                for(const DynamicList<label>& iterCell : cellInds)
                {
                    for(label cellInd : iterCell)
                    {
                        patchFaceToCellMap[faceInd].append({neighborProcess,cellInd});
                    }
                }
            }
        }
    }
    
    /*
    if(Pstream::myProcNo()==0)
    {
        Pout<<"globalHaloCellToIndexMap"<<Foam::endl;
        Pout<<globalHaloCellToIndexMap.size()<<Foam::endl;
        for(label proc=0; proc<globalHaloCellToIndexMap.size(); proc++)
            for(auto iter=globalHaloCellToIndexMap[proc].begin(); iter!=globalHaloCellToIndexMap[proc].end(); iter++)
                Pout<<proc<<" : "<<iter->first<<"->"<<iter->second<<Foam::endl;
            
        Pout<<"globalHaloCellList_Sorted"<<Foam::endl;
        Pout<<globalHaloCellList_Sorted.size()<<Foam::endl;
        for(label proc=0; proc<globalHaloCellList_Sorted.size(); proc++)
            for(label cellInd=0; cellInd<globalHaloCellList_Sorted[proc].size(); cellInd++)
                Pout<<proc<<":"<<cellInd<<":"<<globalHaloCellList_Sorted[proc][cellInd].index<<Foam::endl;
        
        Pout<<"patchFaceToCellMap"<<Foam::endl;
        Pout<<patchFaceToCellMap.size()<<Foam::endl;
        for(auto iter=patchFaceToCellMap.begin(); iter!=patchFaceToCellMap.end(); iter++)
        {
            Pout<<"faceInd:"<<iter->first<<"->"<<Foam::endl;
            for(auto pair : iter->second)
                Pout<<"\t\t"<<"{"<<pair.first<<","<<pair.second<<"}"<<Foam::endl;
        }
    }
    */
}

scalar Foam::Structure::initialSpacingFromMesh
(
    const dynamicRefineFvMesh& mesh
)
{
    Info<<"initialSpacingFromMesh cellNbr:"<<mesh.cells().size()<<Foam::endl;
    if(mesh.cells().size()<1)
        FatalErrorInFunction<<"Mesh not valid"<<exit(FatalError);
    const cell& oneCell = mesh.cells()[0];
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    const pointField& oneCellPoints = oneCell.points(faces,points);
    vector oneCellCentre = oneCell.centre(points,faces);
    scalar minHalfDist = std::numeric_limits<scalar>::max();
    for(vector pnt : oneCellPoints)
    {
        vector conn = pnt-oneCellCentre;
        for(label dim=0; dim<3; dim++)
        {
            scalar connD = std::abs(conn[dim]);
            if(connD<minHalfDist)
                minHalfDist=connD;
        }
    }
    return 2*minHalfDist;
}
