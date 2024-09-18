#include "Structure.H"

Foam::Structure::Structure
(
    const fvMesh& mesh,
    const Time& runTime
):
initialMeshSpacing(initialSpacingFromMesh(mesh)),
runTime(runTime),
runDirectory(runTime.rootPath()),
caseName(runTime.caseName()),
xmlPath(getXMLPath()),
name(getName()),
nR(loadRodsFromXML()),
mesh(mesh),
meshBoundingBox(computeMeshBoundingBox())
{
    collectMeshHaloData(4);
    setupActiveRodMesh();
    setupRodBoundingBoxTree();
}

Foam::Structure::Structure
(
    const fvMesh& mesh,
    const IOdictionary& stuctureDict,
    const Time& runTime
):
initialMeshSpacing(initialSpacingFromMesh(mesh)),
runTime(runTime),
runDirectory(runTime.rootPath()),
caseName(runTime.caseName()),
xmlPath(xmlFromDict(stuctureDict)),
name(getName()),
nR(loadRodsFromXML()),
mesh(mesh),
meshBoundingBox(computeMeshBoundingBox())
{
    collectMeshHaloData(4);
    setupActiveRodMesh();
    setupRodBoundingBoxTree();
}

Foam::Structure::~Structure()
{
    /*
    for (int i = 0; i < nR; i++)
	{
		//delete Geo[i];
		delete Rods[i];
	}
	Geo.clear();
	Rods.clear();
	*/
    cleanupActiveRodMesh();
}

void Foam::Structure::cleanupActiveRodMesh()
{
	for(std::size_t i=0; i<Rods.size(); i++)
        delete Rods[i];
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
        /*
        if(dp->d_name==NULL)
            FatalIOError<<"Reading existing files name as null pointer!"<<exit(FatalIOError);
        */
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

Foam::word Foam::Structure::xmlFromDict(const IOdictionary& stuctureDict)
{
    ITstream rodTypeStream = stuctureDict.lookup("rodFile");
    token rodFileToken;
    rodTypeStream.read(rodFileToken);
    if(!rodFileToken.isString())
        FatalErrorInFunction<<"Invalid entry in structure/structureDict/rodFile -- must be string"<<exit(FatalError);
    word rodFileWord = rodFileToken.stringToken();
    fileName caseDirectory = runDirectory+"/"+caseName+"/structure/";
    word fullPath = caseDirectory+rodFileWord;
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
    //Info<<"loadRodsFromXML"<<Foam::endl;
    //printf("loadRodsFromXML\n");
    std::string rodsXMLFilePath = xmlPath;
    bool importSuccess = ActiveRodMesh::import_xmlCrv(rodsList, rodsXMLFilePath, 3, 1, 0);
    if(!importSuccess)
    {
        FatalIOError<<"Importing of Nurbs into rodMesh failed"<<exit(FatalIOError);
    }
	const int  nR = rodsList.size();
    return nR;
}

void Foam::Structure::createRodScaling()
{
    //Info<<"rodsList.size():"<<rodsList.size()<<endl;

    //Info<<"nR:"<<nR<<endl;
    /*
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
    */

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
    //Info<<"Bounding Box (x:["<<x0<<"-"<<x1<<"], y:["<<y0<<"-"<<y1<<"], z:["<<z0<<"-"<<z1<<"])"<<endl;
    //Info<<"lateScale:"<<latScale<<endl;
    //Info<<"latDir:"<<latDir<<endl;
    latSize << x1 - x0, y1 - y0, z1 - z0;
    
    latScale=1;
    gsVector<double,3> dX;
    dX << -x0, -y0, -z0;
    //std::cout<<"dX:"<<dX<<std::endl;
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

    //#pragma omp parallel for  
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
        
        //std::cout<<"rodsList[i].coeffs:"<<rodsList[i].coefs()<<std::endl;
        //std::cout<<"Rods[i]->m_Curve.coeffs:"<<Rods[i]->m_Curve.coefs()<<std::endl;
        
        Rods[i]->m_Rot.setCoefs(Rods[i]->m_init_Rot.transpose());
        //std::cout<<"Rods[i]->m_Curve.coeffs:"<<Rods[i]->m_Curve.coefs()<<std::endl;

        // Cross-section basis for optimization
        Rods[i]->resetEbasis(ekn, ewts);

        // Gravity 
        if (applyGravity)
            Rods[i]->set_force_lG(loadG);
    }
    
    /*
    for(ActiveRodMesh::rodCosserat* rodPtr : Rods)
        std::cout<<rodPtr->m_Curve.coefs()<<std::endl;
    */
    
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
    //const double bceps = 1e-5;
    //const int dir1 = (latDir + 1) % 3;
    //const int dir2 = (latDir + 2) % 3;
    //gsVector<double, 3> uBCxV;
    //uBCxV.setZero();

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

    //int solveOK;
    //uint64 tg1, tg2;
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
    
    Info<<"setSolverOptions done"<<Foam::endl;
}

void Foam::Structure::checkActiveRodMesh()
{
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        if(rodsList[rodNumber].coefs()!=Rods[rodNumber]->m_Curve.coefs())
            FatalErrorInFunction<<"Mismatch in basis coefficients!"<<exit(FatalError);
    }
}

void Foam::Structure::setupActiveRodMesh()
{
    cleanupActiveRodMesh();
    
    cntOpt.ptsType = 2;
    cntOpt.ptsN = 2;
    cntOpt.csFac = 2.0;
    cntOpt.preg = 0.02 * geoR0;
    cntOpt.kc = 0.2 * geoE0*geoR0*geoR0;
    cntOpt.initOut = 0;

    folder = runDirectory+"/"+caseName;

    createRodScaling();
    createNurbsStructure();
    createNurbsBoundary();
    setSolverOptions();
    checkActiveRodMesh();
    updateRodCoordinateSystem();

    for(label rodI=0; rodI<nR; rodI++)
    {
        const ActiveRodMesh::rodCosserat* rod = Rods[rodI];
        const gsNurbs<scalar>& curve = rod->m_Curve;
        const gsNurbs<scalar>& deformation = rod->m_Def;
        if(curve.domainStart()!=deformation.domainStart())
            FatalErrorInFunction<<"Mismatch in curve and deformation start"<<exit(FatalError);
        if(curve.domainEnd()!=deformation.domainEnd())
            FatalErrorInFunction<<"Mismatch in curve and deformation end"<<exit(FatalError);
    }
    
    //constructCoeffDerivedData();
    
    prevDeformations.clear();
    prevDeformations.resize(nR);
    prevRotations.clear();
    prevRotations.resize(nR);
}

void Foam::Structure::updateRodCoordinateSystem()
{
    for(uint rodIndex=0; rodIndex<myMesh->m_Rods.size(); rodIndex++)
    {
        myMesh->m_Rods[rodIndex]->m_Rot.setCoefs(myMesh->m_Rods[rodIndex]->m_Rot_coefs.transpose());
    }
}

Foam::label Foam::Structure::getNumberRods() const
{
    return nR;
}

Foam::label Foam::Structure::getMaxDegree
(
    const ActiveRodMesh::rodCosserat* oneRod
) const
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

gsNurbs<Foam::scalar> Foam::Structure::createNurbs
(
    std::vector<scalar> knots,
    uint degree,
    std::vector<scalar> weights,
    std::vector<scalar> coefficients
)
{
    gsKnotVector<scalar> cKnots(knots,degree);
    Info<<"Generated gsKnotVector"<<Foam::endl;
    gsMatrix<scalar> cWeight(weights.size(),1);
    for(uint i=0; i<weights.size(); i++)
        cWeight.at(i) = weights[i];
    Info<<"Generated weights"<<Foam::endl;
    gsMatrix<scalar> cCoeff(coefficients.size(),1);
    for(uint i=0; i<coefficients.size(); i++)
        cCoeff.at(i) = coefficients[i];
    Info<<"Generated coefficients"<<Foam::endl;
    return gsNurbs<scalar>(cKnots,cWeight,cCoeff);
}

gsNurbs<Foam::scalar> Foam::Structure::createNurbs
(
    std::vector<scalar> knots,
    uint degree,
    std::vector<scalar> weights,
    std::vector<vector> coefficients
)
{
    gsKnotVector<scalar> cKnots(knots,degree);
    gsMatrix<scalar> cWeight(weights.size(),1);
    for(uint i=0; i<weights.size(); i++)
        cWeight.at(i) = weights[i];
    gsMatrix<scalar> cCoeff(coefficients.size(),3);
    for(uint i=0; i<coefficients.size(); i++)
    {
        for(label d=0;d<3;d++)
        {
            cCoeff(i,d) = coefficients[i][d];
        }
    }            
    return gsNurbs<scalar>(cKnots,cWeight,cCoeff);
}

const std::pair<gsNurbs<Foam::scalar>,Foam::scalar>* Foam::Structure::readPrevRodDeformation
(
    label rodNumber
) const
{
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber given"<<exit(FatalError);
    std::pair<gsNurbs<scalar>,scalar>* prevDef = nullptr;
    if(prevDeformations[rodNumber])
        prevDef = prevDeformations[rodNumber].get();
    return prevDef;
}

const std::pair<gsNurbs<Foam::scalar>,Foam::scalar>* Foam::Structure::readPrevRodRotation
(
    label rodNumber
) const
{
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber given"<<exit(FatalError);
    std::pair<gsNurbs<scalar>,scalar>* prevRot = nullptr;
    if(prevRotations[rodNumber])
        prevRot = prevRotations[rodNumber].get();
    return prevRot;
}

void Foam::Structure::setDeformation
(
    const List<List<vector>>& deformationCoeffs
)
{
    if(deformationCoeffs.size()!=nR)
        FatalErrorInFunction<<"Mismatch in number of deformationCoeffs and number of rods!"<<exit(FatalError);
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        setDeformation(rodNumber,deformationCoeffs[rodNumber]);
    }
    updateRodCoordinateSystem();
    setupRodBoundingBoxTree();
}

void Foam::Structure::setDeformation
(
    int rodNumber,
    const List<vector>& deformationCoeffs
)
{
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber given"<<exit(FatalError);
    gsNurbs<scalar>& def = Rods[rodNumber]->m_Def;
    gsMatrix<scalar>& defCoeffs = def.coefs();
    if(defCoeffs.rows()!=deformationCoeffs.size())
        FatalErrorInFunction<<"Mismatch in deformation coefficient number!"<<exit(FatalError);
    if(defCoeffs.cols()!=3)
        FatalErrorInFunction<<"Invalid coefficient dimension!"<<exit(FatalError);
    
    if(!prevDeformations[rodNumber])
        prevDeformations[rodNumber] = std::make_unique<std::pair<gsNurbs<scalar>,scalar>>();
    prevDeformations[rodNumber]->first = Rods[rodNumber]->m_Def;
    prevDeformations[rodNumber]->second = mesh.time().value();
    
    if(!prevRotations[rodNumber])
        prevRotations[rodNumber] = std::make_unique<std::pair<gsNurbs<scalar>,scalar>>();
    prevRotations[rodNumber]->first = Rods[rodNumber]->m_Rot;
    prevRotations[rodNumber]->second = mesh.time().value();
    
    for(label coeffNum=0; coeffNum<deformationCoeffs.size(); coeffNum++)
    {
        defCoeffs(coeffNum,0) = deformationCoeffs[coeffNum][0];
        defCoeffs(coeffNum,1) = deformationCoeffs[coeffNum][1];
        defCoeffs(coeffNum,2) = deformationCoeffs[coeffNum][2];        
    }
}

void Foam::Structure::store()
{
    std::vector<gsNurbs<scalar>>& timeDefs = storage[mesh.time().value()];
    for(const ActiveRodMesh::rodCosserat* rod : Rods)
    {
        timeDefs.push_back(rod->m_Def);
    }
}

void Foam::Structure::setToTime(scalar time)
{
    if(storage.find(time)==storage.end())
        FatalErrorInFunction<<"No data exists for "<<time<<exit(FatalError);
    
    std::vector<gsNurbs<scalar>>& timeDefs = storage[time];
    for(std::size_t i=0; i<Rods.size(); i++)
    {
        Rods[i]->m_Def = timeDefs[i];
    }
}

Foam::vector Foam::Structure::rodEval
(
    const ActiveRodMesh::rodCosserat* rod,
    scalar parameter
)
{
    vector r;
    rodEval(rod,parameter,r);
    return r;
}

Foam::vector Foam::Structure::rodDerivEval
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

Foam::vector Foam::Structure::rodDeriv2Eval
(
    const ActiveRodMesh::rodCosserat* rod,
    scalar parameter
)
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    gsMatrix<scalar> baseDeriv2;
    const gsNurbs<scalar>& curve = rod->m_Curve;
    curve.deriv2_into(parMat,baseDeriv2);

    gsMatrix<scalar> defDeriv2;
    const gsNurbs<scalar>& def = rod->m_Def;
    def.deriv2_into(parMat,defDeriv2);
    
    gsMatrix<scalar> pnt = baseDeriv2+defDeriv2;    
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
    rodEval(rod->m_Curve,rod->m_Def,rod->m_Rot,parameter,d1,d2,d3,r);
}

void Foam::Structure::rodEval
(
    const ActiveRodMesh::rodCosserat* rod,
    scalar parameter,
    vector& r
)
{    
    rodEval(rod->m_Curve,rod->m_Def,parameter,r);
}

void Foam::Structure::rodEval
(
    const gsNurbs<scalar>& curve,
    const gsNurbs<scalar>& def,
    const gsNurbs<scalar>& rot,
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
    curve.eval_into(parMat,basePnt);

    gsMatrix<scalar> defPnt;
    def.eval_into(parMat,defPnt);
    
    gsMatrix<scalar> pnt = basePnt+defPnt;    
    r = vector(pnt(0,0),pnt(1,0),pnt(2,0));

    gsMatrix<scalar> rotQuat;
    rot.eval_into(parMat,rotQuat);
    gsMatrix<scalar,3,3> R;
    ActiveRodMesh::quat_R(rotQuat,R);
    
    d1 = vector(R(0,0),R(1,0),R(2,0));
    d2 = vector(R(0,1),R(1,1),R(2,1));
    d3 = vector(R(0,2),R(1,2),R(2,2));
}

void Foam::Structure::rodEval
(
    const gsNurbs<scalar>& curve,
    const gsNurbs<scalar>& def,
    scalar parameter,
    vector& r
)
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
        
    gsMatrix<scalar> basePnt;
    curve.eval_into(parMat,basePnt);

    gsMatrix<scalar> defPnt;
    def.eval_into(parMat,defPnt);
    
    gsMatrix<scalar> pnt = basePnt+defPnt;    
    r = vector(pnt(0,0),pnt(1,0),pnt(2,0));
}

void Foam::Structure::constructCoeffDerivedData()
{
    Info<<"constructCoeffDerivedData"<<Foam::endl;
    coeffDerivedCenterline.resize(nR);
    initialRotation.resize(nR);
    coeffDerivedQuaternions.resize(nR);
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        Info<<"rodNumber:"<<rodNumber<<Foam::endl;
        ActiveRodMesh::rodCosserat* rod = Rods[rodNumber];
        const gsNurbs<scalar>& curve = rod->m_Curve;
        std::cout<<"curve.coefs:"<<curve.coefs()<<std::endl;
        const gsNurbs<scalar>& rotation = rod->m_Rot;
        std::cout<<"rotation.coefs:"<<rotation.coefs()<<std::endl;
        int errorCode = rod->bishopFrameS_dQdR0();
        Info<<"-------"<<Foam::endl;
        if(errorCode==0)
            FatalErrorInFunction<<"Failure to compute Quaternion coefficient derivatives"<<exit(FatalError);
        const gsMatrix<scalar>& dQdR = rod->m_Curve_dQdR0;
        std::cout<<"dQdR:"<<dQdR<<std::endl;
        
        coeffDerivedCenterline[rodNumber].resize(numberCoeffs(rodNumber));
        initialRotation[rodNumber] = rotation;
        coeffDerivedQuaternions[rodNumber].resize(numberCoeffs(rodNumber));
        for(label coeffNumber=0; coeffNumber<numberCoeffs(rodNumber); coeffNumber++)
        {
            Info<<"coeffNumber:"<<coeffNumber<<Foam::endl;
            //Create centerline coefficient derivative curve
            coeffDerivedCenterline[rodNumber][coeffNumber].resize(3);
            for(label dim=0; dim<3; dim++)
            {
                coeffDerivedCenterline[rodNumber][coeffNumber][dim] = curve;
                gsNurbs<scalar>& curve = coeffDerivedCenterline[rodNumber][coeffNumber][dim];
                gsMatrix<scalar>& coeffs = curve.coefs();
                for(label col=0; col<coeffs.cols(); col++)
                {
                    for(label row=0; row<coeffs.rows(); row++)
                    {
                        if(coeffNumber==col && row==dim)
                            coeffs(col,row) = 1;
                        else
                            coeffs(col,row) = 0;
                    }
                }

                //Create quaternion coefficient derivative curve
                coeffDerivedQuaternions[rodNumber][coeffNumber].resize(3);
                coeffDerivedQuaternions[rodNumber][coeffNumber][dim] = rotation;
                gsNurbs<scalar>& rotation = coeffDerivedCenterline[rodNumber][coeffNumber][dim];
                gsMatrix<scalar>& qcoeffs = rotation.coefs();
                label coeffCol = coeffNumber*3+dim;
                List<scalar> dq_doneCoeff(dQdR.rows());
                for(label row=0; row<dQdR.rows(); row++)
                    dq_doneCoeff[row] = dQdR(coeffCol,row);
                label listIndex=0;
                for(label col=0; col<qcoeffs.cols(); col++)
                {
                    for(label row=0; row<qcoeffs.rows(); row++)
                    {
                        qcoeffs(col,row) = dq_doneCoeff[listIndex];
                        listIndex++;
                    }
                }
            }
        }
    }
    Info<<"constructCoeffDerivedData done"<<Foam::endl;
    constructedCoeffDerivedData = true;
}

Foam::FixedList<gsMatrix<Foam::scalar>,3> Foam::Structure::compute_dRdq
(
    const gsMatrix<scalar>& q
)
{
    if(q.rows()!=4 || q.cols()!=1)
        FatalErrorInFunction<<"Invalid dimension for q"<<exit(FatalError);
    
    gsMatrix<scalar> dd1dq(3,4);
        dd1dq(0,0)= 0;        dd1dq(0,1)= 0;        dd1dq(0,2)=-4*q(2,0);   dd1dq(0,3)=-4*q(3,0);
        dd1dq(1,0)= 2*q(3,0); dd1dq(1,1)= 2*q(2,0); dd1dq(1,2)= 2*q(1,0);   dd1dq(1,3)= 2*q(0,0);
        dd1dq(2,0)=-2*q(2,0); dd1dq(2,1)= 2*q(3,0); dd1dq(2,2)=-2*q(0,0);   dd1dq(2,3)= 2*q(1,0);
    
    gsMatrix<scalar> dd2dq(3,4);
        dd2dq(0,0)=-2*q(3,0); dd2dq(0,1)= 2*q(2,0); dd2dq(0,2)= 2*q(1,0);   dd2dq(0,3)=-2*q(0,0);
        dd2dq(1,0)= 0;        dd2dq(1,1)=-4*q(1,0); dd2dq(1,2)= 0;          dd2dq(1,3)=-4*q(3,0);
        dd2dq(2,0)= 2*q(1,0); dd2dq(2,1)= 2*q(0,0); dd2dq(2,2)=-2*q(3,0);   dd2dq(2,3)= 2*q(2,0);  
    
    gsMatrix<scalar> dd3dq(3,4);
        dd3dq(0,0)= 2*q(2,0); dd3dq(0,1)= 2*q(3,0); dd3dq(0,2)= 2*q(0,0);   dd3dq(0,3)= 2*q(1,0);
        dd3dq(1,0)=-2*q(1,0); dd3dq(1,1)=-2*q(0,0); dd3dq(1,2)= 2*q(3,0);   dd3dq(1,3)= 2*q(2,0);
        dd3dq(2,0)= 0;        dd3dq(2,1)=-4*q(1,0); dd3dq(2,2)=-4*q(2,0);   dd3dq(2,3)= 0;
        
    return {dd1dq,dd2dq,dd3dq};
}

gsMatrix<Foam::scalar> Foam::Structure::quaternionMultiply
(
    const gsMatrix<scalar>& q1,
    const gsMatrix<scalar>& q2
)
{
    if(q1.rows()!=4 || q1.cols()!=1)
        FatalErrorInFunction<<"Invalid dimension for q1"<<exit(FatalError);
    if(q2.rows()!=4 || q2.cols()!=1)
        FatalErrorInFunction<<"Invalid dimension for q2"<<exit(FatalError);
    gsMatrix<scalar> result = q1;
    result(0,0) = q1(0,0)*q2(0,0) - q1(1,0)*q2(1,0) - q1(2,0)*q2(2,0) - q1(3,0)*q2(3,0);
    result(1,0) = q1(0,0)*q2(1,0) + q1(1,0)*q2(0,0) + q1(2,0)*q2(3,0) - q1(3,0)*q2(2,0);
    result(2,0) = q1(0,0)*q2(2,0) - q1(1,0)*q2(3,0) + q1(2,0)*q2(0,0) + q1(3,0)*q2(1,0);
    result(3,0) = q1(0,0)*q2(3,0) + q1(1,0)*q2(2,0) - q1(2,0)*q2(1,0) + q1(3,0)*q2(0,0);
    return result;
}

gsMatrix<Foam::scalar> Foam::Structure::quaternionInvert
(
    const gsMatrix<scalar>& q
)
{
    if(q.rows()!=4 || q.cols()!=1)
        FatalErrorInFunction<<"Invalid dimension for q"<<exit(FatalError);
    scalar normVal = q(0,0)*q(0,0)+q(1,0)*q(1,0)+q(2,0)*q(2,0)+q(3,0)*q(3,0);
    if(std::abs(normVal)<1e-10)
        FatalErrorInFunction<<"Invalid quaternion length"<<exit(FatalError);
    gsMatrix<scalar> result = q;
    result(0,0) /=  normVal;
    result(1,0) /= -normVal;
    result(2,0) /= -normVal;
    result(3,0) /= -normVal;
    return result;
}

Foam::label Foam::Structure::numberCoeffs
(
    label rodNumber
) const 
{
    /*
    const ActiveRodMesh::rodCosserat* rod = Rods[rodNumber];
    const gsNurbs<scalar>& curve = rod->m_Curve;
    std::cout<<"curve.coefs():"<<curve.coefs()<<std::endl;
    return curve.coefs().rows();
    */
    
    return rodsList[rodNumber].coefs().rows();
}

void Foam::Structure::rodEvalDerivCoeff
(
    label rodNumber,
    label derivCoeffNumber,
    label derivDimension,
    scalar parameter,
    vector& d1dC,
    vector& d2dC,
    vector& d3dC,
    vector& rdC
)
{
    if(!constructedCoeffDerivedData)
        FatalErrorInFunction<<"Data for deriv coeff data not given!"<<exit(FatalError);
        
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    // Derivative of rod center line in respect to coefficient
    const gsNurbs<scalar>& oneCoeffDerivCenterline = coeffDerivedCenterline[rodNumber][derivCoeffNumber][derivDimension];
    gsMatrix<scalar> coeffDerivCenterEval;
    oneCoeffDerivCenterline.eval_into(parMat,coeffDerivCenterEval);    
    rdC = vector(coeffDerivCenterEval(0,0),coeffDerivCenterEval(1,0),coeffDerivCenterEval(2,0));
    
    //Derivative of rod coordinate axis in respect to coefficient
    const gsNurbs<scalar>& oneCoeffDerivQuaternions = coeffDerivedQuaternions[rodNumber][derivCoeffNumber][derivDimension];
    gsMatrix<scalar> coeffDerivQuaternionEval;
    oneCoeffDerivQuaternions.eval_into(parMat,coeffDerivQuaternionEval);
    
    const gsNurbs<scalar>& totalQuaternions = Rods[rodNumber]->m_Rot;
    gsMatrix<scalar> totalQuaternionEval;
    totalQuaternions.eval_into(parMat,totalQuaternionEval);
    
    const gsNurbs<scalar>& initialQuaternions = initialRotation[rodNumber];
    gsMatrix<scalar> initialQuaternionEval;
    initialQuaternions.eval_into(parMat,initialQuaternionEval);
    
    gsMatrix<scalar> deformationQuaternionEval = quaternionMultiply(totalQuaternionEval,
                                                                    quaternionInvert(initialQuaternionEval));
    gsMatrix<scalar> defQuat_inidQuatdCoeff = quaternionMultiply(deformationQuaternionEval,coeffDerivQuaternionEval);
    FixedList<gsMatrix<Foam::scalar>,3> dRdq = compute_dRdq(totalQuaternionEval);
    for(label n=0; n<dRdq.size(); n++)
    {
        dRdq[n] = dRdq[n]*defQuat_inidQuatdCoeff;
        if(dRdq[n].cols()!=1 || dRdq[n].rows()!=3)
            FatalErrorInFunction<<"Invalid d dimension"<<exit(FatalError);
    }
    for(label dim=0; dim<3; dim++)
    {
        d1dC[dim] = dRdq[0](dim,0);
        d2dC[dim] = dRdq[1](dim,0);
        d3dC[dim] = dRdq[2](dim,0);
    }
}

void Foam::Structure::getCurveCoeffs
(
    List<List<vector>>& coeffs
) const
{
    coeffs.resize(nR);
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        List<vector>& rodCoeffs = coeffs[rodNumber];
        rodCoeffs.resize(numberCoeffs(rodNumber));
        const gsMatrix<scalar>& curveCoeffs = rodsList[rodNumber].coefs();
        for(label coeffNumber=0; coeffNumber<numberCoeffs(rodNumber); coeffNumber++)
        {
            rodCoeffs[coeffNumber] = vector(curveCoeffs(coeffNumber,0),curveCoeffs(coeffNumber,1),curveCoeffs(coeffNumber,2));
        }
    }
}

void Foam::Structure::setCurveCoeffs
(
    const List<List<vector>>& coeffs
)
{
    if(coeffs.size()!=nR)
        FatalErrorInFunction<<"Mismatch in rod number to coeffs!"<<exit(FatalError);
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        const List<vector>& rodCoeffs = coeffs[rodNumber];
        gsMatrix<scalar>& curveCoeffs = rodsList[rodNumber].coefs();
        for(label coeffNumber=0; coeffNumber<numberCoeffs(rodNumber); coeffNumber++)
        {
            curveCoeffs(coeffNumber,0) = rodCoeffs[coeffNumber][0];
            curveCoeffs(coeffNumber,1) = rodCoeffs[coeffNumber][1];
            curveCoeffs(coeffNumber,2) = rodCoeffs[coeffNumber][2];
        }
    }
    setupActiveRodMesh();
}

void Foam::Structure::setCurveCoeff
(
    label rodNumber,
    label derivCoeffNumber,
    label dimension,
    scalar value
)
{
    FatalErrorInFunction<<"Insufficient implementation"<<exit(FatalError);
    
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber"<<exit(FatalError);
    ActiveRodMesh::rodCosserat* rod = Rods[rodNumber];
    gsNurbs<scalar>& curve = rod->m_Curve;
    gsMatrix<scalar>& coeffs = curve.coefs();
    if(derivCoeffNumber<0 || derivCoeffNumber>=coeffs.cols())
        FatalErrorInFunction<<"Invalid derivCoeffNumber"<<exit(FatalError);
    if(dimension<0 || dimension>=3)
        FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    coeffs(derivCoeffNumber,dimension) = value;
    
    //std::unique_ptr<std::vector<std::unique_ptr<gsNurbs<scalar>>>>& rodCoeffDerivedCurves = coeffDerivedCurves[rodNumber];
    //rodCoeffDerivedCurves.reset();
}

/*
BoundingBox Foam::Structure::computeBox
(
    const gsMatrix<scalar>& coefs
)
{
    if(coefs.rows()!=3 || coefs.cols()<1)
        FatalErrorInFunction<<"Invalid coefs size"<<exit(FatalError);
    vector lowerCurve,upperCurve;
    lowerCurve = upperCurve = vector(coefs(0,0),coefs(0,1),coefs(0,2));
    for(label col=0; col<coefs.cols(); col++)
    {
        for(label row=0; row<3; row++)
        {
            if(lowerCurve[row]>coefs(col,row))
                lowerCurve[row] = coefs(col,row);
            if(upperCurve[row]<coefs(col,row))
                upperCurve[row] = coefs(col,row);
        }
    }
    return BoundingBox(lowerCurve,upperCurve);
}
*/

Foam::BoundingBox Foam::Structure::computeBox
(
    label rodNumber
)
{
    const ActiveRodMesh::rodCosserat* rod = Rods[rodNumber];
    const gsNurbs<scalar>& curve = rod->m_Curve;
    BoundingBox curve_box = BoundingBox::boundsOfNurbs(curve);
    const gsNurbs<scalar>& deformation = rod->m_Def;
    BoundingBox def_box = BoundingBox::boundsOfNurbs(deformation);
    return curve_box+def_box;
}

Foam::BoundingBox Foam::Structure::computeBox
(
    label rodNumber,
    scalar parStart,
    scalar parEnd
)
{
    const ActiveRodMesh::rodCosserat* rod = Rods[rodNumber];

    if(parStart<rod->m_Curve.domainStart())
        FatalErrorInFunction<<"Out of range"<<exit(FatalError);
    if(parEnd<rod->m_Curve.domainEnd())
        FatalErrorInFunction<<"Out of range"<<exit(FatalError);
    if(parEnd<parStart)
        FatalErrorInFunction<<"End smaller than start"<<exit(FatalError);
    
    const gsNurbs<scalar>& curve = rod->m_Curve;
    const gsNurbs<scalar>& def = rod->m_Def;

    return BoundingBox::boundsOfNurbs(curve,parStart,parEnd)+BoundingBox::boundsOfNurbs(def,parStart,parEnd);
}

Foam::scalar Foam::Structure::characteristicSize
(
    label rodNumber,
    scalar par
)
{
    vector position = rodEval(Rods[rodNumber],par);
    label cell = mesh.findCell(position);
    if(cell!=-1)
        return initialSpacingFromMesh(mesh,cell);
    else
        return initialSpacingFromMesh(mesh);
}

void Foam::Structure::buildTrees()
{
    rodTrees.resize(nR);
    for(label rodI=0; rodI<nR; rodI++)
        buildTreeOnRod(rodI);
    
}

void Foam::Structure::buildTreeOnRod
(
    label rodNumber
)
{
    const ActiveRodMesh::rodCosserat* rod = Rods[rodNumber];
    
    BoundingBoxTree& rodTree = rodTrees[rodNumber];
    std::unique_ptr<BoundingBoxTree::Node>& root = rodTree.getRoot();
    root = std::make_unique<BoundingBoxTree::Node>();
    
    root->key = 0.5*(rod->m_Curve.domainStart()+rod->m_Curve.domainEnd());
    root->value = computeBox(rodNumber);
    root->leftChild = std::make_unique<BoundingBoxTree::Node>();
    root->rightChild = std::make_unique<BoundingBoxTree::Node>();

    std::function<void(std::unique_ptr<BoundingBoxTree::Node>&,std::pair<scalar,scalar>)> recursiveBuildTree =
    [&](std::unique_ptr<BoundingBoxTree::Node>& node, std::pair<scalar,scalar> bound)
    {
        scalar center = 0.5*(bound.first+bound.second);
        node->key = center;
        node->value = computeBox(rodNumber,bound.first,bound.second);
        scalar boxSize = node->value.innerSize();
        scalar charSize = characteristicSize(rodNumber,center);
        if(charSize<boxSize)
        {
            node->leftChild = std::make_unique<BoundingBoxTree::Node>();
            recursiveBuildTree(node->leftChild,{bound.first,center});
            node->rightChild = std::make_unique<BoundingBoxTree::Node>();
            recursiveBuildTree(node->rightChild,{center,bound.second});
        }        
    };
}

void Foam::Structure::setupRodBoundingBoxTree()
{
    buildTrees();
    rodInMesh.resize(nR);
    for(label rodI=0; rodI<nR; rodI++)
    {
        for(const vector& vertice : meshBoundingBox.allVertices())
        {
            std::vector<scalar> parameters;
            rodTrees[rodI].findPointParameters(parameters,vertice);
            if(parameters.empty())
                rodInMesh[rodI] = false;
            else
                rodInMesh[rodI] = true;
        }
    }
}

Foam::BoundingBox Foam::Structure::computeMeshBoundingBox()
{
    vector smaller = mesh.points()[0];
    vector larger = mesh.points()[0];
    for(const vector& pnt : mesh.points())
    {
        for(label d=0; d<3; d++)
        {
            smaller[d] = std::min(pnt[d],smaller[d]);
            larger[d] = std::max(pnt[d],larger[d]);
        }
    }
    return BoundingBox(smaller,larger);
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

Foam::scalar Foam::Structure::supportDomainMinSize
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

/*
template<typename T>
std::unique_ptr<List<List<T>>> Foam::Structure::broadcastHaloFields
(
    const GeometricField<T,fvPatchField,volMesh>& fieldData
)
{
    auto haloFieldPtr = std::make_unique<List<List<T>>>(Pstream::nProcs());
    List<List<T>>& haloField = *haloFieldPtr;
    List<T>& ownHaloField = haloField[Pstream::myProcNo()];
    
    const DynamicList<CellDescription>& ownHaloCells = getHaloCellList(Pstream::myProcNo());    
    ownHaloField.setSize(ownHaloCells.size());
    
    for(label haloCellInd=0; haloCellInd<ownHaloCells.size(); haloCellInd++)
    {
        ownHaloField[haloCellInd] = fieldData[ownHaloCells[haloCellInd].index];
    }
    
    return haloFieldPtr;
}
*/

void Foam::Structure::generateMeshGraph()
{
    const cellList& cells = mesh.cells();
    const labelList& owners = mesh.owner();
    const labelList& neighbours = mesh.neighbour();
    const std::unordered_map<label,List<DynamicList<std::pair<label,label>>>>& patchFaceToCell = getPatchFaceToCellMap();
    
    meshGraph.setSize(cells.size());
    for(label cellInd=0; cellInd<cells.size(); cellInd++)
    {
        const cell& oneCell = cells[cellInd];
        meshGraph[cellInd].setSize(oneCell.size());
        for(label cellFaceInd=0; cellFaceInd<oneCell.size(); cellFaceInd++)
        {
            label faceInd = oneCell[cellFaceInd];            
            Pair<label> neighData = {-1,-1};
            if(owners[faceInd]!=cellInd)
            {
                neighData = {Pstream::myProcNo(),owners[faceInd]};
            }
            else if(faceInd<neighbours.size())
            {
                neighData = {Pstream::myProcNo(),neighbours[faceInd]};
            }
            else if(patchFaceToCell.find(faceInd)!=patchFaceToCell.end())
            {
                auto iter = patchFaceToCell.find(faceInd);
                const List<DynamicList<std::pair<label,label>>>& thisPatchFaceNeighbour = iter->second;
                if(thisPatchFaceNeighbour.size()<1)
                    FatalErrorInFunction<<"Neighbourhood must be larger than one iteration"<<exit(FatalError);
                const DynamicList<std::pair<label,label>>& firstIterNeighbourhood = thisPatchFaceNeighbour[0];
                if(firstIterNeighbourhood.size()!=1)
                    FatalErrorInFunction<<"First iteration neighbourhood must be sized 1"<<exit(FatalError);
                std::pair<label,label> neighbourCell = firstIterNeighbourhood[0];
                neighData = {neighbourCell.first,neighbourCell.second};
            }

            meshGraph[cellInd][cellFaceInd] = neighData;
        }
    }
    
    globalHaloMeshGraph.setSize(Pstream::nProcs());
    globalHaloMeshGraph[Pstream::myProcNo()] = meshGraph;
    
    Pstream::gatherList(globalHaloMeshGraph);
    Pstream::scatterList(globalHaloMeshGraph);    
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
    //[procs]->[patch]
    List<DynamicList<label>> procPatchOwner(Pstream::nProcs());
    List<DynamicList<label>> procPatchNeighbour(Pstream::nProcs());
    List<DynamicList<label>> procPatchIndex(Pstream::nProcs());
    List<DynamicList<label>> procPatchFaceLocalIndex(Pstream::nProcs());
    //[procs]->[patch]->[iteration]->[cell]
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
        
    //Compute halo cell list
    globalHaloCellList_Sorted.resize(Pstream::nProcs());
    globalHaloCellToIndexMap.resize(Pstream::nProcs());
    for(label process=0; process<Pstream::nProcs(); process++)
    {
        const DynamicList<label>& patchOwner = procPatchOwner[process];
        const DynamicList<List<DynamicList<label>>>& patchCellInds = procCellInds[process];
        const DynamicList<List<DynamicList<vector>>>& patchCellCentres = procCellCentres[process];
        const DynamicList<List<DynamicList<scalar>>>& patchCellMags = procCellMags[process];
        
        std::unordered_map<label,std::tuple<label,label,label>> uniqueHaloCells;
        for(label patchFaceInd=0; patchFaceInd<patchOwner.size(); patchFaceInd++)
        {
            label owner = patchOwner[patchFaceInd];
            const List<DynamicList<label>>& cellInds = patchCellInds[patchFaceInd];
            
            if(owner!=process)
                FatalErrorInFunction<<"Invalid data"<<exit(FatalError);
            
            for(label iteration=0; iteration<cellInds.size(); iteration++)
            {
                const DynamicList<label>& iterCellInds = cellInds[iteration];               
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
        //const DynamicList<List<DynamicList<label>>>& patchCellInds = procCellInds[process];
        
        for(label patchFaceInd=0; patchFaceInd<patchOwner.size(); patchFaceInd++)
        {
            label owner = patchOwner[patchFaceInd];
            label neighbour = patchNeighbour[patchFaceInd];
            label index = patchIndex[patchFaceInd];
            label faceLocalIndex = patchFaceLocalIndex[patchFaceInd];
            //const List<DynamicList<label>>& cellInds = patchCellInds[patchFaceInd];
            
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
            //const DynamicList<label>& patchIndex = procPatchIndex[neighborProcess];
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
                //label index = patchIndex[patchFaceInd];
                label faceLocalIndex = patchFaceLocalIndex[patchFaceInd];
                const List<DynamicList<label>>& cellInds = patchCellInds[patchFaceInd];
                
                if(owner!=neighborProcess)
                    FatalErrorInFunction<<"Data mismatch"<<exit(FatalError);
                if(neighbour!=Pstream::myProcNo())
                    FatalErrorInFunction<<"Data mismatch"<<exit(FatalError);
                if(faceLocalIndex!=locFacePatchInd)
                    FatalErrorInFunction<<"Data mismatch"<<exit(FatalError);

                label faceInd = patchStartFace+locFacePatchInd;
                if(cellInds.size()!=iterations)
                    FatalErrorInFunction<<"Size mismatch of iterations and cell inds size"<<exit(FatalError);
                if(patchFaceToCellMap.find(faceInd)!=patchFaceToCellMap.end())
                    FatalErrorInFunction<<"Duplicate halo face to cells"<<exit(FatalError);

                patchFaceToCellMap[faceInd].setSize(cellInds.size());
                for(label iteration=0; iteration<cellInds.size(); iteration++)
                {
                    const DynamicList<label>& iterCell = cellInds[iteration];
                    for(label cellInd : iterCell)
                    {
                        patchFaceToCellMap[faceInd][iteration].append({neighborProcess,cellInd});
                    }                    
                }
            }
        }
    }
    
    generateMeshGraph();
}

Foam::scalar Foam::Structure::initialSpacingFromMesh
(
    const fvMesh& mesh,
    label cellInd
)
{
    if(mesh.cells().size()<1)
        FatalErrorInFunction<<"Mesh not valid"<<exit(FatalError);
    const cell& oneCell = mesh.cells()[cellInd];
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

void Foam::Barrier(bool stop)
{
    labelList test(Pstream::nProcs(),0);
    test[Pstream::myProcNo()] = Pstream::myProcNo();
    Pstream::gatherList(test);
    Pstream::scatterList(test);
    std::cout<<"XXXXXXXXXXXXX---------------Barrier:"<<Pstream::myProcNo()<<"----"<<test[0]<<","<<test[1]<<"---------------XXXXXXXXXXXX"<<std::endl;
    fflush(stdout);
    Pstream::gatherList(test);
    if(stop)
        FatalErrorInFunction<<Pstream::myProcNo()<<"Temp Stop!"<< exit(FatalError);

}
