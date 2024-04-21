#include "Structure.H"

Foam::Structure::Structure
(
    const Time& runTime,
    const dimensionedScalar& alpha,
    const volScalarField& T,
    const volScalarField& p,
    const volVectorField& U,
    dynamicRefineFvMesh& mesh,
    const dimensionedScalar nu
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
nu(nu)
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
    
    //createDeformationCurve();    
    //assignBoundaryFacesToNurbsCurves();
    Info<<"Completed Structure setup"<<Foam::endl;
}

void Foam::Structure::assignBoundaryFacesToNurbsCurves()
{
}

template<typename Tensor_Type>
std::unique_ptr<Foam::List<Foam::List<Tensor_Type>>> Foam::Structure::computeDistributedLoad
(
    const Foam::Field<Tensor_Type> immersedBoundaryField
)
{
}

void Foam::Structure::assignForceOnCurve()
{
}

void Foam::Structure::computeIBHeatFlux()
{
}

Foam::word Foam::Structure::getXMLPath()
{
    //Info<<"lateScale:"<<latScale<<endl;
    //Info<<"latDir:"<<latDir<<endl;

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
    Info<<"Xml path:"<<fullPath<<Foam::endl;
    return fullPath;
}

int Foam::Structure::loadRodsFromXML()
{
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
    forceCurveStorage.resize(nR);
    
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

void Foam::Structure::solveOneStep()
{
    Info<<"Solve Nurbs structure"<<Foam::endl;
    assignForceOnCurve();

    if(Pstream::master())
        myMesh->solve(1.0,solveOpt);
    
    moveNurbs();
}

void Foam::Structure::moveNurbs()
{
    Info<<"Mesh assign Deformation curve"<<Foam::endl;
    label nbrNurbs = myMesh->m_Rods.size();
    std::vector<List<List<vector>>> nurbs_to_new_controlPoints;
    
    //Insert deformation CP_s into data structure in master and resize otherwise
    if(Pstream::master())
    {
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
    }
    else
    {
        nurbs_to_new_controlPoints.resize(nbrNurbs);
    }
    
    //Transfer first dimension of CP list
    labelList controlPntDim1Size(nbrNurbs,0);
    if(Pstream::master())
    {
        for(label nurbsInd=0;nurbsInd<nbrNurbs;nurbsInd++)
        {
            controlPntDim1Size[nurbsInd] = nurbs_to_new_controlPoints[nurbsInd].size();
        }
    }
    Pstream::scatterList(controlPntDim1Size);
    if(!Pstream::master())
    {
        for(label nurbsInd=0;nurbsInd<nbrNurbs;nurbsInd++)
        {
            nurbs_to_new_controlPoints[nurbsInd].setSize(controlPntDim1Size[nurbsInd]);
        }
    }
    
    //Transfer second dimension of CP list
    DynamicList<scalar> CP_data;
    for(label nurbsInd=0;nurbsInd<nbrNurbs;nurbsInd++)
    {
        labelList controlPntDim2Size(nurbs_to_new_controlPoints[nurbsInd].size(),0);
        if(Pstream::master())
        {
            for(label firstDimInd=0;firstDimInd<nurbs_to_new_controlPoints[nurbsInd].size();firstDimInd++)
            {
                controlPntDim2Size[firstDimInd] = nurbs_to_new_controlPoints[nurbsInd][firstDimInd].size();
            }
        }
        Pstream::scatterList(controlPntDim2Size);
        if(!Pstream::master())
        {
            for(label firstDimInd=0;firstDimInd<nurbs_to_new_controlPoints[nurbsInd].size();firstDimInd++)
            {
                nurbs_to_new_controlPoints[nurbsInd][firstDimInd].setSize(controlPntDim2Size[firstDimInd]);
            }
        }
        
        for(label firstDimInd=0;firstDimInd<nurbs_to_new_controlPoints[nurbsInd].size();firstDimInd++)
        {
            for(label secondDimInd=0;
                secondDimInd<nurbs_to_new_controlPoints[nurbsInd][firstDimInd].size();
                secondDimInd++
               )
            {
                for(label vecDimInd=0;vecDimInd<3;vecDimInd++)
                {
                    if(Pstream::master())
                    {
                        CP_data.append(nurbs_to_new_controlPoints[nurbsInd][firstDimInd][secondDimInd][vecDimInd]);
                    }
                    else
                    {
                        CP_data.append(0);
                    }
                }
            }
        }
    }
    Pstream::scatterList(CP_data);
    label CP_data_index=0;
    if(!Pstream::master())
    {
        for(label nurbsInd=0;nurbsInd<nbrNurbs;nurbsInd++)
        {
            for(label firstDimInd=0;firstDimInd<nurbs_to_new_controlPoints[nurbsInd].size();firstDimInd++)
            {
                for(label secondDimInd=0;
                    secondDimInd<nurbs_to_new_controlPoints[nurbsInd][firstDimInd].size();
                    secondDimInd++
                   )
                {
                    for(label vecDimInd=0;vecDimInd<3;vecDimInd++)
                    {
                        if(CP_data_index<CP_data.size())
                            FatalErrorInFunction<<"Invalid index"<< exit(FatalError); 

                        nurbs_to_new_controlPoints[nurbsInd][firstDimInd][secondDimInd][vecDimInd] = CP_data[CP_data_index];
                        CP_data_index++;
                    }
                }
            }
        }
    }    
}

void Foam::Structure::createDeformationCurve()
{
    /*
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
    */
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
