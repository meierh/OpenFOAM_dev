#include "Structure.H"

Foam::Quaternion::Quaternion
(
    scalar x,
    scalar y,
    scalar z,
    scalar w
):
x(x),
y(y),
z(z),
w(w)
{
}

Foam::Quaternion::Quaternion
(
    const gsMatrix<scalar>& gsQuaternion
)
{
    if(gsQuaternion.rows()!=4 || gsQuaternion.cols()!=1)
        FatalErrorInFunction<<"Invalid size of gsQuaternion"<<exit(FatalError);
    
    x = gsQuaternion(0,0);
    y = gsQuaternion(1,0);
    z = gsQuaternion(2,0);
    w = gsQuaternion(3,0);
}
        
Foam::Quaternion Foam::Quaternion::operator*
(
    Quaternion const& q
) const 
{
    Quaternion result;
    result.x = w*q.x - x*q.w - y*q.z - z*q.y;
    result.y = w*q.y - x*q.z - y*q.w - z*q.x;
    result.z = w*q.z - x*q.y - y*q.x - z*q.w;
    result.w = w*q.w - x*q.x - y*q.y - z*q.z;
    return result;
}

Foam::Quaternion Foam::Quaternion::operator/
(
    Quaternion const& q
) const
{
    Quaternion invQ = q.invert();
    return (*this)*invQ;
}

Foam::Quaternion Foam::Quaternion::operator-
(
    Quaternion const& q
) const
{
    Quaternion result;
    for(label i=0; i<4; i++)
        result[i] = (*this)[i]-q[i];
    return result;
}

Foam::Quaternion Foam::Quaternion::invert() const
{
    Quaternion invQ = *this;
    scalar absInvQ = invQ.len();
    absInvQ *= absInvQ;
    if(absInvQ<1e-10)
        FatalErrorInFunction<<"Invalid quaternion length"<<exit(FatalError);
    invQ.w /=  absInvQ;
    invQ.x /= -absInvQ;
    invQ.y /= -absInvQ;
    invQ.z /= -absInvQ;
    return invQ;
}

Foam::scalar Foam::Quaternion::len() const
{
    return std::sqrt(x*x + y*y + z*z + w*w);
}

Foam::scalar Foam::Quaternion::distanceNorm2
(
    Quaternion const& q
) const
{
    Quaternion dq;
    dq.w = w-q.w;
    dq.x = x-q.x;
    dq.y = y-q.y;
    dq.z = z-q.z;
    return dq.len();
}

void Foam::Quaternion::normalize()
{
    scalar len = this->len();
    if(len==0)
    {
        x = 1;
    }
    else
    {
        w /= len;
        x /= len;
        y /= len;
        z /= len;
    }
}

Foam::scalar& Foam::Quaternion::operator[]
(
    uint index
)
{
    switch(index)
    {
        case 0:
            return w;
        case 1:
            return x;
        case 2:
            return y;
        case 3:
            return z;
        default:
            FatalErrorInFunction<<"Invalid index in quaternion"<<exit(FatalError);
            return w;
    }
}

Foam::scalar Foam::Quaternion::operator[]
(
    uint index
) const
{
    switch(index)
    {
        case 0:
            return w;
        case 1:
            return x;
        case 2:
            return y;
        case 3:
            return z;
        default:
            FatalErrorInFunction<<"Invalid index in quaternion"<<exit(FatalError);
            return w;
    }
}

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    Quaternion const& q
)
{
    return os << "[("<<q.x<<","<<q.y<<","<<q.z<<")("<<q.w<<")]";
}


Foam::Rotation::Rotation
(
    vector d1,
    vector d2,
    vector d3
)
{
    T[0] = d1;
    T[1] = d2;
    T[2] = d3;
}

Foam::Rotation::Rotation
(
    const Quaternion& q
)
{
    vector d1 = vector
    (
        q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z,
         2*q.w*q.z + 2*q.x*q.y,
        -2*q.w*q.y + 2*q.x*q.z
    );
    vector d2 = vector
    (
        -2*q.w*q.z + 2*q.x*q.y,
        q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z,
         2*q.w*q.x + 2*q.y*q.z
    );
    vector d3 = vector
    (
         2*q.w*q.y + 2*q.x*q.z,
        -2*q.w*q.x + 2*q.y*q.z,
        q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z
    );
    T = {d1,d2,d3};
}

Foam::Rotation Foam::Rotation::operator-
(
    Rotation const& R
) const
{
    Rotation result;
    for(label d=0; d<3; d++)
    {
        result.T[d] = (T[d]-R.T[d]); 
    }
    return result;
}

Foam::Rotation Foam::Rotation::operator+
(
    Rotation const& R
) const
{
    Rotation result;
    for(label d=0; d<3; d++)
    {
        result.T[d] = (T[d]+R.T[d]); 
    }
    return result;
}

bool Foam::Rotation::operator!=
(
    Rotation const& R
) const
{
    return this->T!=R.T;
}

Foam::Rotation Foam::Rotation::operator/
(
    scalar alpha
) const
{
    Rotation result = *this;
    for(label d=0; d<3; d++)
    {
        result.T[d] /= alpha; 
    }
    return result;
}

Foam::scalar Foam::Rotation::distanceNorm2
(
    Rotation const& R
) const
{
    Rotation diff = *this - R;
    return diff.norm2();
}

Foam::scalar Foam::Rotation::norm2() const
{
    scalar sum = 0;
    for(label i=0; i<3; i++)
        for(label j=0; j<3; j++)
            sum += T[i][j]*T[i][j];
    return std::sqrt(sum);
}

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    Rotation const& m
)
{
    return os << m.T;
}

Foam::Rotation Foam::Rotation::compute_dRdX
(
    const Quaternion& dqdX,
    const Quaternion& q
)
{   
    const FixedList<FixedList<vector,4>,3> dRdq = compute_dRdq(q);
 
    Rotation dRdC;
    std::vector<vector*> ddkdCPtr = {&(dRdC.T[0]),&(dRdC.T[1]),&(dRdC.T[2])};
    for(label dk=0; dk<3; dk++)
    {
        const FixedList<vector,4>& dRkdq = dRdq[dk];
        vector& ddkdC = *(ddkdCPtr[dk]);
        for(label dim=0; dim<3; dim++)
        {
            ddkdC[dim] = 0;
            ddkdC[dim] += dRkdq[0][dim] * dqdX.qw();
            ddkdC[dim] += dRkdq[1][dim] * dqdX.qx();
            ddkdC[dim] += dRkdq[2][dim] * dqdX.qy();
            ddkdC[dim] += dRkdq[3][dim] * dqdX.qz();
        }
    }
    return dRdC;
}

Foam::FixedList<Foam::FixedList<Foam::vector,4>,3> Foam::Rotation::compute_dRdq
(
    const Quaternion& q
)
{
    FixedList<vector,4> dd1dq;
        dd1dq[0][0]= 2*q.qw(); dd1dq[1][0]= 2*q.qx(); dd1dq[2][0]=-2*q.qy(); dd1dq[3][0]=-2*q.qz();
        dd1dq[0][1]= 2*q.qz(); dd1dq[1][1]= 2*q.qy(); dd1dq[2][1]= 2*q.qx(); dd1dq[3][1]= 2*q.qw();
        dd1dq[0][2]=-2*q.qy(); dd1dq[1][2]= 2*q.qz(); dd1dq[2][2]=-2*q.qw(); dd1dq[3][2]= 2*q.qx();
    
    FixedList<vector,4> dd2dq;
        dd2dq[0][0]=-2*q.qz(); dd2dq[1][0]= 2*q.qy(); dd2dq[2][0]= 2*q.qx(); dd2dq[3][0]=-2*q.qw();
        dd2dq[0][1]= 2*q.qw(); dd2dq[1][1]=-2*q.qx(); dd2dq[2][1]= 2*q.qy(); dd2dq[3][1]=-2*q.qz();
        dd2dq[0][2]= 2*q.qx(); dd2dq[1][2]= 2*q.qw(); dd2dq[2][2]= 2*q.qz(); dd2dq[3][2]= 2*q.qy();
        
    FixedList<vector,4> dd3dq;
        dd3dq[0][0]= 2*q.qy(); dd3dq[1][0]= 2*q.qz(); dd3dq[2][0]= 2*q.qw(); dd3dq[3][0]= 2*q.qx();
        dd3dq[0][1]=-2*q.qx(); dd3dq[1][1]=-2*q.qw(); dd3dq[2][1]= 2*q.qz(); dd3dq[3][1]= 2*q.qy();
        dd3dq[0][2]= 2*q.qw(); dd3dq[1][2]=-2*q.qx(); dd3dq[2][2]=-2*q.qy(); dd3dq[3][2]= 2*q.qz();
        
    return {dd1dq,dd2dq,dd3dq};
}

Foam::Rotation Foam::Rotation::compute_d2RdX
(
    const Quaternion& d2qdX,
    const Quaternion& dqdX,
    const Quaternion& q
)
{      
    const FixedList<FixedList<vector,4>,3> d2Rdq = compute_d2Rdq(q);
 
    Rotation d2RdX;
    std::vector<vector*> d2dkdXPtr = {&(d2RdX.T[0]),&(d2RdX.T[1]),&(d2RdX.T[2])};
    for(label dk=0; dk<3; dk++)
    {
        const FixedList<vector,4>& d2Rkdq = d2Rdq[dk];
        vector& d2dkdX = *(d2dkdXPtr[dk]);
        for(label dim=0; dim<3; dim++)
        {
            d2dkdX[dim] = 0;
            d2dkdX[dim] += d2Rkdq[0][dim] * (dqdX.qw()*dqdX.qw());
            d2dkdX[dim] += d2Rkdq[1][dim] * (dqdX.qx()*dqdX.qx());
            d2dkdX[dim] += d2Rkdq[2][dim] * (dqdX.qy()*dqdX.qy());
            d2dkdX[dim] += d2Rkdq[3][dim] * (dqdX.qz()*dqdX.qz());
        }
    }
    
    d2RdX = d2RdX + compute_dRdX(d2qdX,q);
    return d2RdX;
}

Foam::FixedList<Foam::FixedList<Foam::vector,4>,3> Foam::Rotation::compute_d2Rdq
(
    const Quaternion& q
)
{
    FixedList<vector,4> d2d1dq;
        d2d1dq[0][0]= 2; d2d1dq[1][0]= 2; d2d1dq[2][0]=-2; d2d1dq[3][0]=-2;
        d2d1dq[0][1]= 0; d2d1dq[1][1]= 0; d2d1dq[2][1]= 0; d2d1dq[3][1]= 0;
        d2d1dq[0][2]= 0; d2d1dq[1][2]= 0; d2d1dq[2][2]= 0; d2d1dq[3][2]= 0;
    
    FixedList<vector,4> d2d2dq;
        d2d2dq[0][0]= 0; d2d2dq[1][0]= 0; d2d2dq[2][0]= 0; d2d2dq[3][0]= 0;
        d2d2dq[0][1]= 2; d2d2dq[1][1]=-2; d2d2dq[2][1]= 2; d2d2dq[3][1]=-2;
        d2d2dq[0][2]= 0; d2d2dq[1][2]= 0; d2d2dq[2][2]= 0; d2d2dq[3][2]= 0;
        
    FixedList<vector,4> d2d3dq;
        d2d3dq[0][0]= 0; d2d3dq[1][0]= 0; d2d3dq[2][0]= 0; d2d3dq[3][0]= 0;
        d2d3dq[0][1]= 0; d2d3dq[1][1]= 0; d2d3dq[2][1]= 0; d2d3dq[3][1]= 0;
        d2d3dq[0][2]= 2; d2d3dq[1][2]=-2; d2d3dq[2][2]=-2; d2d3dq[3][2]= 2;
        
    return {d2d1dq,d2d2dq,d2d3dq};
}

Foam::Structure::Structure
(
    const fvMesh& mesh,
    const Time& runTime
):
initialMeshSpacing(spacingFromMesh(mesh)),
runTime(runTime),
runDirectory(runTime.rootPath()),
caseName(runTime.caseName()),
xmlPath(getXMLPath()),
name(getName()),
nR(loadRodsFromXML()),
mesh(mesh),
meshBoundingBox(computeMeshBoundingBox())
{
    FatalErrorInFunction<<"Not in use anymore"<<exit(FatalError);
    Info<<"----------------Structure----------------"<<Foam::endl;
    createParallelTopology();
    computeMeshSetup();
    setupActiveRodMesh();
}

Foam::Structure::Structure
(
    const fvMesh& mesh,
    const std::shared_ptr<IOdictionary> structureDict,
    const Time& runTime
):
initialMeshSpacing(spacingFromMesh(mesh)),
runTime(runTime),
runDirectory(runTime.rootPath()),
caseName(runTime.caseName()),
xmlPath(xmlFromDict(*structureDict)),
name(getName()),
nR(loadRodsFromXML()),
mesh(mesh),
structureDict(structureDict),
meshBoundingBox(computeMeshBoundingBox())
{
    Info<<"----------------Structure dir----------------"<<Foam::endl;
    createParallelTopology();
    computeMeshSetup();
    setupActiveRodMesh();
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
    {
        delete Rods[i];
        delete Geo[i];
        delete BasisRef[i];
    }
    Geo.clear();
	Rods.clear();
    myMesh = nullptr;
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
        FatalErrorInFunction<<"Invalid entry in constant/structureDict/rodFile -- must be string"<<exit(FatalError);
    word rodFileWord = rodFileToken.stringToken();
    word caseName = this->caseName;
    if(Pstream::parRun())
    {
        std::size_t sepInd = caseName.find('/');
        if(sepInd==caseName.npos)
            FatalErrorInFunction<<"Invalid caseName:"<<caseName<<exit(FatalError);
        word redCaseName = caseName.substr(0,sepInd);
        caseName = redCaseName;
    }
    fileName caseDirectory = runDirectory+"/"+caseName+"/constant/";
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
    Info<<"Bounding Box (x:["<<x0<<"-"<<x1<<"], y:["<<y0<<"-"<<y1<<"], z:["<<z0<<"-"<<z1<<"])"<<endl;
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
    Info<<"lateScale:"<<latScale<<endl;
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

        /*
        std::cout<<"use_mixed:"<<use_mixed<<std::endl;
        std::cout<<"rodsList[0]:"<<std::endl<<rodsList[0]<<std::endl;
        std::cout<<"*BasisRef[0]:"<<std::endl<<*BasisRef[0]<<std::endl;
        //std::cout<<"*Geo[0]:"<<*Geo[0]<<std::endl;
        std::cout<<"*twist:"<<std::endl<<*twist<<std::endl;
        */

        // Make rod
        if (use_mixed)
            Rods[i] = new ActiveRodMesh::rodCosseratMixed(rodsList[i], *BasisRef[i], *Geo[i], twist, 2, 0);
        else
            Rods[i] = new ActiveRodMesh::rodCosserat(rodsList[i], *BasisRef[i], *Geo[i], twist, 2, 0);
    
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
    Info<<"setupActiveRodMesh"<<Foam::endl;
    for(std::size_t rodI=0; rodI<rodsList.size(); rodI++)
        std::cout<<rodsList[rodI]<<std::endl;
    
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
        const gsNurbs<scalar>& rotation = rod->m_Rot;
        if(curve.domainStart()!=rotation.domainStart())
            FatalErrorInFunction<<"Mismatch in curve and rotation start"<<exit(FatalError);
        if(curve.domainEnd()!=rotation.domainEnd())
            FatalErrorInFunction<<"Mismatch in curve and rotation end"<<exit(FatalError);
        
        
        std::cout<<"curve:"<<curve<<std::endl;
        /*
        std::cout<<"deformation"<<deformation<<std::endl;
        std::cout<<"rotation:"<<rotation<<std::endl;
        */
    }
    constructCoeffDerivedData();
    
    prevDeformations.clear();
    prevDeformations.resize(nR);
    prevRotations.clear();
    prevRotations.resize(nR);
    
    rodEvalBuffer.resize(nR);
    rodDerivEvalBuffer.resize(nR);
    rodDeriv2EvalBuffer.resize(nR);    
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

Foam::scalar Foam::Structure::domainStart
(
    label rodNumber
) const
{
    if(rodNumber<0 || rodNumber>=getNumberRods())
        FatalErrorInFunction<<"Invalid fourierCoeffNumber"<<exit(FatalError);
    return Rods[rodNumber]->m_Curve.domainStart();
}

Foam::scalar Foam::Structure::domainEnd
(
    label rodNumber
) const
{
    if(rodNumber<0 || rodNumber>=getNumberRods())
        FatalErrorInFunction<<"Invalid fourierCoeffNumber"<<exit(FatalError);
    return Rods[rodNumber]->m_Curve.domainEnd();
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

/*
Foam::FixedList<Foam::vector,3> Foam::Structure::quaternionsToRotation
(
    FixedList<scalar,4> quaternions
)
{
    const FixedList<scalar,4>& q = quaternions;
    vector d1 = vector
    (
        q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3],
         2*q[0]*q[3] + 2*q[1]*q[2],
        -2*q[0]*q[2] + 2*q[1]*q[3]
    );
    vector d2 = vector
    (
        -2*q[0]*q[3] + 2*q[1]*q[2],
        q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3],
         2*q[0]*q[1] + 2*q[2]*q[3]
    );
    vector d3 = vector
    (
         2*q[0]*q[2] + 2*q[1]*q[3],
        -2*q[0]*q[1] + 2*q[2]*q[3],
        q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]
    );
    return {d1,d2,d3};
}
*/

Foam::scalar Foam::Structure::matrixDistance
(
    const FixedList<vector,3>& m1,
    const FixedList<vector,3>& m2
)
{
    scalar dist=0;
    for(label j=0; j<3; j++)
    {
        for(label i=0; i<3; i++)
        {
            scalar diff = m1[j][i]-m2[j][i];
            dist += diff*diff;
        }
    }
    return std::sqrt(dist);
}

Foam::scalar Foam::Structure::vectorDistance
(
    const vector& v1,
    const vector& v2
)
{
    scalar dist=0;
    for(label i=0; i<3; i++)
    {
        scalar diff = v1[i]-v2[i];
        dist += diff*diff;
    }
    return std::sqrt(dist);
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
    std::cout<<"cKnots:"<<cKnots<<std::endl;
    std::cout<<"cWeight:"<<cWeight<<std::endl;
    std::cout<<"cCoeff:"<<cCoeff<<std::endl;
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

gsNurbs<Foam::scalar> Foam::Structure::createNurbs
(
    scalar domainStart,
    scalar domainEnd,
    uint degree,
    std::vector<scalar> coefficients
)
{
    std::vector<scalar> knots = computeUniformKnots(degree,coefficients.size(),domainStart,domainEnd);
    std::vector<scalar> weights(coefficients.size());
    std::fill(weights.begin(),weights.end(),1);
    return createNurbs(knots,degree,weights,coefficients);
}

gsNurbs<Foam::scalar> Foam::Structure::createNurbs
(
    scalar domainStart,
    scalar domainEnd,
    uint degree,
    std::vector<vector> coefficients
)
{
    std::vector<scalar> knots = computeUniformKnots(degree,coefficients.size(),domainStart,domainEnd);
    std::vector<scalar> weights(coefficients.size());
    std::fill(weights.begin(),weights.end(),1);
    return createNurbs(knots,degree,weights,coefficients);
}

std::vector<Foam::scalar> Foam::Structure::computeUniformKnots
(
    uint degree,
    uint nbrCoefficients,
    scalar domainStart,
    scalar domainEnd
)
{
    std::vector<scalar> knots;
    scalar deltaSpan = domainEnd-domainStart;
    label centerCoeffs = nbrCoefficients-degree;
    scalar stepsize = deltaSpan/(centerCoeffs);
    
    for(uint p=0; p<=degree; p++)
        knots.push_back(domainStart);
    
    for(label i=0; i<centerCoeffs-1; i++)
        knots.push_back(knots.back()+stepsize);
    
    for(uint p=0; p<=degree; p++)
        knots.push_back(domainEnd);

    return gsKnotVector<Foam::scalar>(knots);
}

void Foam::Structure::fitNurbsCoeffsToPoints
(
    const List<vector>& points,
    const gsNurbs<scalar>& nurbs,
    gsMatrix<scalar>& fittedCoeffs,
    scalar epsilon
)
{
    scalar start = nurbs.domainStart();
    scalar end = nurbs.domainEnd();
    
    List<scalar> parameter;
    List<vector> pnts;
    
    if(points.size()<2)
    {
        parameter.resize(2);
        parameter[0] = start;
        parameter[1] = end;
        pnts.resize(2);
        pnts[0] = points[0];
        pnts[1] = points[0];
    }
    else
    {
        parameter.resize(points.size());
        parameter[0] = start;
        parameter.last() = end;
        scalar spacing = (end-start)/(points.size()-1);
        for(label i=1; i<points.size()-1; i++)
            parameter[i] = parameter[i-1]+spacing;
        pnts = points;
    }
    return fitNurbsCoeffsToPoints(pnts,parameter,nurbs,fittedCoeffs,epsilon);
}

void Foam::Structure::nurbsMinusX
(
    const gsMatrix<scalar>& x,
    const gsMatrix<scalar>& s,
    const gsNurbs<scalar>& nurbs,
    gsMatrix<scalar>& N_x
)
{
    label parameterNbr = x.cols();
    
    if(x.rows()!=3)
        FatalErrorInFunction<<"x row dimension must be 3"<<exit(FatalError);
    if(s.cols()!=parameterNbr)
        FatalErrorInFunction<<"Mismatch x and s dimension"<<exit(FatalError);
    if(s.rows()!=1)
        FatalErrorInFunction<<"s row dimension must be 1"<<exit(FatalError);    
    if(nurbs.targetDim()!=3)
        FatalErrorInFunction<<"Mismatch target dimension"<<exit(FatalError);
    
    if(N_x.rows()!=3 || N_x.cols()!=parameterNbr)
        N_x = gsMatrix<scalar>(3,parameterNbr);
    
    nurbs.eval_into(s,N_x);
    N_x = N_x - x;
}

void Foam::Structure::dNurbsdCoeff
(
    const gsMatrix<scalar>& s,
    const gsNurbs<scalar>& nurbs,
    List<gsMatrix<scalar>>& dNdC
)
{
    label parameterNbr = s.cols();
    label coeffNbr = nurbs.coefs().rows();
    
    if(s.rows()!=1)
        FatalErrorInFunction<<"s row dimension must be 1"<<exit(FatalError);
    if(nurbs.targetDim()!=3)
        FatalErrorInFunction<<"Mismatch target dimension"<<exit(FatalError);
    
    if(dNdC.size()!=coeffNbr)
        dNdC.resize(coeffNbr);
    
    List<gsNurbs<scalar>> dNurbsdC(coeffNbr);
    for(label coeffI=0; coeffI<coeffNbr; coeffI++)
    {
        dNurbsdC[coeffI] = nurbs;
        gsMatrix<scalar>& dNurbsCoefs = dNurbsdC[coeffI].coefs();
        for(label coeffITemp=0; coeffITemp<coeffNbr; coeffITemp++)
        {
            if(coeffITemp==coeffI)
            {
                dNurbsCoefs(coeffITemp,0)=1;
                dNurbsCoefs(coeffITemp,1)=1;
                dNurbsCoefs(coeffITemp,2)=1;
            }
            else
            {
                dNurbsCoefs(coeffITemp,0)=0;
                dNurbsCoefs(coeffITemp,1)=0;
                dNurbsCoefs(coeffITemp,2)=0;
            }
        }
    }
    
    for(label coeffI=0; coeffI<coeffNbr; coeffI++)
        if(dNdC[coeffI].rows()!=3 || dNdC[coeffI].cols()!=parameterNbr)
            dNdC[coeffI] = gsMatrix<scalar>(3,parameterNbr);
    
    for(label coeffI=0; coeffI<coeffNbr; coeffI++)
    {
        dNurbsdC[coeffI].eval_into(s,dNdC[coeffI]);       
    }
}

void Foam::Structure::intNurbsMinusX
(
    const gsMatrix<scalar>& x,
    const gsMatrix<scalar>& s,
    const gsNurbs<scalar>& nurbs,
    vector& intNmX
)
{
    gsMatrix<scalar> N_x(3,x.cols());
    nurbsMinusX(x,s,nurbs,N_x);
    intNmX = vector(0,0,0);
    
    List<scalar> deltaS(s.cols());
    for(label col=0; col<s.cols(); col++)
    {
        if(col==0)
            deltaS[col] = 0.5*(s(0,0)+s(0,1))-s(0,0);
        else if(col==s.cols()-1)
            deltaS[col] = s(0,col)-0.5*(s(0,col)+s(0,col-1));
        else
            deltaS[col] = 0.5*(s(0,col+1)+s(0,col))-0.5*(s(0,col)+s(0,col-1));;
    }
    
    for(label col=0; col<x.cols(); col++)
    {
        intNmX[0] += std::abs(N_x(0,col))*deltaS[col];
        intNmX[1] += std::abs(N_x(1,col))*deltaS[col];
        intNmX[2] += std::abs(N_x(2,col))*deltaS[col];
    }
}

void Foam::Structure::gradIntNurbsMinX
(
    const gsMatrix<scalar>& x,
    const gsMatrix<scalar>& s,
    const gsNurbs<scalar>& nurbs,
    gsMatrix<scalar>& N_x,
    List<gsMatrix<scalar>>& dNdC,
    gsMatrix<scalar>& gradIntN_x
)
{
    nurbsMinusX(x,s,nurbs,N_x);
    dNurbsdCoeff(s,nurbs,dNdC);
    
    for(label row=0; row<N_x.rows(); row++)
    {
        for(label col=0; col<N_x.cols(); col++)
        {
            if(N_x(row,col)>0)
                N_x(row,col) = 1;
            else if(N_x(row,col)<0)
                N_x(row,col) = -1;
            else
                N_x(row,col) = 0;
        }
    }
    
    if(gradIntN_x.rows()!=dNdC.size() || gradIntN_x.cols()!=3)
        gradIntN_x = gsMatrix<scalar>(dNdC.size(),3);
    
    List<scalar> deltaS(s.cols());
    for(label col=0; col<s.cols(); col++)
    {
        if(col==0)
            deltaS[col] = 0.5*(s(0,0)+s(0,1))-s(0,0);
        else if(col==s.cols()-1)
            deltaS[col] = s(0,col)-0.5*(s(0,col)+s(0,col-1));
        else
            deltaS[col] = 0.5*(s(0,col+1)+s(0,col))-0.5*(s(0,col)+s(0,col-1));;
    }
    
    for(label coeffI=0; coeffI<dNdC.size(); coeffI++)
    {
        const gsMatrix<scalar>& dNdCi = dNdC[coeffI];
        vector gradCoeffi(0,0,0);
        for(label dim=0; dim<dNdCi.rows(); dim++)
        {
            for(label paraInd=0; paraInd<dNdCi.cols(); paraInd++)
                gradCoeffi[dim] += N_x(dim,paraInd)*dNdCi(dim,paraInd)*deltaS[paraInd];
            gradIntN_x(coeffI,dim) = gradCoeffi[dim];
        }
    }
}

void Foam::Structure::minGradStep
(
    gsNurbs<scalar>& nurbs,
    const gsMatrix<scalar>& gradient,
    vector stepsize
)
{
    if(nurbs.coefs().rows()!=gradient.rows())
        FatalErrorInFunction<<"Error"<<exit(FatalError);
    if(nurbs.coefs().cols()!=gradient.cols())
        FatalErrorInFunction<<"Error"<<exit(FatalError);
    if(nurbs.coefs().cols()!=3)
        FatalErrorInFunction<<"Error"<<exit(FatalError);
    
    for(label row=0; row<nurbs.coefs().rows(); row++)
    {
        for(label col=0; col<nurbs.coefs().cols(); col++)
        {
            nurbs.coefs()(row,col) = nurbs.coefs()(row,col) - stepsize[col]*gradient(row,col);
        }
    }
}

void Foam::Structure::linearLeastSquare
(
    const List<vector>& points,
    const List<scalar>& parameters,
    const gsNurbs<scalar>& nurbs,
    gsMatrix<scalar>& fittedCoeffs
)
{
    gsMatrix<scalar> s(1,points.size());
    for(label pntInd=0; pntInd<points.size(); pntInd++)
    {
        s(0,pntInd) = parameters[pntInd];
    }
    
    List<gsMatrix<scalar>> R(nurbs.coefs().rows());
    for(label coeffI=0; coeffI<R.size(); coeffI++)
    {
        gsNurbs<scalar> reducedNurbs = nurbs;
        gsMatrix<scalar>& reducedNurbsCoefs = reducedNurbs.coefs();
        for(label coeffITemp=0; coeffITemp<R.size(); coeffITemp++)
        {
            if(coeffITemp==coeffI)
            {
                reducedNurbsCoefs(coeffITemp,0)=1;
                reducedNurbsCoefs(coeffITemp,1)=1;
                reducedNurbsCoefs(coeffITemp,2)=1;
            }
            else
            {
                reducedNurbsCoefs(coeffITemp,0)=0;
                reducedNurbsCoefs(coeffITemp,1)=0;
                reducedNurbsCoefs(coeffITemp,2)=0;
            }
        }
        reducedNurbs.eval_into(s,R[coeffI]);
    }
        
    FixedList<CSR_Matrix_par,3> A;
    FixedList<Vector_par,3> x;
    FixedList<Vector_par,3> c;
    for(label dim=0; dim<3; dim++)
    {
        A[dim] = CSR_Matrix_par(points.size(),0,points.size(),nurbs.coefs().rows(),false);
        x[dim] = Vector_par(points.size(),0,points.size(),false);
        
        for(label parI=0; parI<points.size(); parI++)
        {
            List<scalar> A_row(R.size());
            for(label coeffI=0; coeffI<R.size(); coeffI++)
            {
                A_row[coeffI] = R[coeffI](dim,parI);
            }
            A[dim].addRow(A_row);
            x[dim][parI] = points[parI][dim];
        }
        
        CSR_Matrix_par transpA = A[dim].transpose();
        CSR_Matrix_par transpA_A = transpA*A[dim];
        Vector_par transpA_x = transpA*x[dim];
        
        if(transpA_x.norm2()<1e-30)
        {
            c[dim] = Vector_par(R.size(),0,R.size(),false);
            c[dim].fill(0);
        }
        else
        {        
            BiCGSTAB solver(transpA_A);
            c[dim] = solver.solve(transpA_x);
        }
    }   
    
    fittedCoeffs = nurbs.coefs();
    for(label dim=0; dim<c.size(); dim++)
    {
        for(label coeffI=0; coeffI<c[dim].getGlobalSize(); coeffI++)
        {
            fittedCoeffs(coeffI,dim) = c[dim][coeffI];
        }
    }
}

void Foam::Structure::fitNurbsCoeffsToPoints
(
    const List<vector>& points,
    const List<scalar>& parameters,
    const gsNurbs<scalar>& nurbs,
    gsMatrix<scalar>& fittedCoeffs,
    scalar epsilon
)
{
    if(points.size()!=parameters.size())
        FatalErrorInFunction<<"Mismatch in points and parameter size"<<exit(FatalError);
    
    linearLeastSquare(points,parameters,nurbs,fittedCoeffs);
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
    
    rodEvalBuffer[rodNumber].clear();
    rodDerivEvalBuffer[rodNumber].clear();
    rodDeriv2EvalBuffer[rodNumber].clear();
}

const gsNurbs<Foam::scalar>& Foam::Structure::getDeformation
(
    label rodNumber
) const
{
    if(rodNumber<0 || rodNumber>=getNumberRods())
        FatalErrorInFunction<<"Invalid rodNumber"<<exit(FatalError);
    return Rods[rodNumber]->m_Def;
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
    label rodNumber,
    scalar parameter,
    vector& d1,
    vector& d2,
    vector& d3,
    vector& r
)
{
    std::unordered_map<scalar,std::array<vector,4>>& oneRodEvalBuffer = rodEvalBuffer[rodNumber];
    auto iter = oneRodEvalBuffer.find(parameter);
    if(iter!=oneRodEvalBuffer.end())
    {
        const std::array<vector,4>& bufferedResult = iter->second;
        d1 = bufferedResult[0];
        d2 = bufferedResult[1];
        d3 = bufferedResult[2];
        r  = bufferedResult[3];
    }
    else
    {
        rodEval(Rods[rodNumber],parameter,d1,d2,d3,r);
        oneRodEvalBuffer.insert({parameter,{d1,d2,d3,r}});
    }
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

void Foam::Structure::rodDerivEval
(
    label rodNumber,
    scalar parameter,
    vector& dd1dp,
    vector& dd2dp,
    vector& dd3dp,
    vector& drdp
)
{
    std::unordered_map<scalar,std::array<vector,4>>& oneRodDerivEvalBuffer = rodDerivEvalBuffer[rodNumber];
    auto iter = oneRodDerivEvalBuffer.find(parameter);
    if(iter!=oneRodDerivEvalBuffer.end())
    {
        const std::array<vector,4>& bufferedResult = iter->second;
        dd1dp = bufferedResult[0];
        dd2dp = bufferedResult[1];
        dd3dp = bufferedResult[2];
        drdp  = bufferedResult[3];
    }
    else
    {
        rodDerivEval(Rods[rodNumber],parameter,dd1dp,dd2dp,dd3dp,drdp);
        oneRodDerivEvalBuffer.insert({parameter,{dd1dp,dd2dp,dd3dp,drdp}});
    }
}

void Foam::Structure::rodDerivEval
(
    const ActiveRodMesh::rodCosserat* rod,
    scalar parameter,
    vector& dd1dp,
    vector& dd2dp,
    vector& dd3dp,
    vector& drdp
)
{
    //rodDerivEval(rod->m_Curve,rod->m_Def,rod->m_Rot,parameter,dd1dp,dd2dp,dd3dp,drdp);
    
    const gsNurbs<scalar>& curve = rod->m_Curve;
    const gsNurbs<scalar>& def = rod->m_Def;
    const gsNurbs<scalar>& rot = rod->m_Rot;  
    
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
        
    gsMatrix<scalar> basePnt;
    curve.deriv_into(parMat,basePnt);

    gsMatrix<scalar> defPnt;
    def.deriv_into(parMat,defPnt);
    
    gsMatrix<scalar> pnt = basePnt+defPnt;    
    drdp = vector(pnt(0,0),pnt(1,0),pnt(2,0));

    gsMatrix<scalar> rotQuat;
    rot.eval_into(parMat,rotQuat);
    Quaternion q(rotQuat);
    
    gsMatrix<scalar> drotQuatdp = evalNurbsDeriv(rot,parameter);
    Quaternion dqdp(drotQuatdp);
    
    Rotation dRdp = Rotation::compute_dRdX(dqdp,q);
    
    dd1dp = dRdp.get_d1();
    dd2dp = dRdp.get_d2();
    dd3dp = dRdp.get_d3();
}

void Foam::Structure::rodDeriv2Eval
(
    label rodNumber,
    scalar parameter,
    vector& d2d1dp,
    vector& d2d2dp,
    vector& d2d3dp,
    vector& d2rdp
)
{
    std::unordered_map<scalar,std::array<vector,4>>& oneRodDeriv2EvalBuffer = rodDeriv2EvalBuffer[rodNumber];
    auto iter = oneRodDeriv2EvalBuffer.find(parameter);
    if(iter!=oneRodDeriv2EvalBuffer.end())
    {
        const std::array<vector,4>& bufferedResult = iter->second;
        d2d1dp = bufferedResult[0];
        d2d2dp = bufferedResult[1];
        d2d3dp = bufferedResult[2];
        d2rdp  = bufferedResult[3];
    }
    else
    {
        rodDeriv2Eval(Rods[rodNumber],parameter,d2d1dp,d2d2dp,d2d3dp,d2rdp);
        oneRodDeriv2EvalBuffer.insert({parameter,{d2d1dp,d2d2dp,d2d3dp,d2rdp}});
    }
}

void Foam::Structure::rodDeriv2Eval
(
    const ActiveRodMesh::rodCosserat* rod,
    scalar parameter,
    vector& d2d1dp,
    vector& d2d2dp,
    vector& d2d3dp,
    vector& d2rdp
)
{
    //rodDeriv2Eval(rod->m_Curve,rod->m_Def,rod->m_Rot,parameter,d2d1dp,d2d2dp,d2d3dp,d2rdp);
    
    
    const gsNurbs<scalar>& curve = rod->m_Curve;
    const gsNurbs<scalar>& def = rod->m_Def;
    const gsNurbs<scalar>& rot = rod->m_Rot;  
    
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
        
    gsMatrix<scalar> basePnt = evalNurbsDeriv2(curve,parameter);
    

    gsMatrix<scalar> defPnt = evalNurbsDeriv2(def,parameter);
    
    
    gsMatrix<scalar> pnt = basePnt+defPnt;    
    d2rdp = vector(pnt(0,0),pnt(1,0),pnt(2,0));

    gsMatrix<scalar> rotQuat;
    rot.eval_into(parMat,rotQuat);
    Quaternion q(rotQuat);
    
    
    gsMatrix<scalar> drotQuatdp = evalNurbsDeriv(rot,parameter);
    Quaternion dqdp(drotQuatdp);
    
    
    gsMatrix<scalar> d2rotQuatdp = evalNurbsDeriv2(rot,parameter);
    Quaternion d2qdp(d2rotQuatdp);
    
    
    Rotation d2Rdp = Rotation::compute_d2RdX(d2qdp,dqdp,q);
    
    
    d2d1dp = d2Rdp.get_d1();
    d2d2dp = d2Rdp.get_d2();
    d2d3dp = d2Rdp.get_d3();
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

void Foam::Structure::rodDerivEval
(
    const gsNurbs<scalar>& curve,
    const gsNurbs<scalar>& def,
    const gsNurbs<scalar>& rot,
    scalar parameter,
    vector& dd1dp,
    vector& dd2dp,
    vector& dd3dp,
    vector& drdp
)
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
        
    gsMatrix<scalar> basePnt;
    curve.deriv_into(parMat,basePnt);

    gsMatrix<scalar> defPnt;
    def.deriv_into(parMat,defPnt);
    
    gsMatrix<scalar> pnt = basePnt+defPnt;    
    drdp = vector(pnt(0,0),pnt(1,0),pnt(2,0));

    gsMatrix<scalar> rotQuat;
    rot.eval_into(parMat,rotQuat);
    Quaternion q(rotQuat);
    
    gsMatrix<scalar> drotQuatdp = evalNurbsDeriv(rot,parameter);
    Quaternion dqdp(drotQuatdp);
    
    Rotation dRdp = Rotation::compute_dRdX(dqdp,q);
    
    dd1dp = dRdp.get_d1();
    dd2dp = dRdp.get_d2();
    dd3dp = dRdp.get_d3();
}

void Foam::Structure::rodDeriv2Eval
(
    const gsNurbs<scalar>& curve,
    const gsNurbs<scalar>& def,
    const gsNurbs<scalar>& rot,
    scalar parameter,
    vector& d2d1dp,
    vector& d2d2dp,
    vector& d2d3dp,
    vector& d2rdp
)
{
    auto t0 = std::chrono::system_clock::now();

    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    auto t1 = std::chrono::system_clock::now();
        
    gsMatrix<scalar> basePnt = evalNurbsDeriv2(curve,parameter);
    
    auto t2 = std::chrono::system_clock::now();

    gsMatrix<scalar> defPnt = evalNurbsDeriv2(def,parameter);
    
    auto t3 = std::chrono::system_clock::now();
    
    gsMatrix<scalar> pnt = basePnt+defPnt;    
    d2rdp = vector(pnt(0,0),pnt(1,0),pnt(2,0));

    gsMatrix<scalar> rotQuat;
    rot.eval_into(parMat,rotQuat);
    Quaternion q(rotQuat);
    
    auto t4 = std::chrono::system_clock::now();
    
    gsMatrix<scalar> drotQuatdp = evalNurbsDeriv(rot,parameter);
    Quaternion dqdp(drotQuatdp);
    
    auto t5 = std::chrono::system_clock::now();
    
    gsMatrix<scalar> d2rotQuatdp = evalNurbsDeriv2(rot,parameter);
    Quaternion d2qdp(d2rotQuatdp);
    
    auto t6 = std::chrono::system_clock::now();
    
    Rotation d2Rdp = Rotation::compute_d2RdX(d2qdp,dqdp,q);
    
    auto t7 = std::chrono::system_clock::now();
    
    d2d1dp = d2Rdp.get_d1();
    d2d2dp = d2Rdp.get_d2();
    d2d3dp = d2Rdp.get_d3();
    auto t8 = std::chrono::system_clock::now();
    //std::cout<<"--------------End---------------------"<<std::endl;
    Info<<"\t\t\t rodDeriv2Eval t0-t1 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count()<<Foam::nl;
    Info<<"\t\t\t rodDeriv2Eval t1-t2 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()<<Foam::nl;
    Info<<"\t\t\t rodDeriv2Eval t2-t3 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count()<<Foam::nl;
    Info<<"\t\t\t rodDeriv2Eval t3-t4 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count()<<Foam::nl;
    Info<<"\t\t\t rodDeriv2Eval t4-t5 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count()<<Foam::nl;
    Info<<"\t\t\t rodDeriv2Eval t5-t6 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t5).count()<<Foam::nl;
    Info<<"\t\t\t rodDeriv2Eval t6-t7 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t7-t6).count()<<Foam::nl;
    Info<<"\t\t\t rodDeriv2Eval t7-t8 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t8-t7).count()<<Foam::nl;
    Info<<"\t\t\t sum t0-t8 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t8-t0).count()<<Foam::nl;
}

Foam::Quaternion Foam::Structure::m_Rot_Eval
(
    label rodNumber,
    scalar parameter
)
{
    return m_Rot_Eval(Rods[rodNumber],parameter);
}

Foam::Quaternion Foam::Structure::m_Rot_Eval
(
    const ActiveRodMesh::rodCosserat* rod,
    scalar parameter
)
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
        
    gsMatrix<scalar> rotPnt;
    rod->m_Rot.eval_into(parMat,rotPnt);
    
    return Quaternion(rotPnt);
}

const gsNurbs<Foam::scalar>& Foam::Structure::getQuaternions
(
    label rodNumber
) const
{
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber given"<<exit(FatalError);
    return Rods[rodNumber]->m_Rot;    
}

Foam::List<Foam::Rotation> Foam::Structure::getRotationCurveCoefs
(
    label rodNumber
) const
{
    Info<<"-------------------------------getRotationCurveCoefs---------------------------------------"<<Foam::endl;
    gsNurbs<scalar> quaternions = getQuaternions(rodNumber);
    gsMatrix<scalar> quaternionCoef = quaternions.coefs();
    std::cout<<quaternionCoef<<std::endl;
    List<Rotation> rotations(quaternionCoef.rows());
    for(label coefNbr=0; coefNbr<quaternionCoef.rows(); coefNbr++)
    {
        gsMatrix<scalar> gsQ(1,4);
        for(label qdim=0; qdim<4; qdim++)
            gsQ(0,qdim) = quaternionCoef(coefNbr,qdim);
        rotations[coefNbr] = Rotation(Quaternion(gsQ));
    }
    return rotations;
}

gsMatrix<Foam::scalar> Foam::Structure::evalNurbs
(
    const gsNurbs<scalar>& nurbs,
    scalar parameter
)
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;        
    gsMatrix<scalar> basePnt;
    nurbs.eval_into(parMat,basePnt);
    return basePnt;
}

gsMatrix<Foam::scalar> Foam::Structure::evalNurbsDeriv
(
    const gsNurbs<scalar>& nurbs,
    scalar parameter
)
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;        
    gsMatrix<scalar> basePnt;
    nurbs.deriv_into(parMat,basePnt);
    return basePnt;
}

gsMatrix<Foam::scalar> Foam::Structure::evalNurbsDeriv2
(
    const gsNurbs<scalar>& nurbs,
    scalar parameter
)
{
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;        
    gsMatrix<scalar> basePnt;
    nurbs.deriv2_into(parMat,basePnt);
    return basePnt;
}

void Foam::Structure::constructCoeffDerivedData()
{
    coeffDerivedCenterline.resize(nR);
    initialRotation.resize(nR);
    coeffDerivedQuaternions.resize(nR);
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        ActiveRodMesh::rodCosserat* rod = Rods[rodNumber];
        rod->m_Rot_end0 = false;
        const gsNurbs<scalar>& curve = rod->m_Curve;
        //std::cout<<"curve.coefs:"<<curve.coefs()<<std::endl;
        label numberCurveCoefs = curve.coefs().rows();
        const gsNurbs<scalar>& rotation = rod->m_Rot;
        //std::cout<<"rotation.coefs:"<<rotation.coefs()<<std::endl;
        label numberRotQuaternions = rotation.coefs().rows();
        gsMatrix<scalar>& dQdRreset = rod->m_Curve_dQdR0;
        dQdRreset = gsMatrix<scalar>(numberRotQuaternions*4,numberCurveCoefs*3);
        dQdRreset.fill(0);
        int errorCode = rod->bishopFrameS_dQdR0();
        //Info<<"-------"<<Foam::endl;
        const gsMatrix<scalar>& dQdR = rod->m_Curve_dQdR0;
        if(errorCode==0)
            FatalErrorInFunction<<"Failure to compute Quaternion coefficient derivatives"<<exit(FatalError);
        
        coeffDerivedCenterline[rodNumber].resize(numberCurveCoeffs(rodNumber));
        initialRotation[rodNumber] = rotation;
        coeffDerivedQuaternions[rodNumber].resize(numberCurveCoeffs(rodNumber));
        for(label coeffNumber=0; coeffNumber<numberCurveCoeffs(rodNumber); coeffNumber++)
        {
            
            //Create centerline coefficient derivative curve
            coeffDerivedCenterline[rodNumber][coeffNumber].resize(3);
            for(label dim=0; dim<3; dim++)
            {
                coeffDerivedCenterline[rodNumber][coeffNumber][dim] = curve;
                gsNurbs<scalar>& curve = coeffDerivedCenterline[rodNumber][coeffNumber][dim];
                gsMatrix<scalar>& coeffs = curve.coefs();
                for(label row=0; row<coeffs.rows(); row++)
                {
                    for(label col=0; col<coeffs.cols(); col++)
                    {
                        if(coeffNumber==row && col==dim)
                            coeffs(row,col) = 1;
                        else
                            coeffs(row,col) = 0;
                    }
                }

                //Create quaternion coefficient derivative curve
                coeffDerivedQuaternions[rodNumber][coeffNumber].resize(3);
                coeffDerivedQuaternions[rodNumber][coeffNumber][dim] = rotation;
                gsNurbs<scalar>& rotation = coeffDerivedQuaternions[rodNumber][coeffNumber][dim];
                gsMatrix<scalar>& qcoeffs = rotation.coefs();
                label coeffRow = coeffNumber*3+dim;
                List<scalar> dq_doneCoeff(dQdR.rows());
                for(label row=0; row<dQdR.rows(); row++)
                    dq_doneCoeff[row] = dQdR(row,coeffRow);

                label listIndex=0;
                for(label rotQCoeff=0; rotQCoeff<qcoeffs.rows(); rotQCoeff++)
                {
                    for(label rotQDim=0; rotQDim<qcoeffs.cols(); rotQDim++)
                    {
                        qcoeffs(rotQCoeff,rotQDim) = dq_doneCoeff[listIndex];
                        listIndex++;
                    }
                }
            }
        }
    }
    constructedCoeffDerivedData = true;    
}

Foam::label Foam::Structure::numberCurveCoeffs
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

Foam::vector Foam::Structure::get_drdC
(
    label rodNumber,
    label derivCoeffNumber,
    label derivDimension,
    scalar parameter
) const
{
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber given"<<exit(FatalError);
    if(derivCoeffNumber<0 || derivCoeffNumber>=static_cast<label>(coeffDerivedCenterline[rodNumber].size()))
        FatalErrorInFunction<<"Invalid derivCoeffNumber!"<<exit(FatalError);
    if(derivDimension<0 || derivDimension>=3)
        FatalErrorInFunction<<"Invalid derivDimension!"<<exit(FatalError);
    
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    // Derivative of rod center line in respect to coefficient
    const gsNurbs<scalar>& oneCoeffDerivCenterline = coeffDerivedCenterline[rodNumber][derivCoeffNumber][derivDimension];
    gsMatrix<scalar> coeffDerivCenterEval;
    oneCoeffDerivCenterline.eval_into(parMat,coeffDerivCenterEval);    
    return vector(coeffDerivCenterEval(0,0),coeffDerivCenterEval(1,0),coeffDerivCenterEval(2,0));
}

Foam::Quaternion Foam::Structure::get_iniQuaternions
(
    label rodNumber,
    scalar parameter
) const
{
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber given"<<exit(FatalError);
    
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    const gsNurbs<scalar>& iniQuaternions = initialRotation[rodNumber];
    gsMatrix<scalar> iniQuaternionEval;
    iniQuaternions.eval_into(parMat,iniQuaternionEval);
    
    return Quaternion(iniQuaternionEval);
}

Foam::Quaternion Foam::Structure::get_Quaternions
(
    label rodNumber,
    scalar parameter
) const
{
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber given"<<exit(FatalError);
    
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    const gsNurbs<scalar>& totalQuaternions = Rods[rodNumber]->m_Rot;
    gsMatrix<scalar> totalQuaternionEval;
    totalQuaternions.eval_into(parMat,totalQuaternionEval);
    
    return Quaternion(totalQuaternionEval);
}

Foam::Quaternion Foam::Structure::get_coeffDerivQuaternions
(
    label rodNumber,
    label derivCoeffNumber,
    label derivDimension,
    scalar parameter
) const
{
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber given"<<exit(FatalError);
    if(derivCoeffNumber<0 || derivCoeffNumber>=static_cast<label>(coeffDerivedCenterline[rodNumber].size()))
        FatalErrorInFunction<<"Invalid derivCoeffNumber!"<<exit(FatalError);
    if(derivDimension<0 || derivDimension>=3)
        FatalErrorInFunction<<"Invalid derivDimension!"<<exit(FatalError);
    
    gsMatrix<scalar> parMat(1,1);
    parMat.at(0) = parameter;
    
    const gsNurbs<scalar>& oneCoeffDerivQuaternions = coeffDerivedQuaternions[rodNumber][derivCoeffNumber][derivDimension];
    gsMatrix<scalar> coeffDerivQuaternionEval;
    oneCoeffDerivQuaternions.eval_into(parMat,coeffDerivQuaternionEval);

    return Quaternion(coeffDerivQuaternionEval);
}

const gsNurbs<Foam::scalar>& Foam::Structure::getCoeffDerivedQuaternions
(
    label rodNumber,
    label derivCoeffNumber,
    label derivDimension
) const
{
    if(!constructedCoeffDerivedData)
        FatalErrorInFunction<<"Data for deriv coeff data not given!"<<exit(FatalError);
        
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber given"<<exit(FatalError);
    if(derivCoeffNumber<0 || derivCoeffNumber>=static_cast<label>(coeffDerivedQuaternions[rodNumber].size()))
        FatalErrorInFunction<<"Invalid derivCoeffNumber!"<<exit(FatalError);
    if(derivDimension<0 || derivDimension>=3)
        FatalErrorInFunction<<"Invalid derivDimension!"<<exit(FatalError);
    
    // Derivative of rod center line in respect to coefficient
    return coeffDerivedQuaternions[rodNumber][derivCoeffNumber][derivDimension];
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
    
    // Derivative of rod center line in respect to coefficient    
    rdC = get_drdC(rodNumber,derivCoeffNumber,derivDimension,parameter);
    
    //Derivative of rod coordinate axis in respect to coefficien
    const Quaternion dQdC = get_coeffDerivQuaternions(rodNumber,derivCoeffNumber,derivDimension,parameter);
    
    const Quaternion totalQ = get_Quaternions(rodNumber,parameter);
    const Quaternion iniQ = get_iniQuaternions(rodNumber,parameter);
    
    const Quaternion defQ = totalQ/iniQ;

    const Quaternion defQ_inidQdC = defQ*dQdC;
    
    const Rotation dRdC = Rotation::compute_dRdX(defQ_inidQdC,totalQ);
    
    d1dC = dRdC.get_d1();
    d2dC = dRdC.get_d2();
    d3dC = dRdC.get_d3();    
}

Foam::Quaternion Foam::Structure::m_Rot_Eval_Deriv
(
    label rodNumber,
    label derivCoeffNumber,
    label derivDimension,
    scalar parameter
)
{
    if(!constructedCoeffDerivedData)
        FatalErrorInFunction<<"Data for deriv coeff data not given!"<<exit(FatalError);
        
    const Quaternion dQdC = get_coeffDerivQuaternions(rodNumber,derivCoeffNumber,derivDimension,parameter);
    
    const  Quaternion totalQ = get_Quaternions(rodNumber,parameter);
    const Quaternion iniQ = get_iniQuaternions(rodNumber,parameter);
    
    const Quaternion defQ = totalQ*iniQ;

    const Quaternion defQ_inidQdC = defQ*dQdC;
    
    return defQ_inidQdC;
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
        rodCoeffs.resize(numberCurveCoeffs(rodNumber));
        const gsMatrix<scalar>& curveCoeffs = rodsList[rodNumber].coefs();
        for(label coeffNumber=0; coeffNumber<numberCurveCoeffs(rodNumber); coeffNumber++)
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
    Info<<"--------------- Re-set curve coeffs ----------------"<<Foam::endl;
    
    if(coeffs.size()!=nR)
        FatalErrorInFunction<<"Mismatch in rod number to coeffs!"<<exit(FatalError);
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        const List<vector>& rodCoeffs = coeffs[rodNumber];
        gsMatrix<scalar>& curveCoeffs = rodsList[rodNumber].coefs();
        for(label coeffNumber=0; coeffNumber<numberCurveCoeffs(rodNumber); coeffNumber++)
        {
            curveCoeffs(coeffNumber,0) = rodCoeffs[coeffNumber][0];
            curveCoeffs(coeffNumber,1) = rodCoeffs[coeffNumber][1];
            curveCoeffs(coeffNumber,2) = rodCoeffs[coeffNumber][2];
        }
        rodEvalBuffer[rodNumber].clear();
        rodDerivEvalBuffer[rodNumber].clear();
        rodDeriv2EvalBuffer[rodNumber].clear();
    }
    setupActiveRodMesh();
}

Foam::scalar Foam::Structure::getCurveCoeff
(
    label rodNumber,
    label derivCoeffNumber,
    label dimension
) const
{
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber"<<exit(FatalError);
    ActiveRodMesh::rodCosserat* rod = Rods[rodNumber];
    gsNurbs<scalar>& curve = rod->m_Curve;
    gsMatrix<scalar>& coeffs = curve.coefs();
    if(derivCoeffNumber<0 || derivCoeffNumber>=coeffs.rows())
        FatalErrorInFunction<<"Invalid derivCoeffNumber"<<exit(FatalError);
    if(dimension<0 || dimension>=3)
        FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    return coeffs(derivCoeffNumber,dimension);
}

void Foam::Structure::setCurveCoeff
(
    label rodNumber,
    label derivCoeffNumber,
    label dimension,
    scalar value
)
{
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber"<<exit(FatalError);
    List<List<vector>> rodCoefs;
    getCurveCoeffs(rodCoefs);
    if(derivCoeffNumber<0 || derivCoeffNumber>=rodCoefs[rodNumber].size())
        FatalErrorInFunction<<"Invalid derivCoeffNumber"<<exit(FatalError);
    if(dimension<0 || dimension>=3)
        FatalErrorInFunction<<"Invalid dimension"<<exit(FatalError);
    rodCoefs[rodNumber][derivCoeffNumber][dimension] = value;
    setCurveCoeffs(rodCoefs);  
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

void Foam::Structure::createParallelTopology()
{
    List<Pair<vector>> boundingBoxes(Pstream::nProcs());
    if(meshBoundingBox.isEmpty())
        FatalErrorInFunction<<"meshBoundingBox can not be empty"<<exit(FatalError);
    boundingBoxes[Pstream::myProcNo()] = meshBoundingBox.getBB();
    Pstream::gatherList(boundingBoxes);
    Pstream::scatterList(boundingBoxes);
    
    for(label proc=0; proc<Pstream::nProcs(); proc++)
    {
        BoundingBox procBB(boundingBoxes[proc].first(),boundingBoxes[proc].second());
        for(label nei=proc+1; nei<Pstream::nProcs(); nei++)
        {
            BoundingBox neiBB(boundingBoxes[nei].first(),boundingBoxes[nei].second());
            scalar bbSize = std::max(procBB.innerSize(),neiBB.innerSize());
            if(procBB.distance(neiBB)<=bbSize)
            {
                label comm = Pstream::allocateCommunicator(Pstream::masterNo(),{proc,nei});
                proc_ProcToComm[{proc,nei}] = comm;
                proc_ProcToComm[{nei,proc}] = comm;
            }
        }
    }
    
    std::set<label> orderedComms;
    myProcExistingData.setSize(Pstream::nProcs(),false);
    myProcExistingData[Pstream::myProcNo()] = true;
    for(auto iter=proc_ProcToComm.begin(); iter!=proc_ProcToComm.end(); iter++)
    {
        const Pair<label>& pairProcs = iter->first;
        const label comm = iter->second;
        orderedComms.insert(comm);
        if(pairProcs.first()==Pstream::myProcNo())
        {
            neighbourProcesses.insert(pairProcs.second());
            myProcExistingData[pairProcs.second()] = true;
        }
    }
    for(auto iter=orderedComms.begin(); iter!=orderedComms.end(); iter++)
        this->orderedComms.append(*iter);
}

void Foam::Structure::computePatchFaceToCellMap()
{
    patchFaceToNeighCellMap.clear();
    
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();
    const polyBoundaryMesh& boundaries = mesh.boundaryMesh();
    
    //[proc] -> [...] -> {neighbour,patchFaceInd,cellInd}
    List<DynamicList<Tuple3<label,label,label>>> patchFaceToOwnCell(Pstream::nProcs());
    //[...] -> {neighbour,patchFaceInd,cellInd}
    DynamicList<Tuple3<label,label,label>>& procPatchFaceToOwnCell = patchFaceToOwnCell[Pstream::myProcNo()];
    
    std::unordered_map<label,const processorPolyPatch*> procToPatchPtr;
    
    for(label patchIndex=0; patchIndex<boundaries.size(); patchIndex++)
    {
        const polyPatch& patch = boundaries[patchIndex];
        if(isA<processorPolyPatch>(patch))
        {
            const processorPolyPatch* pPP = dynamic_cast<const processorPolyPatch*>(&patch);
            label neighborProcess = pPP->neighbProcNo();
            label patchStartFace = pPP->start();
            label patchSizeFaces = pPP->size();
            procToPatchPtr[neighborProcess] = pPP;
            
            for(label locFaceInd=0; locFaceInd<patchSizeFaces; locFaceInd++)
            {
                label faceInd = locFaceInd+patchStartFace;
                if(faceInd<neighbour.size())
                    FatalErrorInFunction<<"Boundary face has a neighbour"<<exit(FatalError);
                
                label cellInd = owner[faceInd];
                procPatchFaceToOwnCell.append({neighborProcess,locFaceInd,cellInd});
            }
        }
    }
    
    exchangeBetweenAll(patchFaceToOwnCell);
    
    for(label proc=0; proc<Pstream::nProcs(); proc++)
    {
        const DynamicList<Tuple3<label,label,label>>& procPatchFaceToCell = patchFaceToOwnCell[proc];
        for(const Tuple3<label,label,label>& patchData : procPatchFaceToCell)
        {
            label neighProc = patchData.first();
            label patchFaceInd = patchData.second();
            label cellInd = patchData.third();
            if(neighProc==Pstream::myProcNo())
            {
                auto iterPatch = procToPatchPtr.find(proc);
                if(iterPatch==procToPatchPtr.end())
                    FatalErrorInFunction<<"Missing patch pointer"<<exit(FatalError);
                const processorPolyPatch* pPP = iterPatch->second;
                label neighborProcess = pPP->neighbProcNo();
                if(neighborProcess!=proc)
                {
                    Pout<<"proc:"<<proc<<Foam::nl;
                    Pout<<"patchData:"<<patchData<<Foam::nl;
                    Pout<<"neighborProcess:"<<neighborProcess<<Foam::nl;
                    FatalErrorInFunction<<"Proc mismatch"<<exit(FatalError);
                }
                label patchStartFace = pPP->start();
                label patchSizeFaces = pPP->size();
                if(patchFaceInd>=patchSizeFaces)
                    FatalErrorInFunction<<"Out of range patchFaceInd"<<exit(FatalError);
                label faceInd = patchStartFace+patchFaceInd;                
                auto iter = patchFaceToNeighCellMap.find(faceInd);
                if(iter!=patchFaceToNeighCellMap.end())
                {
                    Pout<<"proc:"<<proc<<Foam::nl;
                    Pout<<"patchData:"<<patchData<<Foam::nl;
                    Pout<<"iter->first:"<<iter->first<<Foam::nl;
                    Pout<<"iter->second:"<<iter->second<<Foam::nl;
                    FatalErrorInFunction<<"Duplicate entry for face"<<exit(FatalError);
                }
                patchFaceToNeighCellMap[faceInd] = {proc,cellInd};
            }
        }
    }
    
    /*
    for(auto iter=procToPatchPtr.begin(); iter!=procToPatchPtr.end(); iter++)
    {
        label proc = iter->first;
        const processorPolyPatch* pPP = iter->second;
        label neighborProcess = pPP->neighbProcNo();
        label patchStartFace = pPP->start();
        label patchSizeFaces = pPP->size();
        
        for(label locInd=0; locInd<patchSizeFaces; locInd++)
        {
            label faceInd = locInd+patchStartFace;
            auto iterFace = patchFaceToNeighCellMap.find(faceInd);
            if(iterFace==patchFaceToNeighCellMap.end())
                FatalErrorInFunction<<"Missing entry"<<exit(FatalError);
            if(iterFace->second.first()!=neighborProcess)
                FatalErrorInFunction<<"Mismatch proc"<<exit(FatalError);
        }
    }
    
    Barrier(true);
    */
}

void Foam::Structure::generateMeshGraph()
{
    /*
    Barrier(false);
    auto t2 = std::chrono::high_resolution_clock::now();
    */
    
    computePatchFaceToCellMap();
    
    /*
    Barrier(false);
    auto t1 = std::chrono::high_resolution_clock::now();
    Info<<"computePatchFaceToCellMap took "<<std::chrono::duration_cast<std::chrono::microseconds>(t1-t2).count()<<" microseconds"<<Foam::nl;    
    */
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const labelList& owners = mesh.owner();
    const labelList& neighbours = mesh.neighbour();
    const pointField& points = mesh.points();
    
    globalMeshGraph.clear();
    globalMeshGraph.setSize(Pstream::nProcs());
    localMeshGraph = &(globalMeshGraph[Pstream::myProcNo()]);
    localMeshGraph->setSize(cells.size());
    List<List<Pair<label>>>& meshGraph = *localMeshGraph;
    
    globalMeshCentres.clear();
    globalMeshCentres.setSize(Pstream::nProcs());
    globalMeshCentres[Pstream::myProcNo()].setSize(cells.size());
    
    globalMeshMagn.clear();
    globalMeshMagn.setSize(Pstream::nProcs());
    globalMeshMagn[Pstream::myProcNo()].setSize(cells.size());
    
    /*
    Barrier(false);
    t2 = std::chrono::high_resolution_clock::now();
    Info<<"Graph setup took "<<std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()<<" microseconds"<<Foam::nl;
    */
    
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
            else if(patchFaceToNeighCellMap.find(faceInd)!=patchFaceToNeighCellMap.end())
            {
                auto iter = patchFaceToNeighCellMap.find(faceInd);
                const Pair<label>& neighCell = iter->second;
                neighData = neighCell;
                if(neighData.first()==Pstream::myProcNo())
                    FatalErrorInFunction<<"Error"<<exit(FatalError);
            }
            meshGraph[cellInd][cellFaceInd] = neighData;
        }
        globalMeshCentres[Pstream::myProcNo()][cellInd] = oneCell.centre(points,faces);
        globalMeshMagn[Pstream::myProcNo()][cellInd] = oneCell.mag(points,faces);
    }
    
    /*
    Barrier(false);
    t1 = std::chrono::high_resolution_clock::now();
    Info<<"Graph collection took "<<std::chrono::duration_cast<std::chrono::microseconds>(t1-t2).count()<<" microseconds"<<Foam::nl;
    */
        
    // Flatten data to transfer
    List<List<label>> cellNeiCount(Pstream::nProcs());
    cellNeiCount[Pstream::myProcNo()].setSize(meshGraph.size());
    List<DynamicList<Pair<label>>> serializedMeshGraph(Pstream::nProcs());
    for(label cellInd=0; cellInd<meshGraph.size(); cellInd++)
    {
        cellNeiCount[Pstream::myProcNo()][cellInd] = meshGraph[cellInd].size();
        for(label cellNeiInd=0; cellNeiInd<meshGraph[cellInd].size(); cellNeiInd++)
        {
            serializedMeshGraph[Pstream::myProcNo()].append(meshGraph[cellInd][cellNeiInd]);
        }
    }
    
    /*
    Barrier(false);
    t2 = std::chrono::high_resolution_clock::now();
    Info<<"Flatten took "<<std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()<<" microseconds"<<Foam::nl;
    */
    
    // Transfer
    exchangeBetweenAll(cellNeiCount);
    /*
    Barrier(false);
    t1 = std::chrono::high_resolution_clock::now();
    Info<<"cellNeiCount took "<<std::chrono::duration_cast<std::chrono::microseconds>(t1-t2).count()<<" microseconds"<<Foam::nl;
    */
    
    exchangeBetweenAll(serializedMeshGraph);
    /*
    Barrier(false);
    t2 = std::chrono::high_resolution_clock::now();
    Info<<"serializedMeshGraph took "<<std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()<<" microseconds"<<Foam::nl;
    */
    
    exchangeBetweenAll(globalMeshCentres);
    /*
    Barrier(false);
    t1 = std::chrono::high_resolution_clock::now();
    Info<<"globalMeshCentres took "<<std::chrono::duration_cast<std::chrono::microseconds>(t1-t2).count()<<" microseconds"<<Foam::nl;
    */
    
    exchangeBetweenAll(globalMeshMagn);
    /*
    Barrier(false);
    t2 = std::chrono::high_resolution_clock::now();
    Info<<"globalMeshMagn took "<<std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()<<" microseconds"<<Foam::nl;
    */
           
    // Restructure
    for(label proc=0; proc<Pstream::nProcs(); proc++)
    {
        if(proc!=Pstream::myProcNo())
        {
            const List<label>& procCellNeiCount = cellNeiCount[proc];
            const DynamicList<Pair<label>>& procSerializedMeshGraph = serializedMeshGraph[proc];
            
            List<List<Pair<label>>>& procGlobalMeshGraph = globalMeshGraph[proc];
            procGlobalMeshGraph.setSize(procCellNeiCount.size());
            
            label index=0;
            for(label cellInd=0; cellInd<procCellNeiCount.size(); cellInd++)
            {
                procGlobalMeshGraph[cellInd].setSize(procCellNeiCount[cellInd]);
                for(label cellNeiInd=0; cellNeiInd<procCellNeiCount[cellInd]; cellNeiInd++)
                {
                    procGlobalMeshGraph[cellInd][cellNeiInd] = procSerializedMeshGraph[index];
                    index++;
                }
            }
        }
    }
    
    /*
    Barrier(false);
    t1 = std::chrono::high_resolution_clock::now();
    Info<<"Restructure took "<<std::chrono::duration_cast<std::chrono::microseconds>(t1-t2).count()<<" microseconds"<<Foam::nl;
    */
}

const Foam::List<Foam::List<Foam::Pair<Foam::label>>>& Foam::Structure::getMeshGraph
(
    label proc
) const
{
    if(proc<0 || proc>=Pstream::nProcs())
        FatalErrorInFunction<<"Out of bounds proc number:"<<proc<<exit(FatalError);
    return globalMeshGraph[proc];
}

void Foam::Structure::generateMeshHaloData
(
    label iterations
)
{   
    if(iterations<1)
        FatalErrorInFunction<<"Must be at least one iterations"<<exit(FatalError);
    
    generateMeshGraph();
    
    //auto t1 = std::chrono::high_resolution_clock::now();
    
    globalHaloCellList_Sorted.clear();
    globalHaloCellList_Sorted.resize(Pstream::nProcs());
    globalHaloCellToIndexMap.clear();
    globalHaloCellToIndexMap.resize(Pstream::nProcs());
    for(label proc=0; proc<Pstream::nProcs(); proc++)
    {
        std::unordered_map<label,bool> procBoundaryCells;
        const List<List<Pair<label>>>& procMeshGraph = globalMeshGraph[proc];
        for(label locCellInd=0; locCellInd<procMeshGraph.size(); locCellInd++)
        {
            const List<Pair<label>>& cellNeis = procMeshGraph[locCellInd];
            bool border=false;
            for(const Pair<label>& nei : cellNeis)
            {
                if(nei.first()!=proc && nei.first()!=-1)
                    border=true;
            }           
            if(border)
                procBoundaryCells.insert({locCellInd,false});
        }
        
        for(label k=1; k<iterations; k++)
        {            
            std::unordered_map<label,bool> addHaloCells;
            for(auto iter=procBoundaryCells.begin(); iter!=procBoundaryCells.end(); iter++)
            {
                const label cellInd = iter->first;
                if(!(iter->second))
                {
                    const List<Pair<label>>& cellGraph = procMeshGraph[cellInd];
                    for(const Pair<label>& cellNei : cellGraph)
                    {
                        if(cellNei.first()==proc)
                        {
                            auto iter = procBoundaryCells.find(cellNei.second());
                            if(iter==procBoundaryCells.end())
                                addHaloCells.insert({cellNei.second(),false});
                        }
                    }
                    iter->second = true;
                }
            }
            procBoundaryCells.insert(addHaloCells.begin(),addHaloCells.end());            
        }
        
        const List<vector>& procMeshCentres = globalMeshCentres[proc];
        const List<scalar>& procMeshMagn = globalMeshMagn[proc];
        
        DynamicList<CellDescription>& haloCellList_Sorted = globalHaloCellList_Sorted[proc];
        for(auto iter=procBoundaryCells.begin(); iter!=procBoundaryCells.end(); iter++)
        {
            label cellInd = iter->first;
            CellDescription oneCellData = {cellInd,procMeshCentres[cellInd],procMeshMagn[cellInd]};
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
            globalHaloCellToIndexMap[proc][haloCellList_Sorted[index].index] = index;
        }
    }
        
    //Ready parallel find Cell
    mesh.findCell(0.5*(meshBoundingBox.getSmaller()+meshBoundingBox.getLarger()));
    

}

void Foam::Structure::selfCheckHalos()
{
    std::unordered_set<Pair<label>,foamPairHash<label>> neighbors;
    
    for(label cellInd=0; cellInd<localMeshGraph->size(); cellInd++)
    {
        const List<Pair<label>>& cellNei = (*localMeshGraph)[cellInd];
        for(const Pair<label>& nei : cellNei)
        {
            if(nei.first()!=Pstream::myProcNo() && nei.first()!=-1)
            {
                const List<List<Pair<label>>>& procMeshGraph = getMeshGraph(nei.first());
                if(nei.second()<0 || nei.second()>=procMeshGraph.size())
                    FatalErrorInFunction<<"Out of range neighbor"<<exit(FatalError);
                const List<Pair<label>> neiProcCell1 = procMeshGraph[nei.second()];
                bool neighbourOrigCell = false;
                for(const Pair<label>& cell : neiProcCell1)
                {
                    if(cell==Pair<label>(Pstream::myProcNo(),cellInd))
                        neighbourOrigCell = true;
                }
                if(!neighbourOrigCell)
                    FatalErrorInFunction<<"Non matching neighbour"<<exit(FatalError);
                                                
                const std::unordered_map<label,label>& procHaloCellToIndexMap = getHaloCellToIndexMap(nei.first());
                auto iter = procHaloCellToIndexMap.find(nei.second());
                if(iter==procHaloCellToIndexMap.end())
                {
                    Pout<<"("<<cellInd<<") -> "<<nei<<Foam::nl;
                    FatalErrorInFunction<<"Out of range halo"<<exit(FatalError);
                }
                
                const DynamicList<CellDescription>& procHaloCellList = getGlobalHaloCellList_Sorted(nei.first());
                label index = iter->second;
                if(index<0 || index>=procHaloCellList.size())
                    FatalErrorInFunction<<"Out of range index"<<exit(FatalError);
                const CellDescription& cellDesc = procHaloCellList[index];
                if(cellDesc.index!=nei.second())
                    FatalErrorInFunction<<"Mismatch cell index"<<exit(FatalError);
                                
                for(const Pair<label>& haloCell1 : neiProcCell1)
                {
                    if(haloCell1.first()!=Pstream::myProcNo() && haloCell1.first()!=-1)
                    {
                        const List<List<Pair<label>>>& procMeshGraph = getMeshGraph(haloCell1.first());
                        if(haloCell1.second()<0 || haloCell1.second()>=procMeshGraph.size())
                            FatalErrorInFunction<<"Out of range neighbor"<<exit(FatalError);
                        const List<Pair<label>> neiProcCell2 = procMeshGraph[haloCell1.second()];
                                                
                        const std::unordered_map<label,label>& procHaloCellToIndexMap = getHaloCellToIndexMap(haloCell1.first());
                        auto iter = procHaloCellToIndexMap.find(haloCell1.second());
                        if(iter==procHaloCellToIndexMap.end())
                        {
                            Pout<<"cellInd:("<<Pstream::myProcNo()<<" "<<cellInd<<")"<<Foam::nl;
                            Pout<<"cellNei:"<<cellNei<<Foam::nl;
                            Pout<<"nei:"<<nei<<Foam::nl;
                            Pout<<"neiProcCell1:"<<neiProcCell1<<Foam::nl;
                            Pout<<"haloCell1:"<<haloCell1<<Foam::nl;
                            FatalErrorInFunction<<"Out of range halo"<<exit(FatalError);
                        }
                        
                        const DynamicList<CellDescription>& procHaloCellList = getGlobalHaloCellList_Sorted(haloCell1.first());
                        label index = iter->second;
                        if(index<0 || index>=procHaloCellList.size())
                            FatalErrorInFunction<<"Out of range index"<<exit(FatalError);
                        const CellDescription& cellDesc = procHaloCellList[index];
                        if(cellDesc.index!=haloCell1.second())
                            FatalErrorInFunction<<"Mismatch cell index"<<exit(FatalError);
                                                
                        for(const Pair<label>& haloCell2 : neiProcCell2)
                        {
                            if(haloCell2.first()!=Pstream::myProcNo() && haloCell2.first()!=-1)
                            {
                                const List<List<Pair<label>>>& procMeshGraph = getMeshGraph(haloCell2.first());
                                if(haloCell2.second()<0 || haloCell2.second()>=procMeshGraph.size())
                                    FatalErrorInFunction<<"Out of range neighbor"<<exit(FatalError);
                                const List<Pair<label>> neiProcCell3 = procMeshGraph[haloCell2.second()];
                                
                                if(Pstream::myProcNo()==0 && cellInd==280)
                                {
                                    Pout<<haloCell2<<"------------neiProcCell3:"<<neiProcCell3<<Foam::nl;
                                    neighbors.insert(neiProcCell3.begin(),neiProcCell3.end());
                                }
                                
                                const std::unordered_map<label,label>& procHaloCellToIndexMap = getHaloCellToIndexMap(haloCell2.first());
                                auto iter = procHaloCellToIndexMap.find(haloCell2.second());
                                if(iter==procHaloCellToIndexMap.end())
                                {
                                    Pout<<"cellInd:("<<Pstream::myProcNo()<<" "<<cellInd<<")"<<Foam::nl;
                                    Pout<<"cellNei:"<<cellNei<<Foam::nl;
                                    Pout<<"nei:"<<nei<<Foam::nl;
                                    Pout<<"neiProcCell1:"<<neiProcCell1<<Foam::nl;
                                    Pout<<"haloCell1:"<<haloCell1<<Foam::nl;
                                    Pout<<"neiProcCell2:"<<neiProcCell2<<Foam::nl;
                                    Pout<<"haloCell2:"<<haloCell2<<Foam::nl;                                    
                                    FatalErrorInFunction<<"Out of range halo"<<exit(FatalError);
                                }
                                
                                const DynamicList<CellDescription>& procHaloCellList = getGlobalHaloCellList_Sorted(haloCell2.first());
                                label index = iter->second;
                                if(index<0 || index>=procHaloCellList.size())
                                    FatalErrorInFunction<<"Out of range index"<<exit(FatalError);
                                const CellDescription& cellDesc = procHaloCellList[index];
                                if(cellDesc.index!=haloCell2.second())
                                    FatalErrorInFunction<<"Mismatch cell index"<<exit(FatalError);
                            }
                        }
                    }
                }
            }
        }
    }
}

const Foam::DynamicList<Foam::Structure::CellDescription>& Foam::Structure::getGlobalHaloCellList_Sorted
(
    label proc
) const
{
    if(proc<0 || proc>=globalHaloCellList_Sorted.size())
        FatalErrorInFunction<<"Out of bounds proc number"<<exit(FatalError);
    return globalHaloCellList_Sorted[proc];
}

const std::unordered_map<Foam::label,Foam::label>& Foam::Structure::getHaloCellToIndexMap
(
    label proc
) const
{
    if(proc<0 || proc>=globalHaloCellToIndexMap.size())
        FatalErrorInFunction<<"Out of bounds proc number"<<exit(FatalError);
    return globalHaloCellToIndexMap[proc];
}

void Foam::Structure::computeMeshSetup()
{
    generateMeshHaloData(haloWidth);
}

Foam::scalar Foam::Structure::spacingFromMesh
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

void Foam::Structure::neighbourCells
(
    const fvMesh& mesh,
    label cellInd,
    DynamicList<label>& neighbours
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const labelListList& pointCells = mesh.pointCells();
    
    if(cellInd<0 || cellInd>=cells.size())
        FatalErrorInFunction<<"Invalid cellInd"<<exit(FatalError);
    
    std::unordered_set<label> neighbourCellSet;
    for(const label vertice : cells[cellInd].labels(faces))
    {
        for(label neighCell : pointCells[vertice])
        {
            neighbourCellSet.insert(neighCell);
        }
    }
    neighbourCellSet.erase(cellInd);
    
    for(auto iter=neighbours.begin(); iter!=neighbours.end(); iter++)
        neighbours.append(*iter);
}

Foam::label Foam::Structure::findCell
(
    vector position,
    label cellHint
) const
{
    label cellInd = mesh.findCell(position);    
    return cellInd;
}

void Foam::Structure::writeCellSizeField
(
    volScalarField& field
) const
{
    for(label cellInd=0; cellInd<mesh.cells().size(); cellInd++)
        field[cellInd] = spacingFromMesh(mesh,cellInd);
}

void Foam::Structure::quaternionCheck()
{
    Quaternion q;
    scalar epsilon = 1e-4;
    scalar step = 0.1;
    for(scalar qw=0; qw<=1; qw+=step)
    {
        for(scalar qx=0; qx<=1; qx+=step)
        {
            for(scalar qy=0; qy<=1; qy+=step)
            {
                for(scalar qz=0; qz<=1; qz+=step)
                {
                    q = {qw,qx,qy,qz};
                    q.normalize();
                    Quaternion invQuaternion = q.invert();
                    Quaternion oneQuaternion1 = invQuaternion*q;
                    if(std::abs(oneQuaternion1.qw()-1)>epsilon && std::abs(oneQuaternion1.len()-1)>epsilon)
                    {
                        Info<<"oneQuaternion1:"<<oneQuaternion1<<Foam::endl;
                        FatalErrorInFunction<<"Error"<<exit(FatalError);
                    }
                    Quaternion oneQuaternion2 = q*invQuaternion;
                    if(std::abs(oneQuaternion2.qw()-1)>epsilon && std::abs(oneQuaternion2.len()-1)>epsilon)
                    {
                        Info<<"oneQuaternion2:"<<oneQuaternion2<<Foam::endl;
                        FatalErrorInFunction<<"Error"<<exit(FatalError);
                    }
                                        
                    FixedList<FixedList<vector,4>,3> dRdq = Rotation::compute_dRdq(q);
                    for(label pertubeDim=0; pertubeDim<4; pertubeDim++)
                    {
                        Rotation pert_dRdq = {dRdq[0][pertubeDim],dRdq[1][pertubeDim],dRdq[2][pertubeDim]};
                        
                        Quaternion pertLowerQuaternion = q;
                        pertLowerQuaternion[pertubeDim]-=epsilon;
                        Rotation lowerR(pertLowerQuaternion);
                        
                        Quaternion pertUpperQuaternion = q;
                        pertUpperQuaternion[pertubeDim]+=epsilon;
                        Rotation upperR(pertUpperQuaternion);
                        
                        Rotation fd_dRdq = (upperR-lowerR)/(2*epsilon);
                        
                        scalar dist = pert_dRdq.distanceNorm2(fd_dRdq);
                        
                        if(dist>epsilon)
                        {
                            Info<<"---------------------------"<<Foam::endl;
                            Info<<"quaternion:"<<q<<Foam::endl;
                            Rotation R(q);
                            Info<<"R:"<<R<<Foam::endl;
                            Info<<"dRdq:"<<dRdq<<Foam::endl;
                            Info<<"\tpertubeDim:"<<pertubeDim<<Foam::endl;
                            Info<<"\tpertLowerQuaternion:"<<pertLowerQuaternion<<Foam::endl;
                            Info<<"\tlowerR:"<<lowerR<<Foam::endl;
                            Info<<"\tpertUpperQuaternion:"<<pertUpperQuaternion<<Foam::endl;
                            Info<<"\tupperR:"<<upperR<<Foam::endl;
                            Info<<"\tfd_dRdq:"<<fd_dRdq<<Foam::endl;
                            Info<<"\tpert_dRdq:"<<pert_dRdq<<Foam::endl;
                            FatalErrorInFunction<<"Error"<<exit(FatalError);
                        }
                    }
                }
            }
        }
    }
}

void Foam::Structure::rotationCheck()
{
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        //scalar domainStart = this->domainStart(rodNumber);
        //scalar domainEnd = this->domainEnd(rodNumber);
        for(label curveCoeffs=0; curveCoeffs<numberCurveCoeffs(rodNumber); curveCoeffs++)
        {
            for(label coefDim=0; coefDim<3; coefDim++)
            {
            }
        }
    }
}

void Foam::Structure::transformationCheck()
{
    OFstream printOutFile(mesh.rootPath()+"/"+mesh.caseName()+"/transformationCheck");
    Info<<printOutFile.name()<<Foam::endl;
    scalar nbrSteps = 20;
    scalar epsilon = 1e-2;
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        printOutFile<<"rodNumber:"<<rodNumber<<Foam::endl;
        scalar domainStart = this->domainStart(rodNumber)+epsilon;
        scalar domainEnd = this->domainEnd(rodNumber)-epsilon;
        scalar delta = domainEnd-domainStart;
        scalar stepsize = delta/nbrSteps;
        for(label curveCoeffs=0; curveCoeffs<numberCurveCoeffs(rodNumber); curveCoeffs++)
        {
            printOutFile<<"\b curveCoeffs:"<<curveCoeffs<<Foam::endl;
            for(label coefDim=0; coefDim<3; coefDim++)
            {
                printOutFile<<"\b\b coefDim:"<<coefDim<<Foam::endl;
                for(scalar parameter=domainStart; parameter<=domainEnd; parameter+=stepsize)
                {
                    // Compute gradients
                    vector drdC = get_drdC(rodNumber,curveCoeffs,coefDim,parameter);
                    Quaternion dqdC = get_coeffDerivQuaternions(rodNumber,curveCoeffs,coefDim,parameter);
                    Rotation dRdC = Rotation::compute_dRdX(dqdC,get_Quaternions(rodNumber,parameter));
                    
                    // Get basic value of coefficient
                    scalar coeffBasicValue = getCurveCoeff(rodNumber,curveCoeffs,coefDim);
                    
                    // Lower value
                    vector lower_r,lower_d1,lower_d2,lower_d3;
                    scalar lowerCoeffValue = coeffBasicValue-epsilon;
                    setCurveCoeff(rodNumber,curveCoeffs,coefDim,lowerCoeffValue);
                    //vector lower_p = rodEval(Rods[rodNumber],parameter);
                    Quaternion lower_Q = m_Rot_Eval(rodNumber,parameter);
                    Rotation lower_R(lower_Q);
                    vector lower_para = rodEval(Rods[rodNumber],parameter-epsilon);
                    vector upper_para = rodEval(Rods[rodNumber],parameter+epsilon);
                    //vector lower_dPara = (upper_para-lower_para)/(2*epsilon);
                    rodEval(Rods[rodNumber],parameter,lower_d1,lower_d2,lower_d3,lower_r);
                    
                    // Upper value
                    vector upper_r,upper_d1,upper_d2,upper_d3;
                    scalar upperCoeffValue = coeffBasicValue+epsilon;
                    setCurveCoeff(rodNumber,curveCoeffs,coefDim,upperCoeffValue);
                    //vector upper_p = rodEval(Rods[rodNumber],parameter);
                    Quaternion upper_Q = m_Rot_Eval(rodNumber,parameter);
                    Rotation upper_R(upper_Q);
                    lower_para = rodEval(Rods[rodNumber],parameter-epsilon);
                    upper_para = rodEval(Rods[rodNumber],parameter+epsilon);
                    //vector upper_dPara = (upper_para-lower_para)/(2*epsilon);
                    rodEval(Rods[rodNumber],parameter,upper_d1,upper_d2,upper_d3,upper_r);

                    vector fd_drdC = (upper_r-lower_r)/(2*epsilon);
                    Rotation fd_dRdC = (upper_R-lower_R)/(2*epsilon);
                        
                    printOutFile<<"\b\b\b parameter:"<<parameter<<" ->  dRdC:"<<dRdC.get_d3()<<" / fd_dRdC:"<<fd_dRdC.get_d3()<<"   dist:"<<vectorDistance(dRdC.get_d3(),fd_dRdC.get_d3())<<Foam::endl;
                    
                    /*
                    if(fd_dRdC.distanceNorm2(dRdC)>epsilon)
                    {
                        Info<<"parameter:"<<parameter<<Foam::endl;
                        Info<<"dRdC:"<<dRdC<<Foam::endl;
                        Info<<"fd_dRdC:"<<fd_dRdC<<Foam::endl;

                        Info<<"--"<<Foam::endl;
                        Info<<"lower_Q:"<<lower_Q<<Foam::endl;
                        Info<<"lower_R:"<<lower_R<<Foam::endl;
                        Info<<"lower_dPara:"<<lower_dPara<<Foam::endl;                        Info<<"lower_r:"<<lower_r<<Foam::endl;
                        Info<<"lower_d1:"<<lower_d1<<Foam::endl;                        Info<<"lower_d2:"<<lower_d2<<Foam::endl;
                        Info<<"lower_d3:"<<lower_d3<<Foam::endl; 
                        Info<<"--"<<Foam::endl;
                        Info<<"upper_Q:"<<upper_Q<<Foam::endl;
                        Info<<"upper_R:"<<upper_R<<Foam::endl;
                        Info<<"upper_dPara:"<<upper_dPara<<Foam::endl;
                        Info<<"upper_r:"<<upper_r<<Foam::endl;
                        Info<<"upper_d1:"<<upper_d1<<Foam::endl;
                        Info<<"upper_d2:"<<upper_d2<<Foam::endl;
                        Info<<"upper_d3:"<<upper_d3<<Foam::endl;
                        FatalErrorInFunction<<"Error"<<exit(FatalError);
                    }
                    */
                    if(vectorDistance(fd_drdC,drdC)>epsilon)
                    {
                        Info<<"parameter:"<<parameter<<Foam::endl;
                        Info<<"drdC:"<<drdC<<Foam::endl;
                        Info<<"lower_r:"<<lower_r<<Foam::endl;
                        Info<<"upper_r:"<<upper_r<<Foam::endl;
                        Info<<"fd_drdC:"<<fd_drdC<<Foam::endl;
                        FatalErrorInFunction<<"Error"<<exit(FatalError);
                    }
                    setCurveCoeff(rodNumber,curveCoeffs,coefDim,coeffBasicValue);
                }
            }
        }
    }
}

void Foam::Structure::parameterGradientCheck()
{
    Info<<"Structure::parameterGradientCheck (rodDerivEval,rodDerivEval2)"<<Foam::endl;

    scalar nbrSteps = 20;
    scalar epsilon = 1e-4;
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        //printOutFile<<"rodNumber:"<<rodNumber<<Foam::endl;
        scalar domainStart = this->domainStart(rodNumber)+5*epsilon;
        scalar domainEnd = this->domainEnd(rodNumber)-5*epsilon;
        scalar delta = domainEnd-domainStart;
        scalar stepsize = delta/nbrSteps;
        for(scalar parameter=domainStart; parameter<=domainEnd; parameter+=stepsize)
        {
            //Info<<"parameter:"<<parameter<<Foam::endl;
            // Compute gradients
            vector dd1dp, dd2dp, dd3dp, drdp;
            rodDerivEval(Rods[rodNumber],parameter,dd1dp,dd2dp,dd3dp,drdp);
            Rotation dRdp(dd1dp,dd2dp,dd3dp);
            vector d2d1dp, d2d2dp, d2d3dp, d2rdp;
            rodDeriv2Eval(Rods[rodNumber],parameter,d2d1dp,d2d2dp,d2d3dp,d2rdp);
            Rotation d2Rdp(dd1dp,dd2dp,dd3dp);
            
            // Lower value
            scalar lower_parameter = parameter-epsilon;
            vector lower_d1,lower_d2,lower_d3,lower_r;
            rodEval(Rods[rodNumber],lower_parameter,lower_d1,lower_d2,lower_d3,lower_r);
            Rotation lower_R(lower_d1,lower_d2,lower_d3);
            vector lower_dd1dp, lower_dd2dp, lower_dd3dp, lower_drdp;
            rodDerivEval(Rods[rodNumber],lower_parameter,lower_dd1dp,lower_dd2dp,lower_dd3dp,lower_drdp);
            Rotation lower_dRdp(lower_dd1dp,lower_dd2dp,lower_dd3dp);
            
            // Upper value
            scalar upper_parameter = parameter+epsilon;
            vector upper_d1,upper_d2,upper_d3,upper_r;
            rodEval(Rods[rodNumber],upper_parameter,upper_d1,upper_d2,upper_d3,upper_r);
            Rotation upper_R(upper_d1,upper_d2,upper_d3);
            vector upper_dd1dp, upper_dd2dp, upper_dd3dp, upper_drdp;
            rodDerivEval(Rods[rodNumber],upper_parameter,upper_dd1dp,upper_dd2dp,upper_dd3dp,upper_drdp);
            Rotation upper_dRdp(upper_dd1dp,upper_dd2dp,upper_dd3dp);
            
            vector fd_drdp = (upper_r-lower_r)/(2*epsilon);
            Rotation fd_dRdp = (upper_R-lower_R)/(2*epsilon);
            vector fd_d2rdp = (upper_drdp-lower_drdp)/(2*epsilon);
            Rotation fd_d2Rdp = (upper_dRdp-lower_dRdp)/(2*epsilon);
            
            scalar error_dRdp = dRdp.distanceNorm2(fd_dRdp);
            scalar error_d2Rdp = d2Rdp.distanceNorm2(fd_d2Rdp);
                        
            if(vectorDistance(fd_drdp,drdp)>epsilon)
            {
                Info<<"parameter:"<<parameter<<Foam::endl;
                Info<<"drdp:"<<drdp<<Foam::endl;
                Info<<"lower_r:"<<lower_r<<Foam::endl;
                Info<<"upper_r:"<<upper_r<<Foam::endl;
                Info<<"fd_drdp:"<<fd_drdp<<Foam::endl;
                FatalErrorInFunction<<"Error"<<exit(FatalError);
            }
            if(vectorDistance(fd_d2rdp,d2rdp)>epsilon)
            {
                Info<<"parameter:"<<parameter;//<<Foam::endl;
                //Info<<"lower_parameter:"<<lower_parameter<<Foam::endl;
                //Info<<"upper_parameter:"<<upper_parameter<<Foam::endl;
                Info<<" d2rdp:"<<d2rdp;//<<Foam::endl;
                //Info<<"lower_drdp:"<<lower_drdp<<Foam::endl;
                //Info<<"upper_drdp:"<<upper_drdp<<Foam::endl;
                Info<<" fd_d2rdp:"<<fd_d2rdp<<Foam::endl;
                //FatalErrorInFunction<<"Error"<<exit(FatalError);
            }
            if(error_dRdp>epsilon)
            {
                Info<<"parameter:"<<parameter<<Foam::endl;
                Info<<"dRdp:"<<dRdp<<Foam::endl;
                Info<<"lower_R:"<<lower_R<<Foam::endl;
                Info<<"upper_R:"<<upper_R<<Foam::endl;
                Info<<"fd_dRdp:"<<fd_dRdp<<Foam::endl;
                FatalErrorInFunction<<"Error"<<exit(FatalError);
            }
            if(error_d2Rdp>epsilon)
            {
                Info<<"parameter:"<<parameter<<Foam::endl;
                Info<<"d2Rdp:"<<d2Rdp<<Foam::endl;
                Info<<"lower_dRdp:"<<lower_dRdp<<Foam::endl;
                Info<<"upper_dRdp:"<<upper_dRdp<<Foam::endl;
                Info<<"fd_d2Rdp:"<<fd_d2Rdp<<Foam::endl;
                FatalErrorInFunction<<"Error"<<exit(FatalError);
            }            
        }
    }
}

void Foam::Structure::selfCheck()
{
    quaternionCheck();
    transformationCheck();
    
    FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
    
    Info<<"-----------Structure derivatives-----------"<<Foam::endl;
    
    std::function<void(std::function<Quaternion(scalar)>,
                       std::function<Quaternion(scalar)>,
                       scalar,scalar,uint)> tComparer =
    [](auto deriv, auto fdDeriv, scalar minPar, scalar maxPar, uint steps)
    {
        scalar deltaPar = maxPar-minPar;
        scalar stepsizePar = deltaPar/steps;
        for(scalar currPara=minPar; currPara<=maxPar; currPara+=stepsizePar)
        {
            Quaternion derivValue = deriv(currPara);
            Quaternion fderivValue = fdDeriv(currPara);
            Quaternion errorVec = derivValue-fderivValue;
            Quaternion error;
            Quaternion derivValueLen;
            Quaternion fderivValueLen;
            Quaternion avgLen;
            Quaternion percError;
            scalar maxError = std::numeric_limits<scalar>::min();
            scalar maxPercError = std::numeric_limits<scalar>::min();
            for(label i=0; i<4; i++)
            {
                error[i] = std::sqrt(errorVec[i]*errorVec[i]);
                derivValueLen[i] = std::sqrt(derivValue[i]*derivValue[i]);
                fderivValueLen[i] = std::sqrt(fderivValue[i]*fderivValue[i]);
                avgLen[i] = 0.5*(derivValueLen[i]+fderivValueLen[i]);
                if(avgLen[i]!=0)
                    percError[i] = error[i]/avgLen[i];
                else
                    percError[i] = error[i];
                maxPercError = std::max(maxPercError,percError[i]);
                maxError = std::max(maxError,error[i]);
            }

            Info<<"("<<currPara<<"): Err:"<<error<<" percErr:"<<percError<<" // grad:"<<derivValue<<" -- fdGrad:"<<fderivValue<<Foam::endl;
            
            FatalErrorInFunction<<"Temp stop!"<<exit(FatalError);
            if(maxPercError>2e-3 && maxError>1e-4)
            {
                Info<<"("<<currPara<<"): Err:"<<error<<" percErr:"<<percError<<" // grad:"<<derivValue<<" -- fdGrad:"<<fderivValue<<Foam::endl;
                FatalErrorInFunction<<"t Comparison failed!"<<exit(FatalError);
            }
        }
    };
    
    std::function<void(std::function<Quaternion(scalar)>,
                       std::function<Quaternion(scalar)>,
                       scalar,scalar,uint)> qComparer =
    [](auto deriv, auto fdDeriv, scalar minPar, scalar maxPar, uint steps)
    {
        scalar deltaPar = maxPar-minPar;
        scalar stepsizePar = deltaPar/steps;
        
        for(scalar currPara=minPar; currPara<=maxPar; currPara+=stepsizePar)
        {
            Quaternion derivValue = deriv(currPara);
            Quaternion fderivValue = fdDeriv(currPara);
            Quaternion error;
            scalar derivValueLen = std::sqrt(derivValue[0]*derivValue[0]+derivValue[1]*derivValue[1]
                                            +derivValue[2]*derivValue[2]+derivValue[3]*derivValue[3]);
            scalar fderivValueLen = std::sqrt(fderivValue[0]*fderivValue[0]+fderivValue[1]*fderivValue[1]
                                             +fderivValue[2]*fderivValue[2]+fderivValue[3]*fderivValue[3]);
            scalar avgLen = 0.5*(derivValueLen+fderivValueLen);
            scalar maxError = std::numeric_limits<scalar>::min();;
            for(label i=0; i<4; i++)
            {
                error[i]= std::abs(derivValue[i]-fderivValue[i]);
                maxError = std::max(maxError,error[i]);
            }
            
            scalar maxPercError;
            if(avgLen!=0)
                maxPercError = maxError / avgLen;
            else
                maxPercError = maxError;
            
            Info<<"("<<currPara<<"): Err:"<<error<<" maxPercError:"<<maxPercError<<" maxError:"<<maxError<<" avgLen:"<<avgLen<<" // grad:"<<derivValue<<" -- fdGrad:"<<fderivValue<<Foam::endl;
            
            FatalErrorInFunction<<"Temp stop!"<<exit(FatalError);
            
            if(maxPercError>2e-3 && maxError>1e-4)
            {
                Info<<"("<<currPara<<"): Err:"<<error<<" maxPercError:"<<maxPercError<<" // grad:"<<derivValue<<" -- fdGrad:"<<fderivValue<<Foam::endl;
                FatalErrorInFunction<<"q Comparison failed!"<<exit(FatalError);
            }
        }
    };
    
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        scalar domainStart = this->domainStart(rodNumber)+0.5;
        scalar domainEnd = this->domainEnd(rodNumber);
        for(label curveCoeffs=0; curveCoeffs<numberCurveCoeffs(rodNumber); curveCoeffs++)
        {
            for(label coefDim=0; coefDim<3; coefDim++)
            {
                std::function<Quaternion(scalar)> derivQ = [&](scalar par)
                {
                    Quaternion deriv = m_Rot_Eval_Deriv(rodNumber,curveCoeffs,coefDim,par);
                    return deriv;
                };
                std::function<Quaternion(scalar)> fdDerivQ = [&](scalar par)
                {
                    scalar epsilon=1e-3;
                    scalar coeffBasicValue = getCurveCoeff(rodNumber,curveCoeffs,coefDim);
                    
                    scalar lowerCoeffValue = coeffBasicValue-epsilon;
                    setCurveCoeff(rodNumber,curveCoeffs,coefDim,lowerCoeffValue);
                    Quaternion lower_q = m_Rot_Eval(rodNumber,par);

                    scalar upperCoeffValue = coeffBasicValue+epsilon;                    
                    setCurveCoeff(rodNumber,curveCoeffs,coefDim,upperCoeffValue);
                    Quaternion upper_q = m_Rot_Eval(rodNumber,par);
                    
                    Quaternion diff = upper_q-lower_q;
                    scalar diffCoeff = upperCoeffValue-lowerCoeffValue;

                    setCurveCoeff(rodNumber,curveCoeffs,coefDim,coeffBasicValue);
                    if(diffCoeff!=0)
                    {
                        diff[0] /= diffCoeff;
                        diff[1] /= diffCoeff;
                        diff[2] /= diffCoeff;
                        diff[3] /= diffCoeff;
                        return diff;
                    }
                    else
                        FatalErrorInFunction<<"No delta in coeff values!"<<exit(FatalError);
                    return diff;
                };
                qComparer(derivQ,fdDerivQ,domainStart,domainEnd,20);
                
                /*
                std::function<FixedList<vector,4>(scalar)> derivT = [&](scalar par)
                {
                    vector d1dC,d2dC,d3dC,rdC;
                    rodEvalDerivCoeff(rodNumber,curveCoeffs,coefDim,par,d1dC,d2dC,d3dC,rdC);
                    FixedList<vector,4> deriv = {d1dC,d2dC,d3dC,rdC};
                    return deriv;
                };
                std::function<FixedList<vector,4>(scalar)> fdDerivT = [&](scalar par)
                {
                    scalar epsilon=0.2;
                    scalar coeffBasicValue = getCurveCoeff(rodNumber,curveCoeffs,coefDim);
                    vector d1,d2,d3,r;
                    rodEval(Rods[rodNumber],par,d1,d2,d3,r);
                    Info<<"d1:"<<d1<<Foam::endl;
                    Info<<"d2:"<<d2<<Foam::endl;
                    Info<<"d3:"<<d3<<Foam::endl;
                    Info<<"r:"<<r<<Foam::endl<<Foam::endl;

                    scalar lowerCoeffValue = coeffBasicValue-epsilon;
                    setCurveCoeff(rodNumber,curveCoeffs,coefDim,lowerCoeffValue);
                    vector lower_d1,lower_d2,lower_d3,lower_r;
                    rodEval(Rods[rodNumber],par,lower_d1,lower_d2,lower_d3,lower_r);

                    scalar upperCoeffValue = coeffBasicValue+epsilon;                    
                    setCurveCoeff(rodNumber,curveCoeffs,coefDim,upperCoeffValue);
                    vector upper_d1,upper_d2,upper_d3,upper_r;
                    rodEval(Rods[rodNumber],par,upper_d1,upper_d2,upper_d3,upper_r);
                    
                    Info<<"lower_d1:"<<lower_d1<<Foam::endl;
                    Info<<"lower_d2:"<<lower_d2<<Foam::endl;
                    Info<<"lower_d3:"<<lower_d3<<Foam::endl;
                    Info<<"lower_r:"<<lower_r<<Foam::endl<<Foam::endl;
                    
                    Info<<"upper_d1:"<<upper_d1<<Foam::endl;
                    Info<<"upper_d2:"<<upper_d2<<Foam::endl;
                    Info<<"upper_d3:"<<upper_d3<<Foam::endl;
                    Info<<"upper_r:"<<upper_r<<Foam::endl<<Foam::endl;


                    vector diff_d1 = upper_d1-lower_d1;
                    vector diff_d2 = upper_d2-lower_d2;
                    vector diff_d3 = upper_d3-lower_d3;
                    vector diff_r = upper_r-lower_r;
                    scalar diffCoeff = upperCoeffValue-lowerCoeffValue;
                    diff_d1 /= diffCoeff;
                    diff_d2 /= diffCoeff;
                    diff_d3 /= diffCoeff;
                    diff_r /= diffCoeff;

                    setCurveCoeff(rodNumber,curveCoeffs,coefDim,coeffBasicValue);
                    if(diffCoeff!=0)
                    {
                        FixedList<vector,4> fdderiv = {diff_d1,diff_d2,diff_d3,diff_r};
                        return fdderiv;
                    }
                    else
                        FatalErrorInFunction<<"No delta in coeff values!"<<exit(FatalError);
                };
                tComparer(derivT,fdDerivT,domainStart,domainEnd,20);
                */
            }
        }
    }
    
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
