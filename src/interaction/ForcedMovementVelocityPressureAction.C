#include "ForcedMovementVelocityPressureAction.H"
#include "codeStream.H"

Foam::ForcedMovementVelocityPressureAction::ForcedMovementVelocityPressureAction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volVectorField& input_U,
    volVectorField& output_Uf,
    const IOdictionary& structureDict,
    std::shared_ptr<MeshRefiner> refinement_,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
VelocityPressureForceInteraction(mesh,structure,input_U,output_Uf,structureDict,refinement_,modusFieldToMarker,modusMarkerToField)
{
    Info<<"Created ForcedMovementVelocityPressureAction"<<Foam::endl;
    constructMovementFunction();
}

Foam::ForcedMovementVelocityPressureAction::~ForcedMovementVelocityPressureAction()
{
    dlclose(rodMoveSo);
}

void Foam::ForcedMovementVelocityPressureAction::preSolveMovement()
{
    structure.pushBackDeformationState();
    std::unique_ptr<List<List<vector>>> allRodsDeformation = readDeformationDict();
    structure.setDeformation(*allRodsDeformation);
    moveMarkers();
}

void Foam::ForcedMovementVelocityPressureAction::preSolveMarkerMeshAdaption()
{
    meshMarkerAdaptation();
};

std::unique_ptr<Foam::List<Foam::List<Foam::vector>>> Foam::ForcedMovementVelocityPressureAction::readDeformationDict()
{   
    const dictionary& rodMovementFieldDict = structureDict.subDict("rodMovementField");
    List<keyType> rodMovementFieldKeys = rodMovementFieldDict.keys();
    
    auto movementListPtr = std::make_unique<List<List<vector>>>(rodMovementFieldKeys.size());
    List<List<vector>>& movementList = *movementListPtr;
    
    if(rodMovementFieldKeys.size()!=structure.getNumberRods())
    {
        Info<<"rodMovementFieldKeys:"<<rodMovementFieldKeys<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in movement field to rod number!"<<exit(FatalError);
    }
    
    for(label rodNumber=0; rodNumber<structure.getNumberRods(); rodNumber++)
    {
        const gsNurbs<scalar>& deformation = structure.getDeformation(rodNumber);
        unsigned int nbrDefParameters = deformation.coefs().rows();
        std::vector<scalar> Px(nbrDefParameters),Py(nbrDefParameters),Pz(nbrDefParameters);
        movementFunction[rodNumber]
        (
            deformation.knots().degree(),deformation.knots().data(),deformation.knots().size(),
            Px.data(),Py.data(),Pz.data(),nbrDefParameters,mesh.time().value()
        );
        List<vector>& movementListRod = movementList[rodNumber];
        movementListRod.resize(nbrDefParameters);
        for(std::size_t index=0; index<nbrDefParameters; index++)
        {
            movementListRod[index] = vector(Px[index],Py[index],Pz[index]);
        }
        Info<<"movementListRod:"<<movementListRod<<Foam::endl;
    }
    return movementListPtr;
}

Foam::vector Foam::ForcedMovementVelocityPressureAction::getVelocity
(
    const LagrangianMarker* marker
)
{
    return marker->getMarkerVelocity();
}

void Foam::ForcedMovementVelocityPressureAction::constructMovementFunction()
{
    Info<<"Foam::ForcedMovementVelocityPressureAction::constructMovementFunction()"<<Foam::endl;
    const dictionary& rodMovementFieldDict = structureDict.subDict("rodMovementField");
    List<keyType> rodMovementFieldKeys = rodMovementFieldDict.keys();
    Info<<"rodMovementFieldKeys:"<<rodMovementFieldKeys<<Foam::endl;
    
    if(rodMovementFieldKeys.size()!=structure.getNumberRods())
    {
        Info<<"rodMovementFieldKeys:"<<rodMovementFieldKeys<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in movement field to rod number!"<<exit(FatalError);
    }
    
    movementFunction.resize(structure.getNumberRods());
    deformationNurbsCodes.resize(structure.getNumberRods());
    
    for(label rodNumber=0; rodNumber<structure.getNumberRods(); rodNumber++)
    {
        keyType oneRodMoveFieldKey = rodMovementFieldKeys[rodNumber];
        Info<<"rodNumber:"<<rodNumber<<" : "<<oneRodMoveFieldKey<<Foam::endl;
        const dictionary& oneRodMovementDict = rodMovementFieldDict.subDict(oneRodMoveFieldKey);
        ITstream moveFunctionStream = oneRodMovementDict.lookup("move");
        token moveFunctionToken;
        moveFunctionStream.read(moveFunctionToken);
        if(!moveFunctionToken.isString())
        {
            FatalErrorInFunction<<"Invalid entry in constant/structureDict/"<<oneRodMoveFieldKey<<"/move -- must be  string"<<exit(FatalError);
        }
        deformationNurbsCodes[rodNumber] = moveFunctionToken.stringToken();
    }
    
    std::ofstream dynCode("constant/rodMovementCode.cpp");
    dynCode<<"#include <math.h>\n"<<"extern \"C\" {\n";
    for(std::size_t rodNumber=0; rodNumber<deformationNurbsCodes.size(); rodNumber++)
    {
        std::string moveFunction = deformationNurbsCodes[rodNumber];
        std::string nameFunction = moveFunction.substr(0,21);
        if(nameFunction!="void deformationNurbs")
            FatalErrorInFunction<<"Function name is wrong:"<<rodNumber<<" -->"<<nameFunction<<"<--"<<exit(FatalError);
        /*
        std::string parameters = moveFunction.substr(23,220);
        std::string::size_type substringStart = moveFunction.find(rodMovementSignature);
        if(substringStart!=0)
        {
            Info<<"nameFunction:"<<nameFunction<<Foam::endl;
            Info<<"parameters:"<<parameters<<Foam::endl;
            FatalErrorInFunction<<"Invalid signature in rod:"<<rodNumber<<Foam::endl<<
            "is: -->"<<moveFunction.substr(0,230)<<"<--"<<Foam::endl<<
            "should be: -->"<<rodMovementSignature<<"<--"<<Foam::endl<<
            exit(FatalError);
        }
        */
        moveFunction.insert(21,std::to_string(rodNumber));        
        dynCode<<moveFunction<<"\n";
    }
    dynCode<<"}";
    dynCode.close();
    int result = system("g++ constant/rodMovementCode.cpp -o constant/rodMovementCode.so -shared -fPIC");
    if(result!=0)
        FatalErrorInFunction<<"Compile of dynamic code failed"<<exit(FatalError);
        
    rodMoveSo = dlopen("constant/rodMovementCode.so", RTLD_NOW);
    if(rodMoveSo==nullptr)
        FatalErrorInFunction<<"Loading of rodMovementCode.so failed"<<exit(FatalError);
    for(std::size_t rodNumber=0; rodNumber<deformationNurbsCodes.size(); rodNumber++)
    {
        void* basisFunctionPtr = dlsym(rodMoveSo, ("deformationNurbs"+std::to_string(rodNumber)).c_str());
        if(basisFunctionPtr==nullptr)
        {
            char *errstr;
            errstr = dlerror();
            if (errstr != NULL)
                printf ("A dynamic linking error occurred: (%s)\n", errstr);
            FatalErrorInFunction<<"Loading of symbol (deformationNurbs) from rodMovementCode.so failed"<<exit(FatalError);
        }
        movementFunction[rodNumber] = nullptr;
        movementFunction[rodNumber] = reinterpret_cast<void(*)(unsigned int,const double*,unsigned int,double*,double*,double*,unsigned int,double)>(basisFunctionPtr);
        if(movementFunction[rodNumber]==nullptr)
        {
            FatalErrorInFunction<<"Casting of function pointer failed"<<exit(FatalError);
        }
    }    
}
