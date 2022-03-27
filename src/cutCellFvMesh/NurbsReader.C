#include "NurbsReader.H"

Foam::NurbsReader::NurbsReader(fileName runDirectory, fileName caseName)
{
    this-> runDirectory = runDirectory;
    this-> caseName = caseName;
    fullXMLPath = getXMLPath();
    nurbsCurves = readOutNurbsFromXML();
}

std::shared_ptr<std::vector<Nurbs1D>> Foam::NurbsReader::getNurbsCurves()
{
    return nurbsCurves;
}

Foam::word Foam::NurbsReader::getXMLPath()
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

std::shared_ptr<std::vector<Nurbs1D>> Foam::NurbsReader::readOutNurbsFromXML()
{
    std::shared_ptr<std::vector<Nurbs1D>> nurbsCurves = make_shared<std::vector<Nurbs1D>>();
    
    XMLDocument nurbsDoc;
	XMLError xmlErrorID = nurbsDoc.LoadFile(fullXMLPath.c_str());
    if(xmlErrorID!=0)
        FatalIOError<<"Reading Nurbs file failed!"<<exit(FatalIOError);
    
    XMLElement* xmlRoot = nurbsDoc.RootElement();
    if(xmlRoot==NULL || word(xmlRoot->Name()).compare("xml")!=0)
            FatalIOError<<"Reading Nurbs file failure! Root name must be xml."<<exit(FatalIOError);
    
    for(XMLNode* geometryNode=xmlRoot->FirstChild();geometryNode!=NULL;geometryNode=geometryNode->NextSibling())
    {
        XMLElement* geometry = geometryNode->ToElement();
        //Second level
        const char* geometryType = geometry->Attribute("type");
        if(word(geometry->Value()).compare("Geometry")!=0 ||
           geometryType==NULL || 
           (word(geometryType).compare("Nurbs")!=0)
           )
            FatalIOError<<"Reading Nurbs file failure! Item name must be Geometry."<<exit(FatalIOError);
        
        //Third level
        XMLElement* basis=geometry->FirstChild()->ToElement();
        const char* basisType = basis->Attribute("type");
        XMLElement* coefs=basis->NextSibling()->ToElement();
        int coefsDim = coefs->IntAttribute("geoDim",-1);
        if(basis==NULL || word(basis->Value()).compare("Basis")!=0)
            FatalIOError<<"Reading Nurbs file failure! Geometry's first child name must be Basis."<<exit(FatalIOError);
        if(basisType==NULL || word(basisType).compare("NurbsBasis")!=0)
            FatalIOError<<"Reading Nurbs file failure! Geometry's basis must be NurbsBasis."<<exit(FatalIOError);
        if(coefs==NULL || word(coefs->Value()).compare("coefs")!=0)
            FatalIOError<<"Reading Nurbs file failure! Geometry's second child name must be coefs."<<exit(FatalIOError);
        if(coefsDim != 3)
            FatalIOError<<"Reading Nurbs file failure! Coefficient dimension must be 3."<<exit(FatalIOError);
        if(coefs->NextSibling()!=NULL)
            FatalIOError<<"Reading Nurbs file failure! Geometry must not have three children."<<exit(FatalIOError);
        
        //Fourth level
        XMLElement* basisbasis=basis->FirstChild()->ToElement();
        const char* basisbasisType = basisbasis->Attribute("type");
        XMLElement* weights=basisbasis->NextSibling()->ToElement();
        const char* coefField = coefs->GetText();
        if(basisbasis==NULL || word(basisbasis->Value()).compare("Basis")!=0)
            FatalIOError<<"Reading Nurbs file failure! Basis first child name must be Basis."<<exit(FatalIOError);
        if(basisbasisType==NULL || word(basisbasisType).compare("BSplineBasis")!=0)
            FatalIOError<<"Reading Nurbs file failure! Geometry's basis must be NurbsBasis."<<exit(FatalIOError);
        if(weights==NULL || word(weights->Value()).compare("weights")!=0)
            FatalIOError<<"Reading Nurbs file failure! Basis first child name must be weights."<<exit(FatalIOError);
        if(weights->NextSibling()!=NULL)
            FatalIOError<<"Reading Nurbs file failure! Basis must not have three children."<<exit(FatalIOError);
        if(coefField==NULL)
            FatalIOError<<"Reading Nurbs file failure! coefs value must not be null."<<exit(FatalIOError);
        
        //Fifth level
        XMLElement* KnotVector=basisbasis->FirstChild()->ToElement();
        int KnotVectorDegree = KnotVector->IntAttribute("degree");
        const char* weightField=weights->GetText();
        if(KnotVector==NULL || word(KnotVector->Value()).compare("KnotVector")!=0)
            FatalIOError<<"Reading Nurbs file failure! BasisBasis name must be KnotVector."<<exit(FatalIOError);
        if(KnotVector->NextSibling()!=NULL)
            FatalIOError<<"Reading Nurbs file failure! BasisBasis must not have two children."<<exit(FatalIOError);        
        if(weightField==NULL)
            FatalIOError<<"Reading Nurbs file failure! weights first child must not be null."<<exit(FatalIOError);
        
        //Sixth level
        const char* knotField=KnotVector->GetText();
        if(knotField==NULL)
            FatalIOError<<"Reading Nurbs file failure! Knots must not be null."<<exit(FatalIOError);    
        
        
        std::stringstream knotString(knotField);
        std::stringstream weightString(weightField);
        std::stringstream coefString(coefField);
        
        
        //Info<<"knotString:"<<knotField<<endl;
        //Info<<"weightString:"<<weightField->Value()<<endl;
        //Info<<"coefString:"<<coefField->Value()<<endl;
        label nurbsDegree = KnotVectorDegree;
        
        DynamicList<scalar> knotList;
        scalar number;
        while(knotString>>number)
        {
            knotList.append(number);
        }
        Info<<"knotList:"<<knotList<<endl;
        
        DynamicList<scalar> weightList;
        while(weightString>>number)
        {
            weightList.append(number);
        }
        Info<<"weightList:"<<weightList<<endl;
        
        DynamicList<vector> coefList;
        DynamicList<scalar> numberList;
        while(coefString>>number)
        {
            numberList.append(number);
        }
        if(numberList.size()%3!=0)
            FatalIOError<<"Reading Nurbs file failure! Coefs must have three dimensions."<<exit(FatalIOError);
        for(int i=0;i<numberList.size();i+=3)
        {
            coefList.append(vector(numberList[i],numberList[i+1],numberList[i+2]));
        }            
        Info<<"coefList:"<<coefList<<endl;
        Info<<"nurbsDegree:"<<nurbsDegree<<endl;        
        
        nurbsCurves->push_back(Nurbs1D(knotList, coefList, weightList, nurbsDegree));
    }
    return nurbsCurves;
}
