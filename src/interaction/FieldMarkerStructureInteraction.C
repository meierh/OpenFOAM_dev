#include "FieldMarkerStructureInteraction.H"

Foam::FieldMarkerStructureInteraction::FieldMarkerStructureInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    const IOdictionary& structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
mesh(mesh),
structure(structure),
structureDict(structureDict),
h(std::cbrt(mesh.cells()[0].mag(mesh.points(),mesh.faces()))),
modusFieldToMarker(modusFieldToMarker),
modusMarkerToField(modusMarkerToField)
{
    Info<<"Completed FieldMarkerStructureInteraction setup"<<Foam::endl;
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
        for(uint i=0;i<knots.size();i++)
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

Foam::FieldMarkerStructureInteraction::MarkerInfoFiles::MarkerInfoFiles
(
    const IOdictionary& structureDict,
    Foam::word dictName,
    bool masterFile,
    List<word> header
):
masterFile(masterFile),
columns(header.size()),
fileActive(false)
{
    if(structureDict.found(dictName))
    {
        ITstream markerInfoFilesStream = structureDict.lookup(dictName);
        token markerInfoFilesToken;
        markerInfoFilesStream.read(markerInfoFilesToken);
        if(!markerInfoFilesToken.isString())
            FatalErrorInFunction<<"Invalid entry in structure/structureDict/"<<dictName<<" -- must be string"<<exit(FatalError);
        fileName = markerInfoFilesToken.stringToken();
        fileActive = true;
        if(masterFile)
        {
            if(Pstream::master())
                filePtr = std::make_unique<std::ofstream>(fileName);
        }
        else
        {
            filePtr = std::make_unique<std::ofstream>(fileName);
        }
    }
    
    if(filePtr)
    {
        for(label headerInd=0; headerInd<header.size(); headerInd++)
        {
            (*filePtr)<<header[headerInd];
            if(headerInd<header.size()-1)
                (*filePtr)<<",";
            else
                (*filePtr)<<std::endl;
        }
    }
}

void Foam::FieldMarkerStructureInteraction::MarkerInfoFiles::write
(
    List<scalar> values
)
{
    if(values.size()!=columns)
        FatalErrorInFunction<<"Mismatch in column number"<<exit(FatalError);
    
    if(filePtr)
    {
        for(label colInd=0; colInd<columns; colInd++)
        {
            (*filePtr)<<values[colInd];
            if(colInd<columns-1)
                (*filePtr)<<",";
            else
                (*filePtr)<<std::endl;
        }
    }
}
