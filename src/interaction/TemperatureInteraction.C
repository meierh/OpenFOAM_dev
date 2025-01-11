#include "TemperatureInteraction.H"

Foam::TemperatureInteraction::TemperatureInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volScalarField& input_T,
    volScalarField& output_Tf,
    const IOdictionary& structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
FieldMarkerStructureInteraction(mesh,structure,structureDict,modusFieldToMarker,modusMarkerToField),
input_T(input_T),
output_Tf(output_Tf),
detailedMarkerTemperatureFileObject(structureDict)
{
}

void Foam::TemperatureInteraction::solve
(
    bool firstIteration,
    scalar time,
    bool finalIteration
)
{
    interpolateTemperatureToMarkers();
    computeCouplingHeatingOnMarkers();
    computeRodHeating();
    interpolateHeatingField();
    
    if(finalIteration)
    {
        detailedMarkerTemperatureFileObject.writeSolution(*this);
    }
}

void Foam::TemperatureInteraction::store()
{
    std::tuple<DynamicList<scalar>,DynamicList<scalar>,DynamicList<scalar>>& markerValues = storage[mesh.time().value()];
    
    std::get<0>(markerValues) = markerFluidTemperature;
    std::get<1>(markerValues) = makerCouplingHeating;
    std::get<2>(markerValues) = rodHeating;
}

void Foam::TemperatureInteraction::setToTime(scalar time)
{
    if(storage.find(time)==storage.end())
        FatalErrorInFunction<<"No data exists for "<<time<<exit(FatalError);
    
    std::tuple<DynamicList<scalar>,DynamicList<scalar>,DynamicList<scalar>>& markerValues = storage[time];
    
    markerFluidTemperature = std::get<0>(markerValues);
    makerCouplingHeating = std::get<1>(markerValues);
    rodHeating = std::get<2>(markerValues);
}

void Foam::TemperatureInteraction::recomputeMarkerValues()
{
    interpolateTemperatureToMarkers();
    computeCouplingHeatingOnMarkers();
    computeRodHeating();
    interpolateHeatingField();
}

void Foam::TemperatureInteraction::interpolateTemperatureToMarkers()
{
    fieldToMarker<scalar>(input_T,markerFluidTemperature);
}

void Foam::TemperatureInteraction::computeCouplingHeatingOnMarkers()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(markerFluidTemperature.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of markerTemperature and markers"<<exit(FatalError);
    
    scalar deltaT = mesh.time().deltaTValue();
    
    makerCouplingHeating.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarkerPtr = markers[markerInd];
        scalar markerTemperature = getTemperature(oneMarkerPtr);
        scalar temperature = markerFluidTemperature[markerInd];
        makerCouplingHeating[markerInd] = (markerTemperature-temperature)/deltaT;
    }
}

void Foam::TemperatureInteraction::interpolateHeatingField()
{
    scalar deltaT = mesh.time().deltaTValue();
    markerToField<scalar>(makerCouplingHeating,output_Tf);
    for(scalar Tf : output_Tf)
        Tf = Tf*deltaT;
}

void Foam::TemperatureInteraction::computeRodHeating()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(makerCouplingHeating.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of makerCouplingHeating and markers"<<exit(FatalError);
    
    rodHeating.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarker = markers[markerInd];
        scalar volume = oneMarker->getMarkerVolume();
        scalar cp = 1.0035;
        scalar rho = 1.225;
        
        rodHeating[markerInd] = rho*cp*volume*makerCouplingHeating[markerInd];
    }
}

Foam::scalar Foam::TemperatureInteraction::sumHeating
(
    std::function<bool(LagrangianMarker)> condition
)
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    if(rodHeating.size()!=static_cast<label>(markers.size()))
    {
        Info<<"rodHeating.size():"<<rodHeating.size()<<Foam::endl;
        Info<<"markers.size():"<<markers.size()<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in size of rodHeating and markers"<<exit(FatalError);
    }
    
    scalar result = Foam::zero();
    for(std::size_t i=0; i<markers.size(); i++)
    {
        LagrangianMarker* oneMarker = markers[i];
        if(condition(*oneMarker))
            result += rodHeating[i];
    }
    return result;
}

void Foam::TemperatureInteraction::DetailedMarkerTemperatureFile::writeSolution
(
    const TemperatureInteraction& interaction
)
{
    const std::vector<LagrangianMarker*>& markers = interaction.getStructure().getCollectedMarkers();
    const DynamicList<scalar>& rodHeating = interaction.getRodHeating();
    const volScalarField& T = interaction.getTemperature();
    
    if(rodHeating.size()!=static_cast<label>(markers.size()))
    {
        Info<<"rodHeating.size():"<<rodHeating.size()<<Foam::endl;
        Info<<"markers.size():"<<markers.size()<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in size of rodHeating and markers"<<exit(FatalError);
    }

    volVectorField gradT = fvc::grad(T);
    
    word interpolationType = "cell";
    autoPtr<interpolation<scalar> > Tinterp = interpolation<scalar>::New(interpolationType,T);
    
    List<DynamicList<List<scalar>>> globalLines(Pstream::nProcs());
    DynamicList<List<scalar>>& lines = globalLines[Pstream::myProcNo()];
    for(std::size_t i=0; i<markers.size(); i++)
    {
        LagrangianMarker* oneMarker = markers[i];
    //{"Time","RodInd","Parameter","Angle","RadiusFrac","Px","Py","Pz","Nx","Ny","Nz","cellSpacing","rodHeating","T_1n","T_2n","T_4n","T_8n","interpolated_T","dTdn","dTdn_1n","dTdn_2n","dTdn_4n","dTdn_8n"})
        
        vector P = oneMarker->getMarkerPosition();
        vector N = oneMarker->getMarkerNormal();
        label cellInd = oneMarker->getMarkerCell();
        scalar interpT = Tinterp->interpolate(P,cellInd);
        vector cellGradT = gradT[cellInd];
        scalar dTdn = cellGradT&N;
        scalar cellSpacing = Structure::spacingFromMesh(interaction.getMesh(),cellInd);
        FixedList<vector,4> steps = {N*cellSpacing,2*N*cellSpacing,4*N*cellSpacing,8*N*cellSpacing};
        FixedList<vector,4> points = {P+steps[0],P+steps[1],P+steps[2],P+steps[3]};
        FixedList<label,4> pointsCell;
        for(label i=0; i<4; i++)
            pointsCell[i] = interaction.getStructure().findCell(points[i],cellInd);
        FixedList<scalar,4> pointsT = {0,0,0,0};
        for(label i=0; i<4; i++)
        {
            if(pointsCell[i]>=0 || pointsCell[i]<T.size())
                pointsT[i] = Tinterp->interpolate(points[i],pointsCell[i]);
        }
        FixedList<scalar,4> pointsdTdn;
        for(label i=0; i<4; i++)
        {
            if(pointsCell[i]>=0 || pointsCell[i]<T.size())
            {
                scalar len = std::sqrt(steps[i]&steps[i]);
                pointsdTdn[i] = (pointsT[i]-interpT)/len;
            }
            else
                pointsdTdn[i] = 0;
        }
        
        List<scalar> line =
        {
            interaction.getMesh().time().value(),
            static_cast<scalar>(oneMarker->getRodNumber()),
            oneMarker->getMarkerParameter(),
            oneMarker->getMarkerAngle(),
            oneMarker->getMarkerRadiusFrac(),
            P[0],
            P[1],
            P[2],
            N[0],
            N[1],
            N[2],
            cellSpacing,
            rodHeating[i],
            pointsT[0],
            pointsT[1],
            pointsT[2],
            pointsT[3],
            interpT,
            dTdn,
            pointsdTdn[0],
            pointsdTdn[1],
            pointsdTdn[2],
            pointsdTdn[3]            
        };
        
        lines.append(line);        
    }
    Pstream::gatherList(globalLines);
    if(Pstream::master())
    {
        for(const DynamicList<List<scalar>>& oneProcLines : globalLines)
        {
            for(const List<scalar>& line : oneProcLines)
                write(line);
        }
    }
}
