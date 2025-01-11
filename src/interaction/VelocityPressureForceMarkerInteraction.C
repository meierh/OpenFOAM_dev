#include "VelocityPressureForceMarkerInteraction.H"


Foam::VelocityPressureForceInteraction::VelocityPressureForceInteraction
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
FieldMarkerStructureInteraction(mesh,structure,structureDict,modusFieldToMarker,modusMarkerToField),
input_U(input_U),
output_Uf(output_Uf),
refinement_(refinement_),
sumMarkerForceFileObject(structureDict),
sumMarkerMomentFileObject(structureDict),
detailedMarkerForceFileObject(structureDict)
{  
    ITstream interiorForcingStream = structureDict.lookup("interiorForcing");
    token interiorForcingToken;
    interiorForcingStream.read(interiorForcingToken);
    if(!interiorForcingToken.isWord())
    {
        FatalErrorInFunction<<"Invalid entry in constant/structureDict/interiorForcingStream -- must be  word"<<exit(FatalError);
    }
    word interiorForcingCmd = interiorForcingToken.wordToken();
    if(interiorForcingCmd=="yes")
        interiorForcingActive = true;
    else if(interiorForcingCmd=="no")
        interiorForcingActive = false;
    else
        FatalErrorInFunction<<"Invalid entry in constant/structureDict/interiorForcingStream -- must be  {yes,no}"<<exit(FatalError);
}

void Foam::VelocityPressureForceInteraction::solve
(
    bool firstIteration,
    scalar time,
    bool finalIteration
)
{   
    interpolateFluidVelocityToMarkers();
    computeCouplingForceOnMarkers();
    computeRodForceMoment();
    interpolateFluidForceField();
    if(interiorForcingActive)
    {
        interiorForcing(time,reconstructInterior(firstIteration,time,finalIteration));
    }
    
    if(finalIteration)
    {
        sumMarkerForceFileObject.writeSolution(*this);
        sumMarkerMomentFileObject.writeSolution(*this);
        detailedMarkerForceFileObject.writeSolution(*this);
    }
}

void Foam::VelocityPressureForceInteraction::store()
{
    std::tuple<DynamicList<vector>,DynamicList<vector>,DynamicList<vector>,DynamicList<vector>>& markerValues = storage[mesh.time().value()];
    
    std::get<0>(markerValues) = markerFluidVelocity;
    std::get<1>(markerValues) = markerCouplingForce;
    std::get<2>(markerValues) = rodForce;
    std::get<3>(markerValues) = rodMoment;
}

void Foam::VelocityPressureForceInteraction::setToTime(scalar time)
{
    if(storage.find(time)==storage.end())
        FatalErrorInFunction<<"No data exists for "<<time<<exit(FatalError);
    
    std::tuple<DynamicList<vector>,DynamicList<vector>,DynamicList<vector>,DynamicList<vector>>& markerValues = storage[time];
    
    markerFluidVelocity = std::get<0>(markerValues);
    markerCouplingForce = std::get<1>(markerValues);
    rodForce = std::get<2>(markerValues);
    rodMoment = std::get<3>(markerValues);
}

void Foam::VelocityPressureForceInteraction::recomputeMarkerValues()
{
    interpolateFluidVelocityToMarkers();
    computeCouplingForceOnMarkers();
    computeRodForceMoment();
    interpolateFluidForceField();
}

void Foam::VelocityPressureForceInteraction::interpolateFluidVelocityToMarkers()
{
    fieldToMarker<vector>(input_U,markerFluidVelocity);
}

void Foam::VelocityPressureForceInteraction::computeCouplingForceOnMarkers()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    
    if(markerFluidVelocity.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of markerFluidVelocity and markers"<<exit(FatalError);
    
    scalar deltaT = mesh.time().deltaTValue();
    Info<<"Foam::VelocityPressureForceInteraction::computeCouplingForceOnMarkers: "<<deltaT<<Foam::nl;
    
    markerCouplingForce.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarkerPtr = markers[markerInd];
        vector markerVelocity = getVelocity(oneMarkerPtr);
        vector velocity = markerFluidVelocity[markerInd];
        markerCouplingForce[markerInd] = (markerVelocity-velocity)/deltaT;
    }
}

void Foam::VelocityPressureForceInteraction::interpolateFluidForceField()
{
    markerToField<vector>(markerCouplingForce,output_Uf);
}

void Foam::VelocityPressureForceInteraction::interiorForcing(scalar time,bool reconstruct)
{
    if(necessaryReconstruct)
    {
        reconstruct = true;
        necessaryReconstruct = false;
    }
    scalar deltaT = mesh.time().deltaTValue();
    const List<std::tuple<label,label,scalar,scalar,scalar>>& interiorCells = structure.getInteriorCells(time,reconstruct);
    for(const std::tuple<label,label,scalar,scalar,scalar>& cell : interiorCells)
    {
        label cellInd = std::get<0>(cell);
        label rodInd = std::get<1>(cell);
        scalar para = std::get<2>(cell);
        scalar angle = std::get<3>(cell);
        scalar radiusFrac = std::get<4>(cell);
        vector cellVelocityGoal = getCellVelocity(rodInd,para,angle,radiusFrac);
        vector cellVelocity = input_U[cellInd];
        vector forcingVelocity = (cellVelocityGoal-cellVelocity)/deltaT;
        output_Uf[cellInd] = forcingVelocity;
    }
    Info<<"interiorForcing:"<<interiorCells.size()<<Foam::endl;
}

void Foam::VelocityPressureForceInteraction::moveMarkers()
{  
    List<bool> prevRodInMesh = structure.getRodInMesh();
    structure.moveMarkersOnRodMovement();
    const List<bool>& rodInMesh = structure.getRodInMesh();
    for(label rodNumber=0; rodNumber<prevRodInMesh.size(); rodNumber++)
    {
        if(!prevRodInMesh[rodNumber])
        {
            if(rodInMesh[rodNumber])
            {
                structure.createMarkersOnRod(rodNumber);
            }
        }
    }
}

Foam::vector Foam::VelocityPressureForceInteraction::summedForces() const
{
    vector sumForces = Foam::zero();
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    for(label cellInd=0; cellInd<mesh.cells().size(); cellInd++)
    {
        scalar volume = cells[cellInd].mag(points,faces);
        sumForces += output_Uf[cellInd]*volume;
    }
    Pstream::gather(sumForces,std::plus<vector>());
    Pstream::scatter(sumForces);
    return sumForces;
}

Foam::vector Foam::VelocityPressureForceInteraction::summedMoments() const
{
    vector sumForces = Foam::zero();
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    for(label cellInd=0; cellInd<mesh.cells().size(); cellInd++)
    {
        scalar volume = cells[cellInd].mag(points,faces);
        vector centre = cells[cellInd].centre(points,faces);
        sumForces += centre ^ (output_Uf[cellInd]*volume);
    }
    Pstream::gather(sumForces,std::plus<vector>());
    Pstream::scatter(sumForces);
    return sumForces;
}

void Foam::VelocityPressureForceInteraction::meshMarkerAdaptation()
{
    if(refinement_)
    {
        Info<<"|||||||||||||||||||||||||Do refinement|||||||||||||||||||||||||"<<Foam::endl;
        //bool meshRefined = refinement_->refineMeshOnStaticMarkers();
        //meshRefined = refinement_->refineMeshAndMarkers(meshRefined);
        refinement_->refineMeshAndMarkers();
        Info<<"||||||||||||||||||||||||Done refinement||||||||||||||||||||||||"<<Foam::endl;
        necessaryReconstruct = true;
    }
    structure.finalizeMarkers(true);
}

void Foam::VelocityPressureForceInteraction::computeRodForceMoment()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();

    if(markerCouplingForce.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of makerCouplingForce and markers"<<exit(FatalError);
    
    rodForce.resize(markers.size());
    rodMoment.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarker = markers[markerInd];
        scalar volume = oneMarker->getMarkerVolume();
        scalar rho = 1.225;
        
        rodForce[markerInd] = rho*volume*markerCouplingForce[markerInd];
        
        vector basePnt;
        Structure::rodEval(oneMarker->getBaseRod(),oneMarker->getMarkerParameter(),basePnt);
        vector vectorToMarker = oneMarker->getMarkerPosition()-basePnt;
        vector momentum = vectorToMarker^markerCouplingForce[markerInd];
        rodMoment[markerInd] = rho*momentum*volume;
    }
}
            
void Foam::VelocityPressureForceInteraction::SumMarkerForceFile::writeSolution
(
    const VelocityPressureForceInteraction& interaction
)
{
    if(fileActive)
    {
        label nR = interaction.getStructure().getNumberRods();
        for(label rodInd=0; rodInd<nR; rodInd++)
        {
            vector summedForces = interaction.summedForces();
            if(Pstream::master())
            {
                write({interaction.getMesh().time().value(),static_cast<scalar>(rodInd),summedForces[0],summedForces[1],summedForces[2]});
            }
        }
    }
}
            
void Foam::VelocityPressureForceInteraction::SumMarkerMomentFile::writeSolution
(
    const VelocityPressureForceInteraction& interaction
)
{
    if(fileActive)
    {
        label nR = interaction.getStructure().getNumberRods();
        for(label rodInd=0; rodInd<nR; rodInd++)
        {
            vector summedMoments = interaction.summedMoments();
            if(Pstream::master())
            {
                write({interaction.getMesh().time().value(),static_cast<scalar>(rodInd),summedMoments[0],summedMoments[1],summedMoments[2]});
            }
        }
    }
}
  
void Foam::VelocityPressureForceInteraction::DetailedMarkerForceFile::writeSolution
(
    const VelocityPressureForceInteraction& interaction
)
{
    if(fileActive)
    {
        const std::vector<LagrangianMarker*>& markers = interaction.getStructure().getCollectedMarkers();
        const DynamicList<vector>& rodForce = interaction.getRodForce();

        if(rodForce.size()!=static_cast<label>(markers.size()))
        {
            Info<<"rodForce.size():"<<rodForce.size()<<Foam::endl;
            Info<<"markers.size():"<<markers.size()<<Foam::endl;
            FatalErrorInFunction<<"Mismatch in size of rodForce and markers"<<exit(FatalError);
        }

        List<DynamicList<List<scalar>>> globalLines(Pstream::nProcs());
        DynamicList<List<scalar>>& lines = globalLines[Pstream::myProcNo()];
        for(std::size_t i=0; i<markers.size(); i++)
        {
            LagrangianMarker* oneMarker = markers[i];
        //{"Time","RodInd","Parameter","Angle","RadiusFrac","Nx","Ny","Nz","Fx","Fy","Fz"})

            vector P = oneMarker->getMarkerPosition();
            vector N = oneMarker->getMarkerNormal();

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
                rodForce[i][0],
                rodForce[i][1],
                rodForce[i][2]
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
}
