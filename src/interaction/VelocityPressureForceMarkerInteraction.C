#include "VelocityPressureForceMarkerInteraction.H"


Foam::VelocityPressureForceInteraction::VelocityPressureForceInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    volVectorField& input_U,
    volVectorField& output_Uf,
    std::shared_ptr<MeshRefiner> refinement_
):
FieldMarkerStructureInteraction(mesh,structure),
input_U(input_U),
output_Uf(output_Uf),
refinement_(refinement_)
{
}

void Foam::VelocityPressureForceInteraction::solve()
{
    interpolateFluidVelocityToMarkers();
    computeCouplingForceOnMarkers();
    computeRodForceMoment();
    interpolateFluidForceField();
}

void Foam::VelocityPressureForceInteraction::store()
{
    std::tuple<DynamicList<vector>,DynamicList<vector>,DynamicList<vector>,DynamicList<vector>>& markerValues = storage[mesh.time().value()];
    
    std::get<0>(markerValues) = markerFluidVelocity;
    std::get<1>(markerValues) = makerCouplingForce;
    std::get<2>(markerValues) = rodForce;
    std::get<3>(markerValues) = rodMoment;
}

void Foam::VelocityPressureForceInteraction::setToTime(scalar time)
{
    if(storage.find(time)==storage.end())
        FatalErrorInFunction<<"No data exists for "<<time<<exit(FatalError);
    
    std::tuple<DynamicList<vector>,DynamicList<vector>,DynamicList<vector>,DynamicList<vector>>& markerValues = storage[time];
    
    markerFluidVelocity = std::get<0>(markerValues);
    markerFluidVelocity = std::get<1>(markerValues);
    rodForce = std::get<2>(markerValues);
    rodMoment = std::get<3>(markerValues);
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
    
    makerCouplingForce.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarkerPtr = markers[markerInd];
        vector markerVelocity = getVelocity(oneMarkerPtr);
        vector velocity = markerFluidVelocity[markerInd];
        makerCouplingForce[markerInd] = (markerVelocity-velocity)/deltaT;
    }
}

void Foam::VelocityPressureForceInteraction::interpolateFluidForceField()
{
    markerToField<vector>(makerCouplingForce,output_Uf);
}

void Foam::VelocityPressureForceInteraction::moveRodsAndMarkers()
{
    
    std::vector<scalar> knotVec = {0,0,0,0.5,1,1,1};
    gsMatrix<scalar> P(4,3);
    for(int r=0; r<P.rows(); r++)
    {
        P(r,0) = r;
        P(r,1) = 0;
        P(r,2) = 0;
    }
    gsMatrix<scalar> w(4,1);
    for(int i=0; i<4; i++)
        w.at(i) = 1;    
    
    gsNurbs<scalar> testRod(knotVec,w,P);
    gsMatrix<scalar> u(1,1);
    u.at(0) = 0.3;
    
    gsMatrix<scalar> res;
    testRod.eval_into(u,res);
    std::cout<<res<<std::endl;
    
    gsMatrix<scalar> fittedCoeffs;
    
    List<vector> points(10);
    for(vector& pnts : points)
        pnts = vector(0,0,0);
    
    std::cout<<"testRod.coefs():"<<testRod.coefs()<<std::endl;

    FatalErrorInFunction<<"Temp stop"<<exit(FatalError);
    Structure::fitNurbsCoeffsToPoints(points,testRod,fittedCoeffs);
    
    
    FatalErrorInFunction<<"Temp stop"<<exit(FatalError);
    
    
    List<bool> prevRodInMesh = structure.getRodInMesh();
    std::unique_ptr<List<List<vector>>> defPtr = getDeformation();
    if(defPtr)
    {
        structure.setDeformation(*defPtr);
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
        if(refinement_)
        {
            refinement_->refineMeshOnStaticMarkers();
            refinement_->refineMeshAndMarkers();
        }
        structure.finalizeMarkers();    
    }
}

Foam::vector Foam::VelocityPressureForceInteraction::sumForces
(
    std::function<bool(LagrangianMarker)> condition
)
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();
    
    if(rodForce.size()!=static_cast<label>(markers.size()))
    {
        Info<<"rodForce.size():"<<rodForce.size()<<Foam::endl;
        Info<<"markers.size():"<<markers.size()<<Foam::endl;
        FatalErrorInFunction<<"Mismatch in size of rodForce and markers"<<exit(FatalError);
    }
    
    vector result = Foam::zero();
    for(std::size_t i=0; i<markers.size(); i++)
    {
        LagrangianMarker* oneMarker = markers[i];
        if(condition(*oneMarker))
            result += rodForce[i];
    }
    Pstream::gather(result,std::plus<vector>());
    Pstream::scatter(result);
    
    return result;
}

void Foam::VelocityPressureForceInteraction::computeRodForceMoment()
{
    const std::vector<LagrangianMarker*>& markers = structure.getCollectedMarkers();

    if(makerCouplingForce.size()!=static_cast<label>(markers.size()))
        FatalErrorInFunction<<"Mismatch in size of makerCouplingForce and markers"<<exit(FatalError);
    
    rodForce.resize(markers.size());
    rodMoment.resize(markers.size());
    for(std::size_t markerInd=0; markerInd<markers.size(); markerInd++)
    {
        LagrangianMarker* oneMarker = markers[markerInd];
        scalar volume = oneMarker->getMarkerVolume();
        scalar rho = 1.225;
        
        rodForce[markerInd] = rho*volume*makerCouplingForce[markerInd];
        
        vector basePnt;
        Structure::rodEval(oneMarker->getBaseRod(),oneMarker->getMarkerParameter(),basePnt);
        vector vectorToMarker = oneMarker->getMarkerPosition()-basePnt;
        vector momentum = vectorToMarker^makerCouplingForce[markerInd];
        rodMoment[markerInd] = rho*momentum*volume;
    }
}
