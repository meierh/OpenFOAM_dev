#include "CrossSectionStructure.H"

Foam::LagrangianMarkerOnCrossSec::LagrangianMarkerOnCrossSec
(
    const LineStructure& structure,
    const fvMesh& mesh,
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* baseRod,
    const scalar markerParameter,
    const CrossSection* baseCrossSec,
    const scalar markerAngle,
    const scalar radiusFrac
):
LagrangianMarker(structure,mesh,rodNumber,baseRod),
markerAngle(markerAngle),
markerRadiusFrac(radiusFrac),
baseCrossSec(baseCrossSec)
{
    this->markerParameter = markerParameter;
    evaluateMarker();
}

void Foam::LagrangianMarkerOnCrossSec::evaluateMarker()
{
    markerPosition = CrossSectionStructure::evaluateRodCircumPos
    (
        baseRod,markerParameter,baseCrossSec,markerAngle,markerRadiusFrac
    );
    markerCell = mesh.findCell(markerPosition);
    computeSupport();
    Pair<vector> h = minMaxNeighbourWidth(directSupport);
    h_plus = h.first();
    h_minus = h.second();
    dilation = dilationFactors(h);
    checkDirectSupport();
    reduceSupport();
    computeCharacLength();
}

Foam::vector Foam::LagrangianMarkerOnCrossSec::getMarkerVelocity() const
{
    return structure.evaluateRodVelocity(rodNumber,markerParameter,markerAngle,markerRadiusFrac);
}

void Foam::LagrangianMarkerOnCrossSec::computeCharacLength()
{
    scalar markerRadiusFrac = this->markerRadiusFrac;
    markerRadiusFrac = std::max(markerRadiusFrac,0.1);
    
    Pair<vector> dMdp_dMdangle = CrossSectionStructure::derivateRodCircumPos(baseRod,markerParameter,baseCrossSec,markerAngle,markerRadiusFrac);
    vector dMdp = dMdp_dMdangle.first();
    vector dMdangle = dMdp_dMdangle.second();
    
    Pair<vector> d2Mdp_d2Mdangle = CrossSectionStructure::derivate2RodCircumPos(baseRod,markerParameter,baseCrossSec,markerAngle,markerRadiusFrac);
    vector d2Mdp = d2Mdp_d2Mdangle.first();
    vector d2Mdangle = d2Mdp_d2Mdangle.second();
    
    scalar len_dMdp = std::sqrt(dMdp & dMdp);
    vector cross_dMdp_d2Mdp = dMdp ^ d2Mdp;
    scalar curvature_dp = std::sqrt(cross_dMdp_d2Mdp & cross_dMdp_d2Mdp)/(len_dMdp*len_dMdp*len_dMdp);
    
    scalar len_dMdangle = std::sqrt(dMdangle & dMdangle);
    vector cross_dMdangle_d2Mdangle = dMdangle ^ d2Mdangle;
    scalar curvature_dangle = std::sqrt(cross_dMdangle_d2Mdangle & cross_dMdangle_d2Mdangle)/(len_dMdangle*len_dMdangle*len_dMdangle);
    
    scalar curvature = std::max(curvature_dangle,curvature_dp);
    
    if(curvature==0)
        FatalErrorInFunction<<"Curvature can not be zero!"<<exit(FatalError);
    
    markerCharLen = 1.0f/curvature;
}
