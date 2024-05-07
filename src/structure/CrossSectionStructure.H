#ifndef CROSSSECTIONSTRUCTURE_H
#define CROSSSECTIONSTRUCTURE_H

#include "LineStructure.H"

namespace Foam
{
class CrossSection
{
    public:
        explicit CrossSection
        (
            scalar radius
        );
        
        explicit CrossSection
        (
            scalar a_0,
            std::vector<scalar> a_k,
            std::vector<scalar> b_k,
            scalar phaseShift=0
        );
        explicit CrossSection
        (
            scalar a_0,
            std::vector<scalar> a_k,
            std::vector<scalar> b_k,
            gsNurbs<scalar> phaseShift
        );
        
        explicit CrossSection
        (
            gsNurbs<scalar> a_0,
            std::vector<gsNurbs<scalar>> a_k,
            std::vector<gsNurbs<scalar>> b_k,
            gsNurbs<scalar> phaseShift
        );
        scalar operator()(scalar parameter, scalar rad) const;
        std::function<scalar(scalar)> getEvalOnPoint(scalar parameter) const;
        
        scalar domainStart(){return domStart;}
        scalar domainEnd(){return domEnd;}

    private:
        std::tuple<gsNurbs<scalar>,std::vector<gsNurbs<scalar>>,std::vector<gsNurbs<scalar>>,gsNurbs<scalar>> constCrossSec
        (
            scalar a_0,
            std::vector<scalar> a_k,
            std::vector<scalar> b_k,
            scalar phaseShift = 0
        );
        std::tuple<gsNurbs<scalar>,std::vector<gsNurbs<scalar>>,std::vector<gsNurbs<scalar>>,gsNurbs<scalar>> twistedConstCrossSec
        (
            scalar a_0,
            std::vector<scalar> a_k,
            std::vector<scalar> b_k,
            gsNurbs<scalar> phaseShift
        );
        void init
        (
            gsNurbs<scalar> a_0,
            std::vector<gsNurbs<scalar>> a_k,
            std::vector<gsNurbs<scalar>> b_k,
            gsNurbs<scalar> phaseShift
        );
        gsNurbs<scalar> createConstNurbs(scalar coeff) const;
        
        gsNurbs<scalar> a_0;
        std::vector<gsNurbs<scalar>> a_k;
        std::vector<gsNurbs<scalar>> b_k;
        gsNurbs<scalar> phaseShift;
        uint numberCoeffs;
        scalar domStart;
        scalar domEnd;
};

class LagrangianMarkerOnCrossSec : public LagrangianMarker
{
    public:
        explicit LagrangianMarkerOnCrossSec
        (
            const scalar markerParameter,
            const vector markerPosition,
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* baseRod,
            const scalar markerAngle,
            const CrossSection baseCrossSec,
            const scalar radiusFrac
        ):
        LagrangianMarker(markerParameter,markerPosition,rodNumber,baseRod),
        markerAngle(markerAngle),
        baseCrossSec(baseCrossSec),
        markerRadiusFrac(radiusFrac)
        {}
        
        const scalar markerAngle;
        const scalar markerRadiusFrac;
        const CrossSection baseCrossSec;        
};

class CrossSectionStructure : public LineStructure
{
    public:
        explicit CrossSectionStructure
        (
            dynamicRefineFvMesh& mesh,
            /*
            const dimensionedScalar& alpha,
            const volScalarField& T,
            const volScalarField& p,
            const volVectorField& U,
            const dimensionedScalar nu,
            */
            std::vector<CrossSection> rodCrossSection
        );
        CrossSectionStructure(const CrossSectionStructure&) = delete;
        void transferMarkers(FieldMarkerStructureInteraction& connector);
        
    protected:
        scalar evaluateCircumArcLen
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameterA,
            scalar parameterB,
            const CrossSection& oneCrossSec,
            scalar angleA,
            scalar angleB,
            scalar radiusFracA,
            scalar radiusFracB
        );
        static vector evaluateRodCircumPos
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection& oneCrossSec,
            scalar angle,
            scalar radiusFrac,
            scalar var_para = 0,
            scalar var_radius = 0
        );
        bool doSubdivision
        (
            const std::list<LagrangianMarkerOnCrossSec>& smallerSide,
            const std::list<LagrangianMarkerOnCrossSec>& largerSide
        );
        bool doSubdivisionCircumferential
        (
            const LagrangianMarkerOnCrossSec& smallerSide,
            const LagrangianMarkerOnCrossSec& largerSide
        );
        
        static scalar restrictAngle(scalar angle);
        static scalar distance(const LagrangianMarkerOnCrossSec& A, const LagrangianMarkerOnCrossSec& B);

    
    private:
        std::vector<CrossSection> rodCrossSection;
        
        std::vector<std::unique_ptr<std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>>> rodMarkers;
        
        LagrangianMarkerOnCrossSec createLagrangianMarkerOnCrossSec
        (
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection& oneCrossSec,
            scalar angle,
            scalar radius
        );
        std::unique_ptr<std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>> constructMarkerSet
        (
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* oneRod,
            const CrossSection& oneCrossSec
        );
        
        void constructMarkerSetCircumferential
        (
            const label rodNumber,
            const ActiveRodMesh::rodCosserat*,
            const CrossSection&,
            std::list<LagrangianMarkerOnCrossSec>&,
            scalar radiusFrac = 1
        );
        
        void constructMarkerSetRadial
        (
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* oneRod,
            const scalar parameter,
            const CrossSection& oneCrossSec,
            std::list<std::list<LagrangianMarkerOnCrossSec>>& radialMarkers
        );
        
        std::list<LagrangianMarkerOnCrossSec> createInitialCircumMarkers
        (
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection& oneCrossSec,
            scalar radiusFrac = 1
        );
};
}
#endif
