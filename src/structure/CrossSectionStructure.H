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
            std::vector<scalar> b_k
        );        
        explicit CrossSection
        (
            gsNurbs<scalar> a_0,
            std::vector<gsNurbs<scalar>> a_k,
            std::vector<gsNurbs<scalar>> b_k
        );
        scalar operator()(scalar parameter, scalar rad) const;
        std::function<scalar(scalar)> getEvalOnPoint(scalar parameter) const;
        
        scalar domainStart(){return domStart;}
        scalar domainEnd(){return domEnd;}

    private:
        std::tuple<gsNurbs<scalar>,std::vector<gsNurbs<scalar>>,std::vector<gsNurbs<scalar>>> constCrossSec
        (
            scalar a_0,
            std::vector<scalar> a_k,
            std::vector<scalar> b_k
        );
        void init
        (
            gsNurbs<scalar> a_0,
            std::vector<gsNurbs<scalar>> a_k,
            std::vector<gsNurbs<scalar>> b_k
        );
        gsNurbs<scalar> a_0;
        std::vector<gsNurbs<scalar>> a_k;
        std::vector<gsNurbs<scalar>> b_k;
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
            const ActiveRodMesh::rodCosserat* baseRod,
            const scalar markerAngle,
            const CrossSection baseCrossSec
        ):
        LagrangianMarker(markerParameter,markerPosition,baseRod),
        markerAngle(markerAngle),
        baseCrossSec(baseCrossSec)
        {}
        
        const scalar markerAngle;
        const CrossSection baseCrossSec;        
};

class CrossSectionStructure : public LineStructure
{
    public:
        explicit CrossSectionStructure
        (
            dynamicRefineFvMesh& mesh,
            const dimensionedScalar& alpha,
            const volScalarField& T,
            const volScalarField& p,
            const volVectorField& U,
            const dimensionedScalar nu,
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
            scalar angleB
        );
        static vector evaluateRodCircumPos
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection& oneCrossSec,
            scalar angle
        );
        bool doSubdivision
        (
            const std::list<LagrangianMarkerOnCrossSec>& smallerSide,
            const std::list<LagrangianMarkerOnCrossSec>& largerSide
        );
        bool doSubdivisionAngle
        (
            const LagrangianMarkerOnCrossSec& smallerSide,
            const LagrangianMarkerOnCrossSec& largerSide
        );
        
        static scalar distance(const LagrangianMarkerOnCrossSec& A, const LagrangianMarkerOnCrossSec& B);

    
    private:
        std::vector<CrossSection> rodCrossSection;
        
        std::vector<std::unique_ptr<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>> rodMarkers;
        
        LagrangianMarkerOnCrossSec createLagrangianMarkerOnCrossSec
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection& oneCrossSec,
            scalar angle
        );
        std::unique_ptr<std::vector<std::vector<LagrangianMarkerOnCrossSec>>> constructMarkerSet
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            const CrossSection& oneCrossSec
        );
        
        void constructMarkerSetCircumferential
        (
            const ActiveRodMesh::rodCosserat*,
            const CrossSection&,
            std::list<LagrangianMarkerOnCrossSec>& 
        );
        
        std::list<LagrangianMarkerOnCrossSec> createInitialCircumMarkers
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection& oneCrossSec
        );
};
}
#endif
