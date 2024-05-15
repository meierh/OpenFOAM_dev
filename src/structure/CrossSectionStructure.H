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
        scalar operator()(scalar parameter, scalar angle) const;
        std::function<scalar(scalar)> getEvalOnPoint(scalar parameter) const;
        
        scalar deriv_angle(scalar parameter, scalar angle) const;
        std::function<scalar(scalar)> getDerivAngleOnPoint(scalar parameter) const;
                
        scalar domainStart(){return domStart;}
        scalar domainEnd(){return domEnd;}
        
        scalar lowerLimitRadius(scalar parameter) const;
        scalar upperLimitRadius(scalar parameter) const;
        
        const gsNurbs<scalar>& get_a_0() const {return a_0;}
        const std::vector<gsNurbs<scalar>>& get_a_k() const {return a_k;}
        const std::vector<gsNurbs<scalar>>& get_b_k() const {return b_k;}
        uint getNumberCoeffs() const {return numberCoeffs;}

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
            const dynamicRefineFvMesh& mesh,
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* baseRod,
            const scalar markerParameter,
            const CrossSection baseCrossSec,
            const scalar markerAngle,
            const scalar radiusFrac
        ):
        LagrangianMarker(mesh,rodNumber,baseRod,markerParameter),
        markerAngle(markerAngle),
        baseCrossSec(baseCrossSec),
        markerRadiusFrac(radiusFrac)
        {}
        
        inline scalar getMarkerAngle() const {return markerAngle;}
        inline scalar getMarkerRadiusFrac() const {return markerRadiusFrac;}
        inline CrossSection getBaseCrossSec() const {return baseCrossSec;}
        
    private:
        scalar markerAngle;
        scalar markerRadiusFrac;
        CrossSection baseCrossSec;        
};

class CrossSectionStructure : public LineStructure
{
    public:
        explicit CrossSectionStructure
        (
            dynamicRefineFvMesh& mesh,
            std::vector<CrossSection> rodCrossSection
        );
        CrossSectionStructure(const CrossSectionStructure&) = delete;
        void transferMarkers(FieldMarkerStructureInteraction& connector);
        
    protected:
        std::unique_ptr<std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>> createSpacedPointsOnCrossSec
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            const CrossSection& crossSec,            
            scalar spacing
        );
        
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
        
        static vector evaluateRodCircumDerivAngle
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection& oneCrossSec,
            scalar angle,
            scalar radiusFrac=1
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
        
        template<typename T>
        static T integrateCircumwise
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection& crossSec,
            std::function<T(scalar)> function
        );
        
        template<typename T>
        static T integrateCircumwise
        (   
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection& crossSec,
            scalar angleStart,
            scalar angleEnd,
            std::function<T(scalar)> function
        );
        
        static scalar distance
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection& crossSec,
            scalar angleStart,
            scalar angleEnd,
            scalar radiusFrac
        );
        
        static scalar distance
        (
            const LagrangianMarkerOnCrossSec& A,
            const LagrangianMarkerOnCrossSec& B
        );

    
    private:
        std::vector<CrossSection> rodCrossSection;
        
        std::vector<std::unique_ptr<std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>>> rodMarkers;
        
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
