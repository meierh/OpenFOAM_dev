#ifndef CROSSSECTIONSTRUCTURE_H
#define CROSSSECTIONSTRUCTURE_H

#include "LineStructure.H"
#include <memory>

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
        
    private:
        std::vector<std::vector<gsNurbs<scalar>*>*> coeffDerivedCurves;
        void delete_coeffDerivedCurves();
        const gsNurbs<scalar>* getCurvePtr(label fourierCoeffNumber) const;
        gsNurbs<scalar>* getCurvePtr(label fourierCoeffNumber);

    public:
        label numberFourierCoeff() const;
        label numberNurbsCoeffs(label fourierCoeffNumber) const;
        scalar evalRadiusDerivCoeff(label fourierCoeffNumber, label derivCoeffNumber, scalar parameter);
        
        void setNurbsCoeff(label fourierCoeffNumber, label derivCoeffNumber, scalar value);
                
        scalar domainStart(){return domStart;}
        scalar domainEnd(){return domEnd;}
        
        scalar lowerLimitRadius(scalar parameter) const;
        scalar upperLimitRadius(scalar parameter) const;
        std::pair<scalar,scalar> nurbsBounds(gsNurbs<scalar> curve) const;
        std::pair<scalar,scalar> nurbsBounds(gsNurbs<scalar> curve, scalar start, scalar end) const;
        std::pair<scalar,scalar> radiusBounds() const;
        std::pair<scalar,scalar> radiusBounds(scalar start, scalar end) const;
        
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
            const Structure& structure,
            const fvMesh& mesh,
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* baseRod,
            const scalar markerParameter,
            const CrossSection* baseCrossSec,
            const scalar markerAngle,
            const scalar radiusFrac
        );
        
        inline scalar getMarkerAngle() const override {return markerAngle;}
        inline scalar getMarkerRadiusFrac() const override {return markerRadiusFrac;}
        inline const CrossSection* getBaseCrossSec() const {return baseCrossSec;}
        void evaluateMarker();
        
        vector getMarkerVelocity() const override;
        
    private:
        scalar markerAngle;
        scalar markerRadiusFrac;
        const CrossSection* baseCrossSec;        
};

class CrossSectionStructure : public LineStructure
{
    public:
        explicit CrossSectionStructure
        (
            const fvMesh& mesh,
            std::vector<CrossSection> rodCrossSection,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        explicit CrossSectionStructure
        (
            const fvMesh& mesh,
            const IOdictionary& stuctureDict,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );

        CrossSectionStructure(const CrossSectionStructure&) = delete;
        //void transferMarkers(FieldMarkerStructureInteraction& connector);
        
    // ----------------- Moving Rod update methods -----------------
        vector evaluateRodVelocity(label rodNumber, scalar parameter, scalar angle, scalar radiusFrac) const override;        
    // ----------------- ------------------------- -----------------
        
        void to_string() override;
        
    // ----------------- Optimization related methods -----------------
        void setCrossSecParameters(label rodNumber, label fourierCoeffNumber, label derivCoeffNumber, scalar value);
        vector dXdParam(const LagrangianMarker* marker,const Parameter& par) override;
        vector dXdParam(label rodNumber,scalar rodParameter,scalar angle,scalar radiusFrac,const Parameter& par) override;
    // ----------------- ---------------------------- -----------------
        
    private:
        std::vector<CrossSection> createCrossSectionsFromDict(const IOdictionary& stuctureDict);

        void check() override;
        
        void createSpacingPoints() override;
        void createSpacedPointsOnRod
        (
            label rodNumber,
            scalar spacing
        ) override;
        void createMarkersFromSpacedPointsOnRod
        (
            label rodNumber
        ) override;
        void createSpacedPointsOnCrossSec
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection* oneCrossSec,
            scalar radFrac,
            scalar spacing,
            std::vector<scalar>& angleData
        );
        
        void refineMarkersOnRod
        (
            label rodNumber,
            std::pair<bool,scalar> forcedSpacing = {false,-1}
        ) override;
        void refineCircumferential
        (
            std::list<LagrangianMarkerOnCrossSec>& circumMarkers,
            std::pair<bool,scalar> refineSpacing = {false,0}
        );
        void refineRadial
        (
            std::list<std::list<LagrangianMarkerOnCrossSec>>& radialMarkers,
            scalar initialSpacing,
            std::pair<bool,scalar> refineSpacing = {false,0}
        );
        void refineTangential
        (
            std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>& tangMarkers,
            scalar initialSpacing,
            std::pair<bool,scalar> refineSpacing = {false,0}
        );
        
        void setMarkerVolumeOnRod
        (
            label rodNumber
        ) override;
        
        void evaluateMarkerMeshRelation() override;
        void evaluateMarkerMeshRelation(std::list<LagrangianMarkerOnCrossSec>& markerList);
        
        void reduceMarkers() override;
        
        BoundingBox computeBox
        (
            label rodNumber
        ) override;
        BoundingBox computeBox
        (
            label rodNumber,
            scalar parStart,
            scalar parEnd
        ) override;
        scalar characteristicSize
        (
            label rodNumber,
            scalar par
        ) override;
        
        void removeOverlapMarkers() override;
        
        //void transferMarkers();
        
        void collectMarkers() override;
        
        /*
        std::unique_ptr<std::vector<LagrangianMarkerOnCrossSec>> constructMarkerSet
        (
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* oneRod,
            const CrossSection* oneCrossSec,
            scalar initialSpacing,
            bool reInitialize = false,
            std::pair<bool,scalar> refineSpacing = {false,0}
        );
        */
        

        /*
        std::unique_ptr<std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>> createSpacedPointsOnCrossSecRod
        (
            label rodNumber,
            const CrossSection* crossSec,            
            scalar spacing
        );
        */
        
        /*
        void computeMarkerVolume
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            const CrossSection* oneCrossSec,
            std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>& markers
        );
        */
        
        std::unique_ptr<std::vector<LagrangianMarkerOnCrossSec>> reduceSurplusMarkers
        (
            std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>& markers
        );
        
        scalar evaluateCircumArcLen
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameterA,
            scalar parameterB,
            const CrossSection* oneCrossSec,
            scalar angleA,
            scalar angleB,
            scalar radiusFracA,
            scalar radiusFracB
        );

    public:
        static vector evaluateRodCircumPos
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection* oneCrossSec,
            scalar angle,
            scalar radiusFrac,
            scalar var_para = 0,
            scalar var_radius = 0
        );
        
        vector evalRodDerivCoeff
        (
            label rodNumber,
            label derivCoeffNumber,
            scalar parameter,
            scalar angle,
            scalar radiusFrac = 1
        );
        
    private:
        static vector evaluateRodCircumDerivAngle
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection* oneCrossSec,
            scalar angle,
            scalar radiusFrac=1
        );
                
        static scalar restrictAngle(scalar angle);
        
        template<typename T>
        static T integrateCircumwise
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection* crossSec,
            std::function<T(scalar)> function
        );
        
        template<typename T>
        static T integrateCircumwise
        (   
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection* crossSec,
            scalar angleStart,
            scalar angleEnd,
            std::function<T(scalar)> function
        );
        
        static scalar distance
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection* crossSec,
            scalar angleStart,
            scalar angleEnd,
            scalar radiusFrac
        );
        
        static scalar distance
        (
            const LagrangianMarkerOnCrossSec& A,
            const LagrangianMarkerOnCrossSec& B
        );

        std::vector<CrossSection> rodCrossSection;
    public:
        const std::vector<CrossSection>& getRodCrossSections() const {return rodCrossSection;}
    
    private:
        //[rod] -> {para , [] -> {rad , [] -> angle}}
        std::vector<std::unique_ptr<std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>>> initialRodPoints;
        std::vector<std::unique_ptr<std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>>> rodMarkersList;

        /*
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
        */
    
        /*
        void constructMarkerSetCircumferential
        (
            const label rodNumber,
            const ActiveRodMesh::rodCosserat*,
            const CrossSection*,
            std::list<LagrangianMarkerOnCrossSec>&,
            scalar radiusFrac = 1
        );
        */
        
        /*
        void constructMarkerSetRadial
        (
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* oneRod,
            const scalar parameter,
            const CrossSection* oneCrossSec,
            std::list<std::list<LagrangianMarkerOnCrossSec>>& radialMarkers
        );
        */
        
        /*
        std::list<LagrangianMarkerOnCrossSec> createInitialCircumMarkers
        (
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter,
            const CrossSection* oneCrossSec,
            scalar radiusFrac = 1
        );
        */
};

class CrossSectionStructureParameters : public LineStructureParameters
{
    public:
        CrossSectionStructureParameters();
        virtual ~CrossSectionStructureParameters() = default;
        void collectParameters(const CrossSectionStructure* structure);
        
    private:
        std::vector<Parameter> parameters;
};

}
#endif
