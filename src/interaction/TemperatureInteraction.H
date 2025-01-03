#ifndef TEMPERATUREINTERACTION_H
#define TEMPERATUREINTERACTION_H
#include "FieldMarkerStructureInteraction.H"
#include "fvc.H"
#include "interpolation.H"

namespace Foam
{
    
class TemperatureInteraction : public FieldMarkerStructureInteraction
{
    public:
        explicit TemperatureInteraction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            volScalarField& input_T,
            volScalarField& output_Tf,
            const IOdictionary& structureDict,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        
        void solve(bool finalIteration) override;
        
        void store() override;
        void setToTime(scalar time) override;
        
        const volScalarField& getTemperatureField() const {return input_T;}
        const volScalarField& getTemperatureForcingField() const {return output_Tf;}
        const DynamicList<scalar>& getMarkerCouplingHeating() const {return makerCouplingHeating;}
        
        void recomputeMarkerValues() override;

    protected:
        void interpolateTemperatureToMarkers();
        void computeCouplingHeatingOnMarkers();
        void computeRodHeating();
        void interpolateHeatingField();
        
        scalar sumHeating
        (
            std::function<bool(LagrangianMarker)> condition = 
                                [](LagrangianMarker marker){return true;}
        );
        
        virtual scalar getTemperature(const LagrangianMarker* marker) = 0;
        
        const DynamicList<scalar>& getRodHeating() const {return rodHeating;}
        const volScalarField& getTemperature() const {return input_T;}

    private:
        volScalarField& input_T;
        volScalarField& output_Tf;
        
        DynamicList<scalar> markerFluidTemperature;
        DynamicList<scalar> makerCouplingHeating;
        DynamicList<scalar> rodHeating;
        
        std::map<scalar,std::tuple<DynamicList<scalar>,DynamicList<scalar>,DynamicList<scalar>>> storage;
        
        //void assignForceOnRod();
        
        class DetailedMarkerTemperatureFile : public MarkerInfoFiles
        {
            public:
                DetailedMarkerTemperatureFile
                (
                    const IOdictionary& structureDict
                ):
                MarkerInfoFiles
                (structureDict,"recordDetailedMarkerTemperature",true,{"Time","RodInd","Parameter","Angle","RadiusFrac","Px","Py","Pz","Nx","Ny","Nz","cellSpacing","rodHeating","T_1n","T_2n","T_4n","T_8n","interpolated_T","dTdn","dTdn_1n","dTdn_2n","dTdn_4n","dTdn_8n"}
                )
                {}
            
                void writeSolution(const TemperatureInteraction& interaction);
        };
        DetailedMarkerTemperatureFile detailedMarkerTemperatureFileObject;
};

}
#endif
