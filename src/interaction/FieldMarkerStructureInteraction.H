#ifndef FIELDMARKERSTRUCTUREINTERACTION_H
#define FIELDMARKERSTRUCTUREINTERACTION_H
#include "LineStructure.H"
#include "pimpleSingleRegionControl.H"

namespace Foam
{

class FieldMarkerStructureInteraction
{
    public:
        enum markerMeshType {Uniform, NonUniform};
        explicit FieldMarkerStructureInteraction
        (
            const fvMesh& mesh,
            LineStructure& structure,
            const IOdictionary& structureDict,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        virtual ~FieldMarkerStructureInteraction()=default;

        virtual void solve(const pimpleSingleRegionControl&) = 0;

        virtual void store() = 0;
        
        virtual void setToTime(scalar time) = 0;
        
        virtual void recomputeMarkerValues() = 0;
        
        const LineStructure& getStructure() const {return structure;}
        const fvMesh& getMesh() const {return mesh;}
        
    protected:
        ActiveRodMesh::rodMesh* rodMesh;        
        const fvMesh& mesh;
        LineStructure& structure;
        const IOdictionary& structureDict;

        //Map field from fielData to markerData according to markers
        template<typename T> void fieldToMarker
        (
            const Field<T>& fieldData,
            DynamicList<T>& markerData
        );
        
        template<typename T> void integrateFieldOverMarker
        (
            const Field<T>& fieldData,
            DynamicList<T>& markerData,
            std::function<T(const LagrangianMarker*, std::function<T(Pair<label>)>)> integration            
        );
                
        template<typename T> void markerToField
        (
            const DynamicList<T>& markerData,
            Field<T>& fieldData
        );
                
        template<typename T> void integrateMarkersOverField
        (
            const DynamicList<T>& markerData,
            GeometricField<T,fvPatchField,volMesh>& fieldData
        );
        
        void scatterNurbs(std::pair<gsNurbs<scalar>,label> in, std::pair<gsNurbs<scalar>,label>& out);
        
        template<typename T> void gatherMarkerData
        (
            std::tuple<List<scalar>,List<List<T>>,label> in,
            std::tuple<List<scalar>,List<List<T>>,label>& out
        );
        
        template<typename T> T integrateField
        (
            const Field<T>& field
        );
        
        template<typename T> T integrateMarkers
        (
            const DynamicList<T>& markerVals
        );
        
        const scalar h;
        const markerMeshType modusFieldToMarker;
        const markerMeshType modusMarkerToField;
        
        class MarkerInfoFiles
        {
            public:
                MarkerInfoFiles
                (
                    const IOdictionary& structureDict,
                    word dictName,
                    bool masterFile,
                    List<word> header
                );

            protected:
                void write(List<scalar> values);
                        
                bool masterFile;
                std::string fileName;
                std::unique_ptr<std::ofstream> filePtr;
                label columns;
                bool fileActive;
        };
};

#include "FieldMarkerStructureInteraction.T.H"

}
#endif
