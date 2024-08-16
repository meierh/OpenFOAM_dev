#ifndef SENSITIVITYINTERACTION_H
#define SENSITIVITYINTERACTION_H
#include "CrossSectionStructure.H"

namespace Foam
{

class SensitivityInteraction : public FieldMarkerStructureInteraction
{
    public:
        explicit SensitivityInteraction
        (
            dynamicRefineFvMesh& mesh,
            LineStructure& structure,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );

    protected:
        template<typename T> void derivedFieldToMarker
        (
            const GeometricField<T,fvPatchField,volMesh>& fieldData,
            DynamicList<T>& markerData,
            const LineStructureParameter& parameter
        );
                
        template<typename T> void derivedMarkerToField
        (
            const DynamicList<T>& markerData,
            GeometricField<T,fvPatchField,volMesh>& fieldData,
            const LineStructureParameter& parameter
        );
        
        template<typename T> void derivedFieldToMarker
        (
            const GeometricField<T,fvPatchField,volMesh>& fieldData,
            DynamicList<T>& markerData,
            const CrossSectionStructureParameter& parameter
        );
                
        template<typename T> void derivedMarkerToField
        (
            const DynamicList<T>& markerData,
            GeometricField<T,fvPatchField,volMesh>& fieldData,
            const CrossSectionStructureParameter& parameter
        );
};
}
#endif
