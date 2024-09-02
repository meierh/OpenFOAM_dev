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
            fvMesh& mesh,
            LineStructure& structure,
            markerMeshType modusFieldToMarker = markerMeshType::NonUniform,
            markerMeshType modusMarkerToField = markerMeshType::NonUniform
        );
        virtual ~SensitivityInteraction()=default;

        virtual scalar computeSensitivity(const Parameter& para) = 0;

    protected:
        template<typename T> void deriveParamFieldToMarker
        (
            const GeometricField<T,fvPatchField,volMesh>& fieldData,
            DynamicList<T>& markerData,
            const Parameter& par
        );
                
        template<typename T> void deriveParamMarkerToField
        (
            const DynamicList<T>& markerData,
            GeometricField<T,fvPatchField,volMesh>& fieldData,
            const Parameter& par
        );
};
#include "SensitivityInteraction.T.H"
}
#endif
