#ifndef MARKERFILTER_H
#define MARKERFILTER_H
#include "LineStructure.H"
#include <unordered_set>

namespace Foam
{

template<typename T>
class MarkerFilter
{
    public:
        explicit MarkerFilter
        (
            const fvMesh& mesh,
            LineStructure& structure,
            VolField<T>& output_Uf,
            const IOdictionary& structureDict,
            std::string filterName
        );
        virtual ~MarkerFilter()=default;

        const LineStructure& getStructure() const {return structure;}
        const fvMesh& getMesh() const {return mesh;}
        
        virtual void record();
        virtual void filter();
        
        scalar compute3dBFrequencyFromEMAParameter(scalar parameter_EMA);
        scalar computeEMAParameterFrom3dBFrequency(scalar frequency_3dB);
        
    protected:
        bool filterActive = false;
        
        void spatialFilter();
        bool spatialFilterActive = false;
        enum spatialFilterWidth {W1x1=0, W3x3=1, W5x5=2};
        spatialFilterWidth width = W3x3;
        enum spatialFilterMethod {Clipping, SmoothLimiting, LogSmoothLimiting, PotentialLogSmoothLimit};
        spatialFilterMethod method = LogSmoothLimiting;
        std::uint64_t meshGraphSyncVersion = 0;
        List<DynamicList<Pair<label>>> cellNeighbourhood;
        const DynamicList<Pair<label>>& getOneCellNeighbourhood(label cellInd){return cellNeighbourhood[cellInd];}
        void computeCellNeighbourhood();
        // Clipping
        scalar clip(scalar referenceForce, scalar actualForce);
        scalar clipLimitFactor = 3; // [2-5]
        // SmoothLimiting
        scalar smoothLimit(scalar referenceForce, scalar actualForce);
        scalar smoothLimitParameter = 0.5; // [0-1] "beta"
        scalar smoothingArea = 0.75; // [0.5-1.0] "gamma"
        // LogSmoothLimiting
        scalar logSmoothLimit(scalar referenceForce, scalar actualForce);
        scalar thresholdParameter = 4;
        scalar dampingHardness = 6;
        // PotentialLogSmoothLimiting
        scalar potentialLogSmoothLimit(scalar referenceForce, scalar actualForce);
        label power = 4;
        
        void temporalFilter();
        bool temporalFilterActive = false;
        scalar parameter_EMA = 0.5; // 
        VolField<T> output_Force_last;
        scalar output_Force_last_update_time = -1;   
        
        VolField<T>& output_Force;     
        const fvMesh& mesh;
        LineStructure& structure;
        const IOdictionary& structureDict;
};
#include "MarkerFilter.T.H"

}
#endif
