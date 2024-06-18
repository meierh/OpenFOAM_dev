#ifndef LINESTRUCTURE_H
#define LINESTRUCTURE_H
#include "Structure.H"
//#include "rodTools.h"

namespace Foam
{

class LagrangianMarker
{
    public:
        explicit LagrangianMarker
        (
            const Structure& structure,
            const dynamicRefineFvMesh& mesh,
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* baseRod,
            const scalar markerParameter
        );
    protected:
        explicit LagrangianMarker
        (    
            const Structure& structure,
            const dynamicRefineFvMesh& mesh,
            const label rodNumber,
            const ActiveRodMesh::rodCosserat* baseRod
        );
        
    public:
        vector getMarkerVelocity();
        scalar getMarkerTemperature();
        scalar getMarkerCellVolume();
        std::pair<scalar,scalar> getMarkerCellSpacing();
        scalar getMarkerCellMinSpacing();
        std::string to_string() const;
        void evaluateMarker();
        
        inline label getRodNumber() const {return rodNumber;}
        inline const ActiveRodMesh::rodCosserat* getBaseRod() const {return baseRod;}
        inline scalar getMarkerParameter() const {return markerParameter;}
        inline vector getMarkerPosition() const {return markerPosition;}
        inline label getMarkerCell() const {return markerCell;}
        inline void setMarkerVolume(scalar volume) {markerVolume=volume;}
        inline scalar getMarkerVolume() const {return markerVolume;}
        inline const DynamicList<std::tuple<bool,label,label>>& getSupportCells() const {return supportCells;}
        inline vector getDilation() const {return dilation;}
        inline const std::array<scalar,10>& getCorrParaB() const {return b;}
        
        static scalar deltaDirac(vector X, vector x, scalar h);
        static scalar deltaDirac(vector X, vector x, vector h);
        static scalar correctedDeltaDirac(vector X, vector x, scalar h, const std::array<scalar,10>& b);
        static scalar correctedDeltaDirac(vector X, vector x, vector h, const std::array<scalar,10>& b);
        static scalar phiFunction(scalar r);
        
        //h=getDilation() b=b;
        scalar correctedDeltaDirac(vector X, vector x) const;
        scalar deltaDirac(vector X, vector x) const;
        
        void compNonUniformCorrWeights() {computeCorrectionWeights();}
        
    protected:        
        const Structure& structure;
        const dynamicRefineFvMesh& mesh;
        label rodNumber;
        const ActiveRodMesh::rodCosserat* baseRod;
        scalar markerParameter;
        vector markerPosition;
        label markerCell;
        scalar markerVolume;
        
        // std::tuple<own processor,processor ID,cell ID>
        DynamicList<std::tuple<bool,label,label>> supportCells;
        Vector<bool> existingDims;
        vector h_plus;
        vector h_minus;
        vector dilation;
        
        void computeSupport(label iterations=1);
        void computeSupportDimensions();
        void minMaxSupportWidth();
        void dilationFactors();
        
        std::array<scalar,10> b;
        
        scalar computeMoment(vector indices) const;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeMomentMatrix3D() const ;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeMomentMatrix2D(label normalDim) const;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeMomentMatrix1D(label dim) const;
        std::unique_ptr<std::array<scalar,10>> rescalingDiagonal3D() const;
        std::unique_ptr<std::array<scalar,6>> rescalingDiagonal2D(label normalDim) const;
        std::unique_ptr<std::array<scalar,3>> rescalingDiagonal1D(label dim) const;
        void computeCorrectionWeights();
                
        scalar computeCorrectedMoment(vector indices) const;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeCorrectedMomentMatrix() const;
};

class FieldMarkerStructureInteraction
{
    public:
        enum markerMeshType {Uniform, NonUniform};
        explicit FieldMarkerStructureInteraction
        (
            dynamicRefineFvMesh& mesh,
            Structure& structure,
            markerMeshType modusFieldToMarker = markerMeshType::Uniform,
            markerMeshType modusMarkerToField = markerMeshType::Uniform
        );
        ActiveRodMesh::rodMesh* rodMesh;
        std::vector<LagrangianMarker*> markers;
        
        //std::vector<std::vector<LagrangianMarker*>> markers;
        dynamicRefineFvMesh& mesh;
        Structure& structure;
        
        /*
        static scalar deltaDirac(vector X, vector x, scalar h);
        static scalar deltaDirac(vector X, vector x, vector h);
        */
        
        /*
        void computeMarkerEpsilon
        (
            const dynamicRefineFvMesh& mesh,
            const DynamicList<LagrangianMarker*>& markers
        );
        */

    protected:
        /*
        static scalar phiFunction(scalar r);
        static scalar correctedDeltaDirac(vector X, vector x, scalar h, std::array<scalar,10> b);
        */
        
        /*
        static scalar computeMoment
        (
            const dynamicRefineFvMesh& mesh,
            const LagrangianMarker& marker,
            vector indices
        );
        static std::unique_ptr<gismo::gsMatrix<scalar>> computeMomentMatrix
        (
            const dynamicRefineFvMesh& mesh,
            const LagrangianMarker& marker
        );
        static std::unique_ptr<std::array<scalar,10>> rescalingDiagonal(const LagrangianMarker& marker);
        std::unique_ptr<std::array<scalar,10>> computeCorrectionWeights
        (
            const dynamicRefineFvMesh& mesh,
            const LagrangianMarker& marker
        );
        */
        
        /*
        static std::unique_ptr<gismo::gsMatrix<scalar>> computeMarkerEpsilonMatrix
        (
            const dynamicRefineFvMesh& mesh,
            const DynamicList<LagrangianMarker*>& markers
        );
        */
        
        template<typename T>
        std::unique_ptr<List<List<T>>> fieldHaloExchange
        (
            const GeometricField<T,fvPatchField,volMesh>& fieldData
        )
        {
            auto globalHaloDataPtr = std::unique_ptr<List<List<T>>>(new List<List<T>>(Pstream::nProcs()));
            List<List<T>>& globalHaloData = *globalHaloDataPtr;
            List<T>& ownHaloData = globalHaloData[Pstream::myProcNo()];
            const DynamicList<Structure::CellDescription>& ownHaloCells = structure.getHaloCellList(Pstream::myProcNo());
            ownHaloData.resize(ownHaloCells);
            for(label haloCellInd=0; haloCellInd<ownHaloCells.size(); haloCellInd++)
            {
                const Structure::CellDescription& oneHaloCell = ownHaloCells[haloCellInd];
                ownHaloData[haloCellInd] = fieldData[oneHaloCell.index];
            }
            Pstream::gatherList(globalHaloData);
            Pstream::scatterList(globalHaloData);
            return globalHaloDataPtr;
        }

        template<typename T>
        void fieldToMarker
        (
            const GeometricField<T,fvPatchField,volMesh>& fieldData,
            DynamicList<T>& markerData
        )
        {
            Info<<"Transfer field to markers"<<Foam::endl;
            const cellList& cells = mesh.cells();
            const faceList& faces = mesh.faces();
            const pointField& points = mesh.points();
            
            std::unique_ptr<List<List<T>>> haloFieldsPtr = fieldHaloExchange(fieldData);
            List<List<T>>& haloFields = *haloFieldsPtr;

            markerData.resize(markers.size());
            scalar factorSumAvg = 0;
            label count=0;
            for(label index=0; index<markerData.size(); index++)
            {
                T value = Foam::zero();
                const LagrangianMarker& oneMarker = *(markers[index]);
                const vector& oneMarkerPos = oneMarker.getMarkerPosition();
                scalar factorSum = 0;
                const DynamicList<std::tuple<bool,label,label>>& supportCells = oneMarker.getSupportCells();
                List<scalar> cellVols(supportCells.size());
                List<vector> cellCentres(supportCells.size());
                List<scalar> cellDistances(supportCells.size());
                List<T> cellValues(supportCells.size());
                for(label suppInd=0; suppInd<supportCells.size(); suppInd++)
                {
                    const std::tuple<bool,label,label>& oneSuppCell = supportCells[suppInd];
                    label proc = std::get<1>(oneSuppCell);
                    label cellInd = std::get<2>(oneSuppCell);
                    if(std::get<0>(oneSuppCell))
                    {
                        cellVols[suppInd] = cells[cellInd].mag(points,faces);
                        cellCentres[suppInd] = cells[cellInd].centre(points,faces);
                        vector markerToCell = cellCentres[suppInd]-oneMarkerPos;
                        cellDistances[suppInd] = std::sqrt(markerToCell&markerToCell);
                        cellValues[suppInd] = fieldData[cellInd];
                    }
                    else
                    {
                        const std::unordered_map<label,label>& procCellToInd = structure.getHaloCellToIndexMap(proc);
                        auto iterIndex = procCellToInd.find(cellInd);
                        if(iterIndex!=procCellToInd.end())
                            FatalErrorInFunction<<"Halo cell does not exist"<<exit(FatalError);
                        label haloIndex = iterIndex->second;
                        const Structure::CellDescription& haloCell = structure.getHaloCellList(proc)[haloIndex];

                        cellVols[suppInd] = haloCell.volume;
                        cellCentres[suppInd] = haloCell.centre;
                        vector markerToCell = cellCentres[suppInd]-oneMarkerPos;
                        cellDistances[suppInd] = std::sqrt(markerToCell&markerToCell);
                        cellValues[suppInd] = haloFields[proc][haloIndex];
                    }
                }
                
                if(modusFieldToMarker==markerMeshType::Uniform)
                {
                    for(label suppCellInd=0; suppCellInd<supportCells.size(); suppCellInd++)
                    {
                        using LM=LagrangianMarker;
                        vector& X = oneMarkerPos;
                        vector& x = cellCentres[suppCellInd];
                        scalar factor = oneMarker.correctedDeltaDirac(X,x);
                        value += factor*cellValues[suppCellInd]*cellVols[suppCellInd];
                    }
                }
                else
                {
                    for(label suppCellInd=0; suppCellInd<supportCells.size(); suppCellInd++)
                    {
                        using LM=LagrangianMarker;
                        vector& X = oneMarkerPos;
                        vector& x = cellCentres[suppCellInd];
                        scalar factor = oneMarker.deltaDirac(X,x);
                        value += factor*cellValues[suppCellInd]*cellVols[suppCellInd];
                    }
                }
                //Info<<"value:"<<value<<Foam::endl;
                markerData[index] = value;
                if(oneMarker.getSupportCells().size()>0)
                {
                    factorSumAvg+=factorSum;
                    count++;
                }
            }
            factorSumAvg/=count;
            Info<<"--------fieldToMarker  factorSumAvg:"<<factorSumAvg<<"------"<<Foam::endl;
            //Info<<markerData<<Foam::endl;
        }
        
        template<typename T>
        std::unique_ptr<List<T>> markerFieldWriteHalo
        (
            const List<List<T>>& globalMarkerFieldData
        )
        {
            // writer proc -> cell proc -> haloCells
            List<List<List<T>>> allProcsGlobalMarkerFieldData(Pstream::nProcs());
            allProcsGlobalMarkerFieldData[Pstream::myProcNo()] = globalMarkerFieldData;
            
            Pstream::gatherList(allProcsGlobalMarkerFieldData);
            Pstream::scatterList(allProcsGlobalMarkerFieldData);
            
            auto ownProcGlobalMarkerFieldDataPtr = std::unique_ptr<List<T>>(new List<T>());
            List<T>& ownProcGlobalMarkerFieldData = *ownProcGlobalMarkerFieldDataPtr;
            ownProcGlobalMarkerFieldData.resize(structure.getHaloCellList(Pstream::myProcNo()).size(),Foam::zero());
            for(label writerProc=0; writerProc<allProcsGlobalMarkerFieldData.size(); writerProc++)
            {
                List<List<T>>& oneProcGlobalMarkerFieldData = allProcsGlobalMarkerFieldData[writerProc];
                List<T>& oneProcOwnProcMarkerFieldData = oneProcGlobalMarkerFieldData[Pstream::myProcNo()];
                for(label haloCellInd=0; haloCellInd<oneProcOwnProcMarkerFieldData.size(); haloCellInd++)
                {
                    ownProcGlobalMarkerFieldData[haloCellInd] += oneProcOwnProcMarkerFieldData[haloCellInd];
                }
            }
            return ownProcGlobalMarkerFieldDataPtr;
        }
        
        template<typename T>
        void markerToField
        (
            const DynamicList<T>& markerData,
            GeometricField<T,fvPatchField,volMesh>& fieldData
        )
        {
            Info<<"Transfer markers to field"<<Foam::endl;
            fieldData = Foam::zero();
            const cellList& cells = mesh.cells();
            const faceList& faces = mesh.faces();
            const pointField& points = mesh.points();
            
            if(markerData.size()!=markers.size())
                FatalErrorInFunction<<"Marker Data size mismatch"<< exit(FatalError);

            List<List<T>>& haloFields(Pstream::nProcs());
            for(label proc=0; proc<Pstream::nProcs(); proc++)
                haloFields[proc].resize(structure.getHaloCellList(proc).size(),Foam::zero());
            
            markerData.resize(markers.size());
            scalar factorSumAvg = 0;
            label count=0;
            for(label index=0; index<markerData.size(); index++)
            {
                T value = Foam::zero();
                const LagrangianMarker& oneMarker = *(markers[index]);
                const vector& oneMarkerPos = oneMarker.getMarkerPosition();
                scalar factorSum = 0;
                const DynamicList<std::tuple<bool,label,label>>& supportCells = oneMarker.getSupportCells();

                for(label suppInd=0; suppInd<supportCells.size(); suppInd++)
                {
                    const std::tuple<bool,label,label>& oneSuppCell = supportCells[suppInd];
                    label proc = std::get<1>(oneSuppCell);
                    label cellInd = std::get<2>(oneSuppCell);
                    if(std::get<0>(oneSuppCell))
                    {
                        vector cellCentre = cells[cellInd].centre(points,faces);
                        scalar factor = oneMarker.deltaDirac(cellCentre,oneMarkerPos);
                        factor *= oneMarker.getMarkerVolume();
                        fieldData[cellInd] += markerData[index]*factor;
                    }
                    else
                    {
                        const std::unordered_map<label,label>& procCellToInd = structure.getHaloCellToIndexMap(proc);
                        auto iterIndex = procCellToInd.find(cellInd);
                        if(iterIndex!=procCellToInd.end())
                            FatalErrorInFunction<<"Halo cell does not exist"<<exit(FatalError);
                        label haloIndex = iterIndex->second;
                        const Structure::CellDescription& haloCell = structure.getHaloCellList(proc)[haloIndex];
                        vector& cellCentre = haloCell.centre;
                        scalar factor = oneMarker.deltaDirac(cellCentre,oneMarkerPos);
                        factor *= oneMarker.getMarkerVolume();
                        haloFields[proc][haloIndex] += markerData[index]*factor;
                    }
                }
            }
            
            std::unique_ptr<List<T>> ownHaloCellValues = markerFieldWriteHalo(haloFields);
            const DynamicList<Structure::CellDescription> ownHaloCells = structure.getHaloCellList(Pstream::myProcNo());
            if(ownHaloCells.size()!=ownHaloCellValues.size())
                FatalErrorInFunction<<"Mismatch in own halo cell size!"<<exit(FatalError);
            
            for(label haloCellInd=0; haloCellInd<ownHaloCells.size(); haloCellInd++)
            {
                fieldData[ownHaloCells[haloCellInd].index] += ownHaloCellValues[haloCellInd];
            }
        }
        
        void scatterNurbs(std::pair<gsNurbs<scalar>,label> in, std::pair<gsNurbs<scalar>,label>& out);
        
        template<typename T>
        void gatherMarkerData
        (
            std::tuple<List<scalar>,List<List<T>>,label> in,
            std::tuple<List<scalar>,List<List<T>>,label>& out
        )
        {
            List<label> rodInIndex(Pstream::nProcs(),-1);
            rodInIndex[Pstream::myProcNo()] = std::get<2>(in);
            Pstream::gatherList(rodInIndex);
            if(Pstream::master())
            {
                for(label ind : rodInIndex)
                    if(ind!=std::get<2>(in))
                        FatalErrorInFunction<<"Inconsistent rod index"<< exit(FatalError);
            }
            
            if(std::get<0>(in).size()!=std::get<1>(in).size())
                FatalErrorInFunction<<"Inconsistent marker array length"<< exit(FatalError);
            
            List<List<scalar>> globalMarkerParam(Pstream::nProcs());
            globalMarkerParam[Pstream::myProcNo()] = std::get<0>(in);
            Pstream::gatherList(globalMarkerParam);
            
            List<List<List<T>>> globalMarkerData(Pstream::nProcs());
            globalMarkerData[Pstream::myProcNo()] = std::get<1>(in);
            Pstream::gatherList(globalMarkerData);
            
            if(Pstream::master())
            {
                std::map<scalar,List<T>> mergedData;
                label elemLen = -1;
                for(label proc=0; proc<globalMarkerParam.size(); proc++)
                {
                    for(label ind=0; ind<globalMarkerParam[proc].size(); ind++)
                    {
                        scalar param = globalMarkerParam[proc][ind];
                        const List<T>& values = globalMarkerData[proc][ind];
                        if(elemLen==-1)
                            elemLen = values.size();
                        else
                            if(elemLen!=values.size())
                                FatalErrorInFunction<<"Inconsistent marker data array length"<< exit(FatalError);
                        auto iter = mergedData.find(param);
                        if(iter==mergedData.end())
                        {
                            mergedData.insert({param,values});
                        }
                        else
                        {
                            for(label i=0; i<values.size(); i++)
                                iter->second[i] += values[i];
                        }
                    }
                }
                
                List<scalar> masterParams(mergedData.size());
                List<List<T>> masterData(mergedData.size());
                label index=0;
                for(auto iter=mergedData.begin(); iter!=mergedData.end(); iter++)
                {
                    masterParams[index] = iter->first;
                    masterData[index] = iter->second;
                    index++;
                }
                
                std::get<0>(out) = masterParams;
                std::get<1>(out) = masterData;
                std::get<2>(out) = std::get<2>(in);
            }
        }
        
        const scalar h;
        const markerMeshType modusFieldToMarker;
        const markerMeshType modusMarkerToField;
};

class LineStructure : public Structure
{
    public:
        enum markerMeshType {Uniform, NonUniform};
        explicit LineStructure
        (
            dynamicRefineFvMesh& mesh,
            const List<scalar> crossSecArea,            
            markerMeshType modusFieldToMarker = markerMeshType::Uniform,
            markerMeshType modusMarkerToField = markerMeshType::Uniform
        );
        LineStructure(const LineStructure&) = delete;
        void connect(FieldMarkerStructureInteraction& connector);
        
        void createMarkers();
        
        void reInitializeMarkers(bool keepMarkers=false,bool keepSeedPoints=false);
        
        static vector evaluateRodPos
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter
        );
        
        void test();
        
    protected:
        
        void createSpacingPoints();
        void createSpacedPointsOnRod
        (
            label rodNumber,
            scalar spacing
        );
        
        void createMarkersFromSpacedPoints();
        void createMarkersFromSpacedPointsOnRod
        (
            label rodNumber
        );
        
        void refineMarkers
        (
            std::pair<bool,scalar> forcedSpacing = {false,-1}
        );
        void refineMarkersOnRod
        (
            label rodNumber,
            std::pair<bool,scalar> forcedSpacing = {false,-1}
        );
        
        void setMarkerVolume();
        void setMarkerVolumeOnRod
        (
            label rodNumber
        );
        
        void evaluateMarkerMeshRelation();
        void evaluateMarkerMeshRelation(std::list<LagrangianMarker>& markerList);
        
        void reduceMarkers();
        template<typename T>
        void reduceMarkers
        (
            const std::vector<std::pair<std::list<T>::iterator,std::list<T>*>>& allMarkers
        )
        {
            using SingleMarkerRef = std::pair<std::list<T>::iterator,std::list<T>*>;

            std::unordered_map<label,DynamicList<SingleMarkerRef>> cellToMarker;
            for(SingleMarkerRef markerData : allMarkers)
            {
                const T& marker = *(markerData.first);
                label cellInd = marker.getMarkerCell();
                cellToMarker[cellInd].append(markerData);
            }
            
            const cellList& cells = mesh.cells();
            const faceList& faces = mesh.faces();
            const pointField& points = mesh.points();
            std::vector<SingleMarkerRef> markersToRemove;
            for(auto iter=cellToMarker.begin(); iter!=cellToMarker.end(); iter++)
            {
                label cellInd = iter->first;
                if(cellInd>=0 && cellInd<mesh.cells().size())
                {

                    const cell& oneCell = cells[cellInd];
                    Info<<"Cell:"<<cellInd<<Foam::endl;
                    vector cellCentre = oneCell.centre(points,faces);
                    std::pair<
                        std::pair<
                            std::pair<
                                DynamicList<SingleMarkerRef>,
                                DynamicList<SingleMarkerRef>
                            >,
                            std::pair<
                                DynamicList<SingleMarkerRef>,
                                DynamicList<SingleMarkerRef>
                            >
                        >,
                        std::pair<
                            std::pair<
                                DynamicList<SingleMarkerRef>,
                                DynamicList<SingleMarkerRef>
                            >,
                            std::pair<
                                DynamicList<SingleMarkerRef>,
                                DynamicList<SingleMarkerRef>
                            >
                        >
                    > xyz_split;
                    for(SingleMarkerRef markerPtr : iter->second)
                    {
                        vector markerPos = markerPtr.first->getMarkerPosition();
                        Info<<"\t"<<"Marker:"<<markerPos<<Foam::endl;
                        
                        std::pair<
                            std::pair<
                                DynamicList<SingleMarkerRef>,
                                DynamicList<SingleMarkerRef>
                            >,
                            std::pair<
                                DynamicList<SingleMarkerRef>,
                                DynamicList<SingleMarkerRef>
                            >
                        >* x_side;
                        if(markerPos[0]<cellCentre[0])
                            x_side = &(xyz_split.first);
                        else
                            x_side = &(xyz_split.second);
                        
                        std::pair<
                            DynamicList<SingleMarkerRef>,
                            DynamicList<SingleMarkerRef>
                        >* xy_side;
                        if(markerPos[1]<cellCentre[1])
                            xy_side = &(x_side->first);
                        else
                            xy_side = &(x_side->second);
                        
                        DynamicList<SingleMarkerRef>* xyz_side;
                        if(markerPos[2]<cellCentre[2])
                            xyz_side = &(xy_side->first);
                        else
                            xyz_side = &(xy_side->second);
                        
                        xyz_side->append(markerPtr);
                    }

                    std::vector<DynamicList<SingleMarkerRef>*> subCells = 
                    {
                        &(xyz_split.first.first.first),
                        &(xyz_split.first.first.second),
                        &(xyz_split.first.second.first),
                        &(xyz_split.first.second.second),
                        &(xyz_split.second.first.first),
                        &(xyz_split.second.first.second),
                        &(xyz_split.second.second.first),
                        &(xyz_split.second.second.first)
                    };
                    
                    for(DynamicList<SingleMarkerRef>* subCell : subCells)
                    {
                        if(subCell->size()>0)
                        {
                            vector averagePosition = vector(0,0,0);
                            for(SingleMarkerRef marker : *subCell)
                                averagePosition += marker.first->getMarkerPosition();
                            averagePosition /= subCell->size();
                            
                            SingleMarkerRef optMarker;
                            scalar optMarkerDistToAvgPos = std::numeric_limits<scalar>::max();
                            for(SingleMarkerRef marker : *subCell)
                            {
                                vector distVec = marker.first->getMarkerPosition() - averagePosition;
                                scalar dist = std::sqrt(distVec&distVec);
                                if(dist<optMarkerDistToAvgPos)
                                {
                                    optMarker = marker;
                                    optMarkerDistToAvgPos = dist;
                                }
                            }
                            for(SingleMarkerRef marker : *subCell)
                                if(marker!=optMarker)
                                    markersToRemove.push_back(optMarker);
                        }
                    }
                }
                else
                {
                    for(SingleMarkerRef markerPtr : iter->second)
                    {
                        markersToRemove.push_back(markerPtr);
                    }
                }
            }
            
            for(SingleMarkerRef markerToRemove : markersToRemove)
            {
                markerToRemove.second->erase(markerToRemove.first);
            }
        }
        
        void collectHaloMarkers();
        
        void exchangeHaloMarkersData();
                
        template<typename T>
        std::unique_ptr<List<List<DynamicList<T>>>> exchangeHaloMarkersValues
        (
            const List<DynamicList<T>>& values
        )
        {
            auto result = std::unique_ptr<List<List<DynamicList<T>>>>(new List<List<DynamicList<T>>>(Pstream::nProcs()));
            List<List<DynamicList<T>>>& globalValues = *result;
            globalValues[Pstream::myProcNo()] = values;
            
            Pstream::gatherList(globalValues);
            Pstream::scatterList(globalValues);
            
            return result;
        }
        
        void computeMarkerCellWeights();
        void computeMarkerCellWeights(const std::vector<LagrangianMarker*>& markers);

        using LinearSystem = std::tuple
        <
            std::vector<LagrangianMarker*>,
            gismo::gsMatrix<scalar>,
            // markerI -> DynList<tuple<procK,haloCellIndK,markerIndK,aIK>>
            std::unordered_map<label,DynamicList<std::tuple<label,label,label,scalar>>>
        >;
        std::unique_ptr<LinearSystem> computeMarkerEpsilonMatrix();
        std::unique_ptr<LinearSystem> computeMarkerEpsilonMatrix(const std::vector<LagrangianMarker*>&);
        void computeMarkerWeights();
        
        
        
        
        
        /*
        std::unique_ptr<std::vector<LagrangianMarker>> constructMarkerSet
        (
            label rodNumber,
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar crossSecArea,
            scalar initialSpacing,
            bool reInitialize,
            std::pair<bool,scalar> refineSpacing = {false,0},
        );
        */
        
        static scalar evaluateRodArcLen
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parA,
            scalar parB
        );
                    
        template<typename T>
        static T integrateRodwise
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            std::function<T(scalar)> function
        )
        {
            return integrateRodwise<T>(oneRod,oneRod->m_Curve.knots().first(),oneRod->m_Curve.knots().last(),function);
        }
        
        template<typename T>
        static T integrateRodwise
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parA,
            scalar parB,
            std::function<T(scalar)> function
        )
        {
            if(parA>parB)
                FatalErrorInFunction<<"Integration from higher to lower value"<<exit(FatalError);
            
            std::set<scalar> knotSet;
            for(scalar knot : oneRod->m_Curve.knots())
                knotSet.insert(knot);
            
            label degree = oneRod->m_Curve.knots().degree();
            auto start = knotSet.insert(parA);
            auto end = knotSet.insert(parB);
            
            label numberOfAbcissa = (degree+1)*(degree+1);
            
            T totalValue = Foam::zero();
            
            std::set<scalar>::iterator currNode = start.first;
            std::set<scalar>::iterator nextNode = ++(start.first);
            for(; currNode!=knotSet.end() && nextNode!=knotSet.end() && currNode!=end.first ; )
            {
                scalar startPar = *currNode;
                scalar endPar = *nextNode;
                scalar dist = endPar-startPar;
                scalar step = dist/numberOfAbcissa;
                scalar initialStep = startPar+step/2;
                T summedValue = Foam::zero();
                for(label i=0;i<numberOfAbcissa;i++)
                {
                    summedValue += function(initialStep)*step;
                    initialStep+=step;
                }
                totalValue+=summedValue;
                
                currNode=nextNode;
                nextNode++;
            }
            return totalValue;
        }
        
        static scalar distance
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parA,
            scalar parB
        );
        
        static scalar distance(const LagrangianMarker& A, const LagrangianMarker& B);
        
        std::vector<FieldMarkerStructureInteraction*> connectedInteractions;
        //const label minRefinement = 1;
        markerMeshType modusFieldToMarker;
        markerMeshType modusMarkerToField;
        
        /*****************************************Collect halo markers******************************************/
        //void extractStructureHaloMarkers(const std::vector<LagrangianMarker*> markers);
        // cell --> [markers]
        //List<DynamicList<LagrangianMarker*>> haloMarkers;
        
        
    private:
        const List<scalar> crossSecArea;
        // [rodNumber] -> ...
    protected:
        std::vector<std::unique_ptr<std::vector<scalar>>> initialRodPoints;

        std::vector<std::unique_ptr<std::list<LagrangianMarker>>> rodMarkersList;
        
        /******************Halo markers data************************/
        // [haloCells] -> [marker of cell] -> ptr        
        std::vector<std::vector<const LagrangianMarker*>> haloCellsRodMarkersList;
        // [nProc] -> [haloCells] -> [markers] -> position        
        List<List<DynamicList<vector>>> globalHaloCellsMarkerPos;
        // [nProc] -> [haloCells] -> [markers] -> volume        
        List<List<DynamicList<scalar>>> globalHaloCellsMarkerVolume;
        // [nProc] -> [haloCells] -> [markers] -> volume        
        List<List<DynamicList<vector>>> globalHaloCellsMarkerDilation;
        // [nProc] -> [haloCells] -> [markers] -> volume
        List<List<DynamicList<DynamicList<vector>>>> globalHaloCellsMarkerSupportCellCentres;
        // [nProc] -> [haloCells] -> [markers] -> volume
        List<List<DynamicList<DynamicList<scalar>>>> globalHaloCellsMarkerSupportCellVolume;
        
        /*
        std::unique_ptr<std::pair<
            std::unordered_multimap<label,std::unordered_multimap<label,DynamicList<std::vector<vector>>>>,
            std::unordered_multimap<label,std::unordered_multimap<label,DynamicList<std::vector<scalar>>>>
        >>shareMarkerData
        (
            const std::vector<LagrangianMarker*>& markers,
            const std::vector<std::vector<vector>>& markersVecValues={},
            const std::vector<std::vector<scalar>>& markersScalValues={}  
        );
        */
        
        scalar initialSpacingFromMesh();
        
        // neighbourProcess -> cell -> markersData
        //std::unordered_multimap<label,std::unordered_multimap<label,DynamicList<std::pair<vector,scalar>>>> globalMarkerData;
        
        /*
        template<typename T>
        std::unique_ptr<List<List<T>>> exchangeHaloFields
        (
            const GeometricField<T,fvPatchField,volMesh>& fieldData
        )
        {
            const DynamicList<label>& selfHalo = getSelfHaloList_Sorted();
            auto results = std::unique_ptr<List<List<T>>>(new List<List<T>>(Pstream::nProcs()));
            List<List<T>>& haloData = *results;
            haloData[Pstream::myProcNo()].resize(selfHalo.size());
            for(label haloCellInd=0; haloCellInd<selfHalo.size(); haloCellInd++)
            {
                label cell = selfHalo[haloCellInd];
                T cellValue = fieldData[cell];
                haloData[haloCellInd] = cellValue;
            }
            Pstream::gatherList(haloData);
            Pstream::scatterList(haloData);
            return results;
        }
        */
        
        void check();

};
}
#endif
