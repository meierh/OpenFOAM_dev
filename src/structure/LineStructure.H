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
        
        static scalar deltaDirac(vector X, vector x, scalar h);
        static scalar deltaDirac(vector X, vector x, vector h);
        static scalar correctedDeltaDirac(vector X, vector x, scalar h, const std::array<scalar,10>& b);
        static scalar correctedDeltaDirac(vector X, vector x, vector h, const std::array<scalar,10>& b);
        static scalar phiFunction(scalar r);
        
        void compNonUniformCorrWeights() {computeCorrectionWeights();}
        
    private:
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
        void minMaxSupportWidth();
        void dilationFactors();
        
        std::array<scalar,10> b;
        
        scalar computeMoment(vector indices) const;
        std::unique_ptr<gismo::gsMatrix<scalar>> computeMomentMatrix() const ;
        std::unique_ptr<std::array<scalar,10>> rescalingDiagonal() const;
        void computeCorrectionWeights();
        
        scalar deltaDirac(vector X, vector x) const;
        scalar correctedDeltaDirac(vector X, vector x) const;
        
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
            markerMeshType modusFieldToMarker = markerMeshType::Uniform,
            markerMeshType modusMarkerToField = markerMeshType::Uniform
        );
        ActiveRodMesh::rodMesh* rodMesh;
        std::vector<std::vector<LagrangianMarker*>> markers;
        dynamicRefineFvMesh& mesh;
        
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
        void fieldToMarker
        (
            const GeometricField<T,fvPatchField,volMesh>& fieldData,
            DynamicList<T>& markerData
        )
        {
            /*
            Info<<"Transfer field to markers"<<Foam::endl;
            const cellList& cells = mesh.cells();
            const faceList& faces = mesh.faces();
            const pointField& points = mesh.points();
            
            markerData.resize(markers.size());
            scalar factorSumAvg = 0;
            label count=0;
            for(label index=0; index<markerData.size(); index++)
            {
                T value = Foam::zero();
                const LagrangianMarker& oneMarker = *(markers[index]);
                const vector& oneMarkerPos = oneMarker.getMarkerPosition();
                scalar factorSum = 0;
                DynamicList<label> supportCells;
                for(auto cellIter=oneMarker.getSupportCells().begin();
                        cellIter!=oneMarker.getSupportCells().end();
                        cellIter++)
                {
                    supportCells.append(std::get<2>(*cellIter));
                }
                List<scalar> cellVols(supportCells.size());
                List<vector> cellCentres(supportCells.size());
                List<scalar> cellDistances(supportCells.size());
                for(label i=0; i<supportCells.size(); i++)
                {
                    label cellInd = supportCells[i];
                    cellVols[i] = cells[cellInd].mag(points,faces);
                    cellCentres[i] = cells[cellInd].centre(points,faces);
                    vector markerToCell = cellCentres[i]-oneMarkerPos;
                    cellDistances[i] = std::sqrt(markerToCell&markerToCell);
                }
                scalar h = std::accumulate(cellDistances.begin(),cellDistances.end(),0.0,std::plus<scalar>());
                if(supportCells.size()>0)
                    h/=supportCells.size();
                
                if(modusFieldToMarker==markerMeshType::Uniform)
                {
                    for(label i=0; i<supportCells.size(); i++)
                    {
                        label cellInd = supportCells[i];
                        scalar factor = deltaDirac(cellCentres[i],oneMarkerPos,h);
                        value += fieldData[cellInd]*factor*cellVols[i];
                        //Info<<cellInd<<":"<<fieldData[cellInd]<<" factor:"<<factor<<Foam::endl;
                        factorSum+=factor;
                    }
                }
                else
                {
                    DynamicList<scalar> epsilon;
                    SymmetricSquareMatrix<scalar> A(oneMarker.getSupportCells().size());
                    for(label i=0; i<supportCells.size(); i++)
                    {
                        
                        for(label j=i+1; j<supportCells.size(); j++)
                        {
                            label cellInd = supportCells[i];
                            scalar factor = deltaDirac(cellCentres[i],oneMarkerPos,h)*cellVols[i];
                            value += fieldData[cellInd]*factor;
                            factorSum+=factor;
                        }
                    }                    
                    for(auto cellIter=oneMarker.getSupportCells().begin();
                        cellIter!=oneMarker.getSupportCells().end();
                        cellIter++)
                    {
                        label cellInd = std::get<2>(*cellIter);
                        const cell& oneCell = cells[cellInd];
                        scalar cellVol = oneCell.mag(points,faces);
                        point cellCentre = oneCell.centre(points,faces);
                        //scalar h = std::cbrt(cellVol);
                        scalar factor = deltaDirac(cellCentre,oneMarkerPos,h)*cellVol;
                        value += fieldData[cellInd]*factor;
                        //Info<<factor<<"  ";
                        factorSum+=factor;
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
            */
        }
        
        template<typename T>
        void markerToField
        (
            const DynamicList<T>& markerData,
            GeometricField<T,fvPatchField,volMesh>& fieldData
        )
        {
            /*
            Info<<"Transfer markers to field"<<Foam::endl;
            fieldData = Foam::zero();
            const cellList& cells = mesh.cells();
            const faceList& faces = mesh.faces();
            const pointField& points = mesh.points();
            
            if(markerData.size()!=markers.size())
                FatalErrorInFunction<<"Marker Data size mismatch"<< exit(FatalError);
            
            List<DynamicList<label>> cellMarkers(cells.size());
            for(label markerInd=0; markerInd<markers.size(); markerInd++)
            {
                const LagrangianMarker& oneMarker = *(markers[markerInd]);
                for(auto cellIter=oneMarker.getSupportCells().begin();
                    cellIter!=oneMarker.getSupportCells().end();
                    cellIter++)
                {
                    label cellInd = std::get<2>(*cellIter);
                    cellMarkers[cellInd].append(markerInd);
                }
            }
            
            scalar factorSumAvg = 0;
            scalar count = 0;
            for(label cellInd=0; cellInd<fieldData.size(); cellInd++)
            {
                const cell& oneCell = cells[cellInd];
                point cellCentre = oneCell.centre(points,faces);
                //scalar cellVol = oneCell.mag(points,faces);
                scalar factorSum = 0;
                if(cellMarkers[cellInd].size()>0)
                {
                    T value = Foam::zero();
                    //Info<<cellInd<<" -- ";
                    for(label markerInd : cellMarkers[cellInd])
                    {
                        const LagrangianMarker& oneMarker = *(markers[markerInd]);
                        const vector& oneMarkerPos = oneMarker.getMarkerPosition();
                        scalar markerVolume = oneMarker.getMarkerVolume();
                        scalar h = std::cbrt(markerVolume);
                        scalar factor = deltaDirac(cellCentre,oneMarkerPos,h)*markerVolume;
                        //Info<<"dd:"<<deltaDirac(cellCentre,oneMarkerPos,h)<<" mV:"<<markerVolume<<" f:"<<factor<<" m:"<<markerData[markerInd]<<"   ";
                        factorSum+=factor;
                        value += markerData[markerInd]*factor;
                    }
                    fieldData[cellInd] += value;
                    //Info<<factorSum<<Foam::endl;
                    count++;
                    factorSumAvg+=factorSum;
                }
            }
            factorSumAvg/=count;
            Info<<"---markerToField factorSumAvg:"<<factorSumAvg<<"--------"<<Foam::endl;
            */
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
        void reduceMarkers
        (
            const std::vector<std::pair<std::list<LagrangianMarker>::iterator,std::list<LagrangianMarker>*>>&
        );
        
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
        const label minRefinement = 1;
        markerMeshType modusFieldToMarker;
        markerMeshType modusMarkerToField;
        
        /*****************************************Collect halo markers******************************************/
        //void extractStructureHaloMarkers(const std::vector<LagrangianMarker*> markers);
        // cell --> [markers]
        List<DynamicList<LagrangianMarker*>> haloMarkers;
        
        
    private:
        const List<scalar> crossSecArea;
        // [rodNumber] -> ...
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
        
        std::unique_ptr<std::pair<
            std::unordered_multimap<label,std::unordered_multimap<label,DynamicList<std::vector<vector>>>>,
            std::unordered_multimap<label,std::unordered_multimap<label,DynamicList<std::vector<scalar>>>>
        >>shareMarkerData
        (
            const std::vector<LagrangianMarker*>& markers,
            const std::vector<std::vector<vector>>& markersVecValues={},
            const std::vector<std::vector<scalar>>& markersScalValues={}  
        );
        scalar initialSpacingFromMesh();
        
        // neighbourProcess -> cell -> markersData
        std::unordered_multimap<label,std::unordered_multimap<label,DynamicList<std::pair<vector,scalar>>>> globalMarkerData;
        
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
