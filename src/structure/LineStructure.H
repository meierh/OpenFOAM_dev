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
        std::string to_string();
        
        inline label getRodNumber() const {return rodNumber;}
        inline const ActiveRodMesh::rodCosserat* getBaseRod() const {return baseRod;}
        inline scalar getMarkerParameter() const {return markerParameter;}
        inline vector getMarkerPosition() const {return markerPosition;}
        inline label getMarkerCell() const {return markerCell;}
        inline void setMarkerVolume(scalar volume) {markerVolume=volume;}
        inline scalar getMarkerVolume() const {return markerVolume;}
        inline const DynamicList<label>& getSupportCells() const {return supportCells;}
        inline vector getDilation() const {return dilation;}

    private:
        const dynamicRefineFvMesh& mesh;
        label rodNumber;
        const ActiveRodMesh::rodCosserat* baseRod;
        scalar markerParameter;
        vector markerPosition;
        label markerCell;
        scalar markerVolume;
        
        DynamicList<label> supportCells;
        vector h_plus;
        vector h_minus;
        vector dilation;
        
        void computeSupport(label iterations=1);
        void minMaxSupportWidth();
        void dilationFactors();
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
        std::vector<LagrangianMarker*> markers;
        dynamicRefineFvMesh& mesh;
        
        static scalar deltaDirac(vector X, vector x, scalar h);
        static scalar deltaDirac(vector X, vector x, vector h);
        
        void computeMarkerEpsilon
        (
            const dynamicRefineFvMesh& mesh,
            const DynamicList<LagrangianMarker*>& markers
        );

    protected:
        static scalar phiFunction(scalar r);
        static scalar correctedDeltaDirac(vector X, vector x, scalar h, std::array<scalar,10> b);
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
        
        static std::unique_ptr<gismo::gsMatrix<scalar>> computeMarkerEpsilonMatrix
        (
            const dynamicRefineFvMesh& mesh,
            const DynamicList<LagrangianMarker*>& markers
        );

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
                    supportCells.append(*cellIter);
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
                    /*
                    if(supportCells.size()>0)
                    {
                        Info<<"cellVols:"<<cellVols<<Foam::endl;
                        Info<<"cellCentres:"<<cellCentres<<Foam::endl;
                        Info<<"cellDistances:"<<cellDistances<<Foam::endl;
                        Info<<"h:"<<h<<Foam::endl;
                        Info<<"oneMarkerPos:"<<oneMarkerPos<<Foam::endl;
                        Info<<"value:"<<value<<Foam::endl;
                        FatalErrorInFunction<<"Temp stop"<<exit(FatalError);
                    }
                    */
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
                        label cellInd = *cellIter;
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
            
            List<DynamicList<label>> cellMarkers(cells.size());
            for(label markerInd=0; markerInd<markers.size(); markerInd++)
            {
                const LagrangianMarker& oneMarker = *(markers[markerInd]);
                for(auto cellIter=oneMarker.getSupportCells().begin();
                    cellIter!=oneMarker.getSupportCells().end();
                    cellIter++)
                {
                    label cellInd = *cellIter;
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
        void transferMarkers(FieldMarkerStructureInteraction& connector);
        
        static vector evaluateRodPos
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parameter
        );
        
        void test();
        
    protected:        
        static std::unique_ptr<std::vector<scalar>> createSpacedPointsOnRod
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar spacing
        );
        
        std::pair<scalar,scalar> minMaxSpan (const cell& thisCell);
        
        void cellDistances
        (
            const cell& thisCell,
            std::vector<scalar>&
        );
        
        //scalar supportDomainMinSize(label cellInd);
        scalar supportDomainMinSize(const DynamicList<label>& supportDomainCells);
        
        static scalar evaluateRodArcLen
        (
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar parA,
            scalar parB
        );
             
        /*
        bool doSubdivision
        (
            const LagrangianMarker& smallerSide,
            const LagrangianMarker& largerSide
        );
        */
        
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
        
        const label minRefinement = 1;
        markerMeshType modusFieldToMarker;
        markerMeshType modusMarkerToField;        
    private:
        const List<scalar> crossSecArea;
        
        std::vector<std::unique_ptr<std::vector<LagrangianMarker>>> rodMarkers;
        std::vector<std::unique_ptr<std::vector<scalar>>> initialRodPoints;
        
        scalar initialSpacing;
        
        std::unique_ptr<std::vector<LagrangianMarker>> constructMarkerSet
        (
            label rodNumber,
            const ActiveRodMesh::rodCosserat* oneRod,
            scalar crossSecArea,
            scalar initialSpacing,
            std::pair<bool,scalar> refineSpacing = {false,0},
            bool reInitialize = false
        );
        
        /*
        std::unique_ptr<gismo::gsMatrix<scalar>> computeMarkerWeightMatrix
        (
            std::vector<LagrangianMarker*> markers
        );
        */
        
        scalar initialSpacingFromMesh();
};
}
#endif
