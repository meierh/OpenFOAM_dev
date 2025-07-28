#include "MarkerFilter.H"
#include <type_traits>

template<typename T>
Foam::MarkerFilter<T>::MarkerFilter
(
    const fvMesh& mesh,
    LineStructure& structure,
    VolField<T>& output_Force,
    const IOdictionary& structureDict,
    std::string filterName
):
output_Force_last(output_Force),
output_Force(output_Force),
mesh(mesh),
structure(structure),
structureDict(structureDict)
{
    if(structureDict.found("markerForceFilter"))
    {
        const dictionary& markerForceFilterDict = structureDict.subDict("markerForceFilter");
        List<keyType> markerForceFilterKeys = markerForceFilterDict.keys();
        for(keyType oneFilterKey : markerForceFilterKeys)
        {
            if(oneFilterKey.match(filterName,true))
            {
                if(filterActive)
                    FatalErrorInFunction<<"Double activation of filter for:"<<filterName<<exit(FatalError);
                filterActive = true;
                
                const dictionary& oneFilterDict = markerForceFilterDict.subDict(filterName);
                List<keyType> oneFilterKeys = oneFilterDict.keys();
                for(keyType oneFilterTypeKey : oneFilterKeys)
                {
                    if(oneFilterTypeKey.match("spatial",false))
                    {
                        if(spatialFilterActive)
                            FatalErrorInFunction<<"Double activation of filter type: spatial"<<exit(FatalError);
                        spatialFilterActive = true;
                        
                        const dictionary& spatialFilterDict = oneFilterDict.subDict("spatial");
                        
                        label width = readFromDictionary<label>("width",spatialFilterDict);
                        switch (width)
                        {
                            case 0:
                                this->width = spatialFilterWidth::W1x1;
                                break;
                            case 1:
                                this->width = spatialFilterWidth::W3x3;
                                break;
                            case 2:
                                this->width = spatialFilterWidth::W5x5;
                                break;
                            default:
                                FatalErrorInFunction<<"Invalid constant/structureDict/"<<filterName<<"/spatial/width; expected {0,1,2} but got "<<width<<exit(FatalError);
                        }

                        word method = readFromDictionary<word>("method",spatialFilterDict);
                        if (method=="Clipping")
                        {
                            this->method = spatialFilterMethod::Clipping;
                            this->clipLimitFactor = readFromDictionary<scalar>("clipFactor",spatialFilterDict);
                            Info<<"clipLimitFactor:"<<clipLimitFactor<<Foam::endl;
                        }
                        else if(method=="SmoothLimiting")
                        {
                            this->method = spatialFilterMethod::SmoothLimiting;
                            this->smoothLimitParameter = readFromDictionary<scalar>("smoothLimitParameter",spatialFilterDict);
                            this->smoothingArea = readFromDictionary<scalar>("smoothingArea",spatialFilterDict);
                            Info<<"smoothLimitParameter:"<<smoothLimitParameter<<Foam::endl;
                            Info<<"smoothingArea:"<<smoothingArea<<Foam::endl;
                        }
                        else if(method=="LogSmoothLimiting")
                        {
                            this->method = spatialFilterMethod::LogSmoothLimiting;
                            this->thresholdParameter = readFromDictionary<scalar>("thresholdFactor",spatialFilterDict);
                            this->dampingHardness = readFromDictionary<scalar>("dampingHardness",spatialFilterDict);
                            Info<<"thresholdParameter:"<<thresholdParameter<<Foam::endl;
                            Info<<"dampingHardness:"<<dampingHardness<<Foam::endl;
                        }
                        else
                            FatalErrorInFunction<<"Invalid constant/structureDict/"<<filterName<<"/spatial/method; expected {Clipping,SmoothLimiting,LogSmoothLimiting}"<<exit(FatalError);
                        
                    }
                    if(oneFilterTypeKey.match("temporal",false))
                    {
                        if(temporalFilterActive)
                            FatalErrorInFunction<<"Double activation of filter type: temporal"<<exit(FatalError);
                        temporalFilterActive = true;
                        
                        const dictionary& temporalFilterDict = oneFilterDict.subDict("temporal");
                        if(temporalFilterDict.found("theta"))
                        {
                            parameter_EMA = readFromDictionary<scalar>("theta",temporalFilterDict);
                        }
                        else if(temporalFilterDict.found("cutOffFrequency"))
                        {
                            scalar cutOffFrequency = readFromDictionary<scalar>("cutOffFrequency",temporalFilterDict);
                            parameter_EMA = computeEMAParameterFrom3dBFrequency(cutOffFrequency);
                        }
                        else 
                            FatalErrorInFunction<<"Missing constant/structureDict/"<<filterName<<"/temporal/...; expected {theta,cutOffFrequency}"<<exit(FatalError);
                        
                        Info<<"cutOffFrequency:"<<compute3dBFrequencyFromEMAParameter(parameter_EMA)<<Foam::endl;
                        Info<<"parameter_EMA:"<<parameter_EMA<<Foam::endl;
                    }
                }
            }
        }
    }
        
    output_Force = Foam::zero();
    output_Force_last = Foam::zero();
}

template<typename T>
void Foam::MarkerFilter<T>::record()
{
    if(filterActive)
    {
        if(temporalFilterActive)
        {
            scalar currentTime = mesh.time().value();
            if(currentTime != output_Force_last_update_time)
            {
                output_Force_last = output_Force;
                output_Force_last_update_time = currentTime;
                Info<<"----------------------------------------------Recorded output_Force before "<<currentTime<<Foam::endl;
            }
        }
    }
}

template<typename T>
void Foam::MarkerFilter<T>::filter()
{
    spatialFilter();
    temporalFilter();
}

template<typename T>
Foam::scalar Foam::MarkerFilter<T>::compute3dBFrequencyFromEMAParameter(scalar parameter_EMA)
{
    return parameter_EMA/(2*Foam::constant::mathematical::pi*(1-parameter_EMA)*mesh.time().deltaTValue());
}

template<typename T>
Foam::scalar Foam::MarkerFilter<T>::computeEMAParameterFrom3dBFrequency(scalar frequency_3dB)
{
    scalar theta = 2 * Foam::constant::mathematical::pi * mesh.time().deltaTValue() * frequency_3dB;
    return theta/(1+theta);
}

template<typename T>
void Foam::MarkerFilter<T>::spatialFilter()
{
    if(filterActive)
    {
        if(spatialFilterActive)
        {
            std::unique_ptr<List<List<T>>> output_Force_halo_Ptr = structure.broadcastHaloFields(output_Force);
            const List<List<T>>& output_Force_halo = *output_Force_halo_Ptr;       
            const List<List<List<Pair<label>>>>& meshGraph = structure.getMeshGraph();
            for(label cellInd = 0; cellInd<output_Force.size(); cellInd++)
            {
                std::unordered_set<Pair<label>,foamPairHash<label>> cellNeighbourhood;
                cellNeighbourhood.insert({Pstream::myProcNo(),cellInd});
                DynamicList<Pair<label>> front;
                front.append({Pstream::myProcNo(),cellInd});
                for(label neigIter=0; neigIter<static_cast<int>(width); neigIter++)
                {
                    std::unordered_set<Pair<label>,foamPairHash<label>> tmp_front;
                    for(const Pair<label>& frontCell : front)
                    {
                        const List<Pair<label>>& cellNode =  meshGraph[frontCell.first()][frontCell.second()];
                        for(const Pair<label>& neiCell : cellNode)
                        {
                            if(cellNeighbourhood.find(neiCell)!=cellNeighbourhood.end())
                            {
                                tmp_front.insert(neiCell);
                            }
                        }
                    }
                    front.clear();
                    for(const Pair<label>& p : tmp_front)
                        front.append(p);
                    cellNeighbourhood.insert(tmp_front.begin(),tmp_front.end());
                    tmp_front.clear();
                }
                
                T cCForce = output_Force[cellInd];
                scalar cCForce_abs;;
                scalar max_Force_abs,min_Force_abs,avg_Force_abs;
                if constexpr (std::is_same<vector,T>::value)
                {
                    cCForce_abs = std::sqrt(cCForce&cCForce);
                    max_Force_abs = cCForce_abs;
                    min_Force_abs = cCForce_abs;
                }
                else
                {
                    cCForce_abs = std::abs(cCForce);
                    max_Force_abs = std::abs(cCForce);
                    min_Force_abs = std::abs(cCForce);
                }
                avg_Force_abs = Foam::zero();
                std::vector<scalar> median_Force_abs_list;
                for(const Pair<label>& neiCell : cellNeighbourhood)
                {
                    T cellVal = Foam::zero();
                    if(neiCell.first()==Pstream::myProcNo())
                    {
                        cellVal = output_Force[neiCell.second()];
                    }
                    else
                    {
                        const std::unordered_map<label,label>& haloCellToIndexMap = structure.getHaloCellToIndexMap(neiCell.first());
                        auto iter = haloCellToIndexMap.find(neiCell.second());
                        if(iter==haloCellToIndexMap.end())
                            FatalErrorInFunction<<"Invalid halo cell index"<<exit(FatalError);
                        label index = iter->second;
                        if(index>=0 || index<output_Force_halo[neiCell.first()].size())
                            FatalErrorInFunction<<"Out of scope error"<<exit(FatalError);
                        cellVal = output_Force_halo[neiCell.first()][index];
                    }
                    scalar cellValLen; 
                    if constexpr (std::is_same<vector,T>::value)
                    {
                        cellValLen = std::sqrt(cellVal&cellVal);
                    }
                    else
                    {
                        cellValLen = std::abs(cellVal);
                    }
                    max_Force_abs = std::max<scalar>(cellValLen,max_Force_abs);
                    min_Force_abs = std::min<scalar>(cellValLen,min_Force_abs);
                    avg_Force_abs += cellValLen;
                    median_Force_abs_list.push_back(cellValLen);
                }
                std::sort(median_Force_abs_list.begin(),median_Force_abs_list.end());
                label med = median_Force_abs_list.size();
                med = med/2;
                scalar median_Force_abs = median_Force_abs_list[med];
                avg_Force_abs /= cellNeighbourhood.size();
                
                scalar scalingFactor;
                switch (method)
                {
                    case Clipping:
                        scalingFactor = clip(median_Force_abs,cCForce_abs); 
                        break;
                    case SmoothLimiting:
                        scalingFactor = smoothLimit(median_Force_abs,cCForce_abs); 
                        break;
                    case LogSmoothLimiting:
                        scalingFactor = logSmoothLimit(median_Force_abs,cCForce_abs);
                        break;
                    default:
                        scalingFactor = 0;
                        FatalErrorInFunction<<"Invalid method"<<exit(FatalError);
                }
                output_Force[cellInd] = scalingFactor*cCForce;
            }
            Info<<"----------------------------------------------Spatially filtered forcing field"<<Foam::endl;
        }
    }
}

template<typename T>
Foam::scalar Foam::MarkerFilter<T>::clip(scalar referenceForce, scalar actualForce)
{
    scalar clippedForce = std::min<scalar>(actualForce,clipLimitFactor*referenceForce);
    if(actualForce==0)
        return 0;
    return clippedForce/actualForce;
}

template<typename T>
Foam::scalar Foam::MarkerFilter<T>::smoothLimit(scalar referenceForce, scalar actualForce)
{
    scalar clippingFactor = 1-smoothLimitParameter*std::tanh(((actualForce/referenceForce)-1)/smoothingArea);
    return clippingFactor;
}


template<typename T>
Foam::scalar Foam::MarkerFilter<T>::logSmoothLimit(scalar referenceForce, scalar actualForce)
{
    if(referenceForce==0)
        return 1;
    scalar s = actualForce/referenceForce;
    if(s>1)
    {       
        return std::log(1+dampingHardness)/std::log(1+dampingHardness*(s/thresholdParameter));
    }
    else
        return 1;
}

template<typename T>
void Foam::MarkerFilter<T>::temporalFilter()
{
    if(temporalFilterActive)
    {
        output_Force = parameter_EMA*output_Force + (1-parameter_EMA)*output_Force_last;
        Info<<"----------------------------------------------Spatially filtered forcing field"<<Foam::endl;
    }
}
