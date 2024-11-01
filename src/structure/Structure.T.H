template<typename T>
std::unique_ptr<List<List<T>>> Foam::Structure::broadcastHaloFields
(
    const GeometricField<T,fvPatchField,volMesh>& fieldData
)
{
    auto haloFieldPtr = std::make_unique<List<List<T>>>(Pstream::nProcs());
    List<List<T>>& haloField = *haloFieldPtr;
    List<T>& ownHaloField = haloField[Pstream::myProcNo()];
    
    const DynamicList<CellDescription>& ownHaloCells = getHaloCellList(Pstream::myProcNo());    
    ownHaloField.setSize(ownHaloCells.size());
    
    for(label haloCellInd=0; haloCellInd<ownHaloCells.size(); haloCellInd++)
    {
        ownHaloField[haloCellInd] = fieldData[ownHaloCells[haloCellInd].index];
    }
    
    Pstream::gatherList(haloField);
    Pstream::scatterList(haloField);
    
    return haloFieldPtr;
}

template<typename T>
void Foam::exchangeBetweenTwo
(
    List<T>& globalData,
    const label comm
)
{
    if(globalData.size()!=Pstream::nProcs())
        FatalErrorInFunction<<"Invalid globalData size"<<exit(FatalError);
    
    if(Pstream::nProcs(comm)!=2)
        FatalErrorInFunction<<"Invalid comm size"<<exit(FatalError);
    
    label worldCommProc0 = Pstream::baseProcNo(comm,0);
    label worldCommProc1 = Pstream::baseProcNo(comm,1);
    
    if(Pstream::myProcNo()==worldCommProc0 || Pstream::myProcNo()==worldCommProc1)
    {
        label localProcNo = Pstream::procNo(comm,Pstream::myProcNo());
        label globalProcNo = Pstream::myProcNo();
        
        label otherLocalProcNo,otherGlobalProcNo;
        if(Pstream::myProcNo()==worldCommProc0)
        {
            otherLocalProcNo = Pstream::procNo(comm,worldCommProc1);
            otherGlobalProcNo = worldCommProc1;
        }
        else
        {
            otherLocalProcNo = Pstream::procNo(comm,worldCommProc0);
            otherGlobalProcNo = worldCommProc0;
        }
        
        List<T> exchangeBuffer(2);
        exchangeBuffer[localProcNo] = globalData[globalProcNo];
    
        Pstream::gatherList(exchangeBuffer,Pstream::msgType(),comm);
        Pstream::scatterList(exchangeBuffer,Pstream::msgType(),comm);
        
        globalData[otherGlobalProcNo] = exchangeBuffer[otherLocalProcNo];
    }
}

template<typename T>
void Foam::Structure::exchangeBetweenAll
(
    List<T>& globalData
)
{
    DynamicList<label> times;
    if( std::is_same<T,DynamicList<List<DynamicList<label>>>>::value ||
        std::is_same<T,DynamicList<List<DynamicList<vector>>>>::value ||
        std::is_same<T,DynamicList<List<DynamicList<scalar>>>>::value
      )
    {
        /*
        List<label> dynSize(globalData.size());
        dynSize[Pstream::myProcNo()] = globalData[Pstream::myProcNo()].size();
        exchangeBetweenAll(dynSize);
        
        List<List<label>> dynListSizes(globalData.size());
        dynListSizes[Pstream::myProcNo()] = List<label>(dynSize[Pstream::myProcNo()]);
        for(label i=0; i<dynListSizes[Pstream::myProcNo()].size(); i++)
            dynListSizes[Pstream::myProcNo()][i] = globalData[Pstream::myProcNo()][i].size();
        exchangeBetweenAll(dynListSizes);
        
        List<List<List<label>>> dynListDynSizes(globalData.size());
        dynListDynSizes[Pstream::myProcNo()] = List<List<label>>(dynSize[Pstream::myProcNo()]);
        for(label i=0; i<dynListDynSizes[Pstream::myProcNo()].size(); i++)
            dynListDynSizes[Pstream::myProcNo()][i].setSize(globalData[Pstream::myProcNo()][i].size());
        
        exchangeBetweenAll(dynListDynSizes);
        */
    }
    else if(std::is_same<T,DynamicList<List<DynamicList<scalar>>>>::value)
    {
        
    }
    else
    {
        if(exchangeSpecific)
        {
            for(label comm : orderedComms)
            {
                auto t1 = std::chrono::high_resolution_clock::now();
                exchangeBetweenTwo<T>(globalData,comm);
                auto t2 = std::chrono::high_resolution_clock::now();
                times.append(std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());
            }
        }
        else
        {
            Pstream::gatherList(globalData);
            Pstream::scatterList(globalData);
        }
    }
    Pout<<"times:"<<times<<Foam::nl;
}
