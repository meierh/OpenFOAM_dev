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
