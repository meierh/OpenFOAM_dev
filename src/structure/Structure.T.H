template<typename T>
std::unique_ptr<List<List<T>>> Foam::Structure::broadcastHaloFields
(
    const Field<T>& fieldData
)
{
    if(fieldData.size()!=mesh.size())
        FatalErrorInFunction<<"Mismatch in cell number"<<exit(FatalError);
    
    auto haloFieldPtr = std::make_unique<List<List<T>>>(Pstream::nProcs());
    List<List<T>>& haloField = *haloFieldPtr;
    List<T>& ownHaloField = haloField[Pstream::myProcNo()];
    
    const DynamicList<CellDescription>& ownHaloCells = getGlobalHaloCellList_Sorted(Pstream::myProcNo());    
    ownHaloField.setSize(ownHaloCells.size());
    
    for(label haloCellInd=0; haloCellInd<ownHaloCells.size(); haloCellInd++)
    {
        ownHaloField[haloCellInd] = fieldData[ownHaloCells[haloCellInd].index];
    }
    
    exchangeBetweenAll(haloField);
    
    return haloFieldPtr;
}

template<typename T>
void Foam::Structure::exchangeBetweenTwo
(
    List<T>& globalData,
    const label comm
)
{
    if(globalData.size()!=Pstream::nProcs())
        FatalErrorInFunction<<"Invalid globalData size"<<exit(FatalError);
    
    if(Pstream::nProcs(comm)!=2)
        FatalErrorInFunction<<"Invalid comm size"<<exit(FatalError);
    
    label globalProcID_0 = Pstream::baseProcNo(comm,0);
    label globalProcID_1 = Pstream::baseProcNo(comm,1);
    label myGlobalProcID = Pstream::myProcNo();
    label otherGlobalProcID = (myGlobalProcID==globalProcID_0)?globalProcID_1:globalProcID_0;

    if(globalProcID_0!=Pstream::myProcNo() && globalProcID_1!=Pstream::myProcNo())
        return;
    
    label localProcID_0 = Pstream::procNo(comm,globalProcID_0);
    label localProcID_1 = Pstream::procNo(comm,globalProcID_1);
    label myLocalProcID = Pstream::procNo(comm,myGlobalProcID);
    label otherLocalProcID = (myLocalProcID==localProcID_0)?localProcID_1:localProcID_0;

    List<T> exchangeBuffer(2);
    exchangeBuffer[myLocalProcID] = globalData[myGlobalProcID];
    
    Pstream::gatherList(exchangeBuffer,Pstream::msgType(),comm);
    Pstream::scatterList(exchangeBuffer,Pstream::msgType(),comm);
    
    globalData[otherGlobalProcID] = exchangeBuffer[otherLocalProcID];
}

template<typename T>
void Foam::Structure::exchangeBetweenAll
(
    List<T>& globalData
) const
{
    //Pout<<"Structure::exchangeBetweenAll"<< Foam::nl;

    if(exchangeSpecific)
    {
        for(label comm : orderedComms)
        {
            exchangeBetweenTwo<T>(globalData,comm);
        }
    }
    else
    {
        Pstream::gatherList(globalData);
        Pstream::scatterList(globalData);
    }
    //Pout<<"Structure::exchangeBetweenAll done"<< Foam::nl;
}
