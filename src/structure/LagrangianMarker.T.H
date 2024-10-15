template<typename VALTYPE>
VALTYPE Foam::LagrangianMarker::convolute
(
    std::function<scalar(vector,vector)> weightFunction,
    std::function<VALTYPE(Pair<label>)> valueFunction
) const
{
    return convolute<VALTYPE>(weightFunction,valueFunction,fullSupport);
}

template<typename VALTYPE>
VALTYPE Foam::LagrangianMarker::convolute
(
    std::function<scalar(vector,vector)> weightFunction,
    std::function<VALTYPE(Pair<label>)> valueFunction,
    const List<Pair<label>>& support
) const
{
    vector X = getMarkerPosition();
    VALTYPE convValue = Foam::zero();
        
    for(const Pair<label>& suppCell : support)
    {
        vector cellCentre;
        scalar cellVolume;
        //auto t1_loop = std::chrono::high_resolution_clock::now();
        getCellData(suppCell,cellCentre,cellVolume);
        //auto t2_loop = std::chrono::high_resolution_clock::now();
        //Info<<support.size()<<"  convolute loop 1-2 took "<<std::chrono::duration_cast<std::chrono::nanoseconds>(t2_loop-t1_loop).count()<<" nanoseconds"<<Foam::nl;
        vector x = cellCentre;

        scalar w = weightFunction(X,x);
        
        /*
        if(w<1e-20)
            Info<<"x:"<<x<<" !w:"<<w<<"  "<<suppCell<<Foam::nl;
        else
            Info<<"x:"<<x<<"  w:"<<w<<"  "<<suppCell<<Foam::nl;
        */
        
        //auto t3_loop = std::chrono::high_resolution_clock::now();
        //Info<<support.size()<<"  convolute loop 2-3 took "<<std::chrono::duration_cast<std::chrono::nanoseconds>(t3_loop-t2_loop).count()<<" nanoseconds"<<Foam::nl;

        VALTYPE f = valueFunction(suppCell);
        //auto t4_loop = std::chrono::high_resolution_clock::now();
        //Info<<support.size()<<"  convolute loop 3-4 took "<<std::chrono::duration_cast<std::chrono::nanoseconds>(t4_loop-t3_loop).count()<<" nanoseconds"<<Foam::nl;
        
        convValue += w*f*cellVolume;
    }    
    return convValue;
}
