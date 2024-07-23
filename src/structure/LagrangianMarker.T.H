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
    //Pout<<"convolute:"<<support<<Foam::endl;
    vector X = getMarkerPosition();
    VALTYPE convValue = Foam::zero();
    for(const Pair<label>& suppCell : support)
    {
        vector cellCentre;
        scalar cellVolume;
        getCellData(suppCell,cellCentre,cellVolume);
        vector x = cellCentre;

        scalar w = weightFunction(X,x);
        VALTYPE f = valueFunction(suppCell);
        
        convValue += w*f*cellVolume;
    }
    //Pout<<"convolute val:"<<convValue<<Foam::endl;
    return convValue;
}

/*
void Foam::solveWeights
(
    std::unique_ptr<Pair<gismo::gsMatrix<scalar>>> system,
    gismo::gsMatrix<scalar>& x
)
{
    linearSolve(system->first(),x,system->second());
}
*/
