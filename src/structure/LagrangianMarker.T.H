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

template<uint dim>
void Foam::LagrangianMarker::solveWeights
(
    std::unique_ptr<std::array<scalar,dim>> Hptr,
    std::unique_ptr<gismo::gsMatrix<scalar>> Mptr,
    gismo::gsMatrix<scalar>& x,
    WeightSystemSolve solutionStrategy
)
{
    std::array<scalar,dim>& H = *Hptr;
    gismo::gsMatrix<scalar>& M = *Mptr;
    x = gismo::gsMatrix<scalar>(dim,1);
    gismo::gsMatrix<scalar> e(dim,1);
    for(uint i=0; i<dim; i++)
    {
        e(i,0) = 0;
        x(i,0) = 0;
    }
    e(0,0) = 1;
    
    gismo::gsMatrix<scalar> P(dim,dim);
    for(uint r=0;r<dim;r++)
        for(uint c=0;c<dim;c++)
            P(r,c)=0;
    if(solutionStrategy==WeightSystemSolve::Raw)
    {
        for(uint d=0; d<dim; d++)
            P(d,d) = 1;
    }
    else if(solutionStrategy==WeightSystemSolve::GeoRescaled)
    {
        for(uint d=0; d<dim; d++)
            P(d,d) = H[d];
    }
    else
    {
        for(uint d=0; d<dim; d++)
            P(d,d) = 1.0/M(d,d);
    }
    gismo::gsMatrix<scalar> A = P*M;
    
    /*
    Pout<<"marker:"<<to_string()<<Foam::endl;
    Pout<<"existingDims:"<<existingDims<<Foam::endl;
    Pout<<"dilation:"<<dilation<<Foam::endl;
    Pout<<"directNeighbours:"<<directNeighbours<<Foam::endl;
    Pout<<"supportCells:"<<supportCells<<Foam::endl;
    */

    std::cout<<"A:"<<A<<std::endl;
    std::cout<<"e:"<<e<<std::endl;
    linearSolve(A,x,e);
    std::cout<<"x:"<<x<<std::endl;
    
    //FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
}
