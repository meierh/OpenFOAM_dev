template<typename VALTYPE>
VALTYPE Foam::LagrangianMarker::convolute
(
    std::function<scalar(vector,vector)> weightFunction,
    std::function<VALTYPE(std::tuple<bool,label,label>)> valueFunction
) const
{
    vector X = getMarkerPosition();
    const DynamicList<std::tuple<bool,label,label>>& supportCells = getSupportCells();
    
    VALTYPE convValue = Foam::zero();
    for(label i=0; i<supportCells.size(); i++)
    {
        const std::tuple<bool,label,label>& suppCellData = supportCells[i];
        vector cellCentre;
        scalar cellVolume;
        getCellData(suppCellData,cellCentre,cellVolume);
        vector x = cellCentre;

        scalar w = weightFunction(X,x);
        VALTYPE f = valueFunction(suppCellData);
        
        convValue += w*f*cellVolume;
    }
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

    linearSolve(A,x,e);
}
