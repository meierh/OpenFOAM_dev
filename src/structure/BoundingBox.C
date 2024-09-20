#include "BoundingBox.H"

Foam::BoundingBox::BoundingBox
(
    vector smaller,
    vector larger
):
smaller(smaller),
larger(larger)
{
    empty = false;
    for(label dim=0; dim<3; dim++)
        if(smaller[dim] >= larger[dim])
            empty = true;
}

Foam::BoundingBox Foam::BoundingBox::operator+
(
    const BoundingBox& rhs
) const
{
    BoundingBox result;
    result.smaller = this->smaller + rhs.smaller;
    result.larger = this->larger + rhs.larger;
    return result;
}

Foam::BoundingBox Foam::BoundingBox::boundingBoxUnion
(
    const BoundingBox& rhs
) const
{
    vector smaller,larger;
    for(label dim=0; dim<3; dim++)
    {
        smaller[dim] = std::min(this->smaller[dim],rhs.smaller[dim]);
        larger[dim] = std::max(this->larger[dim],rhs.larger[dim]);
    }
    BoundingBox result(smaller,larger);
    return result;
}

Foam::BoundingBox Foam::BoundingBox::boundingBoxIntersection
(
    const BoundingBox& rhs
) const
{
    vector smaller,larger;
    for(label dim=0; dim<3; dim++)
    {
        smaller[dim] = std::max(this->smaller[dim],rhs.smaller[dim]);
        larger[dim] = std::min(this->larger[dim],rhs.larger[dim]);
    }
    BoundingBox result(smaller,larger);
    return result;
}

bool Foam::BoundingBox::boundingBoxOverlap
(
    const BoundingBox& rhs
) const
{
    return boundingBoxIntersection(rhs).isEmpty();
}

void Foam::BoundingBox::enlarge
(
    scalar size
)
{
    for(label d=0; d<3; d++)
    {
        larger[d]+=size;
        smaller[d]-=size;
    }
}

bool Foam::BoundingBox::inside
(
    vector point,
    scalar eps    
) const
{
    bool inside = true;
    for(label d=0; d<3; d++)
    {
        if(point[d]>larger[d]+eps)
            inside = false;
        if(point[d]<=smaller[d]+eps)
            inside = false;
    }
    return inside;
}

Foam::scalar Foam::BoundingBox::innerSize()
{
    vector diag = larger-smaller;
    return std::sqrt(diag&diag);
}

Foam::FixedList<Foam::vector,8> Foam::BoundingBox::allVertices()
{
    FixedList<vector,8> allVertices;
    
    allVertices[0][0] = allVertices[1][0] = allVertices[2][0] = allVertices[3][0] = smaller[0];
    allVertices[4][0] = allVertices[5][0] = allVertices[6][0] = allVertices[7][0] = larger[0];
    
    allVertices[0][1] = allVertices[1][1] = allVertices[4][1] = allVertices[5][1] = smaller[1];
    allVertices[2][1] = allVertices[3][1] = allVertices[6][1] = allVertices[7][1] = larger[1];
    
    allVertices[0][2] = allVertices[3][2] = allVertices[4][2] = allVertices[7][2] = smaller[2];
    allVertices[1][2] = allVertices[2][2] = allVertices[5][2] = allVertices[6][2] = larger[2];  
    
    return allVertices;
}

Foam::BoundingBox Foam::BoundingBox::boundsOfCoefficients
(
    const gsMatrix<scalar>& coefs
)
{
    if(coefs.cols()!=3 || coefs.rows()<1)
    {
        std::cout<<"coefs:"<<coefs<<std::endl;
        FatalErrorInFunction<<"Invalid coefs size"<<exit(FatalError);
    }
    vector lowerCurve,upperCurve;
    lowerCurve = upperCurve = vector(coefs(0,0),coefs(0,1),coefs(0,2));
    for(label col=0; col<coefs.cols(); col++)
    {
        for(label row=0; row<coefs.rows(); row++)
        {
            if(lowerCurve[col]>coefs(row,col))
                lowerCurve[col] = coefs(row,col);
            if(upperCurve[col]<coefs(row,col))
                upperCurve[col] = coefs(row,col);
        }
    }
    BoundingBox bb(lowerCurve,upperCurve);
    return bb;
}

Foam::BoundingBox Foam::BoundingBox::boundsOfNurbs
(
    const gsNurbs<scalar>& curve
)
{
    return BoundingBox::boundsOfCoefficients(curve.coefs());
}

Foam::BoundingBox Foam::BoundingBox::boundsOfNurbs
(
    gsNurbs<scalar> curve,
    scalar start,
    scalar end
)
{
    std::unordered_set<label> knotSet;
    for(scalar knot : curve.knots())
        knotSet.insert(knot);

    if(knotSet.find(start)==knotSet.end())
        curve.insertKnot(start);
    if(knotSet.find(end)==knotSet.end())
        curve.insertKnot(end);

    label degree = curve.knots().degree();
    label knot_start = -1;
    label knot_end = -1;
    for(std::size_t knotI=0; knotI<curve.knots().size()-1; knotI++)
    {
        scalar knot_i0 = curve.knots()[knotI];
        scalar knot_i1 = curve.knots()[knotI+1];
        
        if(knot_i0==start && knot_i0!=knot_i1)
        {
            if(knot_start!=-1)
                FatalErrorInFunction<<"Duplicate assigned"<<exit(FatalError);
            knot_start=static_cast<label>(knotI);
        }
        if(knot_i1==end && knot_i0!=knot_i1)
        {
            if(knot_end!=-1)
                FatalErrorInFunction<<"Duplicate assigned"<<exit(FatalError);
            knot_end=static_cast<label>(knotI+1);
        }
    }
    if(knot_start==-1 || knot_end==-1)
        FatalErrorInFunction<<"Not assigned"<<exit(FatalError);

    label coeff_start = knot_start-degree;
    label coeff_end = knot_end-1;

    if(curve.coefs().cols()!=3)
        FatalErrorInFunction<<"Col number must be 3"<<exit(FatalError);
    if(coeff_start<0 || coeff_start>=curve.coefs().rows())
        FatalErrorInFunction<<"Coeff start out of range"<<exit(FatalError);
    if(coeff_end<0 || coeff_end>=curve.coefs().rows())
        FatalErrorInFunction<<"Coeff start out of range"<<exit(FatalError);

    gsMatrix<scalar> coeffs(3,coeff_end-coeff_start+1);
    label ind=0;
    for(label c_s=coeff_start; c_s<coeff_end+1; c_s++,ind++)
    {
        for(label d=0; d<3; d++)
        {
            coeffs(d,ind) = curve.coefs()(d,c_s);
        }
    }
    if(ind!=coeffs.rows())
        FatalErrorInFunction<<"Size mismatch"<<exit(FatalError);
    return BoundingBox::boundsOfCoefficients(coeffs);
}

Foam::BoundingBoxTree& Foam::BoundingBoxTree::operator=
(
    const BoundingBoxTree& rhs
)
{
    root = std::make_unique<Node>(*(rhs.root));
    return *this;
}

void Foam::BoundingBoxTree::findPointParameters
(
    std::vector<scalar>& parameters,
    vector point
) const
{
    std::function<void(const Node*)> recursiveGoDown = [&](const Node* curr)
    {
        if(curr->leftChild && curr->rightChild)
        {
            if(curr->leftChild->value.inside(point))
                recursiveGoDown(curr->leftChild.get());
            if(curr->rightChild->value.inside(point))
                recursiveGoDown(curr->rightChild.get());
        }
        else
        {
            parameters.push_back(curr->key);
            if(curr->leftChild && curr->leftChild->value.inside(point))
                recursiveGoDown(curr->leftChild.get());
            if(curr->rightChild && curr->rightChild->value.inside(point))
                recursiveGoDown(curr->rightChild.get());
        }
    };
    
    if(root)
    {
        if(root->value.inside(point))
        {
            recursiveGoDown(root.get());
        }
    }
}
