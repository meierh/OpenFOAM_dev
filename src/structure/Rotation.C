#include "Structure.H"

Foam::Quaternion::Quaternion
(
    scalar x,
    scalar y,
    scalar z,
    scalar w
):
x(x),
y(y),
z(z),
w(w)
{
}

Foam::Quaternion::Quaternion
(
    const gsMatrix<scalar>& gsQuaternion
)
{
    if(gsQuaternion.rows()!=4 || gsQuaternion.cols()!=1)
        FatalErrorInFunction<<"Invalid size of gsQuaternion"<<exit(FatalError);
    
    x = gsQuaternion(0,0);
    y = gsQuaternion(1,0);
    z = gsQuaternion(2,0);
    w = gsQuaternion(3,0);
}
        
Foam::Quaternion Foam::Quaternion::operator*
(
    Quaternion const& q
) const 
{
    Quaternion result;
    result.x = w*q.x - x*q.w - y*q.z - z*q.y;
    result.y = w*q.y - x*q.z - y*q.w - z*q.x;
    result.z = w*q.z - x*q.y - y*q.x - z*q.w;
    result.w = w*q.w - x*q.x - y*q.y - z*q.z;
    return result;
}

Foam::Quaternion Foam::Quaternion::operator/
(
    Quaternion const& q
) const
{
    Quaternion invQ = q.invert();
    return (*this)*invQ;
}

Foam::Quaternion Foam::Quaternion::operator-
(
    Quaternion const& q
) const
{
    Quaternion result;
    for(label i=0; i<4; i++)
        result[i] = (*this)[i]-q[i];
    return result;
}

Foam::Quaternion Foam::Quaternion::invert() const
{
    Quaternion invQ = *this;
    scalar absInvQ = invQ.len();
    absInvQ *= absInvQ;
    if(absInvQ<1e-10)
        FatalErrorInFunction<<"Invalid quaternion length"<<exit(FatalError);
    invQ.w /=  absInvQ;
    invQ.x /= -absInvQ;
    invQ.y /= -absInvQ;
    invQ.z /= -absInvQ;
    return invQ;
}

Foam::scalar Foam::Quaternion::len() const
{
    return std::sqrt(x*x + y*y + z*z + w*w);
}

Foam::scalar Foam::Quaternion::distanceNorm2
(
    Quaternion const& q
) const
{
    Quaternion dq;
    dq.w = w-q.w;
    dq.x = x-q.x;
    dq.y = y-q.y;
    dq.z = z-q.z;
    return dq.len();
}

void Foam::Quaternion::normalize()
{
    scalar len = this->len();
    if(len==0)
    {
        x = 1;
    }
    else
    {
        w /= len;
        x /= len;
        y /= len;
        z /= len;
    }
}

Foam::scalar& Foam::Quaternion::operator[]
(
    uint index
)
{
    switch(index)
    {
        case 0:
            return w;
        case 1:
            return x;
        case 2:
            return y;
        case 3:
            return z;
        default:
            FatalErrorInFunction<<"Invalid index in quaternion"<<exit(FatalError);
            return w;
    }
}

Foam::scalar Foam::Quaternion::operator[]
(
    uint index
) const
{
    switch(index)
    {
        case 0:
            return w;
        case 1:
            return x;
        case 2:
            return y;
        case 3:
            return z;
        default:
            FatalErrorInFunction<<"Invalid index in quaternion"<<exit(FatalError);
            return w;
    }
}

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    Quaternion const& q
)
{
    return os << "[("<<q.x<<","<<q.y<<","<<q.z<<")("<<q.w<<")]";
}


Foam::Rotation::Rotation
(
    vector d1,
    vector d2,
    vector d3
)
{
    T[0] = d1;
    T[1] = d2;
    T[2] = d3;
}

Foam::Rotation::Rotation
(
    const Quaternion& q
)
{
    vector d1 = vector
    (
        q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z,
         2*q.w*q.z + 2*q.x*q.y,
        -2*q.w*q.y + 2*q.x*q.z
    );
    vector d2 = vector
    (
        -2*q.w*q.z + 2*q.x*q.y,
        q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z,
         2*q.w*q.x + 2*q.y*q.z
    );
    vector d3 = vector
    (
         2*q.w*q.y + 2*q.x*q.z,
        -2*q.w*q.x + 2*q.y*q.z,
        q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z
    );
    T = {d1,d2,d3};
}

Foam::Rotation Foam::Rotation::operator-
(
    Rotation const& R
) const
{
    Rotation result;
    for(label d=0; d<3; d++)
    {
        result.T[d] = (T[d]-R.T[d]); 
    }
    return result;
}

Foam::Rotation Foam::Rotation::operator+
(
    Rotation const& R
) const
{
    Rotation result;
    for(label d=0; d<3; d++)
    {
        result.T[d] = (T[d]+R.T[d]); 
    }
    return result;
}

bool Foam::Rotation::operator!=
(
    Rotation const& R
) const
{
    return this->T!=R.T;
}

Foam::Rotation Foam::Rotation::operator/
(
    scalar alpha
) const
{
    Rotation result = *this;
    for(label d=0; d<3; d++)
    {
        result.T[d] /= alpha; 
    }
    return result;
}

Foam::scalar Foam::Rotation::distanceNorm2
(
    Rotation const& R
) const
{
    Rotation diff = *this - R;
    return diff.norm2();
}

Foam::scalar Foam::Rotation::norm2() const
{
    scalar sum = 0;
    for(label i=0; i<3; i++)
        for(label j=0; j<3; j++)
            sum += T[i][j]*T[i][j];
    return std::sqrt(sum);
}

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    Rotation const& m
)
{
    return os << m.T;
}

Foam::Rotation Foam::Rotation::compute_dRdX
(
    const Quaternion& dqdX,
    const Quaternion& q
)
{   
    const FixedList<FixedList<vector,4>,3> dRdq = compute_dRdq(q);
 
    Rotation dRdC;
    std::vector<vector*> ddkdCPtr = {&(dRdC.T[0]),&(dRdC.T[1]),&(dRdC.T[2])};
    for(label dk=0; dk<3; dk++)
    {
        const FixedList<vector,4>& dRkdq = dRdq[dk];
        vector& ddkdC = *(ddkdCPtr[dk]);
        for(label dim=0; dim<3; dim++)
        {
            ddkdC[dim] = 0;
            ddkdC[dim] += dRkdq[0][dim] * dqdX.qw();
            ddkdC[dim] += dRkdq[1][dim] * dqdX.qx();
            ddkdC[dim] += dRkdq[2][dim] * dqdX.qy();
            ddkdC[dim] += dRkdq[3][dim] * dqdX.qz();
        }
    }
    return dRdC;
}

Foam::FixedList<Foam::FixedList<Foam::vector,4>,3> Foam::Rotation::compute_dRdq
(
    const Quaternion& q
)
{
    FixedList<vector,4> dd1dq;
        dd1dq[0][0]= 2*q.qw(); dd1dq[1][0]= 2*q.qx(); dd1dq[2][0]=-2*q.qy(); dd1dq[3][0]=-2*q.qz();
        dd1dq[0][1]= 2*q.qz(); dd1dq[1][1]= 2*q.qy(); dd1dq[2][1]= 2*q.qx(); dd1dq[3][1]= 2*q.qw();
        dd1dq[0][2]=-2*q.qy(); dd1dq[1][2]= 2*q.qz(); dd1dq[2][2]=-2*q.qw(); dd1dq[3][2]= 2*q.qx();
    
    FixedList<vector,4> dd2dq;
        dd2dq[0][0]=-2*q.qz(); dd2dq[1][0]= 2*q.qy(); dd2dq[2][0]= 2*q.qx(); dd2dq[3][0]=-2*q.qw();
        dd2dq[0][1]= 2*q.qw(); dd2dq[1][1]=-2*q.qx(); dd2dq[2][1]= 2*q.qy(); dd2dq[3][1]=-2*q.qz();
        dd2dq[0][2]= 2*q.qx(); dd2dq[1][2]= 2*q.qw(); dd2dq[2][2]= 2*q.qz(); dd2dq[3][2]= 2*q.qy();
        
    FixedList<vector,4> dd3dq;
        dd3dq[0][0]= 2*q.qy(); dd3dq[1][0]= 2*q.qz(); dd3dq[2][0]= 2*q.qw(); dd3dq[3][0]= 2*q.qx();
        dd3dq[0][1]=-2*q.qx(); dd3dq[1][1]=-2*q.qw(); dd3dq[2][1]= 2*q.qz(); dd3dq[3][1]= 2*q.qy();
        dd3dq[0][2]= 2*q.qw(); dd3dq[1][2]=-2*q.qx(); dd3dq[2][2]=-2*q.qy(); dd3dq[3][2]= 2*q.qz();
        
    return {dd1dq,dd2dq,dd3dq};
}

Foam::Rotation Foam::Rotation::compute_d2RdX
(
    const Quaternion& d2qdX,
    const Quaternion& dqdX,
    const Quaternion& q
)
{      
    const FixedList<FixedList<vector,4>,3> d2Rdq = compute_d2Rdq(q);
 
    Rotation d2RdX;
    std::vector<vector*> d2dkdXPtr = {&(d2RdX.T[0]),&(d2RdX.T[1]),&(d2RdX.T[2])};
    for(label dk=0; dk<3; dk++)
    {
        const FixedList<vector,4>& d2Rkdq = d2Rdq[dk];
        vector& d2dkdX = *(d2dkdXPtr[dk]);
        for(label dim=0; dim<3; dim++)
        {
            d2dkdX[dim] = 0;
            d2dkdX[dim] += d2Rkdq[0][dim] * (dqdX.qw()*dqdX.qw());
            d2dkdX[dim] += d2Rkdq[1][dim] * (dqdX.qx()*dqdX.qx());
            d2dkdX[dim] += d2Rkdq[2][dim] * (dqdX.qy()*dqdX.qy());
            d2dkdX[dim] += d2Rkdq[3][dim] * (dqdX.qz()*dqdX.qz());
        }
    }
    
    d2RdX = d2RdX + compute_dRdX(d2qdX,q);
    return d2RdX;
}

Foam::FixedList<Foam::FixedList<Foam::vector,4>,3> Foam::Rotation::compute_d2Rdq
(
    const Quaternion& q
)
{
    FixedList<vector,4> d2d1dq;
        d2d1dq[0][0]= 2; d2d1dq[1][0]= 2; d2d1dq[2][0]=-2; d2d1dq[3][0]=-2;
        d2d1dq[0][1]= 0; d2d1dq[1][1]= 0; d2d1dq[2][1]= 0; d2d1dq[3][1]= 0;
        d2d1dq[0][2]= 0; d2d1dq[1][2]= 0; d2d1dq[2][2]= 0; d2d1dq[3][2]= 0;
    
    FixedList<vector,4> d2d2dq;
        d2d2dq[0][0]= 0; d2d2dq[1][0]= 0; d2d2dq[2][0]= 0; d2d2dq[3][0]= 0;
        d2d2dq[0][1]= 2; d2d2dq[1][1]=-2; d2d2dq[2][1]= 2; d2d2dq[3][1]=-2;
        d2d2dq[0][2]= 0; d2d2dq[1][2]= 0; d2d2dq[2][2]= 0; d2d2dq[3][2]= 0;
        
    FixedList<vector,4> d2d3dq;
        d2d3dq[0][0]= 0; d2d3dq[1][0]= 0; d2d3dq[2][0]= 0; d2d3dq[3][0]= 0;
        d2d3dq[0][1]= 0; d2d3dq[1][1]= 0; d2d3dq[2][1]= 0; d2d3dq[3][1]= 0;
        d2d3dq[0][2]= 2; d2d3dq[1][2]=-2; d2d3dq[2][2]=-2; d2d3dq[3][2]= 2;
        
    return {d2d1dq,d2d2dq,d2d3dq};
}
