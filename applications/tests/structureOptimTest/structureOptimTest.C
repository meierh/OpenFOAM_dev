/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    structureOptimTest

Description
    Testing of derivative rod methods

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "icoAdjointImmersedBoundary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar RandomFloat(scalar lower, scalar upper)
{
    scalar random = static_cast<scalar>(rand()) / static_cast<scalar>(RAND_MAX);
    scalar diff = upper - lower;
    scalar delta = random * diff;
    return lower + delta;
}

label RandomInt(label lower, label upper)
{
    return static_cast<label>(RandomFloat(lower,upper));
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<<"Start testing"<<Foam::endl;
    
    FixedList<scalar,4> quaternions = {0.000234,-0.000234,-0.707,0.707};
    Info<<Structure::quaternionsToRotation(quaternions)<<Foam::endl;
    
    
    //FatalErrorInFunction<<"Temp stop"<<exit(FatalError);
    
    for(label run=0; run<1; run++)
    {
        scalar a0 = RandomFloat(0,1);
        label k = RandomInt(0,10);
        Info<<"----- Run:"<<run<<" -----:"<<k<<Foam::endl;

        std::vector<scalar> ak(k);
        std::vector<scalar> bk(k);
        for(label i=0; i<k; i++)
        {
            ak[i] = RandomFloat(0,1);
            bk[i] = RandomFloat(0,1);
        }
        scalar phase = RandomFloat(0,6.2);
        std::cout<<"a0:"<<a0<<std::endl;
        std::cout<<"ak: ";
        for(scalar a: ak)
            std::cout<<a<<" ";
        std::cout<<std::endl;
        std::cout<<"bk: ";
        for(scalar b: bk)
            std::cout<<b<<" ";
        std::cout<<std::endl;
        std::cout<<"phase:"<<phase<<std::endl;
        
        List<List<vector>> curveCoeffs(1);
        curveCoeffs[0] = List<vector>(4);
        scalar eps;
        eps = RandomFloat(-1,1);
        curveCoeffs[0][0] = vector(0+eps,0+eps,0+eps);
        eps = RandomFloat(-1,1);
        curveCoeffs[0][1] = vector(10+eps,10+eps,-1+eps);
        eps = RandomFloat(-1,1);
        curveCoeffs[0][2] = vector(20+eps,0+eps,1+eps);
        eps = RandomFloat(-1,1);
        curveCoeffs[0][3] = vector(40+eps,-10+eps,0+eps);
        Info<<"curveCoeffs:"<<curveCoeffs<<Foam::endl;
        
        std::vector<CrossSection> crossSecList = {CrossSection(a0,ak,bk,phase)};
        CrossSectionStructure testStructure(mesh,crossSecList,true);

        gsNurbs<scalar> derivRotation = testStructure.getCoeffDerivedQuaternions(0,1,0);
        List<FixedList<vector,3>> rotationCurveCoefs = testStructure.getRotationCurveCoefs(0);
        Info<<"rotationCurveCoefs:"<<rotationCurveCoefs<<Foam::endl;
                
        scalar epsilon=1e-1;
        scalar coeffBasicValue = testStructure.getCurveCoeff(0,1,0);
                    
        scalar lowerCoeffValue = coeffBasicValue-epsilon;
        testStructure.setCurveCoeff(0,1,0,lowerCoeffValue);
        rotationCurveCoefs = testStructure.getRotationCurveCoefs(0);
        Info<<"rotationCurveCoefs:"<<rotationCurveCoefs<<Foam::endl;
        gsNurbs<scalar> lowerQuaternions = testStructure.getQuaternions(0);
        //std::cout<<lowerQuaternions.coefs()<<std::endl;
        
        scalar upperCoeffValue = coeffBasicValue+epsilon;                    
        testStructure.setCurveCoeff(0,1,0,upperCoeffValue);
        rotationCurveCoefs = testStructure.getRotationCurveCoefs(0);
        Info<<"rotationCurveCoefs:"<<rotationCurveCoefs<<Foam::endl;
        gsNurbs<scalar> upperQuaternions = testStructure.getQuaternions(0);
        //std::cout<<upperQuaternions.coefs()<<std::endl;

        Info<<Foam::endl;
        gsMatrix<scalar> dQuatdeps = (upperQuaternions.coefs()-lowerQuaternions.coefs())/2*epsilon;
        //std::cout<<dQuatdeps<<std::endl;
    
        //testStructure.selfCheck();
        
        /*
        testStructure.printCurves();
        
        curveCoeffs[0][0] = vector(0,0,0);
        curveCoeffs[0][1] = vector(-10,10,-1);
        curveCoeffs[0][2] = vector(-20,0,1);
        curveCoeffs[0][3] = vector(-40,-10,0);
        testStructure.setCurveCoeffs(curveCoeffs);
        testStructure.printCurves();
        
        
        Info<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<Foam::endl;
        Info<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<Foam::endl;
        Info<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<Foam::endl;

        scalar epsilon=1;
        scalar coeffBasicValue = testStructure.getCurveCoeff(0,1,0);
        Info<<"coeffBasicValue:"<<coeffBasicValue<<Foam::endl;
        
        scalar lowerCoeffValue = coeffBasicValue-epsilon;
        Info<<"lowerCoeffValue:"<<lowerCoeffValue<<Foam::endl;
        testStructure.setCurveCoeff(0,1,0,lowerCoeffValue);
        FixedList<scalar,4> lower_q = testStructure.m_Rot_Eval(0,0.3);
        Info<<"lower_q:"<<lower_q<<Foam::endl;

        scalar upperCoeffValue = coeffBasicValue+epsilon;
        Info<<"upperCoeffValue:"<<upperCoeffValue<<Foam::endl;
        testStructure.setCurveCoeff(0,1,0,upperCoeffValue);
        FixedList<scalar,4> upper_q = testStructure.m_Rot_Eval(0,0.3);
        Info<<"upper_q:"<<upper_q<<Foam::endl;
                   
        FixedList<scalar,4> diff = {upper_q[0]-lower_q[0],upper_q[1]-lower_q[1],upper_q[2]-lower_q[2],upper_q[3]-lower_q[3]};
        scalar diffCoeff = upperCoeffValue-lowerCoeffValue;
        
        Info<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<Foam::endl;
        Info<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<Foam::endl;
        Info<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<Foam::endl;
        
        
        testStructure.setCurveCoeff(0,1,0,-7);
        testStructure.printCurves();
        testStructure.setCurveCoeff(0,1,0,-13);
        testStructure.printCurves();
        */
        
        FatalErrorInFunction<<"Temp stop"<<exit(FatalError);
    }
    
    Info<<"Done testing"<<Foam::endl;    
    return 0;
}



// ************************************************************************* //
