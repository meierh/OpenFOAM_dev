#include "CrossSectionStructure.H"
#include <sys/stat.h>

Foam::CrossSectionStructure::CrossSectionStructure
(
    const fvMesh& mesh,
    std::vector<CrossSection> rodCrossSection,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
LineStructure(mesh,modusFieldToMarker,modusMarkerToField),
rodCrossSection(rodCrossSection)
{
    Info<<"CrossSectionStructure"<<Foam::nl;
    initialize();
}

Foam::CrossSectionStructure::CrossSectionStructure
(
    const fvMesh& mesh,
    std::vector<CrossSection> rodCrossSection,
    bool empty
):
LineStructure(mesh),
rodCrossSection(rodCrossSection)
{
    Info<<"CrossSectionStructure empty"<<Foam::nl;
}

Foam::CrossSectionStructure::CrossSectionStructure
(
    const fvMesh& mesh,
    const std::shared_ptr<IOdictionary> structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
LineStructure(mesh,structureDict,true,modusFieldToMarker,modusMarkerToField),
rodCrossSection(createCrossSectionsFromDict(*structureDict))
{
    readRodPntsToMeshSpacingDict(*structureDict);
    initialize();
}

Foam::vector Foam::CrossSectionStructure::evaluateRodVelocity
(
    label rodNumber,
    scalar parameter,
    scalar angle,
    scalar radiusFrac
)
{
    scalar currentTime = mesh.time().value();
    vector currentPosition = evaluateRodCircumPos(Rods[rodNumber],parameter,&(rodCrossSection[rodNumber]),angle,radiusFrac);
    
    const std::pair<gsNurbs<scalar>,scalar>* prevDef = readPrevRodDeformation(rodNumber);
    const std::pair<gsNurbs<scalar>,scalar>* prevRot = readPrevRodRotation(rodNumber);
    
    if(prevDef==nullptr || prevRot==nullptr)
        return vector(0,0,0);
    
    scalar prevTime = prevDef->second;
    if(prevTime!=prevRot->second)
        FatalErrorInFunction<<"Mismatch in prev time stamps"<<exit(FatalError);
    scalar deltaT = currentTime-prevTime;
    
    vector prevR,prevD1,prevD2,prevD3;
    rodEval(Rods[rodNumber]->m_Curve,prevDef->first,prevRot->first,parameter,prevD1,prevD2,prevD3,prevR);
    scalar radius = rodCrossSection[rodNumber](parameter,angle)*radiusFrac;
    vector coordXDir = std::cos(angle)*radius*prevD1;
    vector coordYDir = std::sin(angle)*radius*prevD2;
    vector previousPosition = prevR+coordXDir+coordYDir;
    
    return (currentPosition-previousPosition)/deltaT;
}

void Foam::CrossSectionStructure::to_string()
{
    Info<<"---CrossSectionStructure::Markers---"<<Foam::nl;
    int count=0;
    for(std::unique_ptr<std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>>& oneRodMarkers : rodMarkersList)
    {
        Info<<"Rod "<<count<<Foam::nl;
        if(oneRodMarkers)
        {
            for(std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>& radialMarkers : *oneRodMarkers)
            {
                Info<<"\t-"<<Foam::nl;
                for(std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>& angleMarkers : radialMarkers.second)
                {
                    for(const LagrangianMarkerOnCrossSec& marker : angleMarkers.second)
                    {
                        Info<<"\t"<<marker.to_string()<<Foam::nl;
                    }
                }
            }
        }
    }
    Info<<"------------------------------------"<<Foam::nl;
}

void Foam::CrossSectionStructure::setCrossSecParameters
(
    label rodNumber,
    bool phase,
    label fourierCoeffNumber,
    label derivCoeffNumber,
    scalar value
)
{    
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber"<<exit(FatalError);
    CrossSection& crossSec = rodCrossSection[rodNumber];
    if(phase)
    {
        crossSec.setPhaseNurbsCoeff(derivCoeffNumber,value);
    }
    else
    {
        if(fourierCoeffNumber<0 || fourierCoeffNumber>=crossSec.numberFourierCoeff())
            FatalErrorInFunction<<"Invalid fourierCoeffNumber"<<exit(FatalError);
        crossSec.setFourierCoeffNurbsCoeff(fourierCoeffNumber,derivCoeffNumber,value);
    }
}

Foam::scalar Foam::CrossSectionStructure::getCrossSecParameters
(
    label rodNumber,
    bool phase,
    label fourierCoeffNumber,
    label derivCoeffNumber
)
{    
    if(rodNumber<0 || rodNumber>=nR)
        FatalErrorInFunction<<"Invalid rodNumber"<<exit(FatalError);
    CrossSection& crossSec = rodCrossSection[rodNumber];
    if(phase)
    {
        return crossSec.getPhaseNurbsCoeff(derivCoeffNumber);
    }
    else
    {
        if(fourierCoeffNumber<0 || fourierCoeffNumber>=crossSec.numberFourierCoeff())
            FatalErrorInFunction<<"Invalid fourierCoeffNumber"<<exit(FatalError);
        return crossSec.getFourierCoeffNurbsCoeff(fourierCoeffNumber,derivCoeffNumber);
    }
}

Foam::vector Foam::CrossSectionStructure::dXdParam
(
    const LagrangianMarker* marker,
    const Parameter& par
)
{
    vector rodDerive(0,0,0);
    const LagrangianMarkerOnCrossSec* crossSecMarker;
    if((crossSecMarker = dynamic_cast<const LagrangianMarkerOnCrossSec*>(marker))!=nullptr)
    {
        if(par.getType()!=Parameter::Type::None)
            rodDerive = dXdParam(marker->getRodNumber(),marker->getMarkerParameter(),marker->getMarkerAngle(),marker->getMarkerRadiusFrac(),par);
        else
            FatalErrorInFunction<<"Invalid type of parameter!"<<exit(FatalError);
    }
    else
        FatalErrorInFunction<<"Can not be called with non LagrangianMarkerOnCrossSec marker"<<exit(FatalError);
    return rodDerive;
}

Foam::vector Foam::CrossSectionStructure::dXdParam
(
    label rodNumber,
    scalar rodParameter,
    scalar angle,
    scalar radiusFrac,
    const Parameter& par
)
{
    if(!(par.isValid()))
        FatalErrorInFunction<<"Invalid parameter here!"<<exit(FatalError);
    
    vector d1,d2,d3,r;
    rodEval(Rods[rodNumber],rodParameter,d1,d2,d3,r);
    scalar radius = rodCrossSection[rodNumber](rodParameter,angle)*radiusFrac;

    vector rodDerive(0,0,0);
    if(par.getType()==Parameter::Type::Rod)
    {
        label dimension = par.getDimension();
        const std::vector<NurbsCoeffReference>& nurbsCoeffs = par.getNurbsCoeffs();
        for(const NurbsCoeffReference& ref : nurbsCoeffs)
        {
            if(rodNumber==ref.rodNumber)
            {
                if(dimension!=ref.dimension)
                    FatalErrorInFunction<<"Dimension mismatch!"<<exit(FatalError);
                
                vector d1dC,d2dC,d3dC,rdC;
                rodEvalDerivCoeff(rodNumber,ref.coeffNumber,dimension,rodParameter,d1dC,d2dC,d3dC,rdC);
                
                vector cosd1dC = std::cos(angle)*d1dC;
                vector sind2dC = std::sin(angle)*d2dC;

                rodDerive += ( rdC + (cosd1dC+sind2dC)*radius );
            }
            else
                FatalErrorInFunction<<"RodNumber mismatch!"<<exit(FatalError);
        }
    }
    else if(par.getType()==Parameter::Type::CrossSection)
    {
        vector cosd1 = std::cos(angle)*d1;
        vector sind2 = std::sin(angle)*d2;
        vector radiusDirection = cosd1+sind2;
        scalar derivRadius = 0;
        const std::vector<CrossSectionCoeffReference>& crossSecCoeffs = par.getCrossSecCoeffs();
        for(const CrossSectionCoeffReference& ref : crossSecCoeffs)
        {
            if(rodNumber==ref.rodNumber)
            {
                if(ref.phase)
                {
                    derivRadius +=
                    rodCrossSection[rodNumber].evalRadiusDerivPhaseNurbsCoeff(ref.coeffNumber,
                                                                              rodParameter,angle);
                }
                else
                {
                    derivRadius +=
                    rodCrossSection[rodNumber].evalRadiusDerivFourierCoeffNurbsCoeff(ref.fourierCoeffNumber,
                                                                                     ref.coeffNumber,rodParameter,
                                                                                     angle);
                }
            }
            else
                FatalErrorInFunction<<"RodNumber mismatch!"<<exit(FatalError);
        }
        rodDerive = radiusDirection*derivRadius;
    }
    else
        FatalErrorInFunction<<"Invalid type of parameter here!"<<exit(FatalError);
    rodDerive *= radiusFrac;
    return rodDerive;
}

Foam::List<Foam::scalar> Foam::CrossSectionStructure::getParameterValue
(
    const Parameter& para
)
{
    if(!para.isValid())
        FatalErrorInFunction<<"Invalid parameter!"<<exit(FatalError);
    
    if(para.getType()==Parameter::Type::Rod)
    {
        return LineStructure::getParameterValue(para);
    }
    else if(para.getType()==Parameter::Type::CrossSection)
    {
        List<scalar> values(para.getNurbsCoeffs().size());
        for(label coeffI=0; coeffI<values.size(); coeffI++)
        {
            const NurbsCoeffReference& nurbsCoeff = para.getNurbsCoeffs()[coeffI];
            values[coeffI] = getCurveCoeff(nurbsCoeff.rodNumber,nurbsCoeff.coeffNumber,nurbsCoeff.dimension);
        }
        return values;
    }
    else
        FatalErrorInFunction<<"Invalid parameter type none!"<<exit(FatalError);

    
    List<scalar> values(para.getNurbsCoeffs().size());
    for(label coeffI=0; coeffI<values.size(); coeffI++)
    {
        const NurbsCoeffReference& nurbsCoeff = para.getNurbsCoeffs()[coeffI];
        values[coeffI] = getCurveCoeff(nurbsCoeff.rodNumber,nurbsCoeff.coeffNumber,nurbsCoeff.dimension);
    }
    
    return values;
}

std::vector<Foam::CrossSection> Foam::CrossSectionStructure::createCrossSectionsFromDict
(
    const IOdictionary& structureDict
)
{
    Info<<"createCrossSectionsFromDict"<<Foam::nl;
    
    const dictionary& crossSecDict = structureDict.subDict("crossSections");
    List<keyType> crossSectionsKey = crossSecDict.keys();

    std::vector<CrossSection> crossSec;
    for(keyType oneCrossSec : crossSectionsKey)
    {
        Info<<"Read rod "<<oneCrossSec<<" cross Section"<<Foam::nl;
        const dictionary& oneCrossSecDict = crossSecDict.subDict(oneCrossSec);
        ITstream crossSecTypeStream = oneCrossSecDict.lookup("type");
        token crossSecTypeToken;
        crossSecTypeStream.read(crossSecTypeToken);
        if(!crossSecTypeToken.isWord())
            FatalErrorInFunction<<"Invalid token type in constant/structureDict/crossSections/"<<oneCrossSec<<"/type"<<exit(FatalError);
        word crossSecTypeStr = crossSecTypeToken.wordToken();

        if(crossSecTypeStr=="circle")
        {
            ITstream radiusStream = oneCrossSecDict.lookup("radius");
            token radiusNumber;
            radiusStream.read(radiusNumber);
            if(!radiusNumber.isScalar())
                FatalErrorInFunction<<"Expected scalar but got:"<<radiusNumber<<" at line "<<radiusNumber.lineNumber()<<"in dictionary "<<oneCrossSecDict.name()<<exit(FatalError);
            scalar radius = radiusNumber.scalarToken();
            crossSec.push_back(CrossSection(radius));
        }
        else if(crossSecTypeStr=="cylinder")
        {
            /*
            ITstream a0Stream = oneCrossSecDict.lookup("a0");
            token a0Token;
            a0Stream.read(a0Token);
            if(!a0Token.isScalar())
                FatalErrorInFunction<<"Expected scalar but got:"<<a0Token<<" at line "<<a0Token.lineNumber()<<"in dictionary "<<oneCrossSecDict.name()<<exit(FatalError);
            scalar a0 = a0Token.scalarToken();

            ITstream akStream = oneCrossSecDict.lookup("ak");
            token akToken;
            akStream.read(akToken);
            if(!akToken.isScalar())
                FatalErrorInFunction<<"Expected scalar but got:"<<akToken<<" at line "<<akToken.lineNumber()<<"in dictionary "<<oneCrossSecDict.name()<<exit(FatalError);
            scalar ak = akToken.scalarToken();

            ITstream bkStream = oneCrossSecDict.lookup("bk");
            token bkToken;
            bkStream.read(bkToken);
            if(!bkToken.isScalar())
                FatalErrorInFunction<<"Expected scalar but got:"<<bkToken<<" at line "<<bkToken.lineNumber()<<"in dictionary "<<oneCrossSecDict.name()<<exit(FatalError);
            scalar bk = bkToken.scalarToken();

            bool phaseExists = false;
            scalar phase;
            if(oneCrossSecDict.found("phase"))
            {
                ITstream phaseStream = oneCrossSecDict.lookup("phase");
                token phaseToken;
                phaseStream.read(phaseToken);
                if(!phaseToken.isScalar())
                    FatalErrorInFunction<<"Expected scalar but got:"<<phaseToken<<" at line "<<phaseToken.lineNumber()<<"in dictionary "<<oneCrossSecDict.name()<<exit(FatalError);
                phase = phaseToken.scalarToken();
            }
            */

            FatalErrorInFunction<<"Not yet implemented"<<exit(FatalError);
        }
        else if(crossSecTypeStr=="twistedCylinder")
        {
            FatalErrorInFunction<<"Not yet implemented"<<exit(FatalError);
        }
        else if(crossSecTypeStr=="fullyParam")
        {
            FatalErrorInFunction<<"Not yet implemented"<<exit(FatalError);
        }
        else
            FatalErrorInFunction<<"Invalid CrossSection type:"<<crossSecTypeStr<<" -- {circle,cylinder,twistedCylinder,fullyParam}"<<exit(FatalError);
    }
    if(static_cast<int>(crossSec.size())!=nR)
        FatalErrorInFunction<<"Mismatch in crossSection and rod number"<<exit(FatalError);
    return crossSec;
}

void Foam::CrossSectionStructure::check()
{
    if(!myMesh)
        FatalErrorInFunction<<"Rod Mesh not set!"<<exit(FatalError);
    if(myMesh->m_nR!=static_cast<int>(rodMarkersList.size()))
    {
        rodMarkersList.resize(myMesh->m_nR);
    }   
    if(rodCrossSection.size()!=rodMarkersList.size())
        FatalErrorInFunction<<"Mismatch in size of rodCrossSection and rodMarkersList"<<exit(FatalError);
}

void Foam::CrossSectionStructure::createSpacingPoints()
{    
    CrossSectionStructure::initialRodPoints.resize(myMesh->m_nR);
    LineStructure::createSpacingPoints();
}

void Foam::CrossSectionStructure::createSpacedPointsOnRod
(
    label rodNumber,
    scalar spacing
)
{   
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    CrossSection& crossSec = rodCrossSection[rodNumber];
    
    auto pointsPtr = std::unique_ptr<std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>>
    (
        new std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>()
    );
    std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>& pointsVec = *pointsPtr;
    
    LineStructure::createSpacedPointsOnRod(rodNumber,spacing);
    std::vector<scalar>& spacedPointsOnRod = *(LineStructure::initialRodPoints[rodNumber]);
    std::list<scalar> onRodPoints;
    onRodPoints.insert(onRodPoints.end(),spacedPointsOnRod.begin(),spacedPointsOnRod.end());
    
    //Compute points along the rod
    auto pntIter0 = onRodPoints.begin();
    auto pntIter1 = ++onRodPoints.begin();
    for( ; pntIter1!=onRodPoints.end() ; )
    {
        scalar pnt0Para = *pntIter0;
        scalar pnt1Para = *pntIter1;
        
        vector r,d2,d3;        
        vector pnt0D1;
        rodEval(oneRod,pnt0Para,pnt0D1,d2,d3,r);
        vector pnt1D1;
        rodEval(oneRod,pnt1Para,pnt1D1,d2,d3,r);
        
        if((pnt0D1&pnt1D1) < 0)
            FatalErrorInFunction<<"Vector can not be pointing in different directions"<< exit(FatalError);
        
        scalar angle = std::acos((pnt0D1&pnt1D1)/(std::sqrt(pnt0D1&pnt0D1)*std::sqrt(pnt1D1&pnt1D1)));
        scalar halfAngle = std::abs(angle/2);
        scalar rodDistance = LineStructure::distance(oneRod,pnt0Para,pnt1Para);
        
        scalar maxRadiusPnt0 = crossSec.upperLimitRadius(pnt0Para);
        scalar addedDistPnt0 = maxRadiusPnt0*std::sin(halfAngle);
        scalar maxRadiusPnt1 = crossSec.upperLimitRadius(pnt1Para);
        scalar addedDistPnt1 = maxRadiusPnt1*std::sin(halfAngle);
        
        scalar totalDistance = rodDistance+addedDistPnt0+addedDistPnt1;
        scalar totalDistanceFrac = totalDistance/ (spacing*iniRodPntsDistToMeshSpacing);
        label nbrOfAddIntersec = std::ceil(totalDistanceFrac)-1;
        
        scalar subDistPara = (pnt1Para-pnt0Para)/(nbrOfAddIntersec+1);
        std::vector<scalar> subParas;
        for(label i=1;i<nbrOfAddIntersec+1;i++)
            subParas.push_back(pnt0Para+subDistPara*i);
        onRodPoints.insert(pntIter1,subParas.begin(),subParas.end());
        
        pntIter0 = pntIter1;
        pntIter1++;
    }
    
    //Compute points along the rod end
    pointsVec.resize(onRodPoints.size());
    uint index=0;
    for(auto iter=onRodPoints.begin(); iter!=onRodPoints.end(); iter++,index++)
    {
        scalar parameter = *iter;
        pointsVec[index].first = parameter;
    
        std::vector<std::pair<scalar,std::vector<scalar>>>& radialData = pointsVec[index].second;
                
        if(index==0 || index==onRodPoints.size()-1)
        {
            scalar upperLimitRadius = crossSec.upperLimitRadius(parameter);
            scalar radiusSpacing = upperLimitRadius/(spacing*iniRodPntsDistToMeshSpacing);
            label numberOfFracs = std::ceil(radiusSpacing);
            scalar radFracPart = 1.0/numberOfFracs;
            scalar radFrac = 1.0;
            radialData.resize(numberOfFracs+1);
                        
            for(label r=0; r<numberOfFracs; r++)
            {
                radialData[r].first = radFrac;
                createSpacedPointsOnCrossSec(oneRod,parameter,&crossSec,radFrac,spacing,radialData[r].second);
                //Info<<radFrac<<"  radialData["<<r<<"].second.size():"<<radialData[r].second.size()<<Foam::nl;
                radFrac-=radFracPart;
            }
                        
            radialData[numberOfFracs].first = 0;
            radialData[numberOfFracs].second.resize(1);
            radialData[numberOfFracs].second[0] = 0;
        }
        else
        {
            radialData.resize(1);
            radialData[0].first = 1.0;
            createSpacedPointsOnCrossSec(oneRod,parameter,&crossSec,1.0,spacing,radialData[0].second);
            //Info<<"radialData[0].second.size():"<<radialData[0].second.size()<<Foam::nl;
        }
    }
    initialRodPoints[rodNumber] = std::move(pointsPtr);
}

void Foam::CrossSectionStructure::createSpacedPointsOnCrossSec
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    CrossSection* oneCrossSec,
    scalar radFrac,
    scalar spacing,
    std::vector<scalar>& angleData
)
{    
    std::list<scalar> points;
    points.push_back(0);
    points.push_back(0.5*Foam::constant::mathematical::pi);
    points.push_back(Foam::constant::mathematical::pi);
    points.push_back(1.5*Foam::constant::mathematical::pi);
    points.push_back(2*Foam::constant::mathematical::pi);
    bool refined=true;
    while(refined)
    {
        refined = false;
        auto pntsIter0 = points.begin();
        auto pntsIter1 = ++(points.begin());
        for( ; pntsIter1!=points.end() ; )
        {
            scalar dist = distance(oneRod,parameter,oneCrossSec,*pntsIter0,*pntsIter1,radFrac);
            if(dist > spacing*iniRodPntsDistToMeshSpacing)
            {
                scalar middlePar = 0.5*(*pntsIter0 + *pntsIter1);
                points.insert(pntsIter1,middlePar);
                refined=true;
            }
            pntsIter0 = pntsIter1;
            pntsIter1++;
        }
    }
    for(auto iter=points.begin(); iter!=points.end(); iter++)
    {
        angleData.push_back(*iter);
    }
}

void Foam::CrossSectionStructure::createMarkersFromSpacedPointsOnRod
(
    label rodNumber
)
{  
    if(rodInMesh[rodNumber])
    {
        const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
        CrossSection& crossSec = rodCrossSection[rodNumber];
        
        auto markersPtr = std::unique_ptr<std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>>
        (
            new std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>()
        );
        std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>& oneRodMarkers = *markersPtr;
        
        //[rod] -> {para , [] -> {rad , [] -> angle}}
        const std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>& thisRodIniPnts = *(initialRodPoints[rodNumber]);
        
        oneRodMarkers.resize(thisRodIniPnts.size());
        auto iterPara = oneRodMarkers.begin();
                
        for(uint paraInd=0; paraInd<thisRodIniPnts.size(); paraInd++,iterPara++)
        {
            scalar parameter = thisRodIniPnts[paraInd].first;
            const std::vector<std::pair<scalar,std::vector<scalar>>>& radialPnts = thisRodIniPnts[paraInd].second;
            std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>& radialMarkers = *iterPara;
            radialMarkers.first = parameter;
            radialMarkers.second.resize(radialPnts.size());
            auto iterRad = radialMarkers.second.begin();
            for(uint radInd=0; radInd<radialPnts.size(); radInd++,iterRad++)
            {
                scalar radialFrac = radialPnts[radInd].first;
                const std::vector<scalar>& anglePnts = radialPnts[radInd].second;
                std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>& angleMarkers = *iterRad;
                angleMarkers.first = radialFrac;
                auto iterAngle = angleMarkers.second.begin();
                for(uint angleInd=0; angleInd<anglePnts.size(); angleInd++,iterAngle++)
                {
                    scalar angle = anglePnts[angleInd];
                    angleMarkers.second.push_back
                    (                        
                        LagrangianMarkerOnCrossSec
                        (
                            *this,mesh,rodNumber,oneRod,parameter,&crossSec,angle,radialFrac
                        )
                    );
                }
            }
        }
        rodMarkersList[rodNumber] = std::move(markersPtr);
    }
}

void Foam::CrossSectionStructure::refineMarkersOnRod
(
    label rodNumber,
    bool useMarkerCharLenSpacing,
    std::pair<bool,scalar> forcedSpacing
)
{    
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    CrossSection& crossSec = rodCrossSection[rodNumber];
    if(rodMarkersList[rodNumber])
    {
        std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>& markers = *(rodMarkersList[rodNumber]);
        
        //Info<<"markers.size():"<<markers.size()<<Foam::nl;
        for(auto iterPara=markers.begin(); iterPara!=markers.end(); iterPara++)
        {
            std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>& paraMarkers = *iterPara;
            //Info<<"\t paraMarkers.second.size():"<<paraMarkers.second.size();
            //Info<<"  paraMarkers:"<<paraMarkers.first<<Foam::nl;
            scalar para = paraMarkers.first;
            
            label radialIndex = paraMarkers.second.size()-1;
            for(auto iterRadial=iterPara->second.begin(); iterRadial!=iterPara->second.end(); iterRadial++,radialIndex--)
            {
                std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>& radialMarkers = *iterRadial;
                //Info<<"\t\t radialMarkers.second.size():"<<radialMarkers.second.size();
                //Info<<"  radialMarkers:"<<radialMarkers.first<<Foam::nl;
                scalar radialFrac = radialMarkers.first;
                
                bool isRadialCenter = false;
                if(paraMarkers.second.size()>1 && radialIndex==0)
                    isRadialCenter=true;                
                
                if(isRadialCenter)
                {
                    if(radialMarkers.second.size()==0)
                    {
                        radialMarkers.second.push_back
                        (
                            LagrangianMarkerOnCrossSec
                            (
                                *this,mesh,rodNumber,oneRod,para,&crossSec,0,radialFrac
                            )
                        );
                    }
                }
                else
                {
                    scalar pi = Foam::constant::mathematical::pi;
                    if(radialMarkers.second.size()==0)
                    {
                        radialMarkers.second.push_back
                        (
                            LagrangianMarkerOnCrossSec(*this,mesh,rodNumber,oneRod,para,&crossSec,0,radialFrac)                    
                        );
                        radialMarkers.second.push_back
                        (
                            LagrangianMarkerOnCrossSec(*this,mesh,rodNumber,oneRod,para,&crossSec,0.5*pi,radialFrac)                    
                        );
                        radialMarkers.second.push_back
                         (
                            LagrangianMarkerOnCrossSec(*this,mesh,rodNumber,oneRod,para,&crossSec,pi,radialFrac)                    
                        );
                        radialMarkers.second.push_back
                        (
                            LagrangianMarkerOnCrossSec(*this,mesh,rodNumber,oneRod,para,&crossSec,1.5*pi,radialFrac)                    
                        );
                        radialMarkers.second.push_back
                        (
                            LagrangianMarkerOnCrossSec(*this,mesh,rodNumber,oneRod,para,&crossSec,2*pi,radialFrac)                    
                        );
                    }
                    else if(radialMarkers.second.size()==1)
                    {
                        scalar singleMarkerAngle = radialMarkers.second.front().getMarkerAngle();
                        radialMarkers.second.push_back
                        (
                            LagrangianMarkerOnCrossSec(*this,mesh,rodNumber,oneRod,para,&crossSec,singleMarkerAngle+pi,radialFrac)                    
                        );
                    }
                    refineCircumferential(radialMarkers.second);
                }
            }
        }
        
        // Refine rod endings radially
        refineRadial(markers.front().second,initialMeshSpacing);
        refineRadial(markers.back().second,initialMeshSpacing);
        
        refineTangential(markers,initialMeshSpacing);
        
        /*
        bool refined=true;
        while(refined)
        {
            refined = false;
            bool cond = true;
            auto markersIter0 = markers.begin();
            auto markersIter1 = ++(markers.begin());
            for( ; markersIter1!=markers.end() ; )
            {
                scalar markers0Para = markersIter0->getMarkerParameter();
                scalar markers0CellSpacing = markersIter0->getMarkerCellMinSpacing();
                bool markers0InCell = (markersIter0->getMarkerCell()!=-1);
                
                scalar markers1Para = markersIter0->getMarkerParameter();
                scalar markers1CellSpacing = markersIter0->getMarkerCellMinSpacing();
                bool markers1InCell = (markersIter1->getMarkerCell()!=-1);
                
                scalar dist = distance(oneRod,markers0Para,markers1Para);
                
                bool subdivide = false;
                if(forcedSpacing.first)
                {
                    if(dist>forcedSpacing.second)
                        subdivide=true;
                }
                if(markers0InCell || markers1InCell)
                {
                    scalar minSpacing = std::min(markers0CellSpacing,markers1CellSpacing);
                    if(dist>minSpacing)
                        subdivide=true;
                }
            
                if(subdivide)
                {
                    scalar middlePar = 0.5*(markers0Para+markers1Para);
                    LagrangianMarker middleMarker(*this,mesh,rodNumber,oneRod,middlePar);
                    auto inserted = markers.insert(markersIter1,middleMarker);
                    refined=true;
                }
                markersIter0 = markersIter1;
                markersIter1++;
            }
        }
        */        
    }
}

void Foam::CrossSectionStructure::refineCircumferential
(
    std::list<LagrangianMarkerOnCrossSec>& circumMarkers,
    bool useMarkerCharLenSpacing,
    std::pair<bool,scalar> refineSpacing
)
{
    if(circumMarkers.size()<2)
    {
        Info<<circumMarkers.front().getMarkerRadiusFrac()<<Foam::nl;
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    }
     
    label rodNumber = circumMarkers.front().getRodNumber();
    const ActiveRodMesh::rodCosserat* oneRod = circumMarkers.front().getBaseRod();
    scalar parameter = circumMarkers.front().getMarkerParameter();
    scalar radFrac = circumMarkers.front().getMarkerRadiusFrac();
    CrossSection* oneCrossSec =  &(rodCrossSection[rodNumber]);

    /*
     * ----------------------Initial data collection 1---------------------------- 
     */
    //auto t0 = std::chrono::system_clock::now();
    //auto start = t0;
    /*
    scalar start_avgDeltaAngle=0;
    scalar start_avgDeltaDist=0;
    std::vector<scalar> deltaAngles;
    std::vector<scalar> deltaDist;
    for(auto iter=circumMarkers.begin(); iter!=circumMarkers.end(); iter++)
    {
        std::list<LagrangianMarkerOnCrossSec>::iterator next = iter;
        next++;
        if(next!=circumMarkers.end())
        {
            scalar deltaAngle = next->getMarkerAngle()-iter->getMarkerAngle();
            start_avgDeltaAngle += deltaAngle;
            deltaAngles.push_back(deltaAngle);
            scalar dist = distance
            (
                oneRod,parameter,oneCrossSec,iter->getMarkerAngle(),next->getMarkerAngle(),radFrac
            );
            start_avgDeltaDist+=dist;
            deltaDist.push_back(dist);
        }
        else
            break;
    }
    start_avgDeltaAngle/=circumMarkers.size();
    std::sort(deltaAngles.begin(),deltaAngles.end());
    scalar start_medDeltaAngle = deltaAngles[deltaAngles.size()/2];
    start_avgDeltaDist/=circumMarkers.size();
    std::sort(deltaDist.begin(),deltaDist.end());
    scalar start_medDeltaDist = deltaDist[deltaDist.size()/2];
    label start_size = circumMarkers.size();
    */
    //auto t1 = std::chrono::system_clock::now();
    //Info<<"\t\t IDC1:"<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
    /*
     * -------------------------------------------------------------------------- 
     */
    
    /*
     * ----------------------Initial data collection 2---------------------------- 
     */
    //t0 = std::chrono::system_clock::now();
    //scalar maxMarkerDist = std::numeric_limits<scalar>::min();
    //scalar minMarkerCellSize = std::numeric_limits<scalar>::max();
    //label refineSpanCount=0;
    /*
    for(auto iter=circumMarkers.begin(); iter!=circumMarkers.end(); iter++)
    {
        std::list<LagrangianMarkerOnCrossSec>::iterator next = iter;
        next++;
        if(next==circumMarkers.end())
            next = circumMarkers.begin();
        
        scalar circMarker0Angle = iter->getMarkerAngle();
        scalar circMarker0CellSpacing = iter->getMarkerCellMinSpacing();
        bool circMarker0InCell = (iter->getMarkerCell()!=-1);
            
        scalar circMarker1Angle = next->getMarkerAngle();
        scalar circMarker1CellSpacing = next->getMarkerCellMinSpacing();
        bool circMarker1InCell = (next->getMarkerCell()!=-1);

        
        bool subdivide = false;
        
        scalar dist = lowerBound_distance
        (
            rodNumber,parameter,circMarker0Angle,circMarker1Angle,radFrac
        );
        if(refineSpacing.first)
        {
            if(dist>refineSpacing.second)
                subdivide=true;
        }
        if(circMarker0InCell || circMarker1InCell)
        {
            scalar minSpacing = std::min(circMarker0CellSpacing,circMarker1CellSpacing);
            if(dist>minSpacing)
                subdivide=true;
        }
        if(!subdivide)
        {
            dist = distance
            (
                rodNumber,parameter,circMarker0Angle,circMarker1Angle,radFrac
            );
            if(refineSpacing.first)
            {
                if(dist>refineSpacing.second)
                    subdivide=true;
            }
            if(circMarker0InCell || circMarker1InCell)
            {
                scalar minSpacing = std::min(circMarker0CellSpacing,circMarker1CellSpacing);
                if(dist>minSpacing)
                    subdivide=true;
            }
        }
        
        if(subdivide)
            refineSpanCount++;
                
        maxMarkerDist = std::max(maxMarkerDist,dist);
        minMarkerCellSize = std::min(minMarkerCellSize,circMarker0CellSpacing);
    }
    */
    //t1 = std::chrono::system_clock::now();
    //Info<<" IDC2:"<<std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
    /*
     * -------------------------------------------------------------------------- 
     */
    
    /*
     * -------------------------------- Reset ----------------------------------- 
     */
    //t0 = std::chrono::system_clock::now();
    /*
    if(refineSpanCount>0.1*circumMarkers.size())
    {
        //Info<<" Is-in";
        scalar totalCircumDist = maxMarkerDist*circumMarkers.size();
        label reqMarkerNbr = (totalCircumDist/minMarkerCellSize)+3;
        scalar reqMarkerAngleSpan = 2*Foam::constant::mathematical::pi / reqMarkerNbr;
        circumMarkers.clear();
        for(scalar angle=0; angle<=2*Foam::constant::mathematical::pi; angle+=reqMarkerAngleSpan)
        {
            LagrangianMarkerOnCrossSec newMarker
            (
                *this,mesh,rodNumber,oneRod,parameter,oneCrossSec,angle,radFrac
            );
            circumMarkers.push_back(newMarker);
        }
    }
    */
    //t1 = std::chrono::system_clock::now();
    //Info<<" Reset:"<<std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
    /*
     * -------------------------------------------------------------------------- 
     */
    
    /*
     * ------------------------------- Refine ----------------------------------- 
     */
    //t0 = std::chrono::system_clock::now();
    bool refined=true;
    while(refined)
    {
        refined = false;
        auto circMarkerIter0 = circumMarkers.begin();
        auto circMarkerIter1 = ++(circumMarkers.begin());
        for( ; circMarkerIter1!=circumMarkers.end() ; )
        {
            scalar circMarker0Angle = circMarkerIter0->getMarkerAngle();
            scalar circMarker0CellSpacing = circMarkerIter0->getMarkerCellMinSpacing();
            scalar circMarkers0CharacSpacing = circMarkerIter0->getMarkerCharacLen();
            bool circMarker0InCell = (circMarkerIter0->getMarkerCell()!=-1);
            
            scalar circMarker1Angle = circMarkerIter1->getMarkerAngle();
            scalar circMarker1CellSpacing = circMarkerIter1->getMarkerCellMinSpacing();
            scalar circMarkers1CharacSpacing = circMarkerIter1->getMarkerCharacLen();
            bool circMarker1InCell = (circMarkerIter1->getMarkerCell()!=-1);
                        
            bool subdivide = false;
            scalar dist = lowerBound_distance
            (
                rodNumber,parameter,circMarker0Angle,circMarker1Angle,radFrac
            );
            if(refineSpacing.first)
            {
                if(dist>refineSpacing.second)
                    subdivide=true;
            }
            if(circMarker0InCell || circMarker1InCell)
            {
                scalar minSpacing = std::min(circMarker0CellSpacing,circMarker1CellSpacing);
                if(dist > minSpacing*refnRodMarkersDistToMeshSpacing)
                    subdivide=true;
            }
            if(!subdivide)
            {
                dist = distance
                (
                    rodNumber,parameter,circMarker0Angle,circMarker1Angle,radFrac
                );
                if(refineSpacing.first)
                {
                    if(dist>refineSpacing.second)
                        subdivide=true;
                }
                if(circMarker0InCell || circMarker1InCell)
                {
                    scalar minSpacing = std::min(circMarker0CellSpacing,circMarker1CellSpacing);
                    if(useMarkerCharLenSpacing)
                    {
                        if(circMarker0InCell)
                            minSpacing = std::min(minSpacing,circMarkers0CharacSpacing*rodPntDistToMarkerCharLen);
                        if(circMarker1InCell)
                            minSpacing = std::min(minSpacing,circMarkers1CharacSpacing*rodPntDistToMarkerCharLen);
                    }
                    if(dist > minSpacing*refnRodMarkersDistToMeshSpacing)
                        subdivide=true;
                }
            }
            
            if(subdivide)
            {
                scalar middleAngle = 0.5*(circMarker0Angle + circMarker1Angle);
                LagrangianMarkerOnCrossSec middleAngleMarker
                (
                    *this,mesh,rodNumber,oneRod,parameter,oneCrossSec,middleAngle,radFrac
                );
                circumMarkers.insert(circMarkerIter1,middleAngleMarker);
                refined=true;
            }
            circMarkerIter0 = circMarkerIter1;
            circMarkerIter1++;
        }
    }
    //t1 = std::chrono::system_clock::now();
    //Info<<" Refine:"<<std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
    /*
     * -------------------------------------------------------------------------- 
     */
    
    /*
     * ----------------------Final data collection ---------------------------- 
     */
    /*
    t0 = std::chrono::system_clock::now();
    scalar end_avgDeltaAngle=0;
    scalar end_avgDeltaDist=0;
    deltaAngles.clear();
    deltaDist.clear();
    for(auto iter=circumMarkers.begin(); iter!=circumMarkers.end(); iter++)
    {
        std::list<LagrangianMarkerOnCrossSec>::iterator next = iter;
        next++;
        if(next!=circumMarkers.end())
        {
            scalar deltaAngle = next->getMarkerAngle()-iter->getMarkerAngle();
            end_avgDeltaAngle += deltaAngle;
            deltaAngles.push_back(deltaAngle);
            scalar dist = distance
            (
                oneRod,parameter,oneCrossSec,iter->getMarkerAngle(),next->getMarkerAngle(),radFrac
            );
            end_avgDeltaDist+=dist;
            deltaDist.push_back(dist);
        }
        else
            break;
    }
    end_avgDeltaAngle/=circumMarkers.size();
    std::sort(deltaAngles.begin(),deltaAngles.end());
    scalar end_medDeltaAngle = deltaAngles[deltaAngles.size()/2];
    end_avgDeltaDist/=circumMarkers.size();
    std::sort(deltaDist.begin(),deltaDist.end());
    scalar end_medDeltaDist = deltaDist[deltaDist.size()/2];
    label end_size = circumMarkers.size();
    */
    //t1 = std::chrono::system_clock::now();
    //auto end = t1;
    //Info<<" FDC:"<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
    
    //Info<<" refineCircumferential:"<<std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()<<Foam::nl;
    /*
     * -------------------------------------------------------------------------- 
     */
    
    /*
    Info<<"Circum ("<<start_size<<","<<start_avgDeltaAngle<<","<<start_medDeltaAngle<<","<<start_avgDeltaDist<<","<<start_medDeltaDist<<") -> ("<<end_size<<","<<end_avgDeltaAngle<<","<<end_medDeltaAngle<<","<<end_avgDeltaDist<<","<<end_medDeltaDist<<")"<<Foam::nl;
    */
}

void Foam::CrossSectionStructure::refineRadial
(
    std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>& radialMarkers,
    scalar initialSpacing,
    bool useMarkerCharLenSpacing,
    std::pair<bool,scalar> refineSpacing
)
{
    //auto t0 = std::chrono::system_clock::now();
    //auto start = t0;
    //Info<<"refineRadial"<<Foam::nl;
    
    if(radialMarkers.size()<2)
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    if(radialMarkers.front().second.size()<2)
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    
    label rodNumber = radialMarkers.front().second.front().getRodNumber();
    const ActiveRodMesh::rodCosserat* oneRod = radialMarkers.front().second.front().getBaseRod();
    scalar parameter = radialMarkers.front().second.front().getMarkerParameter();
    CrossSection* oneCrossSec =  &(rodCrossSection[rodNumber]);

    for(auto iter=radialMarkers.begin(); iter!=radialMarkers.end(); iter++)
    {
        if(iter->second.size()<1)
            FatalErrorInFunction<<"Must be at least one marker in list"<< exit(FatalError);
        if(iter->second.front().getMarkerRadiusFrac()!=0)
            refineCircumferential(iter->second,useMarkerCharLenSpacing,refineSpacing);
    }
    //auto t1 = std::chrono::system_clock::now();
    //Info<<"\t IterRadialRefn:"<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();

    //t0 = std::chrono::system_clock::now();
    bool refined=true;
    while(refined)
    {
        //Info<<"\tRefine while"<<Foam::nl;
        refined = false;
        auto radMarkerIter0 = radialMarkers.begin();
        auto radMarkerIter1 = ++(radialMarkers.begin());
        for( ; radMarkerIter1!=radialMarkers.end() ; )
        {
            if(radMarkerIter0->second.size()<1)
            {
                Info<<radMarkerIter0->second.front().getMarkerParameter()<<"/"<<radMarkerIter0->second.front().getMarkerRadiusFrac()<<Foam::nl;
                FatalErrorInFunction<<"Must be at least one markers"<< exit(FatalError);
            }
            scalar radMarker0RadFrac = radMarkerIter0->second.front().getMarkerRadiusFrac();
            scalar radMarker0CellSpacing = std::numeric_limits<scalar>::max();
            scalar radMarker0CharacSpacing = std::numeric_limits<scalar>::max();
            bool radMarker0InCell = false;
            for(auto iter=radMarkerIter0->second.begin(); iter!=radMarkerIter0->second.end(); iter++)
            {
                radMarker0CellSpacing = std::min(radMarker0CellSpacing,iter->getMarkerCellMinSpacing());
                radMarker0CharacSpacing = std::min(radMarker0CharacSpacing,iter->getMarkerCharacLen());
                radMarker0InCell |= (iter->getMarkerCell()!=-1);
            }
            
            if(radMarkerIter1->second.size()<1)
            {
                Info<<radMarkerIter1->second.front().getMarkerParameter()<<"/"<<radMarkerIter1->second.front().getMarkerRadiusFrac()<<Foam::nl;
                FatalErrorInFunction<<"Must be at least one markers"<< exit(FatalError);
            }
            scalar radMarker1RadFrac = radMarkerIter1->second.front().getMarkerRadiusFrac();
            scalar radMarker1CellSpacing = std::numeric_limits<scalar>::max();
            scalar radMarker1CharacSpacing = std::numeric_limits<scalar>::max();
            bool radMarker1InCell = false;
            for(auto iter=radMarkerIter1->second.begin(); iter!=radMarkerIter1->second.end(); iter++)
            {
                radMarker1CellSpacing = std::min(radMarker1CellSpacing,iter->getMarkerCellMinSpacing());
                radMarker1CharacSpacing = std::min(radMarker1CharacSpacing,iter->getMarkerCharacLen());
                radMarker1InCell |= (iter->getMarkerCell()!=-1);
            }
            
            //Info<<"radMarker0RadFrac:"<<radMarker0RadFrac<<"  radMarker1RadFrac:"<<radMarker1RadFrac<<Foam::nl;
            //Info<<"radMarker0CellSpacing:"<<radMarker0CellSpacing<<"  radMarker1CellSpacing:"<<radMarker1CellSpacing<<Foam::nl;
            //Info<<"radMarker0InCell:"<<radMarker0InCell<<"  radMarker1InCell:"<<radMarker1InCell<<Foam::nl;
            
            scalar maxDist = std::numeric_limits<scalar>::max();
            for(auto iter=radMarkerIter0->second.begin(); iter!=radMarkerIter0->second.end(); iter++)
            {
                scalar angle = iter->getMarkerAngle();
                vector radMarker0Pos = evaluateRodCircumPos
                (
                    oneRod,parameter,oneCrossSec,angle,radMarker0RadFrac
                );
                vector radMarker1Pos = evaluateRodCircumPos
                (
                    oneRod,parameter,oneCrossSec,angle,radMarker1RadFrac
                );
                vector distVec = radMarker0Pos-radMarker1Pos;
                maxDist = std::max(maxDist,std::sqrt(distVec&distVec));
            }
            for(auto iter=radMarkerIter1->second.begin(); iter!=radMarkerIter1->second.end(); iter++)
            {
                scalar angle = iter->getMarkerAngle();
                vector radMarker0Pos = evaluateRodCircumPos
                (
                    oneRod,parameter,oneCrossSec,angle,radMarker0RadFrac
                );
                vector radMarker1Pos = evaluateRodCircumPos
                (
                    oneRod,parameter,oneCrossSec,angle,radMarker1RadFrac
                );
                vector distVec = radMarker0Pos-radMarker1Pos;
                maxDist = std::max(maxDist,std::sqrt(distVec&distVec));
            }

            bool subdivide = false;
            if(refineSpacing.first)
            {
                if(maxDist>refineSpacing.second)
                    subdivide=true;
            }
            if(radMarker0InCell || radMarker1InCell)
            {
                scalar minSpacing = std::min(radMarker0CellSpacing,radMarker1CellSpacing);
                if(useMarkerCharLenSpacing)
                {
                    if(radMarker0InCell)
                        minSpacing = std::min(minSpacing,radMarker0CharacSpacing*rodPntDistToMarkerCharLen);
                    if(radMarker1InCell)
                        minSpacing = std::min(minSpacing,radMarker1CharacSpacing*rodPntDistToMarkerCharLen);
                }
                if(maxDist > minSpacing*refnRodMarkersDistToMeshSpacing)
                    subdivide=true;
            }

            if(subdivide)
            {
                scalar middleRadiusFrac = 0.5*(radMarker0RadFrac + radMarker1RadFrac);
                //Info<<"\tSubdivide:"<<"middleRadiusFrac:"<<middleRadiusFrac<<" radMarker0RadFrac:"<<radMarker0RadFrac<<" radMarker1RadFrac:"<<radMarker1RadFrac<<Foam::nl;
                std::vector<scalar> angleData;
                createSpacedPointsOnCrossSec
                (
                    oneRod,parameter,oneCrossSec,middleRadiusFrac,initialSpacing,angleData
                );
                std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>> radialSlice;
                radialSlice.first = middleRadiusFrac;
                std::list<LagrangianMarkerOnCrossSec>& angleMarkers = radialSlice.second;
                for(scalar angle : angleData)
                {
                    LagrangianMarkerOnCrossSec middleRadiusMarker
                    (
                        *this,mesh,rodNumber,oneRod,parameter,oneCrossSec,angle,middleRadiusFrac
                    );
                    angleMarkers.push_back(middleRadiusMarker);
                }
                refineCircumferential(angleMarkers,useMarkerCharLenSpacing,refineSpacing);
                radialMarkers.insert(radMarkerIter1,radialSlice);
                refined=true;
            }
            radMarkerIter0 = radMarkerIter1;
            radMarkerIter1++;
        }
    }    
    //t1 = std::chrono::system_clock::now();
    //Info<<" RadialRefn:"<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
    //auto end = t1;
    //Info<<"  refineRadial:"<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()<<Foam::nl;

}

void Foam::CrossSectionStructure::refineTangential
(
    std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>& tangMarkers,
    scalar initialSpacing,
    bool useMarkerCharLenSpacing,
    std::pair<bool,scalar> refineSpacing
)
{
    //auto t0 = std::chrono::system_clock::now();
    //auto start = t0;
    
    if(tangMarkers.size()<2)
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    if(tangMarkers.front().second.size()<2)
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    if(tangMarkers.front().second.front().second.size()<2)
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    
    label rodNumber = tangMarkers.front().second.front().second.front().getRodNumber();
    const ActiveRodMesh::rodCosserat* oneRod = tangMarkers.front().second.front().second.front().getBaseRod();
    CrossSection* oneCrossSec =  &(rodCrossSection[rodNumber]);

    bool refined=true;
    while(refined)
    {
        refined = false;
        auto tangMarkerIter0 = tangMarkers.begin();
        auto tangMarkerIter1 = ++(tangMarkers.begin());
        for( ; tangMarkerIter1!=tangMarkers.end() ; )
        {
            
            //auto ti0 = std::chrono::system_clock::now();

            if(tangMarkerIter0->second.size()<1)
                FatalErrorInFunction<<"Must be at least one radial marker layer"<< exit(FatalError);
            if(tangMarkerIter0->second.front().second.size()<1)
                FatalErrorInFunction<<"Must be at least one angle marker layer"<< exit(FatalError);
            scalar tangMarker0Para = tangMarkerIter0->second.front().second.front().getMarkerParameter();
            scalar tangMarker0CellSpacing = std::numeric_limits<scalar>::max();
            scalar tangMarker0CharacSpacing = std::numeric_limits<scalar>::max();
            bool tangMarker0InCell = false;
            for(auto iter=tangMarkerIter0->second.front().second.begin(); iter!=tangMarkerIter0->second.front().second.end(); iter++)
            {
                tangMarker0CellSpacing = std::min(tangMarker0CellSpacing,iter->getMarkerCellMinSpacing());
                tangMarker0CharacSpacing = std::min(tangMarker0CharacSpacing,iter->getMarkerCharacLen());
                tangMarker0InCell |= (iter->getMarkerCell()!=-1);
            }
            
            if(tangMarkerIter1->second.size()<1)
                FatalErrorInFunction<<"Must be at least one radial marker layer"<< exit(FatalError);
            if(tangMarkerIter1->second.front().second.size()<1)
                FatalErrorInFunction<<"Must be at least one angle marker layer"<< exit(FatalError);
            scalar tangMarker1Para = tangMarkerIter1->second.front().second.front().getMarkerParameter();
            scalar tangMarker1CellSpacing = std::numeric_limits<scalar>::max();
            scalar tangMarker1CharacSpacing = std::numeric_limits<scalar>::max();
            bool tangMarker1InCell = false;
            for(auto iter=tangMarkerIter1->second.front().second.begin(); iter!=tangMarkerIter1->second.front().second.end(); iter++)
            {
                tangMarker1CellSpacing = std::min(tangMarker1CellSpacing,iter->getMarkerCellMinSpacing());
                tangMarker1CharacSpacing = std::min(tangMarker1CharacSpacing,iter->getMarkerCharacLen());
                tangMarker1InCell |= (iter->getMarkerCell()!=-1);
            }
            
            //auto ti1 = std::chrono::system_clock::now();
            //Info<<" CSC:"<<std::chrono::duration_cast<std::chrono::milliseconds>(ti1-ti0).count();

            
            vector r,d2,d3;        
            vector pnt0D1;
            rodEval(oneRod,tangMarker0Para,pnt0D1,d2,d3,r);
            vector pnt1D1;
            rodEval(oneRod,tangMarker1Para,pnt1D1,d2,d3,r);
            
            if((pnt0D1&pnt1D1) < 0)
                FatalErrorInFunction<<"Vector can not be pointing in different directions"<< exit(FatalError);
            
            scalar angle = std::acos((pnt0D1&pnt1D1)/(std::sqrt(pnt0D1&pnt0D1)*std::sqrt(pnt1D1&pnt1D1)));
            scalar halfAngle = std::abs(angle/2);
            scalar rodDistance = LineStructure::distance(oneRod,tangMarker0Para,tangMarker1Para);
            
            scalar maxRadiusPnt0 = oneCrossSec->upperLimitRadius(tangMarker0Para);
            scalar addedDistPnt0 = maxRadiusPnt0*std::sin(halfAngle);
            scalar maxRadiusPnt1 = oneCrossSec->upperLimitRadius(tangMarker1Para);
            scalar addedDistPnt1 = maxRadiusPnt1*std::sin(halfAngle);
            
            scalar totalDistance = rodDistance+addedDistPnt0+addedDistPnt1;

            bool subdivide = false;
            if(refineSpacing.first)
            {
                if(totalDistance>refineSpacing.second)
                    subdivide=true;
            }
            if(tangMarker0InCell || tangMarker1InCell)
            {
                scalar minSpacing = std::min(tangMarker0CellSpacing,tangMarker1CellSpacing);
                if(useMarkerCharLenSpacing)
                {
                    if(tangMarker0InCell)
                        minSpacing = std::min(minSpacing,tangMarker0CharacSpacing*rodPntDistToMarkerCharLen);
                    if(tangMarker1InCell)
                        minSpacing = std::min(minSpacing,tangMarker1CharacSpacing*rodPntDistToMarkerCharLen);
                }
                if(totalDistance > minSpacing*refnRodMarkersDistToMeshSpacing)
                    subdivide=true;
            }
            
            //auto ti2 = std::chrono::system_clock::now();
            //Info<<" Dec:"<<std::chrono::duration_cast<std::chrono::milliseconds>(ti2-ti1).count();

            if(subdivide)
            {
                scalar middleMarkerPara = 0.5*(tangMarker0Para + tangMarker1Para);
                std::vector<scalar> angleData;
                createSpacedPointsOnCrossSec
                (
                    oneRod,middleMarkerPara,oneCrossSec,1.0,initialSpacing,angleData
                );

                //auto ti3 = std::chrono::system_clock::now();
                //Info<<" CLP:"<<std::chrono::duration_cast<std::chrono::milliseconds>(ti3-ti2).count();

                std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>> angleMarkers;
                angleMarkers.first = 1.0;
                for(scalar angle : angleData)
                {
                    angleMarkers.second.push_back(
                        LagrangianMarkerOnCrossSec
                        (
                            *this,mesh,rodNumber,oneRod,middleMarkerPara,oneCrossSec,angle,1.0
                        ));
                }
                
                //auto ti4 = std::chrono::system_clock::now();
                //Info<<" ToMarkers:"<<std::chrono::duration_cast<std::chrono::milliseconds>(ti4-ti3).count();
                
                refineCircumferential(angleMarkers.second,useMarkerCharLenSpacing,refineSpacing);
                
                //auto ti5 = std::chrono::system_clock::now();
                //Info<<" RefnCirc:"<<std::chrono::duration_cast<std::chrono::milliseconds>(ti5-ti4).count();
                
                std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>> paraSlice;
                paraSlice.first = middleMarkerPara;
                std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>& radialMarkers = paraSlice.second;
                radialMarkers.push_back(angleMarkers);
                tangMarkers.insert(tangMarkerIter1,paraSlice);
                refined=true;
                //auto ti6 = std::chrono::system_clock::now();
                //Info<<" Ins:"<<std::chrono::duration_cast<std::chrono::milliseconds>(ti6-ti5).count();
            }
            tangMarkerIter0 = tangMarkerIter1;
            tangMarkerIter1++;
            
            //auto ti7 = std::chrono::system_clock::now();
            //Info<<" Div:"<<std::chrono::duration_cast<std::chrono::milliseconds>(ti7-ti2).count()<<Foam::nl;
        }
    }
    
    //auto t1 = std::chrono::system_clock::now();
    //auto end = t1;
    //Info<<"  refineTangential:"<<std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count()<<Foam::nl;
}

void Foam::CrossSectionStructure::setMarkerVolumeOnRod
(
    label rodNumber
)
{
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    CrossSection* oneCrossSec = &(rodCrossSection[rodNumber]);
    
    std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>& markersList = *(rodMarkersList[rodNumber]);
    if(rodMarkersList[rodNumber])
    {
        std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec*>>> markers(markersList.size());
        uint P=0;
        for(auto iterPara=markersList.begin(); iterPara!=markersList.end(); iterPara++,P++)
        {
            std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>& radialMarkersList = *iterPara;
            std::vector<std::vector<LagrangianMarkerOnCrossSec*>>& radialMarkers = markers[P];
            radialMarkers.resize(radialMarkersList.second.size());
            uint R=0;
            for(auto iterRadial=radialMarkersList.second.begin(); iterRadial!=radialMarkersList.second.end(); iterRadial++,R++)
            {
                std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>& angleMarkersList = *iterRadial;
                std::vector<LagrangianMarkerOnCrossSec*>& angleMarkers = radialMarkers[R];
                //angleMarkers.resize(radialMarkersList.size());
                for(auto iterAngle=angleMarkersList.second.begin(); iterAngle!=angleMarkersList.second.end(); iterAngle++)
                {
                    angleMarkers.push_back(&(*iterAngle));
                }
            }
        }
        
        std::function<vector(scalar,scalar,scalar,scalar,scalar)> getPosition = 
        [rod=oneRod,crossSec=oneCrossSec]
        (scalar para, scalar radFrac, scalar angle, scalar paraDev, scalar radiusDev)
        {
            return CrossSectionStructure::evaluateRodCircumPos
            (
                rod,para,crossSec,angle,radFrac,paraDev,radiusDev
            );
        };
        
        scalar sumVolume = 0;
        const cellList& cells = mesh.cells();
        const faceList& faces = mesh.faces();
        const pointField& points = mesh.points();
        enum Surface {Circumferential, Radial};
        Surface surfaceType;
        for(uint rodWiseIndex=0; rodWiseIndex<markers.size(); rodWiseIndex++)
        {
            if(rodWiseIndex==0 || rodWiseIndex==markers.size()-1)
                surfaceType = Surface::Radial;
            else
                surfaceType = Surface::Circumferential;
            
            //Info<<"rodWiseIndex:"<<rodWiseIndex<<Foam::nl;
            std::vector<std::vector<LagrangianMarkerOnCrossSec*>>& radialMarkers = markers[rodWiseIndex];
            scalar rodWiseIndexParameter = radialMarkers[0][0]->getMarkerParameter();
            
            scalar prevParameter;
            if(rodWiseIndex==0)
                prevParameter = rodWiseIndexParameter;
            else
                prevParameter = (markers[rodWiseIndex-1][0][0]->getMarkerParameter()+rodWiseIndexParameter)/2;
            
            scalar subseqParameter;
            if(rodWiseIndex==markers.size()-1)
                subseqParameter = rodWiseIndexParameter;
            else
                subseqParameter = (markers[rodWiseIndex+1][0][0]->getMarkerParameter()+rodWiseIndexParameter)/2;
            
            for(uint radialIndex=0; radialIndex<radialMarkers.size(); radialIndex++)
            {
                std::vector<LagrangianMarkerOnCrossSec*>& circumMarkers = radialMarkers[radialIndex];
                scalar radialIndexFrac = circumMarkers[0]->getMarkerRadiusFrac();
                
                scalar outerRadiusFrac, innerRadiusFrac;
                if(radialMarkers.size()==1)
                {
                    outerRadiusFrac = 1;
                    innerRadiusFrac = 0;
                }
                else 
                {
                    if(radialIndex==0)
                        outerRadiusFrac = radialIndexFrac;
                    else
                        outerRadiusFrac = (radialMarkers[radialIndex-1][0]->getMarkerRadiusFrac()+radialIndexFrac)/2;
                    
                    if(radialIndex==radialMarkers.size()-1)
                        innerRadiusFrac = radialIndexFrac;
                    else
                        innerRadiusFrac = (radialMarkers[radialIndex+1][0]->getMarkerRadiusFrac()+radialIndexFrac)/2;
                }
                //Info<<"\tradialIndex:"<<radialIndex<<" - "<<circumMarkers.size()<<Foam::nl;
                for(uint circIndex=0; circIndex<circumMarkers.size(); circIndex++)
                {
                    //Info<<"\t\tcircIndex:"<<circIndex<<"/"<<circumMarkers.size()<<Foam::nl;

                    LagrangianMarkerOnCrossSec* singleMarker = circumMarkers[circIndex];
                    if(singleMarker->getMarkerRadiusFrac()!=1.0 && surfaceType==Surface::Circumferential)
                        FatalErrorInFunction<<"Only end of rod markers can have a radiusFraction != 1"<<exit(FatalError);

                    scalar circIndexAngle = singleMarker->getMarkerAngle();

                    scalar lowerAngle, upperAngle;
                    if(circIndex==0)
                        lowerAngle = circumMarkers.back()->getMarkerAngle()-2*Foam::constant::mathematical::pi;
                    else
                        lowerAngle = circumMarkers[circIndex-1]->getMarkerAngle();
                    lowerAngle += circIndexAngle;
                    lowerAngle /= 2;
                    lowerAngle = CrossSectionStructure::restrictAngle(lowerAngle);
                    
                    if(circIndex==circumMarkers.size()-1)
                        upperAngle = circumMarkers[0]->getMarkerAngle()+2*Foam::constant::mathematical::pi;
                    else
                        upperAngle = circumMarkers[circIndex+1]->getMarkerAngle();
                    
                    upperAngle += circIndexAngle;
                    upperAngle /= 2;
                    upperAngle = CrossSectionStructure::restrictAngle(upperAngle);
                    
                    if(circumMarkers.size()==1)
                    {
                        lowerAngle = 0;
                        upperAngle = Foam::constant::mathematical::pi;
                    }
                    
                    if(singleMarker->getMarkerCell()!=-1)
                    {
                        if(modusMarkerToField == markerMeshType::Uniform)
                        {
                            scalar h = std::cbrt(cells[singleMarker->getMarkerCell()].mag(points,faces));
                            
                            /*
                            Info<<"prevParameter:"<<prevParameter<<"  subseqParameter:"<<subseqParameter<<Foam::nl;
                            Info<<"innerRadiusFrac:"<<innerRadiusFrac<<"  outerRadiusFrac:"<<outerRadiusFrac<<Foam::nl;
                            Info<<"lowerAngle:"<<lowerAngle<<"  upperAngle:"<<upperAngle<<Foam::nl;
                            */
                            
                            std::array<std::array<std::array<label,2>,2>,2> paraRadAnglePnts;
                            pointField points(8);
                            paraRadAnglePnts[0][0][0] = 0;
                            paraRadAnglePnts[0][0][1] = 1;
                            paraRadAnglePnts[0][1][0] = 2;
                            paraRadAnglePnts[0][1][1] = 3;
                            paraRadAnglePnts[1][0][0] = 4;
                            paraRadAnglePnts[1][0][1] = 5;
                            paraRadAnglePnts[1][1][0] = 6;
                            paraRadAnglePnts[1][1][1] = 7;
                            if(surfaceType==Surface::Circumferential)
                            {
                                points[0] = getPosition(prevParameter,radialIndexFrac,lowerAngle,0,-h/2);
                                points[1] = getPosition(prevParameter,radialIndexFrac,upperAngle,0,-h/2);
                                points[2] = getPosition(prevParameter,radialIndexFrac,lowerAngle,0,h/2);
                                points[3] = getPosition(prevParameter,radialIndexFrac,upperAngle,0,h/2);
                                points[4] = getPosition(subseqParameter,radialIndexFrac,lowerAngle,0,-h/2);
                                points[5] = getPosition(subseqParameter,radialIndexFrac,upperAngle,0,-h/2);
                                points[6] = getPosition(subseqParameter,radialIndexFrac,lowerAngle,0,h/2);
                                points[7] = getPosition(subseqParameter,radialIndexFrac,upperAngle,0,h/2);
                            }
                            else
                            {
                                points[0] = getPosition(rodWiseIndexParameter,innerRadiusFrac,lowerAngle,-h/2,0);
                                points[1] = getPosition(rodWiseIndexParameter,innerRadiusFrac,upperAngle,-h/2,0);
                                points[2] = getPosition(rodWiseIndexParameter,outerRadiusFrac,lowerAngle,-h/2,0);
                                points[3] = getPosition(rodWiseIndexParameter,outerRadiusFrac,upperAngle,-h/2,0);
                                points[4] = getPosition(rodWiseIndexParameter,innerRadiusFrac,lowerAngle,h/2,0);
                                points[5] = getPosition(rodWiseIndexParameter,innerRadiusFrac,upperAngle,h/2,0);
                                points[6] = getPosition(rodWiseIndexParameter,outerRadiusFrac,lowerAngle,h/2,0);
                                points[7] = getPosition(rodWiseIndexParameter,outerRadiusFrac,upperAngle,h/2,0);
                            }
                    
                            face prevFace({0,1,3,2});
                            face subseqFace({4,5,7,6});
                            face innerFace({0,1,5,4});
                            face outerFace({2,3,7,6});
                            face lowerFace({0,2,6,4});
                            face upperFace({1,3,7,5});
                            
                            List<face> faces(12);
                            faces[0] = {prevFace[0],prevFace[1],prevFace[3]};
                            faces[1] = {prevFace[1],prevFace[2],prevFace[3]};
                            
                            faces[2] = {subseqFace[0],subseqFace[1],subseqFace[3]};
                            faces[3] = {subseqFace[1],subseqFace[2],subseqFace[3]};
                            
                            faces[4] = {innerFace[0],innerFace[1],innerFace[3]};
                            faces[5] = {innerFace[1],innerFace[2],innerFace[3]};
                            
                            faces[6] = {outerFace[0],outerFace[1],outerFace[3]};
                            faces[7] = {outerFace[1],outerFace[2],outerFace[3]};
                            
                            faces[8] = {lowerFace[0],lowerFace[1],lowerFace[3]};
                            faces[9] = {lowerFace[1],lowerFace[2],lowerFace[3]};
                            
                            faces[10] = {upperFace[0],upperFace[1],upperFace[3]};
                            faces[11] = {upperFace[1],upperFace[2],upperFace[3]};
                            
                            List<label> cellList(12);
                            for(label ind=0; ind<cellList.size(); ind++)
                                cellList[ind] = ind;
                            cell thisCell(cellList);
                            
                            scalar markerVolume = thisCell.mag(points,faces);
                            singleMarker->setMarkerVolume(markerVolume);
                            //Info<<"|||"<<singleMarker.getMarkerVolume()<<"|||"<<singleMarker.getSupportCells().size()<<Foam::nl;
                            sumVolume+=singleMarker->getMarkerVolume();
                        }
                        else
                        {
                            pointField points(4);
                            if(surfaceType==Surface::Circumferential)
                            {                  
                                points[0] = getPosition(prevParameter,radialIndexFrac,lowerAngle,0,0);
                                points[1] = getPosition(prevParameter,radialIndexFrac,upperAngle,0,0);
                                points[2] = getPosition(subseqParameter,radialIndexFrac,upperAngle,0,0);
                                points[3] = getPosition(subseqParameter,radialIndexFrac,lowerAngle,0,0);
                            }
                            else
                            {
                                points[0] = getPosition(rodWiseIndexParameter,innerRadiusFrac,lowerAngle,0,0);
                                points[1] = getPosition(rodWiseIndexParameter,innerRadiusFrac,upperAngle,0,0);
                                points[2] = getPosition(rodWiseIndexParameter,outerRadiusFrac,upperAngle,0,0);
                                points[3] = getPosition(rodWiseIndexParameter,outerRadiusFrac,lowerAngle,0,0);
                            }
                            face thisFace({0,1,2,3});
                            
                            scalar markerVolume = thisFace.mag(points);
                            
                            /*
                            if(markerVolume<1e-10)
                            {
                                Info<<"markerVolume:"<<markerVolume<<"  "<<thisFace<<Foam::nl;
                                Info<<"     points:"<<points<<Foam::nl;
                                Info<<"     prevParameter:"<<prevParameter<<Foam::nl;
                                Info<<"     subseqParameter:"<<subseqParameter<<Foam::nl;
                                Info<<"     innerRadiusFrac:"<<innerRadiusFrac<<Foam::nl;
                                Info<<"     radialIndexFrac:"<<radialIndexFrac<<Foam::nl;
                                Info<<"     outerRadiusFrac:"<<outerRadiusFrac<<Foam::nl;
                                Info<<"     lowerAngle:"<<lowerAngle<<Foam::nl;
                                Info<<"     upperAngle:"<<upperAngle<<Foam::nl;       
                                Info<<"     distAngle:"<<"("<<lowerAngle<<"/"<<circIndexAngle<<"/"<<upperAngle<<"):"<<upperAngle-lowerAngle<<Foam::nl;
                                Info<<"     circumMarkers.size():"<<circumMarkers.size()<<Foam::nl;
                                
                                Info<<"     p1:"<<getPosition(prevParameter,radialIndexFrac,lowerAngle+0.2,0,0)<<Foam::nl;
                                Info<<"     p2:"<<getPosition(prevParameter,radialIndexFrac,lowerAngle+0.4,0,0)<<Foam::nl;
                                Info<<"     p3:"<<getPosition(prevParameter,radialIndexFrac,lowerAngle+0.6,0,0)<<Foam::nl;
                                Info<<"     p4:"<<getPosition(prevParameter,radialIndexFrac,lowerAngle+0.8,0,0)<<Foam::nl;
                                Info<<"     singleMarker:"<<singleMarker->to_string()<<Foam::nl;
                                Info<<"     singleMarkerPara:"<<singleMarker->getMarkerParameter()<<Foam::nl;
                                Info<<"     singleMarkerRadiusFrac:"<<singleMarker->getMarkerRadiusFrac()<<Foam::nl;
                                Info<<"     singleMarkerAnlge:"<<singleMarker->getMarkerAngle()<<Foam::nl;
                            }
                            */
                            
                            singleMarker->setMarkerVolume(markerVolume);
                            //Info<<"|||"<<singleMarker.getMarkerVolume()<<"|||"<<singleMarker.getSupportCells().size()<<Foam::nl;
                            sumVolume+=singleMarker->getMarkerVolume();
                        }
                    }
                }
            }
        }
    }
}

void Foam::CrossSectionStructure::evaluateMarkerMeshRelation()
{
    status.execValid(status.markerMesh);
    for(std::unique_ptr<std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>>& singleRodMarkers :  rodMarkersList)
    {
        if(singleRodMarkers)
        {
            for(std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>& oneParaRodMarkers : *singleRodMarkers)
            {
                for(std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>& radialFracRodMarkers : oneParaRodMarkers.second)
                {
                    evaluateMarkerMeshRelation(radialFracRodMarkers.second);
                }
            }
        }
    }
    status.executed(status.markerMesh);
}

void Foam::CrossSectionStructure::evaluateMarkerMeshRelation
(
    std::list<LagrangianMarkerOnCrossSec>& markerList
)
{
    for(LagrangianMarkerOnCrossSec& marker : markerList)
        marker.evaluateMarker();
}

void Foam::CrossSectionStructure::reduceMarkers()
{   
    std::vector<MarkerReference<LagrangianMarkerOnCrossSec>> allMarkers;
    for(std::unique_ptr<std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>>& singleRodMarkers :  rodMarkersList)
    {
        if(singleRodMarkers)
        {
            for(std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>& oneParaRodMarkers : *singleRodMarkers)
            {
                for(std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>& radialFracRodMarkers : oneParaRodMarkers.second)
                {
                    std::list<LagrangianMarkerOnCrossSec>* singleRadMarkers = &(radialFracRodMarkers.second);
                    for(auto iter=singleRadMarkers->begin(); iter!=singleRadMarkers->end(); iter++)
                    {
                        allMarkers.push_back
                        (
                            MarkerReference<LagrangianMarkerOnCrossSec>(iter,singleRadMarkers)
                        );
                    }
                }
            }
        }
    }
    LineStructure::reduceMarkers(allMarkers);
}

Foam::BoundingBox Foam::CrossSectionStructure::computeBox
(
    label rodNumber
) const
{
    BoundingBox nurbsCurveBox = LineStructure::computeBox(rodNumber);
    std::pair<scalar,scalar> minMax = rodCrossSection[rodNumber].radiusBounds();
    nurbsCurveBox.enlarge(minMax.second);
    return nurbsCurveBox;
}

Foam::BoundingBox Foam::CrossSectionStructure::computeBox
(
    label rodNumber,
    scalar parStart,
    scalar parEnd
) const
{
    BoundingBox nurbsCurveBox = LineStructure::computeBox(rodNumber,parStart,parEnd);
    std::pair<scalar,scalar> minMax = rodCrossSection[rodNumber].radiusBounds(parStart,parEnd);
    nurbsCurveBox.enlarge(minMax.second);
    return nurbsCurveBox;
}

Foam::scalar Foam::CrossSectionStructure::characteristicSize
(
    label rodNumber,
    scalar par
) const
{
    FatalErrorInFunction<<"Not yet implemented"<<exit(FatalError);
    return rodCrossSection[rodNumber].upperLimitRadius(par);
}

void Foam::CrossSectionStructure::removeOverlapMarkers()
{
    std::vector<MarkerReference<LagrangianMarkerOnCrossSec>> deleteMarkers;
    for(std::size_t rodI=0; rodI<rodMarkersList.size(); rodI++)
    {
        if(rodInMesh[rodI])
        {
            std::unique_ptr<std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>>& singleRodMarkers = rodMarkersList[rodI];
            for(std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>& oneParaRodMarkers : *singleRodMarkers)
            {
                for(std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>& radialFracRodMarkers : oneParaRodMarkers.second)
                {
                    std::list<LagrangianMarkerOnCrossSec>* singleRadMarkers = &(radialFracRodMarkers.second);
                    for(auto iter=singleRadMarkers->begin(); iter!=singleRadMarkers->end(); iter++)
                    {
                        LagrangianMarkerOnCrossSec& marker = *iter;
                        vector P = marker.getMarkerPosition();
                        
                        for(std::size_t rodIOther=0; rodIOther<rodMarkersList.size(); rodIOther++)
                        {
                            if(rodInMesh[rodIOther] && rodI!=rodIOther)
                            {
                                const BoundingBoxTree& tree = rodTrees[rodIOther];
                                std::vector<scalar> parameters;
                                tree.findPointParameters(parameters,P);
                                
                                if(parameters.size()>0)
                                {
                                    const ActiveRodMesh::rodCosserat* rod = Rods[rodIOther];
                                    CrossSection& crossSec = rodCrossSection[rodIOther];
                                    std::function<vector(scalar)> C_u=[rod](scalar u){return rodEval(rod,u);};
                                    std::function<vector(scalar)> dC_u=[rod](scalar u){return rodDerivEval(rod,u);};
                                    std::function<vector(scalar)> d2C_u=[rod](scalar u){return rodDeriv2Eval(rod,u);};
                                    
                                    std::function<scalar(scalar)> f_u=[rod,P,dC_u,C_u](scalar u)
                                    {
                                        return dC_u(u)&(C_u(u)-P);                                        
                                    };
                                    std::function<scalar(scalar)> df_u=[rod,P,d2C_u,dC_u,C_u](scalar u)
                                    {
                                        vector dC = dC_u(u);
                                        scalar abs_dC_2 = dC&dC;
                                        return (d2C_u(u)&(C_u(u)-P)) + abs_dC_2;
                                    };
                                    
                                    std::function<scalar(scalar)> f_df_u = [rod,P,f_u,df_u](scalar u)
                                    {
                                        return u-f_u(u)/df_u(u);
                                    };
                                    
                                    scalar u_start = rod->m_Curve.domainStart();
                                    scalar u_end = rod->m_Curve.domainEnd();
                                    
                                    for(scalar& para : parameters)
                                    {
                                        for(label iteration=0; iteration<1000; iteration++)
                                        {
                                            scalar paraNext = f_df_u(para);
                                            paraNext = std::max(paraNext,u_start);
                                            paraNext = std::min(paraNext,u_end);
                                            
                                            scalar orthogonalErr = std::abs(f_u(paraNext));
                                            if(orthogonalErr<1e-8)
                                                break;
                                            
                                            scalar changeMag = std::abs((paraNext-para)*(paraNext-para)*(dC_u(paraNext)&dC_u(paraNext)));
                                            if(changeMag<1e-8)
                                                break;
                                            
                                            para = paraNext;
                                        }
                                    }
                                    scalar minimalPara = parameters[0];
                                    scalar minimalOffOrtho = std::abs(f_u(parameters[0]));
                                    for(scalar& para : parameters)
                                    {
                                        scalar offOrtho = std::abs(f_u(para));
                                        if(offOrtho<minimalOffOrtho)
                                        {
                                            minimalPara = para;
                                            minimalOffOrtho = offOrtho;
                                        }
                                    }
                                    if(minimalOffOrtho>1e-8)
                                    {
                                        if(std::abs(minimalPara-u_start)>1e-4)
                                        {
                                            if(std::abs(minimalPara-u_end)>1e-4)
                                            {
                                                FatalErrorInFunction<<"Non solution found"<<exit(FatalError);
                                            }
                                        }
                                    }
                                    if(u_start==minimalPara || u_end==minimalPara)
                                        continue;
                                        
                                    vector d1,d2,d3,r;
                                    rodEval(rod,minimalPara,d1,d2,d3,r);
                                    
                                    scalar v1 = (d1&P)/(d1&d1);
                                    scalar v2 = (d2&P)/(d2&d2);
                                    
                                    const scalar sign_v1 = v1<0?-1:1;
                                    const scalar sign_v2 = v2<0?-1:1;
                                    
                                    v1 *= sign_v1;
                                    v2 *= sign_v2;
                                    
                                    scalar angle;
                                    if(v1!=0)
                                        angle = std::atan(v2/v1);
                                    else
                                        angle = 0.5*Foam::constant::mathematical::pi;
                                    
                                    if(sign_v1==-1)
                                        angle = Foam::constant::mathematical::pi-angle;
                                    if(sign_v2==-1)
                                        angle = 2*Foam::constant::mathematical::pi-angle;
                                    
                                    vector dist = P-r;
                                    scalar radius = std::sqrt(dist&dist);
                                    
                                    if(radius<crossSec(minimalPara,angle))
                                    {
                                        deleteMarkers.push_back
                                        (
                                            MarkerReference<LagrangianMarkerOnCrossSec>(iter,singleRadMarkers)
                                        );
                                        deleteMarkers.back().deleteMarker();
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void Foam::CrossSectionStructure::collectMarkers()
{
    status.execValid(status.markersCollected);
    collectedMarkers.resize(0);
    for(std::size_t rodNumber=0; rodNumber<rodMarkersList.size(); rodNumber++)
    {
        std::unique_ptr<std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>>& singleRodMarkers = rodMarkersList[rodNumber];
        if(rodInMesh[rodNumber])
        {
            if(singleRodMarkers)
            {
                for(std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>& oneParaRodMarkers : *singleRodMarkers)
                {
                    for(std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>& radialFracRodMarkers : oneParaRodMarkers.second)
                    {
                        for(LagrangianMarkerOnCrossSec& marker : radialFracRodMarkers.second)
                        {
                            collectedMarkers.push_back(&marker);
                        }
                    }
                }
            }
            else
            {
                meshBoundingBox.print();
                Pout<<"rodInMesh["<<rodNumber<<"]:"<<rodInMesh[rodNumber]<<Foam::nl;
                Pout<<"rodTrees["<<rodNumber<<"]:"; rodTrees[rodNumber].printRoot();
                FatalErrorInFunction<<"Rod with no markers given but in mesh"<<exit(FatalError);
            }
        }
        else
        {
            if(singleRodMarkers)
                FatalErrorInFunction<<"Rod out of mesh but markers given"<<exit(FatalError);
        }
    }
    status.executed(status.markersCollected);
}

Foam::scalar Foam::CrossSectionStructure::evaluateCircumArcLen
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameterA,
    scalar parameterB,
    CrossSection* oneCrossSec,
    scalar angleA,
    scalar angleB,
    scalar radiusFracA,
    scalar radiusFracB
)
{
    vector parAVec = CrossSectionStructure::evaluateRodCircumPos(oneRod,parameterA,oneCrossSec,angleA,radiusFracA);
    vector parBVec = CrossSectionStructure::evaluateRodCircumPos(oneRod,parameterB,oneCrossSec,angleB,radiusFracB);
    vector connec = parAVec-parBVec;
    return Foam::mag(connec);
}

Foam::vector Foam::CrossSectionStructure::evaluateRodCircumPos
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    CrossSection* oneCrossSec,
    scalar angle,
    scalar radiusFrac,
    scalar var_para,
    scalar var_radius
)
{
    
    vector d1,d2,d3,r;
    rodEval(oneRod,parameter,d1,d2,d3,r);
    //Info<<"\t\tpara:"<<parameter<<" angle:"<<angle<<" radiusFrac:"<<radiusFrac<<" d1:"<<d1<<" d2:"<<d2<<" d3:"<<d3<<" r:"<<r<<Foam::nl;
    vector tangential = d3;
    scalar tangentialLen = std::sqrt(tangential&tangential);
    tangential /= tangentialLen;
    vector tangentialDev = tangential * var_para;
    //Info<<Foam::nl<<"parameter:"<<parameter<<Foam::nl;
    scalar radius = (*oneCrossSec)(parameter,angle)*radiusFrac;
    radius +=  var_radius;
    //Info<<"radius:"<<radius<<Foam::nl;
    vector coordXDir = std::cos(angle)*radius*d1;
    //Info<<"coordXDir:"<<coordXDir<<"  angle:"<<angle<<"  "<<d2<<Foam::nl;
    vector coordYDir = std::sin(angle)*radius*d2;
    //Info<<"coordYDir:"<<coordYDir<<"  angle:"<<angle<<"  "<<d3<<Foam::nl;
    return (r+tangentialDev)+coordXDir+coordYDir;
}

Foam::vector Foam::CrossSectionStructure::evaluateRodCircumPos
(
    label rodNumber,
    scalar parameter,
    scalar angle,
    scalar radiusFrac,
    scalar var_para,
    scalar var_radius
)
{
    //const ActiveRodMesh::rodCosserat* oneRod = Rods[rodNumber];
    CrossSection* oneCrossSec = &(rodCrossSection[rodNumber]);
    
    vector d1,d2,d3,r;
    rodEval(rodNumber,parameter,d1,d2,d3,r);
    //Info<<"\t\tpara:"<<parameter<<" angle:"<<angle<<" radiusFrac:"<<radiusFrac<<" d1:"<<d1<<" d2:"<<d2<<" d3:"<<d3<<" r:"<<r<<Foam::nl;
    vector tangential = d3;
    scalar tangentialLen = std::sqrt(tangential&tangential);
    tangential /= tangentialLen;
    vector tangentialDev = tangential * var_para;
    //Info<<Foam::nl<<"parameter:"<<parameter<<Foam::nl;
    scalar radius = (*oneCrossSec)(parameter,angle)*radiusFrac;
    radius +=  var_radius;
    //Info<<"radius:"<<radius<<Foam::nl;
    vector coordXDir = std::cos(angle)*radius*d1;
    //Info<<"coordXDir:"<<coordXDir<<"  angle:"<<angle<<"  "<<d2<<Foam::nl;
    vector coordYDir = std::sin(angle)*radius*d2;
    //Info<<"coordYDir:"<<coordYDir<<"  angle:"<<angle<<"  "<<d3<<Foam::nl;
    return (r+tangentialDev)+coordXDir+coordYDir;
}

Foam::Pair<Foam::vector> Foam::CrossSectionStructure::derivateRodCircumPos
(
    label rodNumber,
    scalar parameter,
    scalar angle,
    scalar radiusFrac
)
{
    //const ActiveRodMesh::rodCosserat* oneRod = Rods[rodNumber];
    CrossSection* oneCrossSec = &(rodCrossSection[rodNumber]);
    
    //auto t0 = std::chrono::system_clock::now();
    vector d1,d2,d3,r;
    rodEval(rodNumber,parameter,d1,d2,d3,r);
    
    //auto t1 = std::chrono::system_clock::now();
    
    vector dd1dp,dd2dp,dd3dp,drdp;
    rodDerivEval(rodNumber,parameter,dd1dp,dd2dp,dd3dp,drdp);
    
    //auto t2 = std::chrono::system_clock::now();
    
    scalar radius = (*oneCrossSec)(parameter,angle)*radiusFrac;
    //auto t3 = std::chrono::system_clock::now();
    scalar dradiusdangle = oneCrossSec->deriv_angle(parameter,angle)*radiusFrac;
    //auto t4 = std::chrono::system_clock::now();
    scalar dradiusdp =  oneCrossSec->deriv_para(parameter,angle)*radiusFrac;
    
    //auto t5 = std::chrono::system_clock::now();
    
    vector radVec = std::cos(angle)*d1+std::sin(angle)*d2;
    
    vector dRCPdp = drdp + 
                    dradiusdp*radVec +
                    radius*(std::cos(angle)*dd1dp+std::sin(angle)*dd2dp);
                    
    vector dRCPdangle = dradiusdangle*radVec +
                        radius*(-std::sin(angle)*d1+std::cos(angle)*d2);
                        
    //auto t6 = std::chrono::system_clock::now();
    /*
    Info<<"\t\t d t0-t1 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count()<<Foam::nl;
    Info<<"\t\t d t1-t2 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()<<Foam::nl;
    Info<<"\t\t d t2-t3 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count()<<Foam::nl;
    Info<<"\t\t d t3-t4 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count()<<Foam::nl;
    Info<<"\t\t d t4-t5 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count()<<Foam::nl;
    Info<<"\t\t d t5-t6 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t5).count()<<Foam::nl;
    */
    
    return {dRCPdp,dRCPdangle};
}

Foam::Pair<Foam::vector> Foam::CrossSectionStructure::derivateRodCircumPos
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    CrossSection* oneCrossSec,
    scalar angle,
    scalar radiusFrac
)
{
    //auto t0 = std::chrono::system_clock::now();
    vector d1,d2,d3,r;
    rodEval(oneRod,parameter,d1,d2,d3,r);
    
    //auto t1 = std::chrono::system_clock::now();
    
    vector dd1dp,dd2dp,dd3dp,drdp;
    rodDerivEval(oneRod,parameter,dd1dp,dd2dp,dd3dp,drdp);
    
    //auto t2 = std::chrono::system_clock::now();
    
    scalar radius = (*oneCrossSec)(parameter,angle)*radiusFrac;
    //auto t3 = std::chrono::system_clock::now();
    scalar dradiusdangle = oneCrossSec->deriv_angle(parameter,angle)*radiusFrac;
    //auto t4 = std::chrono::system_clock::now();
    scalar dradiusdp =  oneCrossSec->deriv_para(parameter,angle)*radiusFrac;
    
    //auto t5 = std::chrono::system_clock::now();
    
    vector radVec = std::cos(angle)*d1+std::sin(angle)*d2;
    
    vector dRCPdp = drdp + 
                    dradiusdp*radVec +
                    radius*(std::cos(angle)*dd1dp+std::sin(angle)*dd2dp);
                    
    vector dRCPdangle = dradiusdangle*radVec +
                        radius*(-std::sin(angle)*d1+std::cos(angle)*d2);
                        
    //auto t6 = std::chrono::system_clock::now();
                        
    /*
    Info<<"\t\t d t0-t1 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count()<<Foam::nl;
    Info<<"\t\t d t1-t2 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()<<Foam::nl;
    Info<<"\t\t d t2-t3 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count()<<Foam::nl;
    Info<<"\t\t d t3-t4 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count()<<Foam::nl;
    Info<<"\t\t d t4-t5 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count()<<Foam::nl;
    Info<<"\t\t d t5-t6 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t5).count()<<Foam::nl;
    */
    
    return {dRCPdp,dRCPdangle};
}

Foam::Pair<Foam::vector> Foam::CrossSectionStructure::derivate2RodCircumPos
(
    label rodNumber,
    scalar parameter,
    scalar angle,
    scalar radiusFrac
)
{
    //const ActiveRodMesh::rodCosserat* oneRod = Rods[rodNumber];
    CrossSection* oneCrossSec = &(rodCrossSection[rodNumber]);
    
    //auto t0 = std::chrono::system_clock::now();

    vector d1,d2,d3,r;
    rodEval(rodNumber,parameter,d1,d2,d3,r);
    
    //auto t1 = std::chrono::system_clock::now();
    
    vector dd1dp,dd2dp,dd3dp,drdp;
    rodDerivEval(rodNumber,parameter,dd1dp,dd2dp,dd3dp,drdp);
    
    //auto t2 = std::chrono::system_clock::now();
    
    vector d2d1dp,d2d2dp,d2d3dp,d2rdp;
    rodDeriv2Eval(rodNumber,parameter,d2d1dp,d2d2dp,d2d3dp,d2rdp);
    
    //auto t3 = std::chrono::system_clock::now();

    scalar radius = (*oneCrossSec)(parameter,angle)*radiusFrac;
    scalar dradiusdangle = oneCrossSec->deriv_angle(parameter,angle)*radiusFrac;
    scalar d2radiusdangle = oneCrossSec->deriv2_angle(parameter,angle)*radiusFrac;
    
    //auto t4 = std::chrono::system_clock::now();

    scalar dradiusdp =  oneCrossSec->deriv_para(parameter,angle)*radiusFrac;
    scalar d2radiusdp =  oneCrossSec->deriv2_para(parameter,angle)*radiusFrac;
    
    vector radVec = std::cos(angle)*d1+std::sin(angle)*d2;
    
    //auto t5 = std::chrono::system_clock::now();
    
    vector d2RCPdp = d2rdp + 
                     d2radiusdp*radVec +
                     2*dradiusdp*(std::cos(angle)*dd1dp+std::sin(angle)*dd2dp) +
                     radius*(std::cos(angle)*d2d1dp+std::sin(angle)*d2d2dp);
                    
    vector d2RCPdangle = d2radiusdangle*radVec +
                         2*dradiusdangle*(-std::sin(angle)*d1+std::cos(angle)*d2) + 
                         radius*(-std::cos(angle)*d1-std::sin(angle)*d2);
                         
    //auto t6 = std::chrono::system_clock::now();
    /*
    Info<<"\t\t d2 t0-t1 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count()<<Foam::nl;
    Info<<"\t\t d2 t1-t2 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()<<Foam::nl;
    Info<<"\t\t d2 t2-t3 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count()<<Foam::nl;
    Info<<"\t\t d2 t3-t4 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count()<<Foam::nl;
    Info<<"\t\t d2 t4-t5 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count()<<Foam::nl;
    Info<<"\t\t d2 t5-t6 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t5).count()<<Foam::nl;
    */
    
    return {d2RCPdp,d2RCPdangle};
}

Foam::Pair<Foam::vector> Foam::CrossSectionStructure::derivate2RodCircumPos
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    CrossSection* oneCrossSec,
    scalar angle,
    scalar radiusFrac
)
{
    //auto t0 = std::chrono::system_clock::now();

    vector d1,d2,d3,r;
    rodEval(oneRod,parameter,d1,d2,d3,r);
    
    //auto t1 = std::chrono::system_clock::now();
    
    vector dd1dp,dd2dp,dd3dp,drdp;
    rodDerivEval(oneRod,parameter,dd1dp,dd2dp,dd3dp,drdp);
    
    //auto t2 = std::chrono::system_clock::now();
    
    vector d2d1dp,d2d2dp,d2d3dp,d2rdp;
    rodDeriv2Eval(oneRod,parameter,d2d1dp,d2d2dp,d2d3dp,d2rdp);
    
    //auto t3 = std::chrono::system_clock::now();

    scalar radius = (*oneCrossSec)(parameter,angle)*radiusFrac;
    scalar dradiusdangle = oneCrossSec->deriv_angle(parameter,angle)*radiusFrac;
    scalar d2radiusdangle = oneCrossSec->deriv2_angle(parameter,angle)*radiusFrac;
    
    //auto t4 = std::chrono::system_clock::now();

    scalar dradiusdp =  oneCrossSec->deriv_para(parameter,angle)*radiusFrac;
    scalar d2radiusdp =  oneCrossSec->deriv2_para(parameter,angle)*radiusFrac;
    
    vector radVec = std::cos(angle)*d1+std::sin(angle)*d2;
    
    //auto t5 = std::chrono::system_clock::now();
    
    vector d2RCPdp = d2rdp + 
                     d2radiusdp*radVec +
                     2*dradiusdp*(std::cos(angle)*dd1dp+std::sin(angle)*dd2dp) +
                     radius*(std::cos(angle)*d2d1dp+std::sin(angle)*d2d2dp);
                    
    vector d2RCPdangle = d2radiusdangle*radVec +
                         2*dradiusdangle*(-std::sin(angle)*d1+std::cos(angle)*d2) + 
                         radius*(-std::cos(angle)*d1-std::sin(angle)*d2);
                         
    //auto t6 = std::chrono::system_clock::now();
                         
    /*
    Info<<"\t\t d2 t0-t1 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count()<<Foam::nl;
    Info<<"\t\t d2 t1-t2 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()<<Foam::nl;
    Info<<"\t\t d2 t2-t3 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count()<<Foam::nl;
    Info<<"\t\t d2 t3-t4 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count()<<Foam::nl;
    Info<<"\t\t d2 t4-t5 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count()<<Foam::nl;
    Info<<"\t\t d2 t5-t6 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t5).count()<<Foam::nl;
    */
    
    return {d2RCPdp,d2RCPdangle};
}

Foam::vector Foam::CrossSectionStructure::evalRodDerivCoeff
(
    label rodNumber,
    label derivCoeffNumber,
    scalar parameter,
    scalar angle,
    scalar radiusFrac
)
{
    FatalErrorInFunction<<"Not yet implemented"<<exit(FatalError);
    return vector(0,0,0);
}

Foam::vector Foam::CrossSectionStructure::evaluateRodCircumDerivAngle
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    CrossSection* oneCrossSec,
    scalar angle,
    scalar radiusFrac
)
{
    vector d1,d2,d3,r;
    rodEval(oneRod,parameter,d1,d2,d3,r);
    //Info<<Foam::nl<<"parameter:"<<parameter<<Foam::nl;
    scalar radius = (*oneCrossSec)(parameter,angle)*radiusFrac;
    //Info<<"radius:"<<radius<<Foam::nl;
    scalar dradius_dangle = oneCrossSec->deriv_angle(parameter,angle)*radiusFrac;
    
    //Info<<"radius:"<<radius<<Foam::nl;
    vector coordXDerivAngle = (dradius_dangle*std::cos(angle) - radius*std::sin(angle))*d1;
    //Info<<"coordXDir:"<<coordXDir<<"  angle:"<<angle<<"  "<<d2<<Foam::nl;
    vector coordYDerivAngle = (dradius_dangle*std::sin(angle) + radius*std::cos(angle))*d2;
    //Info<<"coordYDir:"<<coordYDir<<"  angle:"<<angle<<"  "<<d3<<Foam::nl;
    return coordXDerivAngle+coordYDerivAngle;
}

Foam::vector Foam::CrossSectionStructure::evaluateRodCircumDerivAngle
(
    label rodNumber,
    scalar parameter,
    scalar angle,
    scalar radiusFrac
)
{
    CrossSection* oneCrossSec = &(rodCrossSection[rodNumber]);
    
    vector d1,d2,d3,r;
    rodEval(rodNumber,parameter,d1,d2,d3,r);
    //Info<<Foam::nl<<"parameter:"<<parameter<<Foam::nl;
    scalar radius = (*oneCrossSec)(parameter,angle)*radiusFrac;
    //Info<<"radius:"<<radius<<Foam::nl;
    scalar dradius_dangle = oneCrossSec->deriv_angle(parameter,angle)*radiusFrac;
    
    //Info<<"radius:"<<radius<<Foam::nl;
    vector coordXDerivAngle = (dradius_dangle*std::cos(angle) - radius*std::sin(angle))*d1;
    //Info<<"coordXDir:"<<coordXDir<<"  angle:"<<angle<<"  "<<d2<<Foam::nl;
    vector coordYDerivAngle = (dradius_dangle*std::sin(angle) + radius*std::cos(angle))*d2;
    //Info<<"coordYDir:"<<coordYDir<<"  angle:"<<angle<<"  "<<d3<<Foam::nl;
    return coordXDerivAngle+coordYDerivAngle;
}

Foam::scalar Foam::CrossSectionStructure::restrictAngle
(
    scalar angle
)
{
    while(angle<0)
    {
        angle += 2*Foam::constant::mathematical::pi;
    }
    while(angle>2*Foam::constant::mathematical::pi)
    {
        angle -= 2*Foam::constant::mathematical::pi;
    }
    return angle;
}

template<typename T>
T Foam::CrossSectionStructure::integrateCircumwise
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    CrossSection* crossSec,
    std::function<T(scalar)> function
)
{
    return integrateCircumwise<T>(oneRod,parameter,crossSec,0,2*Foam::constant::mathematical::pi,function);
}

template<typename T>
T Foam::CrossSectionStructure::integrateCircumwise
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    CrossSection* crossSec,
    scalar angleStart,
    scalar angleEnd,
    std::function<T(scalar)> function
)
{
    while(angleEnd<angleStart)
        angleEnd += 2*Foam::constant::mathematical::pi;
    
    label highestWaveNumber = (crossSec->getNumberCoeffs()+1)*2;
    
    scalar integrationLen = angleEnd-angleStart;
    scalar segmentLen = integrationLen/highestWaveNumber;
    std::vector<scalar> segmentSet;
    for(label i=0; i<highestWaveNumber; i++)
        segmentSet.push_back(angleStart+segmentLen*i);
    segmentSet.push_back(angleEnd);
    
    label numberOfAbcissa = 20;

    T totalValue = Foam::zero();
    for(uint i=0; i<segmentSet.size()-1;i++)
    {
        scalar startAngle = segmentSet[i];
        scalar endAngle = segmentSet[i+1];
        scalar dist = endAngle-startAngle;
        scalar step = dist/numberOfAbcissa;
        scalar initialStep = startAngle+step/2;
        T summedValue = Foam::zero();
        for(label i=0;i<numberOfAbcissa;i++)
        {
            summedValue += function(initialStep)*step;
            initialStep+=step;
        }
        totalValue+=summedValue;
    }
    return totalValue;
}
        
template<typename T>
T Foam::CrossSectionStructure::integrateCircumwise
(
    label rodNumber,
    scalar parameter,
    scalar angleStart,
    scalar angleEnd,
    std::function<T(scalar)> function
)
{
    return integrateCircumwise<T>(Rods[rodNumber],parameter,&(rodCrossSection[rodNumber]),angleStart,angleEnd,function);
}

Foam::scalar Foam::CrossSectionStructure::distance
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    CrossSection* crossSec,
    scalar angleStart,
    scalar angleEnd,
    scalar radiusFrac
)
{
    std::function<scalar(scalar)> curveLen = 
    [rodPtr=oneRod, para=parameter, crossSecRef=crossSec, radFrac=radiusFrac](scalar angle)
    {
        vector derivCrossSec = evaluateRodCircumDerivAngle(rodPtr,para,crossSecRef,angle,radFrac);
        return std::sqrt(derivCrossSec&derivCrossSec);
    };
    return integrateCircumwise<scalar>(oneRod,parameter,crossSec,angleStart,angleEnd,curveLen);
}

Foam::scalar Foam::CrossSectionStructure::distance
(
    label rodNumber,
    scalar parameter,
    scalar angleStart,
    scalar angleEnd,
    scalar radiusFrac
)
{
    std::function<scalar(scalar)> curveLen = 
    [rodNum=rodNumber, para=parameter, radFrac=radiusFrac, this](scalar angle)
    {
        vector derivCrossSec = evaluateRodCircumDerivAngle(rodNum,para,angle,radFrac);
        return std::sqrt(derivCrossSec&derivCrossSec);
    };
    return integrateCircumwise<scalar>(rodNumber,parameter,angleStart,angleEnd,curveLen);
}

Foam::scalar Foam::CrossSectionStructure::distance
(
    LagrangianMarkerOnCrossSec& A,
    LagrangianMarkerOnCrossSec& B
)
{
    if(A.getBaseRod()!=B.getBaseRod())
        FatalErrorInFunction<<"Distance can not be computed between points on different rods"<<exit(FatalError);
    
    //markerParameter,markerAngle,radiusFrac,markerPosition
    using node = std::tuple<scalar,scalar,scalar,vector>;

    node start = {A.getMarkerParameter(),A.getMarkerAngle(),A.getMarkerRadiusFrac(),A.getMarkerPosition()};
    std::get<1>(start) = restrictAngle(std::get<1>(start));
    node end = {B.getMarkerParameter(),B.getMarkerAngle(),B.getMarkerRadiusFrac(),B.getMarkerPosition()};
    std::get<1>(end) = restrictAngle(std::get<1>(end));
    
    std::list<node> nodes;
    nodes.push_back(start);
    nodes.push_back(end);
    
    scalar distance = std::numeric_limits<scalar>::max();
    bool errorSufficient = false;
    while(!errorSufficient)
    {
        auto node0 = nodes.begin();
        auto node1 = ++nodes.begin();
        for( ; node1!=nodes.end() ; )
        {
            node insertNode;
            std::get<0>(insertNode) = (std::get<0>(*node0)+std::get<0>(*node1))/2;
            std::get<2>(insertNode) = (std::get<2>(*node0)+std::get<2>(*node1))/2;
            
            scalar anglNode0 = std::get<1>(*node0);
            scalar anglNode1 = std::get<1>(*node1);
            scalar anglNodeInter;
            scalar distInner = std::abs(anglNode0-anglNode1);
            if(anglNode0<anglNode1)
            {
                scalar distToUp = 2*Foam::constant::mathematical::pi-anglNode1;
                scalar outerDist = distToUp+anglNode0;
                if(outerDist>distInner)
                    anglNodeInter = (anglNode0+anglNode1)/2;
                else
                {
                    scalar halfOuterDist = outerDist/2;
                    if(distToUp<anglNode0)
                        anglNodeInter = anglNode0-halfOuterDist;
                    else
                        anglNodeInter = anglNode1+halfOuterDist;
                }
            }
            else
            {
                scalar distToUp = 2*Foam::constant::mathematical::pi-anglNode0;
                scalar outerDist = distToUp+anglNode1;
                if(outerDist>distInner)
                    anglNodeInter = (anglNode0+anglNode1)/2;
                else
                {
                    scalar halfOuterDist = outerDist/2;
                    if(distToUp<anglNode1)
                        anglNodeInter = anglNode1-halfOuterDist;
                    else
                        anglNodeInter = anglNode0+halfOuterDist;
                }
            }
            std::get<1>(insertNode) = anglNodeInter;
            std::get<3>(insertNode) = CrossSectionStructure::evaluateRodCircumPos
                                    (A.getBaseRod(),std::get<0>(insertNode),
                                     A.getBaseCrossSec(),std::get<1>(insertNode),
                                     std::get<2>(insertNode));
            nodes.insert(node1,insertNode);
            node0 = node1;
            node1++;
        }
        
        scalar newDist = 0;
        node0 = nodes.begin();
        node1 = ++nodes.begin();
        for( ; node1!=nodes.end() ; )
        {
            vector connec = std::get<3>(*node1) - std::get<3>(*node0);
            newDist += std::sqrt(connec&connec);
            node0++;
            node1++;
        }
        
        scalar change = std::abs(newDist-distance);
        scalar changePerc = change/newDist;
        distance = newDist;
        if(nodes.size()>10 && changePerc<1e-3)
            errorSufficient=true;
    }
    return distance;
}

Foam::scalar Foam::CrossSectionStructure::lowerBound_distance
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    CrossSection* crossSec,
    scalar angleStart,
    scalar angleEnd,
    scalar radiusFrac
)
{
    vector posStart = evaluateRodCircumPos(oneRod,parameter,crossSec,angleStart,radiusFrac);
    vector posEnd = evaluateRodCircumPos(oneRod,parameter,crossSec,angleEnd,radiusFrac);
    vector dist = posEnd-posStart;
    return std::sqrt(dist&dist);
}

Foam::scalar Foam::CrossSectionStructure::lowerBound_distance
(
    label rodNumber,
    scalar parameter,
    scalar angleStart,
    scalar angleEnd,
    scalar radiusFrac
)
{
    vector posStart = evaluateRodCircumPos(rodNumber,parameter,angleStart,radiusFrac);
    vector posEnd = evaluateRodCircumPos(rodNumber,parameter,angleEnd,radiusFrac);
    vector dist = posEnd-posStart;
    return std::sqrt(dist&dist);
}

void Foam::CrossSectionStructure::printMarkerStructure()
{
    for(std::size_t rodNumber=0; rodNumber<rodMarkersList.size(); rodNumber++)
    {
        std::list<std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>>& oneRod = *(rodMarkersList[rodNumber]);
        Info<<"------------------rodNumber:"<<rodNumber<<"--------------------"<<oneRod.size()<<Foam::nl;        
        for(auto iterTang = oneRod.begin(); iterTang!=oneRod.end(); iterTang++)
        {
            std::pair<scalar,std::list<std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>>>& tangRod = *iterTang;
            Info<<"tang:"<<tangRod.first<<"  "<<tangRod.second.size()<<Foam::nl;
            for(auto iterRad = tangRod.second.begin(); iterRad!=tangRod.second.end(); iterRad++)
            {
                std::pair<scalar,std::list<LagrangianMarkerOnCrossSec>>& radRod = *iterRad;
                Info<<"     rad:"<<radRod.first<<"  "<<radRod.second.size()<<Foam::nl;
            }
        }
    }
}

void Foam::CrossSectionStructure::parameterGradientCheck()
{
    LineStructure::parameterGradientCheck();
    Info<<"CrossSectionStructure::parameterGradientCheck (deriv_para,deriv2_para,deriv_angle,deriv2_angle) "<<Foam::nl;
    
    scalar nbrSteps = 20;
    scalar epsilon = 1e-8;
    for(label rodNumber=0; rodNumber<nR; rodNumber++)
    {
        CrossSection* oCS = &(rodCrossSection[rodNumber]);
        
        scalar domainStart = this->domainStart(rodNumber)+5*epsilon;
        scalar domainEnd = this->domainEnd(rodNumber)-5*epsilon;
        scalar delta = domainEnd-domainStart;
        scalar stepsize = delta/nbrSteps;
        for(scalar parameter=domainStart; parameter<=domainEnd; parameter+=stepsize)
        {
            scalar angleStart = 0;
            scalar angleEnd = 2*Foam::constant::mathematical::pi;
            scalar angleDelta = angleEnd-angleStart;
            scalar angleStepsize = angleDelta/nbrSteps;
            for(scalar angle=angleStart; angle<angleEnd; angle+=angleStepsize)
            {
                // Compute gradients
                scalar drdp = oCS->deriv_para(parameter,angle);
                scalar d2rdp = oCS->deriv2_para(parameter,angle);
                scalar drdangle = oCS->deriv_angle(parameter,angle);
                scalar d2rdangle = oCS->deriv2_angle(parameter,angle);
                
                scalar l_para = parameter-epsilon;
                scalar u_para = parameter+epsilon;
                scalar l_angle = angle-epsilon;
                scalar u_angle = angle+epsilon;

                scalar lp_r = (*oCS)(l_para,   angle);
                scalar up_r = (*oCS)(u_para,   angle);
                scalar la_r = (*oCS)(parameter,l_angle);
                scalar ua_r = (*oCS)(parameter,u_angle);
                
                scalar lp_drdp =     oCS->deriv_para(l_para,   angle);
                scalar up_drdp =     oCS->deriv_para(u_para,   angle);
                scalar la_drdangle = oCS->deriv_angle(parameter,l_angle);
                scalar ua_drdangle = oCS->deriv_angle(parameter,u_angle);

                scalar fd_drdp = (up_r-lp_r)/(2*epsilon);
                scalar fd_d2rdp = (up_drdp-lp_drdp)/(2*epsilon);
                scalar fd_drdangle = (ua_r-la_r)/(2*epsilon);
                scalar fd_d2rdangle = (ua_drdangle-la_drdangle)/(2*epsilon);
                
                /*
                Info<<"--------------------------"<<Foam::nl;
                Info<<"parameter:"<<parameter<<Foam::nl;
                Info<<"angle:"<<angle<<Foam::nl;
                Info<<"---"<<Foam::nl;
                Info<<"lp_r:"<<lp_r<<Foam::nl;
                Info<<"up_r:"<<up_r<<Foam::nl;
                Info<<"la_r:"<<la_r<<Foam::nl;
                Info<<"ua_r:"<<ua_r<<Foam::nl;
                Info<<"---"<<Foam::nl;
                Info<<"lp_drdp:"<<lp_drdp<<Foam::nl;
                Info<<"up_drdp:"<<up_drdp<<Foam::nl;
                Info<<"la_drdangle:"<<la_drdangle<<Foam::nl;
                Info<<"ua_drdangle:"<<ua_drdangle<<Foam::nl;
                Info<<"---"<<Foam::nl;
                Info<<"drdp:"<<drdp<<Foam::nl;
                Info<<"d2rdp:"<<d2rdp<<Foam::nl;
                Info<<"drdangle:"<<drdangle<<Foam::nl;
                Info<<"d2rdangle:"<<d2rdangle<<Foam::nl;
                Info<<"---"<<Foam::nl;
                Info<<"fd_drdp:"<<fd_drdp<<Foam::nl;
                Info<<"fd_d2rdp:"<<fd_d2rdp<<Foam::nl;
                Info<<"fd_drdangle:"<<fd_drdangle<<Foam::nl;
                Info<<"fd_d2rdangle:"<<fd_d2rdangle<<Foam::nl;
                Info<<"--------------------------"<<Foam::nl;
                */
                
                scalar err_drdp = std::abs(fd_drdp-drdp);
                scalar denom_drdp = 0.5*(std::abs(fd_drdp)+std::abs(drdp));
                scalar percErr_drdp = (denom_drdp==0)?0:err_drdp/denom_drdp;
                if(percErr_drdp>0.04)
                {
                    Info<<"parameter:"<<parameter<<Foam::nl;
                    Info<<"drdp:"<<drdp<<Foam::nl;
                    Info<<"lp_r:"<<lp_r<<Foam::nl;
                    Info<<"up_r:"<<up_r<<Foam::nl;
                    Info<<"fd_drdp:"<<fd_drdp<<Foam::nl;
                    Info<<"err_drdp:"<<err_drdp<<Foam::nl;
                    Info<<"percErr_drdp:"<<percErr_drdp<<Foam::nl;
                    FatalErrorInFunction<<"Error"<<exit(FatalError);
                }
                scalar err_d2rdp = std::abs(fd_d2rdp-d2rdp);
                scalar denom_d2rdp = 0.5*(std::abs(fd_d2rdp)+std::abs(d2rdp));
                scalar percErr_d2rdp = (denom_d2rdp==0)?0:err_d2rdp/denom_d2rdp;
                if(percErr_d2rdp>0.04)                
                {
                    Info<<"parameter:"<<parameter;
                    Info<<"  d2rdp:"<<d2rdp;
                    //Info<<"lp_drdp:"<<lp_drdp<<Foam::nl;
                    //Info<<"up_drdp:"<<up_drdp<<Foam::nl;
                    Info<<"  fd_d2rdp:"<<fd_d2rdp<<Foam::nl;
                    //Info<<"err_d2rdp:"<<err_d2rdp<<Foam::nl;
                    //Info<<"percErr_d2rdp:"<<percErr_d2rdp<<Foam::nl;
                    //FatalErrorInFunction<<"Error"<<exit(FatalError);
                }
                scalar err_drdangle = std::abs(fd_drdangle-drdangle);
                scalar denom_drdangle = 0.5*(std::abs(fd_drdangle)+std::abs(drdangle));
                scalar percErr_drdangle = (denom_drdangle==0)?0:err_drdangle/denom_drdangle;
                if(percErr_drdangle>0.04)
                {
                    Info<<"angle:"<<angle<<Foam::nl;
                    Info<<"drdangle:"<<drdangle<<Foam::nl;
                    Info<<"la_r:"<<la_r<<Foam::nl;
                    Info<<"ua_r:"<<ua_r<<Foam::nl;
                    Info<<"fd_drdangle:"<<fd_drdangle<<Foam::nl;
                    Info<<"err_drdangle:"<<err_drdangle<<Foam::nl;
                    Info<<"percErr_drdangle:"<<percErr_drdangle<<Foam::nl;
                    FatalErrorInFunction<<"Error"<<exit(FatalError);
                }
                scalar err_d2rdangle = std::abs(fd_d2rdangle-d2rdangle);
                scalar denom_d2rdangle = 0.5*(std::abs(fd_d2rdangle)+std::abs(d2rdangle));
                scalar percErr_d2rdangle = (denom_d2rdangle==0)?0:err_d2rdangle/denom_d2rdangle;
                if(percErr_d2rdangle>0.04)
                {
                    Info<<"angle:"<<angle<<Foam::nl;
                    Info<<"d2rdangle:"<<d2rdangle<<Foam::nl;
                    Info<<"la_drdangle:"<<la_drdangle<<Foam::nl;
                    Info<<"ua_drdangle:"<<ua_drdangle<<Foam::nl;
                    Info<<"fd_d2rdangle:"<<fd_d2rdangle<<Foam::nl;
                    Info<<"err_d2rdangle:"<<err_d2rdangle<<Foam::nl;
                    Info<<"percErr_d2rdangle:"<<percErr_d2rdangle<<Foam::nl;
                    FatalErrorInFunction<<"Error"<<exit(FatalError);
                }
            }
        }
                
        for(scalar parameter=domainStart; parameter<=domainEnd; parameter+=stepsize)
        {
            scalar angleStart = 0;
            scalar angleEnd = 2*Foam::constant::mathematical::pi;
            scalar angleDelta = angleEnd-angleStart;
            scalar angleStepsize = angleDelta/nbrSteps;
            for(scalar angle=angleStart; angle<angleEnd; angle+=angleStepsize)
            {
                scalar radiusStart = 0;
                scalar radiusEnd = 1;
                scalar radiusDelta = radiusEnd-radiusStart;
                scalar radiusStepsize = radiusDelta/nbrSteps;
                for(scalar radFrac=radiusStart; radFrac<radiusEnd; radFrac+=radiusStepsize)
                {
                    // Compute gradients
                    Pair<vector> drd_ = derivateRodCircumPos(Rods[rodNumber],parameter,oCS,angle,radFrac);
                    Pair<vector> d2rd_ = derivate2RodCircumPos(Rods[rodNumber],parameter,oCS,angle,radFrac);
                    
                    scalar l_para = parameter-epsilon;
                    scalar l_angle = angle-epsilon;
                    scalar u_para = parameter+epsilon;
                    scalar u_angle = angle+epsilon;

                    vector lp_r = evaluateRodCircumPos(Rods[rodNumber],l_para,oCS,angle,radFrac);
                    vector la_r = evaluateRodCircumPos(Rods[rodNumber],parameter,oCS,l_angle,radFrac);
                    vector up_r = evaluateRodCircumPos(Rods[rodNumber],u_para,oCS,angle,radFrac);
                    vector ua_r = evaluateRodCircumPos(Rods[rodNumber],parameter,oCS,u_angle,radFrac);
                    
                    vector lp_drdp = derivateRodCircumPos(Rods[rodNumber],l_para,oCS,angle,radFrac).first();
                    vector la_drdangle = derivateRodCircumPos(Rods[rodNumber],parameter,oCS,l_angle,radFrac).second();
                    vector up_drdp = derivateRodCircumPos(Rods[rodNumber],u_para,oCS,angle,radFrac).first();
                    vector ua_drdangle = derivateRodCircumPos(Rods[rodNumber],parameter,oCS,u_angle,radFrac).second();

                    Pair<vector> fd_drd_((up_r-lp_r)/(2*epsilon),(ua_r-la_r)/(2*epsilon));
                    Pair<vector> fd_d2rd_((up_drdp-lp_drdp)/(2*epsilon),(ua_drdangle-la_drdangle)/(2*epsilon));

                    vector err_drdpVec = fd_drd_.first()-drd_.first();
                    scalar err_drdp = std::sqrt(err_drdpVec&err_drdpVec);
                    scalar denom_drdp = 0.5*(std::sqrt(fd_drd_.first()&fd_drd_.first())+std::sqrt(drd_.first()&drd_.first()));
                    scalar percErr_drdp = (denom_drdp==0)?0:err_drdp/denom_drdp;
                    if(percErr_drdp>0.01)
                    {
                        Info<<"parameter:"<<parameter<<Foam::nl;
                        Info<<"drd_.first():"<<drd_.first()<<Foam::nl;
                        Info<<"lp_r:"<<lp_r<<Foam::nl;
                        Info<<"up_r:"<<up_r<<Foam::nl;
                        Info<<"fd_drd_.first():"<<fd_drd_.first()<<Foam::nl;
                        Info<<"err_drdp:"<<err_drdp<<Foam::nl;
                        Info<<"percErr_drdp:"<<percErr_drdp<<Foam::nl;
                        FatalErrorInFunction<<"Error"<<exit(FatalError);
                    }
                    
                    vector err_drdangleVec = fd_drd_.second()-drd_.second();
                    scalar err_drdangle = std::sqrt(err_drdangleVec&err_drdangleVec);
                    scalar denom_drdangle = 0.5*(std::sqrt(fd_drd_.first()&fd_drd_.first())+std::sqrt(drd_.first()&drd_.first()));
                    scalar percErr_drdangle = (denom_drdangle==0)?0:err_drdangle/denom_drdangle;
                    if(percErr_drdangle>0.01)
                    {
                        Info<<"parameter:"<<parameter<<Foam::nl;
                        Info<<"drd_.second():"<<drd_.second()<<Foam::nl;
                        Info<<"la_r:"<<la_r<<Foam::nl;
                        Info<<"ua_r:"<<ua_r<<Foam::nl;
                        Info<<"fd_drd_.second():"<<fd_drd_.second()<<Foam::nl;
                        Info<<"err_drdangle:"<<err_drdangle<<Foam::nl;
                        Info<<"percErr_drdangle:"<<percErr_drdangle<<Foam::nl;
                        FatalErrorInFunction<<"Error"<<exit(FatalError);
                    }
                    
                    vector err_d2rdpVec = fd_d2rd_.first()-d2rd_.first();
                    scalar err_d2rdp = std::sqrt(err_d2rdpVec&err_d2rdpVec);
                    scalar denom_d2rdp = 0.5*(std::sqrt(fd_d2rd_.first()&fd_d2rd_.first())+std::abs(d2rd_.first()&d2rd_.first()));
                    scalar percErr_d2rdp = (denom_d2rdp==0)?0:err_d2rdp/denom_d2rdp;
                    if(percErr_d2rdp>0.01)
                    {
                        Info<<"parameter:"<<parameter;
                        Info<<" d2rdp:"<<fd_d2rd_.first();
                        //Info<<"lp_drdp:"<<lp_drdp<<Foam::nl;
                        //Info<<"up_drdp:"<<up_drdp<<Foam::nl;
                        Info<<" fd_d2rdp:"<<fd_drd_.first();//<<Foam::nl;
                        //Info<<"err_d2rdp:"<<err_d2rdp<<Foam::nl;
                        Info<<" percErr:"<<percErr_d2rdp<<Foam::nl;
                        //FatalErrorInFunction<<"Error"<<exit(FatalError);
                    }
                    
                    vector err_d2rdangleVec = fd_d2rd_.second()-d2rd_.second();
                    scalar err_d2rdangle = std::sqrt(err_d2rdangleVec&err_d2rdangleVec);
                    scalar denom_d2rdangle = 0.5*(std::sqrt(fd_d2rd_.first()&fd_d2rd_.first())+std::sqrt(d2rd_.first()&d2rd_.first()));
                    scalar percErr_d2rdangle = (denom_d2rdangle==0)?0:err_d2rdangle/denom_d2rdangle;
                    if(percErr_d2rdangle>0.01)
                    {
                        Info<<"parameter:"<<parameter;//<<Foam::nl;
                        Info<<" d2rda:"<<d2rd_.second();//<<Foam::nl;
                        //Info<<"la_drdangle:"<<la_drdangle<<Foam::nl;
                        //Info<<"ua_drdangle:"<<ua_drdangle<<Foam::nl;
                        Info<<" fd_d2rda:"<<fd_d2rd_.second();//<<Foam::nl;
                        //Info<<"err_d2rdangle:"<<err_d2rdangle<<Foam::nl;
                        Info<<" percErr:"<<percErr_d2rdangle<<Foam::nl;
                        //FatalErrorInFunction<<"Error"<<exit(FatalError);
                    }
                }
            }
        }
    }
}

void Foam::CrossSectionStructure::selfCheck()
{
    //LineStructure::selfCheck();
    
    Info<<"-----------Check cross Section derivatives-----------"<<Foam::nl;    
    std::function<void(std::function<scalar(scalar,scalar)>,
                       std::function<scalar(scalar,scalar)>,
                       scalar,scalar,scalar,scalar,uint)> crossSecComparer =
    [](auto deriv, auto fdDeriv, scalar minPar, scalar maxPar, scalar minAngle, scalar maxAngle, uint steps)
    {
        scalar deltaPar = maxPar-minPar;
        scalar stepsizePar = deltaPar/steps;
        scalar deltaAngle = maxAngle-minAngle;
        scalar stepsizeAngle = deltaAngle/steps;
        
        for(scalar currPara=minPar; currPara<=maxPar; currPara+=stepsizePar)
        {
            for(scalar currAngle=minAngle; currAngle<=maxAngle; currAngle+=stepsizeAngle)
            {
                scalar derivValue = deriv(currPara,currAngle);
                scalar fderivValue = fdDeriv(currPara,currAngle);
                scalar error = std::abs(derivValue-fderivValue);
                scalar avg = 0.5*(derivValue+fderivValue);
                scalar percError;
                if(avg!=0)
                    percError = error/std::abs(avg);
                else
                    percError = error;
                if(percError>2e-3 && error>1e-4)
                {
                    Info<<"("<<currPara<<" -- "<<currAngle<<"): Err:"<<error<<" percErr:"<<percError<<" // grad:"<<derivValue<<" -- fdGrad:"<<fderivValue<<Foam::nl;
                    FatalErrorInFunction<<"Comparison failed!"<<exit(FatalError);
                }
            }
        }
    };
        
    const std::vector<CrossSection>& crossSec = getRodCrossSections();
    for(std::size_t rodNumber=0; rodNumber<crossSec.size(); rodNumber++)
    {
        //Info<<"rodNumber:"<<rodNumber<<Foam::nl;
        
        CrossSection cpCrossSec = crossSec[rodNumber];
        /*
        scalar domainStart = cpCrossSec.domainStart();
        scalar domainEnd = cpCrossSec.domainEnd();
        */
        
        //Info<<"   numberFourierCoeff:"<<cpCrossSec.numberFourierCoeff()<<Foam::nl;
        
        /*
        Info<<"Changing fourierCoeff"<<Foam::nl;
        for(label fourierCoeff=0; fourierCoeff<cpCrossSec.numberFourierCoeff(); fourierCoeff++)
        {
            //Info<<"      numberNurbsCoeffs:"<<cpCrossSec.numberFourierCoeffNurbsCoeffs(fourierCoeff)<<Foam::nl;
            //std::cout<<"      FourierCoeff:"<<cpCrossSec.getCurve(fourierCoeff)<<std::endl;
            for(label coeffI=0; coeffI<cpCrossSec.numberFourierCoeffNurbsCoeffs(fourierCoeff); coeffI++)
            {
                auto deriv = [&](scalar par, scalar angle)
                {
                    return cpCrossSec.evalRadiusDerivFourierCoeffNurbsCoeff(fourierCoeff,coeffI,par,angle);
                };
                
                auto fdDeriv = [&](scalar par, scalar angle, scalar epsilon=1e-8)
                {
                    CrossSection crossSec = cpCrossSec;
                    
                    scalar coeffBasicValue = cpCrossSec.getFourierCoeffNurbsCoeff(fourierCoeff,coeffI);
                    //Info<<"coeffBasicValue:"<<coeffBasicValue<<Foam::nl;
                    
                    scalar lowerCoeffValue = coeffBasicValue-epsilon;
                    //Info<<"coeffBasicValue:"<<coeffBasicValue<<Foam::nl;
                    crossSec.setFourierCoeffNurbsCoeff(fourierCoeff,coeffI,lowerCoeffValue);
                    scalar lowerR = crossSec(par,angle);
                    //Info<<"lowerR:"<<lowerR<<Foam::nl;
                    
                    scalar upperCoeffValue = coeffBasicValue+epsilon;                    //Info<<"coeffBasicValue:"<<coeffBasicValue<<Foam::nl;
                    crossSec.setFourierCoeffNurbsCoeff(fourierCoeff,coeffI,upperCoeffValue);
                    scalar upperR = crossSec(par,angle);
                    //Info<<"upperR:"<<upperR<<Foam::nl;
                    
                    scalar diffR = upperR-lowerR;
                    //Info<<"diffR:"<<diffR<<Foam::nl;
                    scalar diffCoeff = upperCoeffValue-lowerCoeffValue;
                    //Info<<"diffCoeff:"<<diffCoeff<<Foam::nl;
                    scalar deriv;
                    if(diffCoeff!=0)
                         return diffR/diffCoeff;
                    else
                        FatalErrorInFunction<<"No delta in coeff values!"<<exit(FatalError);
                };
                
                crossSecComparer(deriv,fdDeriv,domainStart,domainEnd,0,6.3,20);
            }
        }
        //std::cout<<"      Phase:"<<cpCrossSec.getPhaseCurve()<<std::endl;
        Info<<"Changing phase"<<Foam::nl;
        for(label phaseCoeff=0; phaseCoeff<cpCrossSec.numberPhaseNurbsCoeffs(); phaseCoeff++)
        {
            auto deriv = [&](scalar par, scalar angle)
            {
                return cpCrossSec.evalRadiusDerivPhaseNurbsCoeff(phaseCoeff,par,angle);
            };
            
            auto fdDeriv = [&](scalar par, scalar angle, scalar epsilon=1e-8)
            {
                CrossSection crossSec = cpCrossSec;
                
                scalar coeffBasicValue = cpCrossSec.getPhaseNurbsCoeff(phaseCoeff);
                
                scalar lowerCoeffValue = coeffBasicValue-epsilon;
                crossSec.setPhaseNurbsCoeff(phaseCoeff,lowerCoeffValue);
                scalar lowerR = crossSec(par,angle);
                
                scalar upperCoeffValue = coeffBasicValue+epsilon;
                crossSec.setPhaseNurbsCoeff(phaseCoeff,upperCoeffValue);
                scalar upperR = crossSec(par,angle);
                
                scalar diffR = upperR-lowerR;
                scalar diffCoeff = upperCoeffValue-lowerCoeffValue;
                scalar deriv;
                if(diffCoeff!=0)
                    return diffR/diffCoeff;
                else
                    FatalErrorInFunction<<"No delta in coeff values!"<<exit(FatalError);
                
                return (upperR-lowerR)/(upperCoeffValue-lowerCoeffValue);
            };
            
            crossSecComparer(deriv,fdDeriv,domainStart,domainEnd,0,6.3,20);
        }
        */
        
        /*
        Info<<"Cross section distance"<<Foam::nl;
        scalar dist = distance(Rods[rodNumber],0.5,&cpCrossSec,0,0.5,1);
        Info<<"dist:"<<dist<<Foam::nl;
        vector angle_0 = evaluateRodCircumPos(Rods[rodNumber],0.5,&cpCrossSec,0,1);
        vector angle_05 = evaluateRodCircumPos(Rods[rodNumber],0.5,&cpCrossSec,0.5,1);
        vector connec = angle_05-angle_0;
        scalar len = std::sqrt(connec&connec);
        
        Info<<"angle_0:"<<angle_0<<Foam::nl;
        Info<<"angle_05:"<<angle_05<<Foam::nl;
        Info<<"len:"<<len<<Foam::nl;
        */
        scalar radius00 = cpCrossSec(0.25,0);
        scalar radius01 = cpCrossSec(0.25,0.01);
        scalar radius02 = cpCrossSec(0.25,0.02);
        scalar radius03 = cpCrossSec(0.25,0.03);
        scalar radius04 = cpCrossSec(0.25,0.04);
        
        Info<<"r000:"<<radius00<<Foam::endl;
        Info<<"r001:"<<radius01<<Foam::endl;
        Info<<"r002:"<<radius02<<Foam::endl;
        Info<<"r003:"<<radius03<<Foam::endl;
        Info<<"r004:"<<radius04<<Foam::endl;
        
        scalar dradius00 = cpCrossSec.deriv_angle(0.25,0);
        scalar dradius01 = cpCrossSec.deriv_angle(0.25,0.01);
        scalar dradius02 = cpCrossSec.deriv_angle(0.25,0.02);
        scalar dradius03 = cpCrossSec.deriv_angle(0.25,0.03);
        scalar dradius04 = cpCrossSec.deriv_angle(0.25,0.04);
        
        Info<<"dr000:"<<dradius00<<Foam::endl;
        Info<<"dr001:"<<dradius01<<Foam::endl;
        Info<<"dr002:"<<dradius02<<Foam::endl;
        Info<<"dr003:"<<dradius03<<Foam::endl;
        Info<<"dr004:"<<dradius04<<Foam::endl;
        
        scalar d2radius00 = cpCrossSec.deriv2_angle(0.25,0);
        scalar d2radius01 = cpCrossSec.deriv2_angle(0.25,0.01);
        scalar d2radius02 = cpCrossSec.deriv2_angle(0.25,0.02);
        scalar d2radius03 = cpCrossSec.deriv2_angle(0.25,0.03);
        scalar d2radius04 = cpCrossSec.deriv2_angle(0.25,0.04);
        
        Info<<"d2r000:"<<d2radius00<<Foam::endl;
        Info<<"d2r001:"<<d2radius01<<Foam::endl;
        Info<<"d2r002:"<<d2radius02<<Foam::endl;
        Info<<"d2r003:"<<d2radius03<<Foam::endl;
        Info<<"d2r004:"<<d2radius04<<Foam::endl;
        
        radius00 = cpCrossSec(0.25,0);
        radius01 = cpCrossSec(0.25,0.01);
        radius02 = cpCrossSec(0.25,0.02);
        radius03 = cpCrossSec(0.25,0.03);
        radius04 = cpCrossSec(0.25,0.04);
        
        Info<<"r000:"<<radius00<<Foam::endl;
        Info<<"r001:"<<radius01<<Foam::endl;
        Info<<"r002:"<<radius02<<Foam::endl;
        Info<<"r003:"<<radius03<<Foam::endl;
        Info<<"r004:"<<radius04<<Foam::endl;
        
        dradius00 = cpCrossSec.deriv_angle(0.25,0);
        dradius01 = cpCrossSec.deriv_angle(0.25,0.01);
        dradius02 = cpCrossSec.deriv_angle(0.25,0.02);
        dradius03 = cpCrossSec.deriv_angle(0.25,0.03);
        dradius04 = cpCrossSec.deriv_angle(0.25,0.04);
        
        Info<<"dr000:"<<dradius00<<Foam::endl;
        Info<<"dr001:"<<dradius01<<Foam::endl;
        Info<<"dr002:"<<dradius02<<Foam::endl;
        Info<<"dr003:"<<dradius03<<Foam::endl;
        Info<<"dr004:"<<dradius04<<Foam::endl;
        
        d2radius00 = cpCrossSec.deriv2_angle(0.25,0);
        d2radius01 = cpCrossSec.deriv2_angle(0.25,0.01);
        d2radius02 = cpCrossSec.deriv2_angle(0.25,0.02);
        d2radius03 = cpCrossSec.deriv2_angle(0.25,0.03);
        d2radius04 = cpCrossSec.deriv2_angle(0.25,0.04);
        
        Info<<"d2r000:"<<d2radius00<<Foam::endl;
        Info<<"d2r001:"<<d2radius01<<Foam::endl;
        Info<<"d2r002:"<<d2radius02<<Foam::endl;
        Info<<"d2r003:"<<d2radius03<<Foam::endl;
        Info<<"d2r004:"<<d2radius04<<Foam::endl;
        
        vector rodPosStatic = evaluateRodCircumPos(Rods[rodNumber],0.25,&cpCrossSec,0.01,1);
        vector rodPos = evaluateRodCircumPos(rodNumber,0.25,0.01,1);
        
        Pair<vector> drodPosStatic = derivateRodCircumPos(Rods[rodNumber],0.25,&cpCrossSec,0.01,1);
        Pair<vector> drodPos = derivateRodCircumPos(rodNumber,0.25,0.01,1);
        
        Pair<vector> d2rodPosStatic = derivate2RodCircumPos(Rods[rodNumber],0.25,&cpCrossSec,0.01,1);
        Pair<vector> d2rodPos = derivate2RodCircumPos(rodNumber,0.25,0.01,1);
        
        Info<<"rodPosStatic:"<<rodPosStatic<<Foam::endl;
        Info<<"rodPos      :"<<rodPos<<Foam::endl;
        Info<<"drodPosStatic:"<<drodPosStatic<<Foam::endl;
        Info<<"drodPos      :"<<drodPos<<Foam::endl;
        Info<<"d2rodPosStatic:"<<d2rodPosStatic<<Foam::endl;
        Info<<"d2rodPos      :"<<d2rodPos<<Foam::endl;
        
        Info<<"rodPosStatic:"<<rodPosStatic<<Foam::endl;
        Info<<"rodPos      :"<<rodPos<<Foam::endl;
        Info<<"drodPosStatic:"<<drodPosStatic<<Foam::endl;
        Info<<"drodPos      :"<<drodPos<<Foam::endl;
        Info<<"d2rodPosStatic:"<<d2rodPosStatic<<Foam::endl;
        Info<<"d2rodPos      :"<<d2rodPos<<Foam::endl;
        
    }
    
    /*
    Info<<"-----------Cross Section structure derivatives-----------"<<Foam::nl;
    for(label rodNumber=0; rodNumber<getNumberRods(); rodNumber++)
    {
        Info<<"rodNumber:"<<rodNumber<<Foam::nl;
        scalar domainStart = this->domainStart(rodNumber);
        scalar domainEnd = this->domainEnd(rodNumber);
        Info<<"   numberCurveCoeffs:"<<numberCurveCoeffs(rodNumber)<<Foam::nl;
        for(label coeffI=0; coeffI<numberCurveCoeffs(rodNumber); coeffI++)
        {
            for(label dim=0; dim<3; dim++)
            {
                auto deriv = [&](scalar par)
                {
                    vector d1,d2,d3,r;
                    rodEvalDerivCoeff(rodNumber,coeffI,dim,par,d1,d2,d3,r);
                    return FixedList<vector,4>({d1,d2,d3,r});
                };
                
                auto fdDeriv = [&](scalar par, scalar epsilon=1e-10)
                {                   
                    scalar coeffBasicValue = getCurveCoeff(rodNumber,coeffI,dim);
                    
                    scalar lowerCoeffValue = coeffBasicValue-epsilon;
                    setCurveCoeff(rodNumber,coeffI,dim,lowerCoeffValue);
                    vector lowerd1,lowerd2,lowerd3,lowerr;
                    rodEval(Rods[rodNumber],par,lowerd1,lowerd2,lowerd3,lowerr);
                    
                    scalar upperCoeffValue = coeffBasicValue+epsilon;
                    setCurveCoeff(rodNumber,coeffI,dim,upperCoeffValue);
                    vector upperd1,upperd2,upperd3,upperr;
                    rodEval(Rods[rodNumber],par,upperd1,upperd2,upperd3,upperr);
                    
                    setCurveCoeff(rodNumber,coeffI,dim,coeffBasicValue);
                        
                    vector d1 = (upperd1-lowerd1)/(upperCoeffValue-lowerCoeffValue);
                    vector d2 = (upperd2-lowerd2)/(upperCoeffValue-lowerCoeffValue);
                    vector d3 = (upperd3-lowerd3)/(upperCoeffValue-lowerCoeffValue);
                    vector r = (upperr-lowerr)/(upperCoeffValue-lowerCoeffValue);
                    return FixedList<vector,4>({d1,d2,d3,r});
                };
            }
        }
    }
    */
    //FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
}
