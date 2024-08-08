#include "CrossSectionStructure.H"

Foam::CrossSectionStructure::CrossSectionStructure
(
    dynamicRefineFvMesh& mesh,
    std::vector<CrossSection> rodCrossSection,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
LineStructure(mesh,modusFieldToMarker,modusMarkerToField),
rodCrossSection(rodCrossSection)
{
    Info<<"CrossSectionStructure"<<Foam::endl;
    initialize();
}

void Foam::CrossSectionStructure::to_string()
{
    Info<<"---CrossSectionStructure::Markers---"<<Foam::endl;
    int count=0;
    for(std::unique_ptr<std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>>& oneRodMarkers : rodMarkersList)
    {
        Info<<"Rod "<<count<<Foam::endl;
        if(oneRodMarkers)
        {
            for(std::list<std::list<LagrangianMarkerOnCrossSec>>& radialMarkers : *oneRodMarkers)
            {
                Info<<"\t-"<<Foam::endl;
                for(std::list<LagrangianMarkerOnCrossSec>& angleMarkers : radialMarkers)
                {
                    for(const LagrangianMarkerOnCrossSec& marker : angleMarkers)
                    {
                        Info<<"\t"<<marker.to_string()<<Foam::endl;
                    }
                }
            }
        }
    }
    Info<<"------------------------------------"<<Foam::endl;
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
    const CrossSection& crossSec = rodCrossSection[rodNumber];
    
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
        scalar totalDistanceFrac = totalDistance/spacing;
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
            scalar radiusSpacing = upperLimitRadius/spacing;
            label numberOfFracs = std::ceil(radiusSpacing);
            scalar radFracPart = 1.0/numberOfFracs;
            scalar radFrac = 1.0;
            radialData.resize(numberOfFracs+1);
                        
            for(label r=0; r<numberOfFracs; r++)
            {
                radialData[r].first = radFrac;
                createSpacedPointsOnCrossSec(oneRod,parameter,&crossSec,radFrac,spacing,radialData[r].second);
                //Info<<radFrac<<"  radialData["<<r<<"].second.size():"<<radialData[r].second.size()<<Foam::endl;
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
            //Info<<"radialData[0].second.size():"<<radialData[0].second.size()<<Foam::endl;
        }
    }
    initialRodPoints[rodNumber] = std::move(pointsPtr);
}

void Foam::CrossSectionStructure::createSpacedPointsOnCrossSec
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    const CrossSection* oneCrossSec,
    scalar radFrac,
    scalar spacing,
    std::vector<scalar>& angleData
)
{    
    std::list<scalar> points;
    points.push_back(0);
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
            if(dist>spacing)
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
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    const CrossSection& crossSec = rodCrossSection[rodNumber];
    
    auto markersPtr = std::unique_ptr<std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>>
    (
        new std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>()
    );
    std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>& oneRodMarkers = *markersPtr;
    
    //[rod] -> {para , [] -> {rad , [] -> angle}}
    const std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>& thisRodIniPnts = *(initialRodPoints[rodNumber]);
    
    oneRodMarkers.resize(thisRodIniPnts.size());
    auto iterPara = oneRodMarkers.begin();
    for(uint paraInd=0; paraInd<thisRodIniPnts.size(); paraInd++,iterPara++)
    {
        scalar parameter = thisRodIniPnts[paraInd].first;
        const std::vector<std::pair<scalar,std::vector<scalar>>>& radialPnts = thisRodIniPnts[paraInd].second;
        std::list<std::list<LagrangianMarkerOnCrossSec>>& radialMarkers = *iterPara;
        radialMarkers.resize(radialPnts.size());
        auto iterRad = radialMarkers.begin();
        for(uint radInd=0; radInd<radialPnts.size(); radInd++,iterRad++)
        {
            scalar radialFrac = radialPnts[radInd].first;
            const std::vector<scalar>& anglePnts = radialPnts[radInd].second;
            std::list<LagrangianMarkerOnCrossSec>& angleMarkers = *iterRad;
            //angleMarkers.resize(anglePnts.size());
            auto iterAngle = angleMarkers.begin();
            for(uint angleInd=0; angleInd<anglePnts.size(); angleInd++,iterAngle++)
            {
                scalar angle = anglePnts[angleInd];
                angleMarkers.push_back
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

void Foam::CrossSectionStructure::refineMarkersOnRod
(
    label rodNumber,
    std::pair<bool,scalar> forcedSpacing
)
{
    //const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>& markers = *(rodMarkersList[rodNumber]);
    
    for(auto iterPara=markers.begin(); iterPara!=markers.end(); iterPara++)
    {
        for(auto iterRadial=iterPara->begin(); iterRadial!=iterPara->end(); iterRadial++)
        {
            if(iterRadial->size()<1)
                FatalErrorInFunction<<"Must be at least one marker in list"<< exit(FatalError);
            if(iterRadial->front().getMarkerRadiusFrac()!=0)
                refineCircumferential(*iterRadial);
        }
    }
    
    // Refine rod endings radially
    refineRadial(markers.front(),initialMeshSpacing);
    refineRadial(markers.back(),initialMeshSpacing);
    
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

void Foam::CrossSectionStructure::refineCircumferential
(
    std::list<LagrangianMarkerOnCrossSec>& circumMarkers,
    std::pair<bool,scalar> refineSpacing
)
{
    if(circumMarkers.size()<2)
    {
        Info<<circumMarkers.front().getMarkerRadiusFrac()<<Foam::endl;
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    }
    
    label rodNumber = circumMarkers.front().getRodNumber();
    const ActiveRodMesh::rodCosserat* oneRod = circumMarkers.front().getBaseRod();
    scalar parameter = circumMarkers.front().getMarkerParameter();
    scalar radFrac = circumMarkers.front().getMarkerRadiusFrac();
    const CrossSection* oneCrossSec = circumMarkers.front().getBaseCrossSec();

    bool refined=true;
    while(refined)
    {
        refined = false;
        auto circMarkerIter0 = circumMarkers.begin();
        auto circMarkerIter1 = ++(circumMarkers.begin());
        for( ; circMarkerIter1!=circumMarkers.end() ; )
        {
            scalar circMarker0Angle = circMarkerIter0->getMarkerParameter();
            scalar circMarker0CellSpacing = circMarkerIter0->getMarkerCellMinSpacing();
            bool circMarker0InCell = (circMarkerIter0->getMarkerCell()!=-1);
            
            scalar circMarker1Angle = circMarkerIter1->getMarkerParameter();
            scalar circMarker1CellSpacing = circMarkerIter1->getMarkerCellMinSpacing();
            bool circMarker1InCell = (circMarkerIter1->getMarkerCell()!=-1);
            
            scalar dist = distance
            (
                oneRod,parameter,oneCrossSec,circMarker0Angle,circMarker1Angle,radFrac
            );
            bool subdivide = false;
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
}

void Foam::CrossSectionStructure::refineRadial
(
    std::list<std::list<LagrangianMarkerOnCrossSec>>& radialMarkers,
    scalar initialSpacing,
    std::pair<bool,scalar> refineSpacing
)
{
    if(radialMarkers.size()<2)
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    if(radialMarkers.front().size()<2)
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    
    label rodNumber = radialMarkers.front().front().getRodNumber();
    const ActiveRodMesh::rodCosserat* oneRod = radialMarkers.front().front().getBaseRod();
    scalar parameter = radialMarkers.front().front().getMarkerParameter();
    const CrossSection* oneCrossSec = radialMarkers.front().front().getBaseCrossSec();
    
    for(auto iter=radialMarkers.begin(); iter!=radialMarkers.end(); iter++)
    {
        if(iter->size()<1)
            FatalErrorInFunction<<"Must be at least one marker in list"<< exit(FatalError);
        if(iter->front().getMarkerRadiusFrac()!=0)
            refineCircumferential(*iter,refineSpacing);
    }

    bool refined=true;
    while(refined)
    {
        refined = false;
        auto radMarkerIter0 = radialMarkers.begin();
        auto radMarkerIter1 = ++(radialMarkers.begin());
        for( ; radMarkerIter1!=radialMarkers.end() ; )
        {
            if(radMarkerIter0->size()<1)
            {
                Info<<radMarkerIter0->front().getMarkerParameter()<<"/"<<radMarkerIter0->front().getMarkerRadiusFrac()<<Foam::endl;
                FatalErrorInFunction<<"Must be at least one markers"<< exit(FatalError);
            }
            scalar radMarker0RadFrac = radMarkerIter0->front().getMarkerRadiusFrac();
            scalar radMarker0CellSpacing = std::numeric_limits<scalar>::max();
            bool radMarker0InCell = false;
            for(auto iter=radMarkerIter0->begin(); iter!=radMarkerIter0->end(); iter++)
            {
                radMarker0CellSpacing = std::min(radMarker0CellSpacing,iter->getMarkerCellMinSpacing());
                radMarker0InCell |= (iter->getMarkerCell()!=-1);
            }
            
            if(radMarkerIter1->size()<1)
            {
                Info<<radMarkerIter1->front().getMarkerParameter()<<"/"<<radMarkerIter1->front().getMarkerRadiusFrac()<<Foam::endl;
                FatalErrorInFunction<<"Must be at least one markers"<< exit(FatalError);
            }
            scalar radMarker1RadFrac = radMarkerIter1->front().getMarkerRadiusFrac();
            scalar radMarker1CellSpacing = std::numeric_limits<scalar>::max();
            bool radMarker1InCell = false;
            for(auto iter=radMarkerIter1->begin(); iter!=radMarkerIter1->end(); iter++)
            {
                radMarker1CellSpacing = std::min(radMarker1CellSpacing,iter->getMarkerCellMinSpacing());
                radMarker1InCell |= (iter->getMarkerCell()!=-1);
            }
            
            scalar maxDist = std::numeric_limits<scalar>::max();
            for(auto iter=radMarkerIter0->begin(); iter!=radMarkerIter0->end(); iter++)
            {
                scalar angle = iter->getMarkerAngle();
                vector radMarker0Pos = evaluateRodCircumPos
                (
                    oneRod,parameter,oneCrossSec,angle,radMarker0CellSpacing
                );
                vector radMarker1Pos = evaluateRodCircumPos
                (
                    oneRod,parameter,oneCrossSec,angle,radMarker1CellSpacing
                );
                vector distVec = radMarker0Pos-radMarker1Pos;
                maxDist = std::max(maxDist,std::sqrt(distVec&distVec));
            }
            for(auto iter=radMarkerIter1->begin(); iter!=radMarkerIter1->end(); iter++)
            {
                scalar angle = iter->getMarkerAngle();
                vector radMarker0Pos = evaluateRodCircumPos
                (
                    oneRod,parameter,oneCrossSec,angle,radMarker0CellSpacing
                );
                vector radMarker1Pos = evaluateRodCircumPos
                (
                    oneRod,parameter,oneCrossSec,angle,radMarker1CellSpacing
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
                if(maxDist>minSpacing)
                    subdivide=true;
            }

            if(subdivide)
            {
                scalar middleRadiusFrac = 0.5*(radMarker0RadFrac + radMarker1RadFrac);
                std::vector<scalar> angleData;
                createSpacedPointsOnCrossSec
                (
                    oneRod,parameter,oneCrossSec,middleRadiusFrac,initialSpacing,angleData
                );
                std::list<LagrangianMarkerOnCrossSec> angleMarkers;
                for(scalar angle : angleData)
                {
                    LagrangianMarkerOnCrossSec middleRadiusMarker
                    (
                        *this,mesh,rodNumber,oneRod,parameter,oneCrossSec,angle,middleRadiusFrac
                    );
                    angleMarkers.push_back(middleRadiusMarker);
                }
                refineCircumferential(angleMarkers,refineSpacing);                
                radialMarkers.insert(radMarkerIter1,angleMarkers);
                refined=true;
            }
            radMarkerIter0 = radMarkerIter1;
            radMarkerIter1++;
        }
    }
}

void Foam::CrossSectionStructure::refineTangential
(
    std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>& tangMarkers,
    scalar initialSpacing,
    std::pair<bool,scalar> refineSpacing
)
{
    if(tangMarkers.size()<2)
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    if(tangMarkers.front().size()<2)
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    if(tangMarkers.front().front().size()<2)
        FatalErrorInFunction<<"Must be at least two markers"<< exit(FatalError);
    
    label rodNumber = tangMarkers.front().front().front().getRodNumber();
    const ActiveRodMesh::rodCosserat* oneRod = tangMarkers.front().front().front().getBaseRod();
    const CrossSection* oneCrossSec = tangMarkers.front().front().front().getBaseCrossSec();
    
    bool refined=true;
    while(refined)
    {
        refined = false;
        auto tangMarkerIter0 = tangMarkers.begin();
        auto tangMarkerIter1 = ++(tangMarkers.begin());
        for( ; tangMarkerIter1!=tangMarkers.end() ; )
        {
            if(tangMarkerIter0->size()<1)
                FatalErrorInFunction<<"Must be at least one radial marker layer"<< exit(FatalError);
            if(tangMarkerIter0->front().size()<1)
                FatalErrorInFunction<<"Must be at least one angle marker layer"<< exit(FatalError);
            scalar tangMarker0Para = tangMarkerIter0->front().front().getMarkerParameter();
            scalar tangMarker0CellSpacing = std::numeric_limits<scalar>::max();
            bool tangMarker0InCell = false;
            for(auto iter=tangMarkerIter0->front().begin(); iter!=tangMarkerIter0->front().end(); iter++)
            {
                tangMarker0CellSpacing = std::min(tangMarker0CellSpacing,iter->getMarkerCellMinSpacing());
                tangMarker0InCell |= (iter->getMarkerCell()!=-1);
            }
            
            if(tangMarkerIter1->size()<1)
                FatalErrorInFunction<<"Must be at least one radial marker layer"<< exit(FatalError);
            if(tangMarkerIter1->front().size()<1)
                FatalErrorInFunction<<"Must be at least one angle marker layer"<< exit(FatalError);
            scalar tangMarker1Para = tangMarkerIter1->front().front().getMarkerParameter();
            scalar tangMarker1CellSpacing = std::numeric_limits<scalar>::max();
            bool tangMarker1InCell = false;
            for(auto iter=tangMarkerIter1->front().begin(); iter!=tangMarkerIter1->front().end(); iter++)
            {
                tangMarker1CellSpacing = std::min(tangMarker1CellSpacing,iter->getMarkerCellMinSpacing());
                tangMarker1InCell |= (iter->getMarkerCell()!=-1);
            }
            
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
                if(totalDistance>minSpacing)
                    subdivide=true;
            }

            if(subdivide)
            {
                scalar middleMarkerPara = 0.5*(tangMarker0Para + tangMarker1Para);
                std::vector<scalar> angleData;
                createSpacedPointsOnCrossSec
                (
                    oneRod,middleMarkerPara,oneCrossSec,1.0,initialSpacing,angleData
                );
                std::list<LagrangianMarkerOnCrossSec> angleMarkers;
                for(scalar angle : angleData)
                    angleMarkers.push_back(
                        LagrangianMarkerOnCrossSec
                        (
                            *this,mesh,rodNumber,oneRod,middleMarkerPara,oneCrossSec,angle,1.0
                        ));
                refineCircumferential(angleMarkers,refineSpacing);
                std::list<std::list<LagrangianMarkerOnCrossSec>> radialMarkers;
                radialMarkers.push_back(angleMarkers);
                tangMarkers.insert(tangMarkerIter1,radialMarkers);
                refined=true;
            }
            tangMarkerIter0 = tangMarkerIter1;
            tangMarkerIter1++;
        }
    }
}

void Foam::CrossSectionStructure::setMarkerVolumeOnRod
(
    label rodNumber
)
{
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    const CrossSection* oneCrossSec = &(rodCrossSection[rodNumber]);
    
    std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>& markersList = *(rodMarkersList[rodNumber]);
    std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec*>>> markers(markersList.size());
    uint P=0;
    for(auto iterPara=markersList.begin(); iterPara!=markersList.end(); iterPara++,P++)
    {
        std::list<std::list<LagrangianMarkerOnCrossSec>>& radialMarkersList = *iterPara;
        std::vector<std::vector<LagrangianMarkerOnCrossSec*>>& radialMarkers = markers[P];
        radialMarkers.resize(radialMarkersList.size());
        uint R=0;
        for(auto iterRadial=radialMarkersList.begin(); iterRadial!=radialMarkersList.end(); iterRadial++,R++)
        {
            std::list<LagrangianMarkerOnCrossSec>& angleMarkersList = *iterRadial;
            std::vector<LagrangianMarkerOnCrossSec*>& angleMarkers = radialMarkers[R];
            //angleMarkers.resize(radialMarkersList.size());
            for(auto iterAngle=angleMarkersList.begin(); iterAngle!=angleMarkersList.end(); iterAngle++)
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
        //Info<<"rodWiseIndex:"<<rodWiseIndex<<Foam::endl;
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
                surfaceType = Surface::Circumferential;
            }
            else 
            {
                surfaceType = Surface::Radial;
                if(radialIndex==0)
                    outerRadiusFrac = radialIndexFrac;
                else
                    outerRadiusFrac = (radialMarkers[radialIndex-1][0]->getMarkerRadiusFrac()+radialIndexFrac)/2;
                
                if(radialIndex==radialMarkers.size()-1)
                    innerRadiusFrac = radialIndexFrac;
                else
                    innerRadiusFrac = (radialMarkers[radialIndex+1][0]->getMarkerRadiusFrac()+radialIndexFrac)/2;
            }
            //Info<<"\tradialIndex:"<<radialIndex<<" - "<<circumMarkers.size()<<Foam::endl;
            for(uint circIndex=0; circIndex<circumMarkers.size(); circIndex++)
            {
                //Info<<"\t\tcircIndex:"<<circIndex<<"/"<<circumMarkers.size()<<Foam::endl;

                LagrangianMarkerOnCrossSec* singleMarker = circumMarkers[circIndex];
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
                //Info<<"\t\tdistAngle:"<<"("<<lowerAngle<<"/"<<circIndexAngle<<"/"<<upperAngle<<"):"<<upperAngle-lowerAngle;
                upperAngle += circIndexAngle;
                upperAngle /= 2;
                upperAngle = CrossSectionStructure::restrictAngle(upperAngle);
                
                if(singleMarker->getMarkerCell()!=-1)
                {
                    if(modusMarkerToField == markerMeshType::Uniform)
                    {
                        scalar h = std::cbrt(cells[singleMarker->getMarkerCell()].mag(points,faces));
                        
                        /*
                        Info<<"prevParameter:"<<prevParameter<<"  subseqParameter:"<<subseqParameter<<Foam::endl;
                        Info<<"innerRadiusFrac:"<<innerRadiusFrac<<"  outerRadiusFrac:"<<outerRadiusFrac<<Foam::endl;
                        Info<<"lowerAngle:"<<lowerAngle<<"  upperAngle:"<<upperAngle<<Foam::endl;
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
                        
                        singleMarker->setMarkerVolume(thisCell.mag(points,faces));
                        //Info<<"|||"<<singleMarker.getMarkerVolume()<<"|||"<<singleMarker.getSupportCells().size()<<Foam::endl;
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
                        
                        singleMarker->setMarkerVolume(thisFace.mag(points));
                        //Info<<"|||"<<singleMarker.getMarkerVolume()<<"|||"<<singleMarker.getSupportCells().size()<<Foam::endl;
                        sumVolume+=singleMarker->getMarkerVolume();
                    }
                }
            }
        }
    }
}

void Foam::CrossSectionStructure::evaluateMarkerMeshRelation()
{
    status.execValid(status.markerMesh);
    for(std::unique_ptr<std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>>& singleRodMarkers :  rodMarkersList)
    {
        for(std::list<std::list<LagrangianMarkerOnCrossSec>>& oneParaRodMarkers : *singleRodMarkers)
        {
            for(std::list<LagrangianMarkerOnCrossSec>& radialFracRodMarkers : oneParaRodMarkers)
            {
                evaluateMarkerMeshRelation(radialFracRodMarkers);
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
    std::cout<<"CrossSectionStructure::reduceMarkers()"<<std::endl;
    
    std::vector<MarkerReference<LagrangianMarkerOnCrossSec>> allMarkers;
    for(std::unique_ptr<std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>>& singleRodMarkers :  rodMarkersList)
    {
        for(std::list<std::list<LagrangianMarkerOnCrossSec>>& oneParaRodMarkers : *singleRodMarkers)
        {
            for(std::list<LagrangianMarkerOnCrossSec>& radialFracRodMarkers : oneParaRodMarkers)
            {
                std::list<LagrangianMarkerOnCrossSec>* singleRadMarkers = &(radialFracRodMarkers);
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
    LineStructure::reduceMarkers(allMarkers);
}

BoundingBox Foam::CrossSectionStructure::computeBox
(
    label rodNumber
)
{
    BoundingBox nurbsCurveBox = Structure::computeBox(rodNumber);
    std::pair<scalar,scalar> minMax = rodCrossSection[rodNumber].radiusBounds();
    nurbsCurveBox.enlarge(minMax.second);
    return nurbsCurveBox;
}

BoundingBox Foam::CrossSectionStructure::computeBox
(
    label rodNumber,
    scalar parStart,
    scalar parEnd
)
{
    BoundingBox nurbsCurveBox = Structure::computeBox(rodNumber,parStart,parEnd);
    std::pair<scalar,scalar> minMax = rodCrossSection[rodNumber].radiusBounds(parStart,parEnd);
    nurbsCurveBox.enlarge(minMax.second);
    return nurbsCurveBox;
}

scalar Foam::CrossSectionStructure::characteristicSize
(
    label rodNumber,
    scalar par
)
{
    return rodCrossSection[rodNumber].upperLimitRadius(par);
}

void Foam::CrossSectionStructure::removeOverlapMarkers()
{
    std::vector<MarkerReference<LagrangianMarkerOnCrossSec>> deleteMarkers;
    for(std::size_t rodI=0; rodI<rodMarkersList.size(); rodI++)
    {
        if(rodInMesh[rodI])
        {
            std::unique_ptr<std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>>& singleRodMarkers = rodMarkersList[rodI];
            for(std::list<std::list<LagrangianMarkerOnCrossSec>>& oneParaRodMarkers : *singleRodMarkers)
            {
                for(std::list<LagrangianMarkerOnCrossSec>& radialFracRodMarkers : oneParaRodMarkers)
                {
                    std::list<LagrangianMarkerOnCrossSec>* singleRadMarkers = &(radialFracRodMarkers);
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
                                    const CrossSection& crossSec = rodCrossSection[rodIOther];
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
    std::cout<<"CrossSectionStructure::collectMarkers"<<std::endl;

    status.execValid(status.markersCollected);
    collectedMarkers.resize(0);    
    for(std::unique_ptr<std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>>& singleRodMarkers :  rodMarkersList)
    {
        if(singleRodMarkers)
        {
            for(std::list<std::list<LagrangianMarkerOnCrossSec>>& oneParaRodMarkers : *singleRodMarkers)
            {
                for(std::list<LagrangianMarkerOnCrossSec>& radialFracRodMarkers : oneParaRodMarkers)
                {
                    for(LagrangianMarkerOnCrossSec& marker : radialFracRodMarkers)
                    {
                        collectedMarkers.push_back(&marker);
                    }
                }
            }
        }
        else
            FatalErrorInFunction<<"Rod with no markers given"<<exit(FatalError);
    }
    status.executed(status.markersCollected);
    std::cout<<"CrossSectionStructure::collectMarkers done"<<std::endl;
}

/*
void Foam::CrossSectionStructure::collectHaloMarkers()
{
    const std::unordered_map<label,label>& selfHaloCellToIndex = getHaloCellToIndexMap(Pstream::myProcNo());
    
    for(std::unique_ptr<std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>>>& singleRodMarkers :  rodMarkersList)
    {
        for(std::list<std::list<LagrangianMarkerOnCrossSec>>& oneParaRodMarkers : *singleRodMarkers)
        {
            for(std::list<LagrangianMarkerOnCrossSec>& radialFracRodMarkers : oneParaRodMarkers)
            {
                for(const LagrangianMarkerOnCrossSec& marker : radialFracRodMarkers)
                {
                    label cellOfMarker = marker.getMarkerCell();
                    auto iter = selfHaloCellToIndex.find(cellOfMarker);
                    if(iter!=selfHaloCellToIndex.end())
                    {
                        label index = iter->second;
                        haloCellsRodMarkersList[index].push_back(&marker);
                    }
                }
            }
        }
    }
}
*/

/*
std::unique_ptr<std::vector<LagrangianMarkerOnCrossSec>> Foam::CrossSectionStructure::constructMarkerSet
(
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* oneRod,
    const CrossSection* oneCrossSec,
    scalar initialSpacing,
    bool reInitialize,
    std::pair<bool,scalar> refineSpacing
)
{
    //Create initial spacing
    std::unique_ptr<std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>>& spacedPoints = initialRodPoints[rodNumber];
    if(reInitialize || !spacedPoints)
    {
        spacedPoints = createSpacedPointsOnCrossSecRod(rodNumber,oneCrossSec,initialSpacing);
    }
        
    // Transfer vector to list
    //(para (radial (angle)))
    std::list<std::list<std::list<LagrangianMarkerOnCrossSec>>> markers;
    for(label paraInd=0; paraInd<spacedPoints->size(); paraInd++)
    {
        scalar para = (*spacedPoints)[paraInd].first;
        const std::vector<std::pair<scalar,std::vector<scalar>>>& radialPoints = (*spacedPoints)[paraInd].second;
        
        std::list<std::list<LagrangianMarkerOnCrossSec>> nextPara;
        markers.push_back(nextPara);
        
        //Info<<paraInd<<" para:"<<para<<"  "<<radialPoints.size()<<Foam::endl;
        
        for(label radialInd=0; radialInd<radialPoints.size(); radialInd++)
        {
            scalar radiusFrac = radialPoints[radialInd].first;
            const std::vector<scalar>& circumPoints = radialPoints[radialInd].second;
            
            std::list<LagrangianMarkerOnCrossSec> nextRadial;
            markers.back().push_back(nextRadial);
            
            //Info<<"     "<<radialInd<<"     radiusFrac:"<<radiusFrac<<"  "<<circumPoints.size()<<Foam::endl;
        
            for(label angleInd=0; angleInd<circumPoints.size(); angleInd++)
            {
                scalar angle = circumPoints[angleInd];
                
                //Info<<"     "<<angleInd<<"             angle:"<<angle<<Foam::endl;
                
                LagrangianMarkerOnCrossSec nextAngleMarker
                (
                    mesh,rodNumber,oneRod,para,oneCrossSec,angle,radiusFrac
                );
                markers.back().back().push_back(nextAngleMarker);
            }
        }
    }
    Info<<"Transfered to vector"<<Foam::endl;
        
    //Refine circumferential
    for(auto iterPara=markers.begin(); iterPara!=markers.end(); iterPara++)
    {
        for(auto iterRadial=iterPara->begin(); iterRadial!=iterPara->end(); iterRadial++)
        {
            if(iterRadial->size()<1)
                FatalErrorInFunction<<"Must be at least one marker in list"<< exit(FatalError);
            if(iterRadial->front().getMarkerRadiusFrac()!=0)
                refineCircumferential(*iterRadial);
        }
    }
    Info<<"Circum refined"<<Foam::endl;
    
    // Refine rod endings radially
    refineRadial(markers.front(),initialSpacing);
    refineRadial(markers.back(),initialSpacing);
    
    refineTangential(markers,initialSpacing);
    
    Info<<"markers.size():"<<markers.size()<<Foam::endl;
    for(auto radial : markers)
    {
        Info<<"     radial.size():"<<radial.size()<<Foam::endl;
        for(auto circum : radial)
        {
            Info<<"         circum.size():"<<circum.size()<<Foam::endl;
        }
    }    

    
    //Transfer list to vector data structure
    std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>> markersVec;
    for(auto iter=markers.begin(); iter!=markers.end(); iter++)
    {
        markersVec.push_back(std::vector<std::vector<LagrangianMarkerOnCrossSec>>());
        for(auto iterSub=iter->begin(); iterSub!=iter->end(); iterSub++)
        {
            markersVec.back().push_back(std::vector<LagrangianMarkerOnCrossSec>());
            for(auto iterSubSub=iterSub->begin(); iterSubSub!=iterSub->end(); iterSubSub++)
            {
                markersVec.back().back().push_back(*iterSubSub);
            }
        }
    }
    
    computeMarkerVolume(oneRod,oneCrossSec,markersVec);

    std::unique_ptr<std::vector<LagrangianMarkerOnCrossSec>> reducedMarkersPtr = reduceSurplusMarkers(markersVec);
    
    FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);


    return reducedMarkersPtr;
}
*/

/*
std::unique_ptr<std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>> Foam::CrossSectionStructure::createSpacedPointsOnCrossSecRod
(
    label rodNumber,
    const CrossSection* crossSec,            
    scalar spacing
)
{
    Info<<"createSpacedPointsOnCrossSec"<<Foam::endl;
    const ActiveRodMesh::rodCosserat* oneRod = myMesh->m_Rods[rodNumber];
    std::unique_ptr<std::vector<scalar>> spacedPointsOnRod;// = createSpacedPointsOnRod(rodNumber,spacing);
    
    std::list<scalar> onRodPoints;
    onRodPoints.insert(onRodPoints.end(),spacedPointsOnRod->begin(),spacedPointsOnRod->end());
    
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
        
        scalar maxRadiusPnt0 = crossSec->upperLimitRadius(pnt0Para);
        scalar addedDistPnt0 = maxRadiusPnt0*std::sin(halfAngle);
        scalar maxRadiusPnt1 = crossSec->upperLimitRadius(pnt1Para);
        scalar addedDistPnt1 = maxRadiusPnt1*std::sin(halfAngle);
        
        scalar totalDistance = rodDistance+addedDistPnt0+addedDistPnt1;
        scalar totalDistanceFrac = totalDistance/spacing;
        label nbrOfAddIntersec = std::ceil(totalDistanceFrac)-1;
        
        scalar subDistPara = (pnt1Para-pnt0Para)/(nbrOfAddIntersec+1);
        std::vector<scalar> subParas;
        for(label i=1;i<nbrOfAddIntersec+1;i++)
            subParas.push_back(pnt0Para+subDistPara*i);
        onRodPoints.insert(pntIter1,subParas.begin(),subParas.end());
        
        pntIter0 = pntIter1;
        pntIter1++;
    }

    Info<<"createSpacedPointsOnRod done"<<Foam::endl;
    Info<<onRodPoints.size()<<Foam::endl;

    auto result = std::unique_ptr<std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>>(new std::vector<std::pair<scalar,std::vector<std::pair<scalar,std::vector<scalar>>>>>(onRodPoints.size()));
    uint index=0;
    Info<<"createSpacedPointsOnCrossSecRod"<<Foam::endl;
    for(auto iter=onRodPoints.begin(); iter!=onRodPoints.end(); iter++,index++)
    {
        scalar parameter = *iter;
        (*result)[index].first = parameter;
    
        std::vector<std::pair<scalar,std::vector<scalar>>>& radialData = (*result)[index].second;
                
        if(index==0 || index==onRodPoints.size()-1)
        {
            scalar upperLimitRadius = crossSec->upperLimitRadius(parameter);
            scalar radiusSpacing = upperLimitRadius/spacing;
            label numberOfFracs = std::ceil(radiusSpacing);
            Info<<"End case:"<<numberOfFracs<<Foam::endl;
            scalar radFracPart = 1.0/numberOfFracs;
            scalar radFrac = 1.0;
            radialData.resize(numberOfFracs+1);
            
            for(label r=0; r<numberOfFracs; r++)
            {
                radialData[r].first = radFrac;
                createSpacedPointsOnCrossSec(oneRod,parameter,crossSec,radFrac,spacing,radialData[r].second);
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
            createSpacedPointsOnCrossSec(oneRod,parameter,crossSec,1.0,spacing,radialData[0].second);
        }
    }
    Info<<"createSpacedPointsOnCrossSecRod done"<<Foam::endl;
    return result;
}
*/

/*
std::unique_ptr<std::vector<LagrangianMarkerOnCrossSec>> Foam::CrossSectionStructure::reduceSurplusMarkers
(
    std::vector<std::vector<std::vector<LagrangianMarkerOnCrossSec>>>& markers
)
{
    Info<<"reduceSurplusMarkers"<<Foam::endl;
    std::unordered_map<label,DynamicList<LagrangianMarkerOnCrossSec*>> cellToMarker;
    for(std::vector<std::vector<LagrangianMarkerOnCrossSec>>& radMarker : markers)
    {
        for(std::vector<LagrangianMarkerOnCrossSec>& circMarker : radMarker)
        {
            for(LagrangianMarkerOnCrossSec& marker : circMarker)
            {
                if(marker.getMarkerCell()!=-1)
                    cellToMarker[marker.getMarkerCell()].append(&marker);
            }
        }
    }
    
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    Info<<"structured"<<Foam::endl;
    
    auto reducededMarkers = std::unique_ptr<std::vector<LagrangianMarkerOnCrossSec>>(new std::vector<LagrangianMarkerOnCrossSec>());
    for(auto iter=cellToMarker.begin(); iter!=cellToMarker.end(); iter++)
    {
        label cellInd = iter->first;
        const cell& oneCell = cells[cellInd];
        vector cellCentre = oneCell.centre(points,faces);
        std::pair<
            std::pair<
                std::pair<
                    DynamicList<LagrangianMarkerOnCrossSec*>,
                    DynamicList<LagrangianMarkerOnCrossSec*>
                >,
                std::pair<
                    DynamicList<LagrangianMarkerOnCrossSec*>,
                    DynamicList<LagrangianMarkerOnCrossSec*>
                >
            >,
            std::pair<
                std::pair<
                    DynamicList<LagrangianMarkerOnCrossSec*>,
                    DynamicList<LagrangianMarkerOnCrossSec*>
                >,
                std::pair<
                    DynamicList<LagrangianMarkerOnCrossSec*>,
                    DynamicList<LagrangianMarkerOnCrossSec*>
                >
            >
        > xyz_split;
        for(LagrangianMarkerOnCrossSec* markerPtr : iter->second)
        {
            vector markerPos = markerPtr->getMarkerPosition();
            
            std::pair<
                std::pair<
                    DynamicList<LagrangianMarkerOnCrossSec*>,
                    DynamicList<LagrangianMarkerOnCrossSec*>
                >,
                std::pair<
                    DynamicList<LagrangianMarkerOnCrossSec*>,
                    DynamicList<LagrangianMarkerOnCrossSec*>
                >
            >* x_side;
            if(markerPos[0]<cellCentre[0])
                x_side = &(xyz_split.first);
            else
                x_side = &(xyz_split.second);
            
            std::pair<
                DynamicList<LagrangianMarkerOnCrossSec*>,
                DynamicList<LagrangianMarkerOnCrossSec*>
            >* xy_side;
            if(markerPos[1]<cellCentre[1])
                xy_side = &(x_side->first);
            else
                xy_side = &(x_side->second);
            
            DynamicList<LagrangianMarkerOnCrossSec*>* xyz_side;
            if(markerPos[2]<cellCentre[2])
                xyz_side = &(xy_side->first);
            else
                xyz_side = &(xy_side->second);
            
            xyz_side->append(markerPtr);
        }

        std::vector<DynamicList<LagrangianMarkerOnCrossSec*>*> subCells = 
        {
            &(xyz_split.first.first.first),
            &(xyz_split.first.first.second),
            &(xyz_split.first.second.first),
            &(xyz_split.first.second.second),
            &(xyz_split.second.first.first),
            &(xyz_split.second.first.second),
            &(xyz_split.second.second.first),
            &(xyz_split.second.second.first)
        };
        
        for(DynamicList<LagrangianMarkerOnCrossSec*>* subCell : subCells)
        {
            if(subCell->size()>0)
            {
                vector averagePosition = vector(0,0,0);
                scalar totalVolume = 0;
                for(LagrangianMarkerOnCrossSec* marker : *subCell)
                {
                    averagePosition += marker->getMarkerPosition();
                    totalVolume += marker->getMarkerVolume();
                }
                averagePosition /= subCell->size();
                
                LagrangianMarkerOnCrossSec* optMarker;
                scalar optMarkerDistToAvgPos = std::numeric_limits<scalar>::max();
                for(LagrangianMarkerOnCrossSec* marker : *subCell)
                {
                    Info<<"markerPtr:"<<marker<<Foam::endl;
                    vector distVec = marker->getMarkerPosition() - averagePosition;
                    Info<<"distVec:"<<distVec<<Foam::endl;
                    scalar dist = std::sqrt(distVec&distVec);
                    Info<<"dist:"<<dist<<Foam::endl;
                    if(dist<optMarkerDistToAvgPos)
                    {
                        optMarker = marker;
                        optMarkerDistToAvgPos = dist;
                    }
                    Info<<"optMarkerDistToAvgPos:"<<optMarkerDistToAvgPos<<Foam::endl;
                }
                reducededMarkers->push_back(*optMarker);
            }
        }
    }
    return reducededMarkers;
}
*/

Foam::scalar Foam::CrossSectionStructure::evaluateCircumArcLen
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameterA,
    scalar parameterB,
    const CrossSection* oneCrossSec,
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
    const CrossSection* oneCrossSec,
    scalar angle,
    scalar radiusFrac,
    scalar var_para,
    scalar var_radius
)
{
    vector d1,d2,d3,r;
    rodEval(oneRod,parameter,d1,d2,d3,r);
    vector tangential = d3;
    scalar tangentialLen = std::sqrt(tangential&tangential);
    tangential /= tangentialLen;
    vector tangentialDev = tangential * var_para;
    //Info<<Foam::endl<<"parameter:"<<parameter<<Foam::endl;
    scalar radius = (*oneCrossSec)(parameter,angle)*radiusFrac;
    radius +=  var_radius;
    //Info<<"radius:"<<radius<<Foam::endl;
    vector coordXDir = std::cos(angle)*radius*d1;
    //Info<<"coordXDir:"<<coordXDir<<"  angle:"<<angle<<"  "<<d2<<Foam::endl;
    vector coordYDir = std::sin(angle)*radius*d2;
    //Info<<"coordYDir:"<<coordYDir<<"  angle:"<<angle<<"  "<<d3<<Foam::endl;
    return (r+tangentialDev)+coordXDir+coordYDir;
}


Foam::vector Foam::CrossSectionStructure::evaluateRodCircumDerivAngle
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    const CrossSection* oneCrossSec,
    scalar angle,
    scalar radiusFrac
)
{
    vector d1,d2,d3,r;
    rodEval(oneRod,parameter,d1,d2,d3,r);
    //Info<<Foam::endl<<"parameter:"<<parameter<<Foam::endl;
    scalar radius = (*oneCrossSec)(parameter,angle)*radiusFrac;
    //Info<<"radius:"<<radius<<Foam::endl;
    scalar dradius_dangle = oneCrossSec->deriv_angle(parameter,angle)*radiusFrac;
    
    //Info<<"radius:"<<radius<<Foam::endl;
    vector coordXDerivAngle = (dradius_dangle*std::cos(angle) - radius*std::sin(angle))*d1;
    //Info<<"coordXDir:"<<coordXDir<<"  angle:"<<angle<<"  "<<d2<<Foam::endl;
    vector coordYDerivAngle = (dradius_dangle*std::sin(angle) + radius*std::cos(angle))*d2;
    //Info<<"coordYDir:"<<coordYDir<<"  angle:"<<angle<<"  "<<d3<<Foam::endl;
    return coordXDerivAngle+coordYDerivAngle;
}

scalar Foam::CrossSectionStructure::restrictAngle
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
    const CrossSection* crossSec,
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
    const CrossSection* crossSec,
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

scalar Foam::CrossSectionStructure::distance
(
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    const CrossSection* crossSec,
    scalar angleStart,
    scalar angleEnd,
    scalar radiusFrac
)
{
    //Info<<"computeDist:"<<angleStart<<"-"<<angleEnd<<Foam::endl;
    std::function<scalar(scalar)> curveLen = 
    [rodPtr=oneRod, para=parameter, crossSecRef=crossSec, radFrac=radiusFrac](scalar angle)
    {
        vector derivCrossSec = evaluateRodCircumDerivAngle(rodPtr,para,crossSecRef,angle,radFrac);
        return std::sqrt(derivCrossSec&derivCrossSec);
    };
    return integrateCircumwise<scalar>(oneRod,parameter,crossSec,angleStart,angleEnd,curveLen);
}

scalar Foam::CrossSectionStructure::distance
(
    const LagrangianMarkerOnCrossSec& A,
    const LagrangianMarkerOnCrossSec& B
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

/*
void Foam::CrossSectionStructure::constructMarkerSetRadial
(
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* oneRod,
    const scalar parameter,
    const CrossSection* oneCrossSec,
    std::list<std::list<LagrangianMarkerOnCrossSec>>& radialMarkers
)
{
    if(!radialMarkers.empty())
        FatalErrorInFunction<<"List must be empty"<< exit(FatalError);
    
    //Info<<"constructMarkerSetRadial:"<<parameter<<Foam::endl;
    
    // Insert outer circum markers
    radialMarkers.push_back(createInitialCircumMarkers(rodNumber,oneRod,parameter,oneCrossSec,1));
    constructMarkerSetCircumferential(rodNumber,oneRod,oneCrossSec,radialMarkers.back());
    
    //Info<<"Done outer cirum markers"<<Foam::endl;

    // Insert inner circum marker
    LagrangianMarkerOnCrossSec marker
    (
        mesh,
        rodNumber,
        oneRod,
        parameter,
        oneCrossSec,
        0,
        0
    );
    radialMarkers.push_back({marker});
    
    //Info<<"Done inner cirum markers"<<Foam::endl;

    
    bool refined=true;
    label refinementCount = 0;
    while(refined)
    {
        refined = false;
        bool cond = true;
        auto markersIter0 = radialMarkers.begin();
        auto markersIter1 = ++radialMarkers.begin();
        for( ; markersIter1!=radialMarkers.end() ; )
        {
            bool subdivide = doSubdivision(*markersIter0, *markersIter1);
            //Info<<"\t\tsubdivide:"<<subdivide<<Foam::endl;
            if(refinementCount<minRefinement)
                subdivide = true;
            if(subdivide)
            {
                scalar middleRadiusFrac = 0.5*(markersIter0->front().getMarkerRadiusFrac()+markersIter1->front().getMarkerRadiusFrac());
                //Info<<"\t\t"<<markersIter0->front().getMarkerRadiusFrac()<<"->"<<markersIter1->front().getMarkerRadiusFrac()<<" | middleRadiusFrac:"<<middleRadiusFrac<<Foam::endl;
                auto inserted = radialMarkers.insert
                (
                    markersIter1,
                    createInitialCircumMarkers(rodNumber,oneRod,parameter,oneCrossSec,middleRadiusFrac)
                );
                constructMarkerSetCircumferential(rodNumber,oneRod,oneCrossSec,*inserted,middleRadiusFrac);
                refined=true;
            }
            markersIter0 = markersIter1;
            markersIter1++;
        }
        refinementCount++;
        //Info<<"\tconstructMarkerSetRadial refinementCount:"<<refinementCount<<Foam::endl;
        if(refinementCount==7)
            FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);    
    }
}
*/

/*
void Foam::CrossSectionStructure::constructMarkerSetCircumferential
(
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* oneRod,
    const CrossSection* oneCrossSec,
    std::list<LagrangianMarkerOnCrossSec>& circleMarkers,
    scalar radiusFrac
)
{    
    //Info<<"constructMarkerSetCircumferential:"<<radiusFrac<<Foam::endl;
    
    if(circleMarkers.size()!=2)
        FatalErrorInFunction<<"Must start with 2 markers"<< exit(FatalError);
    scalar parameter = circleMarkers.front().getMarkerParameter();
    
    bool refined=true;
    label refinementCount = 0;
    while(refined)
    {
        refined = false;
        bool cond = true;
        auto markersIter0 = circleMarkers.begin();
        auto markersIter1 = ++circleMarkers.begin();
        for( ; markersIter1!=circleMarkers.end() ; )
        {
            bool subdivide = doSubdivisionCircumferential(*markersIter0, *markersIter1);
            if(refinementCount<minRefinement)
                subdivide = true;
            if(subdivide)
            {
                scalar middleAngle = 0.5*(markersIter0->getMarkerAngle()+markersIter1->getMarkerAngle());
                //Info<<parameter<<"\t"<<markersIter0->getMarkerAngle()<<markersIter0->markerPosition<<"->"<<markersIter1->getMarkerAngle()<<markersIter1->markerPosition<<" | middleAngle:"<<middleAngle<<"  ||||||||"<<markersIter0->getSupportCells().size()<<"|"<<markersIter1->getSupportCells().size()<<Foam::endl;
                LagrangianMarkerOnCrossSec middleMarker
                    (mesh,rodNumber,oneRod,parameter,oneCrossSec,middleAngle,radiusFrac);
                auto inserted = circleMarkers.insert
                    (
                        markersIter1,
                        middleMarker
                    );
                refined=true;
            }
            markersIter0 = markersIter1;
            markersIter1++;
        }
        bool subdivide = doSubdivisionCircumferential(circleMarkers.back(),circleMarkers.front());
        if(refinementCount<minRefinement)
            subdivide = true;
        if(subdivide)
        {
            scalar middleAngle = 0.5*(circleMarkers.back().getMarkerAngle()+2*Foam::constant::mathematical::pi);
            //Info<<parameter<<"\t"<<circleMarkers.back().getMarkerAngle()<<circleMarkers.back().getSupportCells()<<"->"<<2*Foam::constant::mathematical::pi<<circleMarkers.front().getSupportCells()<<" | middleAngle:"<<middleAngle<<"  ||||||||"<<circleMarkers.back().getSupportCells().size()<<"|"<<circleMarkers.front().getSupportCells().size()<<Foam::endl;
            LagrangianMarkerOnCrossSec middleMarker
                (mesh,rodNumber,oneRod,parameter,oneCrossSec,middleAngle,radiusFrac);
            auto inserted = circleMarkers.insert
                (
                    markersIter1,
                    middleMarker
                );
            refined=true;
        }
        //Info<<" angleRefinementCount:"<<refinementCount<<"  refined:"<<refined<<Foam::endl;
        refinementCount++;
        
        if(refinementCount==6)
        {
            LagrangianMarkerOnCrossSec marker1 = createLagrangianMarkerOnCrossSec(oneRod,0,oneCrossSec,6.18501);
            LagrangianMarkerOnCrossSec marker2 = createLagrangianMarkerOnCrossSec(oneRod,0,oneCrossSec,6.28319);
    
            scalar distOnNurbs = distance(marker1,marker2);
            scalar distDirect = Foam::mag(marker1.getSupportCells()-marker2.getSupportCells());
    
            Info<<"XXXxxxXXX "<<distOnNurbs<<"|"<<distDirect<<"/ "<<marker1.getSupportCells()<<"->"<<marker2.getSupportCells()<<Foam::endl;
            
            
            FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);
        }
    }
}
*/

/*
bool Foam::CrossSectionStructure::doSubdivision
(
    const std::list<LagrangianMarkerOnCrossSec>& smallerSide,
    const std::list<LagrangianMarkerOnCrossSec>& largerSide
)
{
    std::unordered_set<label> smallSideCompleteSupport;
    for(auto smIter=smallerSide.begin(); smIter!=smallerSide.end(); smIter++)
        smallSideCompleteSupport.insert(smIter->getSupportCells().begin(),smIter->getSupportCells().end());
    std::unordered_set<label> smallSideCells;
    for(auto smIter=smallerSide.begin(); smIter!=smallerSide.end(); smIter++)
        smallSideCells.insert(smIter->getMarkerCell());
    smallSideCells.erase(-1);
    
    std::unordered_set<label> largeSideCompleteSupport;
    for(auto laIter=largerSide.begin(); laIter!=largerSide.end(); laIter++)
        largeSideCompleteSupport.insert(laIter->getSupportCells().begin(),laIter->getSupportCells().end());
    std::unordered_set<label> largeSideCells;
    for(auto laIter=largerSide.begin(); laIter!=largerSide.end(); laIter++)
        largeSideCells.insert(laIter->getMarkerCell());
    largeSideCells.erase(-1);

    bool allSmallOverlapLarge = true;
    for(auto laIter=largeSideCells.begin(); laIter!=largeSideCells.end(); laIter++)
    {
        label laCell = *laIter;
        if(smallSideCompleteSupport.find(laCell)==smallSideCompleteSupport.end())
            allSmallOverlapLarge = false;
    }
    //Info<<"\t\t\t"<<smallerSide.size()<<"     allSmallOverlapLarge:"<<allSmallOverlapLarge<<"  smallSideCompleteSupportSize:"<<smallSideCompleteSupport.size()<<"  smallSideCellsSize:"<<smallSideCells.size()<<"  (p:"<<smallerSide.front().getMarkerParameter()<<",r:"<<smallerSide.front().getMarkerRadiusFrac()<<")"<<Foam::endl;
    
    bool allLargeOverlapSmall = true;
    for(auto smIter=smallSideCells.begin(); smIter!=smallSideCells.end(); smIter++)
    {
        label smCell = *smIter;
        if(largeSideCompleteSupport.find(smCell)==largeSideCompleteSupport.end())
            allLargeOverlapSmall = false;
    }
    //Info<<"\t\t\t"<<largerSide.size()<<"     allLargeOverlapSmall:"<<allLargeOverlapSmall<<"  largeSideCompleteSupportSize:"<<largeSideCompleteSupport.size()<<"  largeSideCellsSize:"<<largeSideCells.size()<<"  (p:"<<largerSide.front().getMarkerParameter()<<",r:"<<largerSide.front().getMarkerRadiusFrac()<<")"<<Foam::endl;
        
    if(allSmallOverlapLarge && allLargeOverlapSmall)
        return false;
    if(allSmallOverlapLarge && largeSideCompleteSupport.size()==0)
        return false;
    if(smallSideCompleteSupport.size()==0 && allLargeOverlapSmall)
        return false;
    if(smallSideCompleteSupport.size()==0 && largeSideCompleteSupport.size()==0)
        return false;
    
    return true;
}

bool Foam::CrossSectionStructure::doSubdivisionCircumferential
(
    const LagrangianMarkerOnCrossSec& smallerSide,
    const LagrangianMarkerOnCrossSec& largerSide
)
{
    bool supportDomainOverlap = false;
    for(auto iterSm=smallerSide.getSupportCells().begin();
        iterSm!=smallerSide.getSupportCells().end(); iterSm++)
    {
        for(auto iterLa=largerSide.getSupportCells().begin();
            iterLa!=largerSide.getSupportCells().end(); iterLa++)
        {
            
            if(*iterSm == *iterLa)
                supportDomainOverlap = true;
        }
    }
    bool smallerSideHasSupp = smallerSide.getSupportCells().size()!=0;
    bool largerSideHasSupp = largerSide.getSupportCells().size()!=0;
    bool bothSidesHaveSupp = largerSideHasSupp&&smallerSideHasSupp;
    
    scalar distOnNurbs = distance(smallerSide,largerSide);
    scalar distDirect = Foam::mag(smallerSide.getMarkerPosition()-largerSide.getMarkerPosition());
    
    scalar smallerSideDomainMinSize = std::numeric_limits<scalar>::max();
    if(smallerSideHasSupp)
        smallerSideDomainMinSize = supportDomainMinSize(smallerSide.getSupportCells());
    
    scalar largerSideDomainMinSize = std::numeric_limits<scalar>::max();
    if(largerSideHasSupp)
        largerSideDomainMinSize = supportDomainMinSize(largerSide.getSupportCells());
    
    scalar overallDomainMinSide = std::min(smallerSideDomainMinSize,largerSideDomainMinSize);
    
    //Info<<" doSubdivision:"<<distOnNurbs<<"|"<<distDirect<<"/"<<overallDomainMinSide<<smallerSide.getSupportCells()<<"->"<<largerSide.getSupportCells()<<" -- "<<distance(smallerSide,largerSide)<<Foam::endl;
    
    if(bothSidesHaveSupp)
    {
        if(distOnNurbs<overallDomainMinSide)
        {
            if(!supportDomainOverlap)
            {
                Info<<Foam::endl;
                std::cout<<"smallerSide.getMarkerParameter():"<<smallerSide.getMarkerParameter()<<std::endl;
                std::cout<<"largerSide.getMarkerParameter():"<<largerSide.getMarkerParameter()<<std::endl;
                Info<<"distOnNurbs:"<<distOnNurbs<<Foam::endl;
                Info<<"distDirect:"<<distDirect<<Foam::endl;
                Info<<"smallerSideDomainMinSize:"<<smallerSideDomainMinSize<<Foam::endl;
                Info<<"largerSideDomainMinSize:"<<largerSideDomainMinSize<<Foam::endl;
                Info<<"overallDomainMinSide:"<<overallDomainMinSide<<Foam::endl;
                Info<<"smallerSideHasSupp:"<<smallerSideHasSupp<<Foam::endl;
                Info<<"largerSideHasSupp:"<<largerSideHasSupp<<Foam::endl;
                FatalErrorInFunction<<"Dist matches but no overlap"<<exit(FatalError);
            }
            return false;
        }
        else
            return true;
    }
    else if(smallerSideHasSupp || largerSideHasSupp)
    {
        if(distOnNurbs<overallDomainMinSide)
            return false;
        else
            return true;
    }
    else
    {
        if(distOnNurbs<overallDomainMinSide)
            return false;
        else
            return true;        
    }
}
*/

/*
std::list<Foam::LagrangianMarkerOnCrossSec> Foam::CrossSectionStructure::createInitialCircumMarkers
(
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* oneRod,
    scalar parameter,
    const CrossSection* oneCrossSec,
    scalar radiusFrac
)
{
    std::list<LagrangianMarkerOnCrossSec> initialCircle;
    LagrangianMarkerOnCrossSec zeroAngleMarker
        (mesh,rodNumber,oneRod,parameter,oneCrossSec,0,radiusFrac);
    initialCircle.push_back(zeroAngleMarker);
    
    LagrangianMarkerOnCrossSec piAngleMarker
        (mesh,rodNumber,oneRod,parameter,oneCrossSec,Foam::constant::mathematical::pi,radiusFrac);
    initialCircle.push_back(piAngleMarker);
    
    return initialCircle;
}
*/
