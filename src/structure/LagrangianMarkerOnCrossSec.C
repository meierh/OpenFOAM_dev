#include "CrossSectionStructure.H"

Foam::LagrangianMarkerOnCrossSec::LagrangianMarkerOnCrossSec
(
    LineStructure& structure,
    const fvMesh& mesh,
    const label rodNumber,
    const ActiveRodMesh::rodCosserat* baseRod,
    const scalar markerParameter,
    CrossSection* baseCrossSec,
    const scalar markerAngle,
    const scalar radiusFrac
):
LagrangianMarker(structure,mesh,rodNumber,baseRod),
markerAngle(markerAngle),
markerRadiusFrac(radiusFrac),
baseCrossSec(baseCrossSec)
{
    this->markerParameter = markerParameter;
    evaluateMarker();
    if(markerRadiusFrac<1 && (markerParameter!=0 && markerParameter!=1))
    {
        Info<<"markerRadiusFrac:"<<markerRadiusFrac<<Foam::endl;
        Info<<"markerParameter:"<<markerParameter<<Foam::endl;
        Info<<"markerAngle:"<<markerAngle<<Foam::endl;
        Info<<to_string()<<Foam::endl;
        FatalErrorInFunction<<"Invalid marker"<<exit(FatalError);
    }
}

void Foam::LagrangianMarkerOnCrossSec::evaluateMarker()
{
    //auto t0 = std::chrono::system_clock::now();
    markerPosition = CrossSectionStructure::evaluateRodCircumPos
    (
        baseRod,markerParameter,baseCrossSec,markerAngle,markerRadiusFrac
    );
    //auto t1 = std::chrono::system_clock::now();
    markerCell = mesh.findCell(markerPosition);
    //auto t2 = std::chrono::system_clock::now();
    computeSupport();
    //auto t3 = std::chrono::system_clock::now();
    Pair<vector> h = minMaxNeighbourWidth(directSupport);
    //auto t4 = std::chrono::system_clock::now();
    h_plus = h.first();
    h_minus = h.second();
    dilation = dilationFactors(h);
    //auto t5 = std::chrono::system_clock::now();
    checkDirectSupport();
    //auto t6 = std::chrono::system_clock::now();
    reduceSupport();
    //auto t7 = std::chrono::system_clock::now();
    computeCharacLength();
    //auto t8 = std::chrono::system_clock::now();
    
    /*
    if(paraMarker.find(markerParameter)==paraMarker.end())
        paraMarker.insert(markerParameter);
    else
    {
        if(markerCell!=-1)
        {
            const List<List<Pair<label>>>& graph = structure.getMeshGraph();
            const List<Pair<label>>& cellNeigh = graph[markerCell];
            
            bool noBoundFace = true;
            for(const Pair<label>& faceNei : cellNeigh)
                if(faceNei.second()==-1)
                    noBoundFace = false;
            
            if(noBoundFace)
            {
                Info<<"directSupport:"<<directSupport<<Foam::nl;
                Info<<"fullSupport:"<<fullSupport<<Foam::nl;
                Info<<"cellNeigh:"<<cellNeigh<<Foam::nl;
                
                Info<<"t0-t1 evaluateRodCircumPos:"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count()<<Foam::nl;
                Info<<"t1-t2 findCell:"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()<<Foam::nl;
                Info<<"t2-t3 computeSupport:"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count()<<Foam::nl;
                Info<<"t3-t4 minMaxNeighbourWidth:"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count()<<Foam::nl;
                Info<<"t4-t5:"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count()<<Foam::nl;
                Info<<"t5-t6 checkDirectSupport:"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t5).count()<<Foam::nl;
                Info<<"t6-t7 reduceSupport:"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t7-t6).count()<<Foam::nl;
                Info<<"t7-t8 computeCharacLength:"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t8-t7).count()<<Foam::nl;
            
                Info<<"||t0-t8||"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t8-t0).count()<<Foam::nl;
                FatalErrorInFunction<<"Temp Stop"<<exit(FatalError);
            }
        }
    }
    */
}

Foam::vector Foam::LagrangianMarkerOnCrossSec::getMarkerVelocity() const
{
    return structure.evaluateRodVelocity(rodNumber,markerParameter,markerAngle,markerRadiusFrac);
}

void Foam::LagrangianMarkerOnCrossSec::computeCharacLength()
{
    // Ugly hack
    LineStructure& structure = const_cast<LineStructure&>(this->structure);

    //auto t0 = std::chrono::system_clock::now();

    scalar markerRadiusFrac = this->markerRadiusFrac;
    markerRadiusFrac = std::max(markerRadiusFrac,0.1);
    
    //auto t1 = std::chrono::system_clock::now();
    
    //Pair<vector> dMdp_dMdangle = CrossSectionStructure::derivateRodCircumPos(baseRod,markerParameter,baseCrossSec,markerAngle,markerRadiusFrac);
    Pair<vector> dMdp_dMdangle = structure.derivateRodCircumPos(rodNumber,markerParameter,markerAngle,markerRadiusFrac);
    vector dMdp = dMdp_dMdangle.first();
    vector dMdangle = dMdp_dMdangle.second();
    
    //auto t2 = std::chrono::system_clock::now();
    
    //Pair<vector> d2Mdp_d2Mdangle = CrossSectionStructure::derivate2RodCircumPos(baseRod,markerParameter,baseCrossSec,markerAngle,markerRadiusFrac);
    Pair<vector> d2Mdp_d2Mdangle = structure.derivate2RodCircumPos(rodNumber,markerParameter,markerAngle,markerRadiusFrac);
    vector d2Mdp = d2Mdp_d2Mdangle.first();
    vector d2Mdangle = d2Mdp_d2Mdangle.second();
    
    //auto t3 = std::chrono::system_clock::now();
    
    scalar len_dMdp = std::sqrt(dMdp & dMdp);
    vector cross_dMdp_d2Mdp = dMdp ^ d2Mdp;
    scalar curvature_dp = std::sqrt(cross_dMdp_d2Mdp & cross_dMdp_d2Mdp)/(len_dMdp*len_dMdp*len_dMdp);
    
    //auto t4 = std::chrono::system_clock::now();
    
    scalar len_dMdangle = std::sqrt(dMdangle & dMdangle);
    vector cross_dMdangle_d2Mdangle = dMdangle ^ d2Mdangle;
    scalar curvature_dangle = std::sqrt(cross_dMdangle_d2Mdangle & cross_dMdangle_d2Mdangle)/(len_dMdangle*len_dMdangle*len_dMdangle);
    
    //auto t5 = std::chrono::system_clock::now();
    
    scalar curvature = std::max(curvature_dangle,curvature_dp);
    
    if(curvature==0)
        FatalErrorInFunction<<"Curvature can not be zero!"<<exit(FatalError);
    
    markerCharLen = 1.0f/curvature;
    
    //auto t6 = std::chrono::system_clock::now();
    
    /*
    Info<<"--------------------------"<<Foam::nl;
    Info<<"\t t0-t1 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count()<<Foam::nl;
    Info<<"\t t1-t2 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()<<Foam::nl;
    Info<<"\t t2-t3 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count()<<Foam::nl;
    Info<<"\t t3-t4 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count()<<Foam::nl;
    Info<<"\t t4-t5 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count()<<Foam::nl;
    Info<<"\t t5-t6 :"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t5).count()<<Foam::nl;
    */
}
