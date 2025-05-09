#include "SensitivityInteraction.H"

Foam::SensitivityInteraction::SensitivityInteraction
(
    const fvMesh& mesh,
    LineStructure& structure,
    const IOdictionary& structureDict,
    markerMeshType modusFieldToMarker,
    markerMeshType modusMarkerToField
):
FieldMarkerStructureInteraction(mesh,structure,structureDict,modusFieldToMarker,modusMarkerToField)
{
    Info<<"Completed SensitivityInteraction setup"<<Foam::endl;
}

void Foam::SensitivityInteraction::deltaFunctionParamGradientCheck(std::vector<Parameter> paraList)
{
    /*
    std::vector<scalar> epsilonList = {1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12};
    for(Parameter para : paraList)
    {
        ParameterVariation& variator = ParameterVariation::createParameterVariator(&structure,para);
        List<FixedList<DynamicList<DynamicList<FixedList<scalar,2>>>,2>> f_values(epsilonList.size());
        for(label epsInd=0; epsInd<epsilonList.size(); epsInd++)
        {
            const std::vector<scalar> signs = {-1,1};
            for(label signInd=0; signInd<signs.size(); signInd++)
            {
                DynamicList<DynamicList<FixedList<scalar,2>>>& f_eps_sign = f_values[epsInd][signInd];
                
                // Parameter -/+ eps
                variator.vary(signs[signInd]*eps);
                

                for(const LagrangianMarker* marker : structure.getCollectedMarkers())
                {
                    vector X = marker->getMarkerPosition();
                    //const FixedList<scalar,10>& b = marker->getCorrParaB();
                    f_eps_sign.append(DynamicList<FixedList<scalar,2>>());
                    for(const Pair<label>& suppCell : marker->getSupportCells(LagrangianMarker::SupportType::Full))
                    {
                        vector cellCentre;
                        scalar cellVolume;
                        marker->getCellData(suppCell,cellCentre,cellVolume);
                        vector x = cellCentre;
                        
                        scalar f_dd = marker->deltaDirac(X,x);
                        scalar f_cdd = marker->correctedDeltaDirac(X,x);
                        f_eps_sign.last().append({f_dd,f_cdd});
                    }
                }
                variator.reset();
            }
                        
                        
                        
                        
                        
                        vector dX_dParam = structure.dXdParam(marker,para);
                        vector dcorrectedDeltaDirac_dX = marker->dcorrectedDeltaDirac_dX(X,x);
                        vector ddeltaDirac_dX = marker->ddeltaDirac_dX(X,x);            
                        scalar ddeltaDirac_dParam = ddeltaDirac_dX & dX_dParam;
                        scalar dcorrectedDeltaDirac_dParam = dcorrectedDeltaDirac_dX & dX_dParam;
                        
                        DynamicList<scalar> error_dd;
                        DynamicList<scalar> error_cdd;                
                        for(scalar eps : epsilonList)
                        {                   

                        }
            
            
            
            
            
            
            
            
            

            variator.reset();

            //Parameter + eps
            variator.vary(eps);
            scalar f1_dd = marker->deltaDirac(X,x);
            scalar f1_cdd = marker->correctedDeltaDirac(X,x);
            variator.reset();

            scalar fd_dDeltaDiracdParam = (f1_dd-f0_dd)/(2*eps);
            scalar fd_dcorrectedDeltaDiracdParam = (f1_cdd-f0_cdd)/(2*eps);

            scalar abs_error_dd = ddeltaDirac_dParam-fd_dDeltaDiracdParam;
            error_dd.append(abs_error_dd);

            scalar abs_error_cdd = dcorrectedDeltaDirac_dParam-fd_dcorrectedDeltaDiracdParam;
            error_cdd.append(abs_error_cdd);
            
        }
        
        
        
        for(const LagrangianMarker* marker : structure.getCollectedMarkers())
        {
            vector X = marker->getMarkerPosition();
            const FixedList<scalar,10>& b = marker->getCorrParaB();
            for(const Pair<label>& suppCell : marker->getSupportCells(LagrangianMarker::SupportType::Full))
            {
                vector cellCentre;
                scalar cellVolume;
                marker->getCellData(suppCell,cellCentre,cellVolume);
                vector x = cellCentre;
                
                vector dX_dParam = structure.dXdParam(marker,para);
                vector dcorrectedDeltaDirac_dX = marker->dcorrectedDeltaDirac_dX(X,x);
                vector ddeltaDirac_dX = marker->ddeltaDirac_dX(X,x);            
                scalar ddeltaDirac_dParam = ddeltaDirac_dX & dX_dParam;
                scalar dcorrectedDeltaDirac_dParam = dcorrectedDeltaDirac_dX & dX_dParam;
                
                DynamicList<scalar> error_dd;
                DynamicList<scalar> error_cdd;                
                for(scalar eps : epsilonList)
                {                   

                }
                
                auto minIter_dd = std::min_element(error_dd.begin(),error_dd.end());
                scalar min_error_dd = *minIter_dd;
                scalar min_index_dd = std::distance(error_dd.begin(),minIter_dd);
                auto maxIter_dd = std::max_element(error_dd.begin(),error_dd.end());
                scalar max_error_dd = *maxIter_dd;
                scalar max_index_dd = std::distance(error_dd.begin(),maxIter_dd);
                
                auto minIter_cdd = std::min_element(error_cdd.begin(),error_cdd.end());
                scalar min_error_cdd = *minIter_cdd;
                scalar min_index_cdd = std::distance(error_cdd.begin(),minIter_cdd);
                auto maxIter_cdd = std::max_element(error_cdd.begin(),error_cdd.end());
                scalar max_error_cdd = *maxIter_cdd;
                scalar max_index_cdd = std::distance(error_cdd.begin(),maxIter_cdd);
                
                Info<<"error_dd:"<<error_dd<<Foam::endl;
                Info<<"error_cdd:"<<error_cdd<<Foam::endl;
                Info<<"ddeltaDirac_dParam:"<<ddeltaDirac_dParam<<Foam::endl;
                Info<<"dcorrectedDeltaDirac_dParam:"<<dcorrectedDeltaDirac_dParam<<Foam::endl;
                
                if(min_error_dd > 1e-6 || min_error_cdd > 1e-6)
                {
                    if(ddeltaDirac_dParam!=0 && (min_error_dd/ddeltaDirac_dParam)<1e-6)
                        return;
                    
                    if(dcorrectedDeltaDirac_dParam!=0 && (min_error_cdd/dcorrectedDeltaDirac_dParam)<1e-6)
                        return;
                    
                    Info<<" min error_dd("<<min_index_dd<<"):"<<min_error_dd<<" max error_dd("<<max_index_dd<<"):"<<max_error_dd<<Foam::endl;
                    Info<<" min error_cdd("<<min_index_cdd<<"):"<<min_error_cdd<<" max error_cdd("<<max_index_cdd<<"):"<<max_error_cdd<<Foam::endl;
                    
                    Info<<"ddeltaDirac_dParam:"<<ddeltaDirac_dParam<<" -- "<<(min_error_dd/ddeltaDirac_dParam)<<Foam::endl;
                    Info<<"dcorrectedDeltaDirac_dParam:"<<dcorrectedDeltaDirac_dParam<<" -- "<<(min_error_cdd/dcorrectedDeltaDirac_dParam)<<Foam::endl;
                    
                    Info<<"error_dd:"<<error_dd<<Foam::endl;
                    Info<<"error_cdd:"<<error_cdd<<Foam::endl;
                    Info<<Foam::endl;
                    
                    FatalErrorInFunction<<"Invalid Gradient"<<exit(FatalError);
                }
                
                FatalErrorInFunction<<"Temp stop"<<exit(FatalError);
            }
        }
    }
    */
}
