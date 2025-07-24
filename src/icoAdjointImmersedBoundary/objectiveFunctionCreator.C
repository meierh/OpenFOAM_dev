#include "icoAdjointImmersedBoundary.H"

Foam::solvers::objectiveFunction Foam::solvers::createEmptyObjective()
{
    Foam::solvers::objectiveFunction obj;
    obj.empty = true;
    return obj;
}

Foam::solvers::objectiveFunction Foam::solvers::createTotalPressureLoss()
{
    Foam::solvers::objectiveFunction obj;
    
    // J = int_in (p+0.5*u²) dS - int_out (p+0.5*u²) dS
    
    /* 
     * dJdp = int_in 1 dS - int_out 1 dS
     * dJdu = int_in u dS - int_out u dS
     * dJdT = 0
     */
    
    obj.dJdp_Inlet = [](const icoAdjointVelocityInletBC& bc)
    {
        return Field<scalar>(bc.patch().size(),1);
    };
    obj.dJdp_Wall = [](const icoAdjointVelocityWallBC& bc)
    {
        return Field<scalar>(bc.patch().size(),0);
    };
    obj.dJdu_uOutlet = [](const icoAdjointVelocityOutletBC& bc)
    {
        const fvPatchField<vector>& u = bc.patch().lookupPatchField<volVectorField,vector>("U");
        return Field<vector>(-u);
    };
    obj.dJdu_pOutlet = [](const icoAdjointPressureOutletBC& bc)
    {
        const fvPatchField<vector>& u = bc.patch().lookupPatchField<volVectorField,vector>("U");
        return Field<vector>(-u);
    };
    obj.J = [](Foam::solvers::icoImmersedBoundary const& domain)
    {
        scalar J=0;
        const fvBoundaryMesh& domainBCs = domain.mesh.boundary();
        const volScalarField& p = domain.p;
        const volScalarField::Boundary& p_boundary = p.boundaryField();
        const volVectorField& u = domain.U;
        const volVectorField::Boundary& u_boundary = u.boundaryField();
        
        const label inletPatchInd = domainBCs.findIndex("inlet");
        scalar J_inlet = 0;
        scalar J_inlet_u = 0;
        scalar J_inlet_p = 0;
        if(inletPatchInd!=-1)
        {
            const fvPatch& inletPatch = domainBCs[inletPatchInd];
            const scalarField& inFaceMag = inletPatch.magSf();
            tmp<vectorField> inNormals = -1*inletPatch.nf();
            const fvPatchField<vector>& UInlet = u_boundary[inletPatchInd];

            List<vector> UInlet_data(UInlet.size(),Foam::zero());
            for(int i=0; i<UInlet.size(); i++)
                UInlet_data[i] = UInlet[i];
            //Pout<<"UInlet_data:"<<UInlet_data<<Foam::endl;

            scalarField u_minN_inlet = UInlet & inNormals.ref();
            scalarField abs_u = 0.5*(u_minN_inlet * u_minN_inlet);
            scalarField vol_abs_u = abs_u*inFaceMag;
            for(scalar val : vol_abs_u)
                J_inlet_u += val;
            const fvPatchField<scalar>& pInlet = p_boundary[inletPatchInd];

            List<scalar> pInlet_data(pInlet.size(),Foam::zero());
            for(int i=0; i<pInlet.size(); i++)
                pInlet_data[i] = pInlet[i];
            //Pout<<"pInlet_data:"<<pInlet_data<<Foam::endl;

            scalarField vol_p =  pInlet*inFaceMag;
            for(scalar val : vol_p)
                J_inlet_p += val;
            J_inlet = J_inlet_u+J_inlet_p;
        }
        Pout<<"Inlet J_u:"<<J_inlet_u<<Foam::nl;
        Pout<<"Inlet J_p:"<<J_inlet_p<<Foam::nl;
        Pout<<"Inlet J:"<<J_inlet<<Foam::nl;

        const label outletPatchInd = domainBCs.findIndex("outlet");
        scalar J_outlet = 0;
        scalar J_outlet_u = 0;
        scalar J_outlet_p = 0;
        if(outletPatchInd!=-1)
        {
            const fvPatch& outletPatch = domainBCs[outletPatchInd];
            const scalarField& outFaceMag = outletPatch.magSf();
            tmp<vectorField> outNormals = outletPatch.nf();
            const fvPatchField<vector>& UOutlet = u_boundary[outletPatchInd];

            List<vector> UOutlet_data(UOutlet.size(),Foam::zero());
            for(int i=0; i<UOutlet.size(); i++)
                UOutlet_data[i] = UOutlet[i];
            //Pout<<"UOutlet_data:"<<UOutlet_data<<Foam::endl;

            scalarField u_minN_outlet = UOutlet & outNormals.ref();
            scalarField abs_u = 0.5*(u_minN_outlet * u_minN_outlet);
            scalarField vol_abs_u = abs_u*outFaceMag;
            for(scalar val : vol_abs_u)
                J_outlet_u += val;
            const fvPatchField<scalar>& pOutlet = p_boundary[outletPatchInd];

            List<scalar> pOutlet_data(pOutlet.size(),Foam::zero());
            for(int i=0; i<pOutlet.size(); i++)
                pOutlet_data[i] = pOutlet[i];
            //Pout<<"pOutlet_data:"<<pOutlet_data<<Foam::endl;

            scalarField vol_p = pOutlet*outFaceMag;
            for(scalar val : vol_p)
                J_outlet_p -= val;
            J_outlet = J_outlet_u+J_outlet_p;
        }
        Pout<<"Outlet J_u:"<<J_outlet_u<<Foam::nl;
        Pout<<"Outlet J_p:"<<J_outlet_p<<Foam::nl;
        Pout<<"Outlet J:"<<J_outlet<<Foam::nl;
        
        J = J_inlet-J_outlet;
        Pstream::gather<scalar>(J,std::plus<scalar>());
        Pstream::scatter<scalar>(J);
        
        return J;
    };
    obj.empty = false;
    
    return obj;
}
