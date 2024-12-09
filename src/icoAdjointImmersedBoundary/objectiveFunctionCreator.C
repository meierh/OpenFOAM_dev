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
    obj.J = [](Foam::solvers::icoAdjointImmersedBoundary const& domain)
    {
        scalar J=0;
        const fvBoundaryMesh& domainBCs = domain.mesh.boundary();
        const volScalarField& p = domain.p;
        const volScalarField::Boundary& p_boundary = p.boundaryField();
        const volVectorField& u = domain.U;
        const volVectorField::Boundary& u_boundary = u.boundaryField();
        
        const label inletPatchInd = domainBCs.findIndex("inlet");
        if(inletPatchInd!=-1)
        {
            const fvPatch& inletPatch = domainBCs[inletPatchInd];
            const scalarField& inFaceMag = inletPatch.magSf();
            tmp<vectorField> inNormals = -1*inletPatch.nf();
            const fvPatchField<scalar>& pInlet = p_boundary[inletPatchInd];
            const fvPatchField<vector>& UInlet = u_boundary[inletPatchInd];
            scalarField u_minN_inlet = UInlet & inNormals.ref();
            scalarField abs_u = 0.5*(u_minN_inlet * u_minN_inlet);
            scalarField p_plus_abs_u = abs_u+pInlet;
            scalarField int_p_plus_abs_u = p_plus_abs_u*inFaceMag;
            for(scalar val : int_p_plus_abs_u)
                J += val;
        }

        const label outletPatchInd = domainBCs.findIndex("outlet");
        if(outletPatchInd!=-1)
        {
            const fvPatch& outletPatch = domainBCs[outletPatchInd];
            const scalarField& outFaceMag = outletPatch.magSf();
            tmp<vectorField> outNormals = outletPatch.nf();
            const fvPatchField<scalar>& pOutlet = p_boundary[outletPatchInd];
            const fvPatchField<vector>& UOutlet = u_boundary[outletPatchInd];
            scalarField u_minN_outlet = UOutlet & outNormals.ref();
            scalarField abs_u = 0.5*(u_minN_outlet * u_minN_outlet);
            scalarField p_plus_abs_u = abs_u+pOutlet;
            scalarField int_p_plus_abs_u = p_plus_abs_u*outFaceMag;
            for(scalar val : int_p_plus_abs_u)
                J -= val;
        }
        
        Pstream::gather<scalar>(J,std::plus<scalar>());
        Pstream::scatter<scalar>(J);
        
        return J;
    };
    obj.empty = false;
    
    return obj;
}
