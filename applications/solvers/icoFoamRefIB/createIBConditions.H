dictionary cutCellBoundaryHeatFlux("IB_HeatFlux");
cutCellBoundaryHeatFlux.add("type","wallHeatFlux");
cutCellBoundaryHeatFlux.add("libs","(fieldFunctionObjects)");
cutCellBoundaryHeatFlux.add("patches","(cutCell)");

functionObjects::FSIFluidBoundaryConditions iB_heatFluxFromFluid
(
    "heatFlux",
    runTime,
    cutCellBoundaryHeatFlux    
);
