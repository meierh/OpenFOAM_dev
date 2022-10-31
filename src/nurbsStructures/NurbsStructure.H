#ifndef NURBSSTRUCTURE_H
#define NURBSSTRUCTURE_H

#include <stdio.h>
#include <iostream>
//#include <direct.h>
#include <random>
#include <chrono>

#include "activeRodMesh.h"
#include "rodOpt.h"
#include "rodTools.h"
#include "rodLatticePBC.h"
#include "rodImport.h"
#include "rodCScircle.h"
#include "rodCSrectangle.h"
#include "timing.h"

#include <dirent.h>

#include "cutCellFvMesh.H"

namespace Foam
{
class NurbsStructure
{
    public:
        explicit NurbsStructure(Time& runTime);
        
    private:
        Time& runTime;
        fileName runDirectory;
        fileName caseName;
        fileName xmlPath;
        word getXMLPath();
        
        
        //******************* Parameters *******************//
            const T		latScale = 1.;			// Scaling factor for XML geometry
            const int	latDir = 2;				// Lattice direction (e.g. for gravity load)
        
            // * Parameters for IGA discretization
            const bool	use_mixed = 1;			// Use mixed method?
            const int	p_sim = 6;				// Refine degree of splines to p
            const int	el_sim = 16;			// Refine knot vector of splines to el elements

            // * Parameters for cross-sections
            const T		geoR0 = 0.8;			// Radius [mm]
            const T		geoE0 = 60;				// Young's modulus [MPa]
            const T		geoNu = 0.45;			// Poisson's ratio [.]
            const T		geoRho = 1100e-9;		// Density [kg/mm設
            const T		Temp0 = 25.;			// Ambient temperature [蚓]
            
            // * Parameters for loads and boundary conditions
            const int	bcMode = 0;				// 0:Clamp in latDir
            const bool	bcRot = 1;				// Fix rotations on clamped boundaries
            const bool	applyGravity = 1;		// Apply gravity load
            // Apply non-zero DBC along latDir
            const T		loadStrain = 0.1;		// ... equivalent to an effective strain of ...
            const int	loadDir = +1;			// ... in +1:tension or -1:compression
            const int	loadSteps = 10;			// ... in ... load steps
                
            // * Parameters for simulation control
            const int	omp_nthreads = 16;		// Use multithreading with OpenMP 
            const int	plot_vtk_geo = 1;		// Write VTK file with initial geometry
            const bool	plot_vtk_solid = 0;		// Write VTK files with volume geometry
            const int	plot_vtk = 15;			// Write VTK each ... load steps
            const int	plot_n = 11;			// No. of points along axis in VTK files
            const int	plot_c = 9;				// No. of points around cross-section in VTK files

            // * Parameters for contact
            const bool	isContact = 0;
            ContactOptions cntOpt;
            cntOpt.ptsType = 2;
            cntOpt.ptsN = 2;
            cntOpt.csFac = 2.0;
            cntOpt.preg = 0.02 * geoR0; 
            cntOpt.kc = 0.2 * geoE0*geoR0*geoR0; 
            cntOpt.initOut = 0;
        //******************* Parameters *******************//

};
}
#endif
