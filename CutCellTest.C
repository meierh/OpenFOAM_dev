/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    OFtutorial11_modifyingTheMesh

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cutCellPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // set up the case
    #include "setRootCase.H"

    // create the run time object
    Info<< "Create time\n" << endl;
    Time runTime
    (
        Time::controlDictName,
        args.rootPath(),
        args.caseName()
    );

    // disable post-processing etc.
    runTime.functionObjects().off();

    
    Foam::Info
    << "Create mesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::cutCellPolyMesh basisMesh
    (
        Foam::IOobject
        (
            Foam::polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        ),
        [](const vector coord)->scalar
        {
            scalar center_x = 0;
            scalar center_y = 0;
            scalar radius = 1;
    
            scalar x = coord.x();
            scalar y = coord.y();
            scalar x_c = x-center_x;
            scalar y_c = y-center_y;
            return std::sqrt(x_c*x_c+y_c*y_c)-radius;
        }
    );
    
    basisMesh.printAddedPoints();
    basisMesh.printAddedEdges();
    basisMesh.printAddedFaces();

    // ---
    // Write the grid
    Info << nl << "Writing extruded mesh to time = " << runTime.timeName() << nl << endl;
    Info<< "nPoints: "<<basisMesh.nPoints() <<endl;
    Info<< "nInternalFaces: "<<basisMesh.nInternalFaces()<<endl;
    Info<< "nFaces: "<<basisMesh.nFaces()<<endl;
    Info<< "nCells: "<<basisMesh.nCells()<<endl;
    
    
    /*
    labelList Hello = basisMesh.faces()[0];
    for(int i=0;i<Hello.size();i++)
    {
        Info<<Hello[i]<<endl;
    }
    face thisFace = basisMesh.faces()[0];

    for(int i=0;i<Hello.size();i++)
    {
        Info<<"|"<<thisFace.prevLabel(i)<<"|"<<endl;
    }
    */
    /*
    labelList CutCells = getCutCells(basisMesh);
    pointField newPoints = addedPoints(basisMesh);
    
    volScalarField phi // note that pressure is a scalar field
    (
        IOobject
        (
            "phi", // name of the field
            runTime.timeName(), // name of the current time, i.e. the time folder to read from
            basisMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        basisMesh // initialises the field to match the size of the mesh with default (0) values
    );
    
    for (label cellI=0; cellI<basisMesh.C().size(); cellI++)
    {
        phi[cellI] = 0;
    }
    
    for(int i=0; i<CutCells.size();i++)
    {
        phi[CutCells[i]] = 1;
    }
    */
    
    /*
    addedMeshStructures newItems = addingMeshItems(basisMesh);
    cutCellData cutCells = createCutCells(basisMesh,newItems);
    modifyMesh(basisMesh,newItems,cutCells);
    
    pointField points = basisMesh.points();
    cellList cells = basisMesh.cells();
    faceList faces = basisMesh.faces();
    for(int i=0;i<cells.size();i++)
    {
        Info<<"Cell:"<<i<<" centre:"<<cells[i].centre(points,faces)<<endl;
        for(int k=0;k<cells[i].size();k++)
        {
            Info<<"\t";
            face oneFace = faces[cells[i][k]];
            for(int j=0;j<oneFace.size();j++)
            {
                Info<<points[oneFace[j]]<<"->";
            }
            Info<<endl;
        }
    }
    */
    
    basisMesh.write();
    runTime.loop();
    runTime.write();

    
    // ---
    Info << nl << "To best visualise the results, load the mesh and extract all patches" << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
