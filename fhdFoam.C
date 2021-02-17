/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-201X OpenFOAM Foundation
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
    fhdFoam

Description
    Solves diffusion equation for multispecies
 
    Works with official and extended versions of OpenFOAM
    Currently: OpenFOAM v2012

\*---------------------------------------------------------------------------*/

// common and simpleFoam
#include "fvCFD.H"
#include "faCFD.H"

// OF
//#include "pointPatchField.H"
//#include "syncTools.H"
#include "dynamicFvMesh.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"

// dissolFoam project
//#include "steadyStateControl.H"
//#include "dissolMotionPointPatchVectorField.H"
//#include "pointPatchField.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    //#include "createPimpleControl.H"
    pimpleControl pimple(mesh, "PIMPLE", false);

    // TODO make it general with input from a dictionary
    // create fa mesh using patch
    label patchID=mesh.boundaryMesh().findPatchID("reactive_surface");
    faMesh fam(mesh, patchID);
  
    // local includes
    #include "readDicts.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"

// * * * * *   MAIN LOOP   * * * * * * * * * * * * * * * * * * * * * //

    runTime.functionObjects().execute();     // Execute cntlDict functions

/*##########################################
 *   Time-dependent convection-diffusion solver
 *##########################################*/

    Info << "Time-dependent concentration solver"<< endl;

    Random rand(512462521);

    scalar totVol = 0;
    for(scalar cellVol: mesh.V()) totVol+=cellVol;
    Info<<"total volume: "<<totVol<<nl;

    Info<<"min max Mesh: "<< Foam::pow( min(mesh.V()), 1./3.)
                          << "  "
                          << Foam::pow( max(mesh.V()), 1./3.)
                          << nl;

    scalar dt = runTime.deltaTValue();

    C.correctBoundaryConditions();

    while ( runTime.run() )
    {
        //#include "readDyMControls.H"
        //#include "CourantNo.H"

        ++runTime;
        Info << "Begin cycle: Time = " << runTime.timeName() 
             << "    dt = " << dt
             << nl << endl;


/*###############################################
 *    Mesh motion & relaxation
 *    Control parameters in dynamicMeshDict
 *###############################################*/
        mesh.update();
        //mesh.controlledUpdate();
/*
        Info << "Mesh update: ExecutionTime = " 
             << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl<< endl;
*/
        // Info << " Update curvature"<<nl; 
        curv = fam.faceCurvatures();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())                                                
        {                                                                    

            if (mesh.changing())
            {
                //if (correctPhi)
                {
                    // Calculate absolute flux
                    // from the mapped surface velocity
                    phi = mesh.Sf() & Uf();

                    #include "correctPhi.H"

                    // Make the flux relative to the mesh motion
                    fvc::makeRelative(phi, U);
                }

                //if (checkMeshCourantNo)
                {
                //    #include "meshCourantNo.H"
                }
            }


            #include "UEqn.H"
            #include "CEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                // assuming pimple consistent is true
                #include "pcEqn.H"
            }
        }                                                                    

// *********************************************************
// *    Write Output data
// *********************************************************

        runTime.write();
        //runTime.printExecutionTime(Info);
    }

    Info << "End" << endl;
    return 0;
}

// ********************* End of the solver ************************** //
