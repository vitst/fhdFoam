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
 
    Works with official and extended OpenFOAM
    Currently: OpenFOAM v2012

\*---------------------------------------------------------------------------*/

#include <time.h>
#include <sys/types.h>

// common
#include "fvCFD.H"
#include "faCFD.H"

// OF
#include "mathematicalConstants.H"
#include "dynamicFvMesh.H"
#include "pisoControl.H"
#include "steadyStateControl.H"

static const double Avogadro = 6.02214179e+23; // mol^-1
static const double k_Boltzmann = 1.38064852e-23; // m^2 * kg * s^-2 * K^-1

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createControl.H"

    // TODO overview Info statements for parallel runs.

    // TODO make it general with input from a dictionary
    // create fa mesh using patch
    label patchID=mesh.boundaryMesh().findPatchID("reactive_surface");
    faMesh fam(mesh, patchID);
  
    // local includes
    #include "readDicts.H"
    #include "createFields.H"

    #include "initContinuityErrs.H"

// * * * * *   MAIN LOOP   * * * * * * * * * * * * * * * * * * * * * //

    runTime.functionObjects().execute();     // Execute cntlDict functions

// * * * * *   INITIATING THE FLOW   * * * * * * * * * * * * * * * * * * * * * //
    time_t timer;
    struct tm y2k = {0};
    double seconds;
  
    y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
    y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;
  
    time(&timer);  /* get current time; same as: timer = time(NULL)  */
  
    seconds = difftime(timer,mktime(&y2k));

    Random rand(int(seconds) + 2021*Pstream::myProcNo());
    //Random rand(512462521 + 2021*Pstream::myProcNo());

    scalar dt = runTime.deltaTValue();

    if(runTime.value()==0)
    { 
        // turn nucleation off for equilibration
        if(nucleation)
        {
            growthOn.boundaryFieldRef()[patchID] = 1.0;
        }
        #include "equilibrate.H"
        if(nucleation)
        {
            growthOn.boundaryFieldRef()[patchID] = 0.0;
        }
    }
    {
        #include "equilAfterRestart.H"
    }
    C.correctBoundaryConditions();
    runTime.writeNow();

    int resCMt = mesh.checkTopology(true);
    int resCMg = mesh.checkGeometry(true);
    Info<<"result of check mesh I run:  "<<resCMt<<"  "<<resCMg<<nl;

/*##########################################
 *   Time-dependent convection-diffusion solver
 *##########################################*/
   
    if(nucleation)
    {
        // update nucleation after interpolation
        forAll(growthOn.boundaryField()[patchID], i)
        {
            if(growthOn.boundaryField()[patchID][i]<0.5)
                growthOn.boundaryFieldRef()[patchID][i] = 0.0;
            else
                growthOn.boundaryFieldRef()[patchID][i] = 1.0;
        }
    }
    else
    {
        growthOn.boundaryFieldRef()[patchID] = 1.0;
    }

    Info << "Time-dependent concentration solver"<< endl;

    scalar totVol = 0;
    for(scalar cellVol: mesh.V()) totVol+=cellVol;
    Info<< nl <<"total volume: "<<totVol<<nl;

    Info<<"min max Mesh: "<< Foam::pow( min(mesh.V()), 1./3.)
                          << "  "
                          << Foam::pow( max(mesh.V()), 1./3.)
                          << nl << nl;


    while ( runTime.run() )
    {
        #include "CourantNo.H"

        if(CoNum>10.0)
        {
            SeriousErrorIn("main")
                <<"*** Error. Co Num is too high  "<<CoNum<<"  adjust the mesh."<<nl
                << exit(FatalError);
        }

        ++runTime;
        Info << "Begin cycle: Time = " << runTime.timeName() 
             << "    dt = " << dt
             << nl << endl;


/*###############################################
 *    Mesh motion & relaxation
 *    Control parameters in dynamicMeshDict
 *###############################################*/
        // update nucleation

        if(nucleation)
        {
            const scalarField& faceAreas = mesh.magSf().boundaryField()[patchID]; 
            const scalarField& faceCs = C.boundaryField()[patchID];
            const pointField& faceCentres = mesh.Cf().boundaryField()[patchID];
            //Info<<"min(area): "<<min(faceAreas)<<"  max(faceAreas): "<<max(faceAreas)<<nl;

            point sphCenter(50, 0, 0.05);
            forAll(growthOn.boundaryField()[patchID], i)
            {
                // apply only if the face hasn't been activated
                if(growthOn.boundaryField()[patchID][i] < 0.95)
                {
                    // calculate distance to center at (50 0 0.05)
                    scalar dist = Foam::mag( faceCentres[i] - sphCenter );
                    scalar surfL = faceAreas[i] / 0.8;
                    // here 1.0 is radius of the initial sphere
                    if(dist > surfL+1.0)
                    {
                        growthOn.boundaryFieldRef()[patchID][i] = 1;
                    }
                    else
                    {
                        scalar SI = Foam::log( (theta * theta * faceCs[i] * faceCs[i]).value() );
                        scalar sf = faceAreas[i];
    
                        scalar lnJ = lnA - Bcoef / (SI*SI);
                        scalar J = Foam::exp(lnJ);
    
                        scalar prob = 1 - Foam::exp( - (sf * J * dt * h0 * h0 * td).value() ); 
    
                        // here update
                        scalar tmpZ = rand.sample01<scalar>();
    
                        if(tmpZ < prob)
                        {
                            growthOn.boundaryFieldRef()[patchID][i] = 1;
                        }
                    }
                }
            }
            //Info<<"min(growthOn): "<<min(growthOn.boundaryField()[patchID])<<"  max(growthOn): "<<max(growthOn.boundaryField()[patchID])<<nl;
        }

        mesh.update();

        curv = fam.faceCurvatures();

        resCMg = mesh.checkGeometry(false);
        Info<<"Check mesh (geometry):  "<<resCMg<<nl;

        // Pressure-velocity PISO corrector
        {                                                                    
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }

        #include "CEqn.H"

// *********************************************************
// *    Write Output data
// *********************************************************

        runTime.write();
    }

    Info << "End" << endl;
    return 0;
}

// ********************* End of the solver ************************** //
