
{
    Info << "Steady-state flow solver solver" << endl;
    steadyStateControl equilibrate(mesh, "EQUILIBRATE");
    while ( equilibrate.loop() )
    {
        #include "UEqnEQ.H"        // Navier Stokes
        #include "pEqnEQ.H"
    }


    dictionary convDiffEQ = mesh.solutionDict().subDict("CONVDIFFEQ");
    // convergence criteria for steady state convection-diffusion
    double tolerance = 0;        
    if( !convDiffEQ.readIfPresent<double>("tolerance", tolerance) ){
        SeriousErrorIn("main")
            << "No tolerance parameter in CONVDIFFEQ subdictionary"
            << exit(FatalError);
    }
    // maximum number of iterations
    int maxIter = 0;             
    if( !convDiffEQ.readIfPresent<int>("maxIter", maxIter) ){
        SeriousErrorIn("main")
            << "No maxIter parameter in CONVDIFFEQ subdictionary"
            << exit(FatalError);
    }

    Info << "Steady-state concentration solver"<< endl;

    int  iter = 0;
    while ( iter < maxIter )
    {
        tmp<fvScalarMatrix> tCEqn
        (
            fvm::div(phi, C) - fvm::laplacian(D, C)
        );
        fvScalarMatrix& CEqn = tCEqn.ref();
        CEqn.relax();
        double residual = solve( CEqn ).initialResidual();
  
        iter++;
        Info << " Step " << token::TAB << iter << token::TAB
             << " residual: "<< residual << " > " << tolerance << endl;
  
        if( residual < tolerance )
        {
            Info << "Convection-diffusion: "
                 << "ExecutionTime = " << runTime.elapsedCpuTime() << " s "
                 << "ClockTime = " << runTime.elapsedClockTime() << " s "
                 << nl << "Converged in " << iter << " steps.  Residual = "
                 << residual << nl << endl;
            break;                                      // Done
        }
        else if(iter >= maxIter)
        {
            Warning << "Runtime WARNING:"
                    << "Convection-diffusion solver did not converge." << nl
                    << "Maximum number of iterations"
                    << "  iter: "<< iter << exit(FatalError); // No convergence
        }
    }
}
