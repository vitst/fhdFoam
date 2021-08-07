
{
    tmp<fvScalarMatrix> tCEqn
    (
        fvm::ddt(C) 
        +
        fvm::div(phi, C)
        - 
        fvm::laplacian(D_h, C)
    );

    fvScalarMatrix& CEqn = tCEqn.ref();

    CEqn.relax();

    double residual = solve( CEqn ).initialResidual();

    Info << "CEqn residual: "<< residual << endl;
}
