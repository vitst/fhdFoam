// Solve the Momentum equation

tmp<fvVectorMatrix> tUEqn
(
    fvm::div(phi, U)
    -
    fvm::laplacian(nu_h, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

solve(UEqn == -fvc::grad(p));
