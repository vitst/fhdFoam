// Solve the Momentum equation

tmp<fvVectorMatrix> tUEqn
(
    fvm::div(phi, U)
    -
    fvm::laplacian(nu, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

solve(UEqn == -fvc::grad(p));
