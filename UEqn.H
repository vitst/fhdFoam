
// Solve the Momentum equation with noise Zmom generated outside the loop

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U)
    +
    fvm::div(phi, U)
    -
    fvm::laplacian(nu_h, U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
}
