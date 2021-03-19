
// Generate Gauss noise GWN tensor field
forAll(Zmom, iz)
{
    scalar fV = Foam::sqrt(mesh.V()[iz]);
    tensor tmpZ = rand.GaussNormal<tensor>() / fV;
    Zmom[iz] = tmpZ + tmpZ.T();
}

Zmom = sqrtDTinv * Zu0 * Zmom;

// Solve the Momentum equation

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U)
    +
    fvm::div(phi, U)
    -
    fvm::laplacian(nu_h, U)
    -
    fvc::div(Zmom)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
}
