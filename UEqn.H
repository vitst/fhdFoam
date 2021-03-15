
// Generate Gauss noise GWN tensor field
forAll(Zmom, iz)
{
    tensor tmpZ = rand.GaussNormal<tensor>();
    Zmom[iz] = tmpZ + tmpZ.T();
}

//dimensionedScalar scl = dimensionedScalar("scl", dimensionSet(0, -2, 1, 0, 0, 0, 0), 1.0);

dimensionedScalar momFactor = Foam::sqrt(nu * 1.0);

//Info<<" momFactor: "<<momFactor<<nl;

Zmom = momFactor.value() * Zmom;

// Solve the Momentum equation

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U)
    +
    fvm::div(phi, U)
    -
    fvm::laplacian(nu, U)
    -
    fvc::div(Zmom)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
}
