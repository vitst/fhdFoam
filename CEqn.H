
{
    scalar sqrtDT = Foam::sqrt(dt);
    forAll(RND, i)
    {
        //RND[i] = mesh.V()[i] / totVol * rand.GaussNormal<vector>();
        //RND[i] = totVol / mesh.V()[i] * rand.GaussNormal<vector>();
        //RND[i] = rand.GaussNormal<vector>();
        //RND[i] = 1.0 / std::sqrt(dt * mesh.V()[i]) * rand.GaussNormal<vector>();
        //Info<<"V: "<<mesh.V()[i]<<nl;
        scalar fV = Foam::pow(mesh.V()[i], 1.0/3.0);
        //Info<<fV<<nl;
        RND[i] = rand.GaussNormal<vector>() / fV;
    }
    //Info<<"min(RND): "<<min(RND)<<"  max(RND): "<<max(RND)<<nl;

    tmp<volScalarField> cc = Foam::mag(C * (1-C));
    RND = preFactor * Foam::sqrt( cc ) * sqrtDT * RND;

    tmp<fvScalarMatrix> tCEqn
    (
        fvm::ddt(C) 
        +
        fvm::div(phi, C)
        - 
        fvm::laplacian(D, C)
        +
        fvc::div(RND)
    );
    fvScalarMatrix& CEqn = tCEqn.ref();

    CEqn.relax();

    double residual = solve( CEqn ).initialResidual();
    Info << " residual: "<< residual << endl;
}
