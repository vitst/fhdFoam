
{
    forAll(RND, i)
    {
        scalar fV = Foam::sqrt(mesh.V()[i]);
        RND[i] = rand.GaussNormal<vector>() / fV;
    }

    //tmp<volScalarField> func_c = Foam::sqrt( Foam::mag( (C + ccc) * (1-C) ) );
    //tmp<volScalarField> func_c = Foam::sqrt( Foam::mag( C + ccc ) );
    //tmp<volScalarField> func_c = Foam::sqrt( Foam::mag( C + ccc ) * ( (rho-c_s*M)/(M*(c_in-c_s)) - C) );
    tmp<volScalarField> func_c = Foam::sqrt( C * ( rho/(M*c_in) - C) );

    //Info<<func_c<<nl<<nl;
    //Info<< (rho-c_s*M)/(M*(c_in-c_s)) <<nl;

    //exit(0);

    RND = (sqrtDTinv * Zc0 * func_c) * RND;

    RND.correctBoundaryConditions();

    tmp<fvScalarMatrix> tCEqn
    (
        fvm::ddt(C) 
        +
        fvm::div(phi, C)
        - 
        fvm::laplacian(D_h, C)
        +
        fvc::div(RND)
    );

    fvScalarMatrix& CEqn = tCEqn.ref();

    CEqn.relax();

    double residual = solve( CEqn ).initialResidual();

    Info << "CEqn residual: "<< residual << endl;
}
