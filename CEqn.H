
{
    scalar preFactor = 0.0001 * Foam::sqrt( 2.0 * D.value() );

    scalar sqrtDTinv = 1.0 / Foam::sqrt(dt);
    forAll(RND, i)
    {
        //RND[i] = mesh.V()[i] / totVol * rand.GaussNormal<vector>();
        //RND[i] = totVol / mesh.V()[i] * rand.GaussNormal<vector>();
        //RND[i] = rand.GaussNormal<vector>();
        //RND[i] = 1.0 / std::sqrt(dt * mesh.V()[i]) * rand.GaussNormal<vector>();
        //Info<<"V: "<<mesh.V()[i]<<nl;
        //scalar fV = Foam::pow(mesh.V()[i], 1.0/3.0);
        scalar fV = Foam::sqrt(mesh.V()[i]);
        //Info<<fV<<nl;
        //RND[i] = rand.GaussNormal<vector>() / fV;
        RND[i] = rand.GaussNormal<vector>() / fV;
        //RND[i] = preFactor * Foam::sqrt(Foam::mag(C[i]*(1-C[i]))) * sqrtDTinv * rand.GaussNormal<vector>() / fV;
        //RND[i] = preFactor * sqrtDTinv * rand.GaussNormal<vector>() / fV;
    }
    //Info<<"min(RND): "<<min(RND)<<"  max(RND): "<<max(RND)<<nl;

    tmp<volScalarField> cc = Foam::mag(C * (1-C));
    RND = preFactor * Foam::sqrt( cc ) * sqrtDTinv * RND;

    /*
    vector maxR=vector::zero, minR=vector::one;
    forAll(RND, i)
    {
        if(mesh.C()[i].x()>2.5 && mesh.C()[i].x()<3.5)
        {
            minR = min( minR, RND[i] );
            maxR = max( maxR, RND[i] );
        }
    }
    */

    RND.correctBoundaryConditions();

    //Info<<"min(RND): "<<minR<<"  max(RND): "<<maxR<<"   pre:"<<preFactor*sqrtDTinv<<nl;
    Info<<"min(RND): "<<min(RND)<<"  max(RND): "<<max(RND)<<nl;
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

    //C.min(1.0);
    //C.max(0.0);
    //
    Info<<"min(C): "<<min(C)<<"  max(C): "<<max(C)<<nl;

    Info << "CEqn residual: "<< residual << endl;
}
