
// Generate Gauss noise GWN tensor field

forAll(Zmom, iz)
{
    scalar fV = Foam::sqrt(mesh.V()[iz]);
    tensor tmpZ = rand.GaussNormal<tensor>() / fV;
    Zmom[iz] = tmpZ + tmpZ.T();
    //tmpZ -= tr(tmpZ) * tensor::one / 3.;
}

Zmom = (scaleMomentumNoise * sqrtDTinv * Zu0) * Zmom;

