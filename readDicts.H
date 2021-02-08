  Info << "Reading dissolFoam flags" << nl;
  
  Info << "Reading transportProperties" << nl;

  IOdictionary transportProperties
  (
      IOobject
      (
          "transportProperties",
          runTime.constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      )
  );

  const dimensionSet dimInvL2(0, -2, 0, 0, 0, 0, 0);
  const dimensionSet dimDiffusion(0, 2, -1, 0, 0, 0, 0);
/*
  dimensionedScalar nu
  (
      "nu",
      dimDiffusion,
      transportProperties
  );
*/
  dimensionedScalar D
  (
      "D",
      dimDiffusion,
      transportProperties
  );

  dimensionedScalar lR
  (
      "lR",
      dimLength,
      transportProperties
  );


  Info << "***************************************************" << nl;
//  Info << "transportProperties, nu:    " << nu                  << nl;
  Info << "transportProperties, D:     " << D                   << nl;
  Info << "transportProperties, lR:    " << lR                  << nl;
  Info << "***************************************************" << endl;

