//  Create fields at runTime.timeName()
Info << "Reading field p" << endl;

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info << "Reading field U" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

# include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

Info << "Reading field C" << endl;

volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Creating random component field and surface curvature" << endl;
int numPatches = mesh.boundaryMesh().size();
wordList boundaryTypes(numPatches, "calculated"); // default bc

volVectorField RND
(
    IOobject
    (
        "RND",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "RND", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector(0,0,0)
    ),
    boundaryTypes
);


areaScalarField curv
(
    IOobject
    (
        "curvature",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fam,
    dimensionedScalar(
        "curvature", dimensionSet(0, -1, 0, 0, 0, 0, 0), 1.0
    )
);
curv = fam.faceCurvatures(); 

