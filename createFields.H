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
        "RND", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero
    ),
    boundaryTypes
);

volTensorField Zmom
(
    IOobject
    (
        "Zmom",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "Zmom_dim", dimensionSet(0, 2, -2, 0, 0, 0, 0), tensor::zero
    ),
    boundaryTypes
);

// surface curvature field
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

// by default the growth is ON on all faces
// if nucleation is enabled, then the growth is OFF or read from file
scalar nucleationVal = 1.0;
if(nucleation)
{
    nucleationVal = 0.0;
}

volScalarField growthOn
(
    IOobject
    (
        "growthOn",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("growthOn", dimensionSet(0, 0, 0, 0, 0, 0, 0), nucleationVal)
);



// create time field to test local euler
volScalarField rDeltaT
(
    IOobject
    (
        "rDeltaT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless/dimTime, runTime.deltaTValue()),
    extrapolatedCalculatedFvPatchScalarField::typeName
);

