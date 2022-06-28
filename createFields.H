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
setRefCell(p, piso.dict(), pRefCell, pRefValue);
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

volScalarField chi
(
    IOobject
    (
        "coverage",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(
        "cov", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0
    )
);



volScalarField probNuc
(
    IOobject
    (
        "probNuc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(
        "probNuc", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0
    )
);


volScalarField gamma_ef
(
    IOobject
    (
        "gamma_ef",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(
        "gamma_ef_dim", dimensionSet(0, 0, 0, 0, 0, 0, 0), gamma
    )
);

const pointField& fCs = mesh.Cf().boundaryField()[patchID];
forAll(fCs, i)
{
    if(fCs[i][0]>28 && fCs[i][0]<32 && fCs[i][1]>26 && fCs[i][1]<30)
    {
        gamma_ef.boundaryFieldRef()[patchID][i] = 0.01*gamma;
    }
    if(fCs[i][0]>28 && fCs[i][0]<32 && fCs[i][1]>-30 && fCs[i][1]<-26)
    {
        gamma_ef.boundaryFieldRef()[patchID][i] = 0.01*gamma;
    }
    if(fCs[i][0]>60 && fCs[i][0]<64 && fCs[i][1]>26 && fCs[i][1]<30)
    {
        gamma_ef.boundaryFieldRef()[patchID][i] = 0.01*gamma;
    }
    if(fCs[i][0]>60 && fCs[i][0]<64 && fCs[i][1]>-30 && fCs[i][1]<-26)
    {
        gamma_ef.boundaryFieldRef()[patchID][i] = 0.01*gamma;
    }

}




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

