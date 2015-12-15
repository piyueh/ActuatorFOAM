Info << "\nReading transport properties for region "
     << fluidMesh.name() << nl << endl;

IOdictionary fluidProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        fluidMesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info << "....kinematic viscosity, nu" << endl;
dimensionedScalar nuFluid
(
    "nuFluid",
    dimViscosity,
    fluidProperties.lookup("nu")
);
Info << "    nu = " << nuFluid << endl;

Info << "....relative permitivity, eps" << endl;
dimensionedScalar epsFluid
(
    "epsFluid",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    fluidProperties.lookup("eps")
);
Info << "    eps = " << epsFluid << endl;

Info << "....Debye length, lambda" << endl;
dimensionedScalar lambdaFluid
(
    "lambdaFluid",
    dimensionSet(0, 1, 0, 0, 0, 0, 0),
    fluidProperties.lookup("lambda")
);
Info << "    lambda = " << lambdaFluid << endl;



Info << "\nReading transport properties for region "
     << solidMesh.name() << nl << endl;

IOdictionary solidProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        solidMesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info << "....relative permitivity, eps" << endl;
dimensionedScalar epsSolid
(
    "epsSolid",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    solidProperties.lookup("eps")
);
Info << "    eps = " << epsSolid << endl;


Info << "....Debye length, lambda" << endl;
dimensionedScalar lambdaSolid
(
    "lambdaSolid",
    dimensionSet(0, 1, 0, 0, 0, 0, 0),
    solidProperties.lookup("lambda")
);
Info << "    lambda = " << lambdaSolid << endl;



Info << "\nReading initial field for region "
     << fluidMesh.name() << nl << endl;

Info << "....electric field, ep" << endl;
volScalarField epFluid
(
    IOobject
    (
        "ep",
        runTime.timeName(),
        fluidMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    fluidMesh
);

Info << "....charge density, rhoc" << endl;
volScalarField rhocFluid
(
    IOobject
    (
        "rhoc",
        runTime.timeName(),
        fluidMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    fluidMesh
);

Info << "....flow velocity, U" << endl;
volVectorField UFluid
(
    IOobject
    (
        "U",
        runTime.timeName(),
        fluidMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    fluidMesh
);

Info << "....flow pressure, p" << endl;
volScalarField pFluid
(
    IOobject
    (
        "p",
        runTime.timeName(),
        fluidMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    fluidMesh
);

Info << "....generate eps scalar field, epsField" << endl;
volScalarField epsFieldFluid
(
    IOobject
    (
        "epsField",
        runTime.timeName(),
        fluidMesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluidMesh,
    dimensionedScalar(
        "epsFieldScalar", dimensionSet(0, 0, 0, 0, 0, 0, 0), epsFluid.value())
);

Info << "....reading/calculating face flux field, phi" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        fluidMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(UFluid) & fluidMesh.Sf()
);

Info << "....setting reference pressure cell and value" << endl;

    label   pRefCell = 0;
    scalar  pRefValue = 0.0;

setRefCell(pFluid, 
        fluidMesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
fluidMesh.setFluxRequired(pFluid.name());



Info << "\nReading initial field for region "
     << solidMesh.name() << nl << endl;

Info << "....electric field, ep" << endl;
volScalarField epSolid
(
    IOobject
    (
        "ep",
        runTime.timeName(),
        solidMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    solidMesh
);

Info << "....charge density, rhoc" << endl;
volScalarField rhocSolid
(
    IOobject
    (
        "rhoc",
        runTime.timeName(),
        solidMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    solidMesh
);

Info << "....generate eps scalar field, epsField" << endl;
volScalarField epsFieldSolid
(
    IOobject
    (
        "epsField",
        runTime.timeName(),
        solidMesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    solidMesh,
    dimensionedScalar(
        "epsFieldScalar", dimensionSet(0, 0, 0, 0, 0, 0, 0), epsSolid.value())
);
