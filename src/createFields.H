Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

dimensionedScalar eps0
(
    "eps0",
    dimensionSet(-1, -3, 4, 0, 0, 2, 0),
    transportProperties.lookup("eps0")
);

dimensionedScalar epMax
(
    "epMax",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    transportProperties.lookup("epMax")
);

dimensionedScalar rhocMax
(
    "rhocMax",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    transportProperties.lookup("rhocMax")
);

dimensionedScalar omega
(
    "omega",
    dimensionSet(0, 0, -1, 0, 0, 0, 0),
    transportProperties.lookup("omega")
);



Info << "Reading electric potential ep\n" << endl;

volScalarField ep
(
    IOobject
    (
        "ep",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info << "Reading relative permitivity eps\n" << endl;

volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("eps", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0)
);



Info << "Reading field rhoc\n" << endl;

volScalarField rhoc
(
    IOobject
    (
        "rhoc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading Debye length lambda\n" << endl;

volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("lambda", dimensionSet(0, 1, 0, 0, 0, 0, 0), 1.0)
);


Info<< "Reading field p\n" << endl;
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


Info<< "Reading field U\n" << endl;
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


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
