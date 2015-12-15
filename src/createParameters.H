Info<< "\nReading userDefinedParameters\n" << endl;

IOdictionary userDefinedParameters(
    IOobject
    (
        "userDefinedParameters",
        runTime.constant(),
        defaultMesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info << "....permiticity of free space, eps0" << endl;
dimensionedScalar eps0
(
    "eps0",
    dimensionSet(-1, -3, 4, 0, 0, 2, 0),
    userDefinedParameters.lookup("eps0")
);
Info << "    eps0 = " << eps0 << endl;

Info << "....electric potential coefficient of input sine signal, epMax" << endl;
dimensionedScalar epMax
(
    "epMax",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    userDefinedParameters.lookup("epMax")
);
Info << "    epMax = " << epMax << endl;

Info << "....charge density coefficient of input sine signal, rhocMax" << endl;
dimensionedScalar rhocMax
(
    "rhocMax",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    userDefinedParameters.lookup("rhocMax")
);
Info << "    rhocMax = " << rhocMax << endl;

Info << "....frequency of input sine signal, omega" << endl;
dimensionedScalar omega
(
    "omega",
    dimensionSet(0, 0, -1, 0, 0, 0, 0),
    userDefinedParameters.lookup("omega")
);
Info << "    omega = " << omega << endl;


