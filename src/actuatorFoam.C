/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

# include "fvCFD.H"
# include "pisoControl.H"
# include "regionProperties.H"

# include "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField.H"
# include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    regionProperties    rp(runTime);

    if (rp["fluid"].size() > 1)
    {
        cerr << "This solver currently doesn't "
             << "support more than one fluid region!!" << endl;
        exit(EXIT_FAILURE);
    }

    if (rp["solid"].size() > 1)
    {
        cerr << "This solver currently doesn't "
             << "support more than one fluid region!!" << endl;
        exit(EXIT_FAILURE);
    }

    # include "createRegionMeshes.H"
    # include "createSolverControls.H"
    # include "createParameters.H"
    # include "createRegionFields.H"

    # include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
      
    Info<< "\nCalculating dimensionless electric potential, ep\n" << endl;

    scalar          initRes1 = 1.;

    int             maxI = 0;

    
    while ((initRes1 > 5e-5) && (maxI < 10000))
    {
        maxI ++;
        Info<< "Iteration = " << maxI << nl << endl;

        initRes1 = solve
        (
            fvm::laplacian(epFluid)
        ).initialResidual();

        solve
        (
            fvm::laplacian(epSolid)
        );
    }

    if (maxI == 1000000)
    {
        cerr << "Electric potential solver diverged!" << endl;
        exit(EXIT_FAILURE);
    }


    Info<< "\nCalculating dimensionless charge density, rhoc\n" << endl;

    initRes1 = 1.;
    maxI = 0;

    while ((initRes1 > 1e-7) && (maxI < 10000))
    {
        maxI ++;
        Info<< "Iteration = " << maxI << nl << endl;

        solve
        (
            fvm::laplacian(rhocSolid) 
        );

        initRes1 = solve
        (
            fvm::laplacian(rhocFluid) 
                == fvm::Sp(1. / (lambdaFluid * lambdaFluid * epsFluid), rhocFluid)
        ).initialResidual();

    }

    if (maxI == 1000000)
    {
        cerr << "Charge Density solver diverged!" << endl;
        exit(EXIT_FAILURE);
    }

    volVectorField EFluid
    (
        IOobject
        (
            "EFluid",
            runTime.timeName(),
            fluidMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        - fvc::grad(epFluid)
    );


    Info<< "\nCalculating time-independent force field\n" << endl;

    volVectorField fFluid
    (
        IOobject
        (
            "fFluid",
            runTime.timeName(),
            fluidMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rhocMax * epMax * rhocFluid * EFluid / rhoFluid
    );

    fvMesh  &mesh = fluidMesh;
    volScalarField      &p = pFluid;
    volVectorField      &U = UFluid;
    dimensionedScalar   &nu = nuFluid;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            Info << "Solve momentum predictor" << endl;
            solve
            (
                UEqn == - fvc::grad(p) + 
                    fFluid * std::pow(
                        std::sin(2 * M_PI * omega.value() * runTime.value()), 2)
            ); 
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());

            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn.H();
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                (fvc::interpolate(HbyA) & mesh.Sf())
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    
      
      

    return 0;
}

/*


    */
// ************************************************************************* //
