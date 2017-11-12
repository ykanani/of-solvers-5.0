/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    pisoFoam

Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

	surfaceScalarField uf = fvc::flux(u);
	surfaceScalarField Uf = fvc::flux(U);

	fvVectorMatrix UEqn
	(
	    fvm::ddt(u) + fvm::div(Uf,u)
	   // +fvm::div(uf,U)
	    //+ u&fvc::grad(U) 
	    - fvm::laplacian(nu, u)
	);
	fvScalarMatrix pEqn
	(
	    fvm::ddt(p) 
            //+ fvm::div(Uf,p)
	   // +fvm::div(uf,U)
	    //+ u&fvc::grad(U) 
	    //- fvm::laplacian(nu, p)
	);


	volScalarField M1(UEqn.A());
	//Info << "M1 = "  << M1 << endl;
	Info << "M1 = "  << pEqn << endl;
	Info << "M1 = "  << pEqn+pEqn << endl;
//	Info << "M1 = "  << UEqn*2 << endl;

/*
	fvOptions.constrain(UEqn);

	if (piso.momentumPredictor())
	{
	    solve(UEqn == -fvc::grad(p) );

	    fvOptions.correct(U);
	}

	fvVectorMatrix UhatEqn
	(
	    fvm::ddt(Uhat) - U/D +fvm::Sp(1.0/D,Uhat)
	);

	UhatEqn.relax();
	UhatEqn.solve();

        volScalarField rAU(1.0/UEqn.A());
	volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
	surfaceScalarField phiHbyA
	(
	    "phiHbyA",
	    fvc::flux(HbyA)
	  + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
	);

	MRF.makeRelative(phiHbyA);

	adjustPhi(phiHbyA, U, p);

	// Update the pressure BCs to ensure flux consistency
	constrainPressure(p, U, phiHbyA, rAU, MRF);

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
	fvOptions.correct(U); 
		runTime.write();

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
		    << nl << endl;
	    }
*/
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
