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
//    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //while (runTime.run())
    while (runTime.loop())
    {
        //#include "readTimeControls.H"
        //#include "CourantNo.H"
        //#include "setDeltaT.H"
	//runTime++;
	
	Info<< "Time = " << runTime.timeName() << nl << endl;
        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            #include "UEqn.H"

            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }


        laminarTransport.correct();
        turbulence->correct();
	#include "T1Eqn.H"
        #include "T2Eqn.H"

////calculating other statistics
        const GeometricField<vector, fvPatchField, volMesh>& UMean =
        U.db().objectRegistry::
        lookupObject<GeometricField<vector, fvPatchField, volMesh> >
        (
             "UMean"
        );
        const GeometricField<scalar, fvPatchField, volMesh>& T1Mean =
        U.db().objectRegistry::
        lookupObject<GeometricField<scalar, fvPatchField, volMesh> >
        (
             "T1Mean"
        );
        const GeometricField<scalar, fvPatchField, volMesh>& T2Mean =
        U.db().objectRegistry::
        lookupObject<GeometricField<scalar, fvPatchField, volMesh> >
        (
             "T2Mean"
        );

        const GeometricField<scalar, fvPatchField, volMesh>& pMean =
        U.db().objectRegistry::
        lookupObject<GeometricField<scalar, fvPatchField, volMesh> >
        (
             "pMean"
        );
        const GeometricField<vector, fvPatchField, volMesh>& UT1Mean =
        U.db().objectRegistry::
        lookupObject<GeometricField<vector, fvPatchField, volMesh> >
        (
             "UT1Mean"
        );
        const GeometricField<vector, fvPatchField, volMesh>& UT2Mean =
        U.db().objectRegistry::
        lookupObject<GeometricField<vector, fvPatchField, volMesh> >
        (
             "UT2Mean"
        );

        const GeometricField<vector, fvPatchField, volMesh>& UPMean =
        U.db().objectRegistry::
        lookupObject<GeometricField<vector, fvPatchField, volMesh> >
        (
             "UPMean"
                                       
	);



        Info << "Calculating Velocity-Pressure corrleation" << endl;
        UP = U * p;
        upMean = UPMean - UMean*pMean;

        Info << "Calculating Velocity-Temperature corrleation" << endl;
        UT1 = U * T1;
        UT2 = U * T2;
        ut1Mean = UT1Mean - UMean*T1Mean; 
        ut2Mean = UT2Mean - UMean*T2Mean; 

        Info << "Calculating Dissipation rates and SGS Renolds Stress" << endl;
        volSymmTensorField S(symm(fvc::grad(U)));
        ef = 2 * turbulence->nu()  *(S && S);           //filtered dissipation rate
        esgs = 2 * turbulence->nut() *(S && S);         //rate of production of residual kinetic energy (SGS dissipation)
        Ept = ef + esgs;                                //total epsilon
        epsilon = turbulence->epsilon();

        RSgs = turbulence->R();





        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
