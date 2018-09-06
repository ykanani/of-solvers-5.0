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
#include "IOporosityModelList.H"
#include "tubeBank.H"
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



    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
   
    Info<< "\nStarting time loop\n" << endl;
    Info<< "\nStarting to work\n" << endl;
    volScalarField alpha("alpha",turbulence->nu()/Pr);
    //dimensionedScalar a = ((0.257*pow((D_/lp_), 0.371)*pow((st_/(st_- d_)*Z*Dh_/nu),0.634)*pow(((sl_- lp_)/Dh_),0.134)*pow(((st_- d_)/Dh_),-0.278))/Dh_);//correlation implemented
    dimensionedScalar a = ((0.257*pow((D_/lp_), 0.371)*pow((st_/(st_- d_)*Dh_/nu),0.634)*pow(((sl_- lp_)/Dh_),0.134)*pow(((st_- d_)/Dh_),-0.278))/Dh_);//correlation implemented

    Info << "a.dimensions = " << a.dimensions() << endl;
    Info << "alpha.dimensions = " << alpha.dimensions() << endl;
	Info <<"alpha ="<< alpha <<endl;
    Info<< "a= "<< a << endl;
	dimensionedScalar AOV= ((lp_-D_)*2+(3.1415926535897932384626433832795*D_))/(st_*sl_);
	Info<< "AOV= "<< AOV << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            #include "UEqn.H"
            { 
            //Info<< "\nStarting to run\n" << endl;
            //Info<< "C0= "<< C0 << endl;
            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }
        laminarTransport.correct();
        turbulence->correct();

        
        
        fvScalarMatrix TEqn
        (
         
         //(rho*Z*Cpf)/lambda*fvc::div(T)=fvm::grad(T)+(omega*a)/((Tf-Ts)*Cpf)
        fvm::ddt(Tf) + fvm::div(phi,Tf)==fvm::laplacian(alpha,Tf)-fvm::Sp(a*pow(mag(U.component(vector::Y)),0.634)*alpha*C0*AOV,Tf) + (a*pow(mag(U.component(vector::Y)),0.634)*alpha*C0*AOV*Ts)
         );
         TEqn.relax();
         TEqn.solve();
    }
   
    
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************ //
