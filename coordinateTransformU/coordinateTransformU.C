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
    postChannel

Description
    Post-processes data from channel flow calculations.

    For each time: calculate: txx, txy,tyy, txy,
    eps, prod, vorticity, enstrophy and helicity. Assuming that the mesh
    is periodic in the x and z directions, collapse Umeanx, Umeany, txx,
    txy and tyy to a line and print them as standard output.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "OSspecific.H"
#include "vectorList.H"
#include "timeSelector.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions(true, true);

    #include "setRootCase.H"
    #include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"
    #include "readTransportProperties.H"

    const word& gFormat = runTime.graphFormat();

    IOdictionary pointsDict
                (
                 IOobject
                 (
                  "curve",
                  runTime.constant(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE
                 )
                );

    vectorList pp(pointsDict.lookup("points"));
    scalarList slope(pointsDict.lookup("slope"));
    scalarList s(pointsDict.lookup("s"));

    forAll(pp,i)
        {
                Info << "point " << i << " = " << pp[i] << endl;
                Info << "slope " << i << " = " << slope[i] << endl;
        }
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        // Reading U
	IOobject Uheader
                (
                 "U",
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ
                );
        volVectorField U(Uheader, mesh);

        Info << "reading fields..." << endl;
        volVectorField Utn
        (
        IOobject
        (
        "Utn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
        ),
        U
        );

        volScalarField wallDist
        (
        IOobject
        (
        "wallDist",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
        ),
        U.component(0)*1000.0
        );




        scalar x1(0);
        scalar y1(0);
        scalar x2(0);
        scalar y2(0);
        scalar xc(0),yc(0);
        scalar distMin1(100000);
        scalar i1(0);
        scalar distMin2(100000);
        scalar i2(0);
        scalar a,b,c,m;
        scalar dist(0);
        scalar xi,yi,mi,m1,m2;
	scalar sign(1);
        vector tv,tvnorm;
        vector unitx(1,0,0);
        tensor R(0,0,0,0,0,0,0,0,0);
        forAll( mesh.C(), celli)
        {

                xc = mesh.C()[celli].component(0);
                yc = mesh.C()[celli].component(1);

                //Info << "Cell  (" << xc << ","<< yc << "," <<  mesh.C()[celli].component(2) << endl;
                //Info << "U =" << U[celli] << endl;

                forAll(pp,i)
                {
                        dist = Foam::sqrt(pow((pp[i][0]-mesh.C()[celli].component(0)),2)+pow((pp[i][1]-mesh.C()[celli].component(1)),2));
                        if (dist < distMin1)
                        {
                                distMin2 = distMin1;
                                i2 = i1 ;

                                distMin1 = dist;
                                i1 = i ;
                        }

                        //Info << "slope " << i << " = " << slope[i] << endl;
                }
                if (i1==0)
                {
                        i2=1;
                }
                //Info << "Point1 " << pp[i1] << endl;
                //Info << "Min Dist1 =" << distMin1 << " to point :  " << i1 << endl;
                //Info << "Point2 " << pp[i2] << endl;
                //Info << "Min Dist2 =" << distMin2 << " to point :  " << i2 << endl;
                x1 = pp[i1][0];
                y1 = pp[i1][1];

                x2 = pp[i2][0];
                y2 = pp[i2][1];

                m1 = slope[i1];
                m2 = slope[i2];

                m = (y2-y1)/(x2-x1);
                b = 1;
                a = -m;
                c = -y1+m*x1;

                xi= (b*( b*xc -a*yc)-a*c)/(pow(a,2)+pow(b,2));
                yi= (a*(-b*xc +a*yc)-b*c)/(pow(a,2)+pow(b,2));
                wallDist[celli] = mag(a*xc+b*yc+c)/Foam::sqrt((pow(a,2)+pow(b,2)));
                //Info << "wallDist["<< celli << "]= " << wallDist[celli] << endl;
                mi = m1+ (xi-x1)*(m2-m1)/(x2-x1);

                tv = vector(1,mi,0);
                tvnorm = tv / mag(tv);
		if (mi < 0)
		{
			sign=-1;
		}
		else
		{
			sign=1;
		}
                R= tensor(unitx&tvnorm,sign*mag(tvnorm^unitx),0,-sign*mag(tvnorm^unitx),unitx&tvnorm,0,0,0,1);
                Utn[celli]=R&U[celli];

                distMin1=100000;
                i1=0;
                distMin2=100000;
                i2=0;
                //Info << "m1 =" << m1 << endl;
                //Info << "m2 =" << m2 << endl;
                //Info << "mi =" << mi << endl;

                //Info << "tvnorm =" << tvnorm << endl;
                //Info << "R =" << R << endl;

        }

        Info<< "\nWriting fields...\n" << endl;

         Utn.write();

        wallDist.write();

     }
        //Info << "tvnorm =" << tvnorm << endl;
        //vector vectemp(1,0,0);
        //vectemp=R&vector(0.707107 ,0.707107,0);
        //Info << "R =" << R << endl;
        //Info << "vectemp =" << vectemp << endl;


        Info<< "\nEnd\n" << endl;

    // For each time step read all fields
//    forAll(timeDirs, timeI)
//    {
//        runTime.setTime(timeDirs[timeI], timeI);
//        Info<< "Collapsing fields for time " << runTime.timeName() << endl;

  //      #include "readFields.H"
  //      #include "calculateFields.H"

        // Average fields over channel down to a line
    //    #include "collapse.H"
   // }




    return 0;
}


// ************************************************************************* //
