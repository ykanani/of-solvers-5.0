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

dimensionedScalar Pr(transportProperties.lookup("Pr"));
dimensionedScalar Prt(transportProperties.lookup("Prt"));

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

/*Usage

    div(phi,U)       Gauss localBlended linear upwind; // When blendingfactor is 1 use the first scheme. When blendingfactor is 0 use the second.

*/
//Info<< "Reading field UBlend\n" << endl;
volScalarField UBlend
(
    IOobject
    (
        "UBlend",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
  	IOobject::AUTO_WRITE

    ),
    mesh
);

//Info<< "Transforming field UBlend to surfaceScalarField named: UBlendingFactor\n" << endl;
surfaceScalarField UBlendingFactor
(
    IOobject
    (
        "UBlendingFactor",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
fvc::interpolate(UBlend)  // The interpolation. Please note: for the localBlended scheme to work on 'U' the variable UBlendingFactor should be present. It should be an IOobject and it should be a surfaceScalarField
);
// kinematic turbulent thermal thermal conductivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
 IOobject
 (
  "alphat",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
 ),
 mesh
 );

Info<< "Reading field T1\n" << endl;
volScalarField T1
(
 IOobject
 (
  "T1",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
 ),
 mesh
 );
Info<< "Reading field T2\n" << endl;
volScalarField T2
(
 IOobject
 (
  "T2",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
 ),
 mesh
 );


Info<< "Creating field ef\n" << endl; //Pope, eq 13.122
volScalarField ef
(
 IOobject
 (
  "ef",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar
 (
  "zero",
  dimensionSet( 0, 2, -3, 0, 0, 0, 0), // Or Umean.dimensions()
  0
 )

 );

Info<< "Creating field esgs\n" << endl; //Pope, eq 13.123
volScalarField esgs
(
 IOobject
 (
  "esgs",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar
 (
  "zero",
  dimensionSet( 0, 2, -3, 0, 0, 0, 0), // Or Umean.dimensions()
  0
 )

 );


Info<< "Creating field Ept\n" << endl; // total dissipation, Ept = ef + Pr
volScalarField Ept
(
 IOobject
 (
  "Ept",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar
 (
  "zero",
  dimensionSet( 0, 2, -3, 0, 0, 0, 0), // Or Umean.dimensions()
  0
 )

 );



Info<< "Creating field epsilon\n" << endl; //from SGS model, should be equal to Pr
volScalarField epsilon
(
 IOobject
 (
  "epsilon",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar
 (
  "zero",
  dimensionSet( 0, 2, -3, 0, 0, 0, 0), // Or Umean.dimensions()
  0
 )

 );

Info<< "Creating field RSgs\n" << endl;
volSymmTensorField RSgs
(
 IOobject
 (
  "RSgs",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedSymmTensor
 (
  "zero",
  dimensionSet( 0, 2, -2, 0, 0, 0, 0), // Or Umean.dimensions()
  symmTensor (0,0,0,0,0,0)
 )

 );

Info<< "Creating field UMean\n" << endl;
volVectorField UMean
(
 IOobject
 (
  "UMean",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  U.dimensions(), // Or Umean.dimensions()
  vector (0,0,0)
 )
 );

Info<< "Creating field T1Mean\n" << endl;
volScalarField T1Mean
(
 IOobject
 (
  "T1Mean",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar
 (
  "zero",
  T1.dimensions(), // Or Umean.dimensions()
  0
 )
 );
Info<< "Creating field T2Mean\n" << endl;
volScalarField T2Mean
(
 IOobject
 (
  "T2Mean",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar
 (
  "zero",
  T2.dimensions(), // Or Umean.dimensions()
  0
 )
 );

Info<< "Creating field pMean\n" << endl;
volScalarField pMean
(
 IOobject
 (
  "pMean",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar
 (
  "zero",
  p.dimensions(), // Or Umean.dimensions()
  0
 )
 );
Info<< "Creating field UT1\n" << endl;
volVectorField UT1
(
 IOobject
 (
  "UT1",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimensionSet( 0, 1, -1, 1, 0, 0, 0), // Or Umean.dimensions()
  vector (0,0,0)
 )

 );
Info<< "Creating field UT1Mean\n" << endl;
volVectorField UT1Mean
(
 IOobject
 (
  "UT1Mean",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimensionSet( 0, 1, -1, 1, 0, 0, 0), // Or Umean.dimensions()
  vector (0,0,0)
 )

 );

Info<< "creating field ut1Mean\n" << endl;
volVectorField ut1Mean
(
 IOobject
 (
  "ut1Mean",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimensionSet( 0, 1, -1, 1, 0, 0, 0), // Or Umean.dimensions()
  vector (0,0,0)
 )
 );

Info<< "Creating field UT2\n" << endl;
volVectorField UT2
(
 IOobject
 (
  "UT2",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimensionSet( 0, 1, -1, 1, 0, 0, 0), // Or Umean.dimensions()
  vector (0,0,0)
 )

 );


Info<< "Creating field UT2Mean\n" << endl;
volVectorField UT2Mean
(
 IOobject
 (
  "UT2Mean",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimensionSet( 0, 1, -1, 1, 0, 0, 0), // Or Umean.dimensions()
  vector (0,0,0)
 )

 );

Info<< "creating field ut2Mean\n" << endl;
volVectorField ut2Mean
(
 IOobject
 (
  "ut2Mean",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimensionSet( 0, 1, -1, 1, 0, 0, 0), // Or Umean.dimensions()
  vector (0,0,0)
 )
 );

Info<< "creating field UP\n" << endl;
volVectorField UP
(
 IOobject
 (
  "UP",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimensionSet( 0, 3, -3, 0, 0, 0, 0), // Or Umean.dimensions()
  vector (0,0,0)
 )
 );


Info<< "creating field UPMean\n" << endl;
volVectorField UPMean
(
 IOobject
 (
  "UPMean",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimensionSet( 0, 3, -3, 0, 0, 0, 0), // Or Umean.dimensions()
  vector (0,0,0)
 )
 );


Info<< "creating field upMean\n" << endl;
volVectorField upMean
(
 IOobject
 (
  "upMean",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimensionSet( 0, 3, -3, 0, 0, 0, 0), // Or Umean.dimensions()
  vector (0,0,0)
 )
 );


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

	autoPtr<incompressible::turbulenceModel> turbulence
	(
	 incompressible::turbulenceModel::New(U, phi, laminarTransport)
	);

#include "createMRF.H"
