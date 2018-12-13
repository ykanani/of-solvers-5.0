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
  IOobject::NO_WRITE
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
  IOobject::NO_WRITE
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
  IOobject::NO_WRITE
 ),
 mesh,
 dimensionedSymmTensor
 (
  "zero",
  dimensionSet( 0, 2, -2, 0, 0, 0, 0), // Or Umean.dimensions()
  symmTensor (0,0,0,0,0,0)
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
  IOobject::NO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimVelocity*dimTemperature,
  Zero
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
  IOobject::NO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimVelocity*dimTemperature,
  Zero
 )


 );

Info<< "creating field Up\n" << endl;
volVectorField Up
(
 IOobject
 (
  "Up",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 mesh,
 dimensionedVector
 (
  "zero",
  dimPressure*dimVelocity/dimDensity,
  Zero
 )
 );


