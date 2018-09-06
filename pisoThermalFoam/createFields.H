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


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseThermalTransportModel laminarTransport(U, phi);

	autoPtr<incompressible::turbulenceModel> turbulence
	(
	 incompressible::turbulenceModel::New(U, phi, laminarTransport)
	);

#include "createTurbulenceStatFields.H"

#include "createMRF.H"
