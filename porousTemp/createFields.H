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
Info<< "Reading field Tf\n" << endl;
volScalarField Tf
(
    IOobject
    (
        "Tf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field Ts\n" << endl;
volScalarField Ts
(
    IOobject
    (
        "Ts",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field C0\n" << endl;
volScalarField C0
(
    IOobject
    (
        "C0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
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

Info<< "Reading D\n"<<endl;
dimensionedScalar D_
(
    transportProperties.lookup("D")
);
Info<< "Reading lp\n"<<endl;
dimensionedScalar lp_
(
    transportProperties.lookup("lp")
);
Info<< "Reading st\n"<<endl;
dimensionedScalar st_
(
    transportProperties.lookup("st")
);
Info<< "Reading Dh\n"<<endl;
dimensionedScalar Dh_
(
    transportProperties.lookup("Dh")
);
Info<< "Reading d\n"<<endl;
dimensionedScalar d_
(
    transportProperties.lookup("d")
);
Info<< "Reading sl\n"<<endl;
dimensionedScalar sl_
(
    transportProperties.lookup("sl")
);
Info<< "Reading fluid velocity\n"<<endl;
dimensionedScalar Z
(
    transportProperties.lookup("Z")
);
Info<< "Reading Pr\n"<<endl;
dimensionedScalar Pr
(
    transportProperties.lookup("Pr")
);
Info<< "Reading nu\n"<<endl;
dimensionedScalar nu
(
    transportProperties.lookup("nu")
);
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
