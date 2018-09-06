{
//alphat = turbulence->nut()/Prt;
//alphat.correctBoundaryConditions();


//const volScalarField& alphat = time_.lookupObject<volVectorField>("alphat");
alphat = turbulence->nut()/Prt;
/// const GeometricField<scalar, fvPatchField, volMesh>& alphat1 =
///        U.db().objectRegistry::
///        lookupObject<GeometricField<scalar, fvPatchField, volMesh> >
///    (
///             "alphat1"
///    );

Info <<  max(laminarTransport.alpha()) << endl;
Info <<  min(laminarTransport.alpha()) << endl;

Info <<  max(T1) << endl;
Info <<  min(T1) << endl;

laminarTransport.correct();
volScalarField alphaEff1("alphaEff1", laminarTransport.alpha() + alphat);



fvScalarMatrix T1Eqn
(
fvm::ddt(T1)
+ fvm::div(phi, T1)
- fvm::laplacian(alphaEff1, T1)
);

T1Eqn.relax();
T1Eqn.solve();
}
