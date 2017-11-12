{
//alphat = turbulence->nut()/Prt;
//alphat.correctBoundaryConditions();


//const volScalarField& alphat = time_.lookupObject<volVectorField>("alphat");
alphat = turbulence->nut()/Prt;

/// const GeometricField<scalar, fvPatchField, volMesh>& alphat2 =
///        U.db().objectRegistry::
///        lookupObject<GeometricField<scalar, fvPatchField, volMesh> >
///    (
///             "alphat2"
///    );

//Info << "average(alphat)" << average(alphat) << endl;
volScalarField alphaEff2("alphaEff2", turbulence->nu()/Pr + alphat);



fvScalarMatrix T2Eqn
(
fvm::ddt(T2)
+ fvm::div(phi, T2)
- fvm::laplacian(alphaEff2, T2)
);

T2Eqn.relax();
T2Eqn.solve();
}
