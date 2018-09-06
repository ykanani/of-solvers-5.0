Info << "Calculating Velocity-Pressure corrleation" << endl;
Up = U * p;
Info << "Calculating Velocity-Temperature corrleation" << endl;
UT = U*thermo.T();

Info << "Calculating Dissipation rates and SGS Renolds Stress" << endl;
volSymmTensorField S(symm(fvc::grad(U)));
ef = 2 * turbulence->nu()  *(S && S);           //filtered dissipation rate
esgs = 2 * turbulence->nut() *(S && S);         //rate of production of residual kinetic energy (SGS dissipation)
Ept = ef + esgs;                                //total epsilon
epsilon = turbulence->epsilon();

RSgs = turbulence->R();

