Info << "Calculating local rotation matrix..." << endl;



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

volVectorField Utn
(
IOobject
	(
	"Utn",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
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
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
	U.component(0)*1000.0
);


volTensorField RotTensor
(
IOobject
	(
	"RotTensor",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedTensor("0", dimless, Zero)
);

vectorList pp(pointsDict.lookup("points"));
scalarList slope(pointsDict.lookup("slope"));
scalarList s(pointsDict.lookup("s"));



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
	RotTensor[celli]= tensor(unitx&tvnorm,sign*mag(tvnorm^unitx),0,-sign*mag(tvnorm^unitx),unitx&tvnorm,0,0,0,1);
	Utn[celli]=RotTensor[celli]&U[celli];

	distMin1=100000;
	i1=0;
	distMin2=100000;
	i2=0;

}




