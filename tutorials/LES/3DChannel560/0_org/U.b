/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5.0                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      binary;
    class       volVectorField;
    location    "0";
    object      U.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   #codeStream
{
    codeInclude
    #{
        #include "fvCFD.H"
        #include "Random.H"
    #};

    codeOptions
    #{
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude \
        -I$(LIB_SRC)/OpenFOAM/lnInclude
    #};
    codeLibs
    #{
        -lfiniteVolume \
        -lmeshTools\
        -lOpenFOAM
    #};
    code
    #{
        const IOdictionary& d = static_cast<const IOdictionary&>(dict);
	const fvMesh& mesh = refCast<const fvMesh>(d.db());
    vectorField Ub(mesh.nCells(),vector(0,0,0));
    /*
    IOobject Uheader
    (
        "U.b",
        "0",
        mesh,
        IOobject::MUST_READ
    );
    Info<< "Reading U.b" << endl;
    volVectorField U.b(Uheader, mesh);

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            "../constant/",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );
    dimensionedVector Ubar
    (
        transportProperties.lookup("Ubar")
    );
*/
    const scalar nu(1e-6);
    const vector Ubar(vector (0.5096,0,0));
    const scalar Retau(560);
    const scalar h(0.02);
    const bool setBulk(true);
    const bool perturb(true);
    
    Info<< "nu      = " << nu << endl;
    Info<< "Ubar    = " << Ubar << endl;
    Info<< "Re(tau) = " << Retau << endl;
    const scalar utau = Retau*nu/h;
    Info<< "u(tau)  = " << utau << endl;

    // Streamwise component of flow. 0=x, 1=y, 2=z
    //const scalar streamwise = 0;
    // Streamwise component of flow. 0=x, 1=y, 2=z
    //const scalar height = 1;
    // Spanwise component of flow. 0=x, 1=y, 2=z
    //const scalar spanwise = 2;
    const direction streamDir(0);
    const direction heightDir(1);
    const direction spanDir(2);

    //wall normal circulation
    const scalar duplus = Ubar[streamDir]*0.25/utau;
    //spanwise wavenumber: spacing z+ = 200
    const scalar betaPlus = 2.0*constant::mathematical::pi*(1.0/200.0);
    const scalar sigma = 0.00055;
    //streamwise wave number: spacing x+ = 500
    const scalar alphaPlus = 2.0*constant::mathematical::pi*(1.0/500.0);
    const scalar epsilon = Ubar[streamDir]/200.0;

    // Random number generator
    Random perturbation(1234567);

    const vectorField& centres = mesh.C();

    forAll(centres, celli)
    {
        // add a small (+/-20%) random component to enhance symetry breaking
        scalar deviation=1.0 + 0.2*perturbation.GaussNormal<scalar>();


        const vector& cCentre = centres[celli];

        scalar zplus = cCentre[spanDir]*Retau/h;
        scalar y = min(cCentre[heightDir], 2*h-cCentre[heightDir]);
        scalar yplus = y*Retau/h;
        scalar xplus = cCentre[streamDir]*Retau/h;

        if (setBulk)
        {
            // laminar parabolic profile
            Ub[celli] = vector::zero;

            Ub[celli][streamDir] =
                3.0*Ubar[streamDir] * (y/h - 0.5*sqr(y/h));
        }

        if (perturb)
        {
            // streak streamwise velocity
            Ub[celli][streamDir] +=
                (utau * duplus/2.0) * (yplus/40.0)
                * Foam::exp(-sigma * Foam::sqr(yplus) + 0.5)
                * Foam::cos(betaPlus*zplus)*deviation;

            // streak spanwise perturbation
            Ub[celli][spanDir] =
                epsilon
              * Foam::sin(alphaPlus*xplus)
              * yplus
              * Foam::exp(-sigma*Foam::sqr(yplus))
              * deviation;
        }
    }

    Info<< "Writing modified U field " << endl;
// OpenFoam.com version and Old OpenFoam.org (<= 6)
        Ub.writeEntry("", os);
// New OpenFoam.org version (>=7)
//        writeEntry(os, "", U);
    #};
};


boundaryField
{
    bottomWall
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    topWall
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    sides1_half0
    {
        type            cyclic;
    }
    sides1_half1
    {
        type            cyclic;
    }
    sides2_half0
    {
        type            cyclic;
    }
    sides2_half1
    {
        type            cyclic;
    }
    inout1_half0
    {
        type            cyclic;
    }
    inout1_half1
    {
        type            cyclic;
    }
    inout2_half0
    {
        type            cyclic;
    }
    inout2_half1
    {
        type            cyclic;
    }
}


// ************************************************************************* //
