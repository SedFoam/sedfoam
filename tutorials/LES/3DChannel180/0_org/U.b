/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5                                     |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "1";
    object      Ub;
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
        -I$(LIB_SRC)/meshTools/lnInclude
    #};
    codeLibs
    #{
        -lfiniteVolume \
        -lmeshTools
    #};
    code
    #{
        const IOdictionary& d = static_cast<const IOdictionary&>(dict);
	    const fvMesh& mesh = refCast<const fvMesh>(d.db());
        vectorField Ub(mesh.nCells());
        scalar Ubar = 1.;
        scalar h = 2.;
        scalar retau = 180;
        scalar pi = Foam::constant::mathematical::pi;
        forAll(mesh.C(), i)
        {
            scalar x = mesh.C()[i].x();
            scalar y = mesh.C()[i].y();
            scalar z = mesh.C()[i].z();
            Random perturbation(1234567);
            scalar deviation=1.0 + 0.5*perturbation.GaussNormal<scalar>();
            if (y >= h)
            {
                Ub[i] = vector(3*Ubar*(y/h-0.5*sqr(y/h))+4*Ubar*0.125*(y-2*h)*retau/h/40.0*exp(-0.0001*sqr((y-2*h)*retau/h)+0.5)*cos(2*pi/200*z*retau/h)*deviation, 0, 4*Ubar/200.0*sin(2*pi/500*x*retau/h)*(y-2*h)*retau/h*exp(-0.0002*sqr((y-2*h)*retau/h))*deviation);
            }
            else
            {
                Ub[i] = vector(3*Ubar*(y/h-0.5*sqr(y/h))+4*Ubar*0.125*y*retau/h/40.0*exp(-0.0001*sqr(y*retau/h)+0.5)*cos(2*pi/200*z*retau/h)*deviation, 0, 4*Ubar/200.0*sin(2*pi/500*x*retau/h)*y*retau/h*exp(-0.0002*sqr(y*retau/h))*deviation);
            }
        }
        Ub.writeEntry("", os);
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
    "(left|right)"
    {
        type            cyclic;
    }
    "(inlet|outlet)"
    {
        type            cyclic;
    }
}


// ************************************************************************* //
