/*---------------------------------------------------------------------------*\
Copyright (C) 2015 Cyrille Bonamy, Julien Chauchat, Tian-Jian Hsu
                   and contributors

License
    This file is part of SedFOAM.

    SedFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SedFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SedFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
if (this->runTime_.outputTime())
{
    //write the dissipation term of the specific dissipation
    volScalarField omegad
    (
       IOobject
       (
            "omegad",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       -betaOmega_*sqr(omega_)
    );

    //write the dissipation term of the specific dissipation
    volScalarField omegad2
    (
        IOobject
        (
            "omegad2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        - scalar(0.0708)
         *((scalar(1.0)+scalar(85.0)*XsiOmega)
         /(scalar(1.0)+scalar(100.0)*XsiOmega))*sqr(omega_)
    );

    //write the dissipation term of the specific dissipation
    volScalarField XsiOmega_
    (
        IOobject
        (
            "XsiOmega",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        XsiOmega
    );

    //write the dissipation term of the specific dissipation
    volScalarField cross
    (
        IOobject
        (
            "cross",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    // (min((CDkOmega*pos((1-alpha_)-0.92)),
    //  dimensionedScalar("test", dimensionSet(0, 0, -2, 0, 0), 1.0e-4)))
    // (min(CDkOmega*0.5*(1+tanh(-10*(alpha_-0.1))),
    //  dimensionedScalar("test", dimensionSet(0, 0, -2, 0, 0), 1.0e-2)))
        (CDkOmega)
    );

    //write the dissipation term of the specific dissipation
    volTensorField vort2
    (
        IOobject
        (
            "vort2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        skew(fvc::grad(U)) & skew(fvc::grad(U))
    );

    /*
    //write the total contribution of particle presence on the specific
    //dissipation
    volScalarField omegaTot
    (
        IOobject
        (
            "omegaTot",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    - ESD_*omega_
    + ESD2()*C3om_*KE2_*omega_
    + C4om_*KE4_*(ESD5_*nut_/k_)*omega_
    );

    //write the drag contribution of particle presence on the specific
    //dissipation
    volScalarField omegaDrag
    (
        IOobject
        (
            "omegaDrag",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ESD2()*C3om_*KE2_*omega_
    );

    //write the suspension contribution of particle presence on the specific
    //dissipation
    volScalarField omegaSUSP
    (
        IOobject
        (
            "omegaSUSP",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        C4om_*KE4_*(ESD5_*nut_/k_)*omega_
    );

    //write the drag contribution of particle presence on the specific
    //dissipation
    volScalarField omega1
    (
        IOobject
        (
            "omega1",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        - ESD_*omega_
    );
    */

    volVectorField gradomega
    (
        IOobject
        (
            "gradomega",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Gradomega
    );

    volVectorField gradk
    (
        IOobject
        (
            "gradk",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Gradk
    );

    this->runTime_.write();

}
