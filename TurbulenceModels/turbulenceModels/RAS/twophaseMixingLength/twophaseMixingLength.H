/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RASModels::twophaseMixingLength

Group
    grpRASTurbulence

Description
    Standard high Reynolds-number k-omega turbulence model for
    incompressible flows.

    References:
    \verbatim
        "Turbulence Modeling for CFD"
        D. C. Wilcox,
        DCW Industries, Inc., La Canada,
        California, 1988.

        See also:
        http://www.cfd-online.com/Wiki/Wilcox's_k-omega_model
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        twophaseMixingLengthCoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            alpha       0.52;
            betaOmega   0.072;
            alphak      0.5;
            alphaOmega  0.5;
        }
    \endverbatim

SourceFiles
    twophaseMixingLength.C

\*---------------------------------------------------------------------------*/

#ifndef twophaseMixingLength_H
#define twophaseMixingLength_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class twophaseMixingLength Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class twophaseMixingLength
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        twophaseMixingLength(const twophaseMixingLength&);
        void operator=(const twophaseMixingLength&);


protected:

    // Protected data

        // Dictionary

            IOdictionary twophaseRASProperties_;

        // Model coefficients
            dimensionedScalar expoLM_;
            dimensionedScalar alphaMaxLM_;
            dimensionedScalar kappaLM_;


        // Fields
            volScalarField k_;
            volScalarField epsilon_;

    // Protected Member Functions
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("twophaseMixingLength");

    // Constructors

        //- Construct from components
        twophaseMixingLength
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    #if (OFVERSION < 220 and OFVERSION >= 200) or OFVERSION == 00
        virtual ~twophaseMixingLength() = default;
    #else
        virtual ~twophaseMixingLength()
        {}
    #endif


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();


        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }
        virtual tmp<volScalarField> omega() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "omega",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    epsilon_/(k_),
                    epsilon_.boundaryField().types()
                )
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "twophaseMixingLength.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
