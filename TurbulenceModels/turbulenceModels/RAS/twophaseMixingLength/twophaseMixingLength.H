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

Class
    Foam::RASModels::twophaseMixingLength

Group
    grpRASTurbulence

Description
    Mixing length model for two-phase flow
    (only valid for 1D case with Y the wall-normal direction)

    References:
    \verbatim
        Revil-Baudard, T. and Chauchat, J. (2013)
        A two-phase model for sheet flow regime based on dense granular flow
        rheology.
        Journal of Geophysical Research - Oceans, 118, 619-634
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        twophaseMixingLengthCoeffs
        {
            Cmu         0.09;
            expoLM      1.0;
            alphaMaxLM  0.55;
            kappaLM     0.225;
        }
    \endverbatim

    This file is in the \dir turbulenceModels/RAS/twophaseMixingLength
    subdirectory
    inside the \dir turbulenceModels/RAS subdirectory
    inside the \dir turbulenceModels subdirectory

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
protected:

    // Protected data


    // Model coefficients
    dimensionedScalar Cmu_;
    dimensionedScalar expoLM_;
    dimensionedScalar alphaMaxLM_;
    dimensionedScalar kappaLM_;


    // Fields
    volScalarField k_;
    volScalarField epsilon_;

    // Protected Member Functions

        virtual void correctNut();

        //- No copy construct and assignment
        twophaseMixingLength(const twophaseMixingLength&) = delete;

        //- No copy assignment
        void operator=(const twophaseMixingLength&) = delete;


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
    virtual ~twophaseMixingLength()
    {}


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
