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
    Foam::RASModels::twophaseMixingLengthBerzi

Group
    grpRASTurbulence

Description
    Mixing length model for two-phase flow
    (only valid for 1D case with Y the wall-normal direction)

    References:
    \verbatim
    "Turbulence Locality and Granularlike Fluid Shear Viscosity in Collisional
    Suspensions."
    Berzi, D. & Fraccarollo, L. (2015);
    Physical Review Letters 115 (19), 194501.

    "Stresses and Drag in Turbulent Bed Load From Refractive Index-Matched
    Experiments."
    Ni, W-J & Capart, H. (2018).
    Geophysical Research Letters 45 (14), 7000–7009.
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        twophaseMixingLengthBerziCoeffs
        {
            Cmu         0.09;
            expoLM      1.0;
            alphaMaxLM  0.55;
            kappaLM     0.225;
        }
    \endverbatim

    This file is in the \dir turbulenceModels/RAS/twophaseMixingLengthBerzi
    subdirectory
    inside the \dir turbulenceModels/RAS subdirectory
    inside the \dir turbulenceModels subdirectory

SourceFiles
    twophaseMixingLengthBerzi.C

\*---------------------------------------------------------------------------*/

#ifndef twophaseMixingLengthBerzi_H
#define twophaseMixingLengthBerzi_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class twophaseMixingLengthBerzi Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class twophaseMixingLengthBerzi
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
    dimensionedScalar d_;


    // Fields
    volScalarField k_;
    volScalarField epsilon_;

    // Protected Member Functions

        virtual void correctNut();

        //- No copy construct and assignment
        twophaseMixingLengthBerzi(const twophaseMixingLengthBerzi&) = delete;

        //- No copy assignment
        void operator=(const twophaseMixingLengthBerzi&) = delete;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("twophaseMixingLengthBerzi");

    // Constructors

    //- Construct from components
    twophaseMixingLengthBerzi
    (
        const alphaField& beta,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& betafRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    );


    //- Destructor
    virtual ~twophaseMixingLengthBerzi()
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
    #include "twophaseMixingLengthBerzi.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
