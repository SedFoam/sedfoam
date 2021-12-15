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
    Foam::RASModels::twophaseMixingLengthConst

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
        twophaseMixingLengthConstCoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            alpha       0.52;
            betaOmega   0.072;
            alphak      0.5;
            alphaOmega  0.5;
        }
    \endverbatim

    This file is in the \dir turbulenceModels/RAS/twophaseMixingLengthConst
    subdirectory
    inside the \dir turbulenceModels/RAS subdirectory
    inside the \dir turbulenceModels subdirectory

SourceFiles
    twophaseMixingLengthConst.C

\*---------------------------------------------------------------------------*/

#ifndef twophaseMixingLengthConst_H
#define twophaseMixingLengthConst_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class twophaseMixingLengthConst Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class twophaseMixingLengthConst
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
protected:

    // Fields
    volScalarField k_;
    volScalarField epsilon_;


    // Protected Member Functions

        virtual void correctNut();

        //- No copy construct
        twophaseMixingLengthConst(const twophaseMixingLengthConst&) = delete;

        //- No copy assignment
        void operator=(const twophaseMixingLengthConst&) = delete;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("twophaseMixingLengthConst");

    // Constructors

    //- Construct from components
    twophaseMixingLengthConst
    (
        const alphaField& beta,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    );


    //- Destructor
    virtual ~twophaseMixingLengthConst()
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
    #include "twophaseMixingLengthConst.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
