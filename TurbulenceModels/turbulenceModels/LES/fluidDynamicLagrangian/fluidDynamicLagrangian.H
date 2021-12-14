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
    Foam::LESModels::fluidDynamicLagrangian

Group
    grpLESTurbulence

Description
    Dynamic SGS model with Lagrangian averaging for the fluid phase

    Reference:
    \verbatim
        Meneveau, C., Lund, T. S., & Cabot, W. H. (1996).
        A Lagrangian dynamic subgrid-scale model of turbulence.
        Journal of Fluid Mechanics, 319, 353-385.
    \endverbatim

SourceFiles
    fluidDynamicLagrangian.C

\*---------------------------------------------------------------------------*/

#ifndef fluidDynamicLagrangian_H
#define fluidDynamicLagrangian_H

#include "LESeddyViscosity.H"
#include "simpleFilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                       Class fluidDynamicLagrangian Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class fluidDynamicLagrangian
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        fluidDynamicLagrangian(const fluidDynamicLagrangian&);
        void operator=(const fluidDynamicLagrangian&);


protected:

    // Protected data

        volScalarField flmb_;
        volScalarField fmmb_;
        volScalarField Cyb_;

        dimensionedScalar theta_;
        dimensionedScalar D_;
        dimensionedScalar nutMax_;
        dimensionedScalar CybMax_;

        simpleFilter simpleFilter_;
        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;

        dimensionedScalar flmb0_;
        dimensionedScalar fmmb0_;


    // Protected Member Functions

        //- Update sub-grid eddy-viscosity
        void correctNut(const tmp<volTensorField>& gradU);
        
        volScalarField fDelta
        (
            const tmp<volScalarField>& delta,
            const tmp<volScalarField>& alpha,
            const tmp<volVectorField>& Ua,
            const tmp<volVectorField>& Ub,
            const tmp<volScalarField>& rho,
            const tmp<volScalarField>& K
        );

        volScalarField hAlpha(const tmp<volScalarField>& alpha);
        
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("fluidDynamicLagrangian");


    // Constructors

        //- Construct from components
        fluidDynamicLagrangian
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
    virtual ~fluidDynamicLagrangian()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy of phase a
        tmp<volScalarField> k(const tmp<volTensorField>& gradU) const
        {
            return pow(2.0*flmb_/fmmb_, 2.0/3.0)
                   * pow(this->Ce_, -2.0/3.0)
                   * sqr(this->delta())*magSqr(dev(symm(gradU)));
        }

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(this->U_));
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", this->nut_ + this->nu())
            );
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
        
        //- Return spherical part of the SGS stress tensor
        virtual tmp<volScalarField> spherSigmaSGS();

        virtual tmp<volScalarField> sqrDeltaD();

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fluidDynamicLagrangian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
