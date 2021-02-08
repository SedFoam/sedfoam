/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::kineticTheoryModel

Description

SourceFiles
    kineticTheoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryModel_H
#define kineticTheoryModel_H

#include "dragModel.H"
#include "phaseModel.H"
#include "autoPtr.H"
#include "viscosityModel.H"
#include "conductivityModel.H"
#include "radialModel.H"
#include "granularPressureModel.H"
#include "frictionalStressModel.H"
#include "fvCFD.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class kineticTheoryModel Declaration
\*---------------------------------------------------------------------------*/

class kineticTheoryModel
{
    // Private data

        const phaseModel& phasea_;
        const volVectorField& Ua_;
        const volVectorField& Ub_;
        const volScalarField& alpha_;
        const surfaceScalarField& phia_;

        const dragModel& draga_;

        const dimensionedScalar& rhoa_;
        const dimensionedScalar& da_;
        const dimensionedScalar& nua_;

        //- dictionary holding the modeling info
        IOdictionary kineticTheoryProperties_;

        //- use kinetic theory or not.
        Switch kineticTheory_;

        //- use generation == dissipation
        Switch equilibrium_;

         //- use collisions
        Switch collisions_;

        //- use extended kinetic theory
        Switch extended_;

        //- use production limiter for Theta (to kb)
        Switch limitProduction_;

        autoPtr<kineticTheoryModels::viscosityModel> viscosityModel_;

        autoPtr<conductivityModel> conductivityModel_;

        autoPtr<radialModel> radialModel_;

        autoPtr<granularPressureModel> granularPressureModel_;

        autoPtr<frictionalStressModel> frictionalStressModel_;

        //- coefficient of restitution
        const dimensionedScalar e_;

        //- maximum packing
        const dimensionedScalar alphaMax_;

        //- min value for which the frictional stresses are zero
        const dimensionedScalar DiluteCut_;
        const dimensionedScalar ttzero;
        const dimensionedScalar ttone;
        const dimensionedScalar MaxTheta;

        //- angle of internal friction
        const dimensionedScalar phi_;

        //Kill transfert of temperature from fluid to particles
        const dimensionedScalar killJ1_;

        //- The granular energy/temperature
        volScalarField Theta_;

        // The sqrt granular energy/temperature
        //  volScalarField ThetaSqrt_;

        //- The granular viscosity
        volScalarField mua_;

        volScalarField muf_;

        //- The granular bulk viscosity
        volScalarField lambda_;

        //- The granular pressure
        volScalarField pa_;

        //- frictional stress
        volScalarField pf_;

    //- The granular pressure derivative
        volScalarField ppMagf_;

        //- The granular temperature conductivity
        volScalarField kappa_;

        //- The radial distribution function
        volScalarField gs0_;

    //- The derivative of the radial distribution function
        volScalarField gs0Prime_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        kineticTheoryModel(const kineticTheoryModel&);

        //- Disallow default bitwise assignment
        void operator=(const kineticTheoryModel&);


public:

    // Constructors

        //- Construct from components
        kineticTheoryModel
        (
            const phaseModel& phasea,
            const volVectorField& Ub,
            const dragModel& draga
        );


        //- Destructor
        virtual ~kineticTheoryModel();


    // Member Functions

        void solve
        (
            const volTensorField& gradUat,
            const volScalarField& kb,
            const volScalarField& epsilonb,
            const volScalarField& nuft,
            const dimensionedScalar& B,
            const dimensionedScalar& tt
        );

        void updateRheo
        (
            const volScalarField& kb,
            const volScalarField& epsilonb,
            const dimensionedScalar& B
        );

        bool on() const
        {
            return kineticTheory_;
        }

        const volScalarField& Theta() const
        {
            return Theta_;
        }

        const volScalarField& mua() const
        {
            return mua_;
        }

        const volScalarField& muf() const
        {
            return muf_;
        }

        const volScalarField& pa() const
        {
            return pa_;
        }

        const volScalarField& pf() const
        {
            return pf_;
        }

        const volScalarField& ppMagf() const
        {
            return ppMagf_;
        }

        const volScalarField& lambda() const
        {
            return lambda_;
        }

        const volScalarField& kappa() const
        {
            return kappa_;
        }

        const dimensionedScalar& phi() const
        {
            return phi_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
