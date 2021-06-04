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
    Foam::granularRheologyModel

Description
    Dense granular-flow rheology model for particle phase stress
    
    References:
    \verbatim
        GDRmidi (2004).
        On dense granular flows.
        The European Physical Journal E, 14, 341-365

        Forterre, Y. and Pouliquen, O. (2008).
        Flows of Dense Granular Media.
        Annual Review of Fluid Mechanics, 40, 1-24
    \endverbatim

    This file is in the \dir granularRheologyModels/granularRheologyModel 
    subdirectory
    inside the \dir granularRheologyModels subdirectory

SourceFiles
    granularRheologyModel.C

\*---------------------------------------------------------------------------*/

#ifndef granularRheologyModel_H
#define granularRheologyModel_H

#include "dragModel.H"
#include "phaseModel.H"
#include "autoPtr.H"
#include "FrictionModel.H"
#include "PPressureModel.H"
#include "FluidViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class granularRheologyModel Declaration
\*---------------------------------------------------------------------------*/

class granularRheologyModel
{
    // Private data

    const volScalarField& alpha_;
    const surfaceScalarField& phia_;
    volScalarField pa_new_value;
    const dimensionedScalar& rhoa_;
    const dimensionedScalar& da_;

    const dimensionedScalar& rhob_;
    const dimensionedScalar& nub_;

    //- dictionary holding the modeling info
    IOdictionary granularRheologyProperties_;

    //- use granular Rheology or not.
    Switch granularRheology_;

    //- use hydrostatic particulate pressure.
    Switch granularHydrostaticPressure_;

    //- use granular mixture model.
    Switch granularMixtureModel_;

    //- use dilatancy model or not.
    Switch granularDilatancy_;
    
    autoPtr<granularRheologyModels::FrictionModel> FrictionModel_;
    autoPtr<granularRheologyModels::PPressureModel> PPressureModel_;
    autoPtr<granularRheologyModels::FluidViscosityModel> FluidViscosityModel_;

    //- maximum packing
    const dimensionedScalar alphaMaxG_;

    //- Static friction coefficient
    const dimensionedScalar mus_;

    //- Dynamic friction coefficient
    const dimensionedScalar mu2_;

    //- Empirical rheology coefficient
    const dimensionedScalar I0_;

    //- Empirical Dilatancy coefficient
    const dimensionedScalar Bphi_;

    //- Exponent of the Krieger-Dougherty viscosity
    const dimensionedScalar n_;

    //- Bulk viscosity amplification factor
    const dimensionedScalar BulkFactor_;
    
    // critical volume fraction
    const dimensionedScalar alpha_c_;

    // dilatancy factor
    const dimensionedScalar K_dila_;

    //- Relaxation factor for pa
    const dimensionedScalar relaxPa_;

    //- Minimum particle pressure (regularisation)
    const dimensionedScalar PaMin_;

    //- Minimum relaxation timescale
    const dimensionedScalar tau_inv_min_;
    
    //- The friction coefficient
    volScalarField muI_;

    ////- The  dilatancy angle
    volScalarField delta_;

    //- Total particle  pressure
    volScalarField p_p_total_;
    
    //- The granular viscosity
    volScalarField mua_;

    //- The granular bulk viscosity
    volScalarField lambda_;

    //- The granular pressure
    volScalarField pa_;

    //- Fluid Effective viscosity (viscous contribution)
    volScalarField nuvb_;
    
    //-Viscous number
    volScalarField I_;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    granularRheologyModel(const granularRheologyModel&);

    //- Disallow default bitwise assignment
    void operator=(const granularRheologyModel&);


public:

    // Constructors

    //- Construct from components
    granularRheologyModel
    (
        const phaseModel& phasea,
        const phaseModel& phaseb,
        const volScalarField& p_p_total,
        const dimensionedScalar& Dsmall
    );


    //- Destructor
    virtual ~granularRheologyModel();


    // Member Functions

    void solve
    (
        const volScalarField& magD,
        const volScalarField& pf,
        const dimensionedScalar& alphaSmall,
        const dimensionedScalar& Dsmall
    );

    bool on() const
    {
        return granularRheology_;
    }

    bool mixture() const
    {
        return granularMixtureModel_;
    }
    bool dilatancy() const
    {
        return granularDilatancy_;
    }

    const volScalarField& mua() const
    {
        return mua_;
    }

    const volScalarField& pa() const
    {
        return pa_;
    }

    const dimensionedScalar& alphaMaxG() const
    {
        return alphaMaxG_;
    }
    const volScalarField& p_p_total() const
    {
        return p_p_total_;
    }
    const volScalarField& I() const
    {
        return I_;
    }
    const volScalarField& muI() const
    {
        return muI_;
    }
    const volScalarField& delta() const
    {
        return delta_;
    }
    const volScalarField& lambda() const
    {
        return lambda_;
    }
    const volScalarField& nuvb() const
    {
        return nuvb_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
