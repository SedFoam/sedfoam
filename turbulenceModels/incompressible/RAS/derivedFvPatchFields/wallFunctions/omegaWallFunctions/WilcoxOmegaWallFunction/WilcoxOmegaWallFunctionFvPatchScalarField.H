/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::incompressible::WilcoxOmegaWallFunctionFvPatchScalarField

Group
    grpIcoWallFunctions

Description
    This boundary condition provides a wall function constraint on turbulnce
    specific dissipation, omega.  The values are computed using:

        \f[
            \omega = sqrt(\omega_{vis}^2 + \omega_{log}^2)
        \f]

    where

    \vartable
        \omega_{vis} | omega in viscous region
        \omega_{log} | omega in logarithmic region
    \endvartable

    Model described by Eqs.(14, 15) of:
    \verbatim
        Menter, F., Esch, T.
        "Elements of Industrial Heat Transfer Prediction"
        16th Brazilian Congress of Mechanical Engineering (COBEM),
        Nov. 2001
    \endverbatim

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        Cmu          | model coefficient       | no          | 0.09
        kappa        | Von Karman constant     | no          | 0.41
        E            | model coefficient       | no          | 9.8
        beta1        | model coefficient       | no          | 0.075
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            WilcoxOmegaWallFunction;
    }
    \endverbatim

SourceFiles
    WilcoxOmegaWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef WilcoxOmegaWallFunctionFvPatchScalarField_H
#define WilcoxOmegaWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
             Class WilcoxOmegaWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class WilcoxOmegaWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

        //- Tolerance used in weighted calculations
        static scalar tolerance_;

         //- kn coefficient (roughness)
        scalar kn_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;

        //- beta1 coefficient
        scalar beta1_;

        //- Y+ at the edge of the laminar sublayer
        scalar yPlusLam_;

        //- Local copy of turbulence G field
        scalarField G_;

        //- Local copy of turbulence omega field
        scalarField omega_;

        //- Initialised flag
        bool initialised_;

        //- Master patch ID
        label master_;

        //- List of averaging corner weights
        List<List<scalar> > cornerWeights_;


    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;

        //- Set the master patch - master is responsible for updating all
        //  wall function patches
        virtual void setMaster();

        //- Create the averaging weights for cells which are bounded by
        //  multiple wall function faces
        virtual void createAveragingWeights();

        //- Helper function to return non-const access to an omega patch
        virtual WilcoxOmegaWallFunctionFvPatchScalarField& omegaPatch
        (
            const label patchI
        );

        //- Main driver to calculate the turbulence fields
        virtual void calculateTurbulenceFields
        (
            const turbulenceModel& turbulence,
            scalarField& G0,
            scalarField& omega0
        );

        //- Calculate the omega and G
        virtual void calculate
        (
            const turbulenceModel& turbulence,
            const List<scalar>& cornerWeights,
            const fvPatch& patch,
            scalarField& G,
            scalarField& omega
        );

        //- Return non-const access to the master patch ID
        virtual label& master()
        {
            return master_;
        }


public:

    //- Runtime type information
    TypeName("WilcoxOmegaWallFunction");


    // Constructors

        //- Construct from patch and internal field
        WilcoxOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        WilcoxOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // WilcoxOmegaWallFunctionFvPatchScalarField
        //  onto a new patch
        WilcoxOmegaWallFunctionFvPatchScalarField
        (
            const WilcoxOmegaWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        WilcoxOmegaWallFunctionFvPatchScalarField
        (
            const WilcoxOmegaWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new WilcoxOmegaWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        WilcoxOmegaWallFunctionFvPatchScalarField
        (
            const WilcoxOmegaWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new WilcoxOmegaWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return non-const access to the master's G field
            scalarField& G(bool init = false);

            //- Return non-const access to the master's omega field
            scalarField& omega(bool init = false);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs(const scalarField& weights);

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);

            //- Manipulate matrix with given weights
            virtual void manipulateMatrix
            (
                fvMatrix<scalar>& matrix,
                const scalarField& weights
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
