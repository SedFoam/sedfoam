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

\*---------------------------------------------------------------------------*/

#include "WilcoxOmegaWallFunctionFvPatchScalarField.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "nutkWallFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar WilcoxOmegaWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void WilcoxOmegaWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("WilcoxOmegaWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void WilcoxOmegaWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("kn") << kn_ << token::END_STATEMENT << nl; 
    //os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    //os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    //os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    //os.writeKeyword("beta1") << beta1_ << token::END_STATEMENT << nl;
}


void WilcoxOmegaWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& omega =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchI)
    {
        if (isA<WilcoxOmegaWallFunctionFvPatchScalarField>(bf[patchI]))
        {
            WilcoxOmegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchI);

            if (master == -1)
            {
                master = patchI;
            }

            opf.master() = master;
        }
    }
}


void WilcoxOmegaWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchI)
    {
        if (isA<WilcoxOmegaWallFunctionFvPatchScalarField>(bf[patchI]))
        {
            omegaPatches.append(patchI);

            const labelUList& faceCells = bf[patchI].patch().faceCells();
            forAll(faceCells, i)
            {
                label cellI = faceCells[i];
                weights[cellI]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(omegaPatches, i)
    {
        label patchI = omegaPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchI];
        cornerWeights_[patchI] = 1.0/wf.patchInternalField();
    }

    G_.setSize(dimensionedInternalField().size(), 0.0);
    omega_.setSize(dimensionedInternalField().size(), 0.0);

    initialised_ = true;
}


WilcoxOmegaWallFunctionFvPatchScalarField&
WilcoxOmegaWallFunctionFvPatchScalarField::omegaPatch(const label patchI)
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->dimensionedInternalField());

    const volScalarField::GeometricBoundaryField& bf = omega.boundaryField();

    const WilcoxOmegaWallFunctionFvPatchScalarField& opf =
        refCast<const WilcoxOmegaWallFunctionFvPatchScalarField>(bf[patchI]);

    return const_cast<WilcoxOmegaWallFunctionFvPatchScalarField&>(opf);
}


void WilcoxOmegaWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalarField& G0,
    scalarField& omega0
)
{
    // accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchI)
    {
        if (!cornerWeights_[patchI].empty())
        {
            WilcoxOmegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchI);

            const List<scalar>& w = cornerWeights_[patchI];
//            opf.calculate(turbulence, w, opf.patch(), omega0);
            opf.calculate(turbulence, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchI)
    {
        if (!cornerWeights_[patchI].empty())
        {
            WilcoxOmegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchI);

            opf == scalarField(omega0, opf.patch().faceCells());
        }
    }
}


void WilcoxOmegaWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbulence,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G,
    scalarField& omega
)
{
    const label patchI = patch.index();

    const scalarField& y = turbulence.y()[patchI];

 //   const scalar Cmu25 = pow025(Cmu_);

    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();

    const tmp<volScalarField> tnu = turbulence.nu();
    const scalarField& nuw = tnu().boundaryField()[patchI];

    const tmp<volScalarField> tnut = turbulence.nut();
    const volScalarField& nut = tnut();
    const scalarField& nutw = nut.boundaryField()[patchI];

    const fvPatchVectorField& Uw = turbulence.U().boundaryField()[patchI];

    const scalarField magGradUw(mag(Uw.snGrad()));

    // Set omega and G
    forAll(nutw, faceI)
    {
        label cellI = patch.faceCells()[faceI];

        scalar w = cornerWeights[faceI];

        scalar utau = max(1.0e-6,sqrt((nutw[faceI] + nuw[faceI])*magGradUw[faceI]));
        scalar knplus = kn_*utau/nuw[faceI];
        scalar SR = 0e0;
        //Info<<utau<<endl;
        //Info<<knplus<<endl;
        if (knplus<=25e0)
          { 
            SR = pow(50e0/knplus,2);
          }
        else
          { 
            SR = 100e0/knplus;
          }

        omega[cellI] +=  w*pow(utau,2)*SR/nuw[faceI];

    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WilcoxOmegaWallFunctionFvPatchScalarField::WilcoxOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    kn_(1e-6),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


WilcoxOmegaWallFunctionFvPatchScalarField::WilcoxOmegaWallFunctionFvPatchScalarField
(
    const WilcoxOmegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    kn_(ptf.kn_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


WilcoxOmegaWallFunctionFvPatchScalarField::WilcoxOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    kn_(dict.lookupOrDefault<scalar>("kn", 1e-6)),    
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


WilcoxOmegaWallFunctionFvPatchScalarField::WilcoxOmegaWallFunctionFvPatchScalarField
(
    const WilcoxOmegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    kn_(owfpsf.kn_),   
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


WilcoxOmegaWallFunctionFvPatchScalarField::WilcoxOmegaWallFunctionFvPatchScalarField
(
    const WilcoxOmegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    kn_(owfpsf.kn_),   
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
///*
scalarField& WilcoxOmegaWallFunctionFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }
    return omegaPatch(master_).G();
}
//*/

scalarField& WilcoxOmegaWallFunctionFvPatchScalarField::omega(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void WilcoxOmegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>(turbulenceModel::typeName);

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbulence, G(true), omega(true));
 //      calculateTurbulenceFields(turbulence,omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbulence.GName())
        );

    FieldType& omega = const_cast<FieldType&>(dimensionedInternalField());

    forAll(*this, faceI)
    {
        label cellI = patch().faceCells()[faceI];

        G[cellI] = G0[cellI];
        omega[cellI] = omega0[cellI];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void WilcoxOmegaWallFunctionFvPatchScalarField::updateCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>(turbulenceModel::typeName);

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbulence, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbulence.GName())
        );

    FieldType& omega = const_cast<FieldType&>(dimensionedInternalField());

    scalarField& omegaf = *this;

    // only set the values if the weights are > tolerance
    forAll(weights, faceI)
    {
        scalar w = weights[faceI];

        if (w > tolerance_)
        {
            label cellI = patch().faceCells()[faceI];

            G[cellI] = (1.0 - w)*G[cellI] + w*G0[cellI];
            omega[cellI] = (1.0 - w)*omega[cellI] + w*omega0[cellI];
            omegaf[faceI] = omega[cellI];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void WilcoxOmegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void WilcoxOmegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    // filter weights so that we only apply the constraint where the
    // weight > SMALL
    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintomega(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& omega =
        dimensionedInternalField();

    label nConstrainedCells = 0;


    forAll(weights, faceI)
    {
        // only set the values if the weights are > tolerance
        if (weights[faceI] > tolerance_)
        {
            nConstrainedCells++;

            label cellI = faceCells[faceI];

            constraintCells.append(cellI);
            constraintomega.append(omega[cellI]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << nConstrainedCells
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues
    (
        constraintCells,
        scalarField(constraintomega.xfer())
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void WilcoxOmegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
    writeLocalEntries(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    WilcoxOmegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
