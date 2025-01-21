/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "sixDoFRigidBodyMotionPlaneConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionConstraints
{
    defineTypeNameAndDebug(plane, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionConstraint,
        plane,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::plane::plane
(
    const word& name,
    const dictionary& sDoFRBMCDict,
    const sixDoFRigidBodyMotion& motion
)
:
    sixDoFRigidBodyMotionConstraint(name, sDoFRBMCDict, motion)
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::plane::~plane()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionConstraints::plane::setCentreOfRotation
(
    point& CofR
) const
{
    CofR = centreOfRotation_;
}


void Foam::sixDoFRigidBodyMotionConstraints::plane::constrainTranslation
(
    pointConstraint& pc
) const
{
    pc.applyConstraint(normal_);
}


void Foam::sixDoFRigidBodyMotionConstraints::plane::constrainRotation
(
    pointConstraint& pc
) const
{}


bool Foam::sixDoFRigidBodyMotionConstraints::plane::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyMotionConstraint::read(sDoFRBMCDict);

    centreOfRotation_ = sDoFRBMCCoeffs_.getOrDefault
    (
        "centreOfRotation",
        motion_.initialCentreOfMass()
    );

    sDoFRBMCCoeffs_.readEntry("normal", normal_);

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraints::plane::write
(
    Ostream& os
) const
{
    os.writeEntry("centreOfRotation", centreOfRotation_);
    os.writeEntry("normal", normal_);
}

// ************************************************************************* //
