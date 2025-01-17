/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "linearSpringDamper.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(linearSpringDamper, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        linearSpringDamper,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::linearSpringDamper::linearSpringDamper
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict)
{
    oldRestraintForce_ = Zero;
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::linearSpringDamper::~linearSpringDamper()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::linearSpringDamper::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    if (!anchor_)
    {
        anchor_.reset
        (
            Function1<vector>::New
            (
                "anchor",
                coeffDict(),
                &motion.time()
            )
        );
    }

    scalar t = motion.time().timeOutputValue();

    restraintPosition = motion.transform(refAttachmentPt_);

    // Current axis of the spring
    vector r = restraintPosition - anchor_->value(t);
    vector rDir = r/(mag(r) + VSMALL);

    vector v = motion.velocity(restraintPosition);

    scalar m = motion.mass();

    restraintMoment = Zero;

//     scalar dt =  motion.time().deltaTValue();
//
//     oldError_ = error_;
//     oldErrorIntegral_ = errorIntegral_;
//     error_ = (mag(r) - restLength_)/restLength_;
//     errorIntegral_ =
//         oldErrorIntegral_ + 0.5*(error_ + oldError_);
//
//     scalar errorDifferential = (error_ - oldError_)/dt;

    if (mag(r) > restLength_)
    {

//         factor_ =
//             P_*error_ + I_*errorIntegral_ + D_*errorDifferential;
//
//         factor_ = max(factor_, -1);

        scalar damping = psi_*2*m*wn_/numberOfChains_;
        scalar stiffness = sqr(wn_)*m/numberOfChains_;

        restraintForce =
            frelax_
            *(
                - damping*(rDir & v)*rDir
                - stiffness*(mag(r) - restLength_)*rDir
            )
            + (1-frelax_)*oldRestraintForce_;

        oldRestraintForce_ = restraintForce;
    }
    else
    {
        restraintForce = Zero;
        oldRestraintForce_ = Zero;
    }


    if (motion.report())
    {
        Info<< t << " " << restraintForce.x()   //2
            << " " << restraintForce.y()        //3
            << " " << restraintForce.z()        //4
            << " " << mag(r) - restLength_      //5
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::linearSpringDamper::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.readEntry("refAttachmentPt", refAttachmentPt_);
    psi_ = sDoFRBMRCoeffs_.getOrDefault<scalar>("psi", 1);
    sDoFRBMRCoeffs_.readEntry("wn", wn_);
    sDoFRBMRCoeffs_.readEntry("restLength", restLength_);
    sDoFRBMRCoeffs_.readEntry("numberOfChains", numberOfChains_);
    frelax_ = sDoFRBMRCoeffs_.getOrDefault<scalar>("frelax", 0.8);

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::linearSpringDamper::write
(
    Ostream& os
) const
{
    os.writeEntry("refAttachmentPt", refAttachmentPt_);
    os.writeEntry("psi", psi_);
    os.writeEntry("wn", wn_);
    os.writeEntry("restLength", restLength_);
    os.writeEntry("numberOfChains", numberOfChains_);
    os.writeEntryIfDifferent<scalar>("psi", 1, psi_);
    os.writeEntryIfDifferent<scalar>("frelax", 0.8, frelax_);
}


// ************************************************************************* //
