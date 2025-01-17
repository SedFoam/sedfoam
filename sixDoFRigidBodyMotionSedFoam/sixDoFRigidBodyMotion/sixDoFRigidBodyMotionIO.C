/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "sixDoFRigidBodyMotion.H"
#include "IOstreams.H"
#include "sixDoFSolver.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotion::read(const dictionary& dict)
{
    dict.readEntry("mass", mass_);
    dict.readEntry("momentOfInertia", momentOfInertia_);
    aRelax_ = dict.getOrDefault<scalar>("accelerationRelaxation", 1);
    aDamp_ = dict.getOrDefault<scalar>("accelerationDamping", 1);
    report_ = dict.getOrDefault<Switch>("report", false);

    restraints_.clear();
    addRestraints(dict);

    constraints_.clear();
    addConstraints(dict);

    return true;
}


void Foam::sixDoFRigidBodyMotion::write(Ostream& os) const
{
    motionState_.write(os);

    os.writeEntry("centreOfMass", initialCentreOfMass_);
    os.writeEntry("initialOrientation", initialQ_);
    os.writeEntry("mass", mass_);
    os.writeEntry("momentOfInertia", momentOfInertia_);
    os.writeEntry("accelerationRelaxation", aRelax_);
    os.writeEntry("accelerationDamping", aDamp_);
    os.writeEntry("report", report_);

    if (!restraints_.empty())
    {
        os.beginBlock("restraints");

        forAll(restraints_, rI)
        {
            const word& restraintType(restraints_[rI].type());

            os.beginBlock(restraints_[rI].name());

            os.writeEntry("sixDoFRigidBodyMotionRestraint", restraintType);

            restraints_[rI].write(os);

            os.endBlock();
        }

        os.endBlock();
    }

    if (!constraints_.empty())
    {
        os.beginBlock("constraints");

        forAll(constraints_, rI)
        {
            const word& constraintType(constraints_[rI].type());

            os.beginBlock(constraints_[rI].name());

            os.writeEntry("sixDoFRigidBodyMotionConstraint", constraintType);

            constraints_[rI].sixDoFRigidBodyMotionConstraint::write(os);

            constraints_[rI].write(os);

            os.endBlock();
        }

        os.endBlock();
    }

    if (solver_)
    {
        os  << indent << "solver";
        solver_->write(os);
    }
}


// ************************************************************************* //
