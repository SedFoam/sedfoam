/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::sixDoFRigidBodyMotionRestraints::linearSpringDamper

Description
    sixDoFRigidBodyMotionRestraints model.  Linear Spring-Damper system.

    Spring-damper system for restraint. Acts as a "soft" rope when the
    distance between anchor and refAttachmentPt exceed the restLength,
    a system of spring-damper is activated to simulate a soft rope
    de-acceleration

SourceFiles
    linearSpringDamper.C

\*---------------------------------------------------------------------------*/

#ifndef linearSpringDamper_H
#define linearSpringDamper_H

#include "sixDoFRigidBodyMotionRestraint.H"
#include "point.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class linearSpringDamper Declaration
\*---------------------------------------------------------------------------*/

class linearSpringDamper
:
    public sixDoFRigidBodyMotionRestraint
{
    // Private data

        //- Anchor point, where the spring is attached to an immovable
        //  object
        mutable autoPtr<Function1<vector>> anchor_;

        //- Reference point of attachment to the solid body
        point refAttachmentPt_;

        //- Damping factor
        scalar psi_;

        //- System frequency  [1/s]
        scalar wn_;

        //- Number of chains in the same system
        scalar numberOfChains_;

        //- Rest length - length of spring when no forces are applied to it
        scalar restLength_;

        //- Cache old restrain force
        mutable vector oldRestraintForce_;

        //- Relaxation factor
        scalar frelax_;


        /*
        //- Error between ideal length and actual length
        mutable scalar error_;

         //- Old error between ideal length and actual length
        mutable scalar oldError_;

        //- Integral length error in time
        mutable scalar errorIntegral_;

        //- Old integral length error in time
        mutable scalar oldErrorIntegral_;

        //- Proportional constant
        scalar P_;

        //- Integral constant
        scalar I_;

        //- Differential constant
        scalar D_;

        //- Initial force
        vector f0_;

        mutable scalar wn_;
        */


public:

    //- Runtime type information
    TypeName("linearSpringDamper");


    // Constructors

        //- Construct from components
        linearSpringDamper
        (
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionRestraint>
            (
                new linearSpringDamper(*this)
            );
        }


    //- Destructor
    virtual ~linearSpringDamper();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const sixDoFRigidBodyMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
