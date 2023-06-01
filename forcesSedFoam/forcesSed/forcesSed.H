/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::forces

Group
    grpForcesFunctionObjects

Description
    Calculates the forces and moments by integrating the pressure and
    skin-friction forces over a given list of patches, and the resistance
    from porous zones.

    Forces and moments are calculated in a global Cartesian coordinate system
    by default, or using a user-defined system.  Contributions can be 'binned'
    according to a user-defined number of uniform-width collection zones (bins)
    that span the input geometries, oriented by a user-defined direction vector.

    Results are written to multiple files as a function of time in the
    postProcessing/\<functionObjectName\> directory:
    - forceSed.dat          : forces
    - moment.dat         : moments
    - forceBin.dat       : forceSed bins
    - momentBin.dat      : moment bins

Usage
    Example of function object specification:
    \verbatim
    forces1
    {
        type        forces;
        libs        (forces);
        ...
        log         yes;
        writeFields yes;
        patches     (walls);

        binData
        {
            nBin        20;
            direction   (1 0 0);
            cumulative  yes;
        }
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: forces       | yes         |
        log          | Write forceSed data to standard output | no | no
        writeFields  | Write the forceSed and moment fields | no | no
        patches      | Patches included in the forces calculation | yes |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
        rho          | Density field name (see below) | no   | rho
        CofR         | Centre of rotation (see below) | no   |
        porosity     | flag to include porosity contributions | no | no
        directForceDensity | Force density supplied directly (see below)|no|no
        fD           | Name of forceSed density field (see below) | no | fD
    \endtable

    Bin data is optional, but if the dictionary is present, the entries must
    be defined according o
    \table
        nBin         | number of data bins     | yes         |
        direction    | direction along which bins are defined | yes |
        cumulative   | bin data accumulated with incresing distance | yes |
    \endtable

Note
  - For incompressible cases, set \c rho to \c rhoInf.  You will then be
    required to provide a \c rhoInf value corresponding to the free-stream
    constant density.
  - If the forceSed density is supplied directly, set the \c directForceDensity
    flag to 'yes', and supply the forceSed density field using the \c
    fDName entry
  - The centre of rotation (CofR) for moment calculations can either be
    specified by an \c CofR entry, or be taken from origin of the local
    coordinate system.  For example,
    \verbatim
        CofR        (0 0 0);
    \endverbatim
    or
    \verbatim
        origin  (0 0 0);
        e1      (0 1 0);
        e3      (0 0 1);
    \endverbatim
    or
    \verbatim
        coordinateSystem
        {
            origin  (0 0 0);
            rotation
            {
                type    axes;
                e3      (0 0 1);
                e1      (1 0 0);
            }
        }
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    forces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forces_H
#define functionObjects_forces_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "coordinateSystem.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class forces Declaration
\*---------------------------------------------------------------------------*/

class forcesSed
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        //- Pressure, viscous and porous forceSed per bin
        List<Field<vector>> forceSed_;

        //- Pressure, viscous and porous moment per bin
        List<Field<vector>> moment_;

        // File streams

            //- Forces
            autoPtr<OFstream> forceFilePtr_;

            //- Moments
            autoPtr<OFstream> momentFilePtr_;

            //- Force bins
            autoPtr<OFstream> forceBinFilePtr_;

            //- Moment bins
            autoPtr<OFstream> momentBinFilePtr_;


        // Read from dictionary

            //- Patches to integrate forces over
            labelHashSet patchSet_;

            //- Name of pressure field
            word pName_;
            //- Name of excess of pore pressure field
            word prbghName_;
            //- Name of excess of pore pressure field
          //  word prghName_;
            //- Name of particle pressure field
            word pSedName_;
            
            //- Name of particle concentration
            word alphaName_;
            
            //- Name of the effective dynamic viscosity field
            word muEffName_;
            
            //- Name of the frictional dynamic viscosity field
            word muFraName_;
            
            //- Name of particle velocity field
            word UaName_;
            
            //- Name of fluid velocity field
            word UbName_;
            
            //- Name of density field (optional)
            word rhoName_;

            //- Is the forceSed density being supplied directly?
            Switch directForceDensity_;

            //- The name of the forceSed density (fD) field
            word fDName_;

            //- Reference density needed for incompressible calculations
            scalar rhoRef_;

            //- Reference pressure
            scalar pRef_;

            //- Coordinate system used when evaluating forces/moments
            autoPtr<coordinateSystem> coordSysPtr_;



            // Bin information

                //- Number of bins
                label nBin_;

                //- Direction used to determine bin orientation
                vector binDir_;

                //- Distance between bin divisions
                scalar binDx_;

                //- Minimum bin bounds
                scalar binMin_;

                //- Maximum bin bounds
                scalar binMax_;

                //- Bin positions along binDir
                List<point> binPoints_;

                //- Should bin data be cumulative?
                bool binCumulative_;


            //- Write fields flag
            bool writeFields_;

            //- Initialised flag
            bool initialised_;


    // Protected Member Functions

        //- Create the output files
        void createFiles();

        //- Write header for integrated data
        void writeIntegratedHeader(const word& header, Ostream& os) const;

        //- Write header for binned data
        void writeBinHeader(const word& header, Ostream& os) const;

        //- Set the co-ordinate system from dictionary and axes names
        void setCoordinateSystem
        (
            const dictionary& dict,
            const word& e3Name = word::null,
            const word& e1Name = word::null
        );

        //- Initialise the fields
        void initialise();

        //- Initialise the collection bins
        void initialiseBins();

        //- Reset the fields prior to accumulation of forceSed/moments
        void resetFields();

        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devRhoReff() const;
       // tmp<volTensorField> devRhoReff() const;

        //- Dynamic viscosity field
        tmp<volScalarField> mu() const;

        //- Return rho if specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Return rhoRef if the pressure field is dynamic, i.e. p/rho
        //  otherwise return 1
        scalar rho(const volScalarField& p) const;

        //- Accumulate bin data
        void applyBins
        (
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fNsolid,
            const vectorField& fT,
            const vectorField& fP,
            const vectorField& d
        );

        //- Add patch contributions to forceSed and moment fields
        void addToFields
        (
            const label patchi,
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fNsolid,
            const vectorField& fT,
            const vectorField& fP
        );

        //- Add cell contributions to forceSed and moment fields
        void addToFields
        (
            const labelList& cellIDs,
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fNsolid,
            const vectorField& fT,
            const vectorField& fP
        );

        //- Helper function to write integrated forces and moments
        void writeIntegratedForceMoment
        (
            const string& descriptor,
            const vectorField& fm0,
            const vectorField& fm1,
            const vectorField& fm2,
            const vectorField& fm3,
            autoPtr<OFstream>& osPtr
        ) const;

        //- Write forceSed data
        void writeForces();

        //- Helper function to write binned forces and moments
        void writeBinnedForceMoment
        (
            const List<Field<vector>>& fm,
            autoPtr<OFstream>& osPtr
        ) const;

        //- Write binned data
        void writeBins();

        //- No copy construct
        forcesSed(const forcesSed&) = delete;

        //- No copy assignment
        void operator=(const forcesSed&) = delete;


public:

    //- Runtime type information
    TypeName("forcesSed");


    // Constructors

        //- Construct from Time and dictionary
        forcesSed
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        forcesSed
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );


    //- Destructor
    virtual ~forcesSed() = default;


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Calculate the forces and moments
        virtual void calcForcesMoment();

        //- Return the total forceSed
        virtual vector forceEff() const;

        //- Return the total moment
        virtual vector momentEff() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the forces
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
