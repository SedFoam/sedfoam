/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::functionObjects::forceCoeffsSed

Group
    grpForcesFunctionObjects

Description
    Extends the \c forcesSed functionObject by providing coefficients for:
    - drag, side and lift forcesSed (Cd, Cs, and Cl)
    - roll, pitch and yaw moments (CmRoll, CmPitch, and CmYaw)
    - front and rear axle forceSed contributions (C(f) and C(r)) wherein

    \verbatim
        Cd(f/r) = 0.5*Cd \pm CmRoll
        Cs(f/r) = 0.5*Cs \pm CmYaw
        Cl(f/r) = 0.5*Cl \pm CmPitch
    \endverbatim

    The data can optionally be output into bins, defined in a given direction.

    The binned data provides the total and consitituent components per bin:
    - total coefficient
    - pressure coefficient contribution
    - viscous coefficient contribution
    - porous coefficient contribution

    Data is written into multiple files in the
    postProcessing/\<functionObjectName\> directory:
    - coefficient.dat   : integrated coefficients over all geometries
    - CdBin.dat         : drag coefficient bins
    - CsBin.dat         : side coefficient bins
    - ClBin.dat         : lift coefficient bins
    - CmRollBin.dat     : roll moment coefficient bins
    - CmPitchBin.dat    : pitch moment coefficient bins
    - CmYawBin.dat      : yaw moment coefficient bins

Usage
    Example of function object specification:
    \verbatim
    forceCoeffs1
    {
        type        forceCoeffsSed;
        libs        (forcesSed);
        ...
        log         yes;
        writeFields yes;
        patches     (walls);

        // input keywords for directions of forceSed/moment coefficients
        // refer below for options, and relations

        magUInf     100;
        lRef        3.5;
        Aref        2.2;
        porosity    no;

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
        Property     | Description                          | Required | Default
        type         | Type name: forceCoeffsSed               | yes |
        log          | Write forceSed data to standard output  | no  | no
        writeFields  | Write forceSed,moment coefficient fields | no | no
        patches      | Patches included in the forcesSed calculation | yes |
        magUInf      | Free stream velocity magnitude       | yes |
        rhoInf       | Free stream density | for compressible cases |
        lRef         | Reference length scale for moment calculations | yes |
        Aref         | Reference area                       | yes |
        porosity     | Include porosity contributions       | no  | false
    \endtable

    Bin data is optional, but if the dictionary is present, the entries must
    be defined according to following:
    \table
        nBin         | Number of data bins                    | yes |
        direction    | Direction along which bins are defined | yes |
        cumulative   | Bin data accumulated with incresing distance | yes |
    \endtable

    Input of forceSed/moment coefficient directions:
    - require an origin, and two orthogonal directions; the remaining orthogonal
    direction is determined accordingly.
    - can be added by the three options below.

    \verbatim
        CofR        (0 0 0); // Centre of rotation
        dragDir     (1 0 0);
        liftDir     (0 0 1);
    \endverbatim

    \verbatim
        origin (0 0 0);
        e1     (1 0 0);
        e3     (0 0 1); // combinations: (e1, e2) or (e2, e3) or (e3, e1)
    \endverbatim

    \verbatim
        coordinateSystem
        {
            origin  (0 0 0);
            rotation
            {
                type axes;
                e1 (1 0 0);
                e3 (0 0 1); // combinations: (e1, e2) or (e2, e3) or (e3, e1)
            }
        }
    \endverbatim

    The default direction relations are shown below:

    \table
        Property     | Description           | Alias | Direction
        dragDir      | Drag direction        | e1    | (1 0 0)
        sideDir      | Side forceSed direction  | e2    | (0 1 0)
        liftDir      | Lift direction        | e3    | (0 0 1)
        rollAxis     | Roll axis             | e1    | (1 0 0)
        pitchAxis    | Pitch axis            | e2    | (0 1 0)
        yawAxis      | Yaw axis              | e3    | (0 0 1)
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::timeControl
    Foam::functionObjects::forcesSed

SourceFiles
    forceCoeffsSed.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forceCoeffs_H
#define functionObjects_forceCoeffs_H

#include "forcesSed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class forceCoeffsSed Declaration
\*---------------------------------------------------------------------------*/

class forceCoeffsSed
:
    public forcesSed
{
    // Private data

        // Free-stream conditions

            //- Free-stream velocity magnitude
            scalar magUInf_;


        // Reference scales

            //- Reference length [m]
            scalar lRef_;

            //- Reference area [m^2]
            scalar Aref_;


        // File streams

            //- Integrated coefficients
            autoPtr<OFstream> coeffFilePtr_;

            //- Drag coefficient
            autoPtr<OFstream> CdBinFilePtr_;

            //- Side coefficient
            autoPtr<OFstream> CsBinFilePtr_;

            //- Lift coefficient
            autoPtr<OFstream> ClBinFilePtr_;

            //- Roll moment coefficient
            autoPtr<OFstream> CmRollBinFilePtr_;

            //- Pitch moment coefficient
            autoPtr<OFstream> CmPitchBinFilePtr_;

            //- Yaw moment coefficient
            autoPtr<OFstream> CmYawBinFilePtr_;


    // Private Member Functions

        //- No copy construct
        forceCoeffsSed(const forceCoeffsSed&) = delete;

        //- No copy assignment
        void operator=(const forceCoeffsSed&) = delete;


protected:

    // Protected Member Functions

        //- Create the output files
        void createFiles();

        //- Write header for integrated data
        void writeIntegratedHeader(const word& header, Ostream& os) const;

        //- Write header for binned data
        void writeBinHeader(const word& header, Ostream& os) const;

        //- Write integrated data
        void writeIntegratedData
        (
            const word& title,
            const List<Field<scalar>>& coeff
        ) const;

        //- Write binned data
        void writeBinData(const List<Field<scalar>> coeffs, Ostream& os) const;


public:

    //- Runtime type information
    TypeName("forceCoeffsSed");


    // Constructors

        //- Construct from Time and dictionary
        forceCoeffsSed
        (
            const word& name,
            const Time& runTime,
            const dictionary&,
            const bool readFields = true
        );


    //- Destructor
    virtual ~forceCoeffsSed() = default;


    // Member Functions

        //- Read the forcesSed data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write the forcesSed
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
