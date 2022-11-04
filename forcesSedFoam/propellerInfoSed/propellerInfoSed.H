/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::functionObjects::propellerInfoSed

Group
    grpForcesFunctionObjects

Description
    Calculates propeller performance and wake field properties.

    Controlled by executeControl:
    - Propeller performance
      - Thrust coefficient, Kt
      - Torque coefficient, 10*Kq
      - Advance coefficient, J
      - Open water efficiency, etaO
      - Written to postProcessing/<name>/<time>/propellerPerformance.dat

    Controlled by writeControl:
    - Wake field text file
      - Wake: 1 - UzMean/URef
      - Velocity in cylindrical coordinates at xyz locations
      - Written to postProcessing/<name>/<time>/wake.dat
    - Axial wake field text file
      - 1 - Uz/URef at r/R and angle
      - Written to postProcessing/<name>/<time>/axialWake.dat
    - Velocity surface
      - Written to postProcessing/<name>/surfaces/time>/disk.<fileType>

Usage
    Example of function object specification:
    \verbatim
    propellerInfo1
    {
        type            propellerInfoSed;
        libs            (forces);
        writeControl    writeTime;

        patches         ("propeller.*");

        URef            5; // Function1 type; 'constant' form shown here
        rho             rhoInf; // incompressible
        rhoInf          1.2;

        // Optionally write propeller performance data
        writePropellerPerformance yes;


        // Propeller data:

        // Radius
        radius          0.1;

        rotationMode    specified; // specified | MRF

        // rotationMode = specified:
        origin          (0 -0.1 0);
        n               25.15;
        axis            (0 1 0);

        // Optional reference direction for angle (alpha) = 0
        alphaAxis       (1 0 0);


        //// rotationMode = mrf
        //// MRF             MRFZoneName;
        //// (origin, n and axis retrieved from MRF model)

        // Optionally write wake text files
        // Note: controlled by writeControl
        writeWakeFields yes;

        // Sample plane (disk) properties
        // Note: controlled by writeControl
        sampleDisk
        {
            surfaceWriter   vtk;
            r1              0.05;
            r2              0.2;
            nTheta          36;
            nRadial         10;
            interpolationScheme cellPoint;
            errorOnPointNotFound false;
        }
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property        | Description                   | Required | Deflt value
        type            | Type name: propellerInfoSed      | yes      |
        log             | Write to standard output      | no       | no
        patches         | Patches included in the forces calculation | yes |
        p               | Pressure field name           | no       | p
        U               | Velocity field name           | no       | U
        rho             | Density field name            | no       | rho
        URef            | Reference velocity            | yes      |
        rotationMode    | Rotation mode (see below)     | yes      |
        origin          | Sample disk centre            | no*      |
        n               | Revolutions per second        | no*      |
        axis            | Propeller axis                | no*      |
        alphaAxis       | Axis that defines alpha=0 dir | no       |
        MRF             | Name of MRF zone              | no*      |
        originOffset    | Origin offset for MRF mode    | no       | (0 0 0)
        writePropellerPerformance| Write propeller performance text file | yes |
        writeWakeFields | Write wake field text files   | yes      |
        surfaceWriter   | Sample disk surface writer    | no*      |
        r1              | Sample disk inner radius      | no       | 0
        r2              | Sample disk outer radius      | no*      |
        nTheta          | Divisions in theta direction  | no*      |
        nRadial         | Divisions in radial direction | no*      |
        interpolationScheme | Sampling interpolation scheme | no* | cell
    \endtable


Note
- URef is a scalar Function1 type, i.e. supports constant, table, lookup values
- rotationMode is used to set the origin, axis and revolutions per second
  - if set to 'specified' all 3 entries are required
    - note: origin is the sample disk origin
  - if set to 'MRF' only the MRF entry is required
    - to move the sample disk away from the MRF origin, use the originOffset
- if writePropellerPerformance is set to on|true:
  - propellerPerformance text file will be written
- if writeWakeFields is set to on|true:
  - wake and axialWake text files will be written
- if the surfaceWriter entry is set, the sample disk surface will be written
  - extents set according to the r1 and r2 entries
  - discretised according to the nTheta and nRadial entries

See also
    Foam::functionObject::forces
    Foam::Function1

SourceFiles
    propellerInfoSed.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_propellerInfo_H
#define functionObjects_propellerInfo_H

#include "forcesSed.H"
#include "Enum.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class Function1;

class surfaceWriter;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class propellerInfoSed Declaration
\*---------------------------------------------------------------------------*/

class propellerInfoSed
:
    public forcesSed
{

public:

    enum class rotationMode
    {
        SPECIFIED,
        MRF
    };

    static const Enum<rotationMode> rotationModeNames_;


protected:

    // Protected data

        //- Propeller radius
        scalar radius_;

        //- Reference velocity
        autoPtr<Function1<scalar>> URefPtr_;

        //- Rotation mode
        rotationMode rotationMode_;

        //- Name of MRF zone (if applicable)
        word MRFName_;

        //- Propeller speed (revolutions per second)
        scalar n_;

        //- Flag to write performance data
        bool writePropellerPerformance_;

        //- Propeller performance file
        autoPtr<OFstream> propellerPerformanceFilePtr_;

        //- Flag to write wake fields
        bool writeWakeFields_;


        // Wake field output

            //- Surface writer
            autoPtr<surfaceWriter> surfaceWriterPtr_;

            //- Number of surface divisions in theta direction
            label nTheta_;

            //- Number of surface divisions in radial direction
            label nRadial_;

            //- Surface points
            pointField points_;

            //- Flag to raise an error if the sample point is not found in the
            //- mesh. Default = false to enable. e.g. reduced geometry/symmetric
            //- cases
            bool errorOnPointNotFound_;

            //- Surface faces
            faceList faces_;

            //- Surface point sample cell IDs
            labelList cellIds_;

            //- List of participating points (parallel reduced)
            boolList pointMask_;

            //- Interpolation scheme
            word interpolationScheme_;

            //- Wake field file
            autoPtr<OFstream> wakeFilePtr_;

            //- Axial wake field file
            autoPtr<OFstream> axialWakeFilePtr_;

            //- Default value when a sample point is not found; default =
            //- scalar::min
            scalar nanValue_;


    // Protected Member Functions

        //- Set the coordinate system
        void setCoordinateSystem(const dictionary& dict);

        //- Set the rotational speed
        void setRotationalSpeed();

        //- Create output files
        void createFiles();

        //- Return the velocity field
        const volVectorField& U() const;


        // Propeller performance text file

            //- Write the wake fields
            void writePropellerPerformance();


        // Wake text files

            //- Write the wake text file
            void writeWake(const vectorField& U, const scalar URef);

            //- Write the axial wake text file
            void writeAxialWake(const vectorField& U, const scalar URef);

            //- Write the wake fields
            void writeWakeFields(const scalar URef);


        // Sample surface functions

            //- Set the faces and points for the sample surface
            void setSampleDiskGeometry
            (
                const coordinateSystem& coordSys,
                const scalar r1,
                const scalar r2,
                const scalar nTheta,
                const label nRadius,
                faceList& faces,
                pointField& points
            ) const;

            //- Set the sample surface based on dictionary settings
            void setSampleDiskSurface(const dictionary& dict);

            //- Set the sample cells corresponding to the sample points
            void updateSampleDiskCells();

            //- Return the area average of a field
            scalar meanSampleDiskField(const scalarField& field) const;

            //- Write the sample surface
            void writeSampleDiskSurface
            (
                const vectorField& U,
                const vectorField& Ur,
                const scalar URef
            );

            //- Interpolate from the mesh onto the sample surface
            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>& psi,
                const Type& defaultValue
            ) const;

        //- No copy construct
        propellerInfoSed(const propellerInfoSed&) = delete;

        //- No copy assignment
        void operator=(const propellerInfoSed&) = delete;


public:

    //- Runtime type information
    TypeName("propellerInfoSed");


    // Constructors

        //- Construct from Time and dictionary
        propellerInfoSed
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        propellerInfoSed
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );


    //- Destructor
    virtual ~propellerInfoSed() = default;


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the forces
        virtual bool write();

        void UpdateMesh(const mapPolyMesh& mpm);

        void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
