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

\*---------------------------------------------------------------------------*/

#include "propellerInfoSed.H"
#include "cylindricalCS.H"
#include "fvMesh.H"
#include "IOMRFZoneList.H"
#include "mathematicalConstants.H"
#include "interpolation.H"
#include "Function1.H"
#include "surfaceWriter.H"
#include "treeDataCell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(propellerInfoSed, 0);
    addToRunTimeSelectionTable(functionObject, propellerInfoSed, dictionary);
}
}

const Foam::Enum<Foam::functionObjects::propellerInfoSed::rotationMode>
Foam::functionObjects::propellerInfoSed::rotationModeNames_
({
    { rotationMode::SPECIFIED, "specified" },
    { rotationMode::MRF, "MRF" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::propellerInfoSed::setCoordinateSystem
(
    const dictionary& dict
)
{
    vector origin(Zero);
    vector axis(Zero);

    switch (rotationMode_)
    {
        case rotationMode::SPECIFIED:
        {
            origin = dict.get<vector>("origin");
            axis = dict.get<vector>("axis");
            axis.normalise();

            n_ = dict.get<scalar>("n");
            break;
        }
        case rotationMode::MRF:
        {
            MRFName_ = dict.get<word>("MRF");
            const auto* MRFZones =
                mesh_.cfindObject<IOMRFZoneList>("MRFProperties");

            if (!MRFZones)
            {
                FatalIOErrorInFunction(dict)
                    << "Unable to find MRFProperties in the database. "
                    << "Is this an MRF case?"
                    << exit(FatalIOError);
            }
            const auto& mrf = MRFZones->MRFZoneList::getFromName(MRFName_);
            vector offset = dict.getOrDefault("originOffset", vector::zero);
            origin = offset + mrf.origin();
            axis = mrf.axis();

            // Convert rad/s to revolutions per second
            n_ = (mrf.Omega() & axis)/constant::mathematical::twoPi;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << rotationModeNames_[rotationMode_]
                << abort(FatalError);
        }
    }

    vector alphaAxis;
    if (!dict.readIfPresent("alphaAxis", alphaAxis))
    {
        // Value has not been set - find vector orthogonal to axis

        vector cand(Zero);
        scalar minDot = GREAT;
        for (direction d = 0; d < 3; ++d)
        {
            vector test(Zero);
            test[d] = 1;
            scalar dotp = mag(test & axis);
            if (dotp < minDot)
            {
                minDot = dotp;
                cand = test;
            }
        }

        alphaAxis = axis ^ cand;
    }

    alphaAxis.normalise();

    coordSysPtr_.reset(new coordSystem::cylindrical(origin, axis, alphaAxis));
}


void Foam::functionObjects::propellerInfoSed::setRotationalSpeed()
{
    switch (rotationMode_)
    {
        case rotationMode::SPECIFIED:
        {
            // Set on dictionary re-read
            break;
        }
        case rotationMode::MRF:
        {
            const auto* MRFZones =
                mesh_.cfindObject<IOMRFZoneList>("MRFProperties");

            if (!MRFZones)
            {
                FatalErrorInFunction
                    << "Unable to find MRFProperties in the database. "
                    << "Is this an MRF case?"
                    << exit(FatalError);
            }

            const auto& mrf = MRFZones->MRFZoneList::getFromName(MRFName_);

            // Convert rad/s to revolutions per second
            n_ = (mrf.Omega() & mrf.axis())/constant::mathematical::twoPi;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << rotationModeNames_[rotationMode_]
                << abort(FatalError);
        }
    }
}


void Foam::functionObjects::propellerInfoSed::createFiles()
{
    if (!writeToFile())
    {
        return;
    }

    if (writePropellerPerformance_ && !propellerPerformanceFilePtr_)
    {
        propellerPerformanceFilePtr_ = createFile("propellerPerformance");
        auto& os = propellerPerformanceFilePtr_();

        writeHeader(os, "Propeller performance");
        writeHeaderValue(os, "CofR", coordSysPtr_->origin());
        writeHeaderValue(os, "Radius", radius_);
        writeHeaderValue(os, "Axis", coordSysPtr_->e3());

        writeHeader(os, "");

        writeCommented(os, "Time");
        writeTabbed(os, "n");
        writeTabbed(os, "URef");
        writeTabbed(os, "J");
        writeTabbed(os, "KT");
        writeTabbed(os, "10*KQ");
        writeTabbed(os, "eta0");
        os  << nl;
    }

    if (writeWakeFields_)
    {
        if (!wakeFilePtr_) wakeFilePtr_ = createFile("wake");
        if (!axialWakeFilePtr_) axialWakeFilePtr_ = createFile("axialWake");
    }
}


const Foam::volVectorField& Foam::functionObjects::propellerInfoSed::U() const
{
    const auto* UPtr = mesh_.cfindObject<volVectorField>(UName_);

    if (!UPtr)
    {
        FatalErrorInFunction
            << "Unable to find velocity field " << UName_
            << " . Available vector fields are: "
            << mesh_.names<volVectorField>()
            << exit(FatalError);
    }

    return *UPtr;
}


void Foam::functionObjects::propellerInfoSed::setSampleDiskGeometry
(
    const coordinateSystem& coordSys,
    const scalar r1,
    const scalar r2,
    const scalar nTheta,
    const label nRadius,
    faceList& faces,
    pointField& points
) const
{
    label nPoint = nRadius*nTheta;
    if (r1 < SMALL)
    {
        nPoint += 1; // 1 for origin
    }
    else
    {
        nPoint += nTheta;
    }
    const label nFace = nRadius*nTheta;

    points.setSize(nPoint);
    faces.setSize(nFace);

    const point& origin = coordSys.origin();
    const scalar zCoord = 0;
    label pointi = 0;
    for (int radiusi = 0; radiusi <= nRadius; ++radiusi)
    {
        if (r1 < SMALL && radiusi == 0)
        {
            points[pointi++] = origin;
        }
        else
        {
            const scalar r = r1 + radiusi*(r2 - r1)/nRadius;

            for (label i = 0; i < nTheta; ++i)
            {
                point p
                (
                    r,
                    (i/scalar(nTheta))*constant::mathematical::twoPi,
                    zCoord
                );

                points[pointi++] = coordSys.globalPosition(p);
            }
        }
    }


    const List<label> ptIDs(identity(nTheta));

    // Faces
    label facei = 0;
    label pointOffset0 = 0;
    label radiusOffset = 0;
    DynamicList<label> facePts(4);
    for (int radiusi = 0; radiusi < nRadius; ++radiusi)
    {
        if (r1 < SMALL && radiusi == 0)
        {
            radiusOffset = 1;
            pointOffset0 = 1;

            // Adding faces as a fan about the origin
            for (label thetai = 1; thetai <= nTheta; ++thetai)
            {
                facePts.clear();

                // Append triangle about origin
                facePts.append(0);
                facePts.append(thetai);
                facePts.append(1 + ptIDs.fcIndex(thetai - 1));

                faces[facei++] = face(facePts);
            }
        }
        else
        {
            for (label thetai = 0; thetai < nTheta; ++thetai)
            {
                facePts.clear();

                label offset = pointOffset0 + (radiusi-radiusOffset)*nTheta;

                // Inner
                facePts.append(offset + ptIDs.fcIndex(thetai - 1));
                facePts.append(offset + ptIDs.fcIndex(thetai));

                // Outer
                facePts.append(offset + nTheta + ptIDs.fcIndex(thetai));
                facePts.append(offset + nTheta + ptIDs.fcIndex(thetai - 1));

                faces[facei++] = face(facePts);
            }
        }
    }
}


void Foam::functionObjects::propellerInfoSed::setSampleDiskSurface
(
    const dictionary& dict
)
{
    const dictionary& sampleDiskDict(dict.subDict("sampleDisk"));

    const scalar r1 = sampleDiskDict.getScalar("r1");
    const scalar r2 = sampleDiskDict.getScalar("r2");

    nTheta_ = sampleDiskDict.getLabel("nTheta");
    nRadial_ = sampleDiskDict.getLabel("nRadial");

    setSampleDiskGeometry
    (
        coordSysPtr_(),
        r1,
        r2,
        nTheta_,
        nRadial_,
        faces_,
        points_
    );

    // Surface writer
    word surfWriterType;
    if (sampleDiskDict.readIfPresent("surfaceWriter", surfWriterType))
    {
        const auto writeOptions = sampleDiskDict.subOrEmptyDict("writeOptions");

        surfaceWriterPtr_ = surfaceWriter::New
        (
            surfWriterType,
            writeOptions.subOrEmptyDict(surfWriterType)
        );

        // Use outputDir/TIME/surface-name
        surfaceWriterPtr_->useTimeDir(true);
    }

    errorOnPointNotFound_ =
        sampleDiskDict.getOrDefault("errorOnPointNotFound", false);

    updateSampleDiskCells();
}


void Foam::functionObjects::propellerInfoSed::updateSampleDiskCells()
{
    if (!writeWakeFields_)
    {
        return;
    }

    treeBoundBox bb(points_);
    bb.inflate(0.05);
    DynamicList<label> treeCellIDs(10*points_.size());

    const auto& meshCells = mesh_.cells();
    const auto& meshFaces = mesh_.faces();
    const auto& meshPoints = mesh_.points();

    forAll(meshCells, celli)
    {
        bool found = false;

        for (const label facei : meshCells[celli])
        {
            for (const label fpi : meshFaces[facei])
            {
                if (bb.contains(meshPoints[fpi]))
                {
                    found = true;
                    break;
                }
            }

            if (found)
            {
                treeCellIDs.append(celli);
                break;
            }
        }
    }

    indexedOctree<treeDataCell> tree
    (
        treeDataCell(true, mesh_, std::move(treeCellIDs), polyMesh::CELL_TETS),
        bb,
        10,
        100,
        10
    );

    cellIds_.setSize(points_.size(), -1);
    pointMask_.setSize(points_.size(), false);

    // Kick the tet base points calculation to avoid parallel comms later
    (void)mesh_.tetBasePtIs();

    const auto& cellLabels = tree.shapes().cellLabels();

    forAll(points_, pointi)
    {
        const vector& pos = points_[pointi];

//        label meshCelli = mesh_.findCell(pos);
        label treeCelli = tree.findInside(pos);

        label proci = treeCelli >= 0 ? Pstream::myProcNo() : -1;

        reduce(proci, maxOp<label>());

        pointMask_[pointi] = treeCelli != -1;

        if (proci >= 0)
        {
            cellIds_[pointi] =
                proci == Pstream::myProcNo()
              ? cellLabels[treeCelli]
              : -1;
        }
        else
        {
            if (errorOnPointNotFound_)
            {
                FatalErrorInFunction
                    << "Position " << pos << " not found in mesh"
                    << abort(FatalError);
            }
            else
            {
                DebugInfo
                    << "Position " << pos << " not found in mesh"
                    << endl;
            }
        }
    }

    Pstream::listCombineGather(pointMask_, orEqOp<bool>());
    Pstream::listCombineScatter(pointMask_);
}


Foam::scalar Foam::functionObjects::propellerInfoSed::meanSampleDiskField
(
    const scalarField& field
) const
{
    if (field.size() != points_.size())
    {
        FatalErrorInFunction
            << "Inconsistent field sizes: input:" << field.size()
            << " points:" << points_.size()
            << abort(FatalError);
    }

    scalar sumArea = 0;
    scalar sumFieldArea = 0;
    forAll(faces_, facei)
    {
        const face& f = faces_[facei];

        bool valid = true;
        scalar faceValue = 0;
        for (const label pti : f)
        {
            // Exclude contributions where sample cell for point was not found
            if (!pointMask_[pti])
            {
                valid = false;
                break;
            }
            faceValue += field[pti];
        }

        if (valid)
        {
            scalar area = f.mag(points_);
            sumArea += area;
            sumFieldArea += faceValue/f.size()*area;
        }
    }

    return sumFieldArea/(sumArea + ROOTVSMALL);
}


void Foam::functionObjects::propellerInfoSed::writeSampleDiskSurface
(
    const vectorField& U,
    const vectorField& Ur,
    const scalar URef
)
{
    // Write surface
    if (!surfaceWriterPtr_)
    {
        return;
    }


    // Time-aware, with time spliced into the output path
    surfaceWriterPtr_->isPointData(true);
    surfaceWriterPtr_->beginTime(time_);
    surfaceWriterPtr_->open
    (
        points_,
        faces_,
        (baseFileDir() / name() / "surfaces" / "disk"),
        false  // serial - already merged
    );
    surfaceWriterPtr_->nFields(4); // Legacy VTK
    surfaceWriterPtr_->write("U", U);
    surfaceWriterPtr_->write("Ur", Ur);
    surfaceWriterPtr_->write("UNorm", U/URef);
    surfaceWriterPtr_->write("UrNorm", Ur/URef);
    surfaceWriterPtr_->endTime();
    surfaceWriterPtr_->clear();
}


void Foam::functionObjects::propellerInfoSed::writePropellerPerformance()
{
    if (!writePropellerPerformance_)
    {
        return;
    }

    // Update n_
    setRotationalSpeed();

    const vector sumForce(sum(forceSed_[0]) + sum(forceSed_[1]) + sum(forceSed_[2]));
    const vector sumMoment(sum(moment_[0]) + sum(moment_[1]) + sum(moment_[2]));

    const scalar diameter = 2*radius_;
    const scalar URef = URefPtr_->value(time_.timeOutputValue());
    const scalar j = -URef/mag(n_ + ROOTVSMALL)/diameter;
    const scalar denom = rhoRef_*sqr(n_)*pow4(diameter);
    const scalar kt = (sumForce & coordSysPtr_->e3())/denom;
    const scalar kq =
        -sign(n_)*(sumMoment & coordSysPtr_->e3())/(denom*diameter);
    const scalar etaO = kt*j/(kq*constant::mathematical::twoPi + ROOTVSMALL);

    if (writeToFile())
    {
        auto& os = propellerPerformanceFilePtr_();

        writeCurrentTime(os);
        os  << tab << n_
            << tab << URef
            << tab << j
            << tab << kt
            << tab << 10*kq
            << tab << etaO
            << nl;

        os.flush();
    }

    Log << type() << " " << name() <<  " output:" << nl
        << "    Revolutions per second, n : " << n_ << nl
        << "    Reference velocity, URef  : " << URef << nl
        << "    Advance coefficient, J    : " << j << nl
        << "    Thrust coefficient, Kt    : " << kt << nl
        << "    Torque coefficient, 10*Kq : " << 10*kq << nl
        << "    Efficiency, etaO          : " << etaO << nl
        << nl;


    // Write state/results information
    setResult("n", n_);
    setResult("URef", URef);
    setResult("Kt", kt);
    setResult("Kq", kq);
    setResult("J", j);
    setResult("etaO", etaO);
}


void Foam::functionObjects::propellerInfoSed::writeWake
(
    const vectorField& U,
    const scalar URef
)
{
    if (!Pstream::master()) return;

    // Velocity
    auto& os = wakeFilePtr_();

    const pointField propPoints(coordSysPtr_->localPosition(points_));
    const label offset =
        mag(propPoints[1][0] - propPoints[0][0]) < SMALL ? 0 : 1;
    const scalar rMax = propPoints.last()[0];

    const scalar UzMean = meanSampleDiskField(U.component(2));

    writeHeaderValue(os, "Time", time_.timeOutputValue());
    writeHeaderValue(os, "Reference velocity", URef);
    writeHeaderValue(os, "Direction", coordSysPtr_->e3());
    writeHeaderValue(os, "Wake", 1 - UzMean/URef);
    writeHeader(os, "");
    writeCommented(os, "r/R");
    writeTabbed(os, "alpha");
    writeTabbed(os, "(x y z)");
    writeTabbed(os, "(Ur Utheta Uz)");
    os << nl;

    for (label thetai = 0; thetai < nTheta_; ++thetai)
    {
        const scalar deg = 360*thetai/scalar(nTheta_);

        for (label radiusi = 0; radiusi <= nRadial_; ++radiusi)
        {
            label pointi = radiusi*nTheta_ + thetai + offset;

            if (radiusi == 0 && offset == 1)
            {
                // Only a single point at the centre - repeat for all thetai
                pointi = 0;
            }

            if (pointMask_[pointi])
            {
                const scalar rR = propPoints[radiusi*nTheta_][0]/rMax;

                os  << rR << tab << deg << tab
                    << points_[pointi] << tab << U[pointi] << nl;
            }
        }
    }

    writeBreak(os);

    os  << endl;
}


void Foam::functionObjects::propellerInfoSed::writeAxialWake
(
    const vectorField& U,
    const scalar URef
)
{
    if (!Pstream::master()) return;

    // Alternative common format - axial wake component
    auto& os = axialWakeFilePtr_();

    const pointField propPoints(coordSysPtr_->localPosition(points_));
    const label offset =
        mag(propPoints[1][0] - propPoints[0][0]) < SMALL ? 0 : 1;
    const scalar rMax = propPoints.last()[0];

    writeHeaderValue(os, "Time", time_.timeOutputValue());

    os  << "# angle";
    for (label radiusi = 0; radiusi <= nRadial_; ++radiusi)
    {
        label pointi = radiusi*nTheta_;
        scalar r = propPoints[pointi][0];
        os  << tab << "r/R=" << r/rMax;
    }
    os  << nl;

    for (label thetai = 0; thetai < nTheta_; ++thetai)
    {
        os  << 360*thetai/scalar(nTheta_);

        for (label radiusi = 0; radiusi <= nRadial_; ++radiusi)
        {
            label pointi = radiusi*nTheta_ + thetai + offset;

            if (radiusi == 0 && offset == 1)
            {
                // Only a single point at the centre - repeat for all thetai
                pointi = 0;
            }

            if (pointMask_[pointi])
            {
                os << tab << 1 - U[pointi][2]/URef;
            }
            else
            {
                os << tab << "undefined";
            }
        }

        os  << nl;
    }

    writeBreak(os);

    os  << endl;
}


void Foam::functionObjects::propellerInfoSed::writeWakeFields(const scalar URef)
{
    if (!writeWakeFields_)
    {
        return;
    }

    // Normalised velocity
    const vectorField UDisk(interpolate(U(), vector::uniform(nanValue_))());
    const vectorField UrDisk(coordSysPtr_->localVector(UDisk));

    // Surface field data
    writeSampleDiskSurface(UDisk, UrDisk,  URef);

    // Write wake text files
    writeWake(UrDisk, URef);
    writeAxialWake(UrDisk, URef);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::functionObjects::propellerInfoSed::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    const Type& defaultValue
) const
{
    auto tfield = tmp<Field<Type>>::New(points_.size(), defaultValue);
    auto& field = tfield.ref();

    autoPtr<interpolation<Type>> interpolator
    (
        interpolation<Type>::New(interpolationScheme_, psi)
    );

    forAll(points_, pointi)
    {
        const label celli = cellIds_[pointi];

        if (cellIds_[pointi] != -1)
        {
            const point& position = points_[pointi];
            field[pointi] = interpolator().interpolate(position, celli);
        }
    }

    Pstream::listCombineGather(field, maxEqOp<Type>());
    Pstream::listCombineScatter(field);

    return tfield;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::propellerInfoSed::propellerInfoSed
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    forcesSed(name, runTime, dict, false),
    radius_(0),
    URefPtr_(nullptr),
    rotationMode_(rotationMode::SPECIFIED),
    MRFName_(),
    n_(0),
    writePropellerPerformance_(true),
    propellerPerformanceFilePtr_(nullptr),
    writeWakeFields_(true),
    surfaceWriterPtr_(nullptr),
    nTheta_(0),
    nRadial_(0),
    points_(),
    errorOnPointNotFound_(false),
    faces_(),
    cellIds_(),
    pointMask_(),
    interpolationScheme_("cell"),
    wakeFilePtr_(nullptr),
    axialWakeFilePtr_(nullptr),
    nanValue_(pTraits<scalar>::min)
{
    if (readFields)
    {
        read(dict);
        Log << endl;
    }
}


Foam::functionObjects::propellerInfoSed::propellerInfoSed
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    propellerInfoSed(name, obr.time(), dict, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::propellerInfoSed::read(const dictionary& dict)
{
    if (forcesSed::read(dict))
    {
        radius_ = dict.getScalar("radius");
        URefPtr_.reset(Function1<scalar>::New("URef", dict, &mesh_));
        rotationMode_ = rotationModeNames_.get("rotationMode", dict);

        // Must be set before setting the surface
        setCoordinateSystem(dict);

        writePropellerPerformance_ =
            dict.get<bool>("writePropellerPerformance");

        writeWakeFields_ = dict.get<bool>("writeWakeFields");
        if (writeWakeFields_)
        {
            setSampleDiskSurface(dict);

            dict.readIfPresent("interpolationScheme", interpolationScheme_);

            dict.readIfPresent("nanValue", nanValue_);
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::propellerInfoSed::execute()
{
    calcForcesMoment();

    createFiles();

    if (writeWakeFields_)
    {
        // Only setting mean axial velocity result during execute
        // - wake fields are 'heavy' and controlled separately using the
        //   writeControl
        const vectorField
            UDisk
            (
                coordSysPtr_->localVector
                (
                    interpolate
                    (
                        U(),
                        vector::uniform(nanValue_)
                    )()
                )
            );
        const scalar UzMean = meanSampleDiskField(UDisk.component(2));

        setResult("UzMean", UzMean);
    }

    writePropellerPerformance();

    return true;
}


bool Foam::functionObjects::propellerInfoSed::write()
{
    const scalar URef = URefPtr_->value(time_.timeOutputValue());
    writeWakeFields(URef);

    return true;
}


void Foam::functionObjects::propellerInfoSed::UpdateMesh(const mapPolyMesh& mpm)
{
    updateSampleDiskCells();
}


void Foam::functionObjects::propellerInfoSed::movePoints(const polyMesh& mesh)
{
    updateSampleDiskCells();
}


// ************************************************************************* //
