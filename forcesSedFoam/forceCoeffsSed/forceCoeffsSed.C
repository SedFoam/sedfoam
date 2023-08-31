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

\*---------------------------------------------------------------------------*/

#include "forceCoeffsSed.H"
#include "dictionary.H"
#include "Time.H"
#include "Pstream.H"
#include "IOmanip.H"
#include "fvMesh.H"
#include "dimensionedTypes.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forceCoeffsSed, 0);
    addToRunTimeSelectionTable(functionObject, forceCoeffsSed, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forceCoeffsSed::createFiles()
{
    // Note: Only possible to create bin files after bins have been initialised

    if (writeToFile() && !coeffFilePtr_)
    {
        coeffFilePtr_ = createFile("coefficient");
        writeIntegratedHeader("Coefficients", coeffFilePtr_());

        if (nBin_ > 1)
        {
            CdBinFilePtr_ = createFile("CdBin");
            writeBinHeader("Drag coefficient bins", CdBinFilePtr_());

            CsBinFilePtr_ = createFile("CsBin");
            writeBinHeader("Side coefficient bins", CsBinFilePtr_());

            ClBinFilePtr_ = createFile("ClBin");
            writeBinHeader("Lift coefficient bins", ClBinFilePtr_());

            CmRollBinFilePtr_ = createFile("CmRollBin");
            writeBinHeader("Roll moment coefficient bins", CmRollBinFilePtr_());

            CmPitchBinFilePtr_ = createFile("CmPitchBin");
            writeBinHeader("Moment coefficient bins", CmPitchBinFilePtr_());

            CmYawBinFilePtr_ = createFile("CmYawBin");
            writeBinHeader("Yaw moment coefficient bins", CmYawBinFilePtr_());
        }
    }
}


void Foam::functionObjects::forceCoeffsSed::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    const auto& coordSys = coordSysPtr_();
    writeHeader(os, "Force coefficients in SedFoam version");
    writeHeaderValue(os, "dragDir", coordSys.e1());
    writeHeaderValue(os, "sideDir", coordSys.e2());
    writeHeaderValue(os, "liftDir", coordSys.e3());
    writeHeaderValue(os, "rollAxis", coordSys.e1());
    writeHeaderValue(os, "pitchAxis", coordSys.e2());
    writeHeaderValue(os, "yawAxis", coordSys.e3());
    writeHeaderValue(os, "magUInf", magUInf_);
    writeHeaderValue(os, "lRef", lRef_);
    writeHeaderValue(os, "Aref", Aref_);
    writeHeaderValue(os, "CofR", coordSys.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, "Cd");
    writeTabbed(os, "Cs");
    writeTabbed(os, "Cl");
    writeTabbed(os, "CmRoll");
    writeTabbed(os, "CmPitch");
    writeTabbed(os, "CmYaw");
    writeTabbed(os, "Cd(f)");
    writeTabbed(os, "Cd(r)");
    writeTabbed(os, "Cs(f)");
    writeTabbed(os, "Cs(r)");
    writeTabbed(os, "Cl(f)");
    writeTabbed(os, "Cl(r)");
    os  << endl;
}


void Foam::functionObjects::forceCoeffsSed::writeBinHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeHeaderValue(os, "bins", nBin_);
    writeHeaderValue(os, "start", binMin_);
    writeHeaderValue(os, "delta", binDx_);
    writeHeaderValue(os, "direction", binDir_);

    vectorField binPoints(nBin_);
    writeCommented(os, "x co-ords  :");
    forAll(binPoints, pointi)
    {
        binPoints[pointi] = (binMin_ + (pointi + 1)*binDx_)*binDir_;
        os  << tab << binPoints[pointi].x();
    }
    os  << nl;

    writeCommented(os, "y co-ords  :");
    forAll(binPoints, pointi)
    {
        os  << tab << binPoints[pointi].y();
    }
    os  << nl;

    writeCommented(os, "z co-ords  :");
    forAll(binPoints, pointi)
    {
        os  << tab << binPoints[pointi].z();
    }
    os  << nl;

    writeHeader(os, "");
    writeCommented(os, "Time");

    for (label j = 0; j < nBin_; ++j)
    {
        const word jn(Foam::name(j) + ':');
        writeTabbed(os, jn + "total");
        writeTabbed(os, jn + "pressure");
        writeTabbed(os, jn + "viscous");


    }

    os  << endl;
}


void Foam::functionObjects::forceCoeffsSed::writeIntegratedData
(
    const word& title,
    const List<Field<scalar>>& coeff
) const
{
    if (!log)
    {
        return;
    }

    const scalar pressure = sum(coeff[0]);
    const scalar viscous = sum(coeff[1]);
    const scalar total = pressure + viscous;

    Info<< "        " << title << "       : " << total << token::TAB
        << '('
        << "pressure: " << pressure << token::TAB
        << "viscous: " << viscous;


    Info<< ')' << endl;
}


void Foam::functionObjects::forceCoeffsSed::writeBinData
(
    const List<Field<scalar>> coeffs,
    Ostream& os
) const
{
    writeCurrentTime(os);

    for (label bini = 0; bini < nBin_; ++bini)
    {
        scalar total = coeffs[0][bini] + coeffs[1][bini] + coeffs[2][bini];

        os  << tab << total << tab << coeffs[0][bini] << tab << coeffs[1][bini];


    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffsSed::forceCoeffsSed
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool readFields
)
:
    forcesSed(name, runTime, dict, false),
    magUInf_(Zero),
    lRef_(Zero),
    Aref_(Zero),
    coeffFilePtr_(),
    CdBinFilePtr_(),
    CsBinFilePtr_(),
    ClBinFilePtr_(),
    CmRollBinFilePtr_(),
    CmPitchBinFilePtr_(),
    CmYawBinFilePtr_()
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict, "liftDir", "dragDir");
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forceCoeffsSed::read(const dictionary& dict)
{
    forcesSed::read(dict);

    // Free stream velocity magnitude
    dict.readEntry("magUInf", magUInf_);

    // If case is compressible we must read rhoInf (store in rhoRef_) to
    // calculate the reference dynamic pressure
    // Note: for incompressible, rhoRef_ is already initialised
    if (rhoName_ != "rhoInf")
    {
        dict.readEntry("rhoInf", rhoRef_);
    }

    // Reference length and area scales
    dict.readEntry("lRef", lRef_);
    dict.readEntry("Aref", Aref_);

    if (writeFields_)
    {
        volVectorField* forceCoeffPtr
        (
            new volVectorField
            (
                IOobject
                (
                    scopedName("forceCoeff"),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector(dimless, Zero)
            )
        );

        mesh_.objectRegistry::store(forceCoeffPtr);

        volVectorField* momentCoeffPtr
        (
            new volVectorField
            (
                IOobject
                (
                    scopedName("momentCoeff"),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector(dimless, Zero)
            )
        );

        mesh_.objectRegistry::store(momentCoeffPtr);
    }

    return true;
}


bool Foam::functionObjects::forceCoeffsSed::execute()
{
    forcesSed::calcForcesMoment();

    createFiles();

    // Storage for pressure and viscous contributions to coeffs
    List<Field<scalar>> dragCoeffs(3);
    List<Field<scalar>> sideCoeffs(3);
    List<Field<scalar>> liftCoeffs(3);
    List<Field<scalar>> rollMomentCoeffs(3);
    List<Field<scalar>> pitchMomentCoeffs(3);
    List<Field<scalar>> yawMomentCoeffs(3);

    forAll(liftCoeffs, i)
    {
        dragCoeffs[i].setSize(nBin_);
        sideCoeffs[i].setSize(nBin_);
        liftCoeffs[i].setSize(nBin_);
        rollMomentCoeffs[i].setSize(nBin_);
        pitchMomentCoeffs[i].setSize(nBin_);
        yawMomentCoeffs[i].setSize(nBin_);
    }

    // Calculate coefficients
    scalar CdTot = 0;
    scalar CsTot = 0;
    scalar ClTot = 0;
    scalar CmRollTot = 0;
    scalar CmPitchTot = 0;
    scalar CmYawTot = 0;

    const scalar pDyn = 0.5*rhoRef_*sqr(magUInf_);

    // Avoid divide by zero in 2D cases
    const scalar momentScaling = 1.0/(Aref_*pDyn*lRef_ + SMALL);
    const scalar forceScaling = 1.0/(Aref_*pDyn + SMALL);

    const auto& coordSys = coordSysPtr_();

    forAll(liftCoeffs, i)
    {
        const Field<vector> localForce(coordSys.localVector(forceSed_[i]));
        const Field<vector> localMoment(coordSys.localVector(moment_[i]));

        dragCoeffs[i] = forceScaling*(localForce.component(0));
        sideCoeffs[i] = forceScaling*(localForce.component(1));
        liftCoeffs[i] = forceScaling*(localForce.component(2));
        rollMomentCoeffs[i] = momentScaling*(localMoment.component(0));
        pitchMomentCoeffs[i] = momentScaling*(localMoment.component(1));
        yawMomentCoeffs[i] = momentScaling*(localMoment.component(2));

        CdTot += sum(dragCoeffs[i]);
        CsTot += sum(sideCoeffs[i]);
        ClTot += sum(liftCoeffs[i]);
        CmRollTot += sum(rollMomentCoeffs[i]);
        CmPitchTot += sum(pitchMomentCoeffs[i]);
        CmYawTot += sum(yawMomentCoeffs[i]);
    }

    // Single contributions to the front and rear
    const scalar CdfTot = 0.5*CdTot + CmRollTot;
    const scalar CdrTot = 0.5*CdTot - CmRollTot;
    const scalar CsfTot = 0.5*CsTot + CmYawTot;
    const scalar CsrTot = 0.5*CsTot - CmYawTot;
    const scalar ClfTot = 0.5*ClTot + CmPitchTot;
    const scalar ClrTot = 0.5*ClTot - CmPitchTot;

    Log << type() << " " << name() << " execute:" << nl
        << "    Coefficients" << nl;

    writeIntegratedData("Cd", dragCoeffs);
    writeIntegratedData("Cs", sideCoeffs);
    writeIntegratedData("Cl", liftCoeffs);
    writeIntegratedData("CmRoll", rollMomentCoeffs);
    writeIntegratedData("CmPitch", pitchMomentCoeffs);
    writeIntegratedData("CmYaw", yawMomentCoeffs);

    Log << "        Cd(f)    : " << CdfTot << nl
        << "        Cd(r)    : " << CdrTot << nl;

    Log << "        Cs(f)    : " << CsfTot << nl
        << "        Cs(r)    : " << CsrTot << nl;

    Log << "        Cl(f)    : " << ClfTot << nl
        << "        Cl(r)    : " << ClrTot << nl;

    if (writeToFile())
    {
        writeCurrentTime(coeffFilePtr_());
        coeffFilePtr_()
            << tab << CdTot << tab << CsTot << tab << ClTot
            << tab << CmRollTot << tab << CmPitchTot << tab << CmYawTot
            << tab << CdfTot << tab << CdrTot
            << tab << CsfTot << tab << CsrTot
            << tab << ClfTot << tab << ClrTot << endl;

        if (nBin_ > 1)
        {
            if (binCumulative_)
            {
                forAll(liftCoeffs, i)
                {
                    for (label bini = 1; bini < nBin_; ++bini)
                    {
                        dragCoeffs[i][bini] += dragCoeffs[i][bini-1];
                        sideCoeffs[i][bini] += sideCoeffs[i][bini-1];
                        liftCoeffs[i][bini] += liftCoeffs[i][bini-1];
                        rollMomentCoeffs[i][bini] +=
                            rollMomentCoeffs[i][bini-1];
                        pitchMomentCoeffs[i][bini] +=
                            pitchMomentCoeffs[i][bini-1];
                        yawMomentCoeffs[i][bini] += yawMomentCoeffs[i][bini-1];
                    }
                }
            }

            writeBinData(dragCoeffs, CdBinFilePtr_());
            writeBinData(sideCoeffs, CsBinFilePtr_());
            writeBinData(liftCoeffs, ClBinFilePtr_());
            writeBinData(rollMomentCoeffs, CmRollBinFilePtr_());
            writeBinData(pitchMomentCoeffs, CmPitchBinFilePtr_());
            writeBinData(yawMomentCoeffs, CmYawBinFilePtr_());
        }
    }

    // Write state/results information
    {
        setResult("Cd", CdTot);
        setResult("Cs", CsTot);
        setResult("Cl", ClTot);
        setResult("CmRoll", CmRollTot);
        setResult("CmPitch", CmPitchTot);
        setResult("CmYaw", CmYawTot);
        setResult("Cd(f)", CdfTot);
        setResult("Cd(r)", CdrTot);
        setResult("Cs(f)", CsfTot);
        setResult("Cs(r)", CsrTot);
        setResult("Cl(f)", ClfTot);
        setResult("Cl(r)", ClrTot);
    }

    if (writeFields_)
    {
        const volVectorField& forceSed =
            lookupObject<volVectorField>(scopedName("forceSed"));

        const volVectorField& moment =
            lookupObject<volVectorField>(scopedName("moment"));

        volVectorField& forceCoeff =
            lookupObjectRef<volVectorField>(scopedName("forceCoeff"));

        volVectorField& momentCoeff =
            lookupObjectRef<volVectorField>(scopedName("momentCoeff"));

        dimensionedScalar f0("f0", dimForce, Aref_*pDyn);
        dimensionedScalar m0("m0", dimForce*dimLength, Aref_*lRef_*pDyn);

        forceCoeff == forceSed/f0;
        momentCoeff == moment/m0;
    }

    return true;
}


bool Foam::functionObjects::forceCoeffsSed::write()
{
    if (writeFields_)
    {
        const volVectorField& forceCoeff =
            lookupObject<volVectorField>(scopedName("forceCoeff"));

        const volVectorField& momentCoeff =
            lookupObject<volVectorField>(scopedName("momentCoeff"));

        forceCoeff.write();
        momentCoeff.write();
    }

    return true;
}


// ************************************************************************* //
