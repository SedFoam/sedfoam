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

#include "forcesSed.H"
#include "fvcGrad.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "cartesianCS.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forcesSed, 0);
    addToRunTimeSelectionTable(functionObject, forcesSed, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forcesSed::createFiles()
{
    // Note: Only possible to create bin files after bins have been initialised

    if (writeToFile() && !forceFilePtr_)
    {
        forceFilePtr_ = createFile("forcesSed");
        writeIntegratedHeader("Force", forceFilePtr_());
        momentFilePtr_ = createFile("moment");
        writeIntegratedHeader("Moment", momentFilePtr_());

        if (nBin_ > 1)
        {
            forceBinFilePtr_ = createFile("forceBin");
            writeBinHeader("Force", forceBinFilePtr_());
            momentBinFilePtr_ = createFile("momentBin");
            writeBinHeader("Moment", momentBinFilePtr_());
        }
    }
}


void Foam::functionObjects::forcesSed::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeHeaderValue(os, "CofR", coordSysPtr_->origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, "(total_x total_y total_z)");
    writeTabbed(os, "(pressureFluid_x pressureFluid_y pressureFluid_z)");
    writeTabbed(os, "(pressureParticle_x pressureParticle_y pressureParticle_z)");
    writeTabbed(os, "(viscous_x viscous_y viscous_z)");


    os  << endl;
}


void Foam::functionObjects::forcesSed::writeBinHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header + " bins");
    writeHeaderValue(os, "bins", nBin_);
    writeHeaderValue(os, "start", binMin_);
    writeHeaderValue(os, "end", binMax_);
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

    for (label j = 0; j < nBin_; j++)
    {
        const word jn(Foam::name(j) + ':');
        os  << tab << jn << "(total_x total_y total_z)"
            << tab << jn << "(pressureFluid_x pressureFluid_y pressureFluid_z)"
            << tab << jn << "(pressureParticle_x pressureFluid_y pressureFluid_z)"
            << tab << jn << "(viscous_x viscous_y viscous_z)";

    }

    os << endl;
}


void Foam::functionObjects::forcesSed::setCoordinateSystem
(
    const dictionary& dict,
    const word& e3Name,
    const word& e1Name
)
{
    coordSysPtr_.clear();

    point origin(Zero);
    if (dict.readIfPresent<point>("CofR", origin))
    {
        const vector e3 = e3Name == word::null ?
            vector(0, 0, 1) : dict.get<vector>(e3Name);
        const vector e1 = e1Name == word::null ?
            vector(1, 0, 0) : dict.get<vector>(e1Name);

        coordSysPtr_.reset(new coordSystem::cartesian(origin, e3, e1));
    }
    else
    {
        // The 'coordinateSystem' sub-dictionary is optional,
        // but enforce use of a cartesian system if not found.

        if (dict.found(coordinateSystem::typeName_()))
        {
            // New() for access to indirect (global) coordinate system
            coordSysPtr_ =
                coordinateSystem::New
                (
                    obr_,
                    dict,
                    coordinateSystem::typeName_()
                );
        }
        else
        {
            coordSysPtr_.reset(new coordSystem::cartesian(dict));
        }
    }

}


void Foam::functionObjects::forcesSed::initialise()
{
    if (initialised_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!foundObject<volVectorField>(fDName_))
        {
            FatalErrorInFunction
                << "Could not find " << fDName_ << " in database"
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            !foundObject<volVectorField>(UaName_)
         || !foundObject<volScalarField>(pName_)

        )
        {
            FatalErrorInFunction
                << "Could not find Ua: " << UaName_ << " or p:" << pName_
                << " in database"
                << exit(FatalError);
        }

        if (rhoName_ != "rhoInf" && !foundObject<volScalarField>(rhoName_))
        {
            FatalErrorInFunction
                << "Could not find rho:" << rhoName_
                << exit(FatalError);
        }
    }

    initialiseBins();

    initialised_ = true;
}


void Foam::functionObjects::forcesSed::initialiseBins()
{
    if (nBin_ > 1)
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Determine extents of patches
        scalar geomMin = GREAT;
        scalar geomMax = -GREAT;
        for (const label patchi : patchSet_)
        {
            const polyPatch& pp = pbm[patchi];
            scalarField d(pp.faceCentres() & binDir_);
            geomMin = min(min(d), geomMin);
            geomMax = max(max(d), geomMax);
        }



        reduce(geomMin, minOp<scalar>());
        reduce(geomMax, maxOp<scalar>());

        // Slightly boost max so that region of interest is fully within bounds
        geomMax = 1.0001*(geomMax - geomMin) + geomMin;

        // Use geometry limits if not specified by the user
        if (binMin_ == GREAT)
        {
            binMin_ = geomMin;
        }
        if (binMax_ == GREAT)
        {
            binMax_ = geomMax;
        }

        binDx_ = (binMax_ - binMin_)/scalar(nBin_);

        // Create the bin mid-points used for writing
        binPoints_.setSize(nBin_);
        forAll(binPoints_, i)
        {
            binPoints_[i] = (i + 0.5)*binDir_*binDx_;
        }
    }

    // Allocate storage for forcesSed and moments
    forAll(forceSed_, i)
    {
        forceSed_[i].setSize(nBin_, vector::zero);
        moment_[i].setSize(nBin_, vector::zero);
    }
}


void Foam::functionObjects::forcesSed::resetFields()
{
    forceSed_[0] = Zero;
    forceSed_[1] = Zero;
    forceSed_[2] = Zero;
    forceSed_[3] = Zero;

    moment_[0] = Zero;
    moment_[1] = Zero;
    moment_[2] = Zero;
    moment_[3] = Zero;

    if (writeFields_)
    {
        volVectorField& forceSed =
            lookupObjectRef<volVectorField>(scopedName("forceSed"));

        forceSed == dimensionedVector(forceSed.dimensions(), Zero);

        volVectorField& moment =
            lookupObjectRef<volVectorField>(scopedName("moment"));

        moment == dimensionedVector(moment.dimensions(), Zero);
    }
}


Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::forcesSed::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    const auto& Ua = lookupObject<volVectorField>(UaName_);
    const auto& Ub = lookupObject<volVectorField>(UbName_);


    if (foundObject<dictionary>("transportProperties"))
    {
        //dimensionedScalar nu("nu", dimViscosity, transportProperties);
        //return -rho()*nu*dev(twoSymm(fvc::grad(U)));
        
        const volScalarField& alpha = lookupObject<volScalarField>(alphaName_);
        const volScalarField& muEff = lookupObject<volScalarField>(muEffName_);
        const volScalarField& muFra = lookupObject<volScalarField>(muFraName_);
        return -muEff*dev(twoSymm(fvc::grad(Ub)))-muFra*dev(twoSymm(fvc::grad(Ua)));
        //return -(muEff+muFra)*dev(twoSymm(fvc::grad(Ub*(1-alpha)+Ua*alpha)));
       // return -muEff*dev(twoSymm(fvc::grad(Ub)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forcesSed::mu() const
{

    if (foundObject<dictionary>("transportProperties"))
    {
		
        const volScalarField& muEff = lookupObject<volScalarField>(muEffName_);
        return muEff;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forcesSed::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<volScalarField>::New
        (
            IOobject
            (
                "rho",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("rho", dimDensity, rhoRef_)
        );
    }

    return(lookupObject<volScalarField>(rhoName_));
}


Foam::scalar Foam::functionObjects::forcesSed::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1.0;
    }

    if (rhoName_ != "rhoInf")
    {
        FatalErrorInFunction
            << "Dynamic pressure is expected but kinematic is provided."
            << exit(FatalError);
    }

    return rhoRef_;
}


void Foam::functionObjects::forcesSed::applyBins
(
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fNsolid,
    const vectorField& fT,
    const vectorField& fP,
    const vectorField& d
)
{
    if (nBin_ == 1)
    {
        forceSed_[0][0] += sum(fN);
        forceSed_[1][0] += sum(fNsolid);
        forceSed_[2][0] += sum(fT);
        forceSed_[3][0] += sum(fP);
        moment_[0][0] += sum(Md^fN);
        moment_[1][0] += sum(Md^fNsolid);
        moment_[2][0] += sum(Md^fT);
        moment_[3][0] += sum(Md^fP);
        
    }
    else
    {
        scalarField dd((d & binDir_) - binMin_);

        forAll(dd, i)
        {
            label bini = min(max(floor(dd[i]/binDx_), 0), forceSed_[0].size() - 1);

            forceSed_[0][bini] += fN[i];
            forceSed_[1][bini] += fNsolid[i];
            forceSed_[2][bini] += fT[i];
            forceSed_[3][bini] += fP[i];
            moment_[0][bini] += Md[i]^fN[i];
            moment_[1][bini] += Md[i]^fNsolid[i];
            moment_[2][bini] += Md[i]^fT[i];
            moment_[3][bini] += Md[i]^fP[i];
        }
    }
}


void Foam::functionObjects::forcesSed::addToFields
(
    const label patchi,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fNsolid,
    const vectorField& fT,
    const vectorField& fP
)
{
    if (!writeFields_)
    {
        return;
    }

    auto& forceSed = lookupObjectRef<volVectorField>(scopedName("forceSed"));
    vectorField& pf = forceSed.boundaryFieldRef()[patchi];
    pf += fN+fNsolid + fT + fP;

    auto& moment = lookupObjectRef<volVectorField>(scopedName("moment"));
    vectorField& pm = moment.boundaryFieldRef()[patchi];
    pm = Md^pf;
}


void Foam::functionObjects::forcesSed::addToFields
(
    const labelList& cellIDs,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fNsolid,
    const vectorField& fT,
    const vectorField& fP
)
{
    if (!writeFields_)
    {
        return;
    }

    auto& forceSed = lookupObjectRef<volVectorField>(scopedName("forcesSed"));
    auto& moment = lookupObjectRef<volVectorField>(scopedName("moment"));

    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];
        forceSed[celli] += fN[i] + fNsolid[i] + fT[i] + fP[i];
        moment[celli] = Md[i]^forceSed[celli];
    }
}


void Foam::functionObjects::forcesSed::writeIntegratedForceMoment
(
    const string& descriptor,
    const vectorField& fm0,
    const vectorField& fm1,
    const vectorField& fm2,
    const vectorField& fm3,
    autoPtr<OFstream>& osPtr
) const
{
    vector pressureFluid = sum(fm0);
    vector pressureParticle = sum(fm1);
    vector viscous = sum(fm2);
    vector porous = sum(fm3);
    vector total = pressureFluid+pressureParticle + viscous + porous;

    Log << "   SedFoam Version. Sum of " << descriptor.c_str() << nl
        << "        Total    : " << total << nl
        << "        pressureFluid : " << pressureFluid << nl
        << "        pressureParticle : " << pressureParticle << nl
        << "        Viscous  : " << viscous << nl;


    if (writeToFile())
    {
        Ostream& os = osPtr();

        writeCurrentTime(os);

        os  << tab << total
            << tab << pressureFluid
            << tab << pressureParticle
            << tab << viscous;



        os  << endl;
    }
}


void Foam::functionObjects::forcesSed::writeForces()
{
    Log << type() << " " << name() << " write:" << nl;

    const auto& coordSys = coordSysPtr_();

    writeIntegratedForceMoment
    (
        "forcesSed",
        coordSys.localVector(forceSed_[0]),
        coordSys.localVector(forceSed_[1]),
        coordSys.localVector(forceSed_[2]),
        coordSys.localVector(forceSed_[3]),
        forceFilePtr_
    );

    writeIntegratedForceMoment
    (
        "moments",
        coordSys.localVector(moment_[0]),
        coordSys.localVector(moment_[1]),
        coordSys.localVector(moment_[2]),
        coordSys.localVector(moment_[3]),
        momentFilePtr_
    );

    Log << endl;
}


void Foam::functionObjects::forcesSed::writeBinnedForceMoment
(
    const List<Field<vector>>& fm,
    autoPtr<OFstream>& osPtr
) const
{
    if ((nBin_ == 1) || !writeToFile())
    {
        return;
    }

    List<Field<vector>> f(fm);

    if (binCumulative_)
    {
        for (label i = 1; i < f[0].size(); i++)
        {
            f[0][i] += f[0][i-1];
            f[1][i] += f[1][i-1];
            f[2][i] += f[2][i-1];
            f[3][i] += f[3][i-1];

        }
    }

    Ostream& os = osPtr();

    writeCurrentTime(os);

    forAll(f[0], i)
    {
        vector total = f[0][i] + f[1][i] + f[2][i]+ f[3][i];

        os  << tab << total
            << tab << f[0][i]
            << tab << f[1][i]
            << tab << f[2][i];

    }

    os  << nl;
}


void Foam::functionObjects::forcesSed::writeBins()
{
    const auto& coordSys = coordSysPtr_();

    List<Field<vector>> lf(4);
    List<Field<vector>> lm(4);
    lf[0] = coordSys.localVector(forceSed_[0]);
    lf[1] = coordSys.localVector(forceSed_[1]);
    lf[2] = coordSys.localVector(forceSed_[2]);
    lf[3] = coordSys.localVector(forceSed_[3]);

    lm[0] = coordSys.localVector(moment_[0]);
    lm[1] = coordSys.localVector(moment_[1]);
    lm[2] = coordSys.localVector(moment_[2]);
    lm[3] = coordSys.localVector(moment_[3]);

    writeBinnedForceMoment(lf, forceBinFilePtr_);
    writeBinnedForceMoment(lm, momentBinFilePtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forcesSed::forcesSed
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    forceSed_(4),
    moment_(4),
    forceFilePtr_(),
    momentFilePtr_(),
    forceBinFilePtr_(),
    momentBinFilePtr_(),
    patchSet_(),
    pName_("p"),
    prbghName_("p_rbgh"),
    pSedName_("pS"),
    alphaName_("alpha.a"),
    muEffName_("muEff"),
    muFraName_("muFra"),
    UaName_("U.a"),
    UbName_("U.b"),
    rhoName_("rho"),
    directForceDensity_(false),
    fDName_("fD"),
    rhoRef_(VGREAT),
    pRef_(0),
    coordSysPtr_(nullptr),
    nBin_(1),
    binDir_(Zero),
    binDx_(0),
    binMin_(GREAT),
    binMax_(GREAT),
    binPoints_(),
    binCumulative_(true),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


Foam::functionObjects::forcesSed::forcesSed
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name),
    forceSed_(4),
    moment_(4),
    forceFilePtr_(),
    momentFilePtr_(),
    forceBinFilePtr_(),
    momentBinFilePtr_(),
    patchSet_(),
    pName_("p"),
    prbghName_("p_rbgh"),
    pSedName_("pS"),
    alphaName_("alpha.a"),
    muEffName_("muEff"),
    muFraName_("muFra"),
    UaName_("U.a"),
    UbName_("U.b"),
    rhoName_("rho"),
    directForceDensity_(false),
    fDName_("fD"),
    rhoRef_(VGREAT),
    pRef_(0),
    coordSysPtr_(nullptr),
    nBin_(1),
    binDir_(Zero),
    binDx_(0),
    binMin_(GREAT),
    binMax_(GREAT),
    binPoints_(),
    binCumulative_(true),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }

/*
    // Turn off writing to file
    writeToFile_ = false;

    forAll(forceSed_, i)
    {
        forceSed_[i].setSize(nBin_, vector::zero);
        moment_[i].setSize(nBin_, vector::zero);
    }
*/
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forcesSed::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    initialised_ = false;

    Info<< type() << " " << name() << ":" << nl;

    directForceDensity_ = dict.getOrDefault("directForceDensity", false);

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        );

    if (directForceDensity_)
    {
        // Optional entry for fDName
        if (dict.readIfPresent<word>("fD", fDName_))
        {
            Info<< "    fD: " << fDName_ << endl;
        }
    }
    else
    {
        // Optional field name entries
        if (dict.readIfPresent<word>("p", pName_))
        {
            Info<< "    p: " << pName_ << endl;
        }
        if (dict.readIfPresent<word>("p_rbgh", prbghName_))
        {
            Info<< "    p_rbgh: " << prbghName_ << endl;
        }
        if (dict.readIfPresent<word>("pS", pSedName_))
        {
            Info<< "    pS: " << pSedName_ << endl;
        }
        if (dict.readIfPresent<word>("alpha.a", alphaName_))
        {
            Info<< "    alpha: " << alphaName_ << endl;
        }
        if (dict.readIfPresent<word>("muEff", muEffName_))
        {
            Info<< "    muEff: " << muEffName_ << endl;
        }
        if (dict.readIfPresent<word>("muFra", muFraName_))
        {
            Info<< "    muFra: " << muFraName_ << endl;
        }
        if (dict.readIfPresent<word>("U.a", UaName_))
        {
            Info<< "    Ua: " << UaName_ << endl;
        }
        if (dict.readIfPresent<word>("U.b", UbName_))
        {
            Info<< "    Ub: " << UbName_ << endl;
        }
        if (dict.readIfPresent<word>("rho", rhoName_))
        {
            Info<< "    rho: " << rhoName_ << endl;
        }

        // Reference density needed for incompressible calculations
        if (rhoName_ == "rhoInf")
        {
            rhoRef_ = dict.get<scalar>("rhoInf");
            Info<< "    Freestream density (rhoInf) set to " << rhoRef_ << endl;
        }

        // Reference pressure, 0 by default
        if (dict.readIfPresent<scalar>("pRef", pRef_))
        {
            Info<< "    Reference pressure (pRef) set to " << pRef_ << endl;
        }
    }


    if (dict.found("binData"))
    {
        Info<< "    Activated data bins" << endl;
        const dictionary& binDict(dict.subDict("binData"));
        nBin_ = binDict.get<label>("nBin");

        if (nBin_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "Number of bins (nBin) must be zero or greater"
                << exit(FatalIOError);
        }
        else if (nBin_ == 0)
        {
            // Case of no bins equates to a single bin to collect all data
            nBin_ = 1;
        }
        else
        {
            Info<< "    Employing " << nBin_ << " bins" << endl;
            if (binDict.readIfPresent("min", binMin_))
            {
                Info<< "    - min         : " << binMin_ << endl;
            }
            if (binDict.readIfPresent("max", binMax_))
            {
                Info<< "    - max         : " << binMax_ << endl;
            }

            binCumulative_ = binDict.get<bool>("cumulative");
            Info<< "    - cumuluative : " << binCumulative_ << endl;

            binDir_ = binDict.get<vector>("direction");
            binDir_.normalise();
            Info<< "    - direction   : " << binDir_ << endl;
        }
    }

    writeFields_ = dict.getOrDefault("writeFields", false);

    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;

        volVectorField* forcePtr
        (
            new volVectorField
            (
                IOobject
                (
                    scopedName("forceSed"),
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector(dimForce, Zero)
            )
        );

        mesh_.objectRegistry::store(forcePtr);

        volVectorField* momentPtr
        (
            new volVectorField
            (
                IOobject
                (
                    scopedName("moment"),
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector(dimForce*dimLength, Zero)
            )
        );

        mesh_.objectRegistry::store(momentPtr);
    }

    return true;
}


void Foam::functionObjects::forcesSed::calcForcesMoment()
{
    initialise();

    resetFields();

    const point& origin = coordSysPtr_->origin();

    if (directForceDensity_)
    {
        const volVectorField& fD = lookupObject<volVectorField>(fDName_);

        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        for (const label patchi : patchSet_)
        {
            vectorField Md(mesh_.C().boundaryField()[patchi] - origin);

            scalarField sA(mag(Sfb[patchi]));

            // Normal forceSed = surfaceUnitNormal*(surfaceNormal & forceDensity)
            vectorField fN
            (
                Sfb[patchi]/sA
               *(
                    Sfb[patchi] & fD.boundaryField()[patchi]
                )
            );
            vectorField fNsolid
            (
                Sfb[patchi]/sA
               *(
                    Sfb[patchi] & fD.boundaryField()[patchi]
                )
            );
            // Tangential forceSed (total forceSed minus normal fN)
            vectorField fT(sA*fD.boundaryField()[patchi] - fN);

            // Porous forceSed
            vectorField fP(Md.size(), Zero);

            addToFields(patchi, Md, fN, fNsolid, fT, fP);

            applyBins(Md, fN, fNsolid, fT, fP, mesh_.C().boundaryField()[patchi]);
        }
    }
    else
    {
        const volScalarField& p = lookupObject<volScalarField>(pName_);
        const volScalarField& p_rbgh = lookupObject<volScalarField>(prbghName_);
        const volScalarField& pS = lookupObject<volScalarField>(pSedName_);
        const volScalarField& alpha = lookupObject<volScalarField>(alphaName_);


        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
        const volSymmTensorField::Boundary& devRhoReffb
            = tdevRhoReff().boundaryField();

        // Scale pRef by density for incompressible simulations
        scalar pRef = pRef_/rho(p);

        for (const label patchi : patchSet_)
        {
            vectorField Md(mesh_.C().boundaryField()[patchi] - origin);




            vectorField fN
            (
                rho(p)*Sfb[patchi]*(1.*p_rbgh.boundaryField()[patchi]- pRef)
                
            );
            vectorField fNsolid
            (
                rho(p)*Sfb[patchi]*(1*pS.boundaryField()[patchi] - pRef)
            );
            
           
            
            vectorField fT(Sfb[patchi] & devRhoReffb[patchi]);

            vectorField fP(Md.size(), Zero);

            addToFields(patchi, Md, fN,fNsolid, fT, fP);

            applyBins(Md, fN,fNsolid, fT, fP, mesh_.C().boundaryField()[patchi]);
        }
    }

   

    Pstream::listCombineGather(forceSed_, plusEqOp<vectorField>());
    Pstream::listCombineGather(moment_, plusEqOp<vectorField>());
    Pstream::listCombineScatter(forceSed_);
    Pstream::listCombineScatter(moment_);
}


Foam::vector Foam::functionObjects::forcesSed::forceEff() const
{
    return sum(forceSed_[0]) + sum(forceSed_[1]) + sum(forceSed_[2])+ sum(forceSed_[3]);
}


Foam::vector Foam::functionObjects::forcesSed::momentEff() const
{
    return sum(moment_[0]) + sum(moment_[1]) + sum(moment_[2])+ sum(moment_[3]);
}


bool Foam::functionObjects::forcesSed::execute()
{
    calcForcesMoment();

    if (Pstream::master())
    {
        createFiles();

        writeForces();

        writeBins();

        Log << endl;
    }

    // Write state/results information
    setResult("normalForce", sum(forceSed_[0]+forceSed_[1]));
    setResult("tangentialForce", sum(forceSed_[2]));
    setResult("porousForce", sum(forceSed_[3]));

    setResult("normalMoment", sum(moment_[0]+moment_[1]));
    setResult("tangentialMoment", sum(moment_[2]));
    setResult("porousMoment", sum(moment_[3]));

    return true;
}


bool Foam::functionObjects::forcesSed::write()
{
    if (writeFields_)
    {
        lookupObject<volVectorField>(scopedName("forceSed")).write();
        lookupObject<volVectorField>(scopedName("moment")).write();
    }

    return true;
}


// ************************************************************************* //
