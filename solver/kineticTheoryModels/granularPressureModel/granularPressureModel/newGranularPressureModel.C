/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "granularPressureModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::granularPressureModel> Foam::granularPressureModel::New
(
    const dictionary& dict
)
{
    word granularPressureModelType(dict.get<word>("granularPressureModel"));

    Info<< "Selecting granularPressureModel "
        << granularPressureModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(granularPressureModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "granularPressureModel::New(const dictionary&) : " << endl
            << "    unknown granularPressureModelType type "
            << granularPressureModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid granularPressureModel types are :"  << endl
            <<    dictionaryConstructorTablePtr_->sortedToc() << endl;
        Info << abort(FatalError) << endl;
    }

    return autoPtr<granularPressureModel>(cstrIter()(dict));
}


// ************************************************************************* //
