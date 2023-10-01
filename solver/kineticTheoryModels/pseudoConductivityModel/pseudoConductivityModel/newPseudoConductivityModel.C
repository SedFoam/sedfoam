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

#include "pseudoConductivityModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::pseudoConductivityModel> Foam::pseudoConductivityModel::New
(
    const dictionary& dict
)
{
    word pseudoConductivityModelType
            (dict.getOrDefault<word>("pseudoConductivityModel", "none"));

    Info<< "Selecting pseudoConductivityModel "
        << pseudoConductivityModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(pseudoConductivityModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "pseudoConductivityModel::New(const dictionary&) : " << endl
            << "    unknown pseudoConductivityModelType type "
            << pseudoConductivityModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid pseudoConductivityModelType types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<pseudoConductivityModel>(cstrIter()(dict));
}


// ************************************************************************* //
