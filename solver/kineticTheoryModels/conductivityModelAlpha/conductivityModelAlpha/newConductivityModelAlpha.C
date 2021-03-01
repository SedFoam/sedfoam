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

#include "conductivityModelAlpha.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::conductivityModelAlpha> Foam::conductivityModelAlpha::New
(
    const dictionary& dict
)
{
    word conductivityModelAlphaType(dict.lookup("conductivityModelAlpha"));

    Info<< "Selecting conductivityModelAlpha "
        << conductivityModelAlphaType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(conductivityModelAlphaType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "conductivityModelAlpha::New(const dictionary&) : " << endl
            << "    unknown conductivityModelAlphaType type "
            << conductivityModelAlphaType
            << ", constructor not in hash table" << endl << endl
            << "    Valid conductivityModelAlphaType types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<conductivityModelAlpha>(cstrIter()(dict));
}


// ************************************************************************* //
