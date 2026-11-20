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

#include "vegDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::vegDragModel> Foam::vegDragModel::New
(
    const dictionary& interfaceDict,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
{
    word vegDragModelType
    (
        interfaceDict.get<word>("vegDragModel" + phasea.name())
    );

    Info << "Selecting vegDragModel for phase "
        << phasea.name()
        << ": "
        << vegDragModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(vegDragModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "vegDragModel::New : " << endl
                << "    unknown vegDragModelType type "
                << vegDragModelType
                << ", constructor not in hash table" << endl << endl
                << "    Valid vegDragModel types are : " << endl;
        Info << dictionaryConstructorTablePtr_->sortedToc()
             << abort(FatalError);
    }

    return cstrIter()(interfaceDict, phasea, phaseb);
}


// ************************************************************************* //
