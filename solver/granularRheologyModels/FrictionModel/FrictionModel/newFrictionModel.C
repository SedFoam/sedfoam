/*---------------------------------------------------------------------------*\
Copyright (C) 2015 Cyrille Bonamy, Julien Chauchat, Tian-Jian Hsu
                   and contributors

License
    This file is part of SedFOAM.

    SedFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SedFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SedFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FrictionModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::granularRheologyModels::FrictionModel>
Foam::granularRheologyModels::FrictionModel::New
(
    const dictionary& dict
)
{
    word FrictionModelType(dict.get<word>("FrictionModel"));

    Info<< "Selecting FrictionModel "
        << FrictionModelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(FrictionModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "FrictionModel::New(const dictionary&) : " << endl
            << "    unknown FrictionModelType type "
            << FrictionModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid FrictionModelType types are :"  << endl
            <<    dictionaryConstructorTablePtr_->sortedToc() << endl;
        Info << abort(FatalError) << endl;
    }

    return autoPtr<FrictionModel>(cstrIter()(dict));
}


// ************************************************************************* //
