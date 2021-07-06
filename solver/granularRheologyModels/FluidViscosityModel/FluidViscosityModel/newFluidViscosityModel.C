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

#include "FluidViscosityModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::granularRheologyModels::FluidViscosityModel>
Foam::granularRheologyModels::FluidViscosityModel::New
(
    const dictionary& dict
)
{
    word FluidViscosityModelType(dict.get<word>("FluidViscosityModel"));

    Info<< "Selecting FluidViscosityModel "
        << FluidViscosityModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(FluidViscosityModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "FluidViscosityModel::New(const dictionary&) : " << endl
            << "    unknown FluidViscosityModelType type "
            << FluidViscosityModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid FluidViscosityModelType types are :"  << endl
            <<    dictionaryConstructorTablePtr_->sortedToc() << endl;
        Info << abort(FatalError) << endl;
    }

    return autoPtr<FluidViscosityModel>(cstrIter()(dict));
}


// ************************************************************************* //
