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

//Production of TKE by drag force
ESD6 = (alphaVeg/rhob*Kveg*Ub) & Ub;
//Direct dissipation due to the drag force
//ESD7 = KE7* (alphaVeg/rhob*Kveg);
//sn = min(dc*(sqrt(constant::mathematical::pi/(4*(alphaVeg+alphaSmall)))-1), dp*(pow(alphaMax/(alpha+alphaSmall), 1.0/3.0)-1));
sn = dc*(sqrt(constant::mathematical::pi/(4*(1-beta+alphaSmall)))-1);
ESD7 = 1.0 / dc * pos(sn-1.79*dc) + 1.0 /sn * neg(sn-1.79*dc);
//ESD7 *= alphaVeg/(alphaVeg+alphaSmall) //Dissipation should be zero if alphaVeg=0
// Spectral shortcut
ESD8 = (alphaVeg/rhob*Kveg);
