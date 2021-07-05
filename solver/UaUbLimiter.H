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
// if alpha < alphaSmall then set Ua = Ub
if (ClipUa.value()==1)
{
    surfaceScalarField alphaf(fvc::interpolate(alpha));
    forAll(alpha, cell)
    {
        if (alpha[cell]<alphaSmall.value())
        {
            Ua[cell] = Ub[cell];
        }
    }
    forAll(alphaf, face)
    {
        if (alphaf[face]<alphaSmall.value())
        {
            phia[face] = phib[face];
        }
    }
    phi = alphaf*phia + (1-alphaf)*phib;
}

if (InitFreeze.value()==1)
{
// in this way the bed is immobile in streamwise
    volScalarField temp(Ua & gravity);
// extract vertical component of velocity component,we get ( 0 ,0,w)
    volVectorField Uatemp(temp*gravity/magSqr(gravity));

    scalar indv1((vector(0, 1, 0) & gravity.value()));
    scalar indv2((vector(0, 0, 1) & gravity.value()));

    int ind = 0;
    if (indv1 != 0) {ind = 1;}
    if (indv2 != 0) {ind = 2;}

    forAll(alpha, cell)
    {
// initially, we need to fix the bed level, and let the flow develop
// only cell centers are considered
// 0-->x, 1-->y, 2-->z
        if (mesh.cellCentres()[cell][ind] < fixbeddepth.value())
        {
// exclude the horizontal conponents
            Ua[cell] = Uatemp[cell];
        }
    }
    phia = (fvc::interpolate(Ua) & mesh.Sf());
    phi = fvc::interpolate(alpha)*phia + fvc::interpolate(beta)*phib;
}

if (InitFreeze.value()==2)
{
// in this way the bed is immobile in streamwise
    volScalarField temp(Ua & gravity);
// extract vertical component of velocity component,we get ( 0 ,0,w)
    volVectorField Uatemp(temp*gravity/magSqr(gravity));
// extract spanwise component of velocity, we get ( 0,v,0)
    volVectorField Uatemp2((Ua & vector(0, 1, 0))*vector(0, 1, 0));

    scalar indv1((vector(0, 1, 0) & gravity.value()));
    scalar indv2((vector(0, 0, 1) & gravity.value()));

    int ind = 0;
    if (indv1 != 0) {ind = 1;}
    if (indv2 != 0) {ind = 2;}

    forAll(alpha, cell)
    {
// initially, we need to fix the bed level, and let the flow develop
// only cell centers are considered
// 0-->x, 1-->y, 2-->z
       if (mesh.cellCentres()[cell][ind] < fixbeddepth.value())
       {
// only exclude the streamwise conponents
           Ua[cell] = Uatemp[cell]+Uatemp2[cell];
       }
    }
    phia = (fvc::interpolate(Ua) & mesh.Sf());
    phi = fvc::interpolate(alpha)*phia + fvc::interpolate(beta)*phib;
}
