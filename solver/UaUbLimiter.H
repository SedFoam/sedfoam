// if alpha < alphaSmall then set Ua = Ub
if (ClipUa.value()==1)
{
    surfaceScalarField alphaf = fvc::interpolate(alpha);
    forAll(alpha,cell)
    {	
        if (alpha[cell]<alphaSmall.value())
        {
            Ua[cell] = Ub[cell];
        }
    }
    forAll(alphaf,face)
    {	
        if (alphaf[face]<alphaSmall.value())
        {
            phia[face] = phib[face];
        }
    }
    phi = alphaf*phia + fvc::interpolate(beta)*phib;
}

if (InitFreeze.value()==1)
{
// in this way the bed is immobile in streamwise
    volScalarField temp = Ua & g;
// extract vertical component of velocity component,we get ( 0 ,0,w)
    volVectorField Uatemp = temp*g/magSqr(g);
    
    scalar indv1 = (vector(0,1,0) & g.value());
    scalar indv2 = (vector(0,0,1) & g.value());
    
    int ind = 0;
    if (indv1 !=0) ind = 1;
    if (indv2 !=0) ind = 2;
    
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
    volScalarField temp = Ua & g;
// extract vertical component of velocity component,we get ( 0 ,0,w)
    volVectorField Uatemp = temp*g/magSqr(g);
// extract spanwise component of velocity, we get ( 0,v,0)
    volVectorField Uatemp2 = (Ua & vector(0,1,0))*vector(0,1,0);

    scalar indv1 = (vector(0,1,0) & g.value());
    scalar indv2 = (vector(0,0,1) & g.value());

    int ind = 0;
    if (indv1 !=0) ind = 1;
    if (indv2 !=0) ind = 2;

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

