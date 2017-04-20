// if alpha < alphaSmall then set Ua = Ub
if (ClipUa.value()==1)
{
    forAll(alpha,cell)
    {	
        if (alpha[cell]<alphaSmall.value())
        {
            Ua[cell] = Ub[cell];
        }
    }
}

if (InitFreeze.value()==1)
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
    
    forAll(alpha, celll)
    {
// initially, we need to fix the bed level, and let the flow develop
// only cell centers are considered
// 0-->x, 1-->y, 2-->z
        if (mesh.cellCentres()[celll][ind] < fixbeddepth.value())
        {
// exclude the horizontal conponents 
            Ua[celll] = Uatemp[celll];
        }
    }
    Ua.correctBoundaryConditions();
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

   forAll(alpha, celll)
   {
// initially, we need to fix the bed level, and let the flow develop
// only cell centers are considered
// 0-->x, 1-->y, 2-->z
       if (mesh.cellCentres()[celll][ind] < fixbeddepth.value())
       {
// only exclude the streamwise conponents
           Ua[celll] = Uatemp[celll]+Uatemp2[celll];
       }
    }
    Ua.correctBoundaryConditions();
}

