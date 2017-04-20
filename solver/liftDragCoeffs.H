volVectorField Ur(-Ua + Ub);
volScalarField magUr(mag(Ur));

volScalarField Ka(draga->K(magUr));

K = Ka;

if (dragPhase == "b")
{
    volScalarField Kb(dragb->K(magUr));
    K = Kb;
}
else if (dragPhase == "blended")
{
    volScalarField Kb(dragb->K(magUr));
    K = (beta*Ka + alpha*Kb);
}

volVectorField liftCoeff(Cl*(beta*rhob + alpha*rhoa)*(Ur ^ fvc::curl(U)));
