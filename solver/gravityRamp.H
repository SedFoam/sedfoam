if (gRampTime != 0 && Foam::mag(runTime.value()/gRampTime)<=1.0)
{
    Foam::scalar gRampFactor(min(1, Foam::mag(runTime.value() / gRampTime)));

    Info << "gRampFactor: " << gRampFactor << endl;
    g = gRampFactor * gEnd;
    
    // Update gh
    Info << "Updating gravity field\n" << endl;
    gh = g & mesh.C();
}
