if (Pstream::master())
{
    gradPOSCFile
        << runTime.timeName() << "\t" << gradPOSC.value().x()
        <<"\t"<< gradPOSC.value().y() << "\t" << gradPOSC.value().z() << endl;
}
