if (packingLimiter)
{
    //since the particle pressure near the alphaMax is extremely large
    //we need to limit the particle concentration to alphaMax-smallnumber
    // Calculating exceeding volume fractions
    volScalarField alphaEx(max(alpha - alphaMax+scalar(0.005), scalar(0)));

    // Finding neighbouring cells of the whole domain
    labelListList neighbour = mesh.cellCells();
    scalarField cellVolumes(mesh.cellVolumes());

    forAll(alphaEx, celli)
    {
        // Finding the labels of the neighbouring cells
        labelList neighbourCell = neighbour[celli];

        // Initializing neighbouring cells contribution
        scalar neighboursEx = 0.0;

        forAll(neighbourCell, cellj)
        {
            labelList neighboursNeighbour = neighbour[neighbourCell[cellj]];
            scalar neighboursNeighbourCellVolumes = 0.0;

            forAll(neighboursNeighbour, cellk)
            {
                neighboursNeighbourCellVolumes +=
                    cellVolumes[neighboursNeighbour[cellk]];
            }

            neighboursEx +=
                alphaEx[neighbourCell[cellj]]*cellVolumes[celli]
               /neighboursNeighbourCellVolumes;
        }

        alpha[celli] += neighboursEx - alphaEx[celli];
    }
	
    alpha.correctBoundaryConditions();
}
