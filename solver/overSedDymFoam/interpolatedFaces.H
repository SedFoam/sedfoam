// Interpolation used
interpolationCellPoint<vector> UaInterpolator(HabyA);
interpolationCellPoint<vector> UbInterpolator(HbbyA);

// Determine faces on outside of interpolated cells
bitSet isOwnerInterpolatedFace(mesh.nInternalFaces());
bitSet isNeiInterpolatedFace(mesh.nInternalFaces());

// Determine donor cells
labelListList donorCell(mesh.nInternalFaces());

scalarListList weightCellCells(mesh.nInternalFaces());

// Interpolated HbyA faces
vectorField UaIntFaces(mesh.nInternalFaces(), Zero);
vectorField UbIntFaces(mesh.nInternalFaces(), Zero);

// Determine receptor neighbour cells
labelList receptorNeigCell(mesh.nInternalFaces(), -1);

{
    const cellCellStencilObject& overlap = Stencil::New(mesh);
    const labelList& cellTypes = overlap.cellTypes();
    const labelIOList& zoneID = overlap.zoneID();

    label nZones = gMax(zoneID)+1;
    PtrList<fvMeshSubset> meshParts(nZones);
    labelList nCellsPerZone(nZones, Zero);

    // A mesh subset for each zone
    forAll(meshParts, zonei)
    {
        meshParts.set
        (
            zonei,
            // Select cells where the zoneID == zonei
            new fvMeshSubset(mesh, zonei, zoneID)
        );
    }

    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        label ownType = cellTypes[mesh.faceOwner()[faceI]];
        label neiType = cellTypes[mesh.faceNeighbour()[faceI]];
        if
        (
            ownType == cellCellStencil::INTERPOLATED
            && neiType == cellCellStencil::CALCULATED
        )
        {
            isOwnerInterpolatedFace.set(faceI);

            const vector& fc = mesh.faceCentres()[faceI];

            for (label zoneI = 0; zoneI < nZones; zoneI++)
            {
                if (zoneI != zoneID[mesh.faceOwner()[faceI]])
                {
                    const fvMesh& partMesh = meshParts[zoneI].subMesh();
                    const labelList& cellMap = meshParts[zoneI].cellMap();
                    label cellI = partMesh.findCell(fc);

                    if (cellI != -1)
                    {
                         // Determine weights
                        labelList stencil(partMesh.cellCells()[cellI]);

                        stencil.append(cellI);

                        label st = stencil.size();

                        donorCell[faceI].setSize(st);

                        weightCellCells[faceI].setSize(st);

                        scalarField weights(st);

                        forAll(stencil, i)
                        {
                            scalar d = mag
                            (
                                partMesh.cellCentres()[stencil[i]]
                              - fc
                            );
                            weights[i] = 1.0/d;
                            donorCell[faceI][i] = cellMap[stencil[i]];
                        }
                        weights /= sum(weights);

                        weightCellCells[faceI] = weights;

                        forAll(stencil, i)
                        {
                            UaIntFaces[faceI] +=
                                weightCellCells[faceI][i]
                               *UaInterpolator.interpolate
                                (
                                    fc,
                                    donorCell[faceI][i]
                                );
                            UbIntFaces[faceI] +=
                                weightCellCells[faceI][i]
                               *UbInterpolator.interpolate
                                (
                                    fc,
                                    donorCell[faceI][i]
                                );
                        }

                        break;
                    }
                }
            }

            receptorNeigCell[faceI] = mesh.faceNeighbour()[faceI];
        }
        else if
        (
            ownType == cellCellStencil::CALCULATED
            && neiType == cellCellStencil::INTERPOLATED
        )
        {
            isNeiInterpolatedFace.set(faceI);

            const vector& fc = mesh.faceCentres()[faceI];
            for (label zoneI = 0; zoneI < nZones; zoneI++)
            {
                if (zoneI != zoneID[mesh.faceNeighbour()[faceI]])
                {
                    const fvMesh& partMesh = meshParts[zoneI].subMesh();
                    const labelList& cellMap = meshParts[zoneI].cellMap();
                    label cellI = partMesh.findCell(fc);

                    if (cellI != -1)
                    {
                        // Determine weights
                        labelList stencil(partMesh.cellCells()[cellI]);

                        stencil.append(cellI);

                        label st = stencil.size();

                        donorCell[faceI].setSize(st);

                        weightCellCells[faceI].setSize(st);

                        scalarField weights(st);

                        forAll(stencil, i)
                        {
                            scalar d = mag
                            (
                                partMesh.cellCentres()[stencil[i]]
                              - fc
                            );
                            weights[i] = 1.0/d;
                            donorCell[faceI][i] = cellMap[stencil[i]];
                        }
                        weights /= sum(weights);

                        weightCellCells[faceI] = weights;

                        forAll(stencil, i)
                        {
                            UaIntFaces[faceI] +=
                                weightCellCells[faceI][i]
                               *UaInterpolator.interpolate
                                (
                                    fc,
                                    donorCell[faceI][i]
                                );
                            UbIntFaces[faceI] +=
                                weightCellCells[faceI][i]
                               *UbInterpolator.interpolate
                                (
                                    fc,
                                    donorCell[faceI][i]
                                );
                        }

                        break;
                    }
                }
            }

            receptorNeigCell[faceI] = mesh.faceOwner()[faceI];
        }
    }
}

// contravariant U
vectorField Ua1Contrav(mesh.nInternalFaces(), Zero);
vectorField Ub1Contrav(mesh.nInternalFaces(), Zero);

surfaceVectorField faceNormals(mesh.Sf()/mesh.magSf());

forAll(isNeiInterpolatedFace, faceI)
{
    label cellId = -1;
    if (isNeiInterpolatedFace.test(faceI))
    {
        cellId = mesh.faceNeighbour()[faceI];
    }
    else if (isOwnerInterpolatedFace.test(faceI))
    {
        cellId = mesh.faceOwner()[faceI];
    }

    if (cellId != -1)
    {
        const vector& n = faceNormals[faceI];
        vector n1(Zero);

        // 2-D cases
        if (mesh.nSolutionD() == 2)
        {
            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                if (mesh.geometricD()[cmpt] == -1)
                {
                    switch (cmpt)
                    {
                        case vector::X:
                        {
                            n1 = vector(0, n.z(), -n.y());
                            break;
                        }

                        case vector::Y:
                        {
                            n1 = vector(n.z(), 0, -n.x());
                            break;
                        }

                        case vector::Z:
                        {
                            n1 = vector(n.y(), -n.x(), 0);
                            break;
                        }
                    }
                }
            }
        }
        else if (mesh.nSolutionD() == 3)
        {
            //Determine which is the primary direction
            if (mag(n.x()) > mag(n.y()) && mag(n.x()) > mag(n.z()))
            {
                n1 = vector(n.y(), -n.x(), 0);
            }
            else if (mag(n.y()) > mag(n.z()))
            {
                n1 = vector(0, n.z(), -n.y());
            }
            else
            {
                n1 = vector(-n.z(), 0, n.x());
            }
        }
        n1.normalise();

        const vector n2 = normalised(n ^ n1);

        tensor rot =
            tensor
            (
               n.x() ,n.y(), n.z(),
               n1.x() ,n1.y(), n1.z(),
               n2.x() ,n2.y(), n2.z()
            );

//         tensor rot =
//             tensor
//             (
//                n  & x ,n  & y, n  & z,
//                n1 & x ,n1 & y, n1 & z,
//                n2 & x ,n2 & y, n2 & z
//             );

        Ua1Contrav[faceI].x() =
           2*transform(rot, UaIntFaces[faceI]).x()
           - transform(rot, HabyA[receptorNeigCell[faceI]]).x();

        Ua1Contrav[faceI].y() = transform(rot, HabyA[cellId]).y();

        Ua1Contrav[faceI].z() = transform(rot, HabyA[cellId]).z();

        Ub1Contrav[faceI].x() =
           2*transform(rot, UbIntFaces[faceI]).x()
           - transform(rot, HbbyA[receptorNeigCell[faceI]]).x();

        Ub1Contrav[faceI].y() = transform(rot, HbbyA[cellId]).y();

        Ub1Contrav[faceI].z() = transform(rot, HbbyA[cellId]).z();
        
        HbbyA[cellId] = transform(inv(rot), Ub1Contrav[faceI]);
    }
}
