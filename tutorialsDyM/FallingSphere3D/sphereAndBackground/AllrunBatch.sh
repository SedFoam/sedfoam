#!/bin/bash
#OAR -n sphere_withSedim
#OAR -l /nodes=1/core=4,walltime=24:00:00
#OAR --stdout log.out
#OAR --stderr errors.err

source /etc/profile
module load openfoam/2206plus
export OMPI_MCA_plm_rsh_agent=oar-envsh

# Nombre de processus parallèles
NBCPUS=$(cat ${OAR_NODEFILE} | wc -l)

# create the mesh
#foamCleanPolyMesh
#blockMesh


# create the intial time folder
#cp -r 0_org 0
#module load openfoam/1906plus
#mapFields -sourceTime 202 ../phi580_refined64
# Decompose the case in order to run in parallel 
#funkySetFields -time 0

#decomposePar


# Lancement
mpirun -np ${NBCPUS} -machinefile  ${OAR_NODEFILE} overSedDymFoam_rbgh -parallel 
#overSedDymFoam_rbgh
