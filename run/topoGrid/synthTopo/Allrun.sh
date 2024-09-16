# clear folder
foamCleanCase
rm -rf 0/polymesh

touch case.foam

blockMesh 

# SERIAL
# topoGrid

# PARALLEL
decomposePar
mpirun -np xx topoGrid -parallel
reconstructPar -withZero


cp -r 0/polyMesh constant/

