# clear folder
foamCleanCase
rm -rf 0/polymesh

touch case.foam

blockMesh 
topoGrid

cp -r 0/polyMesh constant/

