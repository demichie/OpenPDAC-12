# clear folder
foamCleanCase
rm -rf 0/polymesh
rm 0/U*

touch case.foam

blockMesh 
topoGrid

cp -r 0/polyMesh constant/

