# module load openFOAM-10
# source $FOAM_BASHRC

foamCleanCase

cd preprocessing
python createSTL.py
cd ..

cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution

blockMesh 
checkMesh -allTopology -allGeometry

snappyHexMesh -overwrite
extrudeMesh
changeDictionary

# sbatch MPIJob_snappy.script
# extrudeMesh
# changeDictionary


rm -rf 0
cp -r org.0 0
mv 0/alpha.air.init 0/alpha.air
mv 0/alpha.particles.init 0/alpha.particles
mv 0/T.air.init 0/T.air
mv 0/T.particles.init 0/T.particles
mv 0/U.air.init 0/U.air
mv 0/U.particles.init 0/U.particles



#FOR PARALLEL RUN:
#sbatch MPIJob_init.script
#squeue

#FOR SCALAR RUN:
foamRun

mv 0/alpha.air.run 0/alpha.air
mv 0/alpha.particles.run 0/alpha.particles
mv 0/T.air.run 0/T.air
mv 0/T.particles.run 0/T.particles
mv 0/U.air.run 0/U.air
mv 0/U.particles.run 0/U.particles

cp ./system/controlDict.run system/controlDict
cp ./system/fvSolution.run system/fvSolution

#FOR PARALLEL RUN:
#sbatch MPIJob_run.script
#squeue

#FOR SCALAR RUN:
foamRun

