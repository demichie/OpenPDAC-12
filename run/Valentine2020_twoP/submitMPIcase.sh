foamCleanCase

cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution

blockMesh 
checkMesh -allTopology -allGeometry

touch case.foam

changeDictionary

rm -rf 0
cp -r org.0 0
mv 0/alpha.air.init 0/alpha.air
mv 0/alpha.particles1.init 0/alpha.particles1
mv 0/alpha.particles2.init 0/alpha.particles2
mv 0/T.air.init 0/T.air
mv 0/T.particles1.init 0/T.particles1
mv 0/T.particles2.init 0/T.particles2
mv 0/U.air.init 0/U.air
mv 0/U.particles1.init 0/U.particles1
mv 0/U.particles2.init 0/U.particles2

#FOR SCALAR RUN:
foamRun


mv 0/alpha.air.run 0/alpha.air
mv 0/alpha.particles1.run 0/alpha.particles1
mv 0/alpha.particles2.run 0/alpha.particles2
mv 0/T.air.run 0/T.air
mv 0/T.particles1.run 0/T.particles1
mv 0/T.particles2.run 0/T.particles2
mv 0/U.air.run 0/U.air
mv 0/U.particles1.run 0/U.particles1
mv 0/U.particles2.run 0/U.particles2

cp ./system/controlDict.run system/controlDict
cp ./system/fvSolution.run system/fvSolution

#FOR PARALLEL RUN:
#sbatch MPIJob_run.script
#squeue

#FOR SCALAR RUN:
foamRun

