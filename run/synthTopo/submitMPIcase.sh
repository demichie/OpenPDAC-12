# clear previous stl surfaces and output of previous runs
rm -rf constant/triSurface
rm -rf constant/extendedFeatureEdgeMesh
rm system/probes
rm system/cuttingPlane
foamCleanCase

touch case.foam

# create the STL surface from the ASC file
cd preprocessing
unzip synthDEM.zip
python3 ASCtoSTL.py
cd ..

# check that the surfaces for initial conditions are closed
surfaceCheck constant/triSurface/surface_crater_closed.stl
surfaceCheck constant/triSurface/surface_conduit_closed.stl
surfaceCheck constant/triSurface/surface_total_closed.stl

# create the initial uniform grid 
cp ./system/controlDict.init ./system/controlDict
blockMesh 
checkMesh -allTopology -allGeometry

# refinement of the grid
surfaceFeatures
snappyHexMesh -overwrite
checkMesh -allTopology -allGeometry

# create VTKs for cells/faces with problems
foamToVTK -cellSet underdeterminedCells -fields '(none)'
foamToVTK -cellSet concaveCells -fields '(none)'
foamToVTK -faceSet concaveFaces -fields '(none)'
foamToVTK -faceSet skewFaces -fields '(none)'

# assign common names to boundary faces 
changeDictionary

# create zones based on stl files of closed surfaces (crater, conduit,...)
topoSet -dict topoSetDict-conduit

# set fvSolution and cloudProperties for initialization of atm. profile and ballistics
cp ./system/controlDict.init ./system/controlDict
cp ./system/fvSolution.init ./system/fvSolution
cp ./constant/cloudProperties.init ./constant/cloudProperties

# remove previous 0 folder and copy not-initialized fields
rm -rf 0
cp -r org.0 0

# (first run for initialization of the solution: fileds and ballistics) 


#FOR PARALLEL RUN:
#sbatch MPIJob_init.script
#squeue

#FOR SCALAR RUN:
foamRun

# set different values in the crater/conduit zone
setFields

# set the run parameters for the actual simulation
cp ./system/controlDict.run system/controlDict
cp ./system/fvSolution.run system/fvSolution
cp ./constant/cloudProperties.run ./constant/cloudProperties

#FOR PARALLEL RUN:
#sbatch MPIJob_run.script
#squeue

#FOR PARALLEL RUN ON PC:
decomposePar
mpirun -np xx foamRun -parallel

#FOR SCALAR RUN:
foamRun

# monitoring of probe value
foamMonitor -ascii postProcessing/probes/0/p

# create symbolic links for VTKs with sorted names for Paraview
python sortVTK.py

# reconstruct alpha.particles and convert the values on terrain to VTK
reconstructPar -newTimes -fields '(alpha.particles)' -lagrangianFields '(none)'
foamToVTK -fields '(alpha.particles)' -noInternal -noFaceZones -excludePatches '(atm top)'
python plotFieldsMap.py

# remove VTK 
rm -rf VTK

# reconstruct and convert ballistic output to VTK
reconstructPar -withZero -fields '(none)' -lagrangianFields '(origId U d rho)'
foamToVTK -fields '()' -noInternal -noFaceZones -excludePatches '(atm top terrain_in terrain_out)'
python plotBallistics.py
