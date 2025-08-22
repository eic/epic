Curved OB geometry test
-----------------------

This branch has been created by Sam Henry to test a new outer barrel geometry with curved silicon surface.

The OB stave structure has replaced the flat silicon surface with a curved surface. The carbon fibre frame has been replaced with a flat layer. The staves are no longer tilted and are positioned alternately at radii separated by 6mm.

In this test, the curved silicon surface for each stave is modelled as a single module - a segment of a cylinder with a radius of 80mm. In the final version we will include the top and bottom layers of the stave with 4 modules for L3 and 8 for L4. 

This is not yet working - eicrecon does not see the outer barrel hits.

**Modified files:**

```src/BarrelTrackerWithFrame_geo.cpp```

```compact/tracking/silicon_barrel.xml```

**Additional files:**

```SimpleCurved_silicon_barrel.xml``` - alternative version of silicon_barrel.xml with a simpler curved model where all staves are segments of a single big cylinder.

```TestOB.sh``` - script to test the geometry, using ```TestOB.C``` and  ```epic_craterlake_tracking_only_cut.xml```


The test script runs a simulation shooting 1000 muons through the vertex and silicon barrels. It then runs eicrecon and prints the mean of CentralCKFTrajectories.nMeasurements - i.e. the number the barrels the tracks go through. If everything is well, this should be close to 5. If eicrecon is not seeing the outer barrels, it will be 3


**Tests**

Control test with epic-main geometry
```
source /opt/detector/epic-main/bin/thisepic.sh
./TestOB.sh
```
Mean nMeasurements: 5.02398

Test new curved model
```
source install/bin/thisepic.sh
./TestOB.sh

Mean nMeasurements: 3
```
Test simple curved model
```
cp SimpleCurved_silicon_barrel.xml install/share/epic/compact/tracking/silicon_barrel.xml 
./TestOB.sh
```
Mean nMeasurements: 4.947
