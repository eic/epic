# Material Map for ACTS
The material map needs to be updated from the default version in calibration/ when __ANY__ geometry or material thickness is changed within the tracking volume, even the change happens on a non-sensitive structure.

## Generate a new map with auto script
Steps:
* geantino scan to record material from dd4hep simulation
* map materials on to selected sets of ACTS surfaces and boundaries. Default: use all entrance and exit surfaces of tracking layers with grid size set in xml file ("layer_material"...)
* result validation and plots.

pre-requiests:
1. more than 10 GB disk space, and two hours of time.
2. set up your ```$DETECTOR_PATH``` and install epic.

to run:
```./run_material_map_validation.sh --nevents 5000 --nparticles 1000 ```
See comments for details.

This takes about two hours, and >10GB disk space.

## Use a local material map with EICrecon
```eicrecon -Pacts:MaterialMap=/your_path/material-maps.cbor```

## Update the official material map
1. You can either generate the map locally as described above, or download the artifact ```material_map``` from a PR CI.
2. Check the generated comparison plots for any outstanding issues. Then upload the cbor file and relevant plots to [gitlab](https://eicweb.phy.anl.gov/EIC/detectors/athena/-/issues/153).
3. Copy the url of your uploaded cbor file, and update the [path](https://github.com/eic/epic/blob/540a9e1e255e276548993449be09bd275cb3ef05/compact/tracking/definitions_craterlake.xml#L203) at the bottom epic/compact/tracking/definitions_craterlake.xml with a PR.



References:
* [ACTS how_to](https://acts.readthedocs.io/en/latest/examples/howto/material_mapping.html)
* [Presentation at Tracking meeting](https://indico.bnl.gov/event/22490/contributions/87822/attachments/52848/90391/Auto%20Script%20for%20ACTS%20Material%20Map%20V30.pdf)
