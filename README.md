Overview
--------
The EPIC Detector at IP6 for Electron-Ion Collider experiment.

**Detector geometry viewer:**
- [Central detector](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/detector_geo_central.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
- [Full Detector (including beamline)](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/detector_geo_full.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY290;ROTZ350;trz0;trr0;ctrl;all)
- [Inner detector (without calorimetry)](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/detector_geo_inner_detector.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
- Subsystem views:
  - [Calorimetry](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/detector_geo_calorimeters.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
  - [PID](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/detector_geo_pid_only.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
    - [dRICH](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/detector_geo_drich_only.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY290;ROTZ350;trz0;trr0;ctrl;all)
    - [pfRICH](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/detector_geo_pfrich_only.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom55;ROTY49;ROTZ350;trz0;trr0;ctrl;all&)
    - [DIRC](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/detector_geo_dirc_only.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
    - [ToF](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/detector_geo_tof_only.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom55;ROTY49;ROTZ350;trz0;trr0;ctrl;all&)
  - [Tracking](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/detector_geo_tracking_only.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
    - [Vertex detector](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/vertex_only_geo.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
  - [Far-forward](https://eic.phy.anl.gov/geoviewer/index.htm?nobrowser&file=https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/geo/detector_geo_ip6.root?job=dump_geometry&item=default;1&opt=clipx;clipy;transp30;zoom40;ROTY290;ROTZ350;trz0;trr0;ctrl;all)

<a href="https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/images/view01.pdf?job=report">
<img src="https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/images/view01.png?job=report" width="400px" />
</a>

<br />
<a href="https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/images/view01_top.pdf?job=report">
<img src="https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/images/view01_top.png?job=report" width="400px" />
</a>

[Browse latest](https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/browse/images?job=report)

[Detector views](https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/doc/dawn_views.md?job=report)

[Constants listing](https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/doc/constants.out?job=report)

Getting Started
---------------

You will  likely want to use this repository along with the IP6 repository:
```bash
git clone https://github.com/eic/epic.git
git clone https://github.com/eic/ip6.git
ln -s ../ip6/ip6 epic/ip6
```

### Compilation

To configure, build, and install the geometry (to the `install` directory), use the following commands:
```bash
cmake -B build -S . -D CMAKE_INSTALL_PREFIX=install
cmake --build build
cmake --install build
```
To load the geometry, you can use the scripts in the `install` directory:
```bash
source install/setup.sh
```
or
```tcsh
source install/setup.csh
```

### Adding/changing detector geometry

Hint: **Use the CI/CD pipelines**.

To avoid dealing with setting up all the dependencies, we recommend using the continuous integration/continuous deployment (CI/CD) pipelines to make changes and assess their effects. Any feedback to help this process is appreciated.

Here is how to begin:

1. Look at existing detector constructions and reuse if possible. Note that "compact detector descriptions" -> xml files, and "detector construction" -> cpp file.
2. Modify xml file or detector construction.
3. Create a WIP (or draft) merge request or pull request and look at the CI output for debugging. Then go to back to 2 if changes are needed.
4. Remove the WIP/Draft part of the merge request if you would like to see your changes merged into the main.

See:

- [Talk at computing round table](https://indico.jlab.org/event/420/#17-automated-workflow-for-end)

### Compiling (avoid it)

First, see if the use case above is best for you. It most likely is and can save a lot of time for newcomers.
To run the simulation locally, we suggest using the singularity image.
More details can be found at the links below:

- https://dd4hep.web.cern.ch/dd4hep/page/beginners-guide/
- https://eic.phy.anl.gov/tutorials/eic_tutorial/
- https://eicweb.phy.anl.gov/containers/eic_container/


Related useful links
--------------------

- [EIC tutorial](https://eic.phy.anl.gov/tutorials/eic_tutorial)
- [DD4hep repository](https://github.com/AIDAsoft/DD4hep)
- [DD4hep user manual](https://dd4hep.web.cern.ch/dd4hep/usermanuals/DD4hepManual/DD4hepManual.pdf)
