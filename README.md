Overview
--------
The EPIC Detector at IP6 for Electron-Ion Collider experiment.

**Detector geometry viewer:**
- [Viewer only](https://eic.github.io/epic/geoviewer)
- Main configurations:
  - [Arches](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/arches.root&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
  - [Brycecanyon](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/brycecanyon.root&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
- Additional viewers:
  - [Inner detector (without calorimetry)](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/inner_detector.root&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
  - [Calorimetry](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/calorimeters.root&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
    - [Sciglass](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/sciglass_only.root&item=default;1&opt=clipx;clipy;transp30;zoom55;ROTY49;ROTZ350;trz0;trr0;ctrl;all)
  - [PID](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/pid_only.root&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
    - [dRICH](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/drich_only.root&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY290;ROTZ350;trz0;trr0;ctrl;all)
    - [mRICH](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/mrich_only.root&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY290;ROTZ350;trz0;trr0;ctrl;all)
    - [pfRICH](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/pfrich_only.root&item=default;1&opt=clipx;clipy;transp30;zoom55;ROTY49;ROTZ350;trz0;trr0;ctrl;all)
    - [DIRC](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/dirc_only.root&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
  - [Tracking](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/tracking_only.root&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
    - [Vertex](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/vertex_only.root&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all)
    - [ToF](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/tof_only.root&item=default;1&opt=clipx;clipy;transp30;zoom55;ROTY49;ROTZ350;trz0;trr0;ctrl;all)
  - [Beamline](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/ip6.root&item=default;1&opt=clipx;clipy;transp30;zoom40;ROTY290;ROTZ350;trz0;trr0;ctrl;all)

[<img title="arches" src="https://eic.github.io/epic/artifacts/epic_arches_views/view1_top.png" width="49%">](https://eic.github.io/epic/artifacts/epic_arches_views/view1_top.pdf) [<img title="brycecanyon" src="https://eic.github.io/epic/artifacts/epic_brycecanyon_views/view1_top.png" width="49%">](https://eic.github.io/epic/artifacts/epic_brycecanyon_views/view1_top.pdf)

[Browse latest](https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/browse/images?job=report)

[Detector views](https://eicweb.phy.anl.gov/EIC/detectors/ecce/-/jobs/artifacts/main/raw/doc/dawn_views.md?job=report)

**Detector parameters:**
- text: [Arches](https://eic.github.io/epic/artifacts/constants/epic_arches.out) [Brycecanyon](https://eic.github.io/epic/artifacts/constants/epic_brycecanyon.out)
- toml: [Arches](https://eic.github.io/epic/artifacts/constants/epic_arches.toml) [Brycecanyon](https://eic.github.io/epic/artifacts/constants/epic_brycecanyon.toml)
- csv: [Arches](https://eic.github.io/epic/artifacts/DetectorParameterTable/epic_arches.csv) [Brycecanyon](https://eic.github.io/epic/artifacts/DetectorParameterTable/epic_brycecanyon.csv)


Getting Started
---------------

Get a copy of the latest version from this repository:
```bash
git clone https://github.com/eic/epic.git
```

### Compilation

To configure, build, and install the geometry (to the `install` directory), use the following commands:
```bash
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install
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
- [ACTS DD4hep plugin documentation](https://acts.readthedocs.io/en/latest/plugins/dd4hep.html)
