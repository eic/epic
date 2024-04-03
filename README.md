[![CI status](https://github.com/eic/epic/actions/workflows/linux-eic-shell.yml/badge.svg)](https://github.com/eic/epic/actions/workflows/linux-eic-shell.yml)

Overview
--------

<a href="https://eic.github.io/epic/artifacts/epic_craterlake_views/view1_top.pdf"><img align="right" alt="craterlake" src="https://eic.github.io/epic/artifacts/epic_craterlake_views/view1_top.png" width="70%"/></a>

**Detector geometry:**
- [Empty viewer](https://eic.github.io/epic/geoviewer)
- Craterlake: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_craterlake.root&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all) [step](https://eic.github.io/epic//artifacts/epic_craterlake_no_bhcal.stp/epic_craterlake_no_bhcal.stp)
- Subsystems:
  - Inner detector: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_inner_detector.root&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all) [tgeo](https://eic.github.io/epic//artifacts/tgeo/epic_inner_detector.root)
  - Calorimetry: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_calorimeters.root&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all) [tgeo](https://eic.github.io/epic//artifacts/tgeo/epic_calorimeters.root)
    - Imaging: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_imaging_only.root&item=default;1&opt=clipx;clipy;transp30;zoom55;ROTY49;ROTZ350;trz0;trr0;ctrl;all) [tgeo](https://eic.github.io/epic//artifacts/tgeo/epic_imaging.root) [step](https://eic.github.io/epic//artifacts/epic_imaging_only.stp/epic_imaging_only.stp)
  - PID: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_pid_only.root&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY320;ROTZ340;trz0;trr0;ctrl;all) [tgeo](https://eic.github.io/epic//artifacts/tgeo/epic_pid_only.root)
    - dRICH: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_drich_only.root&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY290;ROTZ350;trz0;trr0;ctrl;all) [tgeo](https://eic.github.io/epic//artifacts/tgeo/epic_drich_only.root) [step](https://eic.github.io/epic//artifacts/epic_drich_only.stp/epic_drich_only.stp)
    - pfRICH: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_pfrich_only.root&item=default;1&opt=clipx;clipy;transp30;zoom55;ROTY49;ROTZ350;trz0;trr0;ctrl;all) [tgeo](https://eic.github.io/epic//artifacts/tgeo/epic_pfrich_only.root)
    - DIRC: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_dirc_only.root&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all) [tgeo](https://eic.github.io/epic//artifacts/tgeo/epic_dirc_only.root) [step](https://eic.github.io/epic//artifacts/epic_dirc_only.stp/epic_dirc_only.stp)
  - Tracking: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_tracking_only.root&item=default;1&opt=clipx;clipy;transp30;zoom75;ROTY320;ROTZ340;trz0;trr0;ctrl;all) [tgeo](https://eic.github.io/epic//artifacts/tgeo/epic_tracking_only.root) [step](https://eic.github.io/epic//artifacts/epic_craterlake_tracking_only.stp/epic_craterlake_tracking_only.stp)
    - Vertex: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_vertex_only.root&item=default;1&opt=clipx;clipy;transp30;zoom120;ROTY320;ROTZ340;trz0;trr0;ctrl;all) [tgeo](https://eic.github.io/epic//artifacts/tgeo/epic_vertex_only.root)
    - TOF: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_tof_only.root&item=default;1&opt=clipx;clipy;transp30;zoom55;ROTY49;ROTZ350;trz0;trr0;ctrl;all) [tgeo](https://eic.github.io/epic//artifacts/tgeo/epic_tof_only.root)
  - Beamline: [viewer](https://eic.github.io/epic/geoviewer?nobrowser&file=artifacts/tgeo/epic_ip6.root&item=default;1&opt=clipx;clipy;transp30;zoom40;ROTY290;ROTZ350;trz0;trr0;ctrl;all)  [tgeo](https://eic.github.io/epic//artifacts/tgeo/epic_ip6.root) [step](https://eic.github.io/epic//artifacts/epic_ip6.stp/epic_ip6.stp)

**Detector parameters:**
- Craterlake: [text](https://eic.github.io/epic/artifacts/constants/epic_craterlake_constants.out) [toml](https://eic.github.io/epic/artifacts/constants/epic_craterlake_constants.toml) [csv](https://eic.github.io/epic/artifacts/DetectorParameterTable/epic_craterlake.csv) [html](https://eic.github.io/epic/artifacts/DetectorParameterTable/epic_craterlake.html)

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
