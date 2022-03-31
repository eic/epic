# subsytem_views

These compact detector files are copies of the main `athena.xml` compact description file, modified to only show certain subsystems.

To add or update a subsystem:

1. Copy `athena.xml` to this directory (`compact/subsystem_views`) with the appropriate name.
2. Remove all the undesired detectors from this xml file. **Do not modify any of the included xml files** -- only modify the top level xml file here.
3. Check that that you have the desired color scheme (eg. `colors.xml`) and display attributes (eg. `display_detailed.xml`) included.

These compact detector files should not be used for any real simulations or studies.
