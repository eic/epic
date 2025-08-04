Tracking performance benchmark tests
-----------------------------------------

This directory contains scripts to produce plots of the track momentum resolution and other parameters. This is based on scripts written by Shyam Kumar:

https://github.com/eic/detector_benchmarks/tree/master/benchmarks/tracking_performances

and adapted by Tuna Tasali:

https://github.com/Tunat66/epicUK/blob/UK_Contribution

This version is set up to run on the Oxford Particle Physics Linux system and will need to be adapted for non-Oxford users. Oxford users can follow these steps:

Edit lines 15, 16, 19, 21 of ```run_single_momentum.sh``` to point to the startup script, compact file, and material map that you want to use.

Edit ```epicUK.sh``` to point to your eic-shell

Check you are not in the eic-shell, then enter:  

```./run_multi_momentum.sh ```

This will create the ```results_test``` directory and submit ``` run_single_momentum.sh``` to the cluster for all 14 momentum points, which runs ```ddsim``` and ```eicrecon```. The ```.edm4hep.root``` and ```edm4eic.root``` files are saved in the ```results_test``` directory.

Once this is done, enter:

```./analyse.sh results_test```

and it will run ```Tracking_Performances.C``` on all files, and then the ```doCompare-``` macros to produce the plots.

The momentum resolution plots are in ```results_test_control/Final_Results/pi-/mom/```


