


================
Quick analysis
================

Once the dependencies have been installed, the quickest and easiest way to run the analysis is to simply
type::


    $ python <path_to_proj029>/proj029/scripts/run_analysis.py 
             --cgatdir=<path_to_cgat>
             --proj029dir=<path_to_proj029>
             --rdir=<path_to_R_install>
             --datadir=<path_to_data>


This will run through all the analyses in the subsequent pages (except for the computationally intensive ones). That is symbolic
links will be created to data that is required to recreate analyses and figures. It would be advisable to read through the rest of this
documentation however to learn which files relate to which steps.




