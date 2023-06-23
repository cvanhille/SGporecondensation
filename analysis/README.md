## Important notes

The codes analyse_frames_seeds.py and analyse_distances.py will analyse all frames of all different replicas for a particular parameter set (LLPS interaction, membrane wetting, initial pore size and protein concentration. Note that the code is designed (in its current form) to run in the specific workflow used for this set of simulations, so if you want to analyse specific files you will need to change the paths in the main bit of the code (bpath, etc., lines 291 / 384 onwards for each file respectively).

The codes compile_FrData.py and compile_distances.py recompile the output of the previous ones in more manageable data sets for each parameter set. Again, the paths are set to fit the local workflow of this simulation set, as done to generate the data for the article. You will need to tweak these to run it on your local machine according to how you set up the simulations.

Same applies to the jupyter notebooks, some of the paths from which files are read may not be valid where you clone this repository and will need adjusting.
