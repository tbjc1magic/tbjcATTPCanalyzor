# tbjcATTPCanalyzor

########################
##### instruction ######
########################

Though most of the programs are unfinished works, but they should give good insights of what the data look like.

The most complete program is the VertexAnalyzer, which reconstruct the

########################
##### instruction ######
########################

besides built-in Anaconda packages, you will also need opencv2 and seaborn (mainly for visualizing, when debug option is on)

conda install -c menpo opencv
conda install seaborn

########################
##### running ##########
########################

under the main folder, run for the multi-process mode

python Multi.py

the Multi.py will produce a text file contains a list of ranges for reaction length.

Or run the jupyter notebook interactive mode

jupyter notebook

once you in the jupyter browser, open either main.ipynb or MultiProcess.ipynb

#####################
##### test ##########
#####################

This problem has been tested under local workstations, google cloud Ubuntu images, and Regulus HPC node. Ideally, this analysis package should work under all environment.

########################
##### responsibility ###
########################

Well, use at your own risk. I won't be able to maintain the code.
