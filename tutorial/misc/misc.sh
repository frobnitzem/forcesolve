#!/bin/sh

# Example usage of miscallaneous programs included in this distro.
# David M. Rogers
# Dr. Thomas Beck Lab
# University of Cincinnati
# 12/7/2008
# This work was sponsored by a DOE CSGF.

# Please run these commands manually so you know what you are doing.
# also, 
echo "This is not a conventional shell script." && exit 1

# Use spline_data to make a spline out of a rdf file you have laying around:
# This creates two files, out.espl -- a file listing posterior average B-spline parameters
# and out.zla -- a file listing the posterior avg. and dev. of the variance and the spline tension
spline_data.py -d 1 -m 2 ~/proj/energy/methane/full_npt/run.en.xvg out

# Plot the espl function at h/10 intervals from the B-spline paramater file.
plot_bspline.py out.espl out.dat

