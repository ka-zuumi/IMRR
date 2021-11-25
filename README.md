# IMRR

The Interpolating Moving Ridge Regression (IMRR) predicts energy gradients of molecular geometries given some number of configurations with known energy gradients.

The example below of the IMRR is for the bimolecular reaction between HBr+ and CO2.
A number of frames from previous ab initio molecular dynamics (AIMD) trajectories are consolidated and indexed in a "library".
Here, they are indexed by two collective variables (CVs), the distance between the hydrogen and carbon, and the minimum distance between the oxygens and the bromine. The density of frames can be visualized over these two collective variables with a heatmap.
A portion of the total number of frames are displayed below (those with CVs in the range (6,6) to (7,7)).
An example AIMD trajectory is shown as the black curve overlayed on the heatmap below.

By making use of this library, a number Ninterpolation of frames can be used as inputs for the IMRR to predict energy gradients along the AIMD trajectory.

![Alt text](heatmapTrajectoryTrace.png?raw=true "Example Trajectory")

