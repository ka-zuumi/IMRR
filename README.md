# Interpolating Moving Ridge Regression

The Interpolating Moving Ridge Regression (IMRR) predicts energy gradients of molecular geometries given some number of configurations with known energy gradients. With the predictive risk, at any point during an *ab initio* molecular dynamics (AIMD) trajectory, the *ab initio* energy gradient can be substituted by a cheap IMRR predicted energy gradient provided the risk is low enough. An example AIMD trajectory for the bimolecular reaction between HBr<sup>+</sup> and CO<sub>2</sub> is shown below.

![Alt text](hbr+co2traj1.png?raw=true "Example Trajectory")

AIMD studies are rigorous ways of exploring the configuration space of a system for machine learning methods, and here for IMRR.
A number of frames from previous AIMD trajectories are consolidated and indexed in a "library".
Here, they are indexed by two collective variables (CVs), the distance between the hydrogen and carbon, and the minimum distance between the oxygens and the bromine. The density of frames can be visualized over these two collective variables with a heatmap.
A portion of the total number of frames are displayed below (those with CVs in the range (6,6) to (7,7)) as the light-green square.
An example AIMD trajectory is shown as the black curve overlayed on the heatmap below.

![Alt text](heatmapTrajectoryTrace.png?raw=true "Example Trajectory")

By making use of this library, a number *Ninterpolation* of frames can be used as inputs for the IMRR to predict energy gradients along the AIMD trajectory. The manner and number of frames chosen are specified in the ANALYSIS.f90 file. By default, for each step along the trajectory, three previous consecutive frames are used as inputs for the IMRR. All others are taken from the library. Those steps along the trajectory with inputs taken from the library are highlighted in dark green. Two examples are provided; to compile them and run the first example, execute:

```
./build_IMRR.sh
./bin/example1.o
```

All output from the analysis on this system is placed in a folder in the library. In this case, the analysis is placed in "HBrCO2library/expcompareGradientsInGridtoNWChemGradients0\_history3\_20/". The name of the output folder is specified in the ANALYSIS.f90 file. Errors between the predicted and true forces, as well as other variables, are recorded in the file "data/interpolation.dat" within the output folder. In this example, as this is a trajectory, errors and other variables are printed at each step. After running the first example, ignoring the data in the grid, the directory structure looks like this:

```
├── bin
├── build_IMRR.sh
├── HBrCO2library
│   └── expcompareGradientsInGridtoNWChemGradients0_history3_20
│       └── data
│           └── interpolation.dat
│   └── 001
│       └── grid
├── heatmapTrajectoryTrace.gnu
├── heatmapTrajectoryTrace.png
├── README.md
├── readtrajectories.txt
└── VENUSwNWChem-1.out
```

To visualize the data and make a figure like above, you can use the pre-made gnuplotfile supplied after changing the inputfile:
```
gnuplot heatmapTrajectoryTrace.gnu
```
