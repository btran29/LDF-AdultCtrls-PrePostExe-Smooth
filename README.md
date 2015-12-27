# LDF-AdultCtrls-PrePostExe-Smooth

A collection of scripts to batch-analyze laser Doppler flowmetry (LDF) data with an experimentally-induced post-occlusive reactive hyperemia. The main method takes time, pressure, temperature, and perfusion data and organizes it into a summary table. 

![](https://github.com/btran29/LDF-AdultCtrls-PrePostExe-Smooth/blob/master/example/fig1.png)
**Example figure**: Summary figure for a single testing session. The top graph shows the testing session data as a whole, as well as the selected portions of the data to be analyzed. The bottom graph shows the hyperemia in detail, with the variables we used to model the decay in perfusion levels following the occlusion.

![](https://github.com/btran29/LDF-AdultCtrls-PrePostExe-Smooth/blob/master/example/fig2.png)
**Locating the hyperemia**: Using a polynomial fitting tool to determine 'steady-state flow' following an occlusion-induced hyperemia. Steady-state is defined here as the first critical point of the second derivative of a 3rd degree polynomial fit to the data.

To run, add the batch script and its subroutines into your Matlab path. Manually pointing occlusions in a batch file is required for this particular workflow (more info in 'LDF_batch.m'). Study identifier variables are handled via 'csv_to_mat.m.' Be sure to take a look at that file if you're going to change up the naming scheme.

Also, as it was my first foray into programming, it is poorly written for your convenience. I have uploaded the files to github to improve upon them someday.
