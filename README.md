# LDF-AdultCtrls-PrePostExe-Smooth
A collection of scripts to batch-analyze laser Doppler flowmetry (LDF) data with an experimentally-induced post-occlusive reactive hyperemia. These scripts were originally intended for internal use with the Periflux 5000 Laser Doppler Flowmeter system. There is one main method, LDF_batch. Required Matlab functions are located in the Subroutines folder. 

Time, pressure, temperature, and relative perfusion is analyzed. The output is stratified by testing session and type of session (e.g. pre-exercise, post-exercise). Outputs are generated in tabular (Excel) and graphical formats (PNG).

### Example Outputs

![](https://raw.githubusercontent.com/btran29/LDF-PORH-Analysis-Tools/master/example/fig1.png)

**Summary figure**: Summary figure for a single testing session. The top graph shows the testing session data as a whole, as well as the selected portions of the data to be analyzed. The bottom graph shows the hyperemia in detail, with the variables we used to model the decay in perfusion levels following the occlusion.

### Usage information (more information in script comments)
To run, add the batch script and its subroutines into your Matlab path. Manually pointing occlusions in a batch file is required for this particular workflow (more info in 'LDF_batch.m'). Study identifier variables are handled via 'csv_to_mat.m.' Be sure to take a look at that file if you're going to change up the naming scheme.

Warning: These scripts were my first foray into programming, and thus may be poorly written.
