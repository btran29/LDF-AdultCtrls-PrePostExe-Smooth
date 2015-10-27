# LDF-AdultCtrls-PrePostExe-Smooth
Laser Doppler Flowmetry batch processing script. Takes in time, pressure, temp, and perfusion data from a laser doppler flowmeter + post occlusive hyperemia test and organizes it into a neat little workbook. The batch script will hopefully save you some time. 

To run, add the batch script and its subroutines into your Matlab path. Manual finding of occlusions is necessary for this script. There is more info in 'LDF_batch.m'.

Study identifier variables are handled via 'csv_to_mat.m.' Be sure to take a look at that file if you're going to change up the naming scheme.

Also, as it was my first foray into programming, it is poorly written for your convenience. I will attempt to make it better someday.
