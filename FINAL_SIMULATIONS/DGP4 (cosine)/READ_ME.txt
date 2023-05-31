This is the code corresponding to simulation 2 in the thesis. 
This code should run without issues on any 64bit machine with R version 4.1.2 installed
and with internet access (in order to access the Laguerre package from github).

The results are already in the folder, so one can run lines 1 - 53, and then 136 until the 
end of the script. This way one skips the simulation itself, which can take between 
3 and 8 days to complete, depending on the machine it is run on, and uses the simulation
results already generated. If one wishes to reproduce the simulation results, one can run the entire script.
The results obtained should be identical with the results reported on the thesis, due to the 
random seed setup. 

The script produces the sigma plots, the bias plots, and the MSE plots. In addition, this script produces
plots of the median absolute deviation of the estimators, but this is not reported in the thesis. 

There should be no need to modify the code in any way, as the relevant working directories are 
coded in an adaptive way. The folder structure of the file should not be altered in order not
to introduce bugs in the folder dependencies within the code. 