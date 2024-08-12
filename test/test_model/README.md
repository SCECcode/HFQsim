== On Discovery ==
module load gnuplot
module load gcc/12.3.0

./run_test.sh

Should build and try to run program. However, it will error out on the head node.

$ sbatch hfqsim.job


==Verification  Science Details ==

This new version of the code is set to run for 80 years. This version will generate a synthetic catalog, and save time evolution of all model parameters roughly every 30 simulation years. To verify your results, here are a few things you can look at: (I attached a sample figure as well for illustration)

1. Plot magnitude (column 5 of test_model_eqcat.out) vs. time (column 1 of test_model_eqcat.out), same as black circles in the first panel. You should see clear aftershock sequence following each system size events (M>6, red star in the first panel). In my figure I ignore events that are within 7.5km distance from the model boundary, but don’t worry about that for now. 

2. You can also plot average stress vs. time (red line in panel 1). The name of the stress output file should look like the following, with number (in this case 11)showing the time sequence. Since c++ does not have matrix, the output a one-dimension array, which should be divided into 256 x 64 = 16,384 points per time step. The average stress should increase roughly linearly and drop when a large earthquake occurs. 
test_model11_hist_tau.txt 

3. Time step file looks like the following. Each value shows one time step. 
test_model10_timeStep.txt

4. You can also plot standard deviation of stress value at each time step (e.g. panel 2). The value will increase at each large events and gradually decrease with the aftershock sequences. 

5. You can also check temperature/activation energy evolution (panel 3,4) but I don’t think that’s necessary right now. 

## Changes Log ##
## old version of test_model.cpp
04_22_19:
- A linear activation energy with depth, zBD = 10 km, E0 = 50 kJ/mol, Emax = 130 kJ/mol, fluctuation = +/- 5 kJ/mol
- constant Arrhenius amplitude 
- run for 10 years
- no surface event, depth distribution looks reasonable except for z between 8 and 9 km
- one big event (Mw~6.0) at year = 7, however no clear aftershocks
