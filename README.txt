
#############################################################
#THIS README FILE DESCRIBES THE FILES IN THE BARS repository
#############################################################

############
plot_simulation_funcs.R
############

Plots the true survival and hazard curves for the 4 scenarios under which data are generated in the simulation studies.

############
sim_data.R
############

Simulates data according to scenarios 1-4 and is used at the start of the files standard_analysis_scenario1.R etc.

############
standard_analysis_scenario1.R (and scenarios 2,3,4)
############

Simulates data using sim_data.R.
Fits flexible parametric survival models with 0,1,2,...,10 internal knots using flexsurvspline.
Saves parameter estimates, covariance matrices, AICs, knot positions from all models. 
1000 simulated data sets are used.

############
compile_results.R
############

Summarises simulation results in a number of ways:
- calculates absolute and squared difference between estimated and true survival/hazard curves under flex surv models with 0-1 internal knots.
- identifies the best fitting flex surv model in each simulated data set.
- plots estimated survival/hazard curves from 1000 simulations (from the best fitting model) for comparison with the true curves.
- calculates coverage based on the best fitting model in each simulated data set. 

############
optim_0knot.R, optim_1knot.R, ..., option_10knot.R
############

Fits flex surv models 'by hand', i.e. by defining the likelihood and optimizing using 'optim'.
The variance matrix is obtained using fdhess.

############
standard_analysis_scenario1_OPTIM.R (and scenarios 2,3,4)
############

Simulates data using sim_data.R.
Fits flexible parametric survival models with 0,1,2,...,10 internal knots using both flexsurvspline and optim (to do this it refers to the files optim_0knot.R etc).
Saves parameter estimates, covariance matrices, AICs, knot positions from all models. 
1000 simulated data sets are used.
