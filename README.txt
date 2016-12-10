
#############################################################
#THIS README FILE DESCRIBES THE FILES IN THE FOLDER BARS_sim
#############################################################

############
plot_simulation_funcs.R
############

Plots the true survival and hazard curves for the 4 scenarios under which data are generated in the simulation studies.

############
sim_data.R
############

Simulates data according to scenarios 1-4 and is used at the start of the files standa-analysis_scenario1.R etc.

############
standard_analysis_scenario1.R (and scenarios 2,3,4)
############

Simulates data using sim_data.R.
Fits flexible parametric survival models with 0,1,2,...,10 internal knots using flexsurvspline.
Saves parameter estimates, covariance matrices, AICs, knot positions from all models. 
1000 simulated data sets are used.

############
FOLDER sim_results
############

The results from the simulations are stored in this folder. 

############
compile_results.R
############

Summarises simulation results in a number of ways:
- calculates absolute and squared difference between estimated and true survival/hazard curves under flex surv models with 0-1 internal knots.
- identifies the best fitting flex surv model in each simulated data set.
- plots estimated survival/hazard curves from 1000 simulations (from the best fitting model) for comparison with the true curves.
- calculates coverage based on the best fitting model in each simulated data set. 

############
optim_test.R
############

Fits flex surv models 'by hand', i.e. by defining the likelihood and optimizing using 'optim'.
Curretn this only has the likelihood for a flex surv model with 3 knots (as an example).
This was found to give the same parameter estimates as using flexsurvspline, the variances can be slightly different (using fdHess).
This could be used to obatin variances when flexsurvspline complains about the Hessian not being positive definite. 

############
FOLDER plots
############

Contains plots of true survival/hazard curves and plots of simulation results.