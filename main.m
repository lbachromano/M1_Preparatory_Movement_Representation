
%Load data and set parameters
[Rates,param,J,p_tuning,color_map] = Load_data_set_param('Firing_rates.mat');

%Compute the variables \eta_A, \eta_B, \theta_A, \theta_B
Variable=Compute_tuning(Rates,param,p_tuning,1);  %To not generate the plots, set the last argument to 0

%Compute the order parameters r0,r0A,r0B,rA,rB
[OrdParam,Err_OrdParam]=Compute_order_param(Rates,param,Variable,color_map,1);

%%
%Compute the 4-dim distribution of the variabels
[moments,measure,distrib] = Distribution_all_variables(Variable,color_map,1);

%Plot the bifurcation diagram separating the homogeneous phase from the
%bump phase
Plot_bifurcation_diagram;

%% Infer the external fields
% Infer the temporal evolution of the external input parameters
[ExtInputs,costfun] = Infer_inputs(OrdParam,param,measure,J,color_map,1);

% Check the dynamics of the order parameters with the fields we inferred
Reconstructed_ordparam_meanfield;


%% Generate variables for simulations
% Generate \eta_A, \eta_B, \theta_A, \theta_B from the distribution
% inferred from data. Create the network for simulations.
DiscreteVar = Generate_and_print_variables(distrib,1);

%% Run simulations

ext_dir=1; %this variable sets the direction of motion. Goes from 1 to 8.
Simulation_run %Run simulations and prints the results

 
%% THIS PRINTS trajectories for ALL of THE NEURONS
% GENENRATES A VERY LARGE FILE

% fullFileName =  sprintf('RIS_all_dyn_ext%d_N%d.dat',ext_dir,Network);
% writematrix(Y,fullFileName,'Delimiter','tab');

%% THIS PRINTS trajectories ONLY THE 141 NEURONS closest to data in the variable space
% Print_neurons_closest_to_data; 
