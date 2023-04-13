%% Data

%Load data and set parameters
[Rates,param,p_tuning,color_map] = Load_data_set_param('Firing_rates.mat');

%Compute the variables \eta_A, \eta_B, \theta_A, \theta_B
Variable=Compute_tuning(Rates,param,p_tuning,1);  %To not generate the plots, set the last argument to 0

%Compute the order parameters r0,r0A,r0B,rA,rB
[OrdParam,Err_OrdParam]=Compute_order_param(Rates,param,Variable,color_map,1);

%Compute the 4-dim distribution of the variabels
[rho_4d,moments] = Distribution_all_variables(Variable,color_map,1);

%Plot the bifurcation diagram separating the homogeneous phase from the
%bump phase
Plot_bifurcation_diagram;


 
