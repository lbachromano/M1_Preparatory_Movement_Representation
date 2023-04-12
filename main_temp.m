
[Rates,param,p_tuning,color_map] = Load_data_set_param('Firing_rates.mat'); %Load data and parameters

Variable=Compute_tuning(Rates,param,p_tuning,0);  
% Computes the variables \eta_A, \eta_B, \theta_A and \theta_B
% If the last argument is set to 1, plots are generated

[OrdParam,Err_OrdParam]=Compute_order_param(Rates,param,Variable,color_map,1);   
% Computes the order parameters
% If the last argument is set to 1, plots are generated
    
 
