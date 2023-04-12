
function [Rates_fil,param,p_tuning,color_map] = Load_data_set_param(inputfile)

    ld=load(inputfile); % Load firing rates
    Rates_fil=ld.Rates_fil; % dimension = Nr_units x t_end x Nr_trials

    %% Parameters

    param.Nr_units=size(Rates_fil,1);
    param.Nr_target=size(Rates_fil,3);

    param.t_go=ld.t_go;
    param.t_start=ld.t_start;
    param.t_end=ld.t_end;

    p_tuning.window=300; 
    p_tuning.start_P=150;
    p_tuning.start_E=param.t_start-50;

    %% Colors for Plots

    color_map=[254 202 48;...
    2 40 75;...
    177 34 31;...
    232 127 46;...
    0 127 0;...  %
    28 177 175;...
    25 110 180;...
    110 47 140] / 255;

end