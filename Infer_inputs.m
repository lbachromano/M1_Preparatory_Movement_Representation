function [ExtInputs,costfun] = Infer_inputs(OrdParam,param,measure,J,color_map,idplot)      

    % We want to infer the value of the external inputs such that the 
    % network generates the correct dynamics of the order parameter
    First = @(V) V(1);
    % We want to reproduce the dynamics of the 5 order parameters computed
    % from data
    nr_op=length(fieldnames(OrdParam)); 
    nr_inputs=nr_op; %number of external input parameters to infer at each time step
    %it's the same as the number of order parameters

    %For the sake of computational speed, we infer the strength of the external
    %fields every "Delta_inference*param.dt_simulaz" steps. 
    %By default, Delta_inference=10, param.dt_simulaz=0.5, so we infer the 
    %value of the external inputs every 5ms. It can be changed to smaller
    %or larger values
    
    Delta_inference=10; %It can be set to 1 to infer the fields at every time step  
    param.dt_inference=Delta_inference*param.dt_simulaz; %rescale the size of the time step accordingly
    tau=param.tau/Delta_inference; %rescale the integration time constant accordingly
     
    time=(1:param.dt_inference:param.t_end); 
    time_ms=1:1:param.t_end;
    param.t_end_inference=floor(param.t_end/Delta_inference);
    
    
    %interpolate the order parameter for inference:
    OrdParam.rB(OrdParam.rB<0)=0; %does not change the results significantly 
    target_OP=[interp1(time_ms,OrdParam.r_zero',time); interp1(time_ms,OrdParam.rA',time); ...
        interp1(time_ms,OrdParam.rB',time); interp1(time_ms,OrdParam.r0A',time); interp1(time_ms,OrdParam.r0B',time)]'; 
    
    %average order parameters:
    average_ordparam= structfun(@(x) mean(x),OrdParam,'Uni',false); 
    %initialize variables:
    Inferred_fields=zeros(size(time,2),nr_inputs);
    Initial_OP=zeros(size(time,2),nr_op);
    costfun=zeros(size(time,2),1);
    lb=[-Inf; -Inf; -Inf; 0; 0;]'; %parameters \epsilon_A,\epsilon_B are constrained to be positive
    ub=[Inf; Inf; Inf; Inf; Inf;]';
    opts = optimoptions(@fmincon,'Display','off');

    %% Infer external fields at every time step
        time_ind=1;
        obj = @(input) First(To_minimize_t(input,target_OP(time_ind,:),average_ordparam,param, Initial_OP(time_ind,:),measure,J,tau));
       
        [x,fval] = fmincon(obj,Inferred_fields(time_ind,:),[],[],[],[],lb,ub);  
        Inferred_fields(time_ind,:)=x;
        Error(time_ind)=fval;   
        f=To_minimize_t(x,target_OP(time_ind,:),average_ordparam,param,Initial_OP(time_ind,:),measure,J,tau);
        Initial_OP(time_ind,:)=f(3:end)'; %updated order parameters
        costfun(time_ind,1)=f(2);
        
        %%
    for time_ind=2:size(time,2)
        % returns error + regularization for minimization of reconstruction
        % error:
        obj = @(input) First(To_minimize_t(input,target_OP(time_ind,:),average_ordparam,param, Initial_OP(time_ind-1,:),measure,J,tau));
        % perform minimization:
        [x,fval] = fmincon(obj,Inferred_fields(time_ind-1,:),[],[],[],[],lb,ub);     
        Inferred_fields(time_ind,:)=x;       
        Error(time_ind)=fval;   
        f=To_minimize_t(x,target_OP(time_ind,:),average_ordparam,param,Initial_OP(time_ind-1,:),measure,J,tau);
        Initial_OP(time_ind,:)=f(3:end)'; %updated order parameters
        costfun(time_ind,1)=f(2); %sav4 the value of the cost function 
    end

    %% Filter and plot the results
    res_time=0.5.*(param.t_end_inference.*Delta_inference-param.t_end);
    time_simul=1:param.dt_simulaz:param.t_end-res_time;    
    % optional
    % smooth the inferred fields
    kf=3; 
    C_0=hampel(Inferred_fields(:,1),kf);   
    C_A=hampel(Inferred_fields(:,2),kf);   
    C_B=hampel(Inferred_fields(:,3),kf);   
    epsilon_A=hampel(Inferred_fields(:,4),kf);   
    epsilon_B=hampel(Inferred_fields(:,5),kf); 
    %interpolate the inferred fields to the original time 
    %filter the results
    sigma_filter=20;
    ExtInputs.C0=gaussfilt(time_simul,interp1(time,C_0, time_simul),sigma_filter);
    ExtInputs.CA=gaussfilt(time_simul,interp1(time,C_A, time_simul),sigma_filter);
    ExtInputs.CB=gaussfilt(time_simul,interp1(time,C_B, time_simul),sigma_filter);
    ExtInputs.epsA=gaussfilt(time_simul,interp1(time,epsilon_A, time_simul),sigma_filter);
    ExtInputs.epsB=gaussfilt(time_simul,interp1(time,epsilon_B, time_simul),sigma_filter);

    
    %%
    if idplot ==1 
        figure
        plot(time_simul,ExtInputs.C0,'Color',color_map(8, :),'LineWidth',4)
        hold on
        plot(time_simul,ExtInputs.CA,'Color',color_map(4, :),'LineWidth',4)
        hold on
        plot(time_simul,ExtInputs.CB,'Color',color_map(5, :),'LineWidth',4)
        hold on 
        plot(time_simul,ExtInputs.epsA,'Color',color_map(2, :),'LineWidth',4)
        hold on
        plot(time_simul,ExtInputs.epsB,'Color',color_map(3, :),'LineWidth',4)
        xlabel('time (ms)')
        xlim([0 1860])
        legend('C_0','C_A','C_B','\epsilon_A','\epsilon_B')
        set(gca,'fontsize',16)
    end

end

