
function f = To_minimize_t(input,ord_param,average_ordparam,param,InitialCondition,measure,J,tau)


    % ord_param is the target order parameters that we are trying to reconstruct, from data
    % InitialCondition is the value of the order parameters computed at the previous
    % time step
    
    %use the original tau
    temp_y=Integrate_order_parameters(InitialCondition,input,measure,J,tau); % Solve OD
    updated_ord_param=InitialCondition'+param.dt_simulaz.*temp_y;
   
 
    %Error in reconstructing the dynamics of the order parameters:
    err=(updated_ord_param(1)-ord_param(1)).*(updated_ord_param(1)-ord_param(1))./average_ordparam.r_zero...
    +(updated_ord_param(2)-ord_param(2)).*(updated_ord_param(2)-ord_param(2))./average_ordparam.rA ...
    +(updated_ord_param(3)-ord_param(3)).*(updated_ord_param(3)-ord_param(3))./average_ordparam.rB...
    +(updated_ord_param(4)-ord_param(4)).*(updated_ord_param(4)-ord_param(4))./average_ordparam.r0A...
    +(updated_ord_param(5)-ord_param(5)).*(updated_ord_param(5)-ord_param(5))./average_ordparam.r0B;

    f=zeros(7,1);
    % We minimize the reconstruction error and add a L-1 regularization term
    % preventing the fields from attaining unrealistic high positive and negative values 
    % the regularization parameter we use for this data set is 
    reg_parameter=5e-4;
    f(1)=err+reg_parameter*(abs(input(1))+abs(input(2))+abs(input(3))+input(4)+input(5)); %error + L1 regularization
    f(2)=err; 
    f(3:end)=updated_ord_param;

end
