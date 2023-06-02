function dydt = MF_Integrate_order_parameters(t,y,Input,param,ms,J)

    ExtInputs= structfun(@(x) x(t),Input,'Uni',false);

    % Input: : C0, CA, CB, epsilonA,
    % We fix the value of the external angle to 0   
    ext=0; %could be any angle as long as ext=PSIA
    PSIA=0;
    PSIB=0;
    % The parameters PSI_A, PSI_B can also be calculated from equation 15 in the paper.
    % Here, we fixed them. We will later check that from simulations 
    % (by using the external input inferred through this code)
    % recover the correct value of PSI_A and PSI_B

    interp_filed=J.j0.*y(1)+ExtInputs.C0+permute(ms.edges_eta,[1,3,2]).*ExtInputs.CA...
       +permute(ms.edges_eta,[1,3,4,2]).*ExtInputs.CB...
       +permute(ms.edges_eta,[1,3,4,2]).*J.ja.*y(2).*cos(ms.edges_ang-PSIB)...
       +permute(ms.edges_eta,[1,3,4,2]).*J.jse.*y(3).*cos(ms.edges_ang-PSIB)...          
       +permute(ms.edges_eta,[1,3,2])*J.jsp.*y(2).*cos(ms.edges_ang'-PSIA)...    
       +permute(ms.edges_eta,[1,3,4,2])*ExtInputs.epsB.*cos(ms.edges_ang-ext)...
       +permute(ms.edges_eta,[1,3,2])*ExtInputs.epsA.*cos(ms.edges_ang'-ext);
   
    interp_filed(interp_filed<0)=0;  % ReLu non-linearity

    sum0=sum(ms.dx.*ms.rho_4d.*ms.wA.*ms.wA'.*ms.wE.*permute(ms.wE,[1,2,4,3]).*interp_filed,"all");
    sum1=sum(ms.dx.*ms.rho_4d.*ms.wA.*ms.wA'.*ms.wE.*permute(ms.wE,[1,2,4,3]).*permute(ms.edges_eta,[1,3,2]).*cos(ms.edges_ang'-PSIA).*interp_filed,"all");
    sum2=sum(ms.dx.*ms.rho_4d.*ms.wA.*ms.wA'.*ms.wE.*permute(ms.wE,[1,2,4,3]).*permute(ms.edges_eta,[1,3,4,2]).*cos(ms.edges_ang-PSIB).*interp_filed,"all");
    sum3=sum(ms.dx.*ms.rho_4d.*ms.wA.*ms.wA'.*ms.wE.*permute(ms.wE,[1,2,4,3]).*interp_filed.*permute(ms.edges_eta,[1,3,2]),"all");
    sum4=sum(ms.dx.*ms.rho_4d.*ms.wA.*ms.wA'.*ms.wE.*permute(ms.wE,[1,2,4,3]).*interp_filed.*permute(ms.edges_eta,[1,3,4,2]),"all");
    
    
    r0i=(-y(1)+sum0)/param.tau;
    rAi=(-y(2)+sum1)/param.tau;
    rBi=(-y(3)+sum2)/param.tau;
    r0A=(-y(4)+sum3)/param.tau;
    r0B=(-y(5)+sum4)/param.tau;
    dydt = [r0i; rAi; rBi; r0A; r0B];

end
