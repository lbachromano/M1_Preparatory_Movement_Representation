function dOPdt = Integrate_order_parameters(OP,input,ms,J,tau)


    
    ext=0; %could be anOP angle as long as theOP are the same 
    PSIA=0;
    PSIB=0;

    interp_filed=J.j0.*OP(1)+input(1)+permute(ms.edges_eta,[1,3,2]).*input(2)...
       +permute(ms.edges_eta,[1,3,4,2]).*input(3)...
       +permute(ms.edges_eta,[1,3,4,2]).*J.ja.*OP(2).*cos(ms.edges_ang-PSIB)...
       +permute(ms.edges_eta,[1,3,4,2]).*J.jse.*OP(3).*cos(ms.edges_ang-PSIB)...          
       +permute(ms.edges_eta,[1,3,2])*J.jsp.*OP(2).*cos(ms.edges_ang'-PSIA)...    
       +permute(ms.edges_eta,[1,3,4,2])*input(5).*cos(ms.edges_ang-ext)...
       +permute(ms.edges_eta,[1,3,2])*input(4).*cos(ms.edges_ang'-ext);
   
    interp_filed(interp_filed<0)=0;  % ReLu non-linearitOP

    sum0=sum(ms.dx.*ms.rho_4d.*ms.wA.*ms.wA'.*ms.wE.*permute(ms.wE,[1,2,4,3]).*interp_filed,"all");
    sum1=sum(ms.dx.*ms.rho_4d.*ms.wA.*ms.wA'.*ms.wE.*permute(ms.wE,[1,2,4,3]).*permute(ms.edges_eta,[1,3,2]).*cos(ms.edges_ang'-PSIA).*interp_filed,"all");
    sum2=sum(ms.dx.*ms.rho_4d.*ms.wA.*ms.wA'.*ms.wE.*permute(ms.wE,[1,2,4,3]).*permute(ms.edges_eta,[1,3,4,2]).*cos(ms.edges_ang-PSIB).*interp_filed,"all");
    sum3=sum(ms.dx.*ms.rho_4d.*ms.wA.*ms.wA'.*ms.wE.*permute(ms.wE,[1,2,4,3]).*interp_filed.*permute(ms.edges_eta,[1,3,2]),"all");
    sum4=sum(ms.dx.*ms.rho_4d.*ms.wA.*ms.wA'.*ms.wE.*permute(ms.wE,[1,2,4,3]).*interp_filed.*permute(ms.edges_eta,[1,3,4,2]),"all");
    
    
    r0i=(-OP(1)+sum0)/tau;
    rAi=(-OP(2)+sum1)/tau;
    rBi=(-OP(3)+sum2)/tau;
    r0A=(-OP(4)+sum3)/tau;
    r0B=(-OP(5)+sum4)/tau;
    dOPdt = [r0i; rAi; rBi; r0A; r0B];

end
