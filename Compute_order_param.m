
function [Ord_param,Err_OrdParam] = Compute_order_param(Rates_fil,param,Variable,map,plotid)

 %% Compute the order parameters
 
    t=(1:1:param.t_end)';
    ang=pi.*(linspace(0,param.Nr_target-1,param.Nr_target)'...
        +linspace(1,param.Nr_target,param.Nr_target)')./8; %location of the target on the screen
    target(1,1,:)=ang;
    %% Order parameters
    Ord_param.r_zero=mean(mean(Rates_fil(:,1:param.t_end,:),3),1)';
    Ord_param.r0A=mean(mean(Variable.Eta_data_A.*Rates_fil,3),1)';
    Ord_param.r0B=mean(mean(Variable.Eta_data_B.*Rates_fil,3),1)';
    Ord_param.rA=mean(mean(Variable.Eta_data_A.*cos(Variable.Preferred_angle_A-target).*Rates_fil,3),1)';
    Ord_param.rB=mean(mean(Variable.Eta_data_B.*cos(Variable.Preferred_angle_B-target).*Rates_fil,3),1)';

    %% Standard error
    Err_OrdParam.st0=mean(mean(Rates_fil(:,1:param.t_end,:).*Rates_fil(:,1:param.t_end,:),3),1)';
    Err_OrdParam.Err_OrdParam.st0A= mean(mean(Variable.Eta_data_A.*Rates_fil.*Variable.Eta_data_A.*Rates_fil,3),1)';
    Err_OrdParam.Err_OrdParam.st0B= mean(mean(Variable.Eta_data_B.*Rates_fil.*Variable.Eta_data_B.*Rates_fil,3),1)';
    Err_OrdParam.stA=mean(mean(Variable.Eta_data_A.*cos(Variable.Preferred_angle_A-target).*Rates_fil.*Variable.Eta_data_A.*cos(Variable.Preferred_angle_A-target).*Rates_fil,3),1)';
    Err_OrdParam.stB=mean(mean(Variable.Eta_data_B.*cos(Variable.Preferred_angle_B-target).*Rates_fil.*Variable.Eta_data_B.*cos(Variable.Preferred_angle_B-target).*Rates_fil,3),1)';
   
    Err_OrdParam.st0=sqrt(mean(Err_OrdParam.st0,2)-Ord_param.r_zero.*Ord_param.r_zero)./sqrt(param.Nr_target*length(Variable.Eta_data_B));
    Err_OrdParam.stA=sqrt(mean(Err_OrdParam.stA,2)-Ord_param.rA.*Ord_param.rA)./sqrt(param.Nr_target*length(Variable.Eta_data_A));
    Err_OrdParam.stB=sqrt(mean(Err_OrdParam.stB,2)-Ord_param.rB.*Ord_param.rB)./sqrt(param.Nr_target*length(Variable.Eta_data_B));
    Err_OrdParam.Err_OrdParam.st0A=sqrt(mean(Err_OrdParam.Err_OrdParam.st0A,2)-Ord_param.r0A.*Ord_param.r0A)./sqrt(param.Nr_target*length(Variable.Eta_data_A));
    Err_OrdParam.Err_OrdParam.st0B=sqrt(mean(Err_OrdParam.Err_OrdParam.st0B,2)-Ord_param.r0B.*Ord_param.r0B)./sqrt(param.Nr_target*length(Variable.Eta_data_B));

    

    %% Smooth
    sigma_filter=50;
    Ord_param.r_zero=gaussfilt(t,Ord_param.r_zero,sigma_filter);
    Ord_param.rA=gaussfilt(t,Ord_param.rA,sigma_filter);
    Ord_param.rB=gaussfilt(t,Ord_param.rB,sigma_filter);
    Ord_param.r0A=gaussfilt(t,Ord_param.r0A,sigma_filter);
    Ord_param.r0B=gaussfilt(t,Ord_param.r0B,sigma_filter);

    
    %%
 
    %% Plot 

    if plotid==1
 
        figure
        clf
        subplot(3,1,1)
        y = Ord_param.r_zero;
        dy = Err_OrdParam.st0;   
        fill([t;flipud(t)],[y-dy;flipud(y+dy)],[.95 .95 1.0],'linestyle','none');
        line(t,y)
        hold on
        hh=plot(Ord_param.r_zero,'k','LineWidth',3);
        legend(hh,'r_0','Location','SouthEast');
        %legend('botoff')
        set(gca,'FontSize',15)
        xlim([0 param.t_end])



        subplot(3,1,2)
        y = Ord_param.r0A;
        dy = Err_OrdParam.Err_OrdParam.st0A;   
        fill([t;flipud(t)],[y-dy;flipud(y+dy)],[1.0 .9 0.9],'linestyle','none');
        line(t,y)
        hold on
        y = Ord_param.r0B;
        dy = Err_OrdParam.Err_OrdParam.st0B;  
        fill([t;flipud(t)],[y-dy;flipud(y+dy)],[0.8 0.9 0.8],'linestyle','none');
        line(t,y)
        hold on
        h(1)=plot(Ord_param.r0A,'Color',colormap(map(4, :)),'LineWidth',3);
        hold on
        h(2)=plot(Ord_param.r0B,'Color',colormap(map(5, :)),'LineWidth',3);
        legend([h(1) h(2)],'r_{0A}','r_{0B}','Location','SouthEast'); 
        set(gca,'FontSize',15)
        xlim([0 param.t_end])

        subplot(3,1,3)
        y = Ord_param.rA;
        dy = Err_OrdParam.stA;   
        fill([t;flipud(t)],[y-dy;flipud(y+dy)],[.9 .9 1.0],'linestyle','none');
        line(t,y)
        hold on
        y = Ord_param.rB;
        dy = Err_OrdParam.stB;   
        fill([t;flipud(t)],[y-dy;flipud(y+dy)],[1.0 .9 1.0],'linestyle','none');
        line(t,y)
        hold on
        gg(1)=plot(Ord_param.rA,'Color',colormap(map(2, :)),'LineWidth',3);
        hold on
        gg(2)=plot(Ord_param.rB,'Color',colormap(map(3, :)),'LineWidth',3);
        xlabel('time (ms)')
        legend([gg(1) gg(2)],'r_A','r_B','Location','SouthEast');
        xlim([0 param.t_end])
        set(gca,'FontSize',15)
   

    end

end
 
 