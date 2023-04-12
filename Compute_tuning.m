function [Variable]=Compute_tuning(Rates_fil,param,p_tuning, plotid)


    %% Preparatory and execution activity
    Prep_matrix=Rates_fil(:,p_tuning.start_P:p_tuning.start_P+p_tuning.window,:);
    Exec_matrix=Rates_fil(:,p_tuning.start_E:p_tuning.start_E+p_tuning.window,:);

    %% Tuning properties 
    [Variable.Eta_data_A,Variable.Preferred_angle_A]=Prepare_tuning(Prep_matrix);
    [Variable.Eta_data_B,Variable.Preferred_angle_B]=Prepare_tuning(Exec_matrix);

    %% z-scored spatially localized activity
    Plot_P=(mean(Prep_matrix,2)-mean(mean(Prep_matrix,2),3))./max(mean(Prep_matrix,2),[],3);
    Plot_E=(mean(Exec_matrix,2)-mean(mean(Exec_matrix,2),3))./max(mean(Exec_matrix,2),[],3);

    [ll idx_thetaA]=sort(Variable.Preferred_angle_A);
    [ll idx_thetaB]=sort(Variable.Preferred_angle_B);



    if plotid ==1
        %% Plot
        condh=6; %choose one example to show from 1 to Nr_target=8
        idxn=randsample(param.Nr_units,1);  %one neuron shown in orange
        figure
        subplot(2,2,1)
        plot(Variable.Preferred_angle_A(idx_thetaA),Plot_P(idx_thetaA,:,condh),'k.','markersize',30)
        hold on 
        plot(Variable.Preferred_angle_A(idxn),Plot_P(idxn,:,condh),'.','markersize',40)
        xlabel('Preferred direction during preparation')
        ylabel('Preparation')
        set(gca,'fontsize',12)
        subplot(2,2,2)
        plot(Variable.Preferred_angle_B(idx_thetaB),Plot_E(idx_thetaB,:,condh),'k.','markersize',30)
        hold on 
        plot(Variable.Preferred_angle_B(idxn),Plot_E(idxn,:,condh),'.','markersize',40)
        xlabel('Preferred direction during execution')
        set(gca,'fontsize',12)
        subplot(2,2,3)
        plot(Variable.Preferred_angle_A(idx_thetaA),Plot_E(idx_thetaA,:,condh),'k.','markersize',30)
        hold on    
        plot(Variable.Preferred_angle_A(idxn),Plot_E(idxn,:,condh),'.','markersize',40)
        xlabel('Preferred direction during preparation')
        ylabel('Execution')
        set(gca,'fontsize',12)
        subplot(2,2,4)
        plot(Variable.Preferred_angle_B(idx_thetaB),Plot_P(idx_thetaB,:,condh),'k.','markersize',30)
        hold on    
        plot(Variable.Preferred_angle_B(idxn),Plot_P(idxn,:,condh),'.','markersize',40)
        xlabel('Preferred direction during execution')
        set(gca,'fontsize',12)

        %%

        Exampl_neurons=[14 9 6];
        r_plot=Prepare_tuning_plots(Prep_matrix,Exec_matrix,Exampl_neurons);

        %%
        x0=10;
        y0=10;
        width=400;
        height=400

        label_vector={'0','\pi/2' ,'\pi' ,'3 \pi/2' ,'2 \pi'};
        figure
        plot(Variable.Preferred_angle_A,Variable.Preferred_angle_B,'k.','MarkerSize',20)
        set(gca, 'XTick', 0:2*pi/4:2*pi)
        set(gca, 'YTick', 0:2*pi/4:2*pi)
        set(gca, 'XTickLabel', label_vector)
        set(gca, 'YTickLabel', label_vector)
        set(gca, 'fontsize',18)
        xlabel('\theta_A','FontSize',20)
        ylabel('\theta_B','FontSize',20)
        pbaspect([1 1 1])
        set(gcf,'position',[x0,y0,width,height])


        figure
        plot(Variable.Eta_data_A,Variable.Eta_data_B,'k.','MarkerSize',20)
        set(gca, 'fontsize',18)
        xlabel('\eta_A','FontSize',20)
        ylabel('\eta_B','FontSize',20)
        pbaspect([1 1 1])
        set(gcf,'position',[x0,y0,width,height])
    end 

end
 