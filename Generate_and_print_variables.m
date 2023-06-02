function DiscreteVar = Generate_and_print_variables(distrib,plot)

    rng(19); %NOISE SEED. THIS IS TO REPRODUCE THE HISTOGRAMS IN THE PAPER

    N_TH=25; %this creates a network with aroud 16000 units

    %% THIS SECTION USES THE INVERSE METHOD TO SAMPLE POINTS FROM THE JOINT DISTRIBUTION
    % OF THETA_A AND THETA_B, in such a way to have equidistant point along the diagonal theta_A-theta_B=0, which is
    % important for the case when couplings are above the bifurcation line 

    clear Th_Aq
    clear Th_Bq

    tha_B_min=0;
    tha_B_max=4*pi;
    tha_A_min=0;
    tha_A_max=4*pi;

    Th_B_list=linspace(tha_B_min,tha_B_max,N_TH);
    eqn=@(x)0.117047*x +0.0780314*sin(x);  
    startys=eqn(tha_A_min);
    endys=eqn(tha_A_max);
    interval_to_invert=linspace(startys,endys,N_TH);

    indice=0;
    for i=1:N_TH
        for j=1:N_TH
            indice=indice+1;     
            eqn=@(x)0.117047*x +0.0780314*sin(x)- interval_to_invert(i); 
            x=fsolve(eqn,1);  
            Th_Aq(indice)=x;                 
            Th_Bq(indice)=Th_B_list(j);   
        end
    end
    angle_rotation=-pi/4;
    New_thA=cos(angle_rotation).*Th_Aq-sin(angle_rotation).*Th_Bq;
    New_thB=sin(angle_rotation).*Th_Aq+cos(angle_rotation).*Th_Bq;
    New_thA=New_thA.*(2/sqrt(8))-pi;
    New_thB=New_thB.*(2/sqrt(8))+pi;

    thr=0.01;
    in=find(New_thA>=thr & New_thA-2*pi<thr & New_thB>=thr & New_thB-2*pi<thr);
    Th_Aq=New_thA(in);
    Th_Bq=New_thB(in);

    NTOT=length(Th_Aq);


    %% this uses matlab routines to sample from the distribution of eta_A and eta_B that 
    % was inferred earlier

    NTOT=56;
    nr_bins=500;

    clear vals
    clear dumbdist2
    pts1=linspace(0.0,1,nr_bins);
    nlst=linspace(0,2*pi,nr_bins);

    ps_eta_A=interp1(distrib.lst,distrib.distrib_eta_A,pts1);
    ps_eta_B=interp1(distrib.lst,distrib.distrib_eta_B,pts1);



    for i=1:nr_bins
      for j=1:nr_bins
        pro1=ps_eta_A(i);
        pro2=ps_eta_B(j);
        dumbdist2(j,i)=pro1*pro2;
      end
    end
    dumbdist2(isnan(dumbdist2)==1)=0;
    vals=zeros(2,NTOT);
    for i=1:NTOT
    [vals(1,i),vals(2,i)]=pinky(pts1,pts1,dumbdist2);
    end
    Eta_Aq=vals(1,:)';
    Eta_Bq=vals(2,:)';

    [h,XE,YE]=histcounts2(Eta_Aq,Eta_Bq,10);
    matrix_eta_histogram=floor(h);
    sum(sum(matrix_eta_histogram))
    Eta_center_bins=XE(1:end-1);


    %%

    N_effective=length(Th_Aq)*sum(sum(matrix_eta_histogram));
    Eta_A_effective=zeros(N_effective,1);
    Eta_B_effective=zeros(N_effective,1);
    Th_A_effective=zeros(N_effective,1);
    Th_B_effective=zeros(N_effective,1);

    contatore=1;
    
    for i=1:size(matrix_eta_histogram,1)
        for j=1:size(matrix_eta_histogram,1)
            for k=1:length(Th_Aq)


                    n_etas=matrix_eta_histogram(i,j);
                    n_thetas=1;
                    n_h=n_etas*n_thetas;
                    Eta_A_effective(contatore:contatore+n_h-1)=Eta_center_bins(i);
                    Eta_B_effective(contatore:contatore+n_h-1)=Eta_center_bins(j);
                    Th_A_effective(contatore:contatore+n_h-1)=Th_Aq(k);
                    Th_B_effective(contatore:contatore+n_h-1)=Th_Bq(k);
                    contatore=contatore+n_h;


            end
        end
    end


    %%
    DiscreteVar.Eta_Aq=Eta_A_effective;
    DiscreteVar.Eta_Bq=Eta_B_effective;
    DiscreteVar.Th_Aq=Th_A_effective;
    DiscreteVar.Th_Bq=Th_B_effective;
    DiscreteVar.NTOT=length(Th_B_effective);
    DiscreteVar.zero_etas=find(Eta_Aq==0 & Eta_Bq==0);
    DiscreteVar.N_zero=length(DiscreteVar.zero_etas);

    %% Plots
    if plot > 0
        figure
        histogram2(Eta_Aq,Eta_Bq,10);
        xlabel('\eta_A','FontSize',20)
        ylabel('\eta_B','FontSize',20)
        set(gca,'fontsize',15)

        label_vector={'0','\pi/2','\pi','3\pi/2','2 \pi'};
        figure
        scatter(Th_Aq,Th_Bq,'filled')
        xlabel('\theta_A','FontSize',20)
        ylabel('\theta_B','FontSize',20)
        set(gca, 'YTick', 0:2*pi/4:2*pi)
        set(gca, 'YTickLabel', label_vector)
        set(gca, 'XTick', 0:2*pi/4:2*pi)
        set(gca, 'XTickLabel', label_vector)
        set(gca,'fontsize',20)
        ylim([0 2*pi])
        % s.AlphaData = distfromzero;
        % s.MarkerFaceAlpha = 'flat';
    end
    
end

%% Extra plots 

% (NTOT-N_zero)^2+N_zero
% Th_Aq=Th_Aq+sigma_J*randn(size(Th_Aq));
% Th_Bq=Th_Bq+sigma_J*randn(size(Th_Bq));
% le=linspace(0,1,6);
% figure
% histogram(Eta_Aq,le,'Normalization','pdf')
% hold on
% histogram(Eta_data_A,le,'Normalization','pdf')
% 
% figure
% histogram(Eta_Bq,le,'Normalization','pdf')
% hold on
% histogram(Eta_data_B,le,'Normalization','pdf')
% 
% le=linspace(-pi,pi,9);
% figure
% histogram(angdiff(Th_Aq,Th_Bq),le,'Normalization','pdf')
% hold on
% histogram(angdiff(Preferred_angle_A,Preferred_angle_B),le,'Normalization','pdf')
% % 
% figure
% plot(Th_Aq,Th_Bq,'.','MarkerSize',10)

% label_vector={'0','\pi/2','\pi','3\pi/2','2 \pi'};
% figure
% scatter(Th_Aq,Th_Bq,'filled')
% xlabel('\theta_A','FontSize',20)
% ylabel('\theta_B','FontSize',20)
% set(gca, 'YTick', 0:2*pi/4:2*pi)
% set(gca, 'YTickLabel', label_vector)
% set(gca, 'XTick', 0:2*pi/4:2*pi)
% set(gca, 'XTickLabel', label_vector)
% set(gca,'fontsize',20)
% ylim([0 2*pi])
% s.AlphaData = distfromzero;
% s.MarkerFaceAlpha = 'flat';


%     m=linspace(0,2*pi,8);
%     figure
%     histogram2(Th_Aq,Th_Bq,m,m);
