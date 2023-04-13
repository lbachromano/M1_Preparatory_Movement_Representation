function [rho_4d,moments] = Distribution_all_variables(Variable,color_map,plotid)

    % Distribution of eta variables from kernel density estimation
    nr_bins=160;
    lst=linspace(0.0,1,nr_bins);
    [bandwidth,densityA,xmeshA,cdf]=kde(Variable.Eta_data_A,2^5,0,1);
    Interpl_dens_A=interp1(xmeshA,densityA,lst);
    distrib_etaA=gaussfilt(lst,Interpl_dens_A,0.01);
    [bandwidth,densityB,xmeshB,cdf]=kde(Variable.Eta_data_B,2^6,0,1);
    Interpl_dens_B=interp1(xmeshB,densityB,lst);
    distrib_etaB=gaussfilt(lst,Interpl_dens_B,0.01);

    if plotid==1
        x_increm=diff(lst);
        x_i=x_increm(1);
        figure
        subplot(2,1,1)
        plot(lst,distrib_etaA./(x_i*sum(distrib_etaA)),'Color',color_map(2, :),'LineWidth',3);
        xlabel('\eta_A')
        ylabel('\rho_{pA}')
        set(gca,'fontsize',18)
        subplot(2,1,2)
        plot(lst,distrib_etaB./(x_i*sum(distrib_etaB)),'Color',color_map(3, :),'LineWidth',3);
        xlabel('\eta_B')
        ylabel('\rho_{pB}')
        set(gca,'fontsize',18)
    end

    %% Compute the distribution that I will use for the 4-dim integration.
    % I will use the trapezoid method to integrate over the 4-dim distribution
    % of variables eta_A, eta_B, theta_A, theta_B
    % Here I define the measure that I will use later 
    NB_angles_h=8; %This is the minimum number of bins that yields accurate results for the 4-dim integration over dx. 
                   %A larger number of bins will yield more accurate integrals,
                   %but larger integration time
    NB_eta_h=12;
    dx=(2*pi./(NB_angles_h)).^2*(1./NB_eta_h).^2;
    dim=(NB_angles_h+1)^2*(NB_eta_h+1)^2;
    wA=ones(NB_angles_h+1,1); %for trapezoid method
    wE=ones(1,1,1,NB_eta_h+1);
    wA(1)=0.5;
    wA(end)=0.5;
    wE(1,1,1,1)=0.5;
    wE(1,1,1,end)=0.5;
    edges_ang=linspace(0,2*pi,NB_angles_h+1);
    edges_eta=linspace(0,1,NB_eta_h+1);
    rho_eta_A=interp1(lst,distrib_etaA,edges_eta);
    rho_eta_B=interp1(lst,distrib_etaB,edges_eta);

    rho_4d=(1.5+cos(edges_ang - edges_ang')).*permute(rho_eta_A,[1,3,2]).*permute(rho_eta_B,[1,3,4,2]);  
    %Distribution for the angular variables is proportional to (1.5+cos(\theta_A - \theta_B))
    Normalization = sum((rho_4d.*wA.*wA'.*wE.*permute(wE,[1,2,4,3])).*dx,"all");
    rho_4d=bsxfun(@rdivide, rho_4d, Normalization);
    %4-dimensional measure for integration with trapezoid method

    %% I compute here the first and second moments of eta_A and eta_B
    % I will use these to compute the stability
    moments.average_eta_A=sum(dx.*rho_4d.*wA.*wA'.*wE.*permute(wE,[1,2,4,3]).*permute(edges_eta,[1,3,2]),"all");
    moments.average_eta_B=sum(dx.*rho_4d.*wA.*wA'.*wE.*permute(wE,[1,2,4,3]).*permute(edges_eta,[1,3,4,2]),"all");
    moments.average_eta_A_squared=sum(dx.*rho_4d.*wA.*wA'.*wE.*permute(wE,[1,2,4,3]).*permute(edges_eta,[1,3,2]).*permute(edges_eta,[1,3,2]),"all");
    moments.average_eta_B_squared=sum(dx.*rho_4d.*wA.*wA'.*wE.*permute(wE,[1,2,4,3]).*permute(edges_eta,[1,3,4,2]).*permute(edges_eta,[1,3,4,2]),"all");
    moments.average_eta_A_eta_B=sum(dx.*rho_4d.*wA.*wA'.*wE.*permute(wE,[1,2,4,3]).*permute(edges_eta,[1,3,2]).*permute(edges_eta,[1,3,4,2]),"all");
    % Check: the last one is the product of the first two

end