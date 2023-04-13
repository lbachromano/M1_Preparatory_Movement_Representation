 
F0=1;
FA=moments.average_eta_A_squared/2;
FB=moments.average_eta_B_squared/2;
FAB=moments.average_eta_A*moments.average_eta_B/6;

%%

Pts_plots=10000;
vector_ja=[0; 10;20]; %Pick three values of ja to visualize
vector_jB=linspace(0,46,Pts_plots);
figure
for j=1:length(vector_ja)
    ja=vector_ja(j)
    for i=1:Pts_plots
        jsA=vector_jB(i);  
        lam2 = @(jsB) 0.5* (-2 + FAB*ja + FA*jsA + FB*jsB + ...
           sqrt((FA*jsA - FB*jsB)^2 + 2*FAB*ja*(FA*jsA + FB*jsB) + ...
            FAB^2*(ja^2 + 4*jsA*jsB))); 
        sol(i,1) = jsA;
        sol(i,2) =  fzero(lam2, 10);
    end
    plot(sol(:,1),sol(:,2),'LineWidth',3)
    hold on
end
legend('j_a=0','j_a=10','j_a=20','FontSize',15)
ylim([0,max(sol(:,2))+10])
xlabel('j_s^A')
ylabel('j_s^B')
set(gca,'fontsize',20)

%     Color the area 
%     figure
%     a=area(sol(:,1),sol(:,2))
%     a.FaceColor = [0.9290 0.6940 0.1250];
%     xlabel('j_s^A','FontSize',40)
%     ylabel('j_s^B','FontSize',40)
%     ylim([0 120]);
%     xlim([0 60]);
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])


