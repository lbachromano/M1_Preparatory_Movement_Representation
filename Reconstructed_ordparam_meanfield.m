% Compute the reconstructed order parameters from mean-field equations


Opt    = odeset('Events', @myEvent);
integrand=@(t,y) MF_Integrate_order_parameters(t,y,ExtInputs,param,measure,J);
Y=zeros(length(ExtInputs.C0),length(fieldnames(OrdParam)));
IC = structfun(@(x) x(1),OrdParam,'Uni',false);
Y(1,:)=[IC.r_zero; IC.rA; IC.rB; IC.r0A; IC.r0B]';

for time_ind=2:length(ExtInputs.C0)
Y(time_ind,:)=Y(time_ind-1,:)+param.dt_simulaz.*(integrand(time_ind-1,Y(time_ind-1,:)))';
end

%%


x=(1:1:param.t_end)';
time_simul=[0:param.dt_simulaz:param.t_end+1];
figure
clf
subplot(3,1,1)
y = OrdParam.r_zero;
dy = Err_OrdParam.st0;  % made-up error values
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.95 .95 1.0],'linestyle','none');
line(x,y)
hold on
hh=plot(y,'k','LineWidth',3);
hold on
plot(time_simul,Y(:,1),'k','LineWidth',1);
legend(hh,'r_0','Location','SouthEast');
legend('boxoff')
set(gca,'FontSize',15)
subplot(3,1,2)
y = OrdParam.r0A;
dy = Err_OrdParam.st0A;   
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[1.0 .9 0.9],'linestyle','none');
line(x,y)
hold on
h(1)=plot(y,'Color',color_map(4, :),'LineWidth',3);
hold on
y = OrdParam.r0B;
dy = Err_OrdParam.st0B;   
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.8 0.85 0.8],'linestyle','none');
line(x,y)
hold on
h(2)=plot(y,'Color',color_map(5, :),'LineWidth',3);
hold on
plot(time_simul,Y(:,4),'Color',color_map(4, :),'LineWidth',1);
hold on
plot(time_simul,Y(:,5),'Color',color_map(5, :),'LineWidth',1);
legend([h(1) h(2)],'r_{0A}','r_{0B}','Location','SouthEast');
legend('boxoff')
set(gca,'FontSize',15)
subplot(3,1,3)
y = OrdParam.rA;
dy = Err_OrdParam.stA;   
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.9 .9 1.0],'linestyle','none');
line(x,y)
hold on
h(3)=plot(y,'Color',color_map(2, :),'LineWidth',3);
hold on
y = OrdParam.rB;
dy = Err_OrdParam.stB;  
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[1.0 .9 1.0],'linestyle','none');
line(x,y)
hold on
h(4)=plot(y,'Color',color_map(3, :),'LineWidth',3);
hold on
plot(time_simul,Y(:,2),'Color',color_map(2, :),'LineWidth',1);
hold on
plot(time_simul,Y(:,3),'Color',color_map(3, :),'LineWidth',1);
xlabel('time (ms)')
legend([h(3) h(4)],'r_A','r_B','Location','SouthEast');
legend('boxoff')
set(gca,'FontSize',15)

