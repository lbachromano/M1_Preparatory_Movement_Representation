
%% integration parameters

dt_simulaz=0.5;  
tau=25; % membrane time constant. does not change the results significantly
T_end=length(ExtInputs.C0);

ang=pi.*(linspace(0,param.Nr_target-1,param.Nr_target)'...
        +linspace(1,param.Nr_target,param.Nr_target)')./8; %location of the target on the screen
ext_ang=ang(ext_dir);
    
%% These parameters set the level of the noise.
% It was chosen so that PCA plots matcht the data
% Our model is instrinsically 4 dimensional. Extra dimensions are due to noise,
% in the current version of the model
noise_fields=0.05;
tau_fields=2;
sigmaq=3;
dt_noise=300*dt_simulaz;

%% This adds noise to the external fields. It is optional
% ANOTHER OPTION, IS TO ADD A GAUSSIAN NOISE TO THE COUPLINGS MATRIX.
% by adding to the matrix "couplings" defined at line 55
% a noise matrix  
% rand(DiscreteVar.NTOT,DiscreteVar.NTOT).*some_variance

Z=Generate_OU_process(DiscreteVar.NTOT,T_end,noise_fields,tau_fields);
ExtInputs.C0=ExtInputs.C0+Z';

Z=Generate_OU_process(DiscreteVar.NTOT,T_end,noise_fields,tau_fields);
ExtInputs.CA=ExtInputs.CA+Z';

Z=Generate_OU_process(DiscreteVar.NTOT,T_end,noise_fields,tau_fields);
ExtInputs.CB=ExtInputs.CB+Z';

Z=Generate_OU_process(DiscreteVar.NTOT,T_end,noise_fields,tau_fields);
ExtInputs.epsA=ExtInputs.epsA+Z';

Z=Generate_OU_process(DiscreteVar.NTOT,T_end,noise_fields,tau_fields);
ExtInputs.epsB=ExtInputs.epsB+Z';


%% Define the couplings matrix and the external fields

couplings=J.j0+J.jsp.*(DiscreteVar.Eta_Aq.*DiscreteVar.Eta_Aq').*cos(DiscreteVar.Th_Aq-DiscreteVar.Th_Aq')+J.jse.*(DiscreteVar.Eta_Bq.*DiscreteVar.Eta_Bq')...
    .*cos(DiscreteVar.Th_Bq-DiscreteVar.Th_Bq')+J.ja.*(DiscreteVar.Eta_Bq.*DiscreteVar.Eta_Aq').*cos(DiscreteVar.Th_Bq-DiscreteVar.Th_Aq');

%%
field_ext=ExtInputs.C0'+ExtInputs.CA'.*DiscreteVar.Eta_Aq'+ExtInputs.CB'.*DiscreteVar.Eta_Bq'+ExtInputs.epsA'.*(DiscreteVar.Eta_Aq'.*cos(DiscreteVar.Th_Aq'-ext_ang))...
    +ExtInputs.epsB'.*(DiscreteVar.Eta_Bq'.*cos(DiscreteVar.Th_Bq'-ext_ang));

%%

noise=zeros(1,DiscreteVar.NTOT);
Y=zeros(T_end,DiscreteVar.NTOT);

Initial_X=OrdParam.r_zero(1).*ones(1,DiscreteVar.NTOT)+0.2*cos((DiscreteVar.Th_Aq')-ext_ang);
Y(1,:)=Initial_X;

tic
for j = 2:T_end
    funz_here=myODE_simulaz(Y(j-1,:),couplings,field_ext(j,:),tau,noise,DiscreteVar.NTOT);
    Y(j,:) =Y(j-1,:)+dt_simulaz.*funz_here;
    noise=noise+(-noise+sigmaq.*sqrt(dt_noise).*randn(1,DiscreteVar.NTOT))/dt_noise;
end
toc
 


