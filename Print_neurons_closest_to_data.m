
Number_selected=1; %number of closest neurons, for each one in the data
Closest_index=zeros(length(Variable.Eta_data_A),Number_selected);
for i=1:length(Variable.Eta_data_A)
   i
 
   d_ea=(DiscreteVar.Eta_Aq-Variable.Eta_data_A(i)).^2;
   d_eb=(DiscreteVar.Eta_Bq-Variable.Eta_data_B(i)).^2;
   d_tha=(angdiff(DiscreteVar.Th_Aq,ones(size(DiscreteVar.Th_Aq)).*Variable.Preferred_angle_A(i)))'./pi;
   d_thb=(angdiff(DiscreteVar.Th_Bq,ones(size(DiscreteVar.Th_Aq)).*Variable.Preferred_angle_B(i)))'./pi;
   dist=d_ea+d_eb +d_tha.^2+d_thb.^2;

   [mind idx_sort]=sort(dist);
   Closest_index(i,:)=idx_sort(1:Number_selected)';
end

%%
one_dim_idx=reshape(Closest_index',[],1);
selected_traj=Y(:,one_dim_idx); % THESE ARE THE TRAJECTORIES  
selected_ea=DiscreteVar.Eta_Aq(one_dim_idx); % THESE IS ETA_A FOR THE SELECTED NURONS
selected_eb=DiscreteVar.Eta_Bq(one_dim_idx);
selected_tha=DiscreteVar.Th_Aq(one_dim_idx);
selected_theb=DiscreteVar.Th_Bq(one_dim_idx);


fullFileName =  sprintf('Closest_traj%d_N%d.mat',ext_dir,Network);
save(fullFileName,'selected_traj','selected_ea','selected_eb','selected_tha','selected_theb');
