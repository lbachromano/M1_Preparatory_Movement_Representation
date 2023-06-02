
function dydt = myODE_shorter_simulaz(y,couplings,field_ext,tau,rumore,NTOT)

   temp_field_here=rumore+(couplings*y')'./NTOT+field_ext;   
   temp_field_here(temp_field_here <0)=0;
   dydt = (-y+temp_field_here)./tau;
   
end


