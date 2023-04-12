function [Et,Th]=Prepare_tuning(Activity)
    Th=zeros(size(Activity,1),1);
    Et=zeros(size(Activity,1),1);
    for u=1:size(Activity,1)
        Ct=squeeze(mean(Activity(u,:,:),2));
        b0=mean(Ct);
        b1=((Ct(2)+Ct(4)-Ct(6)-Ct(8))/sqrt(2) + Ct(3)-Ct(7))/4;
        b2=((Ct(2)-Ct(4)-Ct(6)+Ct(8))/sqrt(2) + Ct(1)-Ct(5))/4;
        if (b2==0)    
            if(b1>0)
            Inv_t=pi/2;
            elseif (b1==0)
            Inv_t=NaN;
            else
             Inv_t=3*pi/2;
            end
            else
            Inv_t=atan(b1/b2);
        end
        if ( b2>=0 && b1>=0 )
            pt=Inv_t;       
            elseif  (b2>=0 && b1<0)
                pt=Inv_t+2*pi;
            elseif  ( b2<0)
                pt=Inv_t+pi;
        end
        direction_max = mod(pt,2*pi);
        Th(u,1)=direction_max;
        Et(u,1)=2*sqrt(b1*b1+b2*b2);
    end
    Et=Et./max(Et);
    Et(Et==0)=0.001;
    
end
