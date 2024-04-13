function [Xmod_r,v]=GFDM_demodulate_OAMP_complex(yr,sigPow,A,snr,N)
   %------------------------------------------------initialization
   x_B=zeros(N,1);
   v_AB=0;
   x_BA=zeros(N,1);
   v_BA=sigPow; % initialization almost equals to 1 
   sigma=(10.^(-snr/10))*sigPow; 
   I=eye(N);
   target=3; %iteration
   d=10;
   %temp=Xmod.';
   iter=0;   
   prior.N = N;                                                        
        prior.omega = 0;
        prior.rho =1;
        prior.dump_mes = 1e-20;
        prior.mean = 0;
        prior.sg2 = 1;
        prior.learn = 1;
   %------------------------------------------------Loop
   while iter<target
        W=sigma/v_BA*I + A'*A;    
        W=inv(W)*A';           
        
        eta=real(N/trace(W*A));  
        x_AB=x_BA+eta*W*(yr-A*x_BA);
        v_AB=v_BA*(eta-1);  
   
        %[x_B,v_B]=estim_cond_complex(Xmod,x_AB,v_AB);

    
        [x_B,v_B,c] = Fun_ac_RGaussianLearn(v_AB, x_AB, prior);
       
        x_B = x_B(1:N);
        v_B = mean(v_B(1:N));
        v_BA_new=1./v_B-1/v_AB;
        v_BA_new=1./v_BA_new; 
        if v_BA_new<0
            v_BA_new=v_BA;
        end
   
        x_BA_new=v_BA_new.*(x_B./v_B-x_AB./v_AB);
        
        v_BA=damping(v_BA,v_BA_new,0.2);
        x_BA=damping(x_BA,x_BA_new,0.2);
        
        d_new=sum(abs(bsxfun(@minus, x_B, Xmod.')).^2)/N; 
%         if d_new>d || iter>20
% 
%             x_B=temp;
%             break;
%         else
%             d=d_new;
%         end
%         
%         temp=x_B;
        iter=iter+1;
       
   end 
   
   Xmod_r = x_B;
   v = v_AB;
end