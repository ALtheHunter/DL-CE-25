function [X_hat, Diffs_hat,prior] = oamp_solver(Y, op, opT, opSq, opSqT, opt,snr,A,sigPow,N,l,nvar)

% initialisation 

Y = Y(:); 
opt.signal = opt.signal(:);


R_init = ones(opt.N, 1); 
S2_init = R_init;
S2_init(1 : opt.Ntrue) = opt.var_noise(1); 
S2_init(1 + opt.Ntrue : end) = opt.var_noise(end); 
av_mess_init = R_init;
var_mess_init = S2_init; 


prior = PriorScampi(opt.Ntrue, R_init, S2_init, av_mess_init, var_mess_init, opt.omega); %≥ı º≈‰÷√

% Starting main code
t = 1;
x_B=zeros(N+l,1);
v_AB=zeros(N+l,1);
x_BA=zeros(N+l,1);
v_BA=sigPow; % initialization almost equals to 1 
sigma=(10.^(-snr/10))*sigPow; 
I=eye(N+l);
Gamma = zeros(N,1);
while t <= opt.nb_iter
    
    av_mess_temp = prior.av_mess;
    var_mess_temp = prior.var_mess;    
    
    % OAMP for the augmented system  
    %W=sigma/v_BA*I + A'*A;
    Gamma = 1./sigPow;
    Gamma_diag = diag(Gamma);
    
    var_1 = inv(Gamma_diag + A'*A ./nvar);
    var_2 = var_1*( A'*Y ./nvar + x_BA); 
    
    prior.S2 = diag(var_1)./(1.-Gamma_diag*diag(var_1));
    prior.R = diag(prior.S2)*(var_2./diag(var_1)-x_BA);
     
     
    prior = Prior(prior);
    
    x_B = prior.av_mess;
    v_B = prior.var_mess;
    
    Gamma_new=1./v_B-1./prior.S2;
    x_BA_new=(x_B./v_B-prior.R./prior.S2);
    Gamma=damping(Gamma,Gamma_new,0.2);
    x_BA=damping(x_BA,x_BA_new,0.2);
    
    % use the Bethe free energy to learn the noise variance terms    
    if opt.learnNoise;       
        par = (Y - op(prior.av_mess) ).^2;        
        newNoise = .5 .* (par + sqrt(par.^2 + 4 .* opSq(prior.var_mess) .* par) );
        newNoise1 = opt.dump_learn * opt.var_noise(1) + (1 - opt.dump_learn) * mean(newNoise(1 : opt.Mtrue) );    
        newNoise2 = opt.dump_learn * opt.var_noise(end) + (1 - opt.dump_learn) * mean(newNoise(1 + opt.Mtrue : end) );                
        opt.var_noise(1 : opt.Mtrue) = newNoise1;
        opt.var_noise(1 + opt.Mtrue : end) = newNoise2;    
    end     

    % mse and convergence        
    if mod(t, opt.print) == 0
        conv_ = mean(abs(av_mess_temp - prior.av_mess) );    
        if (t > 0) & (conv_ < opt.conv)
            disp('converged'); 
            break; 
        end          
        nsnr_ = 10 * log10(norm(opt.signal(:) - prior.av_mess(1 : opt.Ntrue) ).^2 ./ norm(opt.signal).^2);       
        disp(sprintf('t=%d, convergence=%0.2e, pixel noise var=%0.2e, dual noise var=%0.2e, nsnr=%0.2e', [t, conv_, opt.var_noise(1), opt.var_noise(end), nsnr_] ) );         
    end   

    % show real time reconstruction
    if opt.showDynamics && (mod(t, opt.print) == 0);
        figure(1); 
        imshow(reshape(prior.av_mess(1 : opt.Ntrue), sqrt(opt.Ntrue), sqrt(opt.Ntrue) ) ); 
        drawnow;
    end

    t = t + 1;    
end

X_hat = prior.av_mess(1 : opt.Ntrue);
Diffs_hat = prior.av_mess(opt.Ntrue + 1 : end);

end