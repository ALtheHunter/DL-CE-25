function [ F_a_obtained, F_c_obtained, prior] = Fun_ac_CGaussian(S2, R, prior)  
% (67), (68) in the long version
% rho : fraction of non zero components of the original signal
% S2, R : arguments of f_a
% mean_gauss : expectetion of the gaussian taken in the probability measure
% sg2 : it's variance 

% Gaussian sparse prior : p(x) ~ (1 - rho) * delta(x) + rho / sqrt(2 * pi * var) * exp(-(x - m)^2 / (2 * var) )
R_ = R; 
S2_ = S2; 

rho_ = 1; 
m_ = prior.av_mess(1 : prior.N); 
var_ = prior.var_mess(1 : prior.N);   
N_ = prior.N ;  

% rho = 1 ./ (1 + exp(-.5 * prior.R(1 + prior.N : end).^2 ./ prior.S2(1 + prior.N : end) + prior.omega) );


    
a = exp(-abs(R_).^2 ./ S2_ + abs(m_ - R_).^2 ./ (var_ + S2_) );
c = 1 ./ (pi .* (S2_ + var_) );
Z = (1 - rho_) ./ (pi .* S2_) .* a + rho_ .* c;

F_a_obtained = (rho_ .* c .* (m_ .* S2_ + R_ .* var_) ./ (S2_ + var_) ) ./ Z;
F_b = (rho_ ./ (pi .* S2_ .* var_) .* (S2_.^(-1) + var_.^(-1) ).^(-2) .* (1 + abs(m_ ./ var_ + R_ ./ S2_).^2 ./ (1 ./ S2_ + 1 ./ var_) ) ) ./ Z;
F_c_obtained = max(1e-18, F_b - abs(F_a_obtained).^2);


% if (prior.learn == 1)
%     
% %     est_mean_ = sum(F_a_obtained)/(rho_*N_);
%     est_sg2_ = sum(F_c_obtained + abs(F_a_obtained).^2)/(rho_*N_) - abs(m_)^2;  
%     est_rho_ = sum(((1./S2+1/var_)./(R./S2+m_/var_)).*F_a_obtained) ...
%         / sum( 1./(1-rho_+rho_./(1+var_./S2).*exp( abs(R+m_/var_).^2./(1./S2+1/var_) - abs(m_)^2/var_ ) ) );    
% 
%     est_rho_ = real(est_rho_) ;  
%     
% %     prior.mean = prior.dump_mes*prior.mean + (1-prior.dump_mes)* est_mean_ ;    
%     prior.sg2 = prior.dump_mes*prior.sg2 + (1-prior.dump_mes)* est_sg2_ ;
%     prior.rho = min(prior.dump_mes*prior.rho + (1-prior.dump_mes)* est_rho_, prior.alpha) ;   
%      
% end


end