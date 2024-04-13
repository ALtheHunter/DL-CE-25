function [ F_a_obtained, F_c_obtained, prior] = Fun_ac_RGaussian(S2, R, prior)  
% (67), (68) in the long version
% rho : fraction of non zero components of the original signal
% S2, R : arguments of f_a
% mean_gauss : expectetion of the gaussian taken in the probability measure
% sg2 : it's variance 

% Gaussian sparse prior : p(x) ~ (1 - rho) * delta(x) + rho / sqrt(2 * pi * var) * exp(-(x - m)^2 / (2 * var) )
R_ = R; 
S2_ = S2; 

rho_ = 1 ;
m_ = prior.av_mess(1 : prior.N); 
var_ = prior.var_mess(1 : prior.N);   
N_ = prior.N ; 


a = exp(-R_.^2 ./ (2 .* S2_) + 0.5 .* (m_ - R_).^2 ./ (var_ + S2_) );
c = 1 ./ sqrt(2 .* pi .* (S2_ + var_) );
Z = (1 - rho_) ./ sqrt(2 .* pi .* S2_) .* a + rho_ .* c;
F_a_obtained = (rho_ .* c .* (m_ .* S2_ + R_ .* var_) ./ (S2_ + var_) ) ./ Z;
F_b = (rho_ ./ sqrt(2 .* pi .* S2_ .* var_) .* (S2_.^(-1) + var_.^(-1) ).^(-3 ./ 2) .* (1 + (m_ ./ var_ + R_ ./ S2_).^2 ./ (1 ./ S2_ + 1 ./ var_) ) ) ./ Z;
F_c_obtained = max(1e-18, F_b - F_a_obtained.^2);
 

% if (prior.learn == 1)
%     
%     disp('not support for the parameters learning for the real Gaussian prior') 
%      
% end


end