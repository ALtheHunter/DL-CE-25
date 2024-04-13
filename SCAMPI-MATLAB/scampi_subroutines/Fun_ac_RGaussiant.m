function [ F_a_obtained, F_v_obtained, prior] = Fun_ac_RGaussiant(S2, R, prior)  
% (67), (68) in the long version
% rho : fraction of non zero components of the original signal
% S2, R : arguments of f_a
% mean_gauss : expectetion of the gaussian taken in the probability measure
% sg2 : it's variance 

% Gaussian sparse prior : p(x) ~ (1 - rho) * delta(x) + rho / sqrt(2 * pi * var) * exp(-(x - m)^2 / (2 * var) )
R_ = R; 
S2_ = S2; 


m_ = prior.av_mess(1 : prior.N); 
var_ = prior.var_mess(1 : prior.N);   
N_ = prior.N ; 




F_a_obtained = (m_ .* S2_ + R_ .* var_) ./ (S2_ + var_) ;
F_v_obtained = (S2_ .* var_) ./ (S2_ + var_) ;
 


end