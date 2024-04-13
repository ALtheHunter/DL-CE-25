function [ F_a_obtained, F_c_obtained, prior] = Fun_ac_MixedCGaussian(S2, R, prior)

% Gaussian-mixture prior : p(x) ~ \sum_{i=1}^{L} rho_i / sqrt(2 * pi * var_i) * exp(-|x_i - m_i|^2 / sg2_i)
% rho : weight of the original signal
% sg2 : variance of the original signal
% m : mean of the original signal (not used!!!)
%
% S2, R : arguments of f_a

[in_m, in_n] = size(R) ;


rho_ = zeros(prior.signal_L,1);
var_ = zeros(prior.signal_L,1);
for k = 1:prior.signal_L; rho_(k) = prior.rho(k) ; var_(k) = prior.sg2(k) ;   end;
N_ = prior.N ;


Z = zeros(in_m, in_n);
F_a_obtained = zeros(in_m, in_n); F_b = zeros(in_m, in_n);
a = zeros(in_m, in_n, prior.signal_L); c = zeros(in_m, in_n, prior.signal_L); e = zeros(in_m, in_n, prior.signal_L);
for k = 1:prior.signal_L
    a(:,:,k) = exp(- abs(R).^2 ./ (var_(k) + S2) );
    c(:,:,k) = 1 ./ ( pi .* (S2 + var_(k)) );
    e(:,:,k) = R .* var_(k) ./ (S2 + var_(k));
    
    Z = Z + rho_(k) .* a(:,:,k) .* c(:,:,k) ;
    F_a_obtained = F_a_obtained + rho_(k) .* a(:,:,k) .* c(:,:,k) .* e(:,:,k) ;
    F_b = F_b + rho_(k) ./  (pi .* S2 .* var_(k)) .* (1 ./ var_(k) + 1 ./ S2).^(-2) .* a(:,:,k) .* (1 + abs(R ./ S2).^2 ./ (1 ./ var_(k) + 1 ./ S2) ) ;
end

F_a_obtained = F_a_obtained ./ Z ;
F_b = F_b ./ Z ;
F_c_obtained = max(1e-18, F_b - abs(F_a_obtained).^2);


% if (prior.learn ~= 0)
%     
%     % Using the method by Vila & Schniter ---
%     tho_ = zeros(in_m, in_n, prior.signal_L);
%     for k = 1:prior.signal_L
%         tho_(:,:,k) = ( rho_(k) .* a(:,:,k) .* c(:,:,k) ) ./ Z ;
%         up_rho = sum( tho_(:,:,k) ) ;
%         down_rho = N_ ;
%         prior.rho(k) = prior.dump_mes .* prior.rho(k) + (1 - prior.dump_mes) .* up_rho ./ down_rho;
%         %         if (prior.rho(k) > prior.alpha); prior.rho(k) = prior.alpha ; end;
%     end
%     %     prior.rho(1) = 1 - sum( prior.rho(2:prior.signal_L) ) ;
%     
%     for k = 1:prior.signal_L
%         S2_mixed = ( S2 .* var_(k) ) ./ ( S2 + var_(k) ) ;
%         up_var = sum( tho_(:,:,k) .* (abs(0 - e(:,:,k)).^2 + S2_mixed ) );
%         down_var = sum( tho_(:,:,k) );
%         prior.sg2(k) = prior.dump_mes .* prior.sg2(k) + (1 - prior.dump_mes) .* up_var ./ down_var ;
%     end
    
    %     prior.mean = prior.dump_mes .* prior.mean + (1 - prior.dump_mes) .* 1 ./ (N_ .* rho_) .* sum(F_a_obtained);
    
end

end