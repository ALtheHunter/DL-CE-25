classdef PriorScampi    
    properties
        av_mess; 
        var_mess; 
        R; 
        S2; 
        N;                 
        omega;
        rho;
        dump_mes;
        mean;
        sg2;
        learn;
    end
    
    methods
        
        function prior = PriorScampi(N, R_init, S2_init, av_mess_init, var_mess_init, omega)                                      
            % Constructor function
            prior.R = R_init;
            prior.S2 = S2_init; 
            prior.N = N;             
            prior.av_mess = av_mess_init; 
            prior.var_mess = var_mess_init;                                                               
            prior.omega = omega;
            prior.rho = 1;
            prior.dump_mes = 1e-20;
            prior.mean = 0;
            prior.sg2 = 1;
            prior.learn = 1;
           
        end
        
        function prior = Prior(prior)
            % pixels (Gaussian, i.e. no prior)      
%             prior.av_mess(1 : prior.N) = prior.R(1 : prior.N);                        
%             prior.var_mess(1 : prior.N) = prior.S2(1 : prior.N);            
   
            % pixels (Gaussian)
            [a,b,c] = Fun_ac_RGaussianLearn(prior.S2(1 : prior.N), prior.R(1 : prior.N), prior);
            prior.av_mess(1 : prior.N) = a;
            prior.var_mess(1 : prior.N) = b;
            
            % differences (SNIPE)                                    
            rho = 1 ./ (1 + exp(-.5 * prior.R(1 + prior.N : end).^2 ./ prior.S2(1 + prior.N : end) + prior.omega) );
            prior.av_mess(1 + prior.N : end) = rho .* prior.R(1 + prior.N : end);              
            prior.var_mess(1 + prior.N : end) = rho .* (prior.R(1 + prior.N : end).^2 .* (1 - rho) + prior.S2(1 + prior.N : end) );                                                                       
        end
        
        function prior = PriorU(prior)
            % pixels (Uniform, i.e. no prior)      
             prior.av_mess(1 : prior.N) = prior.R(1 : prior.N);                        
             prior.var_mess(1 : prior.N) = prior.S2(1 : prior.N);            
   
            % pixels (Gaussian)
%             [a,b,c] = Fun_ac_RGaussian(prior.S2(1 : prior.N), prior.R(1 : prior.N), prior);
%             prior.av_mess(1 : prior.N)=a;
%             prior.var_mess(1 : prior.N)=b;
%             
            % differences (SNIPE)                                    
            rho = 1 ./ (1 + exp(-.5 * prior.R(1 + prior.N : end).^2 ./ prior.S2(1 + prior.N : end) + prior.omega) );
            prior.av_mess(1 + prior.N : end) = rho .* prior.R(1 + prior.N : end);              
            prior.var_mess(1 + prior.N : end) = rho .* (prior.R(1 + prior.N : end).^2 .* (1 - rho) + prior.S2(1 + prior.N : end) );                                                                       
        end
    end    
end