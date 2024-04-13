function [H]=test(N1,N2,K,L,azimuth,elevation,D_y,D_z,A)
% n: number of transmit beams (transmit antennas)
% K: number of users
% L: number of paths, L=1 for LoS, L>1 for multipath
lamada = 1; % wavelength
d = lamada/2; % distance between two antennas
H = zeros(N1*N2,K);
if  isempty(azimuth)
    for k = 1 : K
        %%% gain
%         beta(1:L,k) = sqrt(1/2)*(randn(1,L) + 1i*randn(1,L));

  %beta=randn(L,1);
  
   beta(1,k)=1; %LOS
   beta(2:L,k)=randn(L-1,1);%NLOS
  azimuth=unifrnd(-pi/2,pi/2,L,k);
  elevation=unifrnd(-pi/2,pi/2,L,k);
  
%         beta(1,k) = 1;
% %         beta(1,k) = exp(-1i*2*pi*rand(1)); % gain of the LoS
%         beta(2:L,k) = 0.7; % gain of NLoS
%         %% direction
%         azimuth(1,k) = 0.5;
%         elevation(1,k) = 0.5;
%          azimuth(2,k) = 1.5;
%          elevation(2,k) = 1.5;
%           azimuth(3,k) = 1;
%         elevation(3,k) =  1;
%          azimuth(4,k) = 5;
%          elevation(4,k) = 5;
%          azimuth(5,k) = 2.8;
%          elevation(5,k) = 2.1;
        for l = 1 : L
            H(:,k) = H(:,k) + beta(l,k)*LENS(azimuth(l,k),elevation(l,k),N1,N2,D_y,D_z,A);
            %H = H + beta(l,k)*LENS(azimuth(l,k),elevation(l,k),N1,N2,D_y,D_z,A);
        end
    end
else
%     for k = 1 : K
%         %%% gain
%         beta(1,k) = exp(-1i*2*pi*rand(1)); % gain of the LoS
%         beta(2:L,k) = sqrt(0.1*beta(1))*exp(-1i*2*pi*rand(1,L-1)); % gain of NLoS
%         for l = 1 : L
%             H(:,k) = H(:,k) + beta(l,k)*LENS(azimuth(l,k),elevation(l,k),N1,N2,D_y,D_z,A);
%         end
%     end
end
H = sqrt(N1*N2/L)*H;

%nvar = mean(z(z < max(abs(z) ) ).^2) /10^(SNR_dB(ii) / 10)