function [H]=beamspace_channel_UPA(N1,N2,K,L,azimuth,elevation)
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
        beta(1,k) = 1;
%         beta(1,k) = exp(-1i*2*pi*rand(1)); % gain of the LoS
        beta(2:L,k) = sqrt(0.1*beta(1))*exp(-1i*2*pi*rand(1,L-1)); % gain of NLoS
        %%% direction
        azimuth(1,k) = asin(-6.95/16);
        elevation(1,k) = asin(4.5/16);
%         azimuth(1:L,k) = pi*rand(1,L) - pi/2;
%         elevation(1:L,k) = pi*rand(1,L) - pi/2;
        for l = 1 : L
            H(:,k) = H(:,k) + beta(l,k)*UPA(azimuth(l,k),elevation(l,k),N1,N2);
        end
    end
else
    for k = 1 : K
        %%% gain
        beta(1,k) = exp(-1i*2*pi*rand(1)); % gain of the LoS
        beta(2:L,k) = sqrt(0.1*beta(1))*exp(-1i*2*pi*rand(1,L-1)); % gain of NLoS
        for l = 1 : L
            H(:,k) = H(:,k) + beta(l,k)*UPA(azimuth(l,k),elevation(l,k),N1,N2);
        end
    end
end
H = sqrt(N1*N2/L)*H;