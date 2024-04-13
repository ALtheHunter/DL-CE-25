function [H]=beamspace_channel_lens_test(N1,N2,K,L,azimuth,elevation,D_y,D_z,A)
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
        beta(2:L,k) = 0.7; %beta(2:L,k) = sqrt(0.1*beta(1))*exp(-1i*2*pi*rand(1,L-1)); % gain of NLoS
        %%% direction
        azimuth(1,k) = (pi*2*k)/18;
        elevation(1,k) = (pi*2*k)/12;
         azimuth(2,k) = (pi*2*k)/19;
         elevation(2,k) = (pi*2*k)/13;
         azimuth(3:L,k) = (pi*2*k)/20;
         elevation(3:L,k) = (pi*2*k)/14;
        for l = 1 : L
            H(:,k) = H(:,k) + beta(l,k)*LENS(azimuth(l,k),elevation(l,k),N1,N2,D_y,D_z,A);
        end
    end
else
    for k = 1 : K
        %%% gain
        beta(1,k) = exp(-1i*2*pi*rand(1)); % gain of the LoS
        beta(2:L,k) = sqrt(0.1*beta(1))*exp(-1i*2*pi*rand(1,L-1)); % gain of NLoS
        for l = 1 : L
            H(:,k) = H(:,k) + beta(l,k)*LENS(azimuth(l,k),elevation(l,k),N1,N2,D_y,D_z,A);
        end
    end
end
H = sqrt(N1*N2/L)*H;