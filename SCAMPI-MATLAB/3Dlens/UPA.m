function a = UPA(azimuth,elevation,N1,N2)
% N1: number of antennas in horizon
% N2: number of antennas in vertical
lamada = 1;
d1 = lamada/2; % antenna spacing in horizon
d2 = lamada/2;  % antenna spacing in vertical
fre_sp1 = (d1/lamada)*sin(azimuth);
% fre_sp1 = (d1/lamada)*cos(elevation)*sin(azimuth);
fre_sp2 = (d2/lamada)*sin(elevation);
a1 = sqrt(1/N1)*exp(1i*[-(N1-1)/2:1:(N1-1)/2]*2*pi*fre_sp1).';
a2 = sqrt(1/N2)*exp(1i*[-(N2-1)/2:1:(N2-1)/2]*2*pi*fre_sp2).';
a = kron(a1,a2);  