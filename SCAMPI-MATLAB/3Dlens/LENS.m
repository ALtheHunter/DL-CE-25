function a = LENS(azimuth,elevation,N1,N2,D_y,D_z,A)
% N1: number of antennas in horizon
% N2: number of antennas in vertical
lamada = 1;
hat_phi_y = sin(azimuth);
hat_phi_z = sin(elevation);
a=zeros(N1,N2);
N=N1*N2;
for m=1:1:N1
    for n=1:1:N2
a(m,n) = sqrt(A)*sinc(m-D_y*hat_phi_y)*sinc(n-D_z*hat_phi_z);
%a(m,n) = sqrt(A)*sinc(m-5)*sinc(n-5); 
    end 
end
a=reshape(a,[N,1]);