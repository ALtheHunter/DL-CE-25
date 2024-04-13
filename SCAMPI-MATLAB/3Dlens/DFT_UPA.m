function U = DFT_UPA(N1,N2)
% N1: number of antennas in horizon
% N2: number of antennas in vertical
deta1 = 1/N1; % resolution in horizon
deta2 = 1/N2; % resolution in vertical
for i1 = -(N1-1)/2:1:(N1-1)/2
    i = i1 + (N1-1)/2 + 1;
    for i2 = -(N2-1)/2:1:(N2-1)/2
        j = i2 + (N2-1)/2 +1;
        temp1 = sqrt(1/N1)*exp(1i*[-(N1-1)/2:1:(N1-1)/2]*2*pi*i1*deta1).';
        temp2 = sqrt(1/N2)*exp(1i*[-(N2-1)/2:1:(N2-1)/2]*2*pi*i2*deta2).';
        U(:,N2*(i-1)+j) = kron(temp1,temp2);
    end
end

% lamada = 1;
% d1 = lamada/2; % antenna spacing in horizon
% d2 = lamada/2;  % antenna spacing in vertical
% fre_sp1 = (d1/lamada)*cos(elevation)*sin(azimuth);
% fre_sp2 = (d2/lamada)*sin(elevation);
% a1 = sqrt(1/N1)*exp(1i*[-(N1-1)/2:1:(N1-1)/2]*2*pi*fre_sp1).';
% a2 = sqrt(1/N2)*exp(1i*[-(N2-1)/2:1:(N2-1)/2]*2*pi*fre_sp2).';
% a = kron(a1,a2);
% 
%             U = zeros(n,n); 
%             deta = 1/n;
%             for i = -(n-1)/2:1:(n-1)/2
%                 U(:,i+(n+1)/2) = sqrt(1/n)*exp(1i*[0:n-1]*2*pi*deta*i).'; % spatial DFT
%             end