clear all 
clc
SNR_dB=[-20:10:20];  SNR_linear=10.^(SNR_dB/10.);
N_iter=5; 
f=8;
N1 = 16*f;
N2 = 16*f;
n = N1*N2; % number of beams (transmit antennas)
K = 4; % number of users
L = 3; % number of paths per user
lamada = 1; % wavelength
N = K; % number of retained RF chains
d = lamada/2;
D_y = 2.5*f;%equivalent lens length
D_z = 2.5*f;%equivalent lens height
A=20;% length aperture
NMSE = zeros(1,length(SNR_dB));
V1 = 8;
V2 = 8;
% sigma2=1/SNR_linear(5);
%U = DFT_UPA(N1,N2);
Q = 1638;
for i_snr=1:length(SNR_dB) 
    i_snr
    sigma2=1/SNR_linear(i_snr);
%     SNR=SNR_linear(i_snr);
    temp = 0; temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0; temp5 = 0; temp6 = 0;
    error1 = 0; error2 = 0; error3 = 0; error4 = 0; 
    for iter = 1:N_iter;
        iter
%         H = beamspace_channel_UPA(N1,N2,K,L,[],[]);
%         H_beam = U'*H; % beamspace channel
        H_beam = beamspace_channel_lens(N1,N2,K,L,[],[],D_y,D_z,A);    % beamspace channel
        for j = 1 : K
            h = H_beam(:,j);
            B = reshape(h,N1,N2);
            contourf(abs(B));
            % measurement matrix, 3*L = times of measurement
            Phi = (1/sqrt(Q))*(2*(rand(Q,n)<0.5) - 1);
%             Phi = exp(1i*2*pi*rand(96,n));  % analog beamforming
            noise = sqrt(sigma2)*(randn(Q,1)+1i*randn(Q,1))/sqrt(2);
            x = Phi*real(h) + real(noise);
            h_hat1 = OMP_new(x, Phi, 42*L, 42*L);
            [h_hat2 n1 n2 support r] = SSD_simplify(x,Phi,L,N1,N2,h,V1,V2);
            h_wave1 = reshape(h_hat1,N1,N2);
            h_wave2 = reshape(h_hat2,N1,N2);
            error1 = error1 + norm(h - h_hat1,2)^2/norm(h,2)^2;
            error2 = error2 + norm(real(h) - h_hat2,2)^2/norm(real(h),2)^2;
            error3 = error3 + norm(h - r,2)^2/norm(h,2)^2;
%             H_beam_hat1(:,j) = h_hat1;
%             H_beam_hat2(:,j) = h_hat2;
        end
%         
    end
    NMSE1(i_snr) = error1/K/N_iter;
    NMSE2(i_snr) = error2/K/N_iter;
    NMSE3(i_snr) = error3/K/N_iter;
%     C(i_snr) = temp/N_iter;
%     C1(i_snr) = temp1/N_iter;
%     C2(i_snr) = temp2/N_iter;
%     C3(i_snr) = temp3/N_iter;
%     C4(i_snr) = temp4/N_iter;
%     C5(i_snr) = temp5/N_iter;
%     C6(i_snr) = temp6/N_iter;
end
figure(1)
semilogy(SNR_dB,NMSE1,'b','Linewidth',1.5);
hold on 
semilogy(SNR_dB,NMSE2,'r','Linewidth',1.5);
% hold on 
% semilogy(SNR_dB,NMSE3,'y','Linewidth',1.5);
grid on
xlabel('SNR (dB)');
ylabel('NMSE (dB)');
% figure(2)
% semilogy(SNR_dB,C,'k','Linewidth',1.5);
% hold on 
% % semilogy(SNR_dB,C1,'y','Linewidth',1.5);
% hold on 
% semilogy(SNR_dB,C2,'b','Linewidth',1.5);
% hold on 
% % semilogy(SNR_dB,C3,'r','Linewidth',1.5);
% hold on 
% semilogy(SNR_dB,C4,'g','Linewidth',1.5);
% hold on
% % semilogy(SNR_dB,C5,'p','Linewidth',1.5);
% hold on 
% semilogy(SNR_dB,C6,'g','Linewidth',1.5);
% grid on
% xlabel('SNR (dB)');
% ylabel('Achievable sum-rate (bits/s/Hz)');
% % contourf(abs(H_beam'));