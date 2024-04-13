clear all 
clc
SNR_dB=[0:5:10];  SNR_linear=10.^(SNR_dB/10.);
N_iter=100; 
N1 = 21;
N2 = 21;
n = N1*N2; % number of beams (transmit antennas)
K = 1; % number of users
L = 1; % number of paths per user
lamada = 1; % wavelength
N = K; % number of retained RF chains
d = lamada/2;
NMSE = zeros(1,length(SNR_dB));
% sigma2=1/SNR_linear(5);
%U = DFT_UPA(N1,N2);
Q = 256;
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
       H_beam = beamspace_channel_lens(N1,N2,K,L,[],[])    % beamspace channel

        for j = 1 : K
            h = H_beam(:,j);
            B = reshape(h,N1,N2);
            contourf(abs(B));
            % measurement matrix, 3*L = times of measurement
            Phi = (1/sqrt(Q))*(2*(rand(Q,n)<0.5) - 1);
%             Phi = exp(1i*2*pi*rand(96,n));  % analog beamforming
            noise = sqrt(sigma2)*(randn(Q,1)+1i*randn(Q,1))/sqrt(2);
            x = Phi*h + noise;
            h_hat1 = OMP_new(x, Phi, 64*L, 64*L);
            error1 = error1 + norm(h - h_hat1,2)^2/norm(h,2)^2;
            h_hat2 = SSD(x,Phi,L,N1,N2,h);
            error2 = error2 + norm(h - h_hat2,2)^2/norm(h,2)^2;
            
%             H_beam_hat1(:,j) = h_hat1;
%             H_beam_hat2(:,j) = h_hat2;
        end
%         [Hr1,Hr1_e] = IA_BS(H_beam,H_beam);
%         [Hr2,Hr2_e] = MM_S(H_beam,1,K,H_beam);
%         [Hr3,Hr3_e] = IA_BS(H_beam_hat1,H_beam);
%         [Hr4,Hr4_e] = MM_S(H_beam_hat1,2,K,H_beam);
%         [Hr5,Hr5_e] = IA_BS(H_beam_hat2,H_beam);
%         [Hr6,Hr6_e] = MM_S(H_beam_hat2,2,K,H_beam);
%         %%% ZF precoding
%         %%% Full digital
%         F = H_beam*inv(H_beam'*H_beam + 10^(-6)*eye(K,K));
%         beta = sqrt(K/trace(F*F'));
%         H_eq=H_beam'*F;
%         for k=1:K
%             sum_inf=sum(abs(H_eq(:,k)).^2)-abs(H_eq(k,k))^2;
%             temp=temp+log2(1+abs(H_eq(k,k))^2/(sum_inf+K/(SNR*beta^2)));
%         end
%         %%% perfect CSI, IA-BS
%         F1 = Hr1_e*inv(Hr1_e'*Hr1_e + 10^(-6)*eye(K,K));
%         beta1 = sqrt(K/trace(F1*F1'));
%         H_eq1=Hr1'*F1;
%         for k=1:K
%             sum_inf=sum(abs(H_eq1(:,k)).^2)-abs(H_eq1(k,k))^2;
%             temp1=temp1+log2(1+abs(H_eq1(k,k))^2/(sum_inf+K/(SNR*beta1^2)));
%         end
%         %%% perfect CSI, MM-BS
%         F2 = Hr2_e*inv(Hr2_e'*Hr2_e + 10^(-6)*eye(K,K));
%         beta2 = sqrt(K/trace(F2*F2'));
%         H_eq2=Hr2'*F2;
%         for k=1:K
%             sum_inf=sum(abs(H_eq2(:,k)).^2)-abs(H_eq2(k,k))^2;
%             temp2=temp2+log2(1+abs(H_eq2(k,k))^2/(sum_inf+K/(SNR*beta2^2)));
%         end
%         %%% OMP CSI, IA-BS
%         F3 = Hr3_e*inv(Hr3_e'*Hr3_e);
%         beta3 = sqrt(K/trace(F3*F3'));
%         H_eq3=Hr3'*F3;
%         for k=1:K
%             sum_inf=sum(abs(H_eq3(:,k)).^2)-abs(H_eq3(k,k))^2;
%             temp3=temp3+log2(1+abs(H_eq3(k,k))^2/(sum_inf+K/(SNR*beta3^2)));
%         end
%         %%% OMP CSI, MM-BS
%         F4 = Hr4_e*inv(Hr4_e'*Hr4_e);
%         beta4 = sqrt(K/trace(F4*F4'));
%         H_eq4=Hr4'*F4;
%         for k=1:K
%             sum_inf=sum(abs(H_eq4(:,k)).^2)-abs(H_eq4(k,k))^2;
%             temp4=temp4+log2(1+abs(H_eq4(k,k))^2/(sum_inf+K/(SNR*beta4^2)));
%         end
%         %%% SD CSI, IA-BS
%         F5 = Hr5_e*inv(Hr5_e'*Hr5_e);
%         beta5 = sqrt(K/trace(F5*F5'));
%         H_eq5=Hr5'*F5;
%         for k=1:K
%             sum_inf=sum(abs(H_eq5(:,k)).^2)-abs(H_eq5(k,k))^2;
%             temp5=temp5+log2(1+abs(H_eq5(k,k))^2/(sum_inf+K/(SNR*beta5^2)));
%         end 
%         %%% SD CSI, MM-BS
%         F6 = Hr6_e*inv(Hr6_e'*Hr6_e);
%         beta6 = sqrt(K/trace(F6*F6'));
%         H_eq6=Hr6'*F6;
%         for k=1:K
%             sum_inf=sum(abs(H_eq6(:,k)).^2)-abs(H_eq6(k,k))^2;
%             temp6=temp6+log2(1+abs(H_eq6(k,k))^2/(sum_inf+K/(SNR*beta6^2)));
%         end
    end
    NMSE1(i_snr) = error1/K/N_iter;
    NMSE2(i_snr) = error2/K/N_iter;
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