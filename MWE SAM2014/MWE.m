%The references to definitions and equations refer (mostly) to the 
%following publication:
%
%Müller, A., A. Kammoun, E. Björnson, and M. Debbah, "Efficient Linear 
%Precoding for Massive MIMO Systems using Truncated Polynomial Expansion", 
%IEEE 8th Sensor Array and Multichannel Signal Processing Workshop (SAM), 
%Coruna, Spain, 2014.
%Available at: http://www.laneas.com/axel-muller
%
%This is version 1.2. (Last edited: 2016-03-16)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.

clear all
close all
disp(' '),disp(' ');
disp('%%%%%%%%%%%%%%%%%%%');
disp('Minimal Working Example: Avg Rate of UT1 only vs. TPE order vs. mean RX SNR');
disp(' ');

%% Prepare Figure
f1 = figure('Name', 'Avg Rate of UT1 only vs. TPE order vs. mean RX SNR', 'units','pixels');
hold on; grid on;
xlabel('Average Receive SNR (dB)','fontsize',16)
ylabel('Achievable Avg. Rate of UT_1 (nats/s/Hz)','fontsize',16)
set(gca,'FontSize',16)


%% Set system parameters and constants
P=1;
SNR=[0:4:20];				% Transmit SNR

sigma_tab=P*10.^(-SNR/10);	
rho_lin_tab=P./sigma_tab;
K=128;                      % Number of Users
M=512;                      % Number of BS antennas
a=0.1;						% Covariance/Correlation between antennas
tau=0.1;                    % CSI quality
L_vec=[2,3,4];              % TPE order(s)

% Monte-Carlo runs
num_iter=10;               % Choose >1000 for print quality

%% Channel Covariance Matrix
Phi=a.^toeplitz([0:M-1]);
sqrtmPhi = sqrtm(Phi);   

%% Secondary outer loop  (choose varying value)
for L=L_vec
    fprintf('L=%d, M=%d, K=%d\n', L,M,K );
    
    %% Preallocation
    Rate_moy_pol2_tab=nan(1,length(SNR));
    Rate_moy_rzf_tab=nan(1,length(SNR));    
    
    %% Main TxPower Loop
    tic;
    for ii=1:length(sigma_tab)
        Msg = sprintf('Processing %d/%d (SNR=%0.3gdB @tau=%g); last round took %0.2f mins.\n', ii, length(SNR), SNR(ii), tau, toc/60);
        tic; fprintf(Msg);
        
        %% Re-Initialisation
        Rate_moy_pol2=0;
        Rate_moy_rzf=0;
        
        %% Pre-Calculations
        sigma=sigma_tab(ii);
        rho_lin=rho_lin_tab(ii);
        
        % Calculate TPE weights
        gf2=Precoder_polynomial_updated(P,Phi,L,tau,sigma,K);
        csi=fixed_point_solve(Phi,tau,rho_lin,K);
        
        %% Monte-Carlo loop
        for iter=1:num_iter
            %% Create Channel draw with error and covariance
            W=(randn(M,K)+1i*randn(M,K))/sqrt(2*K);
            H=sqrtmPhi*W;
            h_1=H(:,1);
            error_channel=(randn(M,K)+1i*randn(M,K))/sqrt(2*K);
            H_hat = sqrt(1-tau^2)*H + tau*sqrtmPhi*error_channel;
                        
            %% TPE precoder
            G_pol2=zeros(M,K);
            HH=eye(K);
            for l=0:L-1
                G_pol2=G_pol2+gf2(l+1)*H_hat*HH;
                HH=H_hat'*H_hat*HH;
            end
            % Power re-normalization unneeded IFF weights are chosen correctly!
%             frob = sqrt(1/K*trace(G_pol2*G_pol2'));   % frob == 1.00;
%             fprintf('frob=%0.2f\n',frob);
%             G_pol2=G_pol2/frob;                       % (Unneeded) Power normalization polynomial precoder
            
            %% RZF Precoder
            G_rzf=Precoder_regularized_zero_forcing(P,H_hat,csi);	%Function returns version of RZF precoder that is power normalized to P
            
            %% Rate of UT1
            g_1_rzf=G_rzf(:,1);
            G_1_rzf=G_rzf(:,2:end);
            
            g_12=G_pol2(:,1);
            G_12=G_pol2(:,2:end);
            
            SINR_pol2=abs(g_12'*h_1)^2/(norm(G_12'*h_1).^2+sigma);
            SINR_rzf=abs(g_1_rzf'*h_1)^2/(norm(G_1_rzf'*h_1)^2+sigma);
            
            Rate_moy_pol2=Rate_moy_pol2+abs(log(1+SINR_pol2));
            Rate_moy_rzf=Rate_moy_rzf+abs(log(1+SINR_rzf));
            
        end

        %% Averaging Rate 
        Rate_moy_pol2_tab(ii)=Rate_moy_pol2/iter;
        Rate_moy_rzf_tab(ii)=Rate_moy_rzf/iter;
    end %ii - Main TxPwr loop

    
    %% Plot
    if ishandle(f1)
        set(0, 'CurrentFigure', f1)
        hold on; grid on;
        plot(10*log10(P*trace(Phi)/M./sigma_tab),Rate_moy_pol2_tab,'ko-','markersize',10,'linewidth',3, 'DisplayName','TPE')
        plot(10*log10(P*trace(Phi)/M./sigma_tab),Rate_moy_rzf_tab,'r+-','markersize',10,'linewidth',3, 'DisplayName','RZF')
        legend('Location','NorthWest')
        legend show
    end
    drawnow   
    
end %secondary outer loop

