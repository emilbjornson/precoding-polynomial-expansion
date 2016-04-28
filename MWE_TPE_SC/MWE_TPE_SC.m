%The definitions and equations refer (mostly) to the 
%following publication:
%
%Müller, A., A. Kammoun, E. Björnson, and M. Debbah, "Linear Precoding 
%Based on Polynomial Expansion: Reducing Complexity in Massive MIMO", 
%EURASIP Journal on Wireless Communications and Networking, 2016(1), 1-22, 
%DOI: 10.1186/s13638-016-0546-z.
%Pre-print available at: http://www.laneas.com/axel-muller
%
%This is version 1.0. (Last edited: 2016-03-16)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.

clear all
disp(' '),disp(' ');
disp('%%%%%%%%%%%%%%%%%%%');
disp('Minimal Working Example: Rate of UT1 vs. power class vs. SNR');
disp(' ');

%% Prepare Figure
figure(1), clf;
hold on; grid on;
xlabel('SNR','fontsize',16)
ylabel('Achievable Avg. Rate','fontsize',16)
set(gca,'FontSize',16)

%% System and Simulation Dimension
Power = 1;
rho_tab = [-10:5:20];
num_class=4;                % Number of different classes for users
K=64;                       % Number of Users
M=256;                      % Number of BS antennas
a=0.1;						% Covariance/Correlation between antennas
tau=0.1;                    % CSI quality
L=4;                        % TPE Order

% Monte-Carlo runs
num_iter=1;                 % Choose >1000 for print quality

%% Power allocation Matrix
P=diag(floor(num_class*linspace(0,1-0.01,K))+1)/(K);    % Example power allocation with proper scaling
% P=diag(ones(1,K)/(K));                                % Equal power allocation

%% Channel Covariance Matrix
Phi=a.^toeplitz([0:M-1]);
sqrtmPhi = sqrtm(Phi); 

%% Preallocate
sinr_th_tab=nan(K,length(rho_tab));
Rate_moy_tab_pol=nan(K,length(rho_tab));
Rate_moy_tab_rzf=nan(K,length(rho_tab));

rate_th_pol=nan(num_class,length(rho_tab));
rate_emp_pol=nan(num_class,length(rho_tab));
rate_emp_rzf=nan(num_class,length(rho_tab));

for index_rho_tab=1:length(rho_tab)
    fprintf('Current SNR = %0.2f dB\n',rho_tab(index_rho_tab));
    rho=rho_tab(index_rho_tab);
    sigma=Power*10^(-rho/10);
    rho_lin=Power/sigma;

    %% Re-Allocate
    Rate_moy_pol=zeros(K,1);
    Rate_moy_rzf=zeros(K,1);


    %% Calculate RZF regularization value
    csi=fixed_point_solve(Phi,tau,rho_lin/K,K);
    % Remark: Asymptotically optimal regularization term only known/valid for
    % equal power allocation.


    %% Calculate TPE Precoder weights (+DE of SINR)
    [gf,sinr_th,A_bar,B_bar,C_bar]=Precoder_polynomial_power_control(Power,P,Phi,L-1,tau,sigma,K);
    sinr_th_tab(:,index_rho_tab)=sinr_th;


    %% Mont-Carlo Loop
%     randn('state',0);
    for iter=1:num_iter
        %% Create Channel draw with error and covariance
        W=(randn(M,K)+1i*randn(M,K))/sqrt(2);
        H=sqrtmPhi*W;
        error_channel=(randn(M,K)+1i*randn(M,K))/sqrt(2);
        H_hat = sqrt(1-tau^2)*H + tau*sqrtm(Phi)*error_channel;

        %% Find TPE Precoding Matrix
        G_pol=zeros(M,K);
        HH=eye(K);
        for l=0:L-1
            G_pol=G_pol+gf(l+1)*H_hat*HH/sqrt(K);
            HH=H_hat'*H_hat*HH/K;
        end
        % Power re-normalization unneeded IFF weights are chosen correctly!

        % Power allocation TPE
        G_pol=G_pol*sqrtm(P);

        %% Find RZF Precoding Matrix
        G_rzf=Precoder_regularized_zero_forcing(Power/trace(P),H_hat,csi);
        % Power allocation RZF
        G_rzf=G_rzf*sqrtm(P);    

        %% Rate of UT1 of each power class
        for k=1:K
            h_1=H(:,k);
            g_1_pol=G_pol(:,k);
            g_1_rzf=G_rzf(:,k);

            SINR_pol=abs( h_1'*g_1_pol*g_1_pol'*h_1/((h_1'*G_pol*G_pol'*h_1-h_1'*g_1_pol*g_1_pol'*h_1+sigma)) );
            SINR_rzf=abs( h_1'*g_1_rzf*g_1_rzf'*h_1/((h_1'*G_rzf*G_rzf'*h_1-h_1'*g_1_rzf*g_1_rzf'*h_1+sigma)) );
            Rate_moy_pol(k)=Rate_moy_pol(k)+log2(1+SINR_pol);
            Rate_moy_rzf(k)=Rate_moy_rzf(k)+log2(1+SINR_rzf);
        end
    end %end iter

    %% Averaging MC iterations
    Rate_moy_pol=Rate_moy_pol/num_iter;
    Rate_moy_tab_pol(:,index_rho_tab)=Rate_moy_pol;
    Rate_moy_rzf=Rate_moy_rzf/num_iter;
    Rate_moy_tab_rzf(:,index_rho_tab)=Rate_moy_rzf;

end % end rho_tab


%% Calculate Rate Expectations
number_of_users_per_class = K/num_class;
for index_rho_tab=1:length(rho_tab)
for ii=1:num_class
    rate_th_pol(ii,index_rho_tab) = mean(log2(1+sinr_th_tab((ii-1)*number_of_users_per_class+1:ii*number_of_users_per_class,index_rho_tab)));
    rate_emp_pol(ii,index_rho_tab)= mean(Rate_moy_tab_pol((ii-1)*+number_of_users_per_class+1:ii*number_of_users_per_class,index_rho_tab));
    rate_emp_rzf(ii,index_rho_tab)= mean(Rate_moy_tab_rzf((ii-1)*+number_of_users_per_class+1:ii*number_of_users_per_class,index_rho_tab));
end
end

%% Plot results
color_tab_class=['r','g','b','k'];
figure(1)
for k=1:num_class
    plot(rho_tab,rate_th_pol(k,:),'LineWidth',2,'color',color_tab_class(k),'Marker','x', 'DisplayName',['DE TPE pwr class ', num2str(k)]);
    hold on
    plot(rho_tab,rate_emp_pol(k,:),'LineWidth',1,'color',color_tab_class(k),'LineStyle','none','Marker','s','MarkerSize',10, 'DisplayName','MC Avg');
    plot(rho_tab,rate_emp_rzf(k,:),'LineWidth',2,'color',color_tab_class(k),'Marker','none','LineStyle',':', 'DisplayName',['RZF pwr class ', num2str(k)]);
end
legend('Location','NorthWest')
legend show
