function G_rzf=Precoder_regularized_zero_forcing(P,H,csi)
%The references to definitions and equations refer (mostly) to the 
%following publication:
%
%Müller, A., A. Kammoun, E. Björnson, and M. Debbah, "Efficient Linear 
%Precoding for Massive MIMO Systems using Truncated Polynomial Expansion", 
%IEEE 8th Sensor Array and Multichannel Signal Processing Workshop (SAM), 
%Coruna, Spain, 2014.
%Available at: http://www.laneas.com/axel-muller
%
%This is version 1.1. (Last edited: 2016-03-14)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.

[M,K]=size(H);
% Q=inv(H*H'+csi*eye(M));
% beta=sqrt(P/(1/K*trace(H'*Q^2*H)));
% G_rzf=beta*inv(H*H'+csi*eye(M))*H;

G_rzf_tmp=inv(H*H'+csi*eye(M))*H;
beta=sqrt(P/(1/K*trace(G_rzf_tmp*G_rzf_tmp')));
G_rzf = beta*G_rzf_tmp;